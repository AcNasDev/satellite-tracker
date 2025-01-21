#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
}

void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    elements_.epoch = tle.epoch;

    // Основные элементы
    elements_.i = tle.inclination * DE2RA;
    elements_.Omega = tle.right_ascension * DE2RA;
    elements_.e = tle.eccentricity;
    elements_.omega = tle.argument_perigee * DE2RA;
    elements_.M = tle.mean_anomaly * DE2RA;
    elements_.n0 = tle.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;
    elements_.bstar = tle.bstar;

    // Вспомогательные параметры
    elements_.cosio = cos(elements_.i);
    elements_.sinio = sin(elements_.i);
    elements_.theta2 = elements_.cosio * elements_.cosio;
    elements_.x3thm1 = 3.0 * elements_.theta2 - 1.0;
    elements_.x1mth2 = 1.0 - elements_.theta2;
    elements_.x7thm1 = 7.0 * elements_.theta2 - 1.0;
    elements_.beta0 = sqrt(1.0 - elements_.e * elements_.e);
    elements_.eta = elements_.beta0;

    // Вычисление большой полуоси
    double a1 = pow(XKE / elements_.n0, TOTHRD);
    double delta1 = 1.5 * CK2 * elements_.x3thm1 /
                    (a1 * a1 * elements_.beta0 * elements_.beta0 * elements_.beta0);

    double a0 = a1 * (1.0 - delta1/3.0 - delta1 * delta1 -
                      134.0/81.0 * delta1 * delta1 * delta1);

    double delta0 = 1.5 * CK2 * elements_.x3thm1 /
                    (a0 * a0 * elements_.beta0 * elements_.beta0 * elements_.beta0);

    elements_.a = a0 / (1.0 - delta0);
    elements_.n = elements_.n0 / (1.0 + delta0);

    // Вычисление скоростей изменения элементов орбиты
    elements_.xhdot = -1.5 * CK2 * elements_.n * elements_.cosio /
                      (elements_.a * elements_.a * elements_.beta0);

    elements_.xndot = elements_.n - elements_.n0;

    elements_.xnodot = elements_.xhdot;

    elements_.omgdot = -0.75 * CK2 * elements_.n *
                       (2.0 - 3.0 * elements_.theta2) /
                       (elements_.a * elements_.a * elements_.beta0);

    elements_.xmdot = elements_.n + elements_.xndot;

    // Коэффициенты для возмущений
    elements_.coef = elements_.bstar * 2.0 * elements_.a *
                     elements_.n * XKMPER/AE * elements_.eta;

    elements_.c1 = elements_.coef * elements_.x3thm1;

    elements_.c4 = 2.0 * elements_.n * elements_.coef * elements_.a *
                   elements_.eta * (2.0 + elements_.eta * elements_.eta);
}

void SGP4Propagator::updateForTime(double tsince, double& xll, double& omgadf,
                                   double& xnode, double& em, double& xinc,
                                   double& xn) const {
    // Вековые возмущения
    xll = elements_.M + (elements_.xmdot + elements_.xnodot) * tsince;
    omgadf = elements_.omega + elements_.omgdot * tsince;
    xnode = elements_.Omega + elements_.xnodot * tsince;
    em = elements_.e - elements_.coef * tsince;
    xinc = elements_.i;
    xn = elements_.n;

    // Периодические возмущения
    double axn = em * cos(omgadf);
    double temp = 1.0 / (elements_.a * (1.0 - em * em));

    // Возмущения в долготе
    xll += elements_.c1 * temp * axn * tsince;

    // Возмущения в эксцентриситете
    em += -elements_.c4 * temp * axn;
}

void SGP4Propagator::solveKeplerEquation(double& xll, double& e) const {
    double epw = xll;

    for(int i = 0; i < 10; i++) {
        double sin_epw = sin(epw);
        double cos_epw = cos(epw);
        double delta = (epw - e * sin_epw - xll) / (1.0 - e * cos_epw);
        epw -= delta;

        if(fabs(delta) <= E6A) break;
    }

    xll = epw;
}

QVector3D SGP4Propagator::calculatePosVel(double tsince) const {
    // Инициализация переменных
    double xll = elements_.M;
    double omgadf = elements_.omega;
    double xnode = elements_.Omega;
    double em = elements_.e;
    double xinc = elements_.i;
    double xn = elements_.n;

    // Обновление элементов с учетом возмущений
    updateForTime(tsince, xll, omgadf, xnode, em, xinc, xn);

    // Решение уравнения Кеплера
    solveKeplerEquation(xll, em);

    // Вычисление координат в орбитальной плоскости
    double sin_xll = sin(xll);
    double cos_xll = cos(xll);

    double r = elements_.a * (1.0 - em * cos_xll);

    // Преобразование в ECI
    double sin_omgadf = sin(omgadf);
    double cos_omgadf = cos(omgadf);
    double sin_xinc = sin(xinc);
    double cos_xinc = cos(xinc);
    double sin_xnode = sin(xnode);
    double cos_xnode = cos(xnode);

    double x = r * ((cos_xll * cos_omgadf - sin_xll * sin_omgadf) * cos_xnode -
                    (cos_xll * sin_omgadf + sin_xll * cos_omgadf) * cos_xinc * sin_xnode);

    double y = r * ((cos_xll * cos_omgadf - sin_xll * sin_omgadf) * sin_xnode +
                    (cos_xll * sin_omgadf + sin_xll * cos_omgadf) * cos_xinc * cos_xnode);

    double z = r * (cos_xll * sin_omgadf + sin_xll * cos_omgadf) * sin_xinc;

    return QVector3D(x * XKMPER, y * XKMPER, z * XKMPER);
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = elements_.epoch.msecsTo(time) / (1000.0 * 60.0);

    // Расчет положения
    state.position = calculatePosVel(tsince);

    // Вычисление скорости через конечные разности
    double dt = 0.001; // 0.001 минуты для лучшей точности
    QVector3D pos2 = calculatePosVel(tsince + dt);
    state.velocity = (pos2 - state.position) / (dt * 60.0);

    return state;
}
