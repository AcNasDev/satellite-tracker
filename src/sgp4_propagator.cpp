#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
}

void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    elements_.epoch = tle.epoch;

    // Преобразуем элементы орбиты
    elements_.i = tle.inclination * DE2RA;
    elements_.Omega = tle.right_ascension * DE2RA;
    elements_.e = tle.eccentricity;
    elements_.omega = tle.argument_perigee * DE2RA;
    elements_.M = tle.mean_anomaly * DE2RA;
    elements_.n0 = tle.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;
    elements_.bstar = tle.bstar;

    // Вычисляем вспомогательные параметры
    elements_.cosio = cos(elements_.i);
    elements_.sinio = sin(elements_.i);
    elements_.theta2 = elements_.cosio * elements_.cosio;
    elements_.x3thm1 = 3.0 * elements_.theta2 - 1.0;
    elements_.x1mth2 = 1.0 - elements_.theta2;
    elements_.x7thm1 = 7.0 * elements_.theta2 - 1.0;
    elements_.eta = sqrt(1.0 - elements_.e * elements_.e);

    // Вычисляем начальную большую полуось
    double a1 = pow(XKE / elements_.n0, 2.0/3.0);
    double delta1 = 1.5 * CK2 * elements_.x3thm1 /
                    (a1 * a1 * elements_.eta * elements_.eta * elements_.eta);
    elements_.a = a1 * (1.0 - delta1/3.0 - delta1 * delta1 -
                        134.0/81.0 * delta1 * delta1 * delta1);

    // Коррекция среднего движения
    double delta0 = 1.5 * CK2 * elements_.x3thm1 /
                    (elements_.a * elements_.a * elements_.eta * elements_.eta * elements_.eta);
    elements_.n = elements_.n0 / (1.0 + delta0);

    // Коэффициенты для вековых возмущений
    elements_.coef = elements_.bstar * 2.0 * elements_.a * elements_.n0 * XKMPER/AE *
                     (elements_.eta * (1.0 + elements_.eta * cos(elements_.M)));
    elements_.c1 = elements_.coef * elements_.x3thm1;
    elements_.c4 = 2.0 * elements_.n0 * elements_.coef * elements_.a * elements_.eta *
                   (elements_.eta * (2.0 + 0.5 * elements_.eta));
}

void SGP4Propagator::solveKeplerEquation(double& M, double& E) const {
    E = M;
    for(int i = 0; i < 10; i++) {
        double delta = (E - elements_.e * sin(E) - M) / (1.0 - elements_.e * cos(E));
        E -= delta;
        if(fabs(delta) < 1e-12) break;
    }
}

QVector3D SGP4Propagator::calculatePosVel(double tsince) const {
    // Вековые возмущения
    double xmdf = elements_.M + elements_.n * tsince;
    double omgadf = elements_.omega + 0.75 * CK2 * elements_.n * tsince *
                                          elements_.x3thm1 / (elements_.a * elements_.a * elements_.eta);
    double xnode = elements_.Omega - 1.5 * CK2 * elements_.n * tsince *
                                         elements_.cosio / (elements_.a * elements_.a * elements_.eta);

    // Периодические возмущения
    double e = elements_.e;
    double a = elements_.a;
    double xincl = elements_.i;

    // Решаем уравнение Кеплера
    double E;
    solveKeplerEquation(xmdf, E);

    // Вычисляем истинную аномалию
    double sinE = sin(E);
    double cosE = cos(E);
    double nu = atan2(sqrt(1.0 - e * e) * sinE, cosE - e);

    // Расчет расстояния
    double r = a * (1.0 - e * cosE);

    // Вычисляем позицию в орбитальной плоскости
    double u = omgadf + nu;
    double sin_u = sin(u);
    double cos_u = cos(u);
    double sin_xnode = sin(xnode);
    double cos_xnode = cos(xnode);
    double sin_i = sin(xincl);
    double cos_i = cos(xincl);

    // Преобразуем в ECI координаты
    double x = r * (cos_u * cos_xnode - sin_u * cos_i * sin_xnode);
    double y = r * (cos_u * sin_xnode + sin_u * cos_i * cos_xnode);
    double z = r * sin_u * sin_i;

    return QVector3D(x * XKMPER, y * XKMPER, z * XKMPER);
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = elements_.epoch.msecsTo(time) / (1000.0 * 60.0);

    // Вычисляем позицию
    QVector3D pos = calculatePosVel(tsince);
    state.position = pos;

    // Вычисляем скорость через численное дифференцирование
    double dt = 0.1; // 0.1 минуты
    QVector3D pos2 = calculatePosVel(tsince + dt);
    state.velocity = (pos2 - pos) / (dt * 60.0); // переводим в км/с

    return state;
}
