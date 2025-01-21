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
    elements_.eta = sqrt(1.0 - elements_.e * elements_.e);
    elements_.beta0 = elements_.eta;
    elements_.x3thm1 = 3.0 * elements_.cosio * elements_.cosio - 1.0;
    elements_.x1mth2 = 1.0 - elements_.cosio * elements_.cosio;
    elements_.e0 = elements_.e;

    // Расчет большой полуоси
    double temp = 0.75 * J2 * elements_.x3thm1;
    double a1 = pow(XKE / elements_.n0, 2.0/3.0);
    double del1 = temp / (a1 * a1 * elements_.beta0 * elements_.beta0 * elements_.beta0);
    double ao = a1 * (1.0 - del1 * (0.5 * 2.0/3.0 + del1 * (1.0 + 134.0/81.0 * del1)));
    double delo = temp / (ao * ao * elements_.beta0 * elements_.beta0 * elements_.beta0);

    elements_.aodp = ao;
    elements_.a = ao / (1.0 - delo);
    elements_.n = elements_.n0;

    // Коэффициенты для возмущений
    double theta2 = elements_.cosio * elements_.cosio;
    elements_.xmdot = elements_.n + 0.5 * temp * elements_.beta0 * elements_.x3thm1 / elements_.aodp;
    double xhdot1 = -temp * elements_.cosio / elements_.aodp;

    elements_.omgdot = -0.5 * temp * (1.0 - 5.0 * theta2) / (elements_.aodp * elements_.beta0);
    elements_.xnodot = xhdot1;

    // Коэффициенты для атмосферного сопротивления
    elements_.coef = 2.0 * elements_.bstar * CK2 * elements_.aodp / elements_.beta0;
    elements_.C1 = elements_.coef * elements_.x3thm1;

    // Коэффициенты для короткопериодических возмущений
    elements_.C4 = 2.0 * elements_.n * elements_.coef * elements_.aodp * elements_.beta0 *
                   (elements_.eta * (2.0 + 0.5 * elements_.eta));

    elements_.t2cof = 1.5 * elements_.C1;
    elements_.xnodcf = 3.5 * elements_.beta0 * elements_.C1;
}

void SGP4Propagator::updateForSecularEffects(double tsince, double& xll,
                                             double& omgasm, double& xnodes,
                                             double& em, double& xinc, double& xn) const {
    // Вековые возмущения
    xll = elements_.M + elements_.xmdot * tsince;
    omgasm = elements_.omega + elements_.omgdot * tsince;
    xnodes = elements_.Omega + elements_.xnodot * tsince;
    em = elements_.e - elements_.bstar * elements_.C4 * tsince;
    xinc = elements_.i;
    xn = elements_.n;

    // Учет атмосферного сопротивления
    double temp = 1.0 - elements_.C1 * tsince;
    temp *= temp;
    xn = elements_.n0 / temp;
}

void SGP4Propagator::updateForPeriodicEffects(double& em, double& xinc,
                                              double& omgasm, double& xnodes,
                                              double& xll, double tsince) const {
    // Короткопериодические возмущения
    double axn = em * cos(omgasm);
    double temp = 1.0 / (em * (1.0 - em * em));
    double xlcof = 0.125 * elements_.x3thm1 * temp * elements_.C1;
    double aycof = 0.25 * elements_.sinio * temp * elements_.C1;

    double xmam = xll + elements_.omega;
    xll += xlcof * axn;
    omgasm -= aycof * (xmam + omgasm - xnodes);
}

void SGP4Propagator::solveKeplerEquation(double& xll, double e) const {
    double epw = xll;
    for(int i = 0; i < 10; i++) {
        double sin_epw = sin(epw);
        double cos_epw = cos(epw);
        double delta = (epw - e * sin_epw - xll) / (1.0 - e * cos_epw);
        epw -= delta;
        if(fabs(delta) < 1e-12) break;
    }
    xll = epw;
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = elements_.epoch.msecsTo(time) / (1000.0 * 60.0);

    // Вычисление вековых эффектов
    double xll, omgasm, xnodes, em, xinc, xn;
    updateForSecularEffects(tsince, xll, omgasm, xnodes, em, xinc, xn);

    // Вычисление периодических эффектов
    updateForPeriodicEffects(em, xinc, omgasm, xnodes, xll, tsince);

    // Решение уравнения Кеплера
    solveKeplerEquation(xll, em);

    // Вычисление положения в орбитальной плоскости
    double sin_xll = sin(xll);
    double cos_xll = cos(xll);

    double r = elements_.a * (1.0 - em * cos_xll);
    double rdot = xn * elements_.a * em * sin_xll / sqrt(1.0 - em * em);
    double rfdot = xn * elements_.a * sqrt(1.0 - em * em) / r;

    // Ориентация в пространстве
    double sin_u = sin(omgasm);
    double cos_u = cos(omgasm);
    double sin_Omega = sin(xnodes);
    double cos_Omega = cos(xnodes);
    double sin_i = sin(xinc);
    double cos_i = cos(xinc);

    // Позиция в км
    double x = r * (cos_u * cos_Omega - sin_u * sin_Omega * cos_i);
    double y = r * (cos_u * sin_Omega + sin_u * cos_Omega * cos_i);
    double z = r * sin_u * sin_i;

    state.position = QVector3D(x * XKMPER, y * XKMPER, z * XKMPER);

    // Скорость в км/с
    double vx = (rdot * (cos_u * cos_Omega - sin_u * sin_Omega * cos_i) -
                 r * (sin_u * cos_Omega + cos_u * sin_Omega * cos_i) * rfdot) * XKMPER / 60.0;

    double vy = (rdot * (cos_u * sin_Omega + sin_u * cos_Omega * cos_i) +
                 r * (-sin_u * sin_Omega + cos_u * cos_Omega * cos_i) * rfdot) * XKMPER / 60.0;

    double vz = (rdot * sin_u * sin_i +
                 r * cos_u * sin_i * rfdot) * XKMPER / 60.0;

    state.velocity = QVector3D(vx, vy, vz);

    return state;
}
