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
    elements_.x2o3 = TOTHRD;

    // Вычисление большой полуоси
    double a1 = pow(XKE / elements_.n0, elements_.x2o3);
    double delta1 = (3.0 * elements_.x3thm1 / (4.0 * a1 * a1)) * CK2;
    elements_.a = a1 * (1.0 - delta1 * (1.0/3.0 + delta1 * (1.0 + 134.0/81.0 * delta1)));
    elements_.eta = sqrt(1.0 - elements_.e * elements_.e);

    // Долгопериодические возмущения
    elements_.xlcof = 0.125 * CK2 * elements_.sinio *
                      (3.0 + 5.0 * elements_.cosio) / (1.0 + elements_.cosio);
    elements_.aycof = 0.25 * CK2 * elements_.sinio;

    // Атмосферное сопротивление
    double temp = (elements_.eta * elements_.eta *
                   (1.0 - 1.5 * elements_.theta2 +
                    2.0 * elements_.cosio * elements_.cosio));
    elements_.coef = 2.0 * elements_.bstar * CK2 * temp /
                     (elements_.a * (1.0 - elements_.e * elements_.e));

    // Вековые возмущения
    elements_.xmdot = elements_.n0 + 0.5 * elements_.coef * elements_.x3thm1;
    elements_.omgdot = -0.5 * CK2 * (5.0 * elements_.cosio * elements_.cosio - 1.0) /
                       (elements_.a * elements_.a * elements_.eta);
    elements_.xnodot = -CK2 * elements_.cosio /
                       (elements_.a * elements_.a * elements_.eta);

    // Дополнительные параметры
    elements_.c1 = elements_.coef * elements_.bstar * CK2 * elements_.x3thm1;
    elements_.c4 = 2.0 * elements_.n0 * elements_.coef * elements_.a * elements_.eta *
                   ((2.0 + elements_.eta) * (1.0 + elements_.e * cos(elements_.M)));
    elements_.delmo = pow(1.0 + elements_.eta * cos(elements_.M), 3);
    elements_.sinmo = sin(elements_.M);
    elements_.aodp = elements_.a;
}

void SGP4Propagator::calculatePerturbations(double tsince, double& ep, double& xincp,
                                            double& omgadf, double& xnode, double& xmp) const {
    // Учет вековых возмущений
    double xmdf = elements_.M + elements_.xmdot * tsince;
    omgadf = elements_.omega + elements_.omgdot * tsince;
    xnode = elements_.Omega + elements_.xnodot * tsince;
    double omega_k = omgadf;
    ep = elements_.e - elements_.coef * tsince;
    xincp = elements_.i;
    xmp = xmdf;

    // Долгопериодические возмущения
    double axn = ep * cos(omega_k);
    double temp = 1.0 / (elements_.a * (1.0 - ep * ep));
    double xlcof = 0.125 * CK2 * elements_.sinio *
                   (3.0 + 5.0 * elements_.cosio) /
                   (1.0 + elements_.cosio);
    double aycof = 0.25 * CK2 * elements_.sinio;

    double xll = temp * xlcof * axn;
    double aynl = temp * aycof;
    double xlt = xll + xnode;
    double ayn = ep * sin(omega_k) + aynl;

    omega_k = omgadf + (0.5 * temp * CK2 * elements_.x3thm1) * tsince;
    xnode = xnode + xll;

    // Короткопериодические возмущения
    double elsq = axn * axn + ayn * ayn;
    if (elsq >= 1.0) {
        ep = sqrt(elsq);
    }

    double sin_xlt = sin(xlt);
    double cos_xlt = cos(xlt);

    omgadf = omega_k;
    xmp = xmdf + elements_.xlcof * axn * cos_xlt + elements_.aycof * ayn * sin_xlt;
}

void SGP4Propagator::solveKeplerEquation(double xme, double xec, double& eps, double& epw) const {
    epw = xme;
    for(int i = 0; i < 10; i++) {
        eps = epw - xec * sin(epw) - xme;
        double newton_raphson = eps / (1.0 - xec * cos(epw));
        epw = epw - newton_raphson;
        if(fabs(eps) < 1.0e-12) break;
    }
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = elements_.epoch.msecsTo(time) / (1000.0 * 60.0);

    // Вычисляем возмущения
    double ep, xincp, omgadf, xnode, xmp;
    calculatePerturbations(tsince, ep, xincp, omgadf, xnode, xmp);

    // Решаем уравнение Кеплера
    double eps = 0.0, epw = 0.0;
    solveKeplerEquation(xmp, ep, eps, epw);

    // Вычисляем положение в орбитальной плоскости
    double sin_epw = sin(epw);
    double cos_epw = cos(epw);

    double ecosE = 1.0 - ep * cos_epw;
    double esinE = ep * sin_epw;
    double elsq = ep * ep;
    double pl = elements_.a * (1.0 - elsq);

    double r = elements_.a * ecosE;
    double rdot = -elements_.n0 * elements_.a * esinE / ecosE;
    double rvdot = elements_.n0 * elements_.a / r;

    // Преобразование в ECI
    double sin_xnode = sin(xnode);
    double cos_xnode = cos(xnode);
    double sin_xinc = sin(xincp);
    double cos_xinc = cos(xincp);
    double sin_u = sin(epw + omgadf);
    double cos_u = cos(epw + omgadf);

    // Позиция
    double x = r * (cos_u * cos_xnode - sin_u * cos_xinc * sin_xnode);
    double y = r * (cos_u * sin_xnode + sin_u * cos_xinc * cos_xnode);
    double z = r * sin_u * sin_xinc;

    // Скорость
    double vx = rdot * (cos_u * cos_xnode - sin_u * cos_xinc * sin_xnode) -
                rvdot * (sin_u * cos_xnode + cos_u * cos_xinc * sin_xnode);
    double vy = rdot * (cos_u * sin_xnode + sin_u * cos_xinc * cos_xnode) -
                rvdot * (sin_u * sin_xnode - cos_u * cos_xinc * cos_xnode);
    double vz = rdot * sin_u * sin_xinc + rvdot * cos_u * sin_xinc;

    // Перевод в километры и километры в секунду
    state.position = QVector3D(x * XKMPER, y * XKMPER, z * XKMPER);
    state.velocity = QVector3D(vx * XKMPER / 60.0,
                               vy * XKMPER / 60.0,
                               vz * XKMPER / 60.0);

    return state;
}
