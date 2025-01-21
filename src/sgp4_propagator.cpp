#include "sgp4_propagator.h"
#include <QtMath>
#include <cmath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
}

void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    // Базовые элементы
    sgp4_.epoch = tle.epoch;
    sgp4_.i = tle.inclination * DE2RA;
    sgp4_.Omega = tle.right_ascension * DE2RA;
    sgp4_.e = tle.eccentricity;
    sgp4_.omega = tle.argument_perigee * DE2RA;
    sgp4_.M = tle.mean_anomaly * DE2RA;
    sgp4_.n0 = tle.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;
    sgp4_.bstar = tle.bstar;

    // Вычисление производных элементов
    sgp4_.cosio = cos(sgp4_.i);
    sgp4_.sinio = sin(sgp4_.i);
    sgp4_.x3thm1 = 3.0 * sgp4_.cosio * sgp4_.cosio - 1.0;
    sgp4_.x1mth2 = 1.0 - sgp4_.cosio * sgp4_.cosio;
    sgp4_.x7thm1 = 7.0 * sgp4_.cosio * sgp4_.cosio - 1.0;

    // Коррекция для SGP4
    double beta0 = sqrt(1.0 - sgp4_.e * sgp4_.e);
    sgp4_.eta = sqrt(1.0 - sgp4_.e * sgp4_.e);

    // Вычисление a1 и коррекция mean motion
    double a1 = pow(XKE / sgp4_.n0, 2.0/3.0);
    double delta1 = 1.5 * CK2 * sgp4_.x3thm1 / (a1 * a1 * beta0 * beta0 * beta0);
    sgp4_.a0 = a1 * (1.0 - delta1 * (1.0/3.0 + delta1 * (1.0 + 134.0/81.0 * delta1)));
    double delta0 = 1.5 * CK2 * sgp4_.x3thm1 / (sgp4_.a0 * sgp4_.a0 * beta0 * beta0 * beta0);

    // Обновление mean motion с учетом возмущений
    sgp4_.n = sgp4_.n0 / (1.0 + delta0);
    sgp4_.a = sgp4_.a0 * pow((sgp4_.n0/sgp4_.n), 2.0/3.0);

    // Коэффициенты для long-period periodic effects
    sgp4_.c1 = CK2 * sgp4_.x3thm1;
    sgp4_.c4 = 2.0 * sgp4_.n * CK2 * sgp4_.x3thm1 *
               sgp4_.a * beta0 * beta0 / (sgp4_.eta * sgp4_.eta * sgp4_.eta);

    // Расчет производных mean motion
    sgp4_.ndot = -2.0 * (1.0/3.0) * delta0 * sgp4_.n * (1.0 + 4.0 * delta0);
    sgp4_.nddot = -delta0 * delta0 * sgp4_.n * (15.0/6.0 + 17.0 * delta0);
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = sgp4_.epoch.msecsTo(time) / (1000.0 * 60.0);

    // Расчет вековых возмущений
    double xmdf = sgp4_.M + sgp4_.n * tsince;
    double omgadf = sgp4_.omega + sgp4_.ndot * tsince;
    double xnode = sgp4_.Omega + sgp4_.nddot * tsince * tsince / 2.0;

    // Учет долгопериодических возмущений
    calculateSecularEffects(tsince, xmdf, omgadf, xnode);

    // Расчет периодических возмущений
    double ep = sgp4_.e;
    double xincp = sgp4_.i;
    calculatePeriodicEffects(tsince, ep, xincp, omgadf, xnode);

    // Решение уравнения Кеплера
    double xl = xmdf + omgadf + xnode;
    double u = fmod(xl - xnode, 2.0 * M_PI);
    double eo1 = u;

    double tem5 = 1.0;
    int iteration = 0;
    while (fabs(tem5) >= 1e-12 && iteration < 10) {
        double sineo1 = sin(eo1);
        double coseo1 = cos(eo1);
        tem5 = 1.0 - coseo1 * ep - sineo1 * sgp4_.e - u;
        tem5 = tem5 / (sineo1 * ep - coseo1 * sgp4_.e) *
               (1.0 / (1.0 - coseo1 * ep - sineo1 * sgp4_.e));
        eo1 = eo1 - tem5;
        iteration++;
    }

    // Расчет позиции и скорости
    double sinu = sin(u);
    double cosu = cos(u);

    double rk = sgp4_.a * (1.0 - ep * cos(eo1));
    double uk = u + omgadf;
    double xnodek = xnode;
    double xinck = xincp;

    double rdotk = sgp4_.a * ep * sin(eo1) * sgp4_.n / rk;
    double rfdotk = sgp4_.a * sqrt(1.0 - ep * ep) * cos(eo1) * sgp4_.n / rk;

    QVector3D pos_vel = calculatePositionAndVelocity(rk, uk, xnodek, xinck, rdotk, rfdotk);
    state.position = QVector3D(pos_vel.x(), pos_vel.y(), 0);
    state.velocity = QVector3D(0, 0, pos_vel.z());

    return state;
}

void SGP4Propagator::calculateSecularEffects(double tsince, double& xmdf,
                                             double& omgadf, double& xnode) const {
    double xmp = xmdf + sgp4_.c1 * tsince *
                            (1.0 - 3.0 * sgp4_.cosio * sgp4_.cosio) * 0.5;
    double omega_df = omgadf - sgp4_.c1 * tsince *
                                   (1.0 - 5.0 * sgp4_.cosio * sgp4_.cosio) * 0.25;
    xnode = xnode + sgp4_.c1 * tsince *
                        (4.0 - 19.0 * sgp4_.cosio * sgp4_.cosio) * 0.25;

    xmdf = xmp;
    omgadf = omega_df;
}

void SGP4Propagator::calculatePeriodicEffects(double tsince, double& ep, double& xincp,
                                              double& omgadf, double& xnode) const {
    double sin_omega = sin(omgadf);
    double cos_omega = cos(omgadf);

    // Short-period периодические возмущения
    double dx = sgp4_.c4 * sin(2.0 * omgadf);
    double dy = -sgp4_.c4 * cos(2.0 * omgadf);

    ep = sgp4_.e + dx;
    xincp = sgp4_.i + dy;
}

QVector3D SGP4Propagator::calculatePositionAndVelocity(double rk, double uk,
                                                       double xnodek, double xinck,
                                                       double rdotk, double rfdotk) const {
    double sinuk = sin(uk);
    double cosuk = cos(uk);
    double sinik = sin(xinck);
    double cosik = cos(xinck);
    double sinnok = sin(xnodek);
    double cosnok = cos(xnodek);

    double xmx = -sinnok * cosik;
    double xmy = cosnok * cosik;
    double ux = xmx * sinuk + cosnok * cosuk;
    double uy = xmy * sinuk + sinnok * cosuk;
    double uz = sinik * sinuk;

    // Позиция в км
    double x = rk * ux * XKMPER;
    double y = rk * uy * XKMPER;
    double z = rk * uz * XKMPER;

    // Скорость в км/с
    double vx = (rdotk * ux + rfdotk * (-xmx * cosuk + cosnok * sinuk)) * XKMPER / 60.0;
    double vy = (rdotk * uy + rfdotk * (-xmy * cosuk + sinnok * sinuk)) * XKMPER / 60.0;
    double vz = (rdotk * uz + rfdotk * (sinik * cosuk)) * XKMPER / 60.0;

    return QVector3D(x, y, z);
}
