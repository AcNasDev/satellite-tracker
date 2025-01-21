#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
}

void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    elements_.epoch = tle.epoch;

    // Базовые элементы
    elements_.i = tle.inclination * DE2RA;
    elements_.Omega = tle.right_ascension * DE2RA;
    elements_.e = tle.eccentricity;
    elements_.omega = tle.argument_perigee * DE2RA;
    elements_.M = tle.mean_anomaly * DE2RA;
    elements_.n0 = tle.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;
    elements_.bstar = tle.bstar;

    // Вспомогательные элементы
    elements_.cosio = cos(elements_.i);
    elements_.sinio = sin(elements_.i);
    elements_.theta2 = elements_.cosio * elements_.cosio;
    elements_.x3thm1 = 3.0 * elements_.theta2 - 1.0;
    elements_.x7thm1 = 7.0 * elements_.theta2 - 1.0;
    elements_.beta0 = sqrt(1.0 - elements_.e * elements_.e);
    elements_.eta = elements_.beta0;

    // Коррекция для SGP4
    double a1 = pow(XKE / elements_.n0, 2.0/3.0);
    double delta1 = (3.0/2.0) * CK2 * elements_.x3thm1 /
                    (a1 * a1 * elements_.beta0 * elements_.beta0 * elements_.beta0);

    elements_.a0 = a1 * (1.0 - delta1/3.0 - delta1 * delta1 -
                         134.0/81.0 * delta1 * delta1 * delta1);

    double delta0 = (3.0/2.0) * CK2 * elements_.x3thm1 /
                    (elements_.a0 * elements_.a0 * elements_.beta0 * elements_.beta0 * elements_.beta0);

    // Среднее движение с коррекцией
    elements_.n = elements_.n0 / (1.0 + delta0);
    elements_.a = elements_.a0 / (1.0 - delta0);

    // Коэффициенты для вековых возмущений
    double theta4 = elements_.theta2 * elements_.theta2;
    elements_.coef = 2.0 * elements_.bstar * CK2 * elements_.n /
                     (elements_.beta0 * elements_.beta0 * elements_.beta0);

    elements_.c1 = elements_.coef * elements_.x3thm1;
    elements_.c4 = 2.0 * elements_.n * CK2 * elements_.x3thm1 * elements_.a * elements_.beta0 /
                   (elements_.eta * elements_.eta * elements_.eta);
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = elements_.epoch.msecsTo(time) / (1000.0 * 60.0);

    // Вековые эффекты
    double xmdf = elements_.M + elements_.n * tsince;
    double omgadf = elements_.omega;
    double xnode = elements_.Omega;
    double eps1 = 0.0;
    calculateSecularEffects(tsince, xmdf, omgadf, xnode, eps1);

    // Периодические эффекты
    double em = elements_.e;
    double xinc = elements_.i;
    double xll = xmdf;
    calculatePeriodicEffects(tsince, em, xinc, omgadf, xnode, xll);

    // Решение уравнения Кеплера
    double u = xll - xnode;
    double eo1 = u;
    double tem5 = 1.0;

    // Итерационное решение
    for(int i = 0; i < 10 && fabs(tem5) >= 1e-12; i++) {
        double sineo1 = sin(eo1);
        double coseo1 = cos(eo1);
        tem5 = 1.0 - coseo1 * em - sineo1 * em - u;
        tem5 = tem5 / (sineo1 * em - coseo1 * em) /
               (1.0 - coseo1 * em - sineo1 * em);
        eo1 = eo1 - tem5;
    }

    // Расчёт позиции и скорости в орбитальной плоскости
    double coseo1 = cos(eo1);
    double sineo1 = sin(eo1);
    double rk = elements_.a * (1.0 - em * coseo1);
    double rdot = elements_.a * em * sineo1 * elements_.n / rk;
    double rfdot = elements_.a * sqrt(1.0 - em * em) * coseo1 * elements_.n / rk;

    // Ориентация орбитальной плоскости
    double uk = u + omgadf;
    double cosuk = cos(uk);
    double sinuk = sin(uk);

    // Расчёт позиции и скорости в ECI
    QVector3D pos_vel = calculatePosVel(rk, rdot, rfdot, cosuk, sinuk, rk, uk, xnode, xinc);

    state.position = QVector3D(pos_vel.x(), pos_vel.y(), pos_vel.z());

    // Скорость считаем через конечные разности
    const double dt = 0.0001; // малый шаг времени для численного дифференцирования
    OrbitalState state2 = calculateState(time.addMSecs(dt * 1000));
    state.velocity = (state2.position - state.position) / dt;

    return state;
}

void SGP4Propagator::calculateSecularEffects(double tsince, double& xmdf, double& omgadf,
                                             double& xnode, double& eps1) const {
    double xhdot1 = -elements_.c1 * tsince;
    eps1 = xhdot1 * elements_.x3thm1 / 2.0;

    // Вековые возмущения для e, i и Ω
    omgadf += (2.0 * elements_.c1 + 3.0 * elements_.c1 * elements_.cosio) * tsince / 2.0;
    xnode += (-3.0 * elements_.c1 * elements_.sinio) * tsince / 2.0;
    xmdf += elements_.n * tsince + eps1;
}

void SGP4Propagator::calculatePeriodicEffects(double tsince, double& em, double& xinc,
                                              double& omgda, double& xnode, double& xll) const {
    // Периодические возмущения
    double sin2u = sin(2.0 * omgda);
    double cos2u = cos(2.0 * omgda);

    double temp1 = CK2 * (1.0 - 1.5 * elements_.theta2);
    double temp2 = (elements_.a * (1.0 - elements_.e * elements_.e)) /
                   (1.0 + elements_.beta0);

    // Коррекции для e и i
    em += elements_.c4 * tsince * sin2u;
    xinc += -0.5 * temp1 * elements_.sinio * cos2u * tsince;

    // Коррекция для Ω и ω
    double xlamb = xll + omgda + xnode;
    double temp3 = 2.0 * (xlamb - elements_.Omega);
    double sin_temp3 = sin(temp3);
    double cos_temp3 = cos(temp3);

    omgda += 0.5 * temp1 * temp2 * (3.0 * cos2u - 1.0) -
             0.25 * CK2 * temp2 * elements_.x3thm1 * sin2u;
    xnode += -0.5 * temp1 * elements_.cosio * sin2u * tsince;
}

QVector3D SGP4Propagator::calculatePosVel(double r, double rdot, double rfdot,
                                          double cosuk, double sinuk, double rk,
                                          double uk, double xnodek, double xinck) const {
    double xmx = -sin(xnodek) * cos(xinck);
    double xmy = cos(xnodek) * cos(xinck);
    double ux = xmx * sinuk + cos(xnodek) * cosuk;
    double uy = xmy * sinuk + sin(xnodek) * cosuk;
    double uz = sin(xinck) * sinuk;

    double x = rk * ux * XKMPER;
    double y = rk * uy * XKMPER;
    double z = rk * uz * XKMPER;

    return QVector3D(x, y, z);
}
