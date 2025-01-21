// sgp4_propagator.cpp
#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initParameters(tle);
    calculateDerivatives();
}

void SGP4Propagator::initParameters(const TLEParser::TLEData& tle) {
    epoch_ = tle.epoch;

    qDebug() << "\n=== SGP4 Initialization ===";

    // Преобразование основных элементов
    elements_.inclo = tle.inclination * de2ra;
    elements_.nodeo = tle.right_ascension * de2ra;
    elements_.ecco = tle.eccentricity;
    elements_.argpo = tle.argument_perigee * de2ra;
    elements_.mo = tle.mean_anomaly * de2ra;
    elements_.no = tle.mean_motion * 2.0 * M_PI / min_per_day;
    elements_.bstar = tle.bstar;

    qDebug() << "Initial TLE elements:";
    qDebug() << "Inclination (rad):" << elements_.inclo;
    qDebug() << "RAAN (rad):" << elements_.nodeo;
    qDebug() << "Eccentricity:" << elements_.ecco;
    qDebug() << "Arg of perigee (rad):" << elements_.argpo;
    qDebug() << "Mean anomaly (rad):" << elements_.mo;
    qDebug() << "Mean motion (rad/min):" << elements_.no;
    qDebug() << "BSTAR:" << elements_.bstar;

    // Вспомогательные величины
    double theta2 = cos(elements_.inclo) * cos(elements_.inclo);
    double x3thm1 = 3.0 * theta2 - 1.0;
    double eosq = elements_.ecco * elements_.ecco;
    double betao2 = 1.0 - eosq;
    double betao = sqrt(betao2);
    double temp = 1.5 * ck2 * x3thm1 / (betao * betao2);

    // Восстановление оригинального среднего движения
    double a1 = pow(xke / elements_.no, 2.0/3.0);
    elements_.del1 = temp / (a1 * a1);
    double ao = a1 * (1.0 - elements_.del1 * (1.0/3.0 + elements_.del1 *
                                                              (1.0 + 134.0/81.0 * elements_.del1)));
    elements_.del2 = temp / (ao * ao);
    elements_.no = elements_.no / (1.0 + elements_.del1);
    elements_.a = ao / (1.0 - elements_.del1);

    // Вычисление периода
    elements_.alta = elements_.a * (1.0 + elements_.ecco) - 1.0;
    elements_.altp = elements_.a * (1.0 - elements_.ecco) - 1.0;

    qDebug() << "\nDerived elements:";
    qDebug() << "Semi-major axis (ER):" << elements_.a;
    qDebug() << "Corrected mean motion (rad/min):" << elements_.no;
    qDebug() << "Apogee height (ER):" << elements_.alta;
    qDebug() << "Perigee height (ER):" << elements_.altp;
}

void SGP4Propagator::calculateDerivatives() {
    double theta2 = cos(elements_.inclo) * cos(elements_.inclo);
    double x3thm1 = 3.0 * theta2 - 1.0;
    double eosq = elements_.ecco * elements_.ecco;
    double betao2 = 1.0 - eosq;
    double betao = sqrt(betao2);
    double temp1 = 1.5 * ck2 * x3thm1 / (betao * betao2);
    double temp2 = temp1 / elements_.a;

    elements_.ndot = -temp2 * elements_.no * (1.0 + 1.5 * elements_.ecco);
    elements_.nddot = -temp2 * elements_.ndot;

    qDebug() << "Motion derivatives:";
    qDebug() << "ndot (rad/min²):" << elements_.ndot;
    qDebug() << "nddot (rad/min³):" << elements_.nddot;
}

double SGP4Propagator::solveKepler(double mean_anomaly, double ecc) const {
    double E = mean_anomaly;
    for(int i = 0; i < 10; i++) {
        double E_old = E;
        E = mean_anomaly + ecc * sin(E);
        if(fabs(E - E_old) < 1.0e-12) break;
    }
    return E;
}

void SGP4Propagator::propagate(double tsince, QVector3D& pos, QVector3D& vel) {
    qDebug() << "\nPropagating for tsince =" << tsince << "minutes";

    // Учет вековых возмущений
    double xn = elements_.no + elements_.ndot * tsince +
                elements_.nddot * tsince * tsince / 2.0;
    elements_.e = elements_.ecco - elements_.bstar * tsince;
    elements_.e = qMax(1.0e-6, qMin(0.999999, elements_.e));

    double xmp = elements_.mo + xn * tsince;
    double omega = elements_.argpo + 0.75 * ck2 * sin(elements_.inclo) *
                                         xn * tsince / (elements_.no * betao2);
    double xnode = elements_.nodeo - 1.5 * ck2 * cos(elements_.inclo) *
                                         xn * tsince / (elements_.no * betao2);

    qDebug() << "Secular perturbations:";
    qDebug() << "Mean motion (rad/min):" << xn;
    qDebug() << "Eccentricity:" << elements_.e;
    qDebug() << "Mean anomaly (rad):" << xmp;
    qDebug() << "Arg of perigee (rad):" << omega;
    qDebug() << "RAAN (rad):" << xnode;

    // Решение уравнения Кеплера
    double xl = xmp + omega;
    double E = solveKepler(xl - xnode, elements_.e);

    // Вычисление позиции в орбитальной плоскости
    double sinu = sin(E);
    double cosu = cos(E);
    double r = elements_.a * (1.0 - elements_.e * cosu);
    double u = atan2(sqrt(1.0 - elements_.e * elements_.e) * sinu,
                     cosu - elements_.e) + omega;

    // Вычисление скоростей
    double rdot = xke * elements_.e * sinu * sqrt(1.0 / r);
    double rfdot = xke * sqrt((1.0 - elements_.e * elements_.e) / r);

    qDebug() << "Orbital plane values:";
    qDebug() << "r (ER):" << r;
    qDebug() << "rdot (ER/min):" << rdot;
    qDebug() << "rfdot (rad/min):" << rfdot;

    // Преобразование в ECI
    double sin_u = sin(u);
    double cos_u = cos(u);
    double sin_node = sin(xnode);
    double cos_node = cos(xnode);
    double sin_i = sin(elements_.inclo);
    double cos_i = cos(elements_.inclo);

    pos.setX(r * (cos_u * cos_node - sin_u * sin_node * cos_i) * xkmper);
    pos.setY(r * (cos_u * sin_node + sin_u * cos_node * cos_i) * xkmper);
    pos.setZ(r * sin_u * sin_i * xkmper);

    vel.setX(((rdot * cos_u - r * rfdot * sin_u) * cos_node -
              (rdot * sin_u + r * rfdot * cos_u) * sin_node * cos_i) * xkmper / 60.0);
    vel.setY(((rdot * cos_u - r * rfdot * sin_u) * sin_node +
              (rdot * sin_u + r * rfdot * cos_u) * cos_node * cos_i) * xkmper / 60.0);
    vel.setZ((rdot * sin_u + r * rfdot * cos_u) * sin_i * xkmper / 60.0);

    qDebug() << "Final ECI coordinates:";
    qDebug() << "Position (km):" << pos;
    qDebug() << "Velocity (km/s):" << vel;
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;
    double tsince = epoch_.msecsTo(time) / (1000.0 * 60.0);
    propagate(tsince, state.position, state.velocity);
    return state;
}
