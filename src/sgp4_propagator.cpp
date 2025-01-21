#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initParameters(tle);
}

void SGP4Propagator::initParameters(const TLEParser::TLEData& tle) {
    epoch_ = tle.epoch;

    qDebug() << "\n=== SGP4 Initialization ===";

    // Конвертация основных элементов
    elements_.inclo = tle.inclination * de2ra;
    elements_.nodeo = tle.right_ascension * de2ra;
    elements_.ecco = tle.eccentricity;
    elements_.argpo = tle.argument_perigee * de2ra;
    elements_.mo = tle.mean_anomaly * de2ra;
    elements_.no = tle.mean_motion * 2.0 * M_PI / min_per_day;
    elements_.bstar = tle.bstar / ae;  // Важное исправление: масштабирование BSTAR

    qDebug() << "Initial elements:";
    qDebug() << "Inclination (rad):" << elements_.inclo;
    qDebug() << "RAAN (rad):" << elements_.nodeo;
    qDebug() << "Eccentricity:" << elements_.ecco;
    qDebug() << "Arg of perigee (rad):" << elements_.argpo;
    qDebug() << "Mean anomaly (rad):" << elements_.mo;
    qDebug() << "Mean motion (rad/min):" << elements_.no;
    qDebug() << "BSTAR (1/ER):" << elements_.bstar;

    // Вспомогательные величины
    double cosio = cos(elements_.inclo);
    double theta2 = cosio * cosio;
    double x3thm1 = 3.0 * theta2 - 1.0;
    double eosq = elements_.ecco * elements_.ecco;
    double betao2 = 1.0 - eosq;
    double betao = sqrt(betao2);

    // Восстановление оригинального среднего движения
    double ao = pow(xke/elements_.no, 2.0/3.0);
    double delo = 1.5 * ck2 * x3thm1 / (ao * ao * betao * betao2);
    elements_.a = ao * (1.0 - delo/3.0 - delo * delo - 134.0 * delo * delo * delo/81.0);

    // Сохранение среднего движения Козаи
    elements_.no_kozai = elements_.no / (1.0 + delo);

    // Вычисление производных среднего движения
    double temp = 1.5 * ck2 * x3thm1 / (betao * betao2);
    elements_.ndot = -temp * elements_.no_kozai * (1.0 + 4.0 * betao + eosq);
    elements_.nddot = -temp * elements_.ndot;

    // Вычисление апогея и перигея
    elements_.alta = elements_.a * (1.0 + elements_.ecco) - 1.0;
    elements_.altp = elements_.a * (1.0 - elements_.ecco) - 1.0;

    qDebug() << "\nDerived elements:";
    qDebug() << "Semi-major axis (ER):" << elements_.a;
    qDebug() << "Mean motion (Kozai, rad/min):" << elements_.no_kozai;
    qDebug() << "ndot (rad/min²):" << elements_.ndot;
    qDebug() << "nddot (rad/min³):" << elements_.nddot;
    qDebug() << "Apogee height (ER):" << elements_.alta;
    qDebug() << "Perigee height (ER):" << elements_.altp;
}

void SGP4Propagator::propagate(double tsince, QVector3D& pos, QVector3D& vel) {
    qDebug() << "\nPropagating for tsince =" << tsince << "minutes";

    // Вековые возмущения
    double xn = elements_.no_kozai;
    double em = elements_.ecco;
    double xinc = elements_.inclo;

    // Учет вековых возмущений
    double xnoddf = elements_.ndot * tsince * 2.0;
    double omega = elements_.argpo;
    double xnode = elements_.nodeo + xnoddf;
    double e = em - elements_.bstar * tsince;
    e = qMax(1.0e-6, qMin(0.999999, e));
    double xmp = elements_.mo + xn * tsince + xnoddf/2.0;

    qDebug() << "Secular perturbations:";
    qDebug() << "Mean motion (rad/min):" << xn;
    qDebug() << "Eccentricity:" << e;
    qDebug() << "Mean anomaly (rad):" << xmp;
    qDebug() << "Arg of perigee (rad):" << omega;
    qDebug() << "RAAN (rad):" << xnode;

    // Решение уравнения Кеплера
    double U = solveKepler(xmp, e);

    // Вычисление позиции в орбитальной плоскости
    double sin_U = sin(U);
    double cos_U = cos(U);
    double r = elements_.a * (1.0 - e * cos_U);
    double rdot = xke * sqrt(elements_.a) * e * sin_U / r;
    double rfdot = xke * sqrt(elements_.a * (1.0 - e * e)) / r;

    qDebug() << "Orbital plane values:";
    qDebug() << "r (ER):" << r;
    qDebug() << "rdot (ER/min):" << rdot;
    qDebug() << "rfdot (rad/min):" << rfdot;

    // Вычисление координат в ECI
    double sin_i = sin(xinc);
    double cos_i = cos(xinc);
    double sin_node = sin(xnode);
    double cos_node = cos(xnode);
    double sin_u = sin(U + omega);
    double cos_u = cos(U + omega);

    pos.setX(r * (cos_u * cos_node - sin_u * sin_node * cos_i) * xkmper);
    pos.setY(r * (cos_u * sin_node + sin_u * cos_node * cos_i) * xkmper);
    pos.setZ(r * sin_u * sin_i * xkmper);

    vel.setX(((rdot * cos_u - r * rfdot * sin_u) * cos_node -
              (rdot * sin_u + r * rfdot * cos_u) * sin_node * cos_i) * xkmper / 60.0);
    vel.setY(((rdot * cos_u - r * rfdot * sin_u) * sin_node +
              (rdot * sin_u + r * rfdot * cos_u) * cos_node * cos_i) * xkmper / 60.0);
    vel.setZ((rdot * sin_u + r * rfdot * cos_u) * sin_i * xkmper / 60.0);

    qDebug() << "\nFinal ECI coordinates:";
    qDebug() << "Position (km):" << pos;
    qDebug() << "Velocity (km/s):" << vel;
}

double SGP4Propagator::solveKepler(double M, double e) {
    double E = M;
    for(int i = 0; i < 10; i++) {
        double E_old = E;
        E = M + e * sin(E);
        if(fabs(E - E_old) < 1.0e-12) break;
    }
    return E;
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) {
    OrbitalState state;
    state.epoch = time;
    double tsince = epoch_.msecsTo(time) / (1000.0 * 60.0);
    propagate(tsince, state.position, state.velocity);
    return state;
}

