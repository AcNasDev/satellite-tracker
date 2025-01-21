// sgp4_propagator.cpp
#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initParameters(tle);
}

void SGP4Propagator::initParameters(const TLEParser::TLEData& tle) {
    epoch_ = tle.epoch;

    qDebug() << "\n=== SGP4 Initialization ===";

    // Конвертация базовых элементов
    elements_.io = tle.inclination * de2ra;
    elements_.eo = tle.eccentricity;
    elements_.omegao = tle.argument_perigee * de2ra;
    elements_.xnodeo = tle.right_ascension * de2ra;
    elements_.xmo = tle.mean_anomaly * de2ra;
    elements_.no = tle.mean_motion * 2.0 * M_PI / min_per_day;
    elements_.bstar = tle.bstar;

    qDebug() << "Raw TLE elements:";
    qDebug() << "io (rad):" << elements_.io;
    qDebug() << "eo:" << elements_.eo;
    qDebug() << "omegao (rad):" << elements_.omegao;
    qDebug() << "xnodeo (rad):" << elements_.xnodeo;
    qDebug() << "xmo (rad):" << elements_.xmo;
    qDebug() << "no (rad/min):" << elements_.no;
    qDebug() << "bstar:" << elements_.bstar;

    // Вспомогательные величины
    double cosio = cos(elements_.io);
    double theta2 = cosio * cosio;
    double x3thm1 = 3.0 * theta2 - 1.0;
    double eosq = elements_.eo * elements_.eo;
    double betao2 = 1.0 - eosq;
    double betao = sqrt(betao2);

    // Восстановление оригинальной средней движения (no) и полуоси (a)
    double a1 = pow(ke / elements_.no, tothrd);
    double del1 = 1.5 * ck2 * x3thm1 / (a1 * a1 * betao * betao2);
    double ao = a1 * (1.0 - del1 * (0.5 * tothrd + del1 * (1.0 + 134.0/81.0 * del1)));
    double delo = 1.5 * ck2 * x3thm1 / (ao * ao * betao * betao2);
    double xnodp = elements_.no / (1.0 + delo);
    elements_.a = ao;

    // Инициализация для SGP4
    elements_.inclo = elements_.io;
    elements_.nodeo = elements_.xnodeo;
    elements_.argpo = elements_.omegao;
    elements_.mo = elements_.xmo;

    // Вычисление производных среднего движения
    double temp = 1.5 * ck2 * x3thm1 / (betao * betao2);
    elements_.ndot = -temp * xnodp * (1.0 + 4.0 * betao);
    elements_.nddot = -temp * elements_.ndot;

    // Сохраняем среднее движение Козаи
    elements_.no_kozai = xnodp;

    qDebug() << "\nDerived elements:";
    qDebug() << "Semi-major axis (ER):" << elements_.a;
    qDebug() << "Mean motion (Kozai, rad/min):" << elements_.no_kozai;
    qDebug() << "ndot (rad/min²):" << elements_.ndot;
    qDebug() << "nddot (rad/min³):" << elements_.nddot;
}

void SGP4Propagator::propagate(double tsince, QVector3D& pos, QVector3D& vel) const {
    qDebug() << "\nPropagating for tsince =" << tsince << "minutes";

    // Вековые возмущения
    double xmp = elements_.mo + elements_.no_kozai * tsince;
    double omega = elements_.argpo + 0.75 * ck2 * elements_.no_kozai * tsince *
                                         cos(elements_.inclo) / (elements_.a * elements_.a * (1.0 - elements_.eo * elements_.eo));
    double xnode = elements_.nodeo - 1.5 * ck2 * elements_.no_kozai * tsince *
                                         cos(elements_.inclo) / (elements_.a * elements_.a * (1.0 - elements_.eo * elements_.eo));

    // Обновление эксцентриситета
    double e = elements_.eo - elements_.bstar * elements_.no_kozai * tsince *
                                  (1.0 - elements_.eo * elements_.eo);
    e = qMax(1.0e-6, qMin(0.999999, e));

    qDebug() << "Secular perturbations:";
    qDebug() << "Mean anomaly (rad):" << xmp;
    qDebug() << "Argument of perigee (rad):" << omega;
    qDebug() << "RAAN (rad):" << xnode;
    qDebug() << "Eccentricity:" << e;

    // Решение уравнения Кеплера
    double xl = xmp + omega;
    double E = xl;
    for(int i = 0; i < 10; i++) {
        double E_old = E;
        E = xl + e * sin(E);
        if(fabs(E - E_old) < 1.0e-12) break;
    }

    // Вычисление позиции в орбитальной плоскости
    double sinE = sin(E);
    double cosE = cos(E);
    double sinv = (sqrt(1.0 - e * e) * sinE) / (1.0 - e * cosE);
    double cosv = (cosE - e) / (1.0 - e * cosE);
    double v = atan2(sinv, cosv);
    double u = v + omega;

    double r = elements_.a * (1.0 - e * cosE);
    double rdot = ke * sqrt(elements_.a) * e * sinE / r;
    double rfdot = ke * sqrt(elements_.a * (1.0 - e * e)) / r;

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
