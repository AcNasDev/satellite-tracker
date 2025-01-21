// sgp4_propagator.cpp
#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initParameters(tle);
}

void SGP4Propagator::initParameters(const TLEParser::TLEData& tle) {
    epoch_ = tle.epoch;

    qDebug() << "\n=== SGP4 Initialization ===";

    // Базовые элементы
    elements_.i0 = tle.inclination * de2ra;
    elements_.Omega0 = tle.right_ascension * de2ra;
    elements_.e0 = tle.eccentricity;
    elements_.omega0 = tle.argument_perigee * de2ra;
    elements_.M0 = tle.mean_anomaly * de2ra;
    elements_.n0 = tle.mean_motion * 2.0 * M_PI / min_per_day;
    elements_.bstar = tle.bstar;

    qDebug() << "Initial elements:";
    qDebug() << "i0 (deg):" << tle.inclination;
    qDebug() << "Omega0 (deg):" << tle.right_ascension;
    qDebug() << "e0:" << elements_.e0;
    qDebug() << "omega0 (deg):" << tle.argument_perigee;
    qDebug() << "M0 (deg):" << tle.mean_anomaly;
    qDebug() << "n0 (rev/day):" << tle.mean_motion;
    qDebug() << "bstar:" << elements_.bstar;

    // Вычисление производных элементов
    double cosio = cos(elements_.i0);
    double theta2 = cosio * cosio;
    double x3thm1 = 3.0 * theta2 - 1.0;
    double beta02 = 1.0 - elements_.e0 * elements_.e0;
    double beta0 = sqrt(beta02);

    // Восстановление оригинальной большой полуоси
    double a1 = pow(ke / elements_.n0, 2.0/3.0);
    double delta1 = 1.5 * k2 * x3thm1 / (a1 * a1 * beta02);
    elements_.a0 = a1 * (1.0 - delta1 * (1.0/3.0 + delta1 * (1.0 + 134.0/81.0 * delta1)));

    // Вычисление производных среднего движения
    double temp = 1.5 * k2 * x3thm1 / (beta02 * elements_.a0 * elements_.a0);
    elements_.n0dot = -temp * elements_.n0 * (1.0 + 4.0 * elements_.e0);
    elements_.n0ddot = -temp * elements_.n0dot;

    qDebug() << "\nDerived elements:";
    qDebug() << "a0 (ER):" << elements_.a0;
    qDebug() << "n0dot (rad/min²):" << elements_.n0dot;
    qDebug() << "n0ddot (rad/min³):" << elements_.n0ddot;
}

double SGP4Propagator::solveKepler(double M, double e) const {
    double E = M;
    for(int i = 0; i < 10; i++) {
        double E_old = E;
        E = M + e * sin(E);
        if(fabs(E - E_old) < 1.0e-12) break;
    }
    return E;
}

void SGP4Propagator::propagate(double tsince, QVector3D& pos, QVector3D& vel) const {
    qDebug() << "\nPropagating for tsince =" << tsince << "minutes";

    // Обновление среднего движения
    double n = elements_.n0 + elements_.n0dot * tsince +
               elements_.n0ddot * tsince * tsince / 2.0;

    // Обновление эксцентриситета с учетом атмосферного торможения
    double e = elements_.e0 - elements_.bstar * tsince * n *
                                  (1.0 - elements_.e0 * elements_.e0) / elements_.e0;
    e = qMax(1.0e-6, qMin(0.999999, e));

    // Вековые возмущения
    double M = elements_.M0 + n * tsince;
    double omega = elements_.omega0 + 0.75 * k2 * n * tsince * cos(elements_.i0) /
                                          (elements_.a0 * (1.0 - e * e));
    double Omega = elements_.Omega0 - 1.5 * k2 * n * tsince * cos(elements_.i0) /
                                          (elements_.a0 * (1.0 - e * e));

    qDebug() << "Updated elements:";
    qDebug() << "n (rad/min):" << n;
    qDebug() << "e:" << e;
    qDebug() << "M (rad):" << M;
    qDebug() << "omega (rad):" << omega;
    qDebug() << "Omega (rad):" << Omega;

    // Решение уравнения Кеплера
    double E = solveKepler(M, e);

    // Вычисление истинной аномалии
    double sinE = sin(E);
    double cosE = cos(E);
    double sinv = (sqrt(1.0 - e * e) * sinE) / (1.0 - e * cosE);
    double cosv = (cosE - e) / (1.0 - e * cosE);
    double v = atan2(sinv, cosv);

    // Вычисление расстояния
    double r = elements_.a0 * (1.0 - e * cosE);
    double rdot = ke * sqrt(elements_.a0) * e * sinE / r;
    double rfdot = ke * sqrt(elements_.a0 * (1.0 - e * e)) / r;

    qDebug() << "Orbital plane values:";
    qDebug() << "r (ER):" << r;
    qDebug() << "rdot (ER/min):" << rdot;
    qDebug() << "rfdot (rad/min):" << rfdot;

    // Позиция в орбитальной плоскости
    double u = omega + v;
    double sin_u = sin(u);
    double cos_u = cos(u);
    double sin_Omega = sin(Omega);
    double cos_Omega = cos(Omega);
    double sin_i = sin(elements_.i0);
    double cos_i = cos(elements_.i0);

    // Преобразование в ECI (в км)
    pos.setX(r * (cos_u * cos_Omega - sin_u * sin_Omega * cos_i) * xkmper);
    pos.setY(r * (cos_u * sin_Omega + sin_u * cos_Omega * cos_i) * xkmper);
    pos.setZ(r * sin_u * sin_i * xkmper);

    // Скорости в ECI (в км/с)
    vel.setX(((rdot * cos_u - r * sin_u * rfdot) * cos_Omega -
              (rdot * sin_u + r * cos_u * rfdot) * sin_Omega * cos_i) * xkmper / 60.0);
    vel.setY(((rdot * cos_u - r * sin_u * rfdot) * sin_Omega +
              (rdot * sin_u + r * cos_u * rfdot) * cos_Omega * cos_i) * xkmper / 60.0);
    vel.setZ((rdot * sin_u + r * cos_u * rfdot) * sin_i * xkmper / 60.0);

    qDebug() << "\nFinal ECI coordinates:";
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
