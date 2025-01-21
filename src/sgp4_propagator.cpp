#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
}

void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    elements_.epoch = tle.epoch;

    // Преобразуем основные элементы
    elements_.i = tle.inclination * DE2RA;
    elements_.Omega = tle.right_ascension * DE2RA;
    elements_.e = tle.eccentricity;
    elements_.omega = tle.argument_perigee * DE2RA;
    elements_.M = tle.mean_anomaly * DE2RA;
    elements_.n = tle.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;
    elements_.bstar = tle.bstar;

    // Вычисляем большую полуось
    const double mu = 398600.4418;  // Гравитационный параметр Земли (км³/с²)
    elements_.a = pow(mu / (elements_.n * elements_.n), 1.0/3.0);
}

double SGP4Propagator::solveKepler(double meanAnomaly, double eccentricity) const {
    double E = meanAnomaly;
    double delta;

    // Метод Ньютона с фиксированным числом итераций
    for(int i = 0; i < 10; i++) {
        delta = E - eccentricity * sin(E) - meanAnomaly;
        E = E - delta / (1.0 - eccentricity * cos(E));

        if(fabs(delta) < 1e-12) {
            break;
        }
    }

    return E;
}

QVector3D SGP4Propagator::calculatePositionAndVelocity(double tsince) const {
    // Обновляем среднюю аномалию с учетом времени
    double M = elements_.M + elements_.n * tsince;

    // Решаем уравнение Кеплера
    double E = solveKepler(M, elements_.e);

    // Вычисляем истинную аномалию
    double sinE = sin(E);
    double cosE = cos(E);
    double nu = atan2(sqrt(1.0 - elements_.e * elements_.e) * sinE,
                      cosE - elements_.e);

    // Вычисляем расстояние до спутника
    double r = elements_.a * (1.0 - elements_.e * cosE);

    // Вычисляем координаты в орбитальной плоскости
    double xw = r * (cos(nu + elements_.omega));
    double yw = r * (sin(nu + elements_.omega));

    // Преобразуем в ECI
    double cosO = cos(elements_.Omega);
    double sinO = sin(elements_.Omega);
    double cosI = cos(elements_.i);
    double sinI = sin(elements_.i);

    double x = (xw * cosO - yw * cosI * sinO) * XKMPER;
    double y = (xw * sinO + yw * cosI * cosO) * XKMPER;
    double z = (yw * sinI) * XKMPER;

    return QVector3D(x, y, z);
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = elements_.epoch.msecsTo(time) / (1000.0 * 60.0);

    // Вычисляем позицию в ECI
    QVector3D pos_vel = calculatePositionAndVelocity(tsince);
    state.position = pos_vel;

    // Вычисляем скорость численным методом
    double dt = 1.0; // 1 секунда
    QVector3D pos2 = calculatePositionAndVelocity(tsince + dt/60.0);
    state.velocity = (pos2 - pos_vel) / dt;

    return state;
}
