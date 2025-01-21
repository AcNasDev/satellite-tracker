#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
}

void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    // Сохраняем эпоху
    elements_.epoch = tle.epoch;

    // Преобразуем базовые элементы
    elements_.i = tle.inclination * DE2RA;
    elements_.Omega = tle.right_ascension * DE2RA;
    elements_.e = tle.eccentricity;
    elements_.omega = tle.argument_perigee * DE2RA;
    elements_.M = tle.mean_anomaly * DE2RA;
    elements_.bstar = tle.bstar;

    // Вычисляем среднее движение
    elements_.n0 = tle.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;

    // Вспомогательные величины
    elements_.cosio = cos(elements_.i);
    elements_.sinio = sin(elements_.i);
    elements_.eta = sqrt(1.0 - elements_.e * elements_.e);
    elements_.beta0 = sqrt(1.0 - elements_.e * elements_.e);

    // Вычисляем большую полуось
    double temp = 0.75 * J2 * (3.0 * elements_.cosio * elements_.cosio - 1.0) /
                  (elements_.beta0 * elements_.beta0 * elements_.beta0);

    double a1 = pow(XKE / elements_.n0, 2.0/3.0);
    double delta1 = temp / (a1 * a1);
    elements_.a0 = a1 * (1.0 - delta1 * (1.0/3.0 + delta1 * (1.0 + 134.0/81.0 * delta1)));

    // Корректируем среднее движение
    double delta0 = temp / (elements_.a0 * elements_.a0);
    elements_.n = elements_.n0 * (1.0 + delta0);

    // Коэффициенты для вековых возмущений
    elements_.c1 = elements_.bstar * CK2 * 1.5 / (elements_.beta0 * elements_.beta0 * elements_.beta0);

    // Большая полуось
    elements_.a = elements_.a0 * pow((elements_.n0/elements_.n), 2.0/3.0);
}

QVector3D SGP4Propagator::calculatePosVel(const double tsince) const {
    // Вековые возмущения
    double xmdf = elements_.M + elements_.n * tsince;
    double omgdf = elements_.omega + 0.75 * J2 * elements_.n * tsince *
                                         (2.0 - 3.0 * elements_.cosio * elements_.cosio) /
                                         (elements_.a * elements_.a * elements_.beta0 * elements_.beta0);
    double xnode = elements_.Omega - 1.5 * J2 * elements_.n * tsince *
                                         elements_.cosio / (elements_.a * elements_.a * elements_.beta0 * elements_.beta0);

    // Решаем уравнение Кеплера
    double E = solveKeplerEquation(xmdf, elements_.e);

    // Вычисляем истинную аномалию
    double sinE = sin(E);
    double cosE = cos(E);
    double nu = atan2(elements_.eta * sinE, cosE - elements_.e);

    // Вычисляем радиус-вектор
    double r = elements_.a * (1.0 - elements_.e * cosE);

    // Позиция в орбитальной плоскости
    double sin_nu = sin(nu);
    double cos_nu = cos(nu);
    double x = r * cos_nu;
    double y = r * sin_nu;

    // Преобразуем в ECI
    double sin_i = elements_.sinio;
    double cos_i = elements_.cosio;
    double sin_node = sin(xnode);
    double cos_node = cos(xnode);
    double sin_arg = sin(omgdf);
    double cos_arg = cos(omgdf);

    // Позиция в ECI (км)
    double xh = x * (cos_node * cos_arg - sin_node * sin_arg * cos_i);
    double yh = x * (sin_node * cos_arg + cos_node * sin_arg * cos_i);
    double zh = x * sin_arg * sin_i;

    xh = xh * XKMPER;
    yh = yh * XKMPER;
    zh = zh * XKMPER;

    return QVector3D(xh, yh, zh);
}

double SGP4Propagator::solveKeplerEquation(double M, double e) const {
    double E = M;
    for(int i = 0; i < 10; i++) {
        double delta = (M + e * sin(E) - E) / (1.0 - e * cos(E));
        E += delta;
        if(fabs(delta) < 1e-12) break;
    }
    return E;
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = elements_.epoch.msecsTo(time) / (1000.0 * 60.0);

    // Вычисляем положение и скорость
    QVector3D pos_vel = calculatePosVel(tsince);

    state.position = pos_vel;

    // Вычисляем скорость численным методом
    const double dt = 0.0001; // малый шаг времени
    QVector3D pos2 = calculatePosVel(tsince + dt);
    state.velocity = (pos2 - pos_vel) / dt;

    return state;
}
