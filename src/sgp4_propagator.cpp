#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
}

void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    elements_.inclination = tle.inclination * DE2RA;
    elements_.right_ascension = tle.right_ascension * DE2RA;
    elements_.eccentricity = tle.eccentricity;
    elements_.arg_perigee = tle.argument_perigee * DE2RA;
    elements_.mean_anomaly = tle.mean_anomaly * DE2RA;
    elements_.mean_motion = tle.mean_motion;
    elements_.bstar = tle.bstar;
    elements_.epoch = tle.epoch;

    // Преобразование среднего движения
    elements_.n0 = elements_.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;

    // Вычисление большой полуоси
    double a1 = std::pow(KE / elements_.n0, 2.0/3.0);
    elements_.a = a1 * AE; // Масштабирование к километрам

    // Коррекция для J2
    double delta1 = 1.5 * CK2 * (3.0 * std::pow(std::cos(elements_.inclination), 2.0) - 1.0) /
                    std::pow(1.0 - elements_.eccentricity * elements_.eccentricity, 1.5);
    double a0 = a1 * (1.0 - delta1/3.0 - delta1 * delta1 - 134.0 * delta1 * delta1 * delta1 / 81.0);
    elements_.a = a0 * AE;

    // Коррекция среднего движения
    elements_.n = elements_.n0 / (1.0 + delta1);
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = elements_.epoch.secsTo(time) / 60.0;

    // Учет вековых возмущений
    double xndt = elements_.n * tsince;
    double xldot = elements_.n + xndt;

    // Решение уравнения Кеплера
    double xn = elements_.n;
    double xl = elements_.mean_anomaly + xldot * tsince;
    double e = elements_.eccentricity;

    double u = xl - elements_.right_ascension;
    double sin_u = std::sin(u);
    double cos_u = std::cos(u);

    // Расчет позиции в орбитальной плоскости (с учетом масштаба)
    double r = elements_.a * (1.0 - e * std::cos(xl));
    double rdot = elements_.a * elements_.n * e * std::sin(xl) / std::sqrt(1.0 - e * e);

    // Преобразование в ECI с правильным масштабом
    double xmx = -sin_u * std::cos(elements_.inclination);
    double xmy = cos_u * std::cos(elements_.inclination);
    double xmz = std::sin(elements_.inclination) * std::sin(u);

    state.position = QVector3D(
        r * cos_u,
        r * sin_u * std::cos(elements_.inclination),
        r * sin_u * std::sin(elements_.inclination)
        );

    double vm = std::sqrt(MU / (r * AE)) * AE;
    state.velocity = QVector3D(
        vm * (-sin_u),
        vm * (cos_u * std::cos(elements_.inclination)),
        vm * (cos_u * std::sin(elements_.inclination))
        );

    return state;
}

QVector3D SGP4Propagator::getPosition(double r, double u, double i, double omega) const {
    double cosu = qCos(u);
    double sinu = qSin(u);
    double cosi = qCos(i);
    double sini = qSin(i);
    double cosomega = qCos(omega);
    double sinomega = qSin(omega);

    return QVector3D(
        r * (cosu * cosomega - sinu * sinomega * cosi),
        r * (cosu * sinomega + sinu * cosomega * cosi),
        r * sinu * sini
        );
}

QVector3D SGP4Propagator::getVelocity(double r, double rdot, double u,
                                      double rfdot, double i, double omega) const {
    double cosu = qCos(u);
    double sinu = qSin(u);
    double cosi = qCos(i);
    double sini = qSin(i);
    double cosomega = qCos(omega);
    double sinomega = qSin(omega);

    return QVector3D(
        (rdot * cosu - r * sinu * rfdot) * cosomega -
            (rdot * sinu + r * cosu * rfdot) * sinomega * cosi,

        (rdot * cosu - r * sinu * rfdot) * sinomega +
            (rdot * sinu + r * cosu * rfdot) * cosomega * cosi,

        (rdot * sinu + r * cosu * rfdot) * sini
        );
}

