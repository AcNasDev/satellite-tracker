#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
}

void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    // Преобразование базовых элементов
    elements_.epoch = tle.epoch;
    sgp4_.i = tle.inclination * DE2RA;
    sgp4_.Omega = tle.right_ascension * DE2RA;
    sgp4_.e = tle.eccentricity;
    sgp4_.omega = tle.argument_perigee * DE2RA;
    sgp4_.M = tle.mean_anomaly * DE2RA;
    sgp4_.n0 = tle.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;
    sgp4_.bstar = tle.bstar;

    // Вычисление производных элементов
    const double cos_i = cos(sgp4_.i);
    const double theta2 = cos_i * cos_i;
    const double x3thm1 = 3.0 * theta2 - 1.0;
    const double e2 = sgp4_.e * sgp4_.e;
    const double beta02 = 1.0 - e2;
    const double beta0 = sqrt(beta02);

    // Коррекция для SGP4
    const double a1 = pow(XKE / sgp4_.n0, 2.0/3.0);
    const double del1 = 1.5 * CK2 * x3thm1 / (a1 * a1 * beta02);
    sgp4_.a0 = a1 * (1.0 - del1 * (1.0/3.0 + del1 * (1.0 + 134.0/81.0 * del1)));
    const double del0 = 1.5 * CK2 * x3thm1 / (sgp4_.a0 * sgp4_.a0 * beta02);

    // Обновление mean motion
    sgp4_.n = sgp4_.n0 / (1.0 + del0);
    sgp4_.a = sgp4_.a0;

    // Вычисление secular effects
    const double aux = 1.0 - theta2;
    const double sinio = sin(sgp4_.i);
    const double cosio = cos(sgp4_.i);
    const double x7thm1 = 7.0 * theta2 - 1.0;

    const double c1 = CK2 * 1.5;
    const double c4 = 2.0 * (sgp4_.a0 * sgp4_.a0 * beta02 * beta02);
    const double c5 = 2.0 * c1 * sgp4_.a0 * beta02;

    sgp4_.ndot = -c1 * sgp4_.n * x3thm1 / (sgp4_.a0 * beta02);
    sgp4_.nddot = -sgp4_.n * c1 * (sgp4_.n / (sgp4_.a0 * beta02)) * (-7.0/2.0 + 19.0/4.0 * theta2 + 2.0 * theta2 * theta2);
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Правильный расчет времени с эпохи
    const double minutes_per_day = 1440.0;
    const double days_since_epoch = elements_.epoch.daysTo(time) +
                                    elements_.epoch.secsTo(time) % 86400 / 86400.0;
    const double tsince = days_since_epoch * minutes_per_day;

    // Обновление средних элементов
    const double M = sgp4_.M + sgp4_.n * tsince;
    const double omega = sgp4_.omega + 0.75 * CK2 * sgp4_.n * tsince * cos(sgp4_.i) / sqrt(1.0 - sgp4_.e * sgp4_.e);
    const double Omega = sgp4_.Omega - 1.5 * CK2 * sgp4_.n * tsince * cos(sgp4_.i) / sqrt(1.0 - sgp4_.e * sgp4_.e);

    // Решение уравнения Кеплера
    const double E = Kepler(M, sgp4_.e);
    const double sin_E = sin(E);
    const double cos_E = cos(E);

    // True anomaly
    const double sin_v = (sqrt(1.0 - sgp4_.e * sgp4_.e) * sin_E) / (1.0 - sgp4_.e * cos_E);
    const double cos_v = (cos_E - sgp4_.e) / (1.0 - sgp4_.e * cos_E);

    // Position in orbital plane
    const double r = sgp4_.a * (1.0 - sgp4_.e * cos_E);
    const double v = atan2(sin_v, cos_v);

    // Orientation vectors
    const double sin_omega = sin(omega);
    const double cos_omega = cos(omega);
    const double sin_Omega = sin(Omega);
    const double cos_Omega = cos(Omega);
    const double sin_i = sin(sgp4_.i);
    const double cos_i = cos(sgp4_.i);

    // Position vector (in km)
    const double x = r * (cos_omega * cos_v - sin_omega * sin_v);
    const double y = r * (sin_omega * cos_v + cos_omega * sin_v);

    state.position = QVector3D(
        XKMPER * (x * cos_Omega - y * cos_i * sin_Omega),
        XKMPER * (x * sin_Omega + y * cos_i * cos_Omega),
        XKMPER * y * sin_i
        );

    // Velocity calculation
    const double mu = XKE * XKE;
    const double n = sqrt(mu / (r * r * r));
    const double rdot = sqrt(mu / sgp4_.a) * sgp4_.e * sin_E / r;
    const double rfdot = sqrt(mu * sgp4_.a) * sqrt(1 - sgp4_.e * sgp4_.e) / r;

    state.velocity = QVector3D(
        (XKMPER/60.0) * (-sin_Omega * (rdot * cos_v - rfdot * sin_v) -
                           cos_Omega * (rdot * sin_v + rfdot * cos_v) * cos_i),
        (XKMPER/60.0) * (cos_Omega * (rdot * cos_v - rfdot * sin_v) -
                           sin_Omega * (rdot * sin_v + rfdot * cos_v) * cos_i),
        (XKMPER/60.0) * (rdot * sin_v + rfdot * cos_v) * sin_i
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

