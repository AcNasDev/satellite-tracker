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
    elements_.mean_motion = tle.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;
    elements_.bstar = tle.bstar;
    elements_.epoch = tle.epoch;

    // Вычисление производных параметров
    double a1 = std::pow(XKE / elements_.mean_motion, 2.0/3.0);
    double cosio = qCos(elements_.inclination);
    double theta2 = cosio * cosio;
    double x3thm1 = 3.0 * theta2 - 1.0;
    double eosq = elements_.eccentricity * elements_.eccentricity;
    double betao2 = 1.0 - eosq;
    double betao = qSqrt(betao2);

    // Коррекция среднего движения
    double del1 = 1.5 * CK2 * x3thm1 / (betao * betao2);
    elements_.a = a1 * (1.0 - del1 * (1.0/3.0 + del1 * (1.0 + 134.0/81.0 * del1)));
    double ao = elements_.a;
    elements_.n0 = elements_.mean_motion;
    elements_.n = elements_.n0 + del1 * elements_.n0 * (1.0 + del1 * 134.0/81.0);
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = elements_.epoch.secsTo(time) / 60.0;

    // Учет вековых возмущений
    double xmdf = elements_.mean_anomaly + elements_.n * tsince;
    double omgadf = elements_.arg_perigee;
    double xnode = elements_.right_ascension;
    double bstar = elements_.bstar;

    // Расчет возмущений от J2
    double cosio = qCos(elements_.inclination);
    double sinio = qSin(elements_.inclination);
    double eta = qSqrt(1.0 - elements_.eccentricity * elements_.eccentricity);
    double c1 = CK2 * 1.5;
    double a0 = qPow(XKE / elements_.n0, 2.0/3.0);
    double d0 = 3.0 * c1 * qPow(sinio, 2.0) / (2.0 * eta * eta * a0 * a0);
    double d1 = d0 * elements_.n0;

    // Обновление средней аномалии
    double mp = xmdf + d1 * tsince;

    // Решение уравнения Кеплера
    double e = elements_.eccentricity;
    double E = mp;
    for(int i = 0; i < 10; i++) {
        double delta = E - e * qSin(E) - mp;
        E = E - delta / (1.0 - e * qCos(E));
        if(qAbs(delta) < 1e-12) break;
    }

    // Расчет позиции в орбитальной плоскости
    double cosE = qCos(E);
    double sinE = qSin(E);

    double r = elements_.a * (1.0 - e * cosE);
    double rdot = elements_.a * elements_.n * e * sinE / r;

    double cosnu = (cosE - e) / (1.0 - e * cosE);
    double sinnu = qSqrt(1.0 - e * e) * sinE / (1.0 - e * cosE);
    double nu = qAtan2(sinnu, cosnu);

    // Ориентация орбитальной плоскости
    double u = nu + omgadf;
    double rfdot = elements_.a * elements_.n * qSqrt(1.0 - e * e) / r;

    // Преобразование в ECI координаты
    state.position = getPosition(r, u, elements_.inclination, xnode);
    state.velocity = getVelocity(r, rdot, u, rfdot, elements_.inclination, xnode);

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

