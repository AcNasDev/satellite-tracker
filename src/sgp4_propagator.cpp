#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
}

void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    // Базовые элементы
    elements_.inclination = tle.inclination * M_PI / 180.0;
    elements_.right_ascension = tle.right_ascension * M_PI / 180.0;
    elements_.eccentricity = tle.eccentricity;
    elements_.arg_perigee = tle.argument_perigee * M_PI / 180.0;
    elements_.mean_anomaly = tle.mean_anomaly * M_PI / 180.0;
    elements_.mean_motion = tle.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;
    elements_.bstar = tle.bstar;
    elements_.epoch = tle.epoch;

    // Вычисление глубоких элементов
    deep_.cosio = cos(elements_.inclination);
    deep_.sinio = sin(elements_.inclination);
    deep_.theta2 = deep_.cosio * deep_.cosio;

    // Коррекция для эксцентриситета и движения перигея
    deep_.eosq = elements_.eccentricity * elements_.eccentricity;
    deep_.betao2 = 1.0 - deep_.eosq;
    deep_.betao = sqrt(deep_.betao2);

    // Вычисление среднего движения
    const double a1 = pow(XKE / elements_.mean_motion, 2.0/3.0);
    const double del1 = 1.5 * CK2 * (3.0 * deep_.theta2 - 1.0) /
                        (a1 * a1 * deep_.betao * deep_.betao2);

    const double ao = a1 * (1.0 - del1 * (0.5 + del1 * (1.0 + 134.0/81.0 * del1)));
    const double delo = 1.5 * CK2 * (3.0 * deep_.theta2 - 1.0) /
                        (ao * ao * deep_.betao * deep_.betao2);

    // Окончательные значения
    deep_.xnodp = elements_.mean_motion / (1.0 + delo);
    deep_.aodp = ao / (1.0 - delo);

    // Вычисление вековых возмущений
    const double xhdot1 = -1.5 * CK2 * deep_.xnodp * deep_.sinio * deep_.cosio /
                          (deep_.betao * deep_.betao2);

    deep_.xmdot = deep_.xnodp + 0.5 * xhdot1 * elements_.eccentricity *
                                    (13.0 - 78.0 * deep_.theta2 + 137.0 * deep_.theta2 * deep_.theta2);

    const double x1m5th = 1.0 - 5.0 * deep_.theta2;
    const double omgdt = -0.5 * xhdot1 * x1m5th;
    deep_.xnodot = omgdt;
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    const double tsince = elements_.epoch.secsTo(time) / 60.0;

    // Обновление средних элементов
    const double xmp = elements_.mean_anomaly + deep_.xmdot * tsince;
    const double omgasm = elements_.arg_perigee + deep_.omgdot * tsince;
    const double xnodes = elements_.right_ascension + deep_.xnodot * tsince;

    // Решение уравнения Кеплера
    double E = xmp;
    double delE;
    for(int i = 0; i < 10; i++) {
        delE = (E - elements_.eccentricity * sin(E) - xmp) /
               (1.0 - elements_.eccentricity * cos(E));
        E -= delE;
        if(fabs(delE) < 1e-12) break;
    }

    // Позиция в орбитальной плоскости
    const double cosE = cos(E);
    const double sinE = sin(E);

    const double r = deep_.aodp * XKMPER * (1.0 - elements_.eccentricity * cosE);
    double rdot = XKE * deep_.aodp * elements_.eccentricity * sinE / r;
    const double rfdot = XKE * sqrt(deep_.aodp) * sqrt(1.0 - elements_.eccentricity * elements_.eccentricity) / r;

    // Ориентация
    const double cosomg = cos(omgasm);
    const double sinomg = sin(omgasm);
    const double cosnod = cos(xnodes);
    const double sinnod = sin(xnodes);

    const double xmx = r * (cosomg * cosnod - sinomg * sinnod * deep_.cosio);
    const double xmy = r * (cosomg * sinnod + sinomg * cosnod * deep_.cosio);
    const double xmz = r * sinomg * deep_.sinio;

    // Позиция в ECI
    state.position = QVector3D(xmx, xmy, xmz);

    // Скорость в ECI
    const double rdotk = rdot * XKMPER / 60.0;  // km/s
    const double rfdotk = rfdot * XKMPER / 60.0; // km/s

    state.velocity = QVector3D(
        rdotk * (cosomg * cosnod - sinomg * sinnod * deep_.cosio) -
            rfdotk * (sinomg * cosnod + cosomg * sinnod * deep_.cosio),

        rdotk * (cosomg * sinnod + sinomg * cosnod * deep_.cosio) +
            rfdotk * (cosomg * cosnod - sinomg * sinnod * deep_.cosio),

        rdotk * sinomg * deep_.sinio + rfdotk * cosomg * deep_.sinio
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

