#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
}

void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    elements_.inclination = tle.inclination * M_PI / 180.0;
    elements_.right_ascension = tle.right_ascension * M_PI / 180.0;
    elements_.eccentricity = tle.eccentricity;
    elements_.arg_perigee = tle.argument_perigee * M_PI / 180.0;
    elements_.mean_anomaly = tle.mean_anomaly * M_PI / 180.0;
    elements_.mean_motion = tle.mean_motion * TWOPI / XMNPDA;
    elements_.bstar = tle.bstar;
    elements_.epoch = tle.epoch;

    // Инициализация глубоких элементов
    deep_.cosio = cos(elements_.inclination);
    deep_.sinio = sin(elements_.inclination);
    const double theta2 = deep_.cosio * deep_.cosio;
    const double x3thm1 = 3.0 * theta2 - 1.0;
    const double eosq = elements_.eccentricity * elements_.eccentricity;
    const double betao2 = 1.0 - eosq;
    const double betao = sqrt(betao2);

    // Коррекция для SGP4
    const double a1 = pow(XKE / elements_.mean_motion, 2.0/3.0);
    const double del1 = 1.5 * CK2 * x3thm1 / (a1 * a1 * betao * betao2);
    const double ao = a1 * (1.0 - del1 * (0.5 + del1 * (1.0 + 134.0/81.0 * del1)));
    const double delo = 1.5 * CK2 * x3thm1 / (ao * ao * betao * betao2);

    deep_.aodp = ao;
    deep_.xnodp = elements_.mean_motion / (1.0 + delo);

    // Вековые эффекты
    const double xhdot1 = -1.5 * CK2 * deep_.xnodp * deep_.sinio * deep_.cosio / (betao * betao2);
    deep_.omgdot = -0.5 * xhdot1 * (5.0 * theta2 - 1.0);
    deep_.xnodot = xhdot1 * (1.0 - 2.0 * deep_.cosio);
    deep_.xmdot = deep_.xnodp + 0.5 * xhdot1 * elements_.eccentricity * deep_.sinio;
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    const double tsince = elements_.epoch.secsTo(time) / 60.0;

    // Обновление средних элементов
    const double xmdf = elements_.mean_anomaly + deep_.xmdot * tsince;
    const double omgadf = elements_.arg_perigee + deep_.omgdot * tsince;
    const double xnode = elements_.right_ascension + deep_.xnodot * tsince;

    // Решение уравнения Кеплера
    double E = xmdf;
    for(int i = 0; i < 10; i++) {
        const double delE = (E - elements_.eccentricity * sin(E) - xmdf) /
                            (1.0 - elements_.eccentricity * cos(E));
        E -= delE;
        if(fabs(delE) < 1e-12) break;
    }

    // Позиция в орбитальной плоскости
    const double sinE = sin(E);
    const double cosE = cos(E);

    const double r = deep_.aodp * (1.0 - elements_.eccentricity * cosE);
    const double rdot = deep_.aodp * elements_.eccentricity * sinE * deep_.xnodp / sqrt(1.0 - elements_.eccentricity * elements_.eccentricity);
    const double rfdot = deep_.xnodp * deep_.aodp * sqrt(1.0 - elements_.eccentricity * elements_.eccentricity) / r;

    // Ориентация
    const double cosw = cos(omgadf);
    const double sinw = sin(omgadf);
    const double cosnod = cos(xnode);
    const double sinnod = sin(xnode);

    // Позиция в ECI (km)
    const double sinv = (sqrt(1.0 - elements_.eccentricity * elements_.eccentricity) * sinE) / (1.0 - elements_.eccentricity * cosE);
    const double cosv = (cosE - elements_.eccentricity) / (1.0 - elements_.eccentricity * cosE);

    const double rx = r * (cosw * cosv - sinw * sinv);
    const double ry = r * (sinw * cosv + cosw * sinv);
    const double rz = r * sinv * deep_.sinio;

    state.position = QVector3D(
        XKMPER * (rx * cosnod - ry * sinnod * deep_.cosio),
        XKMPER * (rx * sinnod + ry * cosnod * deep_.cosio),
        XKMPER * ry * deep_.sinio
        );

    // Скорость в ECI (km/s)
    const double rdotk = rdot * XKMPER / 60.0;
    const double rfdotk = rfdot * XKMPER / 60.0;

    state.velocity = QVector3D(
        rdotk * (cosw * cosv - sinw * sinv) * cosnod -
            rfdotk * (sinw * cosv + cosw * sinv) * cosnod -
            (rdotk * (cosw * cosv - sinw * sinv) * sinnod -
             rfdotk * (sinw * cosv + cosw * sinv) * sinnod) * deep_.cosio,

        rdotk * (cosw * cosv - sinw * sinv) * sinnod +
            rfdotk * (sinw * cosv + cosw * sinv) * sinnod +
            (rdotk * (cosw * cosv - sinw * sinv) * cosnod -
             rfdotk * (sinw * cosv + cosw * sinv) * cosnod) * deep_.cosio,

        rdotk * sinv * deep_.sinio + rfdotk * cosv * deep_.sinio
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

