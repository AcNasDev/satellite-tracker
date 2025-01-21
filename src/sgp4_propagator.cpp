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
    deep_.theta2 = deep_.cosio * deep_.cosio;
    deep_.x3thm1 = 3.0 * deep_.theta2 - 1.0;
    const double eosq = elements_.eccentricity * elements_.eccentricity;
    const double betao2 = 1.0 - eosq;
    const double betao = sqrt(betao2);

    // Инициализация SGP4
    const double a1 = pow(XKE / elements_.mean_motion, TOTHRD);
    const double del1 = 1.5 * CK2 * deep_.x3thm1 / (a1 * a1 * betao * betao2);
    const double ao = a1 * (1.0 - del1 * (0.5 + del1 * (1.0 + 134.0/81.0 * del1)));
    const double delo = 1.5 * CK2 * deep_.x3thm1 / (ao * ao * betao * betao2);

    deep_.aodp = ao;
    deep_.xnodp = elements_.mean_motion / (1.0 + delo);

    // Расчет секулярных эффектов
    deep_.omgdot_factor = -1.5 * CK2 * deep_.xnodp * deep_.sinio * deep_.cosio / betao2;
    deep_.xnodot_factor = deep_.omgdot_factor * deep_.cosio;
    deep_.xmdot_factor = deep_.xnodp + 0.5 * deep_.omgdot_factor * elements_.eccentricity * deep_.sinio;

    // Инициализация коэффициентов для долгопериодических возмущений
    deep_.xlcof = 0.125 * CK2 * deep_.sinio * (3.0 + 5.0 * deep_.cosio) / (1.0 + deep_.cosio);
    deep_.aycof = 0.25 * CK2 * deep_.sinio;
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    const double tsince = elements_.epoch.secsTo(time) / 60.0;

    // Обновление средних элементов
    const double xmp = elements_.mean_anomaly + deep_.xmdot_factor * tsince;
    const double omgadf = elements_.arg_perigee + deep_.omgdot_factor * tsince;
    const double xnode = elements_.right_ascension + deep_.xnodot_factor * tsince;

    // Решение уравнения Кеплера
    double E = xmp;
    double delE;
    for(int i = 0; i < 10; i++) {
        delE = (E - elements_.eccentricity * sin(E) - xmp) /
               (1.0 - elements_.eccentricity * cos(E));
        E -= delE;
        if(fabs(delE) < 1e-12) break;
    }

    // Вычисление позиции в орбитальной плоскости
    const double sinE = sin(E);
    const double cosE = cos(E);

    const double sinv = (sqrt(1.0 - elements_.eccentricity * elements_.eccentricity) * sinE) /
                        (1.0 - elements_.eccentricity * cosE);
    const double cosv = (cosE - elements_.eccentricity) / (1.0 - elements_.eccentricity * cosE);

    const double r = deep_.aodp * (1.0 - elements_.eccentricity * cosE);
    const double rdot = deep_.xnodp * deep_.aodp * elements_.eccentricity * sinE / sqrt(1.0 - elements_.eccentricity * elements_.eccentricity);
    const double rfdot = deep_.xnodp * deep_.aodp * sqrt(1.0 - elements_.eccentricity * elements_.eccentricity) / r;

    // Ориентация
    const double cosomg = cos(omgadf);
    const double sinomg = sin(omgadf);
    const double cosnod = cos(xnode);
    const double sinnod = sin(xnode);

    // Позиция и скорость в ECI
    const double xmx = r * (cosomg * cosv - sinomg * sinv);
    const double xmy = r * (sinomg * cosv + cosomg * sinv);

    state.position = QVector3D(
        XKMPER * (xmx * cosnod - xmy * sinnod * deep_.cosio),
        XKMPER * (xmx * sinnod + xmy * cosnod * deep_.cosio),
        XKMPER * xmy * deep_.sinio
        );

    // Преобразование скорости в км/с
    const double vx = rdot * (cosomg * cosv - sinomg * sinv) -
                      rfdot * (cosomg * sinv + sinomg * cosv);
    const double vy = rdot * (sinomg * cosv + cosomg * sinv) -
                      rfdot * (sinomg * sinv - cosomg * cosv);

    state.velocity = QVector3D(
        (XKMPER/60.0) * (vx * cosnod - vy * sinnod * deep_.cosio),
        (XKMPER/60.0) * (vx * sinnod + vy * cosnod * deep_.cosio),
        (XKMPER/60.0) * vy * deep_.sinio
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

