#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
}

// sgp4_propagator.cpp
void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    elements_.inclination = tle.inclination * M_PI / 180.0;
    elements_.right_ascension = tle.right_ascension * M_PI / 180.0;
    elements_.eccentricity = tle.eccentricity;
    elements_.arg_perigee = tle.argument_perigee * M_PI / 180.0;
    elements_.mean_anomaly = tle.mean_anomaly * M_PI / 180.0;

    // Convert mean motion to radians per minute
    const double xno = tle.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;
    elements_.mean_motion = xno;
    elements_.bstar = tle.bstar;
    elements_.epoch = tle.epoch;

    // Calculate semi-major axis (in Earth radii)
    const double a1 = pow(XKE / xno, 2.0/3.0);
    const double cosio = cos(elements_.inclination);
    const double theta2 = cosio * cosio;
    const double x3thm1 = 3.0 * theta2 - 1.0;
    const double eosq = elements_.eccentricity * elements_.eccentricity;
    const double betao2 = 1.0 - eosq;
    const double betao = sqrt(betao2);

    // SGP4 initialization
    const double del1 = 1.5 * CK2 * x3thm1 / (a1 * a1 * betao * betao2);
    const double ao = a1 * (1.0 - del1 * (0.5 + del1 * (1.0 + 134.0/81.0 * del1)));
    const double delo = 1.5 * CK2 * x3thm1 / (ao * ao * betao * betao2);

    elements_.a = ao * (1.0 - delo);  // in Earth radii
    elements_.n0 = xno / (1.0 + delo);
    elements_.n = elements_.n0;
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Time since epoch in minutes
    const double tsince = elements_.epoch.secsTo(time) / 60.0;

    // Update for secular gravity and atmospheric drag
    const double xmp = elements_.mean_anomaly + elements_.n * tsince;
    const double omega = elements_.arg_perigee;
    const double xnode = elements_.right_ascension;

    // Solve Kepler's equation
    double E = xmp;
    double delta;
    for(int i = 0; i < 10; i++) {
        delta = E - elements_.eccentricity * sin(E) - xmp;
        E = E - delta / (1.0 - elements_.eccentricity * cos(E));
        if(fabs(delta) < 1e-12) break;
    }

    // Short-period periodic elements
    const double sin_E = sin(E);
    const double cos_E = cos(E);

    // Position and velocity in orbital plane
    const double r = elements_.a * (1.0 - elements_.eccentricity * cos_E);
    const double v = sqrt(EARTH_MU / (elements_.a * AE_TO_KM));

    const double sinv = (sqrt(1.0 - elements_.eccentricity * elements_.eccentricity) * sin_E) / (1.0 - elements_.eccentricity * cos_E);
    const double cosv = (cos_E - elements_.eccentricity) / (1.0 - elements_.eccentricity * cos_E);

    // Orientation vectors
    const double sin_i = sin(elements_.inclination);
    const double cos_i = cos(elements_.inclination);
    const double sin_omega = sin(omega);
    const double cos_omega = cos(omega);
    const double sin_xnode = sin(xnode);
    const double cos_xnode = cos(xnode);

    // Position in ECI (in km)
    const double rx = r * (cos_omega * cosv - sin_omega * sinv);
    const double ry = r * (sin_omega * cosv + cos_omega * sinv);

    state.position = QVector3D(
        AE_TO_KM * (rx * cos_xnode - ry * cos_i * sin_xnode),
        AE_TO_KM * (rx * sin_xnode + ry * cos_i * cos_xnode),
        AE_TO_KM * (ry * sin_i)
        );

    // Velocity in ECI (in km/s)
    const double vfac = sqrt(EARTH_MU / (r * AE_TO_KM));
    const double vr = vfac * elements_.eccentricity * sin_E;
    const double vt = vfac * sqrt(1.0 - elements_.eccentricity * elements_.eccentricity);

    state.velocity = QVector3D(
        (vr * cos_xnode - vt * sin_xnode) / TIME_TO_SEC,
        (vr * sin_xnode + vt * cos_xnode) / TIME_TO_SEC,
        0.0
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

