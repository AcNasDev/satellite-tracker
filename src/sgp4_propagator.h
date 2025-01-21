#ifndef SGP4_PROPAGATOR_H
#define SGP4_PROPAGATOR_H

#pragma once

#include "tle_parser.h"
#include <QVector3D>
#include <QDateTime>

class SGP4Propagator {
public:
    struct OrbitalState {
        QVector3D position;     // км в системе ECI (Earth-Centered Inertial)
        QVector3D velocity;     // км/с в системе ECI
        QDateTime epoch;        // время расчета
    };

    explicit SGP4Propagator(const TLEParser::TLEData& tle);
    OrbitalState calculateState(const QDateTime& time) const;

private:
    // Константы SGP4
    static constexpr double XKE = 7.43669161e-2;        // sqrt(398600.4418 / earthRadius^3)
    static constexpr double CK2 = 5.413080e-4;          // 0.5 * J2 * earthRadius^2
    static constexpr double CK4 = 0.62098875e-6;        // -0.375 * J4 * earthRadius^4
    static constexpr double XKMPER = 6378.135;          // Earth radius km
    static constexpr double AE = 1.0;                   // Distance units/Earth radii
    static constexpr double DE2RA = 0.174532925e-1;     // π/180 (degrees to radians)
    static constexpr double MINUTES_PER_DAY = 1440.0;
    static constexpr double OMEGA_E = 7.29211514670698e-5; // Earth rotation rad/s
    static constexpr double QOMS2T = 1.880279e-09;      // (q0 - s)^4, q0 = 120 км
    static constexpr double S = 1.01222928;             // для расчета атмосферного торможения

    struct SGP4Elements {
        double a;          // Semi-major axis (Earth radii)
        double e;          // Eccentricity
        double i;          // Inclination (radians)
        double omega;      // Argument of perigee (radians)
        double Omega;      // Right ascension (radians)
        double M;          // Mean anomaly (radians)
        double n;          // Mean motion (radians/minute)
        double bstar;      // Drag term
        QDateTime epoch;   // Epoch time

        // Производные элементы
        double n0;         // Original mean motion
        double a0;         // Original semi major axis
        double ndot;       // First time derivative of mean motion
        double nddot;      // Second time derivative of mean motion
        double cosio;      // Cosine of inclination
        double sinio;      // Sine of inclination
        double eta;        // sqrt(1 - e^2)
        double coef;       // Coefficient for secular perturbations
        double c1;         // Coefficient for long-period perturbations
        double c4;         // Coefficient for periodic perturbations
        double x3thm1;     // 3 * cos^2(i) - 1
        double x1mth2;     // 1 - cos^2(i)
        double x7thm1;     // 7 * cos^2(i) - 1
    };

    void initializeParameters(const TLEParser::TLEData& tle);
    void calculateSecularEffects(double tsince, double& xmdf, double& omgadf, double& xnode) const;
    void calculatePeriodicEffects(double tsince, double& ep, double& xincp,
                                  double& omgadf, double& xnode) const;
    QVector3D calculatePositionAndVelocity(double rk, double uk, double xnodek,
                                           double xinck, double rdotk, double rfdotk) const;

    SGP4Elements sgp4_;
};

#endif // SGP4_PROPAGATOR_H
