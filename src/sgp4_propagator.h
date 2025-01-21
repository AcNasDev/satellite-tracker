#ifndef SGP4_PROPAGATOR_H
#define SGP4_PROPAGATOR_H

#pragma once

#include "tle_parser.h"
#include <QVector3D>
#include <QDateTime>

class SGP4Propagator {
public:
    struct OrbitalState {
        QVector3D position;     // км в ECI
        QVector3D velocity;     // км/с в ECI
        QDateTime epoch;        // время расчета
    };

    explicit SGP4Propagator(const TLEParser::TLEData& tle);
    OrbitalState calculateState(const QDateTime& time) const;

private:
    // SGP4 константы
    static constexpr double ae = 1.0;                    // Earth radius in ER units
    static constexpr double tothrd = 2.0 / 3.0;
    static constexpr double xke = 0.0743669161331734049; // sqrt(GM) in ER^(3/2)/min
    static constexpr double ck2 = 5.413080e-4;           // J2/2
    static constexpr double ck4 = 0.62098875e-6;         // -3J4/8
    static constexpr double j3 = -0.253881e-5;           // J3
    static constexpr double xkmper = 6378.137;           // Earth radius in km
    static constexpr double minute_per_day = 1440.0;
    static constexpr double ae_to_km = 6378.137;
    static constexpr double s = 1.012229;                // AE/RE ratio

    struct SGP4Elements {
        double a;          // Semi-major axis (Earth radii)
        double ecco;      // Eccentricity
        double inclo;     // Inclination (radians)
        double nodeo;     // RAAN (radians)
        double argpo;     // Argument of perigee (radians)
        double mo;        // Mean anomaly (radians)
        double no;        // Mean motion (radians/minute)
        double bstar;     // Drag term
        double aycof;

        // Вычисляемые параметры
        double alta, altp, a0, d2, d3, d4, del1, del2, del3;
        double eta, argpdot, omgcof, sinmao, t2cof, t3cof, t4cof, t5cof;
        double x1mth2, x7thm1, xlcof, xmcof, nodecf, nodedot, xnodot;
        double e0, mdot;
    };

    void initParameters(const TLEParser::TLEData& tle);
    void calculateSecularEffects(double tsince, double& xll, double& omgasm,
                                 double& xnodes, double& em, double& xinc,
                                 double& xn) const;
    void calculatePeriodicEffects(double tsince, double& em, double& xinc,
                                  double& omgasm, double& xnodes, double& xll) const;
    void solveKepler(double& xll, double e) const;

    SGP4Elements elements_;
    QDateTime epoch_;
};

#endif // SGP4_PROPAGATOR_H
