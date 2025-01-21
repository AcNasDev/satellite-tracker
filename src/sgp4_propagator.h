#ifndef SGP4_PROPAGATOR_H
#define SGP4_PROPAGATOR_H

#pragma once

#include "tle_parser.h"
#include <QVector3D>
#include <QDateTime>

class SGP4Propagator {
public:
    struct OrbitalState {
        QVector3D position;     // км в системе ECI
        QVector3D velocity;     // км/с в системе ECI
        QDateTime epoch;        // время расчета
    };

    explicit SGP4Propagator(const TLEParser::TLEData& tle);
    OrbitalState calculateState(const QDateTime& time) const;

private:
    // Константы SGP4
    static constexpr double xke = 0.0743669161331734049;
    static constexpr double ck2 = 5.413080e-4;
    static constexpr double ck4 = 0.62098875e-6;
    static constexpr double e6a = 1.0e-6;
    static constexpr double qoms2t = 1.880279e-09;
    static constexpr double xkmper = 6378.137;
    static constexpr double ae = 1.0;
    static constexpr double s = 1.012229;
    static constexpr double tothrd = 2.0 / 3.0;
    static constexpr double pi = M_PI;
    static constexpr double xj3 = -0.253881e-5;
    static constexpr double xke2 = 7.436685e-2;
    static constexpr double minutes_per_day = 1440.0;
    static constexpr double a3ovk2 = -xj3 / ck2;

    struct SGP4Data {
        double a, altp, alta;
        double epochdays, jdsatepoch;
        double ndot, nddot;
        double bstar, rcse;
        double inclo, nodeo, ecco, argpo, mo, no_kozai;
        double method;
        double aycof, con41, cc1, cc4, cc5;
        double d2, d3, d4;
        double delmo, eta;
        double argpdot, omgcof, sinmao, t2cof, t3cof, t4cof, t5cof;
        double x1mth2, x3thm1, x7thm1, mdot, nodedot, xlcof, xmcof; // Добавлен x3thm1
        double nodecf;
        int isimp;
    };

    void sgp4init(const TLEParser::TLEData& tle);
    QVector3D getPosition(double tsince) const;
    QVector3D getVelocity(double tsince) const;

    SGP4Data satrec_;
    QDateTime epoch_;
};

#endif // SGP4_PROPAGATOR_H
