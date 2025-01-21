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
    // Константы
    static constexpr double XKE = 7.43669161e-2;
    static constexpr double CK2 = 5.413080e-4;
    static constexpr double CK4 = 0.62098875e-6;
    static constexpr double XKMPER = 6378.137;          // Радиус Земли (км)
    static constexpr double AE = 1.0;
    static constexpr double DE2RA = M_PI / 180.0;       // Градусы в радианы
    static constexpr double MINUTES_PER_DAY = 1440.0;
    static constexpr double OMEGA_E = 7.29211514670698e-5;
    static constexpr double XMNPDA = 1440.0;            // Минут в день
    static constexpr double MFACTOR = 7.292115E-5;      // Угловая скорость Земли

    struct Elements {
        double a;          // Большая полуось (км)
        double e;          // Эксцентриситет
        double i;          // Наклонение (рад)
        double omega;      // Аргумент перигея (рад)
        double Omega;      // Долгота восходящего узла (рад)
        double M;          // Средняя аномалия (рад)
        double n;          // Среднее движение (рад/мин)
        double bstar;      // Баллистический коэффициент
        QDateTime epoch;   // Эпоха
    };

    void initializeParameters(const TLEParser::TLEData& tle);
    QVector3D calculatePositionAndVelocity(double tsince) const;
    double solveKepler(double meanAnomaly, double eccentricity) const;

    Elements elements_;
};

#endif // SGP4_PROPAGATOR_H
