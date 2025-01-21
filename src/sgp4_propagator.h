// sgp4_propagator.h
#ifndef SGP4_PROPAGATOR_H
#define SGP4_PROPAGATOR_H

#pragma once

#include "tle_parser.h"
#include <QVector3D>
#include <QDateTime>
#include <QDebug>

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
    // Константы SGP4
    static constexpr double ke = 0.0743669161331734049;  // √(GM) ER^(3/2)/min
    static constexpr double k2 = 5.413080e-4;            // J2/2
    static constexpr double k4 = 0.62098875e-6;          // -3J4/8
    static constexpr double xkmper = 6378.137;           // Радиус Земли (км)
    static constexpr double min_per_day = 1440.0;
    static constexpr double ae = 1.0;                    // Радиус Земли в ER
    static constexpr double de2ra = M_PI / 180.0;

    struct Elements {
        // Базовые элементы
        double n0;        // Начальное среднее движение (рад/мин)
        double e0;        // Начальный эксцентриситет
        double i0;        // Начальное наклонение (рад)
        double omega0;    // Начальный аргумент перигея (рад)
        double Omega0;    // Начальный RAAN (рад)
        double M0;        // Начальная средняя аномалия (рад)
        double bstar;     // Баллистический коэффициент

        // Производные элементы
        double a0;        // Большая полуось (в единицах Земных радиусов)
        double n0dot;     // Первая производная среднего движения
        double n0ddot;    // Вторая производная среднего движения
    };

    void initParameters(const TLEParser::TLEData& tle);
    void propagate(double tsince, QVector3D& pos, QVector3D& vel) const;
    double solveKepler(double M, double e) const;

    Elements elements_;
    QDateTime epoch_;
};

#endif // SGP4_PROPAGATOR_H
