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
    static constexpr double XKMPER = 6378.137;
    static constexpr double AE = 1.0;
    static constexpr double DE2RA = M_PI / 180.0;
    static constexpr double MINUTES_PER_DAY = 1440.0;
    static constexpr double OMEGA_E = 7.29211514670698e-5;

    struct Elements {
        double i;          // Наклонение (рад)
        double Omega;      // Долгота восходящего узла (рад)
        double e;          // Эксцентриситет
        double omega;      // Аргумент перигея (рад)
        double M;          // Средняя аномалия (рад)
        double n;          // Среднее движение (рад/мин)
        double bstar;      // Баллистический коэффициент
        QDateTime epoch;   // Эпоха элементов

        // Производные параметры
        double a;          // Большая полуось (км)
        double n0;         // Начальное среднее движение
        double cosio;      // cos(i)
        double sinio;      // sin(i)
    };

    void initializeParameters(const TLEParser::TLEData& tle);
    double solveKeplerEquation(double M, double e, int maxIter = 10) const;
    QVector3D calculatePosition(double xnode, double xinc, double xll, double e, double omega, double xn) const;

    Elements elements_;
};

#endif // SGP4_PROPAGATOR_H
