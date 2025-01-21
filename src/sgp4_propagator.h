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
    // Фундаментальные константы
    static constexpr double XKE = 7.43669161e-2;
    static constexpr double CK2 = 5.413080e-4;
    static constexpr double CK4 = 0.62098875e-6;
    static constexpr double XKMPER = 6378.137;          // Радиус Земли (км)
    static constexpr double AE = 1.0;
    static constexpr double DE2RA = M_PI / 180.0;       // Градусы в радианы
    static constexpr double MINUTES_PER_DAY = 1440.0;
    static constexpr double OMEGA_E = 7.29211514670698e-5; // Угловая скорость вращения Земли (рад/с)
    static constexpr double J2 = 1.082616e-3;           // J2 гармоника
    static constexpr double J3 = -2.53881e-6;           // J3 гармоника
    static constexpr double J4 = -1.65597e-6;           // J4 гармоника

    struct SGP4Elements {
        double a;          // Большая полуось (радиусы Земли)
        double e;          // Эксцентриситет
        double i;          // Наклонение (радианы)
        double omega;      // Аргумент перигея (радианы)
        double Omega;      // Прямое восхождение (радианы)
        double M;          // Средняя аномалия (радианы)
        double n;          // Среднее движение (рад/мин)
        double bstar;      // Баллистический коэффициент
        QDateTime epoch;   // Эпоха

        // Вспомогательные элементы
        double cosio;      // cos(i)
        double sinio;      // sin(i)
        double eta;        // sqrt(1 - e^2)
        double coef;       // Коэффициент для вековых возмущений
        double c1;         // Коэффициент J2
        double a0;         // Начальная большая полуось
        double n0;         // Начальное среднее движение
        double beta0;      // sqrt(1 - e^2)
    };

    void initializeParameters(const TLEParser::TLEData& tle);
    QVector3D calculatePosVel(const double tsince) const;
    double solveKeplerEquation(double M, double e) const;

    SGP4Elements elements_;
};

#endif // SGP4_PROPAGATOR_H
