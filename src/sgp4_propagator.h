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
    static constexpr double XKE = 0.0743669161331734049;
    static constexpr double CK2 = 5.413080e-4;
    static constexpr double CK4 = 0.62098875e-6;
    static constexpr double XKMPER = 6378.137;
    static constexpr double AE = 1.0;
    static constexpr double DE2RA = M_PI / 180.0;
    static constexpr double MINUTES_PER_DAY = 1440.0;
    static constexpr double QOMS2T = 1.880279e-09;
    static constexpr double S = 1.01222928;

    struct Elements {
        // Основные элементы орбиты
        double i;              // Наклонение (рад)
        double Omega;          // Долгота восходящего узла (рад)
        double e;              // Эксцентриситет
        double omega;          // Аргумент перигея (рад)
        double M;              // Средняя аномалия (рад)
        double n;              // Среднее движение (рад/мин)
        double bstar;          // Баллистический коэффициент
        QDateTime epoch;       // Эпоха

        // Производные параметры
        double a;              // Большая полуось (в радиусах Земли)
        double n0;             // Исходное среднее движение
        double cosio;          // cos(i)
        double sinio;          // sin(i)
        double eta;            // sqrt(1 - e^2)
        double coef;           // Коэффициент для драг-эффекта
        double c1;            // Коэффициент для J2
        double c4;            // Коэффициент для периодических возмущений
        double theta2;        // cos²(i)
        double x3thm1;        // 3*cos²(i) - 1
        double x1mth2;        // 1 - cos²(i)
        double x7thm1;        // 7*cos²(i) - 1
    };

    void initializeParameters(const TLEParser::TLEData& tle);
    void solveKeplerEquation(double& meanAnomaly, double& eccentricAnomaly) const;
    QVector3D calculatePosVel(double tsince) const;

    Elements elements_;
};

#endif // SGP4_PROPAGATOR_H
