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
    static constexpr double OMEGA_E = 7.29211514670698e-5;  // Скорость вращения Земли (рад/с)
    static constexpr double QOMS2T = 1.880279e-09;
    static constexpr double S = 1.01222928;
    static constexpr double TOTHRD = 2.0/3.0;

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
        double a;              // Большая полуось (в земных радиусах)
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
        double xlcof;         // Коэффициент для долгопериодических возмущений
        double aycof;         // Коэффициент для годичных возмущений
        double delmo;         // Коэффициент для средней аномалии
        double sinmo;         // sin(M)
        double x2o3;          // 2/3
        double aodp;          // Параметр орбиты
        double xmdot;         // Скорость изменения средней аномалии
        double omgdot;        // Скорость изменения аргумента перигея
        double xnodot;        // Скорость изменения долготы восходящего узла
    };

    void initializeParameters(const TLEParser::TLEData& tle);
    void calculatePerturbations(double tsince, double& ep, double& xincp,
                                double& omgadf, double& xnode, double& xmp) const;
    void solveKeplerEquation(double xme, double xec, double& eps, double& epw) const;

    Elements elements_;
};

#endif // SGP4_PROPAGATOR_H
