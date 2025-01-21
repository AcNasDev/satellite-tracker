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
    static constexpr double J2 = 0.00108262998905;
    static constexpr double J3 = -0.00000253215306;
    static constexpr double J4 = -0.00000161098761;
    static constexpr double KE = 0.0743669161331734049;
    static constexpr double QOMS2T = 1.880279e-09;

    struct Elements {
        // Основные элементы
        double i;              // Наклонение (рад)
        double Omega;          // Долгота восходящего узла (рад)
        double e;              // Эксцентриситет
        double omega;          // Аргумент перигея (рад)
        double M;              // Средняя аномалия (рад)
        double n;              // Среднее движение (рад/мин)
        double bstar;          // Баллистический коэффициент
        QDateTime epoch;       // Эпоха

        // Производные элементы
        double a;              // Большая полуось (земные радиусы)
        double n0;             // Начальное среднее движение
        double e0;             // Начальный эксцентриситет
        double cosio;          // cos(i)
        double sinio;          // sin(i)
        double eta;            // sqrt(1 - e^2)
        double coef;           // Коэффициент для атмосферного сопротивления
        double C1;             // Коэффициент для J2
        double C4;             // Коэффициент для возмущений
        double x3thm1;         // 3 * cosio^2 - 1
        double x1mth2;         // 1 - cosio^2
        double xmdot;          // Скорость изменения средней аномалии
        double omgdot;         // Скорость изменения аргумента перигея
        double xnodot;         // Скорость изменения долготы восходящего узла
        double xnodcf;         // = 3.5 * beta0 * C1
        double t2cof;          // = 1.5 * C1
        double aodp;           // Параметр орбиты
        double beta0;          // = sqrt(1 - e^2)
    };

    void initializeParameters(const TLEParser::TLEData& tle);
    void updateForSecularEffects(double tsince, double& xll,
                                 double& omgasm, double& xnodes,
                                 double& em, double& xinc, double& xn) const;
    void updateForPeriodicEffects(double& em, double& xinc, double& omgasm,
                                  double& xnodes, double& xll, double tsince) const;
    void solveKeplerEquation(double& xll, double e) const;

    Elements elements_;
};

#endif // SGP4_PROPAGATOR_H
