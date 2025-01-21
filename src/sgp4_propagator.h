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
    // Физические константы
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
    static constexpr double A3OVK2 = -J3/CK2;           // Отношение J3/J2
    static constexpr double MU = 398600.4418;           // Гравитационный параметр Земли (км³/с²)

    struct OrbitalElements {
        double a;              // Большая полуось (км)
        double e;              // Эксцентриситет
        double i;              // Наклонение (рад)
        double omega;          // Аргумент перигея (рад)
        double Omega;          // Долгота восходящего узла (рад)
        double M;              // Средняя аномалия (рад)
        double n;              // Среднее движение (рад/мин)
        double bstar;          // Баллистический коэффициент
        QDateTime epoch;       // Эпоха элементов

        // Вспомогательные элементы
        double cosio;          // cos(i)
        double sinio;          // sin(i)
        double eta;            // sqrt(1 - e^2)
        double coef;          // Коэффициент для драг-эффекта
        double c1;            // Коэффициент J2
        double c4;            // Коэффициент для периодических возмущений
        double a0;            // Начальная большая полуось
        double n0;            // Начальное среднее движение
        double beta0;         // sqrt(1 - e^2)
        double theta2;        // cos²(i)
        double x3thm1;        // 3cos²(i) - 1
        double x7thm1;        // 7cos²(i) - 1
    };

    void initializeParameters(const TLEParser::TLEData& tle);
    void calculateSecularEffects(double tsince, double& xmdf, double& omgadf,
                                 double& xnode, double& eps1) const;
    void calculatePeriodicEffects(double tsince, double& em, double& xinc,
                                  double& omgda, double& xnode, double& xll) const;
    QVector3D calculatePosVel(double r, double rdot, double rfdot,
                              double cosuk, double sinuk, double rk,
                              double uk, double xnodek, double xinck) const;

    OrbitalElements elements_;
};

#endif // SGP4_PROPAGATOR_H
