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
    static constexpr double TOTHRD = 2.0/3.0;
    static constexpr double J2 = 1.082616e-3;
    static constexpr double J3 = -2.53881e-6;
    static constexpr double J4 = -1.65597e-6;

    struct Elements {
        double i;              // Наклонение (рад)
        double Omega;          // Долгота восходящего узла (рад)
        double e;              // Эксцентриситет
        double omega;          // Аргумент перигея (рад)
        double M;              // Средняя аномалия (рад)
        double n;              // Среднее движение (рад/мин)
        double bstar;          // Баллистический коэффициент
        QDateTime epoch;       // Эпоха

        // Вспомогательные параметры
        double a;              // Большая полуось
        double cosio;          // cos(i)
        double sinio;          // sin(i)
        double eta;            // sqrt(1 - e^2)
        double con41;          // Коэффициент для периодических возмущений
        double x1mth2;         // 1 - theta^2
        double x3thm1;         // 3 * theta^2 - 1
        double x7thm1;         // 7 * theta^2 - 1
        double xlcof;          // Коэффициент для долгопериодических возмущений
        double aycof;          // Коэффициент для годовых возмущений
    };

    void initializeParameters(const TLEParser::TLEData& tle);
    void calculateStateVectors(double tsince, QVector3D& pos, QVector3D& vel) const;

    Elements elements_;
};

#endif // SGP4_PROPAGATOR_H
