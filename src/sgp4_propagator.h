#ifndef SGP4_PROPAGATOR_H
#define SGP4_PROPAGATOR_H

#pragma once

#include "tle_parser.h"
#include <QVector3D>
#include <QDateTime>

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
    static constexpr double ke = 7.43669161e-2;        // sqrt(GM) ER^(3/2)/min
    static constexpr double k2 = 5.413080e-4;          // J2/2
    static constexpr double k4 = 0.62098875e-6;        // -3J4/8
    static constexpr double a3ovk2 = -1.0 / 3.0;       // J3/(3*J2)
    static constexpr double xkmper = 6378.137;         // Радиус Земли в км
    static constexpr double min_per_day = 1440.0;      // Минут в сутках
    static constexpr double ae = 1.0;                  // Расстояние в радиусах Земли
    static constexpr double de2ra = M_PI / 180.0;      // Градусы в радианы

    struct Elements {
        double inclo;     // Наклонение (рад)
        double nodeo;     // Долгота восходящего узла (рад)
        double ecco;      // Эксцентриситет
        double argpo;     // Аргумент перигея (рад)
        double mo;        // Средняя аномалия (рад)
        double no;        // Среднее движение (рад/мин)
        double bstar;     // Баллистический коэффициент

        // Производные элементы
        double a;         // Большая полуось (ER)
        double ndot;      // Изменение среднего движения (рад/мин^2)
        double nddot;     // Вторая производная среднего движения (рад/мин^3)
        double alta;      // Высота апогея (км)
        double altp;      // Высота перигея (км)
        double del1;      // Первая поправка к большой полуоси
        double del2;      // Вторая поправка к большой полуоси
        double del3;      // Третья поправка к большой полуоси
        double xincl;     // Наклонение с учетом возмущений
        double xnodp;     // Узел с учетом возмущений
        double aodp;      // Большая полуось с учетом возмущений
        double ao, delo;  // Оригинальные элементы
    };

    void initParameters(const TLEParser::TLEData& tle);
    void calculateDerivatives();
    void propagate(double tsince, QVector3D& pos, QVector3D& vel) const;

    Elements elements_;
    QDateTime epoch_;
};

#endif // SGP4_PROPAGATOR_H
