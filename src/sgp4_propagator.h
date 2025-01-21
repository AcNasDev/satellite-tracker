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
    static constexpr double ke = 7.43669161e-2;
    static constexpr double k2 = 5.413080e-4;
    static constexpr double k4 = 0.62098875e-6;
    static constexpr double xkmper = 6378.137;
    static constexpr double min_per_day = 1440.0;
    static constexpr double ae = 1.0;
    static constexpr double de2ra = M_PI / 180.0;

    struct Elements {
        double inclo;     // Наклонение (рад)
        double nodeo;     // Долгота восходящего узла (рад)
        double ecco;      // Эксцентриситет
        double argpo;     // Аргумент перигея (рад)
        double mo;        // Средняя аномалия (рад)
        double no;        // Среднее движение (рад/мин)
        double bstar;     // Баллистический коэффициент

        // Вспомогательные элементы
        double a;         // Большая полуось (км)
        double n;         // Скорректированное среднее движение
        double e;         // Скорректированный эксцентриситет
        double i;         // Скорректированное наклонение
        double omega;     // Скорректированный аргумент перигея
        double Omega;     // Скорректированная долгота узла
    };

    void initParameters(const TLEParser::TLEData& tle);
    void propagate(double tsince, QVector3D& pos, QVector3D& vel) const;
    void debugElements() const;

    Elements elements_;
    QDateTime epoch_;
};

#endif // SGP4_PROPAGATOR_H
