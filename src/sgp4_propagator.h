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
    static constexpr double SQRT_GM = 0.0743669161331734049; // sqrt(GM) in ER^(3/2)/min

    struct OrbitalElements {
        // Основные элементы
        double inclo;    // Наклонение (рад)
        double nodeo;    // RA восходящего узла (рад)
        double ecco;     // Эксцентриситет
        double argpo;    // Аргумент перигея (рад)
        double mo;       // Средняя аномалия (рад)
        double no;       // Среднее движение (рад/мин)
        double bstar;    // Баллистический коэффициент

        // Вспомогательные элементы
        double a;        // Большая полуось (земные радиусы)
        double alta;     // Апогей (км)
        double altp;     // Перигей (км)
        double jdsatepoch; // Юлианская дата эпохи
    };

    void initializeParameters(const TLEParser::TLEData& tle);
    void sgp4init();
    void sgp4(double tsince, QVector3D& pos, QVector3D& vel) const;
    void rv2coe(const QVector3D& pos, const QVector3D& vel, double& a, double& e,
                double& i, double& omega, double& argp, double& nu) const;

    OrbitalElements elements_;
    QDateTime epoch_;
};

#endif // SGP4_PROPAGATOR_H
