#ifndef SGP4_PROPAGATOR_H
#define SGP4_PROPAGATOR_H

#pragma once

#include "tle_parser.h"
#include <QVector3D>
#include <QDateTime>

class SGP4Propagator {
public:
    struct OrbitalState {
        QVector3D position;     // км в системе ECI (Earth-Centered Inertial)
        QVector3D velocity;     // км/с в системе ECI
        QDateTime epoch;        // время расчета
    };
    struct GeodeticCoords {
        double latitude;    // градусы (положительные для северной широты)
        double longitude;   // градусы (положительные для восточной долготы)
        double altitude;    // км над поверхностью Земли
    };

    explicit SGP4Propagator(const TLEParser::TLEData& tle);
    OrbitalState calculateState(const QDateTime& time) const;
    // Преобразование из ECI в ECEF
    static QVector3D eci2ecef(const QVector3D& eci_pos, const QDateTime& time);
    // Преобразование из ECEF в геодезические координаты
    static GeodeticCoords ecef2geodetic(const QVector3D& ecef_pos);
    // Прямое преобразование из ECI в геодезические координаты
    static GeodeticCoords eci2geodetic(const QVector3D& eci_pos, const QDateTime& time);

private:
    // WGS-72 константы
    static constexpr double XKE = 7.43669161E-2;
    static constexpr double CK2 = 5.413080E-4;
    static constexpr double CK4 = 0.62098875E-6;
    static constexpr double QOMS2T = 1.880279E-09;
    static constexpr double S = 1.01222928;
    static constexpr double XKMPER = 6378.135;    // радиус Земли в км
    static constexpr double AE = 1.0;
    static constexpr double DE2RA = 0.174532925E-1;
    static constexpr double MINUTES_PER_DAY = 1440.0;
    static constexpr double MU = 398600.8;        // гравитационный параметр Земли
    static constexpr double EARTH_RADIUS = 6378.137;    // км (WGS-84)
    static constexpr double EARTH_FLATTENING = 1.0/298.257223563;
    static constexpr double EARTH_ROTATION_RATE = 7.2921150e-5; // рад/с

    struct OrbitalElements {
        double inclination;     // радианы
        double right_ascension; // радианы
        double eccentricity;
        double arg_perigee;     // радианы
        double mean_anomaly;    // радианы
        double mean_motion;     // радианы в минуту
        double bstar;
        QDateTime epoch;

        // Производные элементы
        double a;              // большая полуось
        double n0;            // начальное среднее движение
        double n;             // скорректированное среднее движение
        double e0;            // начальный эксцентриситет
    };

    OrbitalElements elements_;

    void initializeParameters(const TLEParser::TLEData& tle);
    QVector3D getPosition(double r, double u, double i, double omega) const;
    QVector3D getVelocity(double r, double rdot, double u, double rfdot, double i, double omega) const;
};

#endif // SGP4_PROPAGATOR_H
