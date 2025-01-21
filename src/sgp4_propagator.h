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
    static constexpr double XKE = 7.43669161e-2;
    static constexpr double CK2 = 5.413080e-4;
    static constexpr double CK4 = 0.62098875e-6;
    static constexpr double XKMPER = 6378.135;          // Earth radius km
    static constexpr double AE = 1.0;
    static constexpr double DE2RA = 0.174532925e-1;     // π/180
    static constexpr double MINUTES_PER_DAY = 1440.0;
    static constexpr double OMEGA_E = 7.29211514670698e-5; // Earth rotation rad/s
    static constexpr double QOMS2T = 1.880279e-09;
    static constexpr double S = 1.01222928;

    struct SGP4Elements {
        double a;          // Semi-major axis (Earth radii)
        double e;          // Eccentricity
        double i;          // Inclination (radians)
        double omega;      // Argument of perigee (radians)
        double Omega;      // Right ascension (radians)
        double M;          // Mean anomaly (radians)
        double n;          // Mean motion (radians/minute)
        double bstar;      // Drag term

        // Вспомогательные элементы
        double n0;         // Original mean motion
        double a0;         // Original semi major axis
        double ndot;       // First time derivative of mean motion
        double nddot;      // Second time derivative of mean motion
    };

    SGP4Elements sgp4_;

    // Вспомогательные функции
    static double Kepler(double MeanAnomaly, double Eccentricity) {
        double E = MeanAnomaly;
        double delta;
        int iter = 0;
        do {
            delta = E - Eccentricity * sin(E) - MeanAnomaly;
            E = E - delta / (1.0 - Eccentricity * cos(E));
            iter++;
        } while (fabs(delta) > 1e-12 && iter < 10);
        return E;
    }

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
