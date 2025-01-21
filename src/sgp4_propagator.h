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
    // Константы
    static constexpr double XKE = 7.43669161e-2;
    static constexpr double CK2 = 5.413080e-4;
    static constexpr double CK4 = 0.62098875e-6;
    static constexpr double XKMPER = 6378.135;         // Earth radius km
    static constexpr double AE = 1.0;
    static constexpr double XMNPDA = 1440.0;          // Minutes per day
    static constexpr double OMEGA_E = 7.29211514670698e-5;  // Earth rotation rate rad/s
    static constexpr double TWOPI = 2.0 * M_PI;
    static constexpr double E6A = 1.0e-6;
    static constexpr double TOTHRD = 2.0/3.0;

    struct DeepSpaceElements {
        double aodp, cosio, sinio, omgdot, xnodot, xmdot;
        double xnodp, theta2, x3thm1, xlcof, aycof;
        double x1mth2, xmdot_factor, omgdot_factor, xnodot_factor;
    };

    DeepSpaceElements deep_;

    // Вспомогательные функции
    static double FMod2p(double x) {
        double rv = fmod(x, TWOPI);
        if (rv < 0.0) rv += TWOPI;
        return rv;
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
