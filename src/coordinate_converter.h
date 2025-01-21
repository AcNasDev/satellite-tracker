#ifndef COORDINATE_CONVERTER_H
#define COORDINATE_CONVERTER_H

#pragma once

#include <QVector3D>
#include <QDateTime>

class CoordinateConverter {
public:
    struct GeodeticCoords {
        double latitude;    // градусы (положительные для северной широты)
        double longitude;   // градусы (положительные для восточной долготы)
        double altitude;    // км над поверхностью Земли
    };

    // Преобразование из ECI в ECEF
    static QVector3D eci2ecef(const QVector3D& eci_pos, const QDateTime& time);

    // Преобразование из ECEF в геодезические координаты
    static GeodeticCoords ecef2geodetic(const QVector3D& ecef_pos);

    // Прямое преобразование из ECI в геодезические координаты
    static GeodeticCoords eci2geodetic(const QVector3D& eci_pos, const QDateTime& time);

private:
    static constexpr double EARTH_RADIUS = 6378.137;    // км (WGS-84)
    static constexpr double EARTH_FLATTENING = 1.0/298.257223563;
    static constexpr double EARTH_ROTATION_RATE = 7.2921150e-5; // рад/с

    // Получение GMST (Greenwich Mean Sidereal Time)
    static double calculateGMST(const QDateTime& time);
};

#endif // COORDINATE_CONVERTER_H
