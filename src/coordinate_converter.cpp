// coordinate_converter.cpp
#include "coordinate_converter.h"
#include <QtMath>

QVector3D CoordinateConverter::eci2ecef(const QVector3D& eci_pos, const QDateTime& time) {
    double gmst = calculateGMST(time);
    double cos_gmst = qCos(gmst);
    double sin_gmst = qSin(gmst);

    return QVector3D(
        eci_pos.x() * cos_gmst + eci_pos.y() * sin_gmst,
        -eci_pos.x() * sin_gmst + eci_pos.y() * cos_gmst,
        eci_pos.z()
        );
}

CoordinateConverter::GeodeticCoords
CoordinateConverter::ecef2geodetic(const QVector3D& ecef_pos) {
    double e2 = 2.0 * EARTH_FLATTENING - EARTH_FLATTENING * EARTH_FLATTENING;
    double r = qSqrt(ecef_pos.x() * ecef_pos.x() + ecef_pos.y() * ecef_pos.y());
    double E2 = EARTH_RADIUS * EARTH_RADIUS - (EARTH_RADIUS * (1.0 - EARTH_FLATTENING)) *
                                                  (EARTH_RADIUS * (1.0 - EARTH_FLATTENING));
    double F = 54.0 * (EARTH_RADIUS * (1.0 - EARTH_FLATTENING)) * (EARTH_RADIUS * (1.0 - EARTH_FLATTENING)) *
               ecef_pos.z() * ecef_pos.z();
    double G = r * r + (1.0 - e2) * ecef_pos.z() * ecef_pos.z() - e2 * E2;
    double c = e2 * e2 * F * r * r / (G * G * G);
    double s = qPow(1.0 + c + qSqrt(c * c + 2.0 * c), 1.0/3.0);
    double P = F / (3.0 * (s + 1.0/s + 1.0) * (s + 1.0/s + 1.0) * G * G);
    double Q = qSqrt(1.0 + 2.0 * e2 * e2 * P);
    double r0 = -P * e2 * r / (1.0 + Q) +
                qSqrt(0.5 * EARTH_RADIUS * EARTH_RADIUS * (1.0 + 1.0/Q) -
                      P * (1.0 - e2) * ecef_pos.z() * ecef_pos.z()/(Q * (1.0 + Q)) - 0.5 * P * r * r);
    double U = qSqrt((r - e2 * r0) * (r - e2 * r0) + ecef_pos.z() * ecef_pos.z());
    double V = qSqrt((r - e2 * r0) * (r - e2 * r0) + (1.0 - e2) * ecef_pos.z() * ecef_pos.z());
    double z0 = (EARTH_RADIUS * (1.0 - EARTH_FLATTENING)) * (EARTH_RADIUS * (1.0 - EARTH_FLATTENING)) *
                ecef_pos.z() / (EARTH_RADIUS * V);

    GeodeticCoords coords;
    coords.latitude = qAtan((ecef_pos.z() + e2 * z0) / r) * 180.0 / M_PI;
    coords.longitude = qAtan2(ecef_pos.y(), ecef_pos.x()) * 180.0 / M_PI;
    coords.altitude = U * (1.0 - (EARTH_RADIUS * (1.0 - EARTH_FLATTENING)) *
                                     (EARTH_RADIUS * (1.0 - EARTH_FLATTENING)) / (EARTH_RADIUS * V));

    return coords;
}

CoordinateConverter::GeodeticCoords
CoordinateConverter::eci2geodetic(const QVector3D& eci_pos, const QDateTime& time) {
    QVector3D ecef_pos = eci2ecef(eci_pos, time);
    return ecef2geodetic(ecef_pos);
}

double CoordinateConverter::calculateGMST(const QDateTime& time) {
    // Расчет GMST (Greenwich Mean Sidereal Time)
    QDateTime epoch(QDate(2000, 1, 1), QTime(12, 0, 0), Qt::UTC);
    double days_since_epoch = epoch.daysTo(time) +
                              epoch.secsTo(time) / 86400.0;

    double T = days_since_epoch / 36525.0;
    double gmst = 280.46061837 + 360.98564736629 * days_since_epoch +
                  0.000387933 * T * T - T * T * T / 38710000.0;

    // Нормализация к диапазону 0-360 градусов
    gmst = fmod(gmst, 360.0);
    if (gmst < 0) gmst += 360.0;

    return gmst * M_PI / 180.0; // конвертация в радианы
}
