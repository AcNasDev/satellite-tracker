#ifndef SATELLITETRACKER_H
#define SATELLITETRACKER_H

#pragma once

#include <QObject>
#include <QTimer>
#include <QDateTime>
#include <QDebug>
#include <memory>

#include "SGP4.h"
#include "CoordGeodetic.h"
#include "CoordTopocentric.h"
#include "Observer.h"
#include "Tle.h"

#include "tle_parser.h"
#include "sgp4_propagator.h"
#include "coordinate_converter.h"

class SatelliteTracker : public QObject
{
    Q_OBJECT

public:
    explicit SatelliteTracker(const QString& line1, const QString& line2, QObject* parent = nullptr)
        : QObject(parent)
    {
        data = TLEParser::parseTLE(line1, line2).value();
        propagator = std::make_unique<SGP4Propagator>(data);

        try {
            tle = std::make_unique<libsgp4::Tle>("Satellite", line1.toStdString(), line2.toStdString());
            //добавь строчку для вывода в консоль всю информацию из структуры  Tle

            sgp4 = std::make_unique<libsgp4::SGP4>(*tle);

            QTimer* timer = new QTimer(this);
            connect(timer, &QTimer::timeout, this, &SatelliteTracker::updatePosition);
            timer->start(1000);

            updatePosition();

        } catch (const libsgp4::TleException& e) {
            qDebug() << "TLE Error:" << e.what();
            throw;
        }
    }

private slots:
    void updatePosition() {
        // Получаем текущее время в UTC
        QDateTime now = QDateTime::currentDateTimeUtc();

        libsgp4::DateTime dt(
            now.date().year(),
            now.date().month(),
            now.date().day(),
            now.time().hour(),
            now.time().minute(),
            now.time().second()
            );

        libsgp4::Eci eci = sgp4->FindPosition(dt);
        libsgp4::CoordGeodetic geo = eci.ToGeodetic();

        qDebug() << "libsgp4 results:";
        qDebug() << "\nCurrent satellite position (UTC):" << now.toString("yyyy-MM-dd HH:mm:ss");
        qDebug() << "ECF coordinates (km):";
        qDebug() << "X:" << QString::number(eci.Position().x, 'f', 3);
        qDebug() << "Y:" << QString::number(eci.Position().y, 'f', 3);
        qDebug() << "Z:" << QString::number(eci.Position().z, 'f', 3);
        qDebug() << "\nGeographic coordinates:";
        qDebug() << "Latitude:" << QString::number(geo.latitude * 180.0 / M_PI, 'f', 6) << "degrees";
        qDebug() << "Longitude:" << QString::number(geo.longitude * 180.0 / M_PI, 'f', 6) << "degrees";
        qDebug() << "Altitude:" << QString::number(geo.altitude, 'f', 3) << "km";
        qDebug() << "\nVelocity (km/s):";
        qDebug() << "X:" << QString::number(eci.Velocity().x, 'f', 6);
        qDebug() << "Y:" << QString::number(eci.Velocity().y, 'f', 6);
        qDebug() << "Z:" << QString::number(eci.Velocity().z, 'f', 6);
        qDebug() << "----------------------------------------";

        qDebug() << "Propagator results:";
        SGP4Propagator::OrbitalState state = propagator->calculateState(now);
        qDebug() << "ECI coordinates:";
        qDebug() << "Position (km):" << state.position;
        qDebug() << "Velocity (km/s):" << state.velocity;

        // Преобразование в ECEF
        QVector3D ecef_pos = CoordinateConverter::eci2ecef(state.position, now);
        qDebug() << "\nECEF coordinates:";
        qDebug() << "Position (km):" << ecef_pos;

        // Преобразование в геодезические координаты
        auto geodetic = CoordinateConverter::ecef2geodetic(ecef_pos);
        qDebug() << "\nGeodetic coordinates:";
        qDebug() << "Latitude:" << geodetic.latitude << "degrees";
        qDebug() << "Longitude:" << geodetic.longitude << "degrees";
        qDebug() << "Altitude:" << geodetic.altitude << "km";

    }

private:
    std::unique_ptr<libsgp4::Tle> tle;
    std::unique_ptr<libsgp4::SGP4> sgp4;

    TLEParser::TLEData data;
    std::unique_ptr<SGP4Propagator> propagator;
};

#endif // SATELLITETRACKER_H
