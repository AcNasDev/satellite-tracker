#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
}

void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    elements_.epoch = tle.epoch;

    // Преобразуем элементы из TLE
    elements_.i = tle.inclination * DE2RA;
    elements_.Omega = tle.right_ascension * DE2RA;
    elements_.e = tle.eccentricity;
    elements_.omega = tle.argument_perigee * DE2RA;
    elements_.M = tle.mean_anomaly * DE2RA;
    elements_.n = tle.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;
    elements_.bstar = tle.bstar;

    // Вычисляем вспомогательные параметры
    elements_.cosio = cos(elements_.i);
    elements_.sinio = sin(elements_.i);
    elements_.x3thm1 = 3.0 * elements_.cosio * elements_.cosio - 1.0;
    elements_.x1mth2 = 1.0 - elements_.cosio * elements_.cosio;
    elements_.x7thm1 = 7.0 * elements_.cosio * elements_.cosio - 1.0;

    double theta2 = elements_.cosio * elements_.cosio;
    double betao = sqrt(1.0 - elements_.e * elements_.e);
    elements_.eta = elements_.e * betao;

    // Вычисляем начальную большую полуось
    double a1 = pow(XKE / elements_.n, TOTHRD);
    double delta1 = 1.5 * CK2 * elements_.x3thm1 / (a1 * a1 * betao * betao * betao);

    double ao = a1 * (1.0 - delta1 * (0.5 * TOTHRD + delta1 *
                                                         (1.0 + 134.0/81.0 * delta1)));
    double delo = 1.5 * CK2 * elements_.x3thm1 / (ao * ao * betao * betao * betao);

    elements_.a = ao * (1.0 - delo);

    if (elements_.x1mth2 != 0.0) {
        elements_.xlcof = 0.125 * J2 * elements_.sinio *
                          (3.0 + 5.0 * elements_.cosio) / (1.0 + elements_.cosio);
    }
    else {
        elements_.xlcof = 0.0;
    }

    elements_.aycof = 0.25 * J2 * elements_.sinio;
    elements_.con41 = 3.0 * elements_.cosio - 1.0;
}

void SGP4Propagator::calculateStateVectors(double tsince, QVector3D& pos, QVector3D& vel) const {
    // Вычисляем вековые возмущения
    double xmdf = elements_.M + elements_.n * tsince;

    double omgadf = elements_.omega +
                    (-3.0 * J2 * elements_.n * elements_.cosio /
                     (2.0 * elements_.a * elements_.a * (1.0 - elements_.e * elements_.e))) * tsince;

    double xnode = elements_.Omega +
                   (-3.0 * J2 * elements_.n * elements_.cosio /
                    (2.0 * elements_.a * elements_.a * (1.0 - elements_.e * elements_.e))) * tsince;

    // Решаем уравнение Кеплера
    double E = xmdf;
    for(int i = 0; i < 10; i++) {
        double deltaE = (E - elements_.e * sin(E) - xmdf) / (1.0 - elements_.e * cos(E));
        E -= deltaE;
        if(fabs(deltaE) < 1e-12) break;
    }

    // Вычисляем истинную аномалию
    double sinE = sin(E);
    double cosE = cos(E);

    double f = atan2(sqrt(1.0 - elements_.e * elements_.e) * sinE,
                     cosE - elements_.e);

    // Вычисляем расстояние до спутника
    double r = elements_.a * (1.0 - elements_.e * cosE);

    // Ориентация в пространстве
    double cosf = cos(f);
    double sinf = sin(f);

    double cosu = cosf * cos(omgadf) - sinf * sin(omgadf);
    double sinu = sinf * cos(omgadf) + cosf * sin(omgadf);

    // Позиция в плоскости орбиты
    double vx = elements_.a * (-sinE * XKE * sqrt(elements_.a));
    double vy = elements_.a * (sqrt(1.0 - elements_.e * elements_.e) * cosE * XKE * sqrt(elements_.a));

    // Преобразование в ECI
    double sin_xnode = sin(xnode);
    double cos_xnode = cos(xnode);

    // Вычисляем координаты
    pos.setX(r * (cos_xnode * cosu - sin_xnode * sinu * elements_.cosio) * XKMPER);
    pos.setY(r * (sin_xnode * cosu + cos_xnode * sinu * elements_.cosio) * XKMPER);
    pos.setZ(r * sinu * elements_.sinio * XKMPER);

    // Вычисляем скорости
    vel.setX((vx * cos_xnode - vy * sin_xnode * elements_.cosio) * XKMPER / 60.0);
    vel.setY((vx * sin_xnode + vy * cos_xnode * elements_.cosio) * XKMPER / 60.0);
    vel.setZ(vy * elements_.sinio * XKMPER / 60.0);
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = elements_.epoch.msecsTo(time) / (1000.0 * 60.0);

    calculateStateVectors(tsince, state.position, state.velocity);

    return state;
}
