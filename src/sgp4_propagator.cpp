#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initParameters(tle);
    debugElements();
}

void SGP4Propagator::debugElements() const {
    qDebug() << "\n=== SGP4 Debug Information ===";
    qDebug() << "Initial orbital elements:";
    qDebug() << "Inclination (deg):" << elements_.inclo / de2ra;
    qDebug() << "RAAN (deg):" << elements_.nodeo / de2ra;
    qDebug() << "Eccentricity:" << elements_.ecco;
    qDebug() << "Arg of Perigee (deg):" << elements_.argpo / de2ra;
    qDebug() << "Mean Anomaly (deg):" << elements_.mo / de2ra;
    qDebug() << "Mean Motion (rev/day):" << elements_.no * min_per_day / (2.0 * M_PI);
    qDebug() << "BSTAR:" << elements_.bstar;
    qDebug() << "\nDerived elements:";
    qDebug() << "Semi-major axis (km):" << elements_.a;
    qDebug() << "Corrected mean motion (rad/min):" << elements_.n;
    qDebug() << "Corrected eccentricity:" << elements_.e;
    qDebug() << "Corrected inclination (deg):" << elements_.i / de2ra;
}

void SGP4Propagator::initParameters(const TLEParser::TLEData& tle) {
    epoch_ = tle.epoch;

    // Конвертация элементов
    elements_.inclo = tle.inclination * de2ra;
    elements_.nodeo = tle.right_ascension * de2ra;
    elements_.ecco = tle.eccentricity;
    elements_.argpo = tle.argument_perigee * de2ra;
    elements_.mo = tle.mean_anomaly * de2ra;
    elements_.no = tle.mean_motion * 2.0 * M_PI / min_per_day;
    elements_.bstar = tle.bstar;

    // Вычисление производных элементов
    double cosio = cos(elements_.inclo);
    double theta2 = cosio * cosio;
    double x3thm1 = 3.0 * theta2 - 1.0;
    double eosq = elements_.ecco * elements_.ecco;
    double betao2 = 1.0 - eosq;
    double betao = sqrt(betao2);

    // Коррекция средней движения
    double a1 = pow(ke / elements_.no, 2.0/3.0);
    double delta1 = 1.5 * k2 * x3thm1 / (a1 * a1 * betao * betao2);
    double ao = a1 * (1.0 - delta1 * (0.5 * 2.0/3.0 + delta1 * (1.0 + 134.0/81.0 * delta1)));
    double delo = 1.5 * k2 * x3thm1 / (ao * ao * betao * betao2);

    elements_.a = ao * xkmper;
    elements_.n = elements_.no / (1.0 + delo);
    elements_.e = elements_.ecco;
    elements_.i = elements_.inclo;
    elements_.omega = elements_.argpo;
    elements_.Omega = elements_.nodeo;

    qDebug() << "\nInitialization parameters:";
    qDebug() << "a1:" << a1;
    qDebug() << "delta1:" << delta1;
    qDebug() << "ao:" << ao;
    qDebug() << "delo:" << delo;
}

void SGP4Propagator::propagate(double tsince, QVector3D& pos, QVector3D& vel) const {
    qDebug() << "\nPropagation at time:" << tsince << "minutes";

    // Вековые возмущения
    double xmp = elements_.mo + elements_.n * tsince;
    double omega = elements_.omega + 0.75 * k2 * elements_.n * tsince * cos(elements_.i);
    double xnode = elements_.Omega - 1.5 * k2 * elements_.n * tsince * cos(elements_.i);

    // Исправлено: правильный учет BSTAR
    double e = elements_.e;
    if (elements_.bstar != 0.0) {
        e = elements_.e - elements_.bstar * ke * elements_.a * sin(xmp) * tsince;
    }

    qDebug() << "Secular perturbations:";
    qDebug() << "Mean anomaly:" << xmp / de2ra;
    qDebug() << "Argument of perigee:" << omega / de2ra;
    qDebug() << "RAAN:" << xnode / de2ra;
    qDebug() << "Eccentricity:" << e;

    // Решение уравнения Кеплера
    double E = xmp;
    for(int i = 0; i < 10; i++) {
        double E_old = E;
        E = xmp + e * sin(E);
        if(fabs(E - E_old) < 1.0e-12) break;
    }

    // Вычисление истинной аномалии
    double sinE = sin(E);
    double cosE = cos(E);
    double sinv = (sqrt(1.0 - e * e) * sinE) / (1.0 - e * cosE);
    double cosv = (cosE - e) / (1.0 - e * cosE);
    double v = atan2(sinv, cosv);

    // Позиция в плоскости орбиты
    double r = elements_.a * (1.0 - e * cosE);
    double rdot = ke * sqrt(elements_.a) * sinE / sqrt(r);
    double rfdot = ke * sqrt(elements_.a * (1.0 - e * e)) / r;

    qDebug() << "\nOrbital plane values:";
    qDebug() << "r (km):" << r;
    qDebug() << "rdot (km/min):" << rdot;
    qDebug() << "rfdot (rad/min):" << rfdot;

    // Преобразование в ECI
    double sinO = sin(xnode);
    double cosO = cos(xnode);
    double sino = sin(omega + v);
    double coso = cos(omega + v);
    double sini = sin(elements_.i);
    double cosi = cos(elements_.i);

    pos.setX(r * (cosO * coso - sinO * sino * cosi));
    pos.setY(r * (sinO * coso + cosO * sino * cosi));
    pos.setZ(r * sino * sini);

    // Скорости в ECI (переводим из км/мин в км/с)
    double vx = (rdot * (cosO * coso - sinO * sino * cosi) -
                 r * rfdot * (cosO * sino + sinO * coso * cosi)) / 60.0;
    double vy = (rdot * (sinO * coso + cosO * sino * cosi) +
                 r * rfdot * (cosO * coso - sinO * sino * cosi)) / 60.0;
    double vz = (rdot * sino * sini + r * rfdot * coso * sini) / 60.0;

    vel.setX(vx);
    vel.setY(vy);
    vel.setZ(vz);

    qDebug() << "\nFinal ECI coordinates:";
    qDebug() << "Position (km):" << pos;
    qDebug() << "Velocity (km/s):" << vel;
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    double tsince = epoch_.msecsTo(time) / (1000.0 * 60.0);
    qDebug() << "\n=== Calculating state for time since epoch:" << tsince << "minutes ===";

    propagate(tsince, state.position, state.velocity);
    return state;
}
