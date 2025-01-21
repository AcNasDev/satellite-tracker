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

    // Конвертация начальных элементов в каноническую систему единиц
    elements_.inclo = tle.inclination * de2ra;
    elements_.nodeo = tle.right_ascension * de2ra;
    elements_.ecco = tle.eccentricity;
    elements_.argpo = tle.argument_perigee * de2ra;
    elements_.mo = tle.mean_anomaly * de2ra;
    elements_.no = tle.mean_motion * 2.0 * M_PI / min_per_day;
    elements_.bstar = tle.bstar / ae; // Важно: масштабирование BSTAR

    qDebug() << "\nInitial orbital elements:";
    qDebug() << "Inclination (rad):" << elements_.inclo;
    qDebug() << "RAAN (rad):" << elements_.nodeo;
    qDebug() << "Eccentricity:" << elements_.ecco;
    qDebug() << "Arg of perigee (rad):" << elements_.argpo;
    qDebug() << "Mean anomaly (rad):" << elements_.mo;
    qDebug() << "Mean motion (rad/min):" << elements_.no;
    qDebug() << "BSTAR:" << elements_.bstar;

    // Расчет вспомогательных величин
    double cosio = cos(elements_.inclo);
    double cosio2 = cosio * cosio;
    double betao2 = 1.0 - elements_.ecco * elements_.ecco;
    double betao = sqrt(betao2);

    // Восстановление оригинального среднего движения
    double ak = pow(ke / elements_.no, 2.0/3.0);
    double d1 = 0.75 * k2 * (3.0 * cosio2 - 1.0) / (betao2 * betao);
    double del1 = d1 / (ak * ak);
    double ao = ak * (1.0 - del1 * (1.0/3.0 + del1 * (1.0 + 134.0/81.0 * del1)));

    elements_.a = ao;
    elements_.n = elements_.no;

    qDebug() << "\nDerived elements:";
    qDebug() << "Semi-major axis (ER):" << elements_.a;
    qDebug() << "Corrected mean motion (rad/min):" << elements_.n;
}

void SGP4Propagator::propagate(double tsince, QVector3D& pos, QVector3D& vel) const {
    qDebug() << "\nPropagation at time:" << tsince << "minutes";

    // Обновление средних элементов с учетом вековых возмущений
    double xndt = elements_.n * 1.5 * k2 * (3.0 * cos(elements_.inclo) * cos(elements_.inclo) - 1.0) /
                  (pow(1.0 - elements_.ecco * elements_.ecco, 1.5));
    double xnddt = elements_.n * xndt * k2 * cos(elements_.inclo) /
                   (pow(1.0 - elements_.ecco * elements_.ecco, 2.0));

    // Обновление среднего движения с учетом торможения
    double xn = elements_.n + xndt * tsince + xnddt * tsince * tsince * 0.5;

    // Обновление эксцентриситета с учетом атмосферного торможения
    double e = elements_.ecco - elements_.bstar * tsince;
    e = std::max(0.0, std::min(0.999999, e)); // Ограничение эксцентриситета

    qDebug() << "Updated elements:";
    qDebug() << "Mean motion (rad/min):" << xn;
    qDebug() << "Eccentricity:" << e;

    // Обновление средней аномалии
    double xmp = elements_.mo + xn * tsince;

    // Вековые возмущения для аргумента перигея и RAAN
    double omega = elements_.argpo + 0.75 * k2 * sin(elements_.inclo) * sin(elements_.inclo) *
                                         xn * tsince / (elements_.n * (1.0 - e * e));
    double xnode = elements_.nodeo - 1.5 * k2 * cos(elements_.inclo) * xn * tsince /
                                         (elements_.n * (1.0 - e * e));

    qDebug() << "Secular perturbations:";
    qDebug() << "Mean anomaly (rad):" << xmp;
    qDebug() << "Arg of perigee (rad):" << omega;
    qDebug() << "RAAN (rad):" << xnode;

    // Решение уравнения Кеплера
    double xl = xmp;
    double u = xl - xnode;
    double converge = false;
    for(int i = 0; i < 10 && !converge; i++) {
        double sinu = sin(u);
        double cosu = cos(u);
        double old_u = u;
        u = u + (xl - xnode - u + e * sinu) / (1.0 - e * cosu);
        converge = fabs(u - old_u) < 1e-12;
    }

    // Вычисление позиции в орбитальной плоскости
    double r = elements_.a * (1.0 - e * cos(u));
    double rdot = ke * sqrt(elements_.a) * e * sin(u) / r;
    double rfdot = ke * sqrt(elements_.a * (1.0 - e * e)) / r;

    qDebug() << "Orbital plane values:";
    qDebug() << "r (ER):" << r;
    qDebug() << "rdot (ER/min):" << rdot;
    qDebug() << "rfdot (rad/min):" << rfdot;

    // Преобразование в ECI
    double cosu = cos(u);
    double sinu = sin(u);
    double cosw = cos(omega);
    double sinw = sin(omega);
    double cosi = cos(elements_.inclo);
    double sini = sin(elements_.inclo);
    double cosnode = cos(xnode);
    double sinnode = sin(xnode);

    // Позиция (в км)
    double x = r * ((cosnode * cosw - sinnode * sinw * cosi) * cosu -
                    (cosnode * sinw + sinnode * cosw * cosi) * sinu);
    double y = r * ((sinnode * cosw + cosnode * sinw * cosi) * cosu -
                    (sinnode * sinw - cosnode * cosw * cosi) * sinu);
    double z = r * (sinw * cosu + cosw * sinu) * sini;

    pos = QVector3D(x, y, z) * xkmper;

    // Скорость (в км/с)
    double vx = ((-cosnode * cosw + sinnode * sinw * cosi) * sinu * rdot +
                 (-cosnode * sinw - sinnode * cosw * cosi) * cosu * rfdot) / 60.0;
    double vy = ((-sinnode * cosw - cosnode * sinw * cosi) * sinu * rdot +
                 (-sinnode * sinw + cosnode * cosw * cosi) * cosu * rfdot) / 60.0;
    double vz = (sinw * sinu * rdot + cosw * cosu * rfdot) * sini / 60.0;

    vel = QVector3D(vx, vy, vz) * xkmper;

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
