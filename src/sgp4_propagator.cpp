#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    sgp4init(tle);
}

void SGP4Propagator::sgp4init(const TLEParser::TLEData& tle) {
    epoch_ = tle.epoch;

    // Инициализация параметров из TLE
    satrec_.inclo = tle.inclination * pi / 180.0;
    satrec_.nodeo = tle.right_ascension * pi / 180.0;
    satrec_.ecco = tle.eccentricity;
    satrec_.argpo = tle.argument_perigee * pi / 180.0;
    satrec_.mo = tle.mean_anomaly * pi / 180.0;
    satrec_.no_kozai = tle.mean_motion * 2.0 * pi / minutes_per_day;
    satrec_.bstar = tle.bstar;

    // Инициализация переменных
    satrec_.isimp = 0;
    satrec_.method = 'n';
    satrec_.aycof = 0.0;
    satrec_.con41 = 0.0;
    satrec_.cc1 = 0.0;
    satrec_.cc4 = 0.0;
    satrec_.cc5 = 0.0;
    satrec_.d2 = 0.0;
    satrec_.d3 = 0.0;
    satrec_.d4 = 0.0;
    satrec_.delmo = 0.0;
    satrec_.eta = 0.0;
    satrec_.argpdot = 0.0;
    satrec_.omgcof = 0.0;
    satrec_.sinmao = 0.0;
    satrec_.t2cof = 0.0;
    satrec_.t3cof = 0.0;
    satrec_.t4cof = 0.0;
    satrec_.t5cof = 0.0;

    // Вычисление вспомогательных величин
    double cosio = cos(satrec_.inclo);
    double sinio = sin(satrec_.inclo);
    double cosio2 = cosio * cosio;

    satrec_.x3thm1 = 3.0 * cosio2 - 1.0;
    satrec_.x1mth2 = 1.0 - cosio2;
    satrec_.x7thm1 = 7.0 * cosio2 - 1.0;

    // Вычисление большой полуоси
    double a1 = pow(xke / satrec_.no_kozai, tothrd);
    satrec_.rcse = 1.0 / (1.0 - satrec_.ecco);

    // Коррекция для SGP4
    double delta1 = 1.5 * ck2 * satrec_.x3thm1 / (a1 * a1);
    double a0 = a1 * (1.0 - delta1 * (0.5 * tothrd + delta1 * (1.0 + 134.0/81.0 * delta1)));
    double delta0 = 1.5 * ck2 * satrec_.x3thm1 / (a0 * a0);

    satrec_.a = a0 / (1.0 - delta0);
    satrec_.alta = satrec_.a * (1.0 + satrec_.ecco) - 1.0;
    satrec_.altp = satrec_.a * (1.0 - satrec_.ecco) - 1.0;

    // Вычисление производных величин
    double temp1 = 1.5 * ck2 * satrec_.x3thm1;
    double temp2 = 0.5 * temp1;

    satrec_.mdot = satrec_.no_kozai + 0.5 * temp1 * satrec_.rcse / sqrt(satrec_.a);
    satrec_.argpdot = -temp2 * cosio / satrec_.a;
    satrec_.nodedot = -temp2 * cosio / satrec_.a;

    // Дополнительные коэффициенты
    satrec_.xlcof = 0.125 * a3ovk2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
    satrec_.aycof = 0.25 * a3ovk2 * sinio;
}

QVector3D SGP4Propagator::getPosition(double tsince) const {
    // Вычисление средней аномалии
    double xmdf = satrec_.mo + satrec_.mdot * tsince;
    double omgadf = satrec_.argpo + satrec_.argpdot * tsince;
    double xnode = satrec_.nodeo + satrec_.nodedot * tsince;

    // Решение уравнения Кеплера
    double e = satrec_.ecco;
    double xl = xmdf;
    for(int i = 0; i < 10; i++) {
        double xll = xl;
        xl = xmdf + e * sin(xl);
        if(fabs(xl - xll) < 1.0e-12) break;
    }

    // Вычисление позиции в орбитальной плоскости
    double u = xl - xnode;
    double sin_u = sin(u);
    double cos_u = cos(u);

    double r = satrec_.a * (1.0 - e * cos(xl));
    double xmx = r * (cos_u * cos(xnode) - sin_u * sin(xnode) * cos(satrec_.inclo));
    double xmy = r * (cos_u * sin(xnode) + sin_u * cos(xnode) * cos(satrec_.inclo));
    double xmz = r * sin_u * sin(satrec_.inclo);

    return QVector3D(xmx, xmy, xmz) * xkmper;
}

QVector3D SGP4Propagator::getVelocity(double tsince) const {
    const double dt = 0.0001;
    QVector3D r1 = getPosition(tsince - dt);
    QVector3D r2 = getPosition(tsince + dt);
    return (r2 - r1) / (2.0 * dt * 60.0); // переводим в км/с
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    double tsince = epoch_.msecsTo(time) / (1000.0 * 60.0);
    state.position = getPosition(tsince);
    state.velocity = getVelocity(tsince);

    return state;
}
