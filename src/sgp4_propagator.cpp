#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initParameters(tle);
    calculateDerivatives();
}

void SGP4Propagator::initParameters(const TLEParser::TLEData& tle) {
    epoch_ = tle.epoch;

    // Преобразование элементов из TLE
    elements_.inclo = tle.inclination * de2ra;
    elements_.nodeo = tle.right_ascension * de2ra;
    elements_.ecco = tle.eccentricity;
    elements_.argpo = tle.argument_perigee * de2ra;
    elements_.mo = tle.mean_anomaly * de2ra;
    elements_.no = tle.mean_motion * 2.0 * M_PI / min_per_day;
    elements_.bstar = tle.bstar;

    // Предварительные вычисления
    double theta2 = cos(elements_.inclo) * cos(elements_.inclo);
    double x3thm1 = 3.0 * theta2 - 1.0;
    double eosq = elements_.ecco * elements_.ecco;
    double betao2 = 1.0 - eosq;
    double betao = sqrt(betao2);

    // Восстановление оригинальной средней движения
    double a1 = pow(ke / elements_.no, 2.0/3.0);
    elements_.del1 = 1.5 * k2 * x3thm1 / (a1 * a1 * betao * betao2);
    elements_.ao = a1 * (1.0 - elements_.del1 * (0.5 * 2.0/3.0 + elements_.del1 * (1.0 + 134.0/81.0 * elements_.del1)));
    elements_.delo = 1.5 * k2 * x3thm1 / (elements_.ao * elements_.ao * betao * betao2);
    elements_.xnodp = elements_.no / (1.0 + elements_.delo);
    elements_.aodp = elements_.ao / (1.0 - elements_.delo);

    // Вычисление орбитальных периодов
    elements_.a = elements_.aodp * xkmper;
    elements_.alta = elements_.a * (1.0 + elements_.ecco) - xkmper;
    elements_.altp = elements_.a * (1.0 - elements_.ecco) - xkmper;
}

void SGP4Propagator::calculateDerivatives() {
    // Вычисление производных элементов
    double cosio = cos(elements_.inclo);
    double sinio = sin(elements_.inclo);
    double cosio2 = cosio * cosio;
    double x3thm1 = 3.0 * cosio2 - 1.0;
    double x1mth2 = 1.0 - cosio2;
    double x7thm1 = 7.0 * cosio2 - 1.0;
    double betao2 = 1.0 - elements_.ecco * elements_.ecco;
    double betao = sqrt(betao2);

    // Вековые возмущения
    elements_.xincl = elements_.inclo;
    double temp1 = 3.0 * k2 * x3thm1 / (4.0 * elements_.aodp * elements_.aodp * betao * betao2);
    double temp2 = temp1 / elements_.aodp;
    elements_.ndot = -temp1 * elements_.xnodp * (1.0 + 1.5 * elements_.ecco);
    elements_.nddot = -temp2 * elements_.xnodp * (1.0 + 4.0 * elements_.ecco);
}

void SGP4Propagator::propagate(double tsince, QVector3D& pos, QVector3D& vel) const {
    // Вековые возмущения
    double betao2 = 1.0 - elements_.ecco * elements_.ecco;
    double betao = sqrt(betao2);
    double xmp = elements_.mo + elements_.xnodp * tsince;
    double omega = elements_.argpo + 0.75 * k2 * tsince * elements_.xnodp * cos(elements_.inclo) /
                                         (elements_.aodp * elements_.aodp * betao * betao2);
    double xnode = elements_.nodeo + elements_.xnodp * tsince *
                                         (-1.5 * k2 * cos(elements_.inclo) / (elements_.aodp * elements_.aodp * betao * betao2));
    double e = elements_.ecco - elements_.bstar * tsince;

    // Решение уравнения Кеплера
    double ep = e;
    double axn = e * cos(omega);
    double ayn = e * sin(omega);
    double xl = xmp + omega;

    double u = xl;
    for(int i = 0; i < 10; i++) {
        double sin_u = sin(u);
        double cos_u = cos(u);
        u = xl + e * sin_u;
        if(fabs(u - xl) < 1.0e-12) break;
        xl = u;
    }

    // Позиция в плоскости орбиты
    double sin_u = sin(u);
    double cos_u = cos(u);
    double r = elements_.aodp * (1.0 - ep * cos_u);
    double rdot = elements_.xnodp * elements_.aodp * ep * sin_u / sqrt(1.0 - ep * ep);
    double rfdot = elements_.xnodp * elements_.aodp * sqrt(1.0 - ep * ep) / r;

    // Ориентация в пространстве
    double sin_i = sin(elements_.xincl);
    double cos_i = cos(elements_.xincl);
    double sin_node = sin(xnode);
    double cos_node = cos(xnode);

    // Позиция (в км)
    pos.setX(xkmper * r * (cos_u * cos_node - sin_u * sin_node * cos_i));
    pos.setY(xkmper * r * (cos_u * sin_node + sin_u * cos_node * cos_i));
    pos.setZ(xkmper * r * sin_u * sin_i);

    // Скорость (в км/с)
    double temp = rdot * cos_u + r * sin_u * rfdot;
    vel.setX((xkmper/60.0) * ((-sin_node * temp) + (cos_node * ((rdot * sin_u - r * cos_u * rfdot) * cos_i))));
    vel.setY((xkmper/60.0) * ((cos_node * temp) + (sin_node * ((rdot * sin_u - r * cos_u * rfdot) * cos_i))));
    vel.setZ((xkmper/60.0) * ((rdot * sin_u - r * cos_u * rfdot) * sin_i));
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = epoch_.msecsTo(time) / (1000.0 * 60.0);

    // Вычисление позиции и скорости
    propagate(tsince, state.position, state.velocity);

    return state;
}
