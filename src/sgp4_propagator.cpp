#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
}

void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    // Сохраняем основные элементы
    elements_.epoch = tle.epoch;
    elements_.i = tle.inclination * DE2RA;
    elements_.Omega = tle.right_ascension * DE2RA;
    elements_.e = tle.eccentricity;
    elements_.omega = tle.argument_perigee * DE2RA;
    elements_.M = tle.mean_anomaly * DE2RA;
    elements_.n0 = tle.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;
    elements_.bstar = tle.bstar;

    // Вычисляем производные параметры
    elements_.cosio = cos(elements_.i);
    elements_.sinio = sin(elements_.i);

    // Вычисляем большую полуось
    const double theta2 = elements_.cosio * elements_.cosio;
    const double x3thm1 = 3.0 * theta2 - 1.0;
    const double beta0 = sqrt(1.0 - elements_.e * elements_.e);

    double a1 = pow(XKE / elements_.n0, 2.0/3.0);
    double delta1 = 1.5 * CK2 * x3thm1 / (a1 * a1 * beta0 * beta0 * beta0);

    elements_.a = a1 * (1.0 - delta1/3.0 - delta1 * delta1 - 134.0/81.0 * delta1 * delta1 * delta1);
    elements_.n = elements_.n0;
}

double SGP4Propagator::solveKeplerEquation(double M, double e, int maxIter) const {
    // Начальное приближение
    double E = M;

    // Метод Ньютона с ограничением итераций
    for(int i = 0; i < maxIter; i++) {
        double delta = (E - e * sin(E) - M) / (1.0 - e * cos(E));
        E -= delta;

        if(fabs(delta) < 1e-12) break;
    }

    return E;
}

QVector3D SGP4Propagator::calculatePosition(double xnode, double xinc,
                                            double xll, double e, double omega,
                                            double xn) const {
    // Решаем уравнение Кеплера
    double xl = xll + omega;
    double E = solveKeplerEquation(xl - xnode, e);

    // Вычисляем истинную аномалию
    double cosE = cos(E);
    double sinE = sin(E);
    double beta = sqrt(1.0 - e * e);
    double f = atan2(beta * sinE, cosE - e);

    // Вычисляем расстояние
    double r = elements_.a * (1.0 - e * cosE);

    // Вычисляем координаты в орбитальной плоскости
    double cos_f = cos(f);
    double sin_f = sin(f);
    double u = omega + f;
    double cos_u = cos(u);
    double sin_u = sin(u);

    // Преобразование в ECI
    double cos_node = cos(xnode);
    double sin_node = sin(xnode);
    double cos_i = cos(xinc);
    double sin_i = sin(xinc);

    double x = r * (cos_node * cos_u - sin_node * sin_u * cos_i);
    double y = r * (sin_node * cos_u + cos_node * sin_u * cos_i);
    double z = r * sin_u * sin_i;

    return QVector3D(x, y, z);
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = elements_.epoch.msecsTo(time) / (1000.0 * 60.0);

    // Вековые возмущения
    double xnodot = -6.0 * CK2 * elements_.n * elements_.cosio /
                    (elements_.a * elements_.a * pow(1.0 - elements_.e * elements_.e, 3.5));

    double omegadot = -3.0 * CK2 * elements_.n * (1.0 - 5.0 * elements_.cosio * elements_.cosio) /
                      (2.0 * elements_.a * elements_.a * pow(1.0 - elements_.e * elements_.e, 2));

    // Обновляем элементы
    double xnode = elements_.Omega + xnodot * tsince;
    double omega = elements_.omega + omegadot * tsince;
    double xn = elements_.n;
    double xl = elements_.M + xn * tsince;

    // Вычисляем позицию
    state.position = calculatePosition(xnode, elements_.i, xl, elements_.e, omega, xn);

    // Вычисляем скорость через конечные разности
    const double dt = 1.0; // 1 секунда
    double tsince2 = tsince + dt/60.0;

    double xnode2 = elements_.Omega + xnodot * tsince2;
    double omega2 = elements_.omega + omegadot * tsince2;
    double xl2 = elements_.M + xn * tsince2;

    QVector3D pos2 = calculatePosition(xnode2, elements_.i, xl2, elements_.e, omega2, xn);
    state.velocity = (pos2 - state.position) / dt;

    return state;
}
