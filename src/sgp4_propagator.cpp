#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initParameters(tle);
}

void SGP4Propagator::initParameters(const TLEParser::TLEData& tle) {
    epoch_ = tle.epoch;

    // Конвертация основных элементов
    elements_.inclo = tle.inclination * M_PI / 180.0;
    elements_.nodeo = tle.right_ascension * M_PI / 180.0;
    elements_.ecco = tle.eccentricity;
    elements_.argpo = tle.argument_perigee * M_PI / 180.0;
    elements_.mo = tle.mean_anomaly * M_PI / 180.0;
    elements_.no = tle.mean_motion * 2.0 * M_PI / minute_per_day;
    elements_.bstar = tle.bstar;

    // Вычисление вспомогательных величин
    double cosio = cos(elements_.inclo);
    double theta2 = cosio * cosio;
    elements_.x1mth2 = 1.0 - theta2;
    elements_.x7thm1 = 7.0 * theta2 - 1.0;
    double sinio = sin(elements_.inclo);

    // Вычисление большой полуоси
    double a1 = pow(xke / elements_.no, tothrd);
    double betao = sqrt(1.0 - elements_.ecco * elements_.ecco);
    double betao2 = betao * betao;
    double temp0 = 1.5 * ck2 * elements_.x1mth2;
    double del1 = temp0 / (a1 * a1 * betao * betao2);
    double ao = a1 * (1.0 - del1 * (0.5 * tothrd + del1 * (1.0 + 134.0/81.0 * del1)));
    double delo = temp0 / (ao * ao * betao * betao2);

    elements_.a0 = ao;
    elements_.alta = ao * (1.0 + elements_.ecco) - 1.0;
    elements_.altp = ao * (1.0 - elements_.ecco) - 1.0;

    // Вековые возмущения
    elements_.mdot = elements_.no * (1.0 + 1.5 * ck2 * elements_.x1mth2 / (ao * ao * betao * betao2));
    elements_.argpdot = -0.5 * temp0 * elements_.no / (ao * betao * betao2);
    elements_.nodedot = -elements_.argpdot * cosio;

    // Долгопериодические возмущения
    elements_.xlcof = 0.125 * j3 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
    elements_.aycof = 0.25 * j3 * sinio;

    // Коэффициенты для расчета возмущений
    elements_.xmcof = -tothrd * elements_.bstar * xke / elements_.ecco;
    temp0 = 1.0 / (ao - s);
    elements_.t2cof = 1.5 * temp0;
    elements_.t3cof = temp0 * temp0 * temp0;
    elements_.t4cof = elements_.t2cof * temp0;
    elements_.t5cof = elements_.t3cof * temp0;

    elements_.sinmao = sin(elements_.mo);
    elements_.eta = elements_.ecco * elements_.sinmao;
}

void SGP4Propagator::calculateSecularEffects(double tsince, double& xll, double& omgasm,
                                             double& xnodes, double& em, double& xinc,
                                             double& xn) const {
    // Вековые возмущения
    xll = elements_.mo + elements_.mdot * tsince;
    omgasm = elements_.argpo + elements_.argpdot * tsince;
    xnodes = elements_.nodeo + elements_.nodedot * tsince;
    em = elements_.ecco + elements_.xmcof * tsince;
    xinc = elements_.inclo;
    xn = elements_.no;

    // Атмосферное сопротивление
    double temp = 1.0 - elements_.t2cof * tsince;
    temp *= temp;
    xn = elements_.no / temp;
}

void SGP4Propagator::calculatePeriodicEffects(double tsince, double& em, double& xinc,
                                              double& omgasm, double& xnodes, double& xll) const {
    // Долгопериодические возмущения
    double axn = em * cos(omgasm);
    double temp = 1.0 / (em * (1.0 - em * em));
    double xlcof = elements_.xlcof * temp;
    double aycof = elements_.aycof;

    xll += xlcof * axn;
    omgasm -= aycof * (xll + omgasm - xnodes);
}

void SGP4Propagator::solveKepler(double& xll, double e) const {
    for(int i = 0; i < 10; i++) {
        double xl = xll;
        double temp = 1.0 - e * cos(xl);
        xll = xl + (xll - e * sin(xl) - xl) / temp;
        if(fabs(xll - xl) < 1.0e-12) break;
    }
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Время с эпохи в минутах
    double tsince = epoch_.msecsTo(time) / (1000.0 * 60.0);

    // Вычисление возмущенных элементов
    double xll, omgasm, xnodes, em, xinc, xn;
    calculateSecularEffects(tsince, xll, omgasm, xnodes, em, xinc, xn);
    calculatePeriodicEffects(tsince, em, xinc, omgasm, xnodes, xll);

    // Решение уравнения Кеплера
    solveKepler(xll, em);

    // Вычисление позиции в орбитальной плоскости
    double sin_xll = sin(xll);
    double cos_xll = cos(xll);
    double ecosE = 1.0 - em * cos_xll;
    double r = elements_.a0 * (1.0 - em * cos_xll);

    // Вычисление скоростей
    double rdot = xke * sqrt(elements_.a0) * em * sin_xll / r;
    double rfdot = xke * sqrt(elements_.a0 * (1.0 - em * em)) / r;

    // Ориентация в пространстве
    double sin_u = sin(omgasm);
    double cos_u = cos(omgasm);
    double sin_Omega = sin(xnodes);
    double cos_Omega = cos(xnodes);
    double sin_i = sin(xinc);
    double cos_i = cos(xinc);

    // Позиция в ECI (в км)
    double x = r * (cos_u * cos_Omega - sin_u * sin_Omega * cos_i);
    double y = r * (cos_u * sin_Omega + sin_u * cos_Omega * cos_i);
    double z = r * sin_u * sin_i;

    state.position = QVector3D(x, y, z) * xkmper;

    // Скорость в ECI (в км/с)
    double vx = (rdot * (cos_u * cos_Omega - sin_u * sin_Omega * cos_i) -
                 r * (sin_u * cos_Omega + cos_u * sin_Omega * cos_i) * rfdot) * xkmper / 60.0;

    double vy = (rdot * (cos_u * sin_Omega + sin_u * cos_Omega * cos_i) +
                 r * (-sin_u * sin_Omega + cos_u * cos_Omega * cos_i) * rfdot) * xkmper / 60.0;

    double vz = (rdot * sin_u * sin_i + r * cos_u * sin_i * rfdot) * xkmper / 60.0;

    state.velocity = QVector3D(vx, vy, vz);

    return state;
}
