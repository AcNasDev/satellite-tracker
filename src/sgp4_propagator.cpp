#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initParameters(tle);
}

void SGP4Propagator::initParameters(const TLEParser::TLEData& tle) {
    epoch_ = tle.epoch;

    qDebug() << "\n=== SGP4 Initialization ===";

    // Инициализация базовых элементов
    elements_.inclo = tle.inclination * de2ra;
    elements_.nodeo = tle.right_ascension * de2ra;
    elements_.ecco = tle.eccentricity;
    elements_.argpo = tle.argument_perigee * de2ra;
    elements_.mo = tle.mean_anomaly * de2ra;
    elements_.no = tle.mean_motion * 2.0 * M_PI / min_per_day;
    elements_.bstar = tle.bstar;

    qDebug() << "Initial elements:";
    qDebug() << "Inclination (rad):" << elements_.inclo;
    qDebug() << "RAAN (rad):" << elements_.nodeo;
    qDebug() << "Eccentricity:" << elements_.ecco;
    qDebug() << "Arg of perigee (rad):" << elements_.argpo;
    qDebug() << "Mean anomaly (rad):" << elements_.mo;
    qDebug() << "Mean motion (rad/min):" << elements_.no;
    qDebug() << "BSTAR (1/ER):" << elements_.bstar;

    // Вспомогательные параметры
    double cosio = cos(elements_.inclo);
    double theta2 = cosio * cosio;
    double x3thm1 = 3.0 * theta2 - 1.0;
    double eosq = elements_.ecco * elements_.ecco;
    double betao2 = 1.0 - eosq;
    double betao = sqrt(betao2);
    double del1 = 1.5 * ck2 * x3thm1 / (betao * betao2);

    // Восстановление оригинального среднего движения
    elements_.aodp = pow(xke/elements_.no, 2.0/3.0);
    elements_.a = elements_.aodp * (1.0 - del1 * ((1.0/3.0) + del1 * (1.0 + 134.0/81.0 * del1)));
    double delo = 1.5 * ck2 * x3thm1 / (elements_.a * betao * betao2);
    elements_.no_kozai = elements_.no / (1.0 + delo);

    // Дополнительные коэффициенты
    double sinio = sin(elements_.inclo);
    elements_.aycof = -0.5 * xj3 / ck2 * sinio;
    double x7thm1 = 7.0 * theta2 - 1.0;
    if (fabs(cosio + 1.0) > 1.5e-12) {
        elements_.xlcof = 0.125 * xj3 / ck2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
    } else {
        elements_.xlcof = 0.125 * xj3 / ck2 * sinio * (3.0 + 5.0 * cosio) / 1.5e-12;
    }

    // Производные
    double temp = 1.5 * ck2 * x3thm1 / (betao * betao2);
    elements_.ndot = -temp * elements_.no_kozai * (1.0 + 4.0 * betao + eosq);
    elements_.nddot = -temp * elements_.ndot;

    // Апогей и перигей
    elements_.alta = elements_.a * (1.0 + elements_.ecco) - 1.0;
    elements_.altp = elements_.a * (1.0 - elements_.ecco) - 1.0;

    qDebug() << "\nDerived elements:";
    qDebug() << "Semi-major axis (ER):" << elements_.a;
    qDebug() << "Mean motion (Kozai, rad/min):" << elements_.no_kozai;
    qDebug() << "ndot (rad/min²):" << elements_.ndot;
    qDebug() << "nddot (rad/min³):" << elements_.nddot;
    qDebug() << "Apogee height (ER):" << elements_.alta;
    qDebug() << "Perigee height (ER):" << elements_.altp;
}

void SGP4Propagator::propagate(double tsince, QVector3D& pos, QVector3D& vel)  {
    qDebug() << "\nPropagating for tsince =" << tsince << "minutes";

    // Обновление среднего движения с учетом производных
    double xn = elements_.no_kozai;
    double tempa = 1.0 - elements_.bstar * tsince * elements_.ndot / (2.0 * xn);
    double temp = pow(tempa, 3.5);
    xn = xn * temp;

    // Вековые возмущения
    double xll = elements_.mo + xn * tsince;
    double omega = elements_.argpo;
    double xnode = elements_.nodeo + elements_.ndot * tsince * tsince * 1.5;
    double e = elements_.ecco - elements_.bstar * tsince;
    e = qMax(1.0e-6, e);

    qDebug() << "Secular perturbations:";
    qDebug() << "Mean motion (rad/min):" << xn;
    qDebug() << "Eccentricity:" << e;
    qDebug() << "Mean anomaly (rad):" << xll;
    qDebug() << "Arg of perigee (rad):" << omega;
    qDebug() << "RAAN (rad):" << xnode;

    // Учет долгопериодических возмущений
    double axn = e * cos(omega);
    double temp1 = 1.0 / (elements_.a * (1.0 - e * e));
    double xlcof = 0.125 * xj3 / ck2 * elements_.inclo * (3.0 + 5.0 * cos(elements_.inclo)) /
                   (1.0 + cos(elements_.inclo));
    double aynl = e * sin(omega) + temp1 * elements_.aycof;
    double xlt = xll + temp1 * xlcof * axn;

    double u = solveKepler(xlt - xnode, e);
    double sin_u = sin(u);
    double cos_u = cos(u);
    double r = elements_.a * (1.0 - e * cos_u);
    double rdot = xke * sqrt(elements_.a) * e * sin_u / r;
    double rfdot = xke * sqrt(elements_.a * (1.0 - e * e)) / r;

    qDebug() << "Orbital plane values:";
    qDebug() << "r (ER):" << r;
    qDebug() << "rdot (ER/min):" << rdot;
    qDebug() << "rfdot (rad/min):" << rfdot;

    // Расчет координат в ECI
    double sin_i = sin(elements_.inclo);
    double cos_i = cos(elements_.inclo);
    double sin_node = sin(xnode);
    double cos_node = cos(xnode);
    double sin_u_omega = sin(u + omega);
    double cos_u_omega = cos(u + omega);

    pos.setX(r * (cos_u_omega * cos_node - sin_u_omega * sin_node * cos_i) * xkmper);
    pos.setY(r * (cos_u_omega * sin_node + sin_u_omega * cos_node * cos_i) * xkmper);
    pos.setZ(r * sin_u_omega * sin_i * xkmper);

    vel.setX(((rdot * cos_u_omega - r * rfdot * sin_u_omega) * cos_node -
              (rdot * sin_u_omega + r * rfdot * cos_u_omega) * sin_node * cos_i) * xkmper / 60.0);
    vel.setY(((rdot * cos_u_omega - r * rfdot * sin_u_omega) * sin_node +
              (rdot * sin_u_omega + r * rfdot * cos_u_omega) * cos_node * cos_i) * xkmper / 60.0);
    vel.setZ((rdot * sin_u_omega + r * rfdot * cos_u_omega) * sin_i * xkmper / 60.0);

    qDebug() << "\nFinal ECI coordinates:";
    qDebug() << "Position (km):" << pos;
    qDebug() << "Velocity (km/s):" << vel;
}

double SGP4Propagator::solveKepler(double M, double e) {
    double E = M;
    for(int i = 0; i < 10; i++) {
        double E_old = E;
        E = M + e * sin(E);
        if(fabs(E - E_old) < 1.0e-12) break;
    }
    return E;
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) {
    OrbitalState state;
    state.epoch = time;
    double tsince = epoch_.msecsTo(time) / (1000.0 * 60.0);
    propagate(tsince, state.position, state.velocity);
    return state;
}

