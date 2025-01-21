#include "sgp4_propagator.h"
#include <QtMath>

SGP4Propagator::SGP4Propagator(const TLEParser::TLEData& tle) {
    initializeParameters(tle);
    sgp4init();
}

void SGP4Propagator::initializeParameters(const TLEParser::TLEData& tle) {
    epoch_ = tle.epoch;

    // Преобразование базовых элементов из TLE
    elements_.inclo = tle.inclination * DE2RA;
    elements_.nodeo = tle.right_ascension * DE2RA;
    elements_.ecco = tle.eccentricity;
    elements_.argpo = tle.argument_perigee * DE2RA;
    elements_.mo = tle.mean_anomaly * DE2RA;
    elements_.no = tle.mean_motion * 2.0 * M_PI / MINUTES_PER_DAY;
    elements_.bstar = tle.bstar;

    // Расчет юлианской даты эпохи
    QDate baseDate(1949, 12, 31);
    double days = baseDate.daysTo(tle.epoch.date());
    double fractionalDay = (tle.epoch.time().msecsSinceStartOfDay() / 1000.0) / 86400.0;
    elements_.jdsatepoch = 2433281.5 + days + fractionalDay;
}

void SGP4Propagator::sgp4init() {
    // Вычисление вспомогательных параметров
    const double cosio = cos(elements_.inclo);
    const double theta2 = cosio * cosio;
    const double x3thm1 = 3.0 * theta2 - 1.0;
    const double eosq = elements_.ecco * elements_.ecco;
    const double betao2 = 1.0 - eosq;
    const double betao = sqrt(betao2);

    // Вычисление большой полуоси
    double ao, sinio, po, con42, con41, ainv, einv, posq;

    // Восстановление исходной средней движения
    double temp = elements_.no / XKE;
    ao = pow(temp * temp, -1.0 / 3.0);
    sinio = sin(elements_.inclo);
    po = ao * betao2;
    con42 = 1.0 - 5.0 * theta2;
    con41 = -con42 - theta2 - theta2;
    ainv = 1.0 / ao;
    einv = 1.0 / (elements_.ecco - 1.0);
    posq = po * po;

    // Для J2 возмущений
    double temp1 = 3.0 * CK2 * ainv * ainv;
    double temp2 = temp1 * CK2 * ainv;
    double temp3 = 1.25 * CK4 * ainv * ainv * ainv * ainv;

    elements_.a = ao / (1.0 - temp1 * betao * x3thm1 / 2.0 -
                        temp2 * betao * (13.0 - 78.0 * theta2 + 137.0 * theta2 * theta2) / 16.0);

    double delta1 = temp1 * x3thm1 / 2.0 +
                    temp2 * (13.0 - 78.0 * theta2 + 137.0 * theta2 * theta2) / 16.0;

    elements_.no = elements_.no / (1.0 + delta1);

    // Рассчитываем апогей и перигей
    elements_.alta = elements_.a * (1.0 + elements_.ecco) - 1.0;
    elements_.altp = elements_.a * (1.0 - elements_.ecco) - 1.0;

    // Преобразуем в километры
    elements_.alta *= XKMPER;
    elements_.altp *= XKMPER;
}

void SGP4Propagator::sgp4(double tsince, QVector3D& pos, QVector3D& vel) const {
    // Вычисляем возмущенные элементы орбиты
    double xmdf = elements_.mo + elements_.no * tsince;
    double omega = elements_.argpo;
    double xnode = elements_.nodeo;
    double bstar = elements_.bstar;

    // Учет J2 возмущений
    const double cosio = cos(elements_.inclo);
    const double theta2 = cosio * cosio;
    const double x3thm1 = 3.0 * theta2 - 1.0;

    double omgdot = -1.5 * CK2 * elements_.no * x3thm1 / (elements_.a * elements_.a);
    double xnodot = -1.5 * CK2 * elements_.no * cosio / (elements_.a * elements_.a);

    omega += omgdot * tsince;
    xnode += xnodot * tsince;

    // Решение уравнения Кеплера
    double xl = xmdf;
    double u = xl - xnode;
    double e = elements_.ecco;

    double iter = 0;
    double sineo, coseo;
    double eo;
    do {
        eo = xl;
        sineo = sin(xl);
        coseo = cos(xl);
        xl = u + e * sineo;
        iter++;
    } while (fabs(xl - eo) >= 1.0e-12 && iter < 10);

    // Вычисление позиции и скорости в орбитальной плоскости
    double sinuk = sin(u);
    double cosuk = cos(u);

    double r = elements_.a * (1.0 - e * coseo);
    double rdot = sqrt(elements_.a) * e * sineo / r;
    double rfdot = sqrt(elements_.a * (1.0 - e * e)) * coseo / r;

    // Преобразование в ECI координаты
    double sini = sin(elements_.inclo);
    double cosi = cos(elements_.inclo);
    double sinnok = sin(xnode);
    double cosnok = cos(xnode);
    double sinomega = sin(omega);
    double cosomega = cos(omega);

    // Позиция (в км)
    double xmx = r * ((cosuk * cosnok - sinuk * sinnok * cosi) * cosomega -
                      (sinuk * cosnok + cosuk * sinnok * cosi) * sinomega);
    double xmy = r * ((cosuk * sinnok + sinuk * cosnok * cosi) * cosomega -
                      (sinuk * sinnok - cosuk * cosnok * cosi) * sinomega);
    double xmz = r * (sinuk * sini + cosuk * cosi * sinomega);

    pos = QVector3D(xmx, xmy, xmz) * XKMPER;

    // Скорость (в км/с)
    double vx = XKMPER * ((rdot * (cosuk * cosnok - sinuk * sinnok * cosi) -
                           rfdot * (sinuk * cosnok + cosuk * sinnok * cosi)) * cosomega -
                          (rdot * (sinuk * cosnok + cosuk * sinnok * cosi) +
                           rfdot * (cosuk * cosnok - sinuk * sinnok * cosi)) * sinomega) / 60.0;

    double vy = XKMPER * ((rdot * (cosuk * sinnok + sinuk * cosnok * cosi) -
                           rfdot * (sinuk * sinnok - cosuk * cosnok * cosi)) * cosomega -
                          (rdot * (sinuk * sinnok - cosuk * cosnok * cosi) +
                           rfdot * (cosuk * sinnok + sinuk * cosnok * cosi)) * sinomega) / 60.0;

    double vz = XKMPER * (rdot * sinuk * sini + rfdot * cosuk * sini +
                          (rdot * cosuk * cosi - rfdot * sinuk * cosi) * sinomega) / 60.0;

    vel = QVector3D(vx, vy, vz);
}

SGP4Propagator::OrbitalState SGP4Propagator::calculateState(const QDateTime& time) const {
    OrbitalState state;
    state.epoch = time;

    // Расчет времени с эпохи в минутах
    double tsince = epoch_.msecsTo(time) / (1000.0 * 60.0);

    // Вычисление позиции и скорости
    sgp4(tsince, state.position, state.velocity);

    return state;
}
