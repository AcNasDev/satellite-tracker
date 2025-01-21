// sgp4_propagator.h
#ifndef SGP4_PROPAGATOR_H
#define SGP4_PROPAGATOR_H

#include "tle_parser.h"
#include <QVector3D>
#include <QDateTime>
#include <QDebug>

class SGP4Propagator {
public:
    struct OrbitalState {
        QVector3D position;
        QVector3D velocity;
        QDateTime epoch;
    };

    explicit SGP4Propagator(const TLEParser::TLEData& tle);
    OrbitalState calculateState(const QDateTime& time);

private:
    // Константы SGP4
    static constexpr double xke = 0.0743669161331734049; // √(GM) ER^(3/2)/min
    static constexpr double xj2 = 1.082616e-3;           // J2
    static constexpr double xj3 = -0.253881e-5;          // J3
    static constexpr double xj4 = -1.65597e-6;           // J4
    static constexpr double xkmper = 6378.137;           // Радиус Земли (км)
    static constexpr double ae = 1.0;                    // Радиус Земли (ER)
    static constexpr double de2ra = M_PI / 180.0;        // Градусы в радианы
    static constexpr double min_per_day = 1440.0;        // Минут в сутках
    static constexpr double ck2 = xj2/2.0;
    static constexpr double ck4 = -3.0*xj4/8.0;
    static constexpr double qoms2t = 1.880279e-09;       // (q0-s)^4 при ae=1
    static constexpr double s = ae + 78.0/xkmper;        // s при ae=1

    struct Elements {
        // Исходные элементы
        double no;       // Исходное среднее движение (рад/мин)
        double ecco;     // Исходный эксцентриситет
        double inclo;    // Исходное наклонение (рад)
        double nodeo;    // Исходная долгота восходящего узла (рад)
        double argpo;    // Исходный аргумент перигея (рад)
        double mo;       // Исходняя средняя аномалия (рад)
        double bstar;    // Баллистический коэффициент

        // Рабочие элементы
        double ndot;     // Производная среднего движения
        double nddot;    // Вторая производная среднего движения
        double alta;     // Высота апогея (км)
        double altp;     // Высота перигея (км)
        double a;        // Большая полуось (ER)
        double del1;     // Первая поправка к n
        double del2;     // Вторая поправка к n
        double del3;     // Третья поправка к n
        double e;        // Текущий эксцентриситет
        double n;        // Текущее среднее движение
    };

    void initParameters(const TLEParser::TLEData& tle);
    void calculateDerivatives();
    void propagate(double tsince, QVector3D& pos, QVector3D& vel);
    double solveKepler(double mean_anomaly, double ecc);

    Elements elements_;
    QDateTime epoch_;
};

#endif // SGP4_PROPAGATOR_H
