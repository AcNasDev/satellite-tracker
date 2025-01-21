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
    OrbitalState calculateState(const QDateTime& time) const;

private:
    // SGP4 константы
    static constexpr double ae = 1.0;
    static constexpr double tothrd = 2.0/3.0;
    static constexpr double xkmper = 6378.137;
    static constexpr double f = 1.0/298.257223563;
    static constexpr double j2 = 1.082616e-3;
    static constexpr double j3 = -2.53881e-6;
    static constexpr double j4 = -1.65597e-6;
    static constexpr double ke = 7.43669161e-2;
    static constexpr double ck2 = j2/2.0;
    static constexpr double ck4 = -3.0*j4/8.0;
    static constexpr double xj3 = j3;
    static constexpr double qo = ae + 120.0/xkmper;
    static constexpr double s = ae + 78.0/xkmper;
    static constexpr double min_per_day = 1440.0;
    static constexpr double de2ra = M_PI/180.0;

    struct Elements {
        // Исходные элементы
        double no;     // Начальное среднее движение [rad/min]
        double eo;     // Начальный эксцентриситет
        double io;     // Начальное наклонение [rad]
        double omegao; // Начальный аргумент перигея [rad]
        double xnodeo; // Начальная долгота восходящего узла [rad]
        double xmo;    // Начальная средняя аномалия [rad]
        double bstar;  // Баллистический коэффициент [1/earth radii]

        // Производные элементы
        double a;      // Большая полуось [earth radii]
        double ndot;   // Скорость изменения среднего движения [rad/min^2]
        double nddot;  // Ускорение среднего движения [rad/min^3]
        double inclo;  // Наклонение [rad]
        double nodeo;  // Долгота восходящего узла [rad]
        double argpo;  // Аргумент перигея [rad]
        double mo;     // Средняя аномалия [rad]
        double no_kozai; // Среднее движение по Козаи [rad/min]
    };

    void initParameters(const TLEParser::TLEData& tle);
    void deepSpaceInitialize();
    void propagate(double tsince, QVector3D& pos, QVector3D& vel) const;

    Elements elements_;
    QDateTime epoch_;
    bool is_deep_space_;
};

#endif // SGP4_PROPAGATOR_H
