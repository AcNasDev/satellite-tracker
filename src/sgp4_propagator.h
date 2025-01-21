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
    static constexpr double xke = 0.0743669161331734049;
    static constexpr double xj2 = 1.082616e-3;
    static constexpr double xj3 = -2.53881e-6;
    static constexpr double xj4 = -1.65597e-6;
    static constexpr double xkmper = 6378.137;
    static constexpr double ae = 1.0;
    static constexpr double de2ra = M_PI/180.0;
    static constexpr double min_per_day = 1440.0;
    static constexpr double ck2 = xj2/2.0;
    static constexpr double ck4 = -3.0*xj4/8.0;
    static constexpr double q0 = 120.0/xkmper;
    static constexpr double s0 = 78.0/xkmper;

    struct Elements {
        double no;         // Среднее движение [рад/мин]
        double ecco;       // Эксцентриситет
        double inclo;      // Наклонение [рад]
        double nodeo;      // Долгота восходящего узла [рад]
        double argpo;      // Аргумент перигея [рад]
        double mo;         // Средняя аномалия [рад]
        double bstar;      // Баллистический коэффициент [1/ER]

        // Производные элементы
        double a;          // Большая полуось [ER]
        double ndot;       // Первая производная среднего движения [рад/мин^2]
        double nddot;      // Вторая производная среднего движения [рад/мин^3]
        double alta;       // Апогей [ER]
        double altp;       // Перигей [ER]
        double no_kozai;   // Среднее движение Козаи [рад/мин]

        // Дополнительные параметры
        double aodp;       // Первоначальная большая полуось [ER]
        double aycof;      // Коэффициент для y
        double xlcof;      // Коэффициент для l
    };

    void initParameters(const TLEParser::TLEData& tle);
    void propagate(double tsince, QVector3D& pos, QVector3D& vel) ;
    double solveKepler(double M, double e) ;

    Elements elements_;
    QDateTime epoch_;
};
#endif
