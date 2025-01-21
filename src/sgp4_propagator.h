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
    // Константы SGP4
    static constexpr double xke = 0.0743669161331734049;  // sqrt(398600.8) / 6378.137^(3/2)
    static constexpr double xj2 = 1.082616e-3;           // J2 гармоника
    static constexpr double xj3 = -2.53881e-6;           // J3 гармоника
    static constexpr double xj4 = -1.65597e-6;           // J4 гармоника
    static constexpr double xkmper = 6378.137;           // Радиус Земли (км)
    static constexpr double ae = 1.0;                    // Радиус Земли (единицы Земли)
    static constexpr double de2ra = M_PI/180.0;          // Градусы в радианы
    static constexpr double min_per_day = 1440.0;        // Минут в сутках
    static constexpr double ck2 = xj2/2.0;
    static constexpr double ck4 = -3.0*xj4/8.0;

    struct Elements {
        // Начальные элементы орбиты
        double no;           // Среднее движение [рад/мин]
        double ecco;        // Эксцентриситет
        double inclo;       // Наклонение [рад]
        double nodeo;       // Долгота восходящего узла [рад]
        double argpo;       // Аргумент перигея [рад]
        double mo;          // Средняя аномалия [рад]
        double bstar;       // Баллистический коэффициент [1/радиус_земли]

        // Вычисленные параметры
        double a;           // Большая полуось [радиусы_земли]
        double ndot;        // Первая производная среднего движения [рад/мин^2]
        double nddot;       // Вторая производная среднего движения [рад/мин^3]
        double alta;        // Апогей [радиусы_земли]
        double altp;        // Перигей [радиусы_земли]
        double no_kozai;    // Среднее движение по Козаи [рад/мин]
    };

    void initParameters(const TLEParser::TLEData& tle);
    void propagate(double tsince, QVector3D& pos, QVector3D& vel) const;
    double solveKepler(double M, double e) const;

    Elements elements_;
    QDateTime epoch_;
};
#endif
