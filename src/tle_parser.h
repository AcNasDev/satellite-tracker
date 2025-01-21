#ifndef TLE_PARSER_H
#define TLE_PARSER_H

#pragma once

#include <QString>
#include <QDateTime>

class TLEParser {
public:
    // Структура для хранения всех данных из TLE
    struct TLEData {
        // Данные из первой строки
        int satellite_number;        // Номер спутника (столбцы 3-7)
        char classification;         // Классификация (столбец 8)
        QString international_id;    // Международный идентификатор (столбцы 10-17)
        QDateTime epoch;            // Эпоха (столбцы 19-32)
        double first_derivative;     // Первая производная среднего движения (столбцы 34-43)
        double second_derivative;    // Вторая производная среднего движения (столбцы 45-52)
        double bstar;               // BSTAR drag term (столбцы 54-61)
        int ephemeris_type;         // Тип эфемерид (столбец 63)
        int element_number;         // Номер набора элементов (столбцы 65-68)
        int checksum1;              // Контрольная сумма строки 1 (столбец 69)

        // Данные из второй строки
        double inclination;         // Наклонение (столбцы 9-16)
        double right_ascension;     // Прямое восхождение узла (столбцы 18-25)
        double eccentricity;        // Эксцентриситет (столбцы 27-33)
        double argument_perigee;    // Аргумент перигея (столбцы 35-42)
        double mean_anomaly;        // Средняя аномалия (столбцы 44-51)
        double mean_motion;         // Среднее движение (столбцы 53-63)
        int revolution_number;      // Номер витка (столбцы 64-68)
        int checksum2;              // Контрольная сумма строки 2 (столбец 69)
    };

    static std::optional<TLEData> parseTLE(const QString& line1,
                                           const QString& line2,
                                           const QString& name = QString());

private:
    static bool validateChecksum(const QString& line);
    static QDateTime parseEpoch(const QString& epochStr);
    static double parseDecimal(const QString& str, int impliedDecimalPoint);
    static double parseExponential(const QString& mantissa, const QString& exponent);
};

#endif // TLE_PARSER_H
