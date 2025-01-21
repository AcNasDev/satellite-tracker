#include "tle_parser.h"

// tle_parser.hpp остается таким же, как в предыдущем примере

// Добавим в tle_parser.cpp детальный вывод отладочной информации:

std::optional<TLEParser::TLEData> TLEParser::parseTLE(
    const QString& line1, const QString& line2, const QString& name) {

    qDebug() << "\n=== TLE Parser Debug Information ===";
    qDebug() << "Input lines:";
    qDebug() << "Line 1:" << line1;
    qDebug() << "Line 2:" << line2;
    qDebug() << "Name:" << name;

    if (line1.length() != 69 || line2.length() != 69) {
        qDebug() << "ERROR: Invalid line length. Expected 69 characters.";
        qDebug() << "Line 1 length:" << line1.length();
        qDebug() << "Line 2 length:" << line2.length();
        return std::nullopt;
    }

    if (!validateChecksum(line1) || !validateChecksum(line2)) {
        qDebug() << "ERROR: Checksum validation failed";
        return std::nullopt;
    }

    TLEData data;

    // Парсинг первой строки с отладкой
    qDebug() << "\n=== Line 1 Parsing ===";

    qDebug() << "Line number:" << line1.at(0) << "(raw:" << line1.mid(0, 1) << ")";

    data.satellite_number = line1.mid(2, 5).toInt();
    qDebug() << "Satellite Number:" << data.satellite_number << "(raw:" << line1.mid(2, 5) << ")";

    data.classification = line1.at(7).toLatin1();
    qDebug() << "Classification:" << data.classification << "(raw:" << line1.mid(7, 1) << ")";

    data.international_id = line1.mid(9, 8).trimmed();
    qDebug() << "International ID:" << data.international_id << "(raw:" << line1.mid(9, 8) << ")";

    data.epoch = parseEpoch(line1.mid(18, 14));
    qDebug() << "Epoch:" << data.epoch.toString("yyyy-MM-dd HH:mm:ss.zzz")
             << "(raw:" << line1.mid(18, 14) << ")";

    data.first_derivative = parseDecimal(line1.mid(33, 10), 8);
    qDebug() << "First Derivative:" << data.first_derivative
             << "(raw:" << line1.mid(33, 10) << ")";

    QString second_derivative = line1.mid(44, 6).trimmed();
    QString second_derivative_exp = line1.mid(50, 2);
    if (second_derivative == "00000" && second_derivative_exp == "-0") {
        data.second_derivative = 0.0;
        qDebug() << "Second Derivative explicitly set to 0 (special case)";
    } else {
        data.second_derivative = parseExponential(second_derivative, second_derivative_exp);
    }
    qDebug() << "Second Derivative:" << data.second_derivative
             << "(raw:" << second_derivative << "E" << second_derivative_exp << ")";

    QString bstar = line1.mid(53, 6).trimmed();
    QString bstar_exp = line1.mid(59, 2);
    data.bstar = parseExponential(bstar, bstar_exp);
    qDebug() << "BSTAR:" << data.bstar
             << "(raw:" << bstar << "E" << bstar_exp << ")";

    data.ephemeris_type = line1.mid(62, 1).toInt();
    qDebug() << "Ephemeris Type:" << data.ephemeris_type
             << "(raw:" << line1.mid(62, 1) << ")";

    data.element_number = line1.mid(64, 4).toInt();
    qDebug() << "Element Number:" << data.element_number
             << "(raw:" << line1.mid(64, 4) << ")";

    data.checksum1 = line1.mid(68, 1).toInt();
    qDebug() << "Checksum:" << data.checksum1
             << "(raw:" << line1.mid(68, 1) << ")";

    // Парсинг второй строки с отладкой
    qDebug() << "\n=== Line 2 Parsing ===";

    qDebug() << "Line number:" << line2.at(0) << "(raw:" << line2.mid(0, 1) << ")";

    data.inclination = line2.mid(8, 8).toDouble();
    qDebug() << "Inclination:" << data.inclination
             << "(raw:" << line2.mid(8, 8) << ")";

    data.right_ascension = line2.mid(17, 8).toDouble();
    qDebug() << "Right Ascension:" << data.right_ascension
             << "(raw:" << line2.mid(17, 8) << ")";

    data.eccentricity = parseDecimal(line2.mid(26, 7), 7);
    qDebug() << "Eccentricity:" << data.eccentricity
             << "(raw:" << line2.mid(26, 7) << ")";

    data.argument_perigee = line2.mid(34, 8).toDouble();
    qDebug() << "Argument of Perigee:" << data.argument_perigee
             << "(raw:" << line2.mid(34, 8) << ")";

    data.mean_anomaly = line2.mid(43, 8).toDouble();
    qDebug() << "Mean Anomaly:" << data.mean_anomaly
             << "(raw:" << line2.mid(43, 8) << ")";

    data.mean_motion = line2.mid(52, 11).toDouble();
    qDebug() << "Mean Motion:" << data.mean_motion
             << "(raw:" << line2.mid(52, 11) << ")";

    data.revolution_number = line2.mid(63, 5).toInt();
    qDebug() << "Revolution Number:" << data.revolution_number
             << "(raw:" << line2.mid(63, 5) << ")";

    data.checksum2 = line2.mid(68, 1).toInt();
    qDebug() << "Checksum:" << data.checksum2
             << "(raw:" << line2.mid(68, 1) << ")";

    qDebug() << "\n=== Checksum Validation ===";
    qDebug() << "Line 1 checksum valid:" << validateChecksum(line1);
    qDebug() << "Line 2 checksum valid:" << validateChecksum(line2);

    qDebug() << "\n=== Parsing Complete ===\n";

    return data;
}

// Добавим отладочную информацию в helper-функции:

bool TLEParser::validateChecksum(const QString& line) {
    int checksum = 0;
    for(int i = 0; i < line.length() - 1; i++) {
        QChar c = line[i];
        if(c.isDigit()) {
            checksum += c.digitValue();
        } else if(c == '-') {
            checksum += 1;
        }
    }
    checksum %= 10;

    int expected = line.right(1).toInt();
    qDebug() << "Checksum calculation for line:" << line;
    qDebug() << "Calculated checksum:" << checksum;
    qDebug() << "Expected checksum:" << expected;

    return checksum == expected;
}

double TLEParser::parseExponential(const QString& mantissa, const QString& exponent) {
    bool ok1, ok2;

    QString trimmedMantissa = mantissa.trimmed();
    QString trimmedExponent = exponent.trimmed();

    // Проверка на пустое значение или все нули
    if (trimmedMantissa.isEmpty() ||
        trimmedMantissa.count('0') == trimmedMantissa.length()) {
        qDebug() << "Zero or empty mantissa detected";
        return 0.0;
    }

    // Обработка знака в мантиссе
    bool negative = trimmedMantissa.startsWith('-');
    if (negative) {
        trimmedMantissa = trimmedMantissa.mid(1);
    }

    double value = trimmedMantissa.toDouble(&ok1);
    if (negative) {
        value = -value;
    }

    int exp = trimmedExponent.toInt(&ok2);

    qDebug() << "Parsing exponential number:";
    qDebug() << "Mantissa (trimmed):" << trimmedMantissa
             << "(negative:" << negative << ", valid:" << ok1 << ")";
    qDebug() << "Exponent (trimmed):" << trimmedExponent
             << "(valid:" << ok2 << ")";

    if(ok1 && ok2) {
        double result = value * std::pow(10.0, exp);
        qDebug() << "Result:" << result;
        return result;
    }

    qDebug() << "Parse failed, returning 0.0";
    return 0.0;
}

QDateTime TLEParser::parseEpoch(const QString& epochStr) {
    int year = epochStr.left(2).toInt();
    double dayOfYear = epochStr.mid(2).toDouble();

    year += (year < 57) ? 2000 : 1900;

    QDate date(year, 1, 1);
    date = date.addDays(int(dayOfYear) - 1);

    int milliseconds = int((dayOfYear - int(dayOfYear)) * 24 * 60 * 60 * 1000);
    QTime time = QTime(0, 0).addMSecs(milliseconds);

    return QDateTime(date, time, Qt::UTC);
}

double TLEParser::parseDecimal(const QString& str, int impliedDecimalPoint) {
    bool ok;
    QString trimmed = str.trimmed();

    // Обработка десятичных чисел с явной десятичной точкой
    if (trimmed.contains('.')) {
        return trimmed.toDouble(&ok);
    }

    // Обработка чисел с подразумеваемой десятичной точкой
    double value = trimmed.toDouble(&ok);
    if (ok) {
        return value * std::pow(10.0, -impliedDecimalPoint);
    }

    qDebug() << "Failed to parse decimal:" << str;
    return 0.0;
}
