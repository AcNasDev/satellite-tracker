#include <QCoreApplication>
#include "SatelliteTracker.h"

int main(int argc, char *argv[])
{
    QCoreApplication app(argc, argv);

    QString tle1 = "1 57731U 23130E   25021.50910661  .00030604  00000-0  12673-2 0  9994";
    QString tle2 = "2 57731  34.9946 106.4394 0005335  72.4066 287.7231 15.22972654 77656";
    try {
        SatelliteTracker tracker(tle1, tle2);
        return app.exec();
    } catch (const std::exception& e) {
        qDebug() << "Error:" << e.what();
        return 1;
    }
}
