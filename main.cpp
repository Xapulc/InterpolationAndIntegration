#include <iostream>

#include "Interpolation.h"
#include "DrawGridWindow.h"

double f(double x, double y) {
    return (4*x-y*y) * (1 - sin(x)) * cos(y);
}

int main(int argc, char* argv[]) {
    Matrix fValues;
    auto wideGridN1 = atoi(argv[1]);
    auto wideGridN2 = atoi(argv[2]);
    int wideGridmaxIter = atoi(argv[3]);
    double wideGridomega = atof(argv[4]);
    auto denseGridN1 = atoi(argv[5]);
    auto denseGridN2 = atoi(argv[6]);
    int denseGridmaxIter = atoi(argv[7]);
    double denseGridomega = atof(argv[8]);

    if (!GWindow::initX()) {
        std::cout << "Нет подключения к X-серверу." << std::endl;
        return -1;
    }

    double aspect = (double) GWindow::screenMaxX() /
                    (double) GWindow::screenMaxY();
    double width = 30.;
    int height = width / aspect;

    // Строим густую сетку
    auto denseGridX = Matrix(denseGridN1, denseGridN2);
    auto denseGridY = Matrix(denseGridN1, denseGridN2);
    auto denseGridIh = Matrix(denseGridN1, denseGridN2);
    Interpolation::buildGrid(denseGridX, denseGridY, denseGridIh, denseGridmaxIter, denseGridomega);
    std::cout << "Интеграл по густой сетке равен "
              << Interpolation::calculateIntegral(f, denseGridX, denseGridY, denseGridIh)
              << "."
              << std::endl;

    // Рисуем густую сетку

    DrawGridWindow denseGridWindow;
    denseGridWindow.createWindow(
            I2Rectangle(
                    I2Point(10, 10),
                    GWindow::screenMaxX()/3,
                    GWindow::screenMaxY()/3
            ),
            R2Rectangle(
                    R2Point(-width/2., -height/2.),
                    width, height
            ),
            "Густая сетка"
    );
    denseGridWindow.setGrid(denseGridX, denseGridY, denseGridIh);
    fValues = Matrix(denseGridN1, denseGridN2);
    for (int i = 0; i < denseGridN1; i++) {
        for (int j = 0; j < denseGridN2; j++) {
            fValues(i, j) = f(denseGridX(i, j), denseGridY(i, j));
        }
    }
    denseGridWindow.setFunction(fValues);

    // Строим редкую сетку
    auto wideGridX = Matrix(wideGridN1, wideGridN2);
    auto wideGridY = Matrix(wideGridN1, wideGridN2);
    auto wideGridIh = Matrix(wideGridN1, wideGridN2);
    Interpolation::buildGrid(wideGridX, wideGridY, wideGridIh, wideGridmaxIter, wideGridomega);
    std::cout << "Интеграл по редкой сетке равен "
              << Interpolation::calculateIntegral(f, wideGridX, wideGridY, wideGridIh)
              << "."
              << std::endl;

    // Рисуем редкую сетку

    DrawGridWindow wideGridWindow;
    wideGridWindow.createWindow(
            I2Rectangle(
                    I2Point(10, 10),
                    GWindow::screenMaxX()/3,
                    GWindow::screenMaxY()/3
            ),
            R2Rectangle(
                    R2Point(-width/2., -height/2.),
                    width, height
            ),
            "Редкая сетка"
    );
    wideGridWindow.setGrid(wideGridX, wideGridY, wideGridIh);
    fValues = Matrix(wideGridN1, wideGridN2);
    for (int i = 0; i < wideGridN1; i++) {
        for (int j = 0; j < wideGridN2; j++) {
            fValues(i, j) = f(wideGridX(i, j), wideGridY(i, j));
        }
    }
    wideGridWindow.setFunction(fValues);

    // Интерполяция с редкой сетки на густую
    fValues = Interpolation::gridInterpolation(f, wideGridX, wideGridY,
                                               denseGridX, denseGridY,
                                               wideGridIh, denseGridIh);
    std::cout << "Интеграл после интерполяции равен "
              << Interpolation::calculateIntegral(fValues, denseGridX, denseGridY, denseGridIh)
              << "."
              << std::endl;

    DrawGridWindow wideInterWindow;
    wideInterWindow.createWindow(
            I2Rectangle(
                    I2Point(10, 10),
                    GWindow::screenMaxX()/3,
                    GWindow::screenMaxY()/3
            ),
            R2Rectangle(
                    R2Point(-width/2., -height/2.),
                    width, height
            ),
            "Интерполяция с редкой сетки на густую"
    );
    wideInterWindow.setGrid(denseGridX, denseGridY, denseGridIh);
    wideInterWindow.setFunction(fValues);

    GWindow::messageLoop();

    GWindow::closeX();
}
