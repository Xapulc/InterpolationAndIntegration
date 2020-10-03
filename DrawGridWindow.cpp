//
// Created by xapulc on 01.10.2020.
//

#include "DrawGridWindow.h"

I2Point DrawGridWindow::project(R2Point pt) {
    return I2Point(round(winSize.x * (pt.x - this->xmin) / (this->xmax - this->xmin)),
                   round(winSize.y * (this->ymax - pt.y) / (this->ymax - this->ymin)));
}

void DrawGridWindow::onExpose(XEvent&) {
    winSize.x = getWindowRect().width();
    winSize.y = getWindowRect().height();
    setForeground(getBackground());
    fillRectangle(m_RWinRect);
    draw();
}

void DrawGridWindow::draw() {
    setBackground("LightGray");
    drawGrid(true);
}

void DrawGridWindow::setGrid(Matrix &matrixX, Matrix &matrixY, Matrix &ih) {
    this->matrixX = matrixX;
    this->matrixY = matrixY;
    this->ih = ih;
    this->xmin = matrixX.min();
    this->xmax = matrixX.max();
    this->ymin = matrixY.min();
    this->ymax = matrixY.max();
}

void DrawGridWindow::setFunction(Matrix &fValues) {
    this->fValues = fValues;
}


void DrawGridWindow::setGradientColor(double min, double max, double value) {
    double t = (value - min) / (max - min);
    auto red = 105 + (int) (150 * t);
    auto green = 100 + (int) (100 * t);
    auto blue = 55 + (int) (200 * (1-t));
    setForeground(red, green, blue);
}

void DrawGridWindow::drawGrid(bool withIsolines) {
    auto N1 = matrixX.getShape(0);
    auto N2 = matrixX.getShape(1);
    auto fmin = withIsolines ? fValues.min() : 0;
    auto fmax = withIsolines ? fValues.max() : 0;

    setForeground("black");
    drawLine(project(R2Point(0, 0)),
             project(R2Point(10, 0)));
    drawLine(project(R2Point(1, -0.25)),
             project(R2Point(1, 0.25)));
    drawLine(project(R2Point(0, 0)),
             project(R2Point(0, 10)));
    drawLine(project(R2Point(-0.25, 1)),
             project(R2Point(0.25, 1)));
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            if (withIsolines) {
                setGradientColor(fmin, fmax, fValues(i, j));
            }
            if (i == 0) {
                if (j != 0) {
                    drawLine(project(R2Point(matrixX(i, j), matrixY(i, j))),
                             project(R2Point(matrixX(i, j-1), matrixY(i, j-1))));
                }
            } else {
                if (j == 0) {
                    drawLine(project(R2Point(matrixX(i, j), matrixY(i, j))),
                             project(R2Point(matrixX(i-1, j), matrixY(i-1, j))));
                } else {
                    if ((ih(i, j) > 0) && (ih(i-1, j) > 0)) {
                        drawLine(project(R2Point(matrixX(i, j), matrixY(i, j))),
                                 project(R2Point(matrixX(i-1, j), matrixY(i-1, j))));
                    }
                    if ((ih(i, j) > 0) && (ih(i, j-1) > 0)) {
                        drawLine(project(R2Point(matrixX(i, j), matrixY(i, j))),
                                 project(R2Point(matrixX(i, j-1), matrixY(i, j-1))));
                    }
                }
            }
        }
    }
}
