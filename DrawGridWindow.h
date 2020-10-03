//
// Created by xapulc on 01.10.2020.
//

#ifndef INTEGRATION_DRAWGRIDWINDOW_H
#define INTEGRATION_DRAWGRIDWINDOW_H

#include "Matrix.h"
#include "gwindow.h"

class DrawGridWindow: public GWindow {
private:
    Matrix matrixX;
    Matrix matrixY;
    Matrix ih;
    Matrix fValues;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    I2Point winSize = I2Point(getWindowRect().width(), getWindowRect().height());
    void setGradientColor(double min, double max, double value);

public:
    I2Point project(R2Point pt);
    virtual void onExpose(XEvent& event);
    void setGrid(Matrix& matrixX, Matrix& matrixY, Matrix& ih);
    void draw();
    void drawGrid(bool withIsolines=false);
    void setFunction(Matrix& fValues);
};


#endif //INTEGRATION_DRAWGRIDWINDOW_H
