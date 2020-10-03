//
// Created by xapulc on 01.10.2020.
//

#ifndef INTEGRATION_INTERPOLATION_H
#define INTEGRATION_INTERPOLATION_H

#include "Point.h"
#include "Matrix.h"

class Interpolation {
public:
    static void Laplasiation(Matrix& matrixX, Matrix& matrixY, Matrix& ih, int max_iter, double omega);
    static void buildGrid(Matrix& matrixX, Matrix& matrixY, Matrix& ih, int max_iter, double omega);
    static double calculateIntegral(double (*f)(double, double), Matrix& matrixX, Matrix& matrixY, Matrix& ih);
    static Matrix gridInterpolation(double (*f)(double, double),
                                    Matrix& wideGridX, Matrix& wideGridY,
                                    Matrix& denseGridX, Matrix& denseGridY,
                                    Matrix& wideGridIh, Matrix& denseGridIh);
    static double calculateIntegral(Matrix& fValues, Matrix& matrixX, Matrix& matrixY, Matrix& ih);
private:
    static void fill_arc(Matrix& matrixX, Matrix& matrixY, int i1, int i3, int j1, int j3, Point p1, Point p2, Point p3);
    static void contour(Matrix& matrixX, Matrix& matrixY);
    static void linearInterpolation(Matrix& matrixX, Matrix& matrixY);
    static double quadraticElement(double (*f)(double, double), Point a, Point b, Point c, Point d);
    static double interpolationElement(Point a, Point b, Point c, Point d,
                                       double fA, double fB, double fC, double fD);
    static double f_interp(Point a, Point b, Point c, Point d,
                           double fA, double fB, double fC, double fD);
};


#endif //INTEGRATION_INTERPOLATION_H
