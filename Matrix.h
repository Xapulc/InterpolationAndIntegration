//
// Created by xapulc on 22.09.19.
//

#ifndef LINEARSOLVERS_MATRIX_H
#define LINEARSOLVERS_MATRIX_H

#include <vector>
#include <iostream>
#include <cmath>
#include <thread>
#include <future>


class Matrix {
    typedef std::vector<double> DoubleVector;
public:
    void createElems(int len);

    Matrix() = default;
    Matrix(int len);
    Matrix(int h, int w);
    static Matrix eye(int n);

    double& operator[](int i);
    double operator[](int i) const;
    double& operator()(int i, int j);
    double operator()(int i, int j) const;
    Matrix& operator=(const Matrix& other);
    Matrix& operator+=(const Matrix& other);
    Matrix operator+(const Matrix& other) const;
    Matrix& operator-=(const Matrix& other);
    Matrix operator-(const Matrix& other) const;
    Matrix& operator*=(const Matrix& other);
    Matrix operator*(const Matrix& other) const;
    Matrix& operator/=(const Matrix& other);
    Matrix operator/(const Matrix& other) const;
    Matrix& operator*=(double a);
    Matrix operator*(double a) const;
    Matrix& operator/=(double a);
    Matrix operator/(double a) const;
    Matrix dot(const Matrix& other) const;

    const int getShape(const int axis) const { return (axis == 0) ? h : w; }
    void print() const;
    double residual() const;
    Matrix t() const;
    double trace() const;
    double min() const;
    double max() const;
private:
    int h{0};
    int w{0};
    DoubleVector elems;
    static void checkSize(int size);
    static void checkShape(const Matrix& left, const Matrix& right);
};


#endif //LINEARSOLVERS_MATRIX_H
