//
// Created by xapulc on 22.09.19.
//

#include "Matrix.h"

void Matrix::checkSize(const int size) {
    if (size <= 0) {
        std::cout << "Wrong size: " << size << std::endl;
        exit(-1);
    }
}

void Matrix::createElems(const int len) {
    if ((this->elems = DoubleVector(static_cast<unsigned long>(len))).empty()) {
        std::cout << "Cannot allocate enough memory. Len: " << len << std::endl;
        exit(-1);
    }
    for(int i = 0; i < len; i++)
        elems[i] = 0;
}

Matrix::Matrix(const int len) {
    checkSize(len);
    this->h = len;
    this->w = 1;
    createElems(len);
}

Matrix::Matrix(const int h, const int w) {
    checkSize(h);
    checkSize(w);
    this->h = h;
    this->w = w;
    createElems(h*w);
}

Matrix &Matrix::operator=(const Matrix &other) {
    h = other.h;
    w = other.w;
    createElems(h*w);
    for(int i = 0; i < h * w; i++)
        elems[i] = other.elems[i];
    return *this;
}

double &Matrix::operator[](const int i) {
    return elems[i];
}

double Matrix::operator[](const int i) const {
    return elems[i];
}

double &Matrix::operator()(const int i, const int j) {
    return elems[i*w + j];
}

double Matrix::operator()(const int i, const int j) const {
    return elems[i*w + j];
}

void Matrix::checkShape(const Matrix& left, const Matrix& right) {
    if ((left.h != right.h) || (left.w != right.w)) {
        std::cout << "Wrong shapes: (" << left.h << ", " << left.w <<
                  ") and (" << right.h << ", " << right.w << ")" << std::endl;
        exit(-1);
    }
}

Matrix& Matrix::operator+=(const Matrix &other) {
    checkShape(*this, other);
    for(int i = 0; i < h*w; i++)
        elems[i] += other.elems[i];
    return *this;
}

Matrix Matrix::operator+(const Matrix &other) const {
    auto res = *this;
    res += other;
    return res;
}

Matrix& Matrix::operator-=(const Matrix &other) {
    checkShape(*this, other);
    for(int i = 0; i < h*w; i++)
        elems[i] -= other.elems[i];
    return *this;
}

Matrix Matrix::operator-(const Matrix &other) const {
    auto res = *this;
    res -= other;
    return res;
}

Matrix& Matrix::operator*=(const Matrix &other) {
    checkShape(*this, other);
    for(int i = 0; i < h*w; i++)
        elems[i] *= other.elems[i];
    return *this;
}

Matrix Matrix::operator*(const Matrix &other) const {
    auto res = *this;
    res *= other;
    return res;
}

Matrix& Matrix::operator/=(const Matrix &other) {
    checkShape(*this, other);
    for(int i = 0; i < h*w; i++)
        elems[i] /= other.elems[i];
    return *this;
}

Matrix Matrix::operator/(const Matrix &other) const {
    auto res = *this;
    res /= other;
    return res;
}

Matrix& Matrix::operator*=(const double a) {
    for(int i = 0; i < h*w; i++)
        elems[i] *= a;
    return *this;
}

Matrix Matrix::operator*(const double a) const {
    auto res = *this;
    res *= a;
    return res;
}

Matrix& Matrix::operator/=(const double a) {
    for(int i = 0; i < h*w; i++)
        elems[i] /= a;
    return *this;
}

Matrix Matrix::operator/(const double a) const {
    auto res = *this;
    res /= a;
    return res;
}

Matrix Matrix::dot(const Matrix &other) const {
    if (this->w != other.h) {
        std::cout << "Wrong shapes: (" << this->h << ", " << this->w << ") and (" <<
                  other.h << ", " << other.w << ")" << std::endl;
        exit(-1);
    }
    auto res = Matrix(this->h, other.w);

    for(int i = 0; i < this->h; i++)
        for(int j = 0; j < other.w; j++)
            for(int k = 0; k < this->w; k++)
                res(i, j) += (*this)(i, k) * other(k, j);
    return res;
}

Matrix Matrix::eye(const int n) {
    auto res = Matrix(n, n);
    for(int i = 0; i < n; i++)
        res(i, i) = 1;
    return res;
}

void Matrix::print() const {
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++)
            std::cout << (*this)(i, j) << " ";
        std::cout << std::endl;
    }
}

double Matrix::residual() const {
    double sum = 0;
    for(int i = 0; i < h*w; i++)
        sum += elems[i] * elems[i];
    return sum;
}

Matrix Matrix::t() const {
    Matrix res = Matrix(w, h);
    for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++)
            res(j, i) = (*this)(i, j);
    return res;
}

double Matrix::trace() const {
    double trace = 0;
    for(int i = 0; i < ((h > w) ? w : h); i++)
        trace += (*this)(i, i);
    return trace;
}

double Matrix::min() const {
    double min = (*this)(0, 0);
    for (int i = 0; i < this->getShape(0); i++) {
        for (int j = 0; j < this->getShape(1); j++) {
            if ((*this)(i, j) < min) {
                min = (*this)(i, j);
            }
        }
    }
    return min;
}

double Matrix::max() const {
    double max = (*this)(0, 0);
    for (int i = 0; i < this->getShape(0); i++) {
        for (int j = 0; j < this->getShape(1); j++) {
            if ((*this)(i, j) > max) {
                max = (*this)(i, j);
            }
        }
    }
    return max;
}
