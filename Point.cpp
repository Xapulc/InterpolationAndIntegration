//
// Created by xapulc on 29.09.2020.
//

#include "Point.h"

Point::Point(double x, double y) {
    this->x = x;
    this->y = y;
}

Point &Point::operator-=(const Point &other) {
    x -= other.x;
    y -= other.y;
    return *this;
}

Point Point::operator-(const Point &other) const {
    auto res = *this;
    res -= other;
    return res;
}
