//
// Created by xapulc on 29.09.2020.
//

#ifndef INTEGRATION_POINT_H
#define INTEGRATION_POINT_H


class Point {
public:
    Point() = default;
    Point(double x, double y);

    Point& operator-=(const Point& other);
    Point operator-(const Point& other) const;

    double x;
    double y;
};


#endif //INTEGRATION_POINT_H
