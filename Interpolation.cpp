//
// Created by xapulc on 01.10.2020.
//

#include "Interpolation.h"


void Interpolation::fillArc(Matrix& matrixX, Matrix& matrixY, int i1, int i3, int j1, int j3, Point p1, Point p2, Point p3) {
    // По трём точкам, не лежащим на одной прямой, строится дуга A окружности omega,
    // проходящей через p1, p2, p3, причём дуга A соединяет точки p1 и p3.
    // Вектор X-координат дуги A заполняет матрицу matrixX X-координат
    // с (i1, j1) координаты до (i3, j3) координаты.
    // Аналогично, вектор Y-координат дуги A заполняет матрицу matrixY Y-координат
    // с (i1, j1) координаты до (i3, j3) координаты.
    // Заметим, что согласно описанию либо i1 = i3, либо j1 = j3.

    if (((i1 != i3) && (j1 != j3))
        || matrixX.getShape(0) != matrixY.getShape(0)
        || matrixX.getShape(1) != matrixY.getShape(1)
        || i1 < 0 || matrixX.getShape(0) < i1
        || i3 < 0 || matrixX.getShape(0) < i3
        || j1 < 0 || matrixX.getShape(1) < j1
        || j3 < 0 || matrixX.getShape(1) < j3) {
        std::cout << "Ошибка в размерностях или индексах" << std::endl;
        return;
    }

    const double Pi = atan(1) * 4;

    // Заполняем крайние точки.
    matrixX(i1-1, j1-1) = p1.x;
    matrixY(i1-1, j1-1) = p1.y;
    matrixX(i3-1, j3-1) = p3.x;
    matrixX(i3-1, j3-1) = p3.y;

    auto diff21 = p2 - p1;
    auto diff32 = p3 - p2;

    auto f1 = -0.5 * p1.x * p1.x - 0.5 * p1.y * p1.y
              + 0.5 * p2.x * p2.x + 0.5 * p2.y * p2.y;
    auto f2 = -0.5 * p2.x * p2.x - 0.5 * p2.y * p2.y
              + 0.5 * p3.x * p3.x + 0.5 * p3.y * p3.y;
    auto det = diff21.x * diff32.y - diff21.y * diff32.x;

    if (det == 0) {
        std::cout << "Детерминант преобразования равен нулю. Проверьте входные данные." << std::endl;
        return;
    }

    // Находим центр окружности omega и её радиус.
    auto p0 = Point((f1 * diff32.y - f2 * diff21.y) / det,
                    (f2 * diff21.x - f1 * diff32.x) / det);
    auto R = sqrt((p1.x - p0.x) * (p1.x - p0.x)
                  + (p1.y - p0.y) * (p1.y - p0.y));

    // Угол от центра окружности omega до p1.
    double alpha0;
    if ((0 <= p1.x - p0.x) && (0 <= p1.y - p0.y)) {
        alpha0 = acos((p1.x - p0.x) / R);
    } else if ((0 <= p1.x - p0.x) && (p1.y - p0.y <= 0)) {
        alpha0 = 2 * Pi - acos((p1.x - p0.x) / R);
    } else if ((p1.x - p0.x <= 0) && (0 <= p1.y - p0.y)) {
        alpha0 = acos((p1.x - p0.x) / R);
    } else {
        alpha0 = Pi + acos(-(p1.x - p0.x) / R);
    }

    auto alpha1 = acos(((p1.x - p0.x) * (p2.x - p0.x)
                        + (p1.y - p0.y) * (p2.y - p0.y)) / (R * R));
    auto alpha2 = acos(((p2.x - p0.x) * (p3.x - p0.x)
                        + (p2.y - p0.y) * (p3.y - p0.y)) / (R * R));
    auto d = (p1.x - p0.x) * (p2.y - p0.y)
             - (p1.y - p0.y) * (p2.x - p0.x);
    // Угол от центра окружности omega до p3.
    auto alpha3 = (0 <= d) ? (alpha1 + alpha2) : (- alpha1 - alpha2);

    // Генерируем равномерно точки на дуге A и записываем их в матрицы.
    if (j1 == j3) {
        auto current_i = i1;
        auto point_count = abs(i3 - i1);
        auto delta_alpha = alpha3 / point_count;

        for(int i = 0; i < point_count; i++) {
            current_i += (i3 - i1 < 0) ? -1 : 1;
            matrixX(current_i-1, j1-1) = p0.x + R * cos(alpha0 + delta_alpha * abs(current_i - i1));
            matrixY(current_i-1, j1-1) = p0.y + R * sin(alpha0 + delta_alpha * abs(current_i - i1));
        }
    } else {
        auto current_j = j1;
        auto point_count = abs(j3 - j1);
        auto delta_alpha = alpha3 / point_count;

        for(int i = 0; i < point_count; i++) {
            current_j += (j3 - j1 < 0) ? -1 : 1;
            matrixX(i1-1, current_j-1) = p0.x + R * cos(alpha0 + delta_alpha * abs(current_j - j1));
            matrixY(i1-1, current_j-1) = p0.y + R * sin(alpha0 + delta_alpha * abs(current_j - j1));
        }
    }
}

void Interpolation::contour(Matrix& matrixX, Matrix& matrixY) {
    auto N1 = matrixX.getShape(0);
    auto N2 = matrixX.getShape(1);

    auto a = 4.2;
    auto b = 4.8;
    auto c = Point(0.03, 0.03);

    fillArc(matrixX, matrixY, 1, N1, N2, N2,
            Point(-a + c.x, a + c.y), Point(0 + c.x, b + c.y), Point(a + c.x, a + c.y));
    fillArc(matrixX, matrixY, 1, N1, 1, 1,
            Point(-a + c.x, -a + c.y), Point(0 + c.x, -b + c.y), Point(a + c.x, -a + c.y));
    fillArc(matrixX, matrixY, 1, 1, 1, N2,
            Point(-a + c.x, -a + c.y), Point(-b + c.x, 0 + c.y), Point(-a + c.x, a + c.y));
    fillArc(matrixX, matrixY, N1, N1, 1, N2,
            Point(a + c.x, -a + c.y), Point(b + c.x, 0 + c.y), Point(a + c.x, a + c.y));
}

void Interpolation::linearInterpolation(Matrix& matrixX, Matrix& matrixY) {
    auto N1 = matrixX.getShape(0);
    auto N2 = matrixX.getShape(1);

    auto r1 = 1.0 / (N1-1);
    auto r2 = 1.0 / (N2-1);
    for (int i = 1; i < N1-1; i++) {
        for (int j = 1; j < N2-1; j++) {
            matrixX(i, j) = (1 - i*r1) * matrixX(0, j) + (1 - j*r2) * matrixX(i, 0)
                            + i * r1 * matrixX(N1-1, j) + j * r2 * matrixX(i, N2-1)
                            - (1 - i*r1) * (1 - j*r2) * matrixX(0, 0)
                            - (1 - i*r1) * j * r2 * matrixX(0, N2-1)
                            - i * r1 * (1 - j*r2) * matrixX(N1-1, 0)
                            - i * r1 * j * r2 * matrixX(N1-1, N2-1);
            matrixY(i, j) = (1 - i*r1) * matrixY(0, j) + (1 - j*r2) * matrixY(i, 0)
                            + i * r1 * matrixY(N1-1, j) + j * r2 * matrixY(i, N2-1)
                            - (1 - i*r1) * (1 - j*r2) * matrixY(0, 0)
                            - (1 - i*r1) * j * r2 * matrixY(0, N2-1)
                            - i * r1 * (1 - j*r2) * matrixY(N1-1, 0)
                            - i * r1 * j * r2 * matrixY(N1-1, N2-1);
        }
    }
}

void Interpolation::Laplasiation(Matrix& matrixX, Matrix& matrixY, Matrix& ih, int max_iter, double omega) {
    auto N1 = matrixX.getShape(0);
    auto N2 = matrixX.getShape(1);

    auto x11 = (double *) malloc(sizeof(double)*(N1*N2));
    auto x12 = (double *) malloc(sizeof(double)*(N1*N2));
    auto x21 = (double *) malloc(sizeof(double)*(N1*N2));
    auto x22 = (double *) malloc(sizeof(double)*(N1*N2));
    auto yac = (double *) malloc(sizeof(double)*(N1*N2));

    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            x11[i*N2 + j] = 0;
            x12[i*N2 + j] = 0;
            x21[i*N2 + j] = 0;
            x22[i*N2 + j] = 0;
            yac[i*N2 + j] = 0;
        }
    }

    for(int iteration = 0; iteration < max_iter; iteration++) {
        for (int i = 1; i < N1-1; i++) {
            for (int j = 1; j < N2-1; j++) {
                if (ih(i, j) < 3) {
                    continue;
                }
                x11[i*N2 + j] = (matrixX(i+1, j) - matrixX(i-1, j)) / 2;
                x12[i*N2 + j] = (matrixX(i, j+1) - matrixX(i, j-1)) / 2;
                x21[i*N2 + j] = (matrixY(i+1, j) - matrixY(i-1, j)) / 2;
                x22[i*N2 + j] = (matrixY(i, j+1) - matrixY(i, j-1)) / 2;
            }
        }

        double errorX = 0;
        double errorY = 0;
        double sumError = 0;
        for (int i = 0; i < N1-1; i++) {
            for (int j = 0; j < N2-1; j++) {
                if (ih(i, j) < 3) {
                    continue;
                }
                double ccc = 0.2E1 * pow(x22[i*N2 + j], 2)
                             + 0.2E1 * pow(x12[i*N2 + j], 2)
                             + 0.2E1 * pow(x21[i*N2 + j], 2)
                             + 0.2E1 * pow(x11[i*N2 + j], 2);
                if (fabs(ccc) == 0) {
                    ccc = 0.1E-4;
                }

                errorX = matrixX(i, j);
                matrixX(i, j) = -omega * matrixX(i, j)
                                + (1+omega) * ((pow(x22[i*N2 + j], 2) + pow(x12[i*N2 + j], 2))
                                               * (matrixX(i+1, j) + matrixX(i-1, j))
                                               - 0.5 * (x22[i*N2 + j] * x21[i*N2 + j] + x11[i*N2 + j] * x12[i*N2 + j])
                                                 * (matrixX(i+1, j+1) - matrixX(i+1, j-1) - matrixX(i-1, j+1) + matrixX(i-1, j-1))
                                               + (pow(x21[i*N2 + j], 2) + pow(x11[i*N2 + j], 2))
                                                 * (matrixX(i, j+1) + matrixX(i, j-1))) / ccc;
                errorX = pow(errorX - matrixX(i, j), 2);

                errorY = matrixY(i, j);
                matrixY(i, j) = -omega * matrixY(i, j)
                                + (1+omega) * ((pow(x22[i*N2 + j], 2) + pow(x12[i*N2 + j], 2))
                                               * (matrixY(i+1, j) + matrixY(i-1, j))
                                               - 0.5 * (x22[i*N2 + j] * x21[i*N2 + j] + x11[i*N2 + j] * x12[i*N2 + j])
                                                 * (matrixY(i+1, j+1) - matrixY(i+1, j-1) - matrixY(i-1, j+1) + matrixY(i-1, j-1))
                                               + (pow(x21[i*N2 + j], 2) + pow(x11[i*N2 + j], 2))
                                                 * (matrixY(i, j+1) + matrixY(i, j-1))) / ccc;
                errorY = pow(errorY - matrixY(i, j), 2);

                sumError += sqrt(errorX + errorY);
            }
        }
        sumError /= N1 * N2;
    }

    double yacmin;
    double yacmax;

    for(int i = 0; i < N1-1; i++) {
        for(int j = 0; j < N2-1; j++) {
            if (ih(i, j) < 2) {
                continue;
            }
            yac[i*N2 + j] = 0.5 * (matrixX(i+1, j) - matrixX(i, j)) * (matrixY(i, j+1) - matrixY(i, j))
                            - 0.5 * (matrixX(i, j+1) - matrixX(i, j)) * (matrixY(i+1, j) - matrixY(i, j))
                            + 0.5 * (matrixX(i+1, j+1) - matrixX(i, j+1)) * (matrixY(i+1, j+1) - matrixY(i+1, j))
                            - 0.5 * (matrixX(i+1, j+1) - matrixX(i+1, j)) * (matrixY(i+1, j+1) - matrixY(i, j+1));
            if ((i == 0) && (j == 0)) {
                yacmin = yac[i*N2 + j];
                yacmax = yac[i*N2 + j];
            } else {
                if (yac[i*N2 + j] < yacmin) {
                    yacmin = yac[i*N2 + j];
                }

                if (yac[i*N2 + j] > yacmax) {
                    yacmax = yac[i*N2 + j];
                }
            }
        }
    }

    yacmin *= N1 * N2;
    yacmax *= N1 * N2;
    std::cout << "Минимум якобиана равен " << yacmin
              << ", максимум якобиана равен " << yacmax << "." << std::endl;

    free(x11);
    free(x12);
    free(x21);
    free(x22);
}

void Interpolation::buildGrid(Matrix& matrixX, Matrix& matrixY, Matrix& ih, int max_iter, double omega) {
    auto N1 = matrixX.getShape(0);
    auto N2 = matrixX.getShape(1);

    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            ih(i, j) = 1;
        }
    }

    for (int i = 1 + 2*(N1-1) / 5; i < 3*(N1-1) / 5; i++) {
        for (int j = 1 + 2*(N2-1) / 5; j < 3*(N2-1) / 5; j++) {
            ih(i, j) = 0;
        }
    }

    for (int i = 0; i < N1-1; i++) {
        for (int j = 0; j < N2-1; j++) {
            if (ih(i, j) == 0) {
                continue;
            }

            if (i == 0) {
                if (ih(0, j+1) == 1) {
                    ih(i, j) = 2;
                }
                continue;
            }

            if (j == 0) {
                if (ih(i+1, 0) == 1) {
                    ih(i, j) = 2;
                }
                continue;
            }

            double ihx = ih(i+1, j) * ih(i, j+1) * ih(i+1, j+1);
            if (ihx == 0) {
                continue;
            }
            ih(i, j) = 2;

            ihx *= ih(i-1, j) * ih(i, j-1) * ih(i-1, j-1)
                   * ih(i+1, j-1) * ih(i-1, j-1) * ih(i+1, j+1);
            if (ihx == 0) {
                continue;
            }
            ih(i, j) = 3;
        }
    }

    contour(matrixX, matrixY);
    linearInterpolation(matrixX, matrixY);

    const double Pi = atan(1) * 4;

    auto k = 2.3;
    for (int i = 2*(N1-1) / 5; i <= 3*(N1-1) / 5; i++) {
        int j = 3*(N2-1) / 5 - 1;
        double t = Pi / 2 + atan(2)
                   - 2 * atan(2) * (i - 2*(N1-1) / 5) / (3*(N1-1) / 5 - 2*(N1-1) / 5);
        matrixX(i, j+1) = cos(t) / k;
        matrixY(i, j+1) = sin(t);
    }

    for (int j = 2*(N2-1) / 5; j <= 3*(N2-1) / 5; j++) {
        int i = 2*(N1-1) / 5 - 1;
        double t = Pi + atan(2)
                   - 2 * atan(2) * (j - 2*(N2-1) / 5) / (3*(N2-1) / 5 - 2*(N2-1) / 5);
        matrixX(i+1, j) = 0 + cos(t);
        matrixY(i+1, j) = 0 + sin(t) / k;
    }


    for (int i = 3*(N1-1) / 5; 2*(N1-1) / 5 <= i; i--) {
        int j = 2*(N2-1) / 5 - 1;
        double t = -Pi / 2 + atan(2)
                   + 2 * atan(2) * (i - 3*(N1-1) / 5) / (3*(N1-1) / 5 - 2*(N1-1) / 5);
        matrixX(i, j+1) = 0 + cos(t) / k;
        matrixY(i, j+1) = 0 + sin(t);
    }

    for (int j = 3*(N2-1) / 5; 2*(N2-1) / 5 <= j; j--) {
        int i = 3*(N1-1) / 5 - 1;
        double t = atan(2)
                   + 2 * atan(2) * (j - 3*(N2-1) / 5) / (3*(N2-1) / 5 - 2*(N2-1) / 5);
        matrixX(i+1, j) = 0 + cos(t);
        matrixY(i+1, j) = 0 + sin(t) / k;
    }

    Laplasiation(matrixX, matrixY, ih, max_iter, omega);
}

double Interpolation::quadraticElement(double (*f)(double, double), Point a, Point b, Point c, Point d) {
    auto Copt = -c.x*d.y/2.0+c.x*b.y/2.0+a.x*d.y/2.0-a.x*b.y/2.0
                -d.x*a.y/2.0+d.x*c.y/2.0+b.x*a.y/2.0-b.x*c.y/2.0;
    auto E1 = (a.x*a.x*b.y-a.x*a.x*d.y-a.x*a.y*b.x+a.x*a.y*d.x
               +a.x*b.x*b.y-a.x*d.x*d.y-a.y*b.x*b.x+a.y*d.x*d.x
               +b.x*b.x*c.y-b.x*b.y*c.x+b.x*c.x*c.y-b.y*c.x*c.x
               +c.x*c.x*d.y-c.x*c.y*d.x+c.x*d.x*d.y-c.y*d.x*d.x)
              / (a.x*b.y-a.x*d.y-b.x*a.y+d.x*a.y+b.x*c.y
                 -c.x*b.y+c.x*d.y-d.x*c.y);
    auto E2 = (a.x*a.y*b.y-a.x*a.y*d.y+a.x*b.y*b.y-a.x*d.y*d.y
               -a.y*a.y*b.x+a.y*a.y*d.x-a.y*b.x*b.y+a.y*d.x*d.y
               +b.x*b.y*c.y+b.x*c.y*c.y-b.y*b.y*c.x-b.y*c.x*c.y
               +c.x*c.y*d.y+c.x*d.y*d.y-c.y*c.y*d.x-c.y*d.x*d.y)
              /(a.x*b.y-a.x*d.y-b.x*a.y+d.x*a.y+b.x*c.y-c.x*b.y+c.x*d.y-d.x*c.y);

    return Copt*f(E1 / 3,E2 / 3);
}

Matrix Interpolation::gridInterpolation(double (*f)(double, double),
                                        Matrix& wideGridX, Matrix& wideGridY,
                                        Matrix& denseGridX, Matrix& denseGridY,
                                        Matrix& wideGridIh, Matrix& denseGridIh) {
    int n1 = wideGridX.getShape(0);
    int n2 = wideGridX.getShape(1);
    int N1 = denseGridX.getShape(0);
    int N2 = denseGridX.getShape(1);

    Matrix fValues = Matrix(N1, N2);
    Matrix ihmni = Matrix(N1, N2);
    Matrix ihmnj = Matrix(N1, N2);
    Matrix ihmnkey = Matrix(N1, N2);

    double Gamma;
    double Ha;
    double Hb;
    double Hc;
    double S2abc;
    double a1;
    double a2;
    int abc;
    int abd;
    int acd;
    double alpha;
    double b1;
    double b2;
    int bcd;
    double beta;
    int breake;
    double c1;
    double c2;
    double d1;
    double d2;
    double differ_abc;
    double differ_abd;
    double differ_acd;
    double differ_bcd;
    double differ_min;
    double e1;
    double e2;
    int i;
    int ibegin;
    int iend;
    int j;
    int jbegin;
    int jend;
    int key;
    int m;
    int n;
    int none;

    none = 0;
    abc = 1;
    bcd = 2;
    acd = 3;
    abd = 4;
    ibegin = 1;
    iend = n1;
    jbegin = 1;
    jend = n2;
    for(m = 1;m <= N1;m++)
    {
        for (n = 1;n <= N2;n++)
        {
            if (denseGridIh(-1+m, -1+n) == 0)
            {
                ihmni(-1+m, -1+n) = 0;
                ihmnj(-1+m, -1+n)= 0;
                ihmnkey(-1+m, -1+n) = 0;
                continue;
            }
            e1 = denseGridX(-1+m, -1+n);
            e2 = denseGridY(-1+m, -1+n);
            differ_min = 100000.0;
            breake = 0;
            for (i = ibegin;breake < 1 && i <= iend;i++)
                for (j = 1;breake < 1 && j <= n2;j++) {
                    if (wideGridIh(-1+i, j-1) < 2) {
                        continue;
                    }
                    a1 = wideGridX(-1+i, j-1);
                    a2 = wideGridY(-1+i, j-1);
                    b1 = wideGridX(-1+i, j);
                    b2 = wideGridY(-1+i, j);
                    c1 = wideGridX(i, j);
                    c2 = wideGridY(i, j);
                    d1 = wideGridX(i, j-1);
                    d2 = wideGridY(i, j-1);
                    differ_abc = fabs(fabs(-0.5*a1*b2+0.5*c2*a1-0.5*c1*a2+0.5*a2*b1-0.5
                                                                                    *b1*c2+0.5*b2*c1)-fabs(0.5*a1*e2-0.5*a1*b2+0.5*a2*b1-0.5*a2*e1+0.5*e1*b2-0.5*e2
                                                                                                                                                             *b1)-fabs(0.5*e2*b1-0.5*b1*c2+0.5*b2*c1-0.5*e1*b2+0.5*e1*c2-0.5*e2*c1)-fabs(
                            -0.5*e2*c1+0.5*c1*a2-0.5*c2*a1+0.5*e1*c2-0.5*a2*e1+0.5*a1*e2));
                    differ_abd = fabs(fabs(0.5*a1*b2-0.5*d2*a1+0.5*d1*a2-0.5*a2*b1+0.5*
                                                                                   b1*d2-0.5*b2*d1)-fabs(0.5*a1*e2-0.5*a1*b2+0.5*a2*b1-0.5*a2*e1+0.5*e1*b2-0.5*e2*
                                                                                                                                                           b1)-fabs(0.5*e2*b1-0.5*b1*d2+0.5*b2*d1-0.5*e1*b2+0.5*e1*d2-0.5*e2*d1)-fabs(-0.5
                                                                                                                                                                                                                                      *e2*d1+0.5*d1*a2-0.5*d2*a1+0.5*e1*d2-0.5*a2*e1+0.5*a1*e2));
                    differ_acd = fabs(-fabs(-0.5*d2*a1+0.5*c2*a1-0.5*c1*a2+0.5*d1*a2
                                            -0.5*c2*d1+0.5*c1*d2)+fabs(-0.5*e2*c1+0.5*c1*a2-0.5*c2*a1+0.5*e1*c2-0.5*a2*e1+
                                                                       0.5*a1*e2)+fabs(0.5*e2*c1-0.5*c1*d2+0.5*c2*d1-0.5*e1*c2+0.5*e1*d2-0.5*e2*d1)+
                                      fabs(-0.5*e2*d1+0.5*d1*a2-0.5*d2*a1+0.5*e1*d2-0.5*a2*e1+0.5*a1*e2));
                    differ_bcd = fabs(-fabs(0.5*b1*c2-0.5*b1*d2+0.5*b2*d1-0.5*b2*c1+0.5
                                                                                    *c1*d2-0.5*c2*d1)+fabs(0.5*e2*b1-0.5*b1*c2+0.5*b2*c1-0.5*e1*b2+0.5*e1*c2-0.5*e2
                                                                                                                                                             *c1)+fabs(0.5*e2*c1-0.5*c1*d2+0.5*c2*d1-0.5*e1*c2+0.5*e1*d2-0.5*e2*d1)+fabs(0.5
                                                                                                                                                                                                                                         *e2*b1-0.5*b1*d2+0.5*b2*d1-0.5*e1*b2+0.5*e1*d2-0.5*e2*d1));
                    if ((differ_abc < differ_min)
                         && (0.1E-4 < fabs(-0.5*a1*b2+0.5*c2*a1-0.5*c1*a2+0.5*a2*b1-0.5*b1*c2+0.5*b2*c1))
                         && ((wideGridIh(-1+i, -1+j) == 3)
                             || (wideGridIh(-1+i, j) == 3)
                             || (wideGridIh(i, j) == 3))) {
                        ihmni(-1+m, -1+n) = i;
                        ihmnj(-1+m, -1+n) = j;
                        ihmnkey(-1+m, -1+n) = abc;
                        differ_min = differ_abc;
                        if (differ_min < 0.1E-6) {
                            differ_min = 0;
                        }
                    }
                    if ((differ_acd < differ_min)
                        && (0.1E-4 < fabs(-0.5*d2*a1+0.5*c2*a1-0.5*c1*a2+0.5*d1*a2-0.5*c2*d1+0.5*c1*d2))
                        && ((wideGridIh(-1+i, -1+j) == 3)
                            || (wideGridIh(i, j) == 3)
                            || (wideGridIh(i, -1+j) == 3))) {
                        ihmni(-1+m, -1+n) = i;
                        ihmnj(-1+m, -1+n) = j;
                        ihmnkey(-1+m, -1+n) = acd;
                        differ_min = differ_acd;
                        if (differ_min < 0.1E-6) {
                            differ_min = 0;
                        }
                    }
                    if ((differ_bcd < differ_min)
                        && (0.1E-4 < fabs(0.5*b1*c2-0.5*b1*d2+0.5*b2*d1-0.5*b2*c1+0.5*c1*d2-0.5*c2*d1))
                        && ((wideGridIh(i, -1+j) == 3)
                            || (wideGridIh(i, j) == 3))) {
                        ihmni(-1+m, -1+n) = i;
                        ihmnj(-1+m, -1+n) = j;
                        ihmnkey(-1+m, -1+n) = bcd;
                        differ_min = differ_bcd;
                        if (differ_min < 0.1E-6) {
                            differ_min = 0;
                        }
                    }
                    if ((differ_abd < differ_min)
                        && (0.1E-4 < fabs(0.5*a1*b2-0.5*d2*a1+0.5*d1*a2-0.5*a2*b1+0.5*b1*d2-0.5*b2*d1))
                        && ((wideGridIh(-1+i, -1+j) == 3)
                            || (wideGridIh(-1+i, j) == 3)
                            || (wideGridIh(i, -1+j) == 3))) {
                        ihmni(-1+m, -1+n) = i;
                        ihmnj(-1+m, -1+n) = j;
                        ihmnkey(-1+m, -1+n) = abd;
                        differ_min = differ_abd;
                        if (differ_min < 0.1E-4) {
                            differ_min = 0;
                        }
                    }

                    if (differ_min <= 0.1E-4) {
                        breake = 1;
                        ibegin = (-2 + i > 1) ? -2+i : 1;
                        iend = (2 + i <= n1) ? 2+i : n1;
                        jbegin = (-2 + j > 1) ? -2+j : 1;
                        jend = (2 + j<=n2) ? 2+j : n2;
                    }
                }
        }
    }

    for(m = 1;m <= N1;m++)
        for(n = 1;n <= N2;n++)
        {
            if (fabs(denseGridIh(-1+m, -1+n)) < 1e-14) {
                fValues(-1+m, -1+n) = 0;
                continue;
            } else {
                i = (int) ihmni(-1+m, -1+n);
                j = (int) ihmnj(-1+m, -1+n);
                key = (int) ihmnkey(-1+m, -1+n);
                e1 = denseGridX(-1+m, -1+n);
                e2 = denseGridY(-1+m, -1+n);
                if( key == abc )
                {
                    a1 = wideGridX(-1+i, j-1);
                    a2 = wideGridY(-1+i, j-1);
                    b1 = wideGridX(-1+i, j);
                    b2 = wideGridY(-1+i, j);
                    c1 = wideGridX(i, j);
                    c2 = wideGridY(i, j);
                    Ha = f(a1, a2);
                    Hb = f(b1, b2);
                    Hc = f(c1, c2);
                }
                if( key == acd )
                {
                    a1 = wideGridX(-1+i, j-1);
                    a2 = wideGridY(-1+i, j-1);
                    b1 = wideGridX(i, j);
                    b2 = wideGridY(i, j);
                    c1 = wideGridX(i, j-1);
                    c2 = wideGridY(i, j-1);
                    Ha = f(a1, a2);
                    Hb = f(b1, b2);
                    Hc = f(c1, c2);
                }
                if( key == abd )
                {
                    a1 = wideGridX(-1+i, j-1);
                    a2 = wideGridY(-1+i, j-1);
                    b1 = wideGridX(i-1, j);
                    b2 = wideGridY(i-1, j);
                    c1 = wideGridX(i, j-1);
                    c2 = wideGridY(i, j-1);
                    Ha = f(a1, a2);
                    Hb = f(b1, b2);
                    Hc = f(c1, c2);
                }
                if( key == bcd )
                {
                    a1 = wideGridX(-1+i, j);
                    a2 = wideGridY(-1+i, j);
                    b1 = wideGridX(i, j);
                    b2 = wideGridY(i, j);
                    c1 = wideGridX(i, j-1);
                    c2 = wideGridY(i, j-1);
                    Ha = f(a1, a2);
                    Hb = f(b1, b2);
                    Hc = f(c1, c2);
                }

                S2abc = fabs(-b1*c2+c2*a1-a1*b2+b2*c1-c1*a2+a2*b1);
                alpha = (b2*Hc-a2*Hc+a2*Hb+Ha*c2-Ha*b2-Hb*c2)/(-b1*c2+c2*a1-a1*b2+
                                                               b2*c1-c1*a2+a2*b1);
                Gamma = (Hb*c2*a1-a2*c1*Hb+a2*b1*Hc-Ha*b1*c2-b2*a1*Hc+b2*c1*Ha)/(-b1*c2+c2*a1-a1*b2+b2*c1-c1*a2+a2*b1);
                beta = -1/(-b1*c2+c2*a1-a1*b2+b2*c1-c1*a2+a2*b1)*(b1*Hc-c1*Hb-a1*Hc
                                                                  +a1*Hb+c1*Ha-b1*Ha);

                fValues(m-1, n-1) = alpha*denseGridX(m-1, n-1)+beta*denseGridY(m-1, n-1)+Gamma;
            }
        }
    return fValues;
}

double Interpolation::interpolationElement(Point a, Point b, Point c, Point d,
                                           double fA, double fB, double fC, double fD) {
    auto Copt = -c.x*d.y/2.0+c.x*b.y/2.0+a.x*d.y/2.0-a.x*b.y/2.0
                -d.x*a.y/2.0+d.x*c.y/2.0+b.x*a.y/2.0-b.x*c.y/2.0;
    return Copt * fInterp(a, b, c, d, fA, fB, fC, fD);
}

double Interpolation::fInterp(Point a, Point b, Point c, Point d,
                              double fA, double fB, double fC, double fD) {
    auto X1 = (a.x*a.x*b.y-a.x*a.x*d.y-a.x*a.y*b.x+a.x*a.y*d.x+a.x*b.x*b.y-a.x*d.x*d.y-a.y*b.x*b.x+a.y*d.x*
               d.x+b.x*b.x*c.y-b.x*b.y*c.x+b.x*c.x*c.y-b.y*c.x*c.x+c.x*c.x*d.y-c.x*c.y*d.x+c.x*d.x*d.y-c.y*d.x*d.x)
              /(a.x*b.y-a.x*d.y-a.y*b.x+a.y*d.x+b.x*c.y-b.y*c.x+c.x*d.y-c.y*d.x)/3.0;
    auto X2 = (a.x*a.y*b.y-a.x*a.y*d.y+a.x*b.y*b.y-a.x*d.y*d.y-a.y*a.y*b.x+a.y*a.y*d.x-a.y*b.x*b.y+a.y*d.x*
               d.y+b.x*b.y*c.y+b.x*c.y*c.y-b.y*b.y*c.x-b.y*c.x*c.y+c.x*c.y*d.y+c.x*d.y*d.y-c.y*c.y*d.x-c.y*d.x*d.y)
              /(a.x*b.y-a.x*d.y-a.y*b.x+a.y*d.x+b.x*c.y-b.y*c.x+c.x*d.y-c.y*d.x)/3.0;
    return((((((fC-fD)*b.y+(-fB+fD)*c.y+d.y*(fB-fC))*a.x+((-fC+fD)*b.x+(fB-fD)*c.x
           -d.x*(fB-fC))*a.y+((fA-fD)*c.y-d.y*(fA-fC))*b.x+((-fA+fD)*c.x+d.x*(fA-fC))*b.y
           +(c.x*d.y-c.y*d.x)*(fA-fB))*X2+a.y*((-fC+fD)*b.y+(fB-fD)*c.y-d.y*(fB-fC))*a.x+((fC-fD)*b.y*b.x
           -(fB-fD)*c.y*c.x+d.x*d.y*(fB-fC))*a.y-((fA-fD)*c.y-d.y*(fA-fC))*b.y*b.x+((fA-fD)*c.y*c.x
           -d.x*d.y*(fA-fC))*b.y-c.y*d.y*(c.x-d.x)*(fA-fB))*X1+((((fC-fD)*b.x+(-fB+fD)*c.x+d.x*(fB-fC))*a.y
           -(fC-fD)*b.y*b.x+(fB-fD)*c.y*c.x-d.x*d.y*(fB-fC))*a.x+(((fA-fD)*c.x-d.x*(fA-fC))*b.y
           -(fA-fD)*c.y*c.x+d.x*d.y*(fA-fC))*b.x+c.x*d.x*(-d.y+c.y)*(fA-fB))*X2+(((-fC*d.y+fD*c.y)*b.x
           +(fC*d.x-fD*c.x)*b.y+fB*(c.x*d.y-c.y*d.x))*a.y+b.y*(fC*d.y-fD*c.y)*b.x+(-fC*d.x*d.y+fD*c.x*c.y)*b.y
           -fB*c.y*d.y*(c.x-d.x))*a.x+(((-fC*d.x+fD*c.x)*b.y+fC*d.x*d.y-fD*c.x*c.y)*b.x+fB*c.x*d.x*(-d.y+c.y))*a.y
           -(((c.x*d.y-c.y*d.x)*b.y-c.y*d.y*(c.x-d.x))*b.x+b.y*c.x*d.x*(-d.y+c.y))*fA)/((((-d.y+c.y)*b.x
           +(-c.x+d.x)*b.y+c.x*d.y-c.y*d.x)*a.y-b.y*(-d.y+c.y)*b.x+(c.x*c.y-d.x*d.y)*b.y-c.y*d.y*(c.x-d.x))*a.x
           +(((c.x-d.x)*b.y-c.x*c.y+d.x*d.y)*b.x+c.x*d.x*(-d.y+c.y))*a.y+((-c.x*d.y+c.y*d.x)*b.y+c.y*d.y*(c.x-d.x))*b.x-b.y*c.x*d.x*(-d.y+c.y)));
}

double Interpolation::calculateIntegral(double (*f)(double, double), Matrix &matrixX, Matrix &matrixY, Matrix &ih) {
    double integral = 0;
    for (int i = 0; i < matrixX.getShape(0)-1; i++) {
        for (int j = 0; j < matrixX.getShape(1)-1; j++) {
            if (ih(i, j) < 2) {
                continue;
            }
            integral += quadraticElement(f, Point(matrixX(i, j), matrixY(i, j)),
                                         Point(matrixX(i, j+1), matrixY(i, j+1)),
                                         Point(matrixX(i+1, j+1), matrixY(i+1, j+1)),
                                         Point(matrixX(i+1, j), matrixY(i+1, j)));
        }
    }
    return integral;
}

double Interpolation::calculateIntegral(Matrix& fValues, Matrix &matrixX, Matrix &matrixY, Matrix &ih) {
    double integral = 0;
    for (int i = 0; i < matrixX.getShape(0)-1; i++) {
        for (int j = 0; j < matrixX.getShape(1)-1; j++) {
            if (ih(i, j) < 2) {
                continue;
            }
            integral += interpolationElement(Point(matrixX(i, j), matrixY(i, j)),
                                             Point(matrixX(i, j+1), matrixY(i, j+1)),
                                             Point(matrixX(i+1, j+1), matrixY(i+1, j+1)),
                                             Point(matrixX(i+1, j), matrixY(i+1, j)),
                                             fValues(i, j), fValues(i, j+1),
                                             fValues(i+1, j+1), fValues(i+1, j));
        }
    }
    return integral;
}
