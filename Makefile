CC=g++
ADD=-lX11 -std=c++17

main: clean R2G Point Matr GW Inter Grid
	$(CC) -o main main.cpp R2Graph.o Matrix.o Point.o gwindow.o DrawGridWindow.o Interpolation.o $(ADD) -lm


R2G: R2Graph.cpp
	$(CC) -c R2Graph.cpp $(ADD)

Matr: Matrix.cpp
	$(CC) -c Matrix.cpp $(ADD)

Inter: Interpolation.cpp
	$(CC) -c Interpolation.cpp $(ADD)

Point: Point.cpp
	$(CC) -c Point.cpp $(ADD)

Grid: DrawGridWindow.cpp
	$(CC) -c DrawGridWindow.cpp $(ADD)

GW: gwindow.cpp
	$(CC) -c gwindow.cpp $(ADD)

clean:
	rm -f *.o
