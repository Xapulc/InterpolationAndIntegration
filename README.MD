**Описание**

Это семестровый проект практикума на ЭВМ 4 курса 1 семестра
мехмата МГУ.
Написан на C++, использует визуализацию X11.

Программа [main](main.cpp) строит редкую и густую сетки 
на заданной подобласти плоскости с дыркой,
вычисчисляет функции на сетках, 
производит интерполяцию с редкой сетки на густую,
вычисляет интегралы,
визуализирует сетку и изолинии.

**Как это запустить**

Инструкция написана для пользователей Linux (Kubuntu в моём случае).
Возможно, на MacOs тоже будет работать, 
а вот на Windows придётся переписать визуализацию.

Перед запуском программы вам стоит убедиться,
что у вас стоит пакет для визуализации через X11 сервер.
Его можно установить командой
```
sudo apt install libx11-dev
```

В [Makefile](Makefile) содержатся команды для сборки проекта,
поэтому сначала стоит в корне проекта выполнить команду

```
make
```

У вас тем самым создастся файл исполняемый `main`.
Он принимает 8 обязательных аргументов:
* `wideGridN1` - размер по оси `x` редкой решётки;
* `wideGridN2` - размер по оси `y` редкой решётки;
* `wideGridmaxIter` - максимальное количество итераций построения редкой сетки;
* `wideGridomega` - значения параметра `omega` функции построения редкой сетки;
* `denseGridN1` - размер по оси `x` густой решётки;
* `denseGridN2` - размер по оси `y` густой решётки;
* `denseGridmaxIter` - максимальное количество итераций построения густой сетки;
* `denseGridomega` - значения параметра `omega` функции построения густой сетки.

Приемлемый результат получается с параметрами
```
./main 10 10 0.5 60 250 250 0.8 300
```

**Детали**

Внешний контур области задаётся дугами окружностей в функции [contour](Interpolation.cpp).
Эти дуги вычисляются функцией [fillArc](Interpolation.cpp)
и в матрицах `matrixX` и `matrixY` заполняется часть столбца
или часть строки координатами этих дуг.

Сетка строится функцией [buildGrid](Interpolation.cpp).
Матрица `ih` выступает в роли индикатора границ графика, дырки, внутренней части.
С границ графика, определённых [contour](Interpolation.cpp), 
сетка строится линейной интерполяцией (учитывая дырку).
Затем итерационным методом [Laplasiation](Interpolation.cpp)
уточняется структура сетки.
На выходе в матрицах `matrixX` и `matrixY`
мы получаем координаты сетки,
а в `ih` - параметры структуры сетки.

По построенной сетке можно вычислить интеграл 
двухпараметрической вещественной функции `f`.
Функция [calculateIntegral](Interpolation.cpp) вычисляет интеграл 
функции `f` на сетке
с помощью квадратур [quadraticElement](Interpolation.cpp).

Если имеется редкая и густая сетки, можно интерполировать значения функции `f`
с редкой сетки на густую, что делается функцией [gridInterpolation](Interpolation.cpp).
Интеграл от значений проинтерполированной функции на густой сетке 
можно посчитать всё той же функцией [calculateIntegral](Interpolation.cpp).
Только последняя будет использовать 
интерполяционные элементы [interpolationElement](Interpolation.cpp)
для вычисления интеграла.

Итак, у нас есть сетки и значения функций.
Их визуализацией занимается класс [DrawGridWindow](DrawGridWindow.cpp).
После задания решётки функцией [setGrid](DrawGridWindow.cpp)
и значений функции [setFunction](DrawGridWindow.cpp)
можно строить решётку в X11 сервере функцией [draw](DrawGridWindow.cpp).
В результате построится сетка,
раскрашенная по изолиниям значений функции.

**Ссылки**
* Проект реализован по [материалам](https://drive.google.com/drive/folders/1Nxj7Sdg6P2Fannckt0P748msy5j57VYe?usp=sharing) 
[Соколова Александра Германовича](http://numa.math.msu.su/staff/index.html).
* Класс матриц взят из проекта [MatrixSolvers](https://github.com/Xapulc/MatrixSolvers).
* Визуализация выполнена обёрткой над проектом [Gwindow](http://mech.math.msu.su/~nap/2/GWindow/),
как в проекте [RotationIntersectionOfCube](https://github.com/Xapulc/RotationIntersectionOfCube).
