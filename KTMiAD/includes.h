#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <sstream>

const double BIG_VALUE = 1e15;

using namespace std;

template <typename T>
std::string toString(T val)
{
	std::ostringstream oss;
	oss << val;
	return oss.str();
}

typedef // создаем новый прототип (в данном случае указатель на функцию)
double // возвращаемое значение (такое же как в функциях)
(*basis_func) // имя прототипа (в коде употребляется без звездочки)
(double); // список параметров (такое же как в функциях)

typedef // создаем новый прототип (в данном случае указатель на функцию)
double // возвращаемое значение (такое же как в функциях)
(*func_3coord_to_1) // имя прототипа (в коде употребляется без звездочки)
(double, double, double); // список параметров (такое же как в функциях)

typedef // создаем новый прототип (в данном случае указатель на функцию)
double // возвращаемое значение (такое же как в функциях)
(*func_2coord_to_1) // имя прототипа (в коде употребляется без звездочки)
(double, double); // список параметров (такое же как в функциях)

// соответствие локальных узлов/производных 2Д элементам-граням
const int face_to_node[6][4] = { {0, 2, 4, 6}, { 1,3,5,7 }, { 0,1,4,5 }, { 2,3,6,7 }, { 0,1,2,3 }, { 4,5,6,7 } };
const int face_to_derivative[6][2] = { {1, 2}, { 1, 2 }, { 0, 2 }, { 0, 2 }, { 0, 1 }, { 0, 1 } };

//---------------------------линейные одномерные-------------------------------------------
inline double linear1(double ksi) // от -1 до 1
{
	return (1 - ksi) / 2;
}
inline double linear2(double ksi) // от -1 до 1
{
	return (ksi + 1) / 2;
}
inline double qudratic(double ksi) // от -1 до 1
{
	return 1.0 - ksi * ksi;
}

inline double dlinear1_dksi(double ksi)
{
	return -1. / 2;
}
inline double dlinear2_dksi(double ksi)
{
	return 1. / 2;
}
//--------------------------------------------------------------------------------


//---------------одномерные квадратичные функции, где ksi = (x - x_(2k-1)) / (x_(2k+1) - x_(2k-1)))-------------------------------------------------------
inline double quadratic1(double ksi)
{
	return 2 * (ksi - 0.5) * (ksi - 1);
}
inline double quadratic2(double ksi)
{
	return -4 * ksi * (ksi - 1);
}
inline double quadratic3(double ksi)
{
	return 2 * (ksi - 0.5) * ksi;
}
inline double getKsi(double x, double x1, double x2)
{
	return 2.0 * (x - x1) / (x2 - x1) - 1.0;
}
// c. 331
inline int mu(int i)
{
	return (i % 2);
}
inline int nu(int i)
{
	return (i / 2) % 2;
}
inline int teta(int i)
{
	return (i / 4) % 2;
}
//--------------------------------------------------------------------------------

//--------------------------test--------------------------------------------------
inline double f_test(double x, double y, double z)
{
	return// 0.4*(5.0 + 0.2 * x + y + 30.0 * z + 0.5 * x * y + x * z + 10.0 * y * z + x * y * z)
		//2*cos(x + z)
		//x + y + z
		//- 6 + x * x + y * y + z * z
		// -6*x-6*y-6*z + x * x * x + y * y * y + z * z * z
		-12 * (x * x + y * y + z * z) + x * x * x * x + y * y * y * y + z * z * z * z
		;
}

inline double u_test(double x, double y, double z)
{
	return //5.0 + 0.2 * x + y + 30.0 * z + 0.5 * x * y + x * z + 10.0 * y * z + x * y * z
		// cos(x+z)
		//cos(x + y + z)
		//x + y + z
		//x * x + y * y + z * z
		// x * x * x + y * y * y + z * z * z
		x * x * x * x + y * y * y * y + z * z * z * z

		;
}
//--------------------------------------------------------------------------------

struct Node3D
{
	double x;
	double y;
	double z;

	Node3D() { x = 0; y = 0; z = 0; }
	Node3D(double _x, double _y, double _z) { x = _x; y = _y; z = _z; }
};

struct Node2D
{
	double x;
	double y;

	Node2D() { x = 0; y = 0; }
	Node2D(double _x, double _y) { x = _x; y = _y; }
};

const vector<vector<double>> G_linear = {
	{1, -1},
	{-1, 1} }; // * 1/h
const vector<vector<double>> M_linear = {
	{1.0 / 3.0, 1.0 / 6.0},
	{1.0 / 6.0, 1.0 / 3.0} }; // * h


