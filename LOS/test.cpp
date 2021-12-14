#include <iostream>
#include "SimplexMethod.h"

using std::cout;
using std::endl;

int main() {
	// 2 +
	Vec c{ 1, 0, 0, 0 };
	Vec b{ -2, 2 };
	Matrix A{ 
		{-1, 1, 1, 0}, 
		{1, 1, 0, -1}
	};

	// -2 +
	Vec c1{ 0, -1, 0, 0 };
	Vec b1{ 20, 15 };
	Matrix A1{
		{2, 5, -1, 0},
		{1, 5, 0, 1}
	};

	// 0 +
	Vec c2{ -1, 1, 0, 0 };
	Vec b2{ 0, 1 };
	Matrix A2{
		{1, -1, 1, 0},
		{0, 1, 0, 1}
	};

	// -13 +
	Vec c3{ -2, -1, 0, 0, 0 };
	Vec b3{ 32, 17, 5 };
	Matrix A3{
		{3, 4, 1, 0, 0},
		{3, 1, 0, 1, 0},
		{1, 0, 0, 0, 1}
	};

	// -4 +
	Vec c4{ 0, -2, 1, -1 };
	Vec b4{ 6, 6 };
	Matrix A4{
		{2, 1, 1, 0},
		{1, 2, 0, 1}
	};

	// -2 (?)
	Vec c5{ 0, 0, 0, 0, 200, 175, -1100, -2 };
	Vec b5{ 0, 0, 0, 1 };
	Matrix A5{
		{0, 1, 0, 0, -3, -5/4, 7, 1/50},
		{1, 0, 0, 0, 1/3, 1/6, -1, -1/150},
		{0, 0, 1, 0, 75/2, -25/4, 175/2, 1/4},
		{0, 0, 0, 1, 0, 0, 0, 1}
	};

	// 1 +
	Vec c6{ 1, 1, 0, 5 };
	Vec b6{ 1, 1, 5 };
	Matrix A6{
		{1, 1, -1, 3},
		{1, -2, 3, -1},
		{5, -4, 7, 3}
	};

	// -22/7 +
	Vec c7{ 1, -2, -1, -1 };
	Vec b7{ 0, 2 };
	Matrix A7{
		{3, -3, 4, 2},
		{1, 1, 1, 3}
	};

	SimplexMethod sm(c5, b5, A5);
	sm.solve();
	cout << (int)sm.get_state() << endl;
	cout << sm.get_solution() << endl;
	return 0;
}