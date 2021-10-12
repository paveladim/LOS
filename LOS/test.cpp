#include "SimplexMethod.h"

int main() {
	Vec c{ -2, -1, 0, 0, 0 };
	Vec b{ 32, 17, 5 };
	Matrix A{ {3, 4, 1, 0, 0}, {3, 1, 0, 1, 0}, {1, 0, 0, 0, 1} };
	SimplexMethod sm(c, b, A);
	sm.solve();
	return 0;
}