#include "SimplexMethod.h"

int main() {
	Vec c{ -1, -1, 0, -5 };
	Vec b{ 1, 1, 5 };
	Matrix A{ {1, 1, -1, 3}, 
			  {1, -2, 3, -1}, 
			  {5, -4, 7, 3} 
			};
	SimplexMethod sm(c, b, A);
	sm.solve();
	return 0;
}