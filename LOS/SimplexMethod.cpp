#include "SimplexMethod.h"

SimplexMethod::SimplexMethod(const Vec& c, const Vec& b, const Matrix& A) :
	_n(c.size()), 
	_m(b.size()),
	_state(Solution::notfound) {

	// making the simplex table
	_simplex_table.resize(_m + 1);
	for (auto& row : _simplex_table) row.resize(_n + 1);

	_simplex_table[0][0] = 0.0;
	for (int j = 1; j < _n + 1; ++j) _simplex_table[0][j] = -c[j - 1];
	for (int i = 1; i < _m + 1; ++i) _simplex_table[i][0] = b[i - 1];
	for (int i = 1; i < _m + 1; ++i)
		for (int j = 1; j < _n + 1; ++j)
			_simplex_table[i][j] = A[i - 1][j - 1];

	// basis is not detected
	_basis.resize(_m);
	for (auto& basis_col_number : _basis) basis_col_number = -1;
}

int SimplexMethod::detect_leading_column() {
	int leading_column{ -1 };
	double max{ 0.0 };
	
	for (int j = 1; j < _n + 1; ++j)
		if (_simplex_table[0][j] > max) {
			max = _simplex_table[0][j];
			leading_column = j;
		}

	if (leading_column > 0) return leading_column;
	else return -1;
}

int SimplexMethod::detect_leading_row(const int& leading_column) {
	int leading_row{ -1 };
	double ratio{ std::numeric_limits<double>::max() };
	double temp_ratio = 0.0;

	for (int i = 1; i < _m + 1; ++i) {
		if (_simplex_table[i][leading_column] != 0.0) {
			temp_ratio = _simplex_table[i][0] / _simplex_table[i][leading_column];
			if ((temp_ratio < ratio) && (temp_ratio > 0)) {
				ratio = temp_ratio;
				leading_row = i;
			}
		}
	}

	if (leading_row > 0) return leading_row;
	else return -1;
}

// lr means leading row, lc means leading column
void SimplexMethod::basis_transform(const int& lr, const int& lc) {
	bool column_detected = false;
	for (int i = 0; ((i < _basis.size()) && (!column_detected)); ++i)
		if (_simplex_table[lr][_basis[i]] == 1) {
			_basis[i] = lc;
			column_detected = true;
		}
}

// lr means leading row, lc means leading column
void SimplexMethod::gauss_transform(const int& lr, const int& lc) {
	for (int j = 0; j < _n + 1; ++j) 
		_simplex_table[lr][j] / _simplex_table[lr][lc];

	for (int i = 0; i < lr; ++i) {
		double multiplicator = _simplex_table[i][lc];
		for (int j = 0; j < _n + 1; ++j)
			_simplex_table[i][j] -= _simplex_table[lr][j] * multiplicator;
	}

	for (int i = lr + 1; i < _m + 1; ++i) {
		double multiplicator = _simplex_table[i][lc];
		for (int j = 0; j < _n + 1; ++j)
			_simplex_table[i][j] -= _simplex_table[lr][j] * multiplicator;
	}
}

void SimplexMethod::iterate() {
	int leading_column = detect_leading_column();
	if (leading_column == -1) {
		_state = Solution::optimal; 
		return;
	}

	int leading_row = detect_leading_row(leading_column);
	if (leading_row == -1) {
		_state = Solution::unlimited;
		return;
	}

	//basis_transform(leading_row, leading_column);
	gauss_transform(leading_row, leading_column);
}

void SimplexMethod::solve() {
	while (_state == Solution::notfound) iterate();
}
