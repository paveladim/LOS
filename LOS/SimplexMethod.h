#ifndef SIMPLEX_METHOD_H
#define SIMPLEX_METHOD_H

#include <vector>

using Vec = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;
using Basis = std::vector<int>;

/*
*		this is linear optimization solver. it uses simplex method.
* 
*		it solves the following task:
*		<c, x> -> max
*		Ax = b
*		x >= 0
* 
*		make sure that you transformed all inequalities to equalities
*/

enum class Solution { optimal, unlimited, inconsistent, notfound };

class SimplexMethod {
private: // Fields
	Matrix _simplex_table;
	// size of objective vector
	int _n;
	// size of constraints
	int _m;
	// numbers of basis columns
	Basis _basis;
	// solution's state
	Solution _state;
	// marker for M-method
	bool _is_m_method;
private: // Methods
	int find_basis();
	int detect_leading_column();
	int detect_leading_row(const int& leading_column);
	void basis_transform(const int& lr, const int& lc);
	void gauss_transform(const int& lr, const int& lc);
	void iterate();
	void solve_imit(const int& quantity_imit);
	void gauss_tranform_imit();
	double get_solution() { return _simplex_table[0][0]; }
public:
	SimplexMethod(const Vec& c, const Vec& b, const Matrix& A, const bool& m = false);
	void solve();
};

#endif // SIMPLEX_METHOD_H

