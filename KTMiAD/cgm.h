#pragma once
#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

int Solve_Gauss(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double>& q);

class cgm_solver
{
private:
	int n = 0, total_iter = 0, max_iter = 0;
	double eps = 0., pr_norm = 0.;

	std::vector<int> ig, jg;
	std::vector<double> ggu, & ggl = ggu, diag, pr;

	std::vector<double> r, z, rr_vec, az, s, Ll, Ld;


	void llt();

	void matrix_dot_vector(const std::vector<double>& x, std::vector<double>& y);

	void matrix_dot_vector_llt(const std::vector<double>& x, std::vector<double>& y);

	void vec_diff(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& res);

	double scalar_product(const std::vector<double>& x, const std::vector<double>& y);

	double norm(const std::vector<double>& x);

	void solve_auxiliary_system(const std::vector<double>& f, std::vector<double>& x);

	double relative_residual(const std::vector<double>& x);

public:
	cgm_solver(const std::string& path, const bool& ind_from_zero = true);
	cgm_solver(std::vector<int>& ia, std::vector<int> &ja,
		std::vector<double> &di, std::vector<double> &al, std::vector<double> &b, double _eps = 1e-35);

	int get_n();

	int get_total_iter();

	void no_preconditioning(std::vector<double>& x0);

	void diag_preconditioning(std::vector<double>& x0);

	void llt_preconditioning(std::vector<double>& x0);

	void reset();

	void clear_all();

};
