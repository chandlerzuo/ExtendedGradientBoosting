#pragma once

#include <vector>
#include "tree.h"

class Gradient {
public:
	static double min_prob;
	int J;
	std::vector<std::vector<double>> h;
	std::vector<std::vector<double>> g;
	std::vector<double> f;
	std::vector<double> * y;
	std::vector<std::vector<double>> decent_h;
	std::vector<std::vector<double>> decent_g;
	int type; // 0 for binary, 1 for numeric, 2 for multi-class
	double loglik;

	Gradient(int);
	Gradient(int, int);
	void LoadData(std::vector<double> & y);
	void Initialize();
	void Initialize(double);
	void Initialize(double, int);
	void ComputeGradient_g(int);
	void ComputeGradient_h(int);
	void Update_g(int, std::vector<double> &);
	void Update_h(int, std::vector<double> &);
	void Update_g(int, std::vector<double> &, double);
	void Update_h(int, std::vector<double> &, double);
	void RevUpdate_g(int, std::vector<double> &, double);
	void RevUpdate_h(int, std::vector<double> &, double);
	void UpdateLikelihood();
	double UpdateTreeFit(Node *, std::vector<std::vector<double>> &, double, int, bool, bool);
	double ComputeStepSize(std::vector<double> &, int, bool);
	double trim(double);
	double LossWithUpdate(int, std::vector<double> &, double, bool);

private:
	std::vector<std::vector<double>> d;
	std::vector<std::vector<double>> p;
	std::vector<std::vector<double>> dh;
	std::vector<std::vector<double>> d2h;
	std::vector<std::vector<double>> dg;
	std::vector<std::vector<double>> d2g;
	std::vector<double> df;
	std::vector<double> d2f;
	std::vector<std::vector<double>> fdh;
	std::vector<std::vector<double>> fd2h;
	std::vector<std::vector<double>> fdg;
	std::vector<std::vector<double>> fd2g;

	double dLossWithUpdate(int, std::vector<double> &, double, bool);
};