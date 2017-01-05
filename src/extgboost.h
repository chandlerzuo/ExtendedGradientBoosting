#pragma once
#include "tree.h"

class Extgboost {
public:
	static double prop_tune;
	static int n_itr;
	static double alpha;
	static int n_comp;
	std::vector<double> y_test_pred;

	void LoadData(std::vector<std::vector<double>> &, std::vector<double> &, std::vector<std::vector<double>> &, std::vector<bool> &);
	void Fit();
	void Clean();

private:
	std::vector<Node *> vec_h;
	std::vector<Node *> vec_g;
	std::vector<double> y_train_est;
	std::vector<double> y_tune_est;

	double loss_tune;
	double loss_train;

	std::vector<bool> is_numerics;
	std::vector<std::vector<double>> x_train;
	std::vector<std::vector<double>> x_tune;
	std::vector<std::vector<double>> * x_test;
	std::vector<std::vector<double>> x_train_t;
	std::vector<std::vector<double>> x_tune_t;
	std::vector<std::vector<double>> x_test_t;
	std::vector<double> y_train;
	std::vector<double> y_tune;

	void ComputeLoss();
};