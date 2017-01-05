#include "extgboost.h"
#include "utility.h"
#include "gradient.h"
#include <iostream>

void Extgboost::LoadData(std::vector<std::vector<double>> & x, std::vector<double> & y, std::vector<std::vector<double>> & x_test, std::vector<bool> & is_numerics) {
	this->is_numerics = is_numerics;
	this->x_test = &x_test;
	this->x_test_t = Transpose(x_test);

	std::vector<int> tune_ids = Sample(x.size(), Extgboost::prop_tune);
	std::map<int, bool> tune_id_set;
	for (int i = 0; i < tune_ids.size(); ++i) {
		tune_id_set.insert(std::pair<int, bool>(tune_ids[i], true));
	}
	for (int i = 0; i < x.size(); ++i) {
		if (tune_id_set.find(i) != tune_id_set.end()) {
			this->x_tune.push_back(x[i]);
			this->y_tune.push_back(y[i]);
		}
		else {
			this->x_train.push_back(x[i]);
			this->y_train.push_back(y[i]);
		}
	}
	this->x_tune_t = Transpose(this->x_tune);
	this->x_train_t = Transpose(this->x_train);
}

void Extgboost::Clean() {
	for (int i = 0; i < this->vec_h.size(); ++i) {
		delete this->vec_h[i];
	}
	for (int i = 0; i < this->vec_g.size(); ++i) {
		delete this->vec_g[i];
	}
	this->vec_h.clear();
	this->vec_g.clear();

	this->y_test_pred.clear();
	this->x_test = nullptr;
	this->x_test_t.clear();
	
	this->y_train.clear();
	this->y_train_est.clear();
	this->x_train.clear();
	this->x_train_t.clear();
	
	this->y_tune.clear();
	this->y_tune_est.clear();
	this->x_tune.clear();
	this->x_tune_t.clear();
}

void Extgboost::Fit() {
	Gradient grad = Gradient(Extgboost::n_comp);
	Gradient grad_tune = Gradient(Extgboost::n_comp);
	Gradient grad_test = Gradient(Extgboost::n_comp);
	grad.LoadData(this->y_train);
	grad.Initialize();
	grad.UpdateLikelihood();
	std::cout << "Initial loglik = " << grad.loglik << "\n";
	grad_tune.Initialize(grad.f[0], this->y_tune.size());
	grad_test.Initialize(grad.f[0], this->x_test->size());

	for (int i = 0; i < this->y_tune.size(); ++i) {
		this->y_tune_est.push_back(grad.f[0]);
	}
	this->y_train_est = grad.f;

	double alpha = Extgboost::alpha;
	double step_size = 1;
	for (int itr = 0; itr < Extgboost::n_itr; ++itr) {
		for (int j = 0; j < grad.J; ++j) {
			grad.ComputeGradient_g(j);

			/*
			grad.Update_g(j, grad.decent_g[j]);
			grad.UpdateLikelihood();
			double loglik1 = grad.loglik;
			//std::cout << "Loglik = " << loglik1 << "\n";

			
			grad.RevUpdate_g(j, grad.decent_g[j], 1);
			grad.Update_g(j, grad.decent_g[j], alpha);
			grad.UpdateLikelihood();
			double loglik2 = grad.loglik;
			grad.RevUpdate_g(j, grad.decent_g[j], alpha);
			*/
			Node * tree = new Node();
			tree->LoadData(grad.decent_g[j], this->x_train, this->x_train_t, this->is_numerics);
			tree->Build();
			this->vec_g.push_back(tree);
			step_size = grad.UpdateTreeFit(tree, this->x_train, alpha, j, false, true);
			grad.UpdateLikelihood();
			this->y_train_est = grad.f;
			grad_tune.UpdateTreeFit(tree, this->x_tune, alpha * step_size, j, false, false);
			this->y_tune_est = grad_tune.f;
			grad_test.UpdateTreeFit(tree, * this->x_test, alpha * step_size, j, false, false);
			this->y_test_pred = grad_test.f;
			this->ComputeLoss();
			std::cout << "Iteration " << itr << " updated g, training set loss = " << this->loss_train << ", tuning set loss = " << this->loss_tune << "\n";// " target loglik = " << loglik1 << " with step " << loglik2 << " tree loss = " << tree->loss << "\n";
			
			grad.ComputeGradient_h(j);

			/*
			grad.Update_h(j, grad.decent_h[j]);
			grad.UpdateLikelihood();
			loglik1 = grad.loglik;
			//std::cout << "loglik " << loglik1 << "\n";
			
			grad.RevUpdate_h(j, grad.decent_h[j], 1);
			grad.UpdateLikelihood();
			loglik2 = grad.loglik;
			grad.Update_h(j, grad.decent_h[j], alpha);
			grad.UpdateLikelihood();
			double loglik3 = grad.loglik;
			grad.RevUpdate_h(j, grad.decent_h[j], alpha);
			grad.UpdateLikelihood();
			double loglik4 = grad.loglik;
			*/
			tree = new Node();
			tree->LoadData(grad.decent_h[j], this->x_train, this->x_train_t, this->is_numerics);
			tree->Build();
			this->vec_h.push_back(tree);

			step_size = grad.UpdateTreeFit(tree, this->x_train, alpha, j, true, true);
			grad.UpdateLikelihood();
			this->y_train_est = grad.f;
			grad_tune.UpdateTreeFit(tree, this->x_tune, alpha * step_size, j, true, false);
			this->y_tune_est = grad_tune.f;
			grad_test.UpdateTreeFit(tree, * this->x_test, alpha * step_size, j, true, false);
			this->y_test_pred = grad_test.f;
			this->ComputeLoss();
			std::cout << "Iteration " << itr << " updated h, training set loss = " << this->loss_train << ", tuning set loss = " << this->loss_tune << "\n";// " target loglik = " << loglik1 << " with step " << loglik3 << " tree loss = " << tree->loss << " previous loglik ?= " << loglik2 << " ?= " << loglik4 << "\n";
			
		}
		std::cout << "Iteration " << itr << " training set loss = " << this->loss_train << ", tuning set loss = " << this->loss_tune << "\n";
		//alpha *= Extgboost::alpha;
	}
}

void Extgboost::ComputeLoss() {
	this->loss_tune = 0;
	this->loss_train = 0;

	for (int i = 0; i < this->y_train.size(); ++i) {
		if (y_train[i] == 1) {
			this->loss_train += 1 - this->y_train_est[i];
		}
		else {
			this->loss_train += this->y_train_est[i];
		}
	}
	this->loss_train /= this->y_train.size();

	for (int i = 0; i < this->y_tune.size(); ++i) {
		if (y_tune[i] == 1) {
			this->loss_tune += 1 - this->y_tune_est[i];
		}
		else {
			this->loss_tune += this->y_tune_est[i];
		}
	}
	this->loss_tune /= this->y_tune.size();

}

