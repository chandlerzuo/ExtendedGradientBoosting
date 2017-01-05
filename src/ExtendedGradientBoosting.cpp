// ExtendedGradientBoosting.cpp : Defines the entry point for the console application.
//

/*
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define DEBUG_NEW new(_NORMAL_BLOCK, THIS_FILE, __LINE__)

#ifdef _DEBUG
#ifndef DBG_NEW
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
#define new DBG_NEW
#endif
#endif  // _DEBUG
*/

#include "stdafx.h"
#include "utility.h"
#include "tree.h"
#include "csvreader.h"
#include "gradient.h"
#include "extgboost.h"
#include <iostream>

void testSubsample();
void testBuildRuleSet();
void testBuildTree();
void testGradientComputation();
void testEgb();
void testWriteCsv();

int RuleSet::max_rules = 3;
int Node::max_depth = 2;
int Node::min_leaf_size = 10;
double Node::col_subsample = 0.2;
double Node::subsample = 1;
double Gradient::min_prob = 1e-10;
int Extgboost::n_itr = 10;
int Extgboost::n_comp = 3;
double Extgboost::prop_tune = 0.5;
double Extgboost::alpha = 0.1;

int main()
{
//	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
//	_CrtMemState s1, s2, s3;
//	_CrtMemCheckpoint(&s1);
	// memory allocations take place here

	std::cout << "Start program\n";
	//testBuildRuleSet();

	//testSubsample();
	//testBuildTree();
	//testGradientComputation();

	testEgb();

	//testWriteCsv();

//	_CrtDumpMemoryLeaks();
//	_CrtMemCheckpoint(&s2);
//	if (_CrtMemDifference(&s3, &s1, &s2))
//		_CrtMemDumpStatistics(&s3);
}

void testSubsample() {
	std::vector<int> ret = Sample(100, 0.1);
	std::cout << "Subsample results: ";
	for (int i = 0; i < ret.size(); i++) {
		std::cout << ret[i] << " ";
	}
	std::cout << "\n";
	ret = Sample(20, 1);
	std::cout << "Subsample results: ";
	for (int i = 0; i < ret.size(); i++) {
		std::cout << ret[i] << " ";
	}
	std::cout << "\n";
}

void testBuildRuleSet() {
	std::vector<double> y{ 0,0,0,0,0,0,1,1,1,1,1,1 };
	std::vector<double> x1{ 1, -1, 4, 5, 6, 7, 8, 9, 9, 2, 3, 3 };
	std::vector<double> x2{ 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	std::vector<double> x3{ 11, 11, 11, 11, 1, 1, 2, 2, 3, 3, 11, 11 };
	std::vector<bool> types{ true, false, false };
	std::vector<int> members{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

	std::vector<std::vector<double>> x_t{ x1, x2, x3 };
	std::vector<std::vector<double>> x = Transpose(x_t);
	std::vector<double>* yptr = &y;
	std::vector<std::vector<double>>* xptr = &x;
	std::vector<std::vector<double>>* xtptr = &x_t;
	std::vector<int>*mbrptr = &members;
	std::vector<bool>* typeptr = &types;
	RuleSet r;
	r.LoadData(*yptr, *xptr, *xtptr, *typeptr, *mbrptr);
	r.Build();
	r.Print();
	r.Clean();

	r = RuleSet();
	y = std::vector<double>{ 0,0,0,0,0,0,1,1,1,1,1,1 };
	x1 = std::vector<double>{ 11, 11, 11, 11, 1, 1, 1, 1, 1, 1, 1, 1 };
	x2 = std::vector<double>{ 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12 };
	x3 = std::vector<double>{ 1, -1, 4, 5, 6, 7, 8, 9, 9, 2, 3, 3 };
	types = std::vector<bool>{false, true, true };
	x_t = std::vector<std::vector<double>>{ x1, x2, x3 };
	x = Transpose(x_t);
	yptr = &y;
	xptr = &x;
	xtptr = &x_t;
	typeptr = &types;
	r.LoadData(*yptr, *xptr, *xtptr, *typeptr, *mbrptr);
	r.Build();
	r.Print();
	r.Clean();

	r = RuleSet();
	y = std::vector<double>{ 0,0,0,0,0,0,1,1,1,1,1,1 };
	x1 = std::vector<double>{ 6, 5, 10, 11, 4, 2, 3, 1, 20, 15, 12, 9};
	x2 = std::vector<double>{ -2, 2, 2, -1, 0, 0, 1, 1, 1, 1, 2, -1 };
	x3 = std::vector<double>{ 0, -1, 1, 2, 3, 1, 0, -1, 1, 2, 0, -3 };
	types = std::vector<bool>{ true, false, true };
	x_t = std::vector<std::vector<double>>{ x1, x2, x3 };
	x = Transpose(x_t);
	yptr = &y;
	xptr = &x;
	xtptr = &x_t;
	typeptr = &types;
	r.LoadData(*yptr, *xptr, *xtptr, *typeptr, *mbrptr);
	r.Build();
	r.Print();
	r.Clean();

	r = RuleSet();
	y = std::vector<double>{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1 };
	x1 = std::vector<double>{ 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0 };
	x2 = std::vector<double>{ 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0 };
	x3 = std::vector<double>{ 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0 };
	members = std::vector<int> { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
	std::vector<double> x4 = std::vector<double>{ 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1 };
	types = std::vector<bool>{ true, false, true, false };
	x_t = std::vector<std::vector<double>>{ x1, x2, x3, x4 };
	x = Transpose(x_t);
	yptr = &y;
	xptr = &x;
	xtptr = &x_t;
	typeptr = &types;
	mbrptr = &members;
	r.LoadData(*yptr, *xptr, *xtptr, *typeptr, *mbrptr);
	r.Build();
	r.Print();
	r.Clean();
}

void testBuildTree() {
	std::vector<double> * y;
	std::vector<std::vector<double>> * x;
	y = new std::vector<double>();
	x = new std::vector<std::vector<double>>();
	ReadCsv("c:\\Users\\Chandler\\Kaggle\\BNP\\data\\train_data_subsample.csv", 0, y, x, true);

	std::vector<double> * y_new = new std::vector<double>();
	std::vector<std::vector<double>> * x_new = new std::vector<std::vector<double>>();
	ReadCsv("c:\\Users\\Chandler\\Kaggle\\BNP\\data\\tune_data_subsample.csv", 0, y_new, x_new, true);
	std::cout << "Loaded data.\n";

	std::vector<std::vector<double>> * x_t = new std::vector<std::vector<double>>();
	std::vector<std::vector<double>> * x_t_new = new std::vector<std::vector<double>>();

	(*x_t) = Transpose(*x);
	(*x_t_new) = Transpose(*x_new);

	std::vector<bool> is_numerics;
	for (int i = 0; i < x_t->size(); ++i) {
		is_numerics.push_back(true);
	}

	Node root;

	root.LoadData(*y, *x, *x_t, is_numerics);

	root.Build();

	root.Print();

	std::vector<double> y_pred = root.Predict(*x_new);
	double train_error = 0;
	for (int i = 0; i < y_pred.size(); i++) {
		train_error += (y_pred[i] - y_new->at(i)) * (y_pred[i] - y_new->at(i));
	}
	train_error /= y_pred.size();
	std::cout << "Tuning set error " << train_error << "\n";

	std::cout << "Loglik = " << root.loss << "\n";

	delete x, y, x_new, x_t, x_t_new;

}

void testGradientComputation() {
	std::vector<double> y{ 1,0,1,1,0,0,1,0,1 };
	int J = 4;
	Gradient grad = Gradient(J);
	grad.LoadData(y);
	grad.Initialize();

	grad.UpdateLikelihood();
	std::cout << "Initial loglik = " << grad.loglik << "\n";
	for (int k = 0; k < 10; k++) {
		for (int j = 0; j < J; j++) {
			grad.ComputeGradient_h(j);
			for (int i = 0; i < grad.y->size(); ++i) {
				//			std::cout << "g = " << grad.g[j][i] << " h = " << grad.h[j][i] << " p = " << grad.p[j][i] << " d = " << grad.d[j][i] << " f = " << grad.f[i] << " fdg = " << grad.fdg[j][i] << " fd2g = " << grad.fd2g[j][i] << " fdh = " << grad.fdh[j][i] << " fd2h = " << grad.fd2h[j][i] << " dg = " << grad.fdg[j][i] << " d2g = " << grad.d2g[j][i] << " dh = " << grad.dh[j][i] << " d2h = " << grad.d2h[j][i] << "\n";
			}
			grad.Update_h(j, grad.decent_h[j]);
			grad.UpdateLikelihood();
			std::cout << "Iteration " << k << " update h, loglik = " << grad.loglik << "\n";
			grad.ComputeGradient_g(j);
			for (int i = 0; i < grad.y->size(); ++i) {
				//std::cout << "g = " << grad.g[j][i] << " h = " << grad.h[j][i] << " p = " << grad.p[j][i] << " d = " << grad.d[j][i] << " f = " << grad.f[i] << " fdg = " << grad.fdg[j][i] << " fd2g = " << grad.fd2g[j][i] << " fdh = " << grad.fdh[j][i] << " fd2h = " << grad.fd2h[j][i] << " dg = " << grad.fdg[j][i] << " d2g = " << grad.d2g[j][i] << " dh = " << grad.dh[j][i] << " d2h = " << grad.d2h[j][i] << "\n";
			}
			grad.Update_g(j, grad.decent_g[j]);
			grad.UpdateLikelihood();
			std::cout << "Iteration " << k << " update g, loglik = " << grad.loglik << "\n";
		}
	}
	/*
	for (int i = 0; i < y.size(); i++) {
		std::cout << "y = " << y[i];
		for (int j = 0; j < J; j++) {
			std::cout << " p[" << j << "] = " << grad.p[j][i] << " d[" << j << "] = " << grad.d[j][i] << "\n";
		}
	}
	*/
}

void testEgb() {
	std::vector<double> * y;
	std::vector<std::vector<double>> * x;
	y = new std::vector<double>();
	x = new std::vector<std::vector<double>>();
	ReadCsv("c:\\Users\\Chandler\\Kaggle\\BNP\\data\\train_data_subsample.csv", 0, y, x, true);

	std::vector<double> * y_new = new std::vector<double>();
	std::vector<std::vector<double>> * x_new = new std::vector<std::vector<double>>();
	ReadCsv("c:\\Users\\Chandler\\Kaggle\\BNP\\data\\tune_data_subsample.csv", 0, y_new, x_new, true);
	std::cout << "Loaded data.\n";

	std::vector<bool> is_numerics;
	for (int i = 0; i < x->at(0).size(); ++i) {
		is_numerics.push_back(true);
	}

	Extgboost model;
	model.LoadData(*x, *y, *x_new, is_numerics);
	model.Fit();

	double train_error = 0;
	for (int i = 0; i < y_new->size(); i++) {
		train_error += (model.y_test_pred[i] - y_new->at(i)) * (model.y_test_pred[i] - y_new->at(i));
	}
	train_error /= model.y_test_pred.size();
	std::cout << "Tuning set error " << train_error << "\n";

	model.Clean();
	delete x, y, x_new, y_new;

}

void testWriteCsv() {
	std::vector<double> y;
	for (int i = 0; i < 1e6; ++i) {
		y.push_back(i);
	}
	WriteCsv("c:\\Users\\Chandler\\Kaggle\\BNP\\data\\test_cpp_output.csv", y);
}