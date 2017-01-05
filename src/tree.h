#pragma once
#include "ruleset.h"

class Node {
public:
	static int max_depth;
	static int min_leaf_size;
	static double subsample;
	static double col_subsample;

	Node* parent;
	RuleSet rule_set;
	std::vector<Node *> children;
	double value;
	int level;
	std::vector<int> members;
	std::vector<int> varlist;
	double loss;

	Node();
	~Node();
	double Predict(std::vector<double> & );
	std::vector<double> Predict(std::vector<std::vector<double>> &);
	void LoadData(std::vector<double>& , std::vector<std::vector<double>>& , std::vector<std::vector<double>>& , std::vector<bool> & );
	void Build();
	void UnloadData();
	void Print();

private:
	std::vector<double>* y;
	std::vector<std::vector<double>>* x;
	std::vector<std::vector<double>> *x_t;
	std::vector<bool> *is_numerics;

	void SetValue();
	void AllocateChildrenMembers();
	void GenerateChildren();
};