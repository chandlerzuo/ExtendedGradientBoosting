#pragma once
#include "rule.h"
#include <vector>

struct AddRuleResult{
	Rule* rule;
	bool pos_add;
	bool negate;
	double loss;
};

class RuleSet {
public:
	static int max_rules;
	int n_true;
	int n_rules;
	std::vector<Rule*> rule_set;
	std::map<bool, double> predict_values;

	RuleSet();
	~RuleSet();
	void Clean();
	void LoadData(std::vector<double> &, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &, std::vector<bool> &, std::vector<int> &, std::vector<int> & varlist);
	void LoadData(std::vector<double> &, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &, std::vector<bool> &, std::vector<int> &);
	void UnloadData();
	bool Check(std::vector<double>);
	void Check();
	double Predict(std::vector<double>);
	void Build();
	AddRuleResult AddRule(int);
	void Print();
	void Print(int);

private:
	std::vector<double>* y;
	std::vector<std::vector<double>>* x;
	std::vector<std::vector<double>>* x_t;
	std::vector<bool>* is_numerics;
	std::vector<int> varlist;
	std::vector<int>* members;
	std::vector<bool> ruleset_check;

	AddRuleResult AddRuleForNumeric(int);
	AddRuleResult AddRuleForFactor(int);
	AddRuleResult FastBuild(std::map<double, int>, int);
	AddRuleResult FullBuild(std::map<double, int>, int);
};