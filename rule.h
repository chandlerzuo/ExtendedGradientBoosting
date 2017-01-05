#pragma once
#include <map>

class Rule {
public:
	int var_id;
	bool is_numeric;
	bool negate;
	double threshold;
	std::map<double, int> value_set;

	Rule(int, bool);
	~Rule();
	void SetNegate(bool);
	bool Check(double);
	void Print();
};