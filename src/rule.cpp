/*
#ifdef _DEBUG
#ifndef DBG_NEW
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
#define new DBG_NEW
#endif
#endif  // _DEBUG
*/
#include "stdafx.h"
#include "rule.h"
#include <iostream>

//constructor
Rule::Rule(int var_id, bool is_numeric) {
	this->var_id = var_id;
	this->is_numeric = is_numeric;
	this->negate = false;
}

//destructor
Rule::~Rule() {
//	this->value_set.clear();
}

bool Rule::Check(double x) {
	bool ret;
	if (this->is_numeric) {
		ret = x > this->threshold;
	}
	else {
		ret = this->value_set.find(x) != this->value_set.end();
	}
	return this->negate ? !ret : ret;
}

void Rule::Print() {
	if (this->is_numeric) {
		std::cout << "Numeric variable X" << this->var_id << (this->negate ? "<=" : ">" ) << this->threshold << "\n";
	}
	else {
		std::cout << "Categorical variable X" << this->var_id <<  ( this->negate ? " not" : "") << " in set {";
		for (std::map<double, int>::iterator it = this->value_set.begin(); it != this->value_set.end(); ++it) {
			std::cout << " " << it->first;
		}
		std::cout << " }\n";
	}
}

void Rule::SetNegate(bool negate) {
	this->negate = negate;
}