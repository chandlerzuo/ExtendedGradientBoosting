/*
#ifdef _DEBUG
#ifndef DBG_NEW
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
#define new DBG_NEW
#endif
#endif  // _DEBUG
*/
#include "stdafx.h"
#include "RuleSet.h"
#include <limits.h>
#include <iostream>
#include <ctime>
#include<algorithm>
#include "utility.h"

//constructor
RuleSet::RuleSet() {
	this->n_true = 0;
	this->n_rules = 0;
}

//destructor
RuleSet::~RuleSet() {
	this->Clean();
}

void RuleSet::Clean() {
	for (int i = 0; i < this->rule_set.size(); ++i) {
		delete this->rule_set[i];
	}
    this->rule_set.clear();
}

//Load Data
// Notice that we pass a copy of the data set.
void RuleSet::LoadData(std::vector<double>& y, std::vector<std::vector<double>>& x, std::vector<std::vector<double>> & x_t, std::vector<bool>& is_numerics, std::vector<int> & members, std::vector<int> & varlist) {
	this->LoadData(y, x, x_t, is_numerics, members);
	this->varlist = varlist;
	//std::cout << "Loaded data x length " << (*this->x).size() << "\n";
}

void RuleSet::LoadData(std::vector<double>& y, std::vector<std::vector<double>>& x, std::vector<std::vector<double>> & x_t, std::vector<bool>& is_numerics, std::vector<int> & members) {
	for (int i = 0; i < x_t.size(); ++i) {
		if (x_t[i].size() != y.size()) {
			std::cout << "Error: data size on dimension " << i << " is " << x_t[i].size() << " unequal to y " << y.size() << "\n";
			exit;
		}
	}
	if (x_t.size() != is_numerics.size()) {
		std::cout << "Error: dimension of is_numerics " << is_numerics.size() << " <> dimension of x " << x_t.size() << "\n";
		exit;
	}
	this->y = &y;
	this->x = &x;
	this->x_t = &x_t;
	this->is_numerics = &is_numerics;
	this->members = &members;
	this->varlist.clear();
	for (int i = 0; i < x_t.size(); ++i) {
		this->varlist.push_back(i);
	}
	//std::cout << "Loaded data x length " << (*this->x).size() << "\n";
}

// Clean the data
void RuleSet::UnloadData() {
	this->y = nullptr;
	this->x = nullptr;
	this->x_t = nullptr;
	this->ruleset_check.clear();
	this->members = nullptr;
}

// Check for a single input
bool RuleSet::Check(std::vector<double> x) {
	if (x.size() != (*this->x_t).size()) {
		std::cout << "Size of input is " << x.size() << " while the number of variables is " << (*this->x_t).size() << ".\n";
		exit;
	}
	int n_true = 0;
	for (int i = 0; i < this->n_rules; ++i) {
		int var_id = this->rule_set[i]->var_id;
		n_true += this->rule_set[i]->Check(x[var_id]);
	}
	return n_true >= this->n_true;
}

// Check for all observations in the loaded data
void RuleSet::Check() {
	this->ruleset_check.clear();
	//std::cout << "Check result : {";
	for (int i = 0; i < (*this->x).size(); ++i) {
		this->ruleset_check.push_back(this->Check((*this->x)[i]));
		//std::cout << " " << this->ruleset_check[i];
	}
	//std::cout << "}\n";
}

//predict
double RuleSet::Predict(std::vector<double> vals) {
	int n_true = 0;
	for (int i = 0; i < this->n_rules; ++i) {
		n_true += this->rule_set[i]->Check(vals[i]);
	}
	return n_true >= this->n_true ? predict_values[1] : predict_values[0];
}

void RuleSet::Build() {
	// indicator for whether each variable was already used in the rule set
	std::vector<bool> used;
	int n_vars = this->varlist.size();
	for (int i = 0; i < (*this->is_numerics).size(); ++i) {
		used.push_back(false);
	}
	this->Check();
	//recursively add rules to the rule set
	double min_loss = std::numeric_limits<double>::max();
	AddRuleResult add_rule_result;
	while(this->n_rules < RuleSet::max_rules) {
		double best_loss = std::numeric_limits<double>::max();
		AddRuleResult best_add_rule_result;
		best_add_rule_result.rule = new Rule(-1, false);// dummy rule; used to be deleted
		//iterate through all variables and add a rule for each one
		for (int i = 0; i < this->varlist.size(); ++i) {
			if (used[this->varlist[i]]) {
				continue;
			}
			//std::cout << "Start add rule for " << i << "\n";
			AddRuleResult new_add_rule_result = this->AddRule(this->varlist[i]);
			//std::cout << "Best loss " << best_loss << " new loss " << new_add_rule_result.loss << "\n";
			if (new_add_rule_result.loss <= best_loss) {
				delete best_add_rule_result.rule;
				best_add_rule_result = new_add_rule_result;
				best_loss = new_add_rule_result.loss;
			}
			else {
				delete new_add_rule_result.rule;
			}
		}
		//Add the best rule
		if (best_loss <= min_loss) {
			this->rule_set.push_back(best_add_rule_result.rule);
			this->n_rules++;
			//best_add_rule_result.rule->Print();
			//std::cout << (best_add_rule_result.rule->negate?"true":"false") << (this->rule_set[0]->negate?"true\n":"false\n");
			this->n_true += best_add_rule_result.pos_add;
			min_loss = best_loss;
			//Check the result of the current rule set
			this->Check();
			used[best_add_rule_result.rule->var_id] = true;
		}
		else {
			// break if the loss cannot be improved
			delete best_add_rule_result.rule;
			break;
		}
		// if there is no additional variable to be used
		if (this->n_rules >= n_vars || min_loss == 0) {
			break;
		}
	}
	double _pred[2];
	double _n[2];
	_pred[0] = 0;
	_pred[1] = 0;
	_n[0] = 0;
	_n[1] = 0;
	for (int i = 0; i < (*this->members).size(); ++i) {
		int obsid = (*this->members)[i];
		if (this->ruleset_check[obsid]) {
			_pred[1] += (*this->y)[obsid];
			_n[1] ++;
		}
		else {
			_pred[0] += (*this->y)[obsid];
			_n[0] ++;
		}
	}
	if (_n[0] == 0) {
		_pred[0] = (_pred[0] + _pred[1]) / (_n[0] + _n[1]);
	}
	else {
		_pred[0] /= _n[0];
	}
	if (_n[1] == 0) {
		_pred[1] = (_pred[0] + _pred[1]) / (_n[0] + _n[1]);
	}
	else {
		_pred[1] /= _n[1];
	}
	this->predict_values.insert(std::pair<bool, double>(false, _pred[0]));
	this->predict_values.insert(std::pair<bool, double>(true, _pred[1]));
	return;
}

AddRuleResult RuleSet::AddRule(int var_id) {
	if ((*this->is_numerics)[var_id]) {
		//std::cout << "Add rule for numeric variable " << var_id << "\n";
		return this->AddRuleForNumeric(var_id);
	}
	else {
		//std::cout << "Add rule for categorical variable " << var_id << "\n";
		return this->AddRuleForFactor(var_id);
	}
}

// Build RuleSet for numeric predictor
AddRuleResult RuleSet::AddRuleForNumeric(int var_id) {
	// candidate threshold
	std::map<double, int> cand_threshold;
	double threshold;
	double min_loss = std::numeric_limits<double>::max();
	bool pos_add, negate;

	// default
	threshold = this->x_t->at(var_id).at(this->members->at(0));
	pos_add = false;
	negate = false;

	/*
	for (int i = 0; i < this->members->size(); ++i) {
		int obsi = this->members->at(i);
		if (cand_threshold.find(this->x_t->at(var_id)[obsi]) != cand_threshold.end()) {
			cand_threshold.insert(std::pair<double, int>(this->x_t->at(var_id).at(obsi), 1));
		}
	}
	*/

	/*
	for (double i = 0.1; i < 1; i+=0.1 ) {
		if (cand_threshold.find(x_copy[(int)(i * x_copy.size())]) == cand_threshold.end()) {
			cand_threshold.insert(std::pair<double, int>(x_copy[(int)(i * x_copy.size())], 1));
			//std::cout << "Added candidate split " << x_copy[(int)(i * x_copy.size())] << "\n";
		}
	}

	const clock_t begin_time = clock();
	for (int addflag = 0; addflag < 4; addflag++) {
		for (std::map<double, int>::iterator it = cand_threshold.begin(); it != cand_threshold.end(); ++ it) {
			//std::cout << "Check threshold " << (*this->x_t)[var_id][i] << " min loss " << min_loss << "\n";
			// skip if this value has already been checked
			// std::cout << "Started " << i << " " << float(clock() - begin_time) / CLOCKS_PER_SEC << " seconds elapsed.\n";
			
			double cand_split = it->first;
			
			// compute loss
			double _pred[2];
			double _n[2];
			_pred[0] = 0;
			_pred[1] = 0;
			_n[0] = 0;
			_n[1] = 0;
			for (int j = 0; j < this->members->size(); j++) {
				int obsj = this->members->at(j);
				if (CombineConditions(this->ruleset_check[obsj], this->x_t->at(var_id).at(obsj) > cand_split, addflag)) {
					_pred[1] += this->y->at(obsj);
					_n[1]++;
				}
				else {
					_pred[0] += this->y->at(obsj);
					_n[0]++;
				}
			}
			// do not compute loss if trivial
			if (_n[0] == 0 || _n[1] == 0) {
				continue;
			}
			// compute the sum of squared error
			double loss = 0;
			_pred[0] /= _n[0];
			_pred[1] /= _n[1];
			for (int j = 0; j < this->members->size(); j++) {
				int obsj = this->members->at(j);
				if (CombineConditions(this->ruleset_check[obsj], this->x_t->at(var_id).at(obsj) > cand_split, addflag)) {
					loss += (this->y->at(obsj) - _pred[1]) * (this->y->at(obsj) - _pred[1]);
				}
				else {
					loss += (this->y->at(obsj) - _pred[0]) * (this->y->at(obsj) - _pred[0]);
				}
				if (loss >= min_loss) {
					break;
				}
			}
			if (loss < min_loss) {
				threshold = cand_split;
				min_loss = loss;
				//	std::cout << "Updated minimum loss " << min_loss << "\n";
				pos_add = (addflag == 1) || (addflag == 3);
				negate = (addflag == 2) || (addflag == 3);
			}
		}
	}
    //std::cout << "Finished iterations, " << float(clock() - begin_time) / CLOCKS_PER_SEC << " seconds elapsed.\n";
	*/

	const clock_t begin_time = clock();

	std::map<double, std::vector<double>> xy_pairs0, xy_pairs1;
	double y2_sum = 0;
	for (int i = 0; i < this->members->size(); ++i) {
		int obsi = this->members->at(i);
		y2_sum += this->y->at(obsi) * this->y->at(obsi);
		if (this->ruleset_check.at(i)) {
			if (xy_pairs1.find(this->x_t->at(var_id)[obsi]) == xy_pairs1.end()) {
				xy_pairs1[this->x_t->at(var_id)[obsi]];
				xy_pairs1[this->x_t->at(var_id)[obsi]].push_back(this->y->at(obsi));
				xy_pairs1[this->x_t->at(var_id)[obsi]].push_back(1);
			}
			else {
				xy_pairs1[this->x_t->at(var_id)[obsi]][0] += this->y->at(obsi);
				xy_pairs1[this->x_t->at(var_id)[obsi]][1] ++;
			}
		}
		else {
			if (xy_pairs0.find(this->x_t->at(var_id)[obsi]) == xy_pairs0.end()) {
				xy_pairs0[this->x_t->at(var_id)[obsi]];
				xy_pairs0[this->x_t->at(var_id)[obsi]].push_back(this->y->at(obsi));
				xy_pairs0[this->x_t->at(var_id)[obsi]].push_back(1);
			}
			else {
				xy_pairs0[this->x_t->at(var_id)[obsi]][0] += this->y->at(obsi);
				xy_pairs0[this->x_t->at(var_id)[obsi]][1] ++;
			}
		}
	}

	std::vector<std::vector<double>> xy0, xy1;
	for (std::map<double, std::vector<double>>::iterator it = xy_pairs0.begin(); it != xy_pairs0.end(); ++it) {
		xy0.push_back(std::vector<double>{it->first, it->second[0], it->second[1]});
	}
	for (std::map<double, std::vector<double>>::iterator it = xy_pairs1.begin(); it != xy_pairs1.end(); ++it) {
		xy1.push_back(std::vector<double>{it->first, it->second[0], it->second[1]});
	}

	std::sort(xy0.begin(), xy0.end(), less_than_first_coord());
	std::sort(xy1.begin(), xy1.end(), less_than_first_coord());

	// compute partial sums and partial sum of squares
	std::vector<double> y0_cumsum, y1_cumsum, x0, x1, n0_cumsum, n1_cumsum;
	double n1_sum = 0, n0_sum = 0;
	double y_sum = 0, y0_sum = 0, y1_sum = 0;
	for (int i = 0; i < xy0.size(); i++) {
		y0_sum += xy0[i][1];
		n0_sum += xy0[i][2];
		x0.push_back(xy0[i][0]);
		y0_cumsum.push_back(y0_sum);
		n0_cumsum.push_back(n0_sum);
	}
	for (int i = 0; i < xy1.size(); i++) {
		y1_sum += xy1[i][1];
		n1_sum += xy1[i][2];
		x1.push_back(xy1[i][0]);
		y1_cumsum.push_back(y1_sum);
		n1_cumsum.push_back(n1_sum);
	}

	//std::cout << "Finished sorting x, " << float(clock() - begin_time) / CLOCKS_PER_SEC << " seconds elapsed.\n";
	/*
	std::cout << "Neg: \n";
	for (int i = 0; i < x0.size(); ++i) {
		std::cout << x0[i] << " " << y0_cumsum[i] << "\n";
	}
	std::cout << "Pos: \n";
	for (int i = 0; i < x1.size(); ++i) {
		std::cout << x1[i] << " " << y1_cumsum[i] << "\n";
	}
	*/

	for (int addflag = 0; addflag < 4; addflag += 2) {
		if (y0_cumsum.size() <= 1) {
			continue;
		}
		for (int i = 0; i < y0_cumsum.size() - 1; i ++ ) {
			//std::cout << "Check threshold " << (*this->x_t)[var_id][i] << " min loss " << min_loss << "\n";
			// skip if this value has already been checked
			// std::cout << "Started " << i << " " << float(clock() - begin_time) / CLOCKS_PER_SEC << " seconds elapsed.\n";
			double cand_split = x0[i];
			// compute loss
			double _pred[2];
			double _n[2];
			if (addflag == 0) {
				// A || B
				_pred[0] = y0_cumsum[i];
				_n[0] = n0_cumsum[i];
			}
			else {
				// A || ^B
				_pred[0] = y0_sum - y0_cumsum[i];
				_n[0] = n0_sum - n0_cumsum[i];
			}
			_pred[1] = y0_sum + y1_sum - _pred[0];
			_n[1] = n0_sum + n1_sum - _n[0];
			if (_n[1] == 0 || _n[0] == 0) {
				continue;
			}
			// compute the sum of squared error
			_pred[0] /= _n[0];
			_pred[1] /= _n[1];
			double loss = y2_sum - _n[0] * _pred[0] * _pred[0] - _n[1] * _pred[1] * _pred[1];
			if (loss < min_loss) {
				threshold = cand_split;
				min_loss = loss;
				//	std::cout << "Updated minimum loss " << min_loss << "\n";
				pos_add = (addflag == 1) || (addflag == 3);
				negate = (addflag == 2) || (addflag == 3);
			}
		}
	}

	for (int addflag = 1; addflag < 4; addflag += 2) {
		if (y1_cumsum.size() <= 1) {
			continue;
		}
		for (int i = 0; i < y1_cumsum.size() - 1; ++i) {
			//std::cout << "Check threshold " << (*this->x_t)[var_id][i] << " min loss " << min_loss << "\n";
			// skip if this value has already been checked
			// std::cout << "Started " << i << " " << float(clock() - begin_time) / CLOCKS_PER_SEC << " seconds elapsed.\n";

			double cand_split = x1[i];
			// compute loss
			double _pred[2];
			double _n[2];
			if (addflag == 1) {
				// mode: A & B
				_pred[1] = y1_sum - y1_cumsum[i];
				_n[1] = n1_sum - n1_cumsum[i];
			}
			else {
				// mode: A & ^B
				_pred[1] = y1_cumsum[i];
				_n[1] = n1_cumsum[i];
			}
			_pred[0] = y0_sum + y1_sum - _pred[1];
			_n[0] = n0_sum + n1_sum - _n[1];
			if (_n[1] == 0 || _n[0] == 0) {
				continue;
			}
			// compute the sum of squared error
			_pred[0] /= _n[0];
			_pred[1] /= _n[1];
			double loss = y2_sum - _n[0] * _pred[0] * _pred[0] - _n[1] * _pred[1] * _pred[1];
			if (loss < min_loss) {
				threshold = cand_split;
				min_loss = loss;
				//	std::cout << "Updated minimum loss " << min_loss << "\n";
				pos_add = (addflag == 1) || (addflag == 3);
				negate = (addflag == 2) || (addflag == 3);
				// std::cout << "Update loss to " << min_loss << " threshold " << cand_split << " pred " << _pred[0] << " " << _pred[1] << "\n";
			}
		}
	}

	// put the best prediction
	AddRuleResult result;
	result.loss = min_loss;
	result.pos_add = pos_add;
	result.rule = new Rule(var_id, true);
	result.rule->threshold = threshold;
	result.rule->SetNegate(negate);
	return result;
}

// Build RuleSet for categorical predictor
AddRuleResult RuleSet::AddRuleForFactor(int var_id) {
	// default
	AddRuleResult result;		
	result.rule = new Rule(var_id, false);

	result.rule->value_set.insert(std::pair<double, int>((*this->x_t)[var_id][(*this->members)[0]], 1));
	result.loss = std::numeric_limits<double>::max();
	result.pos_add = false;

	std::map<double, std::vector<double>> factor_levels1, factor_levels0;
	double y2_sum = 0;
	for (int i = 0; i < (*this->members).size(); i++) {
		int obsi = (*this->members)[i];
		y2_sum += this->y->at(obsi) * this->y->at(obsi);
		if (this->ruleset_check[obsi]) {
			if (factor_levels1.find((*this->x_t)[var_id][obsi]) == factor_levels1.end()) {
				factor_levels1[this->x_t->at(var_id)[obsi]];
				factor_levels1[this->x_t->at(var_id)[obsi]].push_back(this->y->at(obsi));
				factor_levels1[this->x_t->at(var_id)[obsi]].push_back(1);
			}
			else {
				factor_levels1[this->x_t->at(var_id)[obsi]][0] += this->y->at(obsi);
				factor_levels1[this->x_t->at(var_id)[obsi]][1] ++;
			}
		}
		else {
			if (factor_levels0.find((*this->x_t)[var_id][obsi]) == factor_levels0.end()) {
				factor_levels0[this->x_t->at(var_id)[obsi]];
				factor_levels0[this->x_t->at(var_id)[obsi]].push_back(this->y->at(obsi));
				factor_levels0[this->x_t->at(var_id)[obsi]].push_back(1);
			}
			else {
				factor_levels0[this->x_t->at(var_id)[obsi]][0] += this->y->at(obsi);
				factor_levels0[this->x_t->at(var_id)[obsi]][1] ++;
			}
		}
	}

	// insert into set
	// each vector is {mean of y, x, size of bucket}
	std::vector<std::vector<double>> xy_pairs0, xy_pairs1;
	for (std::map<double, std::vector<double>>::iterator it = factor_levels1.begin(); it != factor_levels1.end(); ++it) {
		xy_pairs1.push_back(std::vector<double>{it->second[0] / it->second[1], it->first, it->second[1]});
	}
	for (std::map<double, std::vector<double>>::iterator it = factor_levels0.begin(); it != factor_levels0.end(); ++it) {
		xy_pairs0.push_back(std::vector<double>{it->second[0] / it->second[1], it->first, it->second[1]});
	}
	std::sort(xy_pairs0.begin(), xy_pairs0.end(), less_than_first_coord());
	std::sort(xy_pairs1.begin(), xy_pairs1.end(), less_than_first_coord());

	// compute partial sums of each set
	std::vector<double> y0_cumsum, y1_cumsum, n0_cumsum, n1_cumsum;
	double y0_sum = 0, y1_sum = 0, n0_sum = 0, n1_sum = 0;
	for (int i = 0; i < xy_pairs0.size(); ++i) {
		y0_sum += xy_pairs0[i][0] * xy_pairs0[i][2];
		n0_sum += xy_pairs0[i][2];
		y0_cumsum.push_back(y0_sum);
		n0_cumsum.push_back(n0_sum);
	}
	for (int i = 0; i < xy_pairs1.size(); ++i) {
		y1_sum += xy_pairs1[i][0] * xy_pairs1[i][2];
		n1_sum += xy_pairs1[i][2];
		y1_cumsum.push_back(y1_sum);
		n1_cumsum.push_back(n1_sum);
	}

	/*
	std::cout << "Neg: \n";
	for (int i = 0; i < xy_pairs0.size(); ++i) {
		std::cout << xy_pairs0[i][1] << " " << y0_cumsum[i] << " " << n0_cumsum[i] << "\n";
	}
	std::cout << "Pos: \n";
	for (int i = 0; i < xy_pairs1.size(); ++i) {
		std::cout << xy_pairs1[i][1] << " " << y1_cumsum[i] << " " << n1_cumsum[i] << "\n";
	}
	*/

	for (int addflag = 0; addflag < 4; addflag += 2) {
		if (y0_cumsum.size() <= 1) {
			continue;
		}
		for (int i = 0; i < y0_cumsum.size() - 1; ++i) {
			double _pred[2], _n[2];
			if (addflag == 0) {
				// A|B
				_pred[0] = y0_sum - y0_cumsum[i];
				_n[0] = n0_sum - n0_cumsum[i];
			}
			else {
				// A|^B
				_pred[0] = y0_cumsum[i];
				_n[0] = n0_cumsum[i];
			}
			_pred[1] = y0_sum + y1_sum - _pred[0];
			_n[1] = n0_sum + n1_sum - _n[0];
			if (_n[0] < 1 || _n[1] < 1) {
				continue;
			}
			_pred[0] /= _n[0];
			_pred[1] /= _n[1];
			double loss = y2_sum - _pred[0] * _pred[0] * _n[0] - _pred[1] * _pred[1] * _n[1];
			if (loss < result.loss) {
				result.loss = loss;
				result.rule->value_set.clear();
				for (int j = 0; j <= i; j++) {
					result.rule->value_set.insert(std::pair<double, int>(xy_pairs0[j][1], 1));
				}
				result.pos_add = (addflag == 1) || (addflag == 3);
				result.rule->SetNegate((addflag == 2) || (addflag == 3));
			}
		}
	}

	for (int addflag = 1; addflag < 4; addflag += 2) {
		if (y1_cumsum.size() <= 1) {
			continue;
		}
		for (int i = 0; i < y1_cumsum.size() - 1; ++i) {
			double _pred[2], _n[2];
			if (addflag == 1) {
				// A&B
				_pred[1] = y1_cumsum[i];
				_n[1] = n1_cumsum[i];
			}
			else {
				// A&^B
				_pred[1] = y1_sum - y1_cumsum[i];
				_n[1] = n1_sum - n1_cumsum[i];
			}
			_pred[0] = y0_sum + y1_sum - _pred[1];
			_n[0] = n0_sum + n1_sum - _n[1];
			if (_n[0] < 1 || _n[1] < 1) {
				continue;
			}
			_pred[0] /= _n[0];
			_pred[1] /= _n[1];
			double loss = y2_sum - _pred[0] * _pred[0] * _n[0] - _pred[1] * _pred[1] * _n[1];
			if (loss < result.loss) {
				result.loss = loss;
				result.rule->value_set.clear();
				for (int j = 0; j <= i; j++) {
					result.rule->value_set.insert(std::pair<double, int>(xy_pairs1[j][1], 1));
				}
				result.pos_add = (addflag == 1) || (addflag == 3);
				result.rule->SetNegate((addflag == 2) || (addflag == 3));
			}
		}
	}

	return result;

	/*
	std::map<double, int> factor_levels;
	for (int i = 0; i < (*this->members).size(); i++) {
	int obsi = (*this->members)[i];
	if (factor_levels.find((*this->x_t)[var_id][obsi]) == factor_levels.end()) {
	factor_levels.insert(std::pair<double, int>((*this->x_t)[var_id][obsi], 1));
	}
	}
	if (factor_levels.size() > 10) {
		std::cout << "Factor levels " << factor_levels.size() << " invoke fast build" << "\n";
		return this->FastBuild(factor_levels, var_id);
	}
	else {
		std::cout << "Factor levels " << factor_levels.size() << " invoke full build" << "\n";
		return this->FullBuild(factor_levels, var_id);
	}
	*/
}

//Build the RuleSet for categorical variable with forward selection
AddRuleResult RuleSet::FastBuild(std::map<double, int> factor_levels, int var_id) {
	AddRuleResult result;
	result.rule = new Rule(var_id, false);
	result.loss = std::numeric_limits<double>::max();
	for (int addflag = 0; addflag < 4; ++addflag) {
		for (std::map<double, int>::iterator it = factor_levels.begin(); it != factor_levels.end(); ++it) {
			double _pred[2];
			double _n[2];
			double loss = 0;
			_pred[0] = 0;
			_pred[1] = 0;
			_n[0] = 0;
			_n[1] = 0;
			for (int i = 0; i < (*this->members).size(); i++) {
				int obsi = (*this->members)[i];
				if (CombineConditions(this->ruleset_check[obsi], ((result.rule->value_set.find((*this->x_t)[var_id][obsi]) != result.rule->value_set.end()) || (*this->x_t)[var_id][obsi] == it->first), addflag)) {
					_pred[1] += (*this->y)[obsi];
					_n[1] ++;
				}
				else {
					_pred[0] += (*this->y)[obsi];
					_n[0] ++;
				}
			}
		// compute loss after adding this level
			if (_n[0] == 0 || _n[1] == 0) {
				loss = std::numeric_limits<double>::max();
			}
			else {
				_pred[0] /= _n[0];
				_pred[1] /= _n[1];
				for (int i = 0; i < (*this->members).size(); i++) {
					int obsi = (*this->members)[i];
					if (CombineConditions(this->ruleset_check[obsi], ((result.rule->value_set.find((*this->x_t)[var_id][obsi]) != result.rule->value_set.end()) || (*this->x_t)[var_id][obsi] == it->first), addflag)) {
						loss += ((*this->y)[obsi] - _pred[1]) * ((*this->y)[obsi] - _pred[1]);
					}
					else {
						loss += ((*this->y)[obsi] - _pred[0]) * ((*this->y)[obsi] - _pred[0]);
					}
				}
			}
			//decide whether add this level
			if (loss < result.loss) {
				result.rule->value_set.insert(std::pair<double, int>(it->first, 1));
				result.loss = loss;
				result.pos_add = (addflag == 1) || (addflag == 3);
				result.rule->SetNegate((addflag == 2) || (addflag == 3));
			}
		}
	}
	return result;
}

AddRuleResult RuleSet::FullBuild(std::map<double, int> factor_levels, int var_id) {
	std::vector<double> _factor_levels;
	std::vector<bool> _include;
	AddRuleResult result;
	result.rule = new Rule(var_id, false);
	result.loss = std::numeric_limits<double>::max();
	for (std::map<double, int>::iterator it = factor_levels.begin(); it != factor_levels.end(); ++it) {
		_factor_levels.push_back(it->first);
	}
	for (int addflag = 0; addflag < 4; ++addflag) {
		_include.clear();
		for (std::map<double, int>::iterator it = factor_levels.begin(); it != factor_levels.end(); ++it) {
			_include.push_back(false);
		}
		int pos = 0;
		double loss;
		while (pos < _include.size() - 1 || !_include[pos]) {
			// This is enumerating all subsets of factor levels
			_include[pos] = true;
			for (int i = 0; i < pos; i++) {
				_include[i] = false;
			}
			std::map<double, int> _value_set;
			for (int i = 0; i < _include.size(); i++) {
				if (_include[i]) {
					_value_set.insert(std::pair<double, int>(_factor_levels[i], 1));
				}
			}
			loss = 0;
			double _pred[2];
			double _n[2];
			_n[0] = 0;
			_n[1] = 0;
			_pred[0] = 0;
			_pred[1] = 0;
			for (int i = 0; i < (*this->members).size(); i++) {
				int obsi = (*this->members)[i];
				if (CombineConditions(this->ruleset_check[obsi],
					_value_set.find((*this->x_t)[var_id][obsi]) != _value_set.end(),
					addflag)) {
					_n[1] ++;
					_pred[1] += (*this->y)[obsi];
				}
				else {
					_n[0] ++;
					_pred[0] += (*this->y)[obsi];
				}
			}

			if (_n[0] == 0 || _n[1] == 0) {
				loss = std::numeric_limits<double>::max();
			}
			else {
				_pred[0] /= _n[0];
				_pred[1] /= _n[1];
				for (int i = 0; i < (*this->members).size(); i++) {
					int obsi = (*this->members)[i];
					if (CombineConditions(this->ruleset_check[obsi],
						_value_set.find((*this->x_t)[var_id][obsi]) != _value_set.end(),
						addflag)) {
						loss += ((*this->y)[obsi] - _pred[1]) * ((*this->y)[obsi] - _pred[1]);
					}
					else {
						loss += ((*this->y)[obsi] - _pred[0]) * ((*this->y)[obsi] - _pred[0]);
					}
				}
			}
			if (loss < result.loss) {
				result.rule->value_set = _value_set;
				result.loss = loss;
				result.pos_add = (addflag == 1) || (addflag == 3);
				result.rule->SetNegate((addflag == 2) || (addflag == 3));
			}

			pos = 0;
			while (pos < _include.size() - 1 && _include[pos]) {
				pos++;
			}
		}
	}
	return result;
}

void RuleSet::Print() {
	this->Print(0);
}

void RuleSet::Print(int n) {
	PrintLeadingSpaces(n);
	std::cout << "Print RuleSet. " << this->n_true << " of the following has to be true: \n";
	for (int i = 0; i < this->n_rules; i++) {
		PrintLeadingSpaces(n);
		this->rule_set[i]->Print();
	}
	PrintLeadingSpaces(n);
	std::cout << "Predicted values: false:" << predict_values[false] << ", true: " << predict_values[true] << "\n";
}