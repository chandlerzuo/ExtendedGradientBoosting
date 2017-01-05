#include "stdafx.h"
#include "tree.h"
#include "utility.h"
#include <limits.h>
#include <iostream>

Node::Node() {
	this->value = 0;
	this->level = 0;
	this->parent = nullptr;
}

Node::~Node() {
	this->parent = nullptr;
}

//Load Data
// Notice that we pass a copy of the data set.
void Node::LoadData(std::vector<double>& y, std::vector<std::vector<double>>& x, std::vector<std::vector<double>> & x_t, std::vector<bool>& is_numerics) {
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
	// subsample
	this->varlist = Sample(x_t.size(), Node::col_subsample);
	//std::cout << "Loaded data x length " << (*this->x).size() << "\n";
}

// Clean the data
void Node::UnloadData() {
	this->y = nullptr;
	this->x = nullptr;
	this->x_t = nullptr;
}

// Predict for one single observation.
double Node::Predict(std::vector<double> & x) {
	if (this->children.size() == 0) {
		return this->value;
	}
	return this->rule_set.Check(x) ? this->children[1]->Predict(x) : this->children[0]->Predict(x);
}

// Predict for a list of observations.
std::vector<double> Node::Predict(std::vector<std::vector<double>> & x) {
	std::vector<double> ret;
	for (int i = 0; i < x.size(); ++i) {
		ret.push_back(this->Predict(x[i]));
	}
	return ret;
}

//Fit for a tree
// This function is supposed to be called at root.
void Node::Build() {
	//error if there is no data
	if (this->y->size() == 0) {
		std::cout << "Error: no data loaded." << "\n";
	}
	// this is the root
	this->members = Sample(this->x->size(), Node::subsample);
	this->level = 0;
	this->SetValue();
	this->GenerateChildren();

	std::vector<double> pred = this->Predict(*this->x);
	this->loss = 0;
	double y2_sum = 0, y_sum = 0;
	for (int i = 0; i < this->x->size(); i++) {
		this->loss += (pred[i] - this->y->at(i)) * (pred[i] - this->y->at(i));
		y_sum += this->y->at(i);
		y2_sum += this->y->at(i) * this->y->at(i);
	}
	double var = y2_sum - y_sum * y_sum / (double) pred.size();
	if (var == 0) {
		this->loss = this->loss == 0 ? 0 : 100;
	}
	else {
		this->loss /= var;
	}
}

void Node::AllocateChildrenMembers() {
	this->children[0]->members.clear();
	this->children[1]->members.clear();
	for (int i = 0; i < this->members.size(); ++i) {
		if (this->rule_set.Check((*this->x)[this->members[i]])) {
			this->children[1]->members.push_back(this->members[i]);
		}
		else {
			this->children[0]->members.push_back(this->members[i]);
		}
	}
}

void Node::SetValue() {
	// if this node has size 0, then return the parent value
	if (this->members.size() == 0) {
		if (this->parent == nullptr) {
			this->value = 0;
			return;
		}
		this->value = this->parent->value;
		return;
	}
	this->value = 0;
	for (int i = 0; i < this->members.size(); ++i) {
		this->value += (*this->y)[this->members[i]];
	}
	this->value /= this->members.size();
	//std::cout << "Value at this node: " << this->value << "\n";
}

// Recursively generate the children nodes
// Note: before calling this function the children vector must be empty
void Node::GenerateChildren() {
	this->SetValue();
	if (this->members.size() <= Node::min_leaf_size || this->level >= Node::max_depth) {
		return;
	}
	// check if all values on this node are the same
	bool single_value = true;
	for (int i = 1; i < this->members.size(); ++i) {
		if (this->y->at(this->members[i]) != this->y->at(this->members[0])) {
			single_value = false;
			break;
		}
	}
	if (single_value) {
		return;
	}
	//Build the rule
	this->rule_set = RuleSet();
	this->rule_set.LoadData(* this->y, * this->x, * this->x_t, * this->is_numerics, this->members, this->varlist);
	//std::cout << "Build rule at node on level " << this->level << "\n";
	this->rule_set.Build();
	//std::cout << "Value at splits: " << this->rule_set.predict_values[0] << " " << this->rule_set.predict_values[1] << "\n";
	//Set children nodes
	this->children.push_back(new Node());
	this->children.push_back(new Node());
	for (int i = 0; i < 2; ++i) {
		this->children[i]->LoadData(* this->y, * this->x, * this->x_t, * this->is_numerics);
		this->children[i]->parent = this;
		this->children[i]->level = this->level + 1;
	}
	// Set members to the children nodes
	this->AllocateChildrenMembers();
	// Fit the rule of each children
	for (int i = 0; i < 2; i++) {
		this->children[i]->GenerateChildren();
	}
}

// Print the tree result
void Node::Print() {
	PrintLeadingSpaces(this->level);
	std::cout << "Value: " << this->value << "\n";
	this->rule_set.Print(this->level);
	for (int i = 0; i < this->children.size(); ++i) {
		this->children[i]->Print();
	}
}