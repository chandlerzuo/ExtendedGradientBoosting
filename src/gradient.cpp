#include "gradient.h"
#include <math.h>
#include <cmath>
#include <iostream>

Gradient::Gradient(int J) {
	this->type = 0;
	this->J = J;
	this->loglik = 0;
}

// Initialize with type
Gradient::Gradient(int J, int type) {
	this->type = type;
	this->J = J;
	this->loglik = 0;
}

void Gradient::LoadData(std::vector<double> & y) {
	this->y = &y;
}

double Gradient::trim(double val) {
	if (val > 1 - Gradient::min_prob) {
		val = 1 - Gradient::min_prob;
	}
	else if (val < Gradient::min_prob) {
		val = Gradient::min_prob;
	}
	return val;
}

void Gradient::Initialize() {
	double avg_y = 0;
	// numeric or binary response
	if (this->type == 0 || this->type == 1) {
		//std::cout << "Begin initializing ... \n";
		for (int i = 0; i < this->y->size(); ++i) {
			avg_y += this->y->at(i);
		}
	}
	avg_y /= (double) this->y->size();
	this->Initialize(avg_y);
}

void Gradient::Initialize(double avg_y) {
	this->Initialize(avg_y, this->y->size());
	// Initialize gradient with f
	for (int i = 0; i < this->y->size(); ++i) {
		if (this->type == 0) {
			// only for binary
			if (this->y->at(i) == 0) {
				this->df.push_back(-1 / (1 - this->f[i]));
				this->d2f.push_back(-1 / (1 - this->f[i]) / (1 - this->f[i]));
			}
			else if (this->y->at(i) == 1) {
				this->df.push_back(1 / this->f[i]);
				this->d2f.push_back(-1 / this->f[i] / this->f[i]);
			}
		}
	}
}

void Gradient::Initialize(double avg_y, int n) {
	avg_y = this->trim(avg_y);
	// Initialize g, h, d, p
	for (int j = 0; j < this->J; j++) {
		//std::cout << "Initialize for j = " << j << "\n";
		std::vector<double> _g, _h, _d, _p, _dg, _d2g, _dh, _d2h, _fdh, _fd2h, _fdg, _fd2g, _decent_h, _decent_g;
		for (int i = 0; i < n; i++) {
			_g.push_back(0);
			_h.push_back(log(avg_y / (1 - avg_y)));
			_d.push_back(avg_y);
			_p.push_back(1 / ((double) this->J));
			_dg.push_back(0);
			_dh.push_back(0);
			_d2g.push_back(1);
			_d2h.push_back(1);
			_fdh.push_back(0);
			_fd2h.push_back(1);
			_fdg.push_back(0);
			_fd2g.push_back(1);
			_decent_h.push_back(0);
			_decent_g.push_back(0);
			if (j == 0) {
				this->f.push_back(avg_y);
			}
		}
		this->d.push_back(_d);
		this->p.push_back(_p);
		this->g.push_back(_g);
		this->h.push_back(_h);
		this->dg.push_back(_dg);
		this->dh.push_back(_dh);
		this->d2g.push_back(_d2g);
		this->d2h.push_back(_d2h);
		this->fdh.push_back(_fdh);
		this->fd2h.push_back(_fd2h);
		this->fdg.push_back(_fdg);
		this->fd2g.push_back(_fd2g);
		this->decent_h.push_back(_decent_h);
		this->decent_g.push_back(_decent_g);
	}
}

void Gradient::ComputeGradient_h(int j) {
	for (int i = 0; i < this->f.size(); ++i) {
		this->fdh[j][i] = this->p[j][i] * this->d[j][i] * (1 - this->d[j][i]);
		this->fd2h[j][i] = this->fdh[j][i] * (1 - 2 * this->d[j][i]);
		this->dh[j][i] = this->df[i] * this->fdh[j][i];
		this->d2h[j][i] = this->d2f[i] * this->fdh[j][i] * this->fdh[j][i] + this->df[i] * this->fd2h[j][i];
		if (this->d2h[j][i] < Gradient::min_prob && this->d2h[j][i] >= 0) {
			this->d2h[j][i] = Gradient::min_prob;
		}
		else if (this->d2h[j][i] > - Gradient::min_prob && this->d2h[j][i] < 0) {
			this->d2h[j][i] = - Gradient::min_prob;
		}
		this->decent_h[j][i] = - this->dh[j][i] / this->d2h[j][i];
	}
}

void Gradient::ComputeGradient_g(int j) {
	if (j == 0) {
		return;
	}
	for (int i = 0; i < this->f.size(); ++i) {
		this->fdg[j][i] = (this->d[j][i] - this->f[i]) * this->p[j][i];
		this->fd2g[j][i] = (this->d[j][i] - this->d[j][i] * this->p[j][i] - this->f[i]) * this->fdg[j][i];
		this->dg[j][i] = this->df[i] * this->fdg[j][i];
		this->d2g[j][i] = this->d2f[i] * this->fdg[j][i] * this->fdg[j][i] + this->df[i] * this->fd2g[j][i];
		if (this->d2g[j][i] < Gradient::min_prob && this->d2g[j][i] >= 0) {
			this->d2g[j][i] = Gradient::min_prob;
		}
		else if (this->d2g[j][i] > -Gradient::min_prob && this->d2g[j][i] < 0) {
			this->d2g[j][i] = -Gradient::min_prob;
		}
		this->decent_g[j][i] = - this->dg[j][i] / this->d2g[j][i];
	}
}

void Gradient::Update_g(int j, std::vector<double> & dg) {
	if (dg.size() != this->f.size()) {
		std::cout << "Dimension error: updated value has length " << dg.size() << " while it should be " << this->f.size() << "\n";
	}
	for (int i = 0; i < dg.size(); ++i) {
		this->g[j][i] += dg[i];
		double total = 0;
		for (int j = 0; j < this->J; j++) {
			total += exp(this->g[j][i]);
		}
		this->f[i] = 0;
		for (int j = 0; j < this->J; j++) {
			this->p[j][i] = this->trim(exp(this->g[j][i]) / total);
			this->f[i] += this->p[j][i] * this->d[j][i];
		}
		// only for binary
		if (this->y != nullptr) {
			if (this->y->at(i) == 0) {
				this->df[i] = -1 / (1 - this->f[i]);
				this->d2f[i] = -1 / (1 - this->f[i]) / (1 - this->f[i]);
			}
			else if (this->y->at(i) == 1) {
				this->df[i] = 1 / this->f[i];
				this->d2f[i] = -1 / this->f[i] / this->f[i];
			}
		}
	}
}

void Gradient::Update_h(int j, std::vector<double> & dh) {
	if (dh.size() != this->f.size()) {
		std::cout << "Dimension error: updated value has length " << dh.size() << " while it should be " << this->f.size() << "\n";
	}
	for (int i = 0; i < dh.size(); ++i) {
		this->h[j][i] += dh[i];
		this->f[i] -= this->p[j][i] * this->d[j][i];
		this->d[j][i] = this->trim(exp(this->h[j][i]) / (1 + exp(this->h[j][i])));
		this->f[i] += this->p[j][i] * this->d[j][i];
		// only for binary
		if (this->y != nullptr) {
			if (this->y->at(i) == 0) {
				this->df[i] = -1 / (1 - this->f[i]);
				this->d2f[i] = -1 / (1 - this->f[i]) / (1 - this->f[i]);
			}
			else if (this->y->at(i) == 1) {
				this->df[i] = 1 / this->f[i];
				this->d2f[i] = -1 / this->f[i] / this->f[i];
			}
		}
	}
}

void Gradient::UpdateLikelihood() {
	this->loglik = 0;
	if (this->y == nullptr) {
		return;
	}
	for (int i = 0; i < this->y->size(); ++i) {
		this->loglik += this->y->at(i) == 1 ? log(this->f[i]) : log(1 - this->f[i]);
	}
}

double Gradient::LossWithUpdate(int j, std::vector<double> & grad, double step, bool is_h) {
	double loss;
	if (is_h) {
		this->Update_h(j, grad, step);
		this->UpdateLikelihood();
		loss = this->loglik;
		this->RevUpdate_h(j, grad, step);
	}
	else {
		this->Update_g(j, grad, step);
		this->UpdateLikelihood();
		loss = this->loglik;
		this->RevUpdate_g(j, grad, step);
	}
	return loss;
}

double Gradient::dLossWithUpdate(int j, std::vector<double> & grad, double step, bool is_h) {
	double loss = this->LossWithUpdate(j, grad, step + 0.0001, is_h) - this->LossWithUpdate(j, grad, step - 0.0001, is_h);
	return loss / 0.0002;
}

double Gradient::ComputeStepSize(std::vector<double> & decent, int j, bool update_h) {
	double step_min = 0.001, step_max = 9.999, step_med = 5;
	int itr = 0;
	double d1, d2, d3;
	while (step_max - step_min > 0.001 && itr < 20) {
		d1 = this->dLossWithUpdate(j, decent, step_min, update_h);
		d2 = this->dLossWithUpdate(j, decent, step_max, update_h);
		if (d1 <= 0) {
			step_max = step_min;
			step_min /= 2;
			//break;
		}
		// no need to interpolate
		if (d1 * d2 >= 0) {
			//if (d1 * d1 > d2 * d2) {
			//	step_min = step_max;
				step_max += step_max - step_min;//aggresive search
			//}
			//else {
			//	step_max = step_min;
			//}
			//break;
		}
		step_med = (step_min + step_max) / 2;
		d3 = this->dLossWithUpdate(j, decent, step_med, update_h);
		if (d1 * d3 < 0) {
			step_max = step_med;
			//d2 = d3;
		}
		else {
			step_min = step_med;
			//d1 = d3;
		}
		++itr;
	}
	std::cout << "optimal step size = " << (step_min + step_max) / 2;
	std::cout << "Loss: " << this->LossWithUpdate(j, decent, 1, update_h) << " --> " << this->LossWithUpdate(j, decent, step_med, update_h) << "\n";
	return (step_min + step_max) / 2;
}

double Gradient::UpdateTreeFit(Node * tree, std::vector<std::vector<double>> & x, double alpha, int j, bool update_h, bool line_search) {
	std::vector<double> decent = tree->Predict(x);
	double step_size = 1;
	if (line_search) {
		step_size = this->ComputeStepSize(decent, j, update_h);
	}
	for (int i = 0; i < decent.size(); ++i) {
		decent[i] *= alpha * step_size;
	}
	update_h ? this->Update_h(j, decent) : this->Update_g(j, decent);
	return step_size;
}

void Gradient::RevUpdate_g(int j, std::vector<double> & dg, double alpha) {
	std::vector<double> rev_d;
	for (int i = 0; i < dg.size(); ++i) {
		rev_d.push_back(-dg[i] * alpha);
	}
	this->Update_g(j, rev_d);
}

void Gradient::RevUpdate_h(int j, std::vector<double> & dh, double alpha) {
	std::vector<double> rev_d;
	for (int i = 0; i < dh.size(); ++i) {
		rev_d.push_back(-dh[i] * alpha);
	}
	this->Update_h(j, rev_d);
}

void Gradient::Update_g(int j, std::vector<double> & dg, double alpha) {
	std::vector<double> dg1;
	for (int i = 0; i < dg.size(); ++i) {
		dg1.push_back(dg[i] * alpha);
	}
	this->Update_g(j, dg1);
}

void Gradient::Update_h(int j, std::vector<double> & dh, double alpha) {
	std::vector<double> dh1;
	for (int i = 0; i < dh.size(); ++i) {
		dh1.push_back(dh[i] * alpha);
	}
	this->Update_h(j, dh1);
}