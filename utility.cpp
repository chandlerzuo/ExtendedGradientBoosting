#include "stdafx.h"
#include "utility.h"
#include <iostream>

#include <fstream>
#include <sstream>
#include <stdlib.h>
#include<stdio.h>
#include<time.h>
#include<boost/iostreams/device/mapped_file.hpp> // for mmap
#include<boost/interprocess/file_mapping.hpp>
#include<boost/interprocess/mapped_region.hpp>
#include<boost/tokenizer.hpp>

bool CombineConditions(bool r1, bool r2, int mode) {
	switch (mode) {
	case 0:
		return r1 || r2;
	case 1:
		return r1 && r2;
	case 2:
		return r1 || !r2;
	case 3:
		return r1 && !r2;
	default:
		exit;
	}
}

std::vector<std::vector<double>> Transpose(std::vector<std::vector<double>> mat) {
	for (int i = 0; i < mat.size(); ++i) {
		if (mat[i].size() != mat[0].size()) {
			std::cout << "Error: dimension of row " << i << " is " << mat[i].size() << " while dimension of row 0 is " << mat[0].size() << "\n";
			exit;
		}
	}
	std::vector<std::vector<double>> result;
	for (int i = 0; i < mat[0].size(); ++i) {
		std::vector<double> tmp;// = new std::vector<double>();
		for (int j = 0; j < mat.size(); ++j) {
			(tmp).push_back(mat[j][i]);
		}
		result.push_back(tmp);
	}
	return result;
}

void PrintLeadingSpaces(int n) {
	for (int i = 0; i < n; i++) {
		std::cout << "  ";
	}
}

std::vector<int> Sample(int n, double prop) {
	srand(time(NULL));
	int m = (int)(n * prop);
	if (m < 1) {
		m = 1;
	}
	else if (m > n) {
		m = n;
	}
	std::vector<int> arr;
	for (int i = 0; i < n; i++) {
		arr.push_back(i);
	}
	int start_id = 0;
	while (start_id < m) {
		int k = rand() % (n - start_id) + start_id;
		int tmp = arr[start_id];
		arr[start_id] = arr[k];
		arr[k] = tmp;
		start_id++;
	}
	arr.erase(arr.begin() + m, arr.end());
	if (arr.size() != m) {
		std::cout << "Error: the subsampled result vector length " << arr.size() << " does not equal to target length " << m << "\n";
	}
	return arr;
}

void ReadData(std::string filename, int index_response, std::vector<double> *& y, std::vector<std::vector<double>> *& x, bool skipheader) {
	boost::interprocess::file_mapping* fm = new boost::interprocess::file_mapping(filename.c_str(), boost::interprocess::read_only);
	boost::interprocess::mapped_region* region = new boost::interprocess::mapped_region(*fm, boost::interprocess::read_only);
	char* bytes = static_cast<char*>(region->get_address());
	const std::string stringBuffer(bytes, region->get_size());
	std::istringstream data(stringBuffer);
	typedef boost::tokenizer< boost::escaped_list_separator<char> > Tokenizer;
	std::string line;
	std::vector<std::string> parsed;
	int lineno = 0;
	if (!data) {
		std::cout << "File " << filename << " does not exist.\n";
	}
	while (std::getline(data, line)) {
		if (skipheader) {
			skipheader = false;
			continue;
		}
		std::string line_value;
		std::vector<double> line_values;
		std::stringstream ss(line);
		int i = 0;
		while (std::getline(ss, line_value, ',')) {
			double val = std::atof(line_value.c_str());
			if (i == index_response) {
				y->push_back(val);
			}
			else {
				line_values.push_back(val);
			}
			++i;
		}
		x->emplace_back(line_values);
		if (lineno < 10) {
			std::cout << "Record " << lineno << ":";
			if (y->size() > 0) {
				std::cout << "y = " << (*y)[lineno];
			}
			std::cout << " x =";
			for (int j = 0; j < line_values.size(); ++j) {
				std::cout << " " << (*x)[lineno][j];
			}
			std::cout << "\n";
		}
		lineno++;
		if (lineno % 1000 == 0) {
			std::cout << "Read " << lineno << " lines.\n";
		}
	}
	std::cout << "Read " << lineno << " lines of data." << "\n";
}
