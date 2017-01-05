#pragma once
#include <vector>

bool CombineConditions(bool, bool, int);

std::vector<std::vector<double>> Transpose(std::vector<std::vector<double>>);

void PrintLeadingSpaces(int);

void ReadData(std::string, int, std::vector<double> *&, std::vector<std::vector<double>> *&, bool);

std::vector<int> Sample(int, double);

struct less_than_first_coord
{
	inline bool operator() (const std::vector<double> & vec1, const std::vector<double> & vec2)
	{
		return (vec1[0] < vec2[0]);
	}
};
