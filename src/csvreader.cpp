#include "stdafx.h"
#include "csvreader.h"

#include<iostream>
#include<fstream>
#include<ctime>

void ReadCsv(std::string filename, int index_response, std::vector<double> *& y, std::vector<std::vector<double>> *& x, bool skipheader) {
	const clock_t begin_time = clock();

	boost::iostreams::mapped_file_source csv(filename);

	std::vector<std::vector<double>> parsed;
	CsvParser<const char*> grammar;
	std::cout << "Size of file is " << csv.size() << ", " << float(clock() - begin_time) / CLOCKS_PER_SEC << " seconds elapsed.\n";
	if (qi::parse(csv.data(), csv.data() + csv.size(), grammar, parsed)) {
		int nlines = parsed.size();
		int ncols = parsed[0].size();
		int offset = 0;
		if (skipheader) {
			nlines--;
			offset = 1;
		}
		std::cout << "Parsed data of " << nlines << " records and " << ncols << " columns, " << float(clock() - begin_time) / CLOCKS_PER_SEC << " seconds elapsed.\n";
		y->reserve(nlines);
		(*x) = parsed;

/*		x->reserve(nlines);
		for (int i = 0; i < parsed.size(); ++i) {
			x->push_back(parsed[i]);
			if (i % 1000 == 0) {
				std::cout << "Parsed " << i << " lines ..., " << float(clock() - begin_time) / CLOCKS_PER_SEC << " seconds elapsed.\n";
			}
		}*/
		if (index_response >= 0 && index_response < ncols) {
			for (int i = 0; i < parsed.size(); i++) {
				y->push_back((*x)[i][index_response]);
				x->at(i).erase(x->at(i).begin() + index_response);
				if (i > 0 && i % 50000 == 0) {
					std::cout << "Parsed " << i << " lines ..., " << float(clock() - begin_time) / CLOCKS_PER_SEC << " seconds elapsed.\n";
				}
			}
		}
	}
}

void WriteCsv(std::string filename, std::vector<double> & pred) {
	std::ofstream output(filename);
	std::string outstr = "";
	for (int i = 0; i < pred.size(); i++) {
		if (i > 0) {
			outstr += "\n";
		}
		outstr += boost::lexical_cast<std::string>(pred[i]);
		if (i % 10000) {
			output << outstr;
			outstr = "";
		}
	}
	output << outstr;
	output.close();
}