#pragma once
#define BOOST_SPIRIT_USE_PHOENIX_V3

#include<vector>
#include<string>
#include<memory>
#include<boost/iostreams/device/mapped_file.hpp> // for mmap
#include<boost/utility/string_ref.hpp>
#include<boost/spirit/include/qi.hpp>
#include<boost/spirit/include/phoenix.hpp>
#include<boost/lexical_cast.hpp>

//typedef boost::string_ref CsvField;
typedef double CsvField;
typedef std::vector<CsvField> CsvLine;
typedef std::vector<CsvLine> CsvFile;  
namespace qi = boost::spirit::qi;

template <typename T> struct CsvParser : qi::grammar<T, CsvFile()> {
	CsvParser() : CsvParser::base_type(lines) {
		using namespace qi;
		using boost::phoenix::construct;
		using boost::phoenix::size;
		using boost::phoenix::begin;
		using boost::spirit::qi::float_;

		//field = raw[*~char_(",\r\n")][_val = construct<CsvField>(begin(qi::_1), size(qi::_1))]; // semantic action
		//field = raw[*~char_(",\r\n")][_val = boost::lexical_cast<double>(construct<boost::string_ref>(begin(qi::_1), size(qi::_1)))]; // semantic action
																								//field = +(~char_(",\r\n"));
		field %= qi::double_;
		line = field % ',';
		lines = line  % eol;
	}
	// declare: line, field, fields
	qi::rule<T, std::vector<std::vector<double>>()> lines;
	qi::rule<T, std::vector<double>()> line;
	qi::rule<T, double()> field;
};

void ReadCsv(std::string, int, std::vector<double> *&, std::vector<std::vector<double>> *&, bool);

void WriteCsv(std::string, std::vector<double> &);

