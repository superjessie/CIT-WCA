// =====================================================================================
//
//       Filename:  Options.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  10/27/2014 10:28:19 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Jinkun Lin, jkunlin@gmail.com
//   Organization:  School of EECS, Peking University
//
// =====================================================================================

#ifndef OPTIONS_H
#define OPTIONS_H

#include <vector>
#include <iostream>
class Options {
public:
	Options () {}
	unsigned initialize(const std::vector<unsigned> &values);
	unsigned size() const { return cumulativeValueCounts.size(); }
	unsigned option(const unsigned symbol) const { return owingOptions[symbol]; };
	unsigned firstSymbol(const unsigned option) const {
		return option ? cumulativeValueCounts[option -1] : 0;
	}
	unsigned lastSymbol(const unsigned option) const {
		return cumulativeValueCounts[option] - 1;
	}
	unsigned symbolCount(const unsigned option) const {
		return option ? cumulativeValueCounts[option] - cumulativeValueCounts[option - 1] : cumulativeValueCounts[option];
	}
#ifndef NDEBUG 
	void print();
#endif
private:
	std::vector<unsigned> cumulativeValueCounts;
	std::vector<unsigned> owingOptions;
};

#endif /* end of include guard: OPTIONS_H */

