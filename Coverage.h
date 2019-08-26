// =====================================================================================
//
//       Filename:  Coverage.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  10/27/2014 09:19:01 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Jinkun Lin, jkunlin@gmail.com
//   Organization:  School of EECS, Peking University
//
// =====================================================================================

#ifndef COVERAGE_H
#define COVERAGE_H


#include <vector>
#include <algorithm>

#include "SpecificationFile.h"
#include "PascalTriangle.h"
#include "Combinadic.h"
#include "SAT.H"

class Coverage {
public:
	Coverage (const SpecificationFile &specificationFile);
	void unconstrained_initialize();
	void set_zero_invalid();
	void initialize(SATSolver &satSovler);
	int coverCount(const unsigned encode) const { return contents[encode]; }
	void cover(const unsigned encode) { ++contents[encode]; }
	void uncover(const unsigned encode) { --contents[encode]; }
	unsigned encode(const std::vector<unsigned> &sortedColumns,
			const std::vector<unsigned> &sortedSubset);
	const std::vector<unsigned>& getColumns(const unsigned encode) const;
	const std::vector<unsigned>& getTuple(unsigned encode) const { return tuples[encode]; }
	size_t tupleSize() { return contents.size(); }
	bool allIsCovered();

	void addWeight(unsigned encode,unsigned step){weights[encode]+= step;}
	unsigned getWeight(unsigned encode){return weights[encode];}
	void resetWeight()
	{
		weights.clear();
		weights.resize(contents.size(),1);
	}


	unsigned tupleCount() const { return contents.size(); }
	void error() {
		for (unsigned i = 0; i < contents.size(); ++i) {
			if (contents[i] < 0) {
				std::cerr << "encode: " << i << "count: " << contents[i] << std::endl; 
				abort();
			}
		}
	}
	void print();
private:
	const SpecificationFile &specificationFile;
	std::vector<unsigned> offsets;
	std::vector<std::vector<unsigned>> columns;
	std::vector<int> contents;
	std::vector<std::vector<unsigned>> tuples;
	std::vector<int> weights;
};

#endif /* end of include guard: COVERAGE_H */
