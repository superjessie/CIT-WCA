// =====================================================================================
//
//       Filename:  TupleSet.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  10/28/2014 09:36:40 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Jinkun Lin, jkunlin@gmail.com
//   Organization:  School of EECS, Peking University
//
// =====================================================================================

#ifndef TUPLESET_H
#define TUPLESET_H

#include <vector>

#include "SpecificationFile.h"
#include "Combinadic.h"
#include "Coverage.h"


class TupleSet {
public:
	TupleSet() {};
	void initialize(const SpecificationFile &specificationFile,
			const Coverage &coverage, bool fullfill = false);
	void pop(const unsigned encode);
	void push(const unsigned encode);
	unsigned encode(const unsigned index) { return tupleSet[index]; }
	unsigned size() const { return tupleSet.size(); }
	std::vector<unsigned>::const_iterator begin() const { return tupleSet.begin(); }
	std::vector<unsigned>::const_iterator end() const { return tupleSet.end(); }
#ifndef NDEBUG 
	void print();
#endif
private:
	std::vector<unsigned> tupleSet; //contents of tuple encode
	std::vector<std::vector<unsigned>::size_type> mapping; //encode -> index in tupleSet
};

#endif /* end of include guard: TUPLESET_H */
