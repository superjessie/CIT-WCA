// =====================================================================================
//
//       Filename:  TupleSet.cc
//
//    Description:  
//
//        Version:  1.0
//        Created:  10/28/2014 09:50:07 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Jinkun Lin, jkunlin@gmail.com
//   Organization:  School of EECS, Peking University
//
// =====================================================================================

#include "TupleSet.h"

void TupleSet::initialize(const SpecificationFile &specificationFile,
		const Coverage &coverage, bool fullfill) {

	const Options &options = specificationFile.getOptions();
	unsigned strength = specificationFile.getStrenth();

	unsigned MaxSize = 0;
	for (std::vector<unsigned> columns = combinadic.begin(strength);
			columns[strength - 1] < options.size(); combinadic.next(columns)) {
		unsigned blockSize = 1;
		for (unsigned i = 0; i < strength; ++i) {
			blockSize *= options.symbolCount(columns[i]);
		}
		MaxSize += blockSize;
	}
	mapping.resize(MaxSize);
	if (fullfill) {
		for (unsigned encode = 0, i = 0; encode < coverage.tupleCount(); ++encode) {
			if (coverage.coverCount(encode) != -1) {
				mapping[encode] = i++;
				tupleSet.push_back(encode);
			}
		}
	}
#ifndef NDEBUG 
	std::cout << "********Debuging TupleSet********" << std::endl;
	print();
	std::cout << "********End of debuging TupleSet********" << std::endl;
#endif
}

void TupleSet::push(const unsigned encode) {
	mapping[encode] = tupleSet.size();
	tupleSet.push_back(encode);
}

void TupleSet::pop(const unsigned encode) {
	tupleSet[mapping[encode]] = tupleSet[tupleSet.size() - 1];
	mapping[tupleSet[tupleSet.size() - 1]] = mapping[encode];
	tupleSet.pop_back();
}

#ifndef NDEBUG 
void TupleSet::print() {
	std::cout << "mapping : " << std::endl;
	for (auto i : mapping) {
		std::cout << i << ' ';
	}
	std::cout << std::endl << "mapping size = " << mapping.size() << std::endl;
	std::cout << "tupleSet : " << std::endl;
	for (auto x : tupleSet) {
		std::cout << x << ' ';
	}
	std::cout << std::endl << "tupleSet size = " << tupleSet.size() << std::endl;
}
#endif
