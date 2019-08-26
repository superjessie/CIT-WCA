// =====================================================================================
//
//       Filename:  CoveringArray.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  10/27/2014 09:18:19 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Jinkun Lin, jkunlin@gmail.com
//   Organization:  School of EECS, Peking University
//
// =====================================================================================

#include <vector>
#include <set>
#include <queue>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

#include "ConstraintFile.H"
#include "Coverage.h"
#include "TupleSet.h"
#include "mersenne.h"
#include "Tabu.h"
#include "SAT.H"

class CoveringArray {
public:
	CoveringArray (const SpecificationFile &specificationFile,
			const ConstraintFile &constraintFile, unsigned long long maxT, int seed);
	void greedyConstraintInitialize();
	void actsInitialize(const std::string file_name);
	void optimize();


private:
	Valid::Validater validater;
	SATSolver satSolver;
	std::vector<bool> option_constrained_indicator;
	Mersenne mersenne;
	const SpecificationFile &specificationFile;
	std::vector<std::vector<unsigned>> array;
	Coverage coverage;
	TupleSet uncoveredTuples;
	std::set<unsigned> varInUncovertuples;
	Tabu<Entry> entryTabu;
	unsigned stepIndex;
	std::vector<std::vector<long long>> Score;
	unsigned valueCount;


	unsigned long long maxTime;
	clock_t clock_start_a;
	
	void printTuple(unsigned Encode);
	void printArray();
	void printScore();
	void printUncovered();
	void scoreInitialze();
	void scoreInitialzeQ();
	void cover(const unsigned encode);
	void uncover(const unsigned encode);
	void updateWeight();
	void updateScoreofTuple(unsigned encode,bool tobecover,unsigned changeline,std::vector<unsigned>& changedoption);
	void updateScoreofOldRow(unsigned line,unsigned diffVar,unsigned oldVar);
	void updateScoreofNewRow(unsigned line,unsigned diffVar,unsigned oldVar);
	//produce one row at least cover one uncovered tuple.
	//Producing the row without update coverage
	void produceSatRow(std::vector<unsigned> &newLine, const unsigned encode);
	//greedily produce one row at least cover one uncovered tuple.
	//producing the row AND updating coverage
	void mostGreedySatRow(std::vector<unsigned> &newLine, const unsigned encode);
	void replaceRow(const unsigned lineIndex, const unsigned encode);
	void randomWalk(unsigned random_size);
	void removeUselessRows();
	void removeOneRow();
	void removeOneRowRandom();
	long long varScoreOfRow(const unsigned var, const unsigned lineIndex,bool flag);
	void replace(const unsigned var, const unsigned lineIndex);
	

	bool tabuStep();


	void tmpPrint(long long tabucount);
	void tmpPrint();
	bool verify(const std::vector<std::vector<unsigned>> &resultArray);
#ifndef NDEBUG 
	void print();
	
#endif
	void t();
};

