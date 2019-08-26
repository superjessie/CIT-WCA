// =====================================================================================
//
//       Filename:  CoveringArray.cc
//
//    Description:
//
//        Version:  1.0
//        Created:  10/27/2014 10:03:43 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Jinkun Lin, jkunlin@gmail.com
//   Organization:  School of EECS, Peking University
//
// =====================================================================================

#include "CoveringArray.h"

CoveringArray::CoveringArray(const SpecificationFile &specificationFile,
							 const ConstraintFile &constraintFile, unsigned long long maxT, int seed) : validater(specificationFile), satSolver(constraintFile.isEmpty()), specificationFile(specificationFile),
																										coverage(specificationFile), entryTabu(4), stepIndex(0),
																										maxTime(maxT)
{

	clock_start_a = clock();
	const Options &options = specificationFile.getOptions();
	//add constraint into satSolver
	const std::vector<InputClause> &clauses = constraintFile.getClauses();
	for (unsigned i = 0; i < clauses.size(); ++i)
	{
		satSolver.addClause(const_cast<InputClause &>(clauses[i]));
	}

	const Valid::Formula &formula = constraintFile.getFormula();
#ifndef NVISIBLE
	formula.print();

#endif
	for (auto &c : formula)
	{
		validater.addClause(c);
	}

	option_constrained_indicator.clear();
	option_constrained_indicator.resize(options.size(), false);
	for (auto &c : formula)
	{
		for (auto &lit : c)
		{
			option_constrained_indicator[options.option(lit.variable())] = true;
		}
	}

	for (unsigned option = 0; option < options.size(); ++option)
	{
		if (!option_constrained_indicator[option])
		{
			continue;
		}
		InputClause atLeast;
		for (unsigned j = options.firstSymbol(option), limit = options.lastSymbol(option);
			 j <= limit; ++j)
		{
			atLeast.append(InputTerm(false, j));
		}
		satSolver.addClause(atLeast);
		for (unsigned j = options.firstSymbol(option), limit = options.lastSymbol(option);
			 j <= limit; ++j)
		{
			for (unsigned k = j + 1; k <= limit; ++k)
			{
				InputClause atMost;
				atMost.append(InputTerm(true, j));
				atMost.append(InputTerm(true, k));
				satSolver.addClause(atMost);
			}
		}
	}

	//coverage.initialize(satSolver);
	//uncoveredTuples.initialize(specificationFile, coverage, true);
	coverage.unconstrained_initialize();
	uncoveredTuples.initialize(specificationFile, coverage);
	valueCount = specificationFile.valueCounts();

	mersenne.seed(seed);
}

void CoveringArray::actsInitialize(const std::string file_name)
{
	const Options &options = specificationFile.getOptions();
	const unsigned &strenth = specificationFile.getStrenth();
	std::ifstream res_file(file_name);
	if (!res_file.is_open())
	{
		std::cerr << "file open failed" << std::endl;
		exit(0);
	}
	std::string line;
	while (getline(res_file, line))
	{
		if (line.find("Test Cases") != std::string::npos)
		{
			break;
		}
	}
	while (true)
	{
		bool begin = false;
		while (getline(res_file, line))
		{
			if (line[0] == '1')
			{
				begin = true;
				break;
			}
		}
		if (!begin)
		{
			break;
		}
		array.push_back(std::vector<unsigned>(options.size()));
		std::vector<unsigned> &newRow = *array.rbegin();
		for (unsigned option = 0; option < options.size(); ++option)
		{
			unsigned value;
			std::string value_str(
				line.substr(line.find_last_of('=') + 1, line.size() - 1));
			value = atoi(value_str.c_str());
			newRow[option] = value + options.firstSymbol(option);
			getline(res_file, line);
		}
	}
	res_file.close();

	clock_t start = clock();
	std::vector<size_t> coverByLineIndex(coverage.tupleSize());
	std::vector<unsigned> tuple(strenth);
	for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex)
	{
		auto &line = array[lineIndex];
		for (std::vector<unsigned> columns = combinadic.begin(strenth);
			 columns[strenth - 1] < line.size(); combinadic.next(columns))
		{
			for (unsigned i = 0; i < strenth; ++i)
			{
				tuple[i] = line[columns[i]];
			}
			//			cover(coverage.encode(columns, tuple));
			unsigned encode = coverage.encode(columns, tuple);
			coverByLineIndex[encode] = lineIndex;
			coverage.cover(encode);
		}
	}
	coverage.set_zero_invalid();
	std::cout << "actsInitialize: " << double(clock() - start) / CLOCKS_PER_SEC
			  << std::endl;
	entryTabu.initialize(Entry(array.size(), array.size()));
	validater.initialize(array);
	scoreInitialzeQ();
	tmpPrint();
#ifndef NDEBUG
	int i = 0;
	for (auto &line : array)
	{
		std::cout << i++ << ": ";
		for (auto v : line)
		{
			std::cout << v << ' ';
		}
		std::cout << std::endl;
	}
#endif
}

void CoveringArray::scoreInitialze()
{
	const Options &options = specificationFile.getOptions();
	const unsigned strength = specificationFile.getStrenth();
	unsigned arraySize = array.size();
	Score.resize(arraySize);
	for (unsigned lineIndex = 0; lineIndex < arraySize; lineIndex++)
	{
		Score[lineIndex].resize(valueCount);
		std::vector<unsigned> &line = array[lineIndex];
		unsigned lineSize = line.size();
		for (unsigned i = 0; i < lineSize; i++)
		{
			Score[lineIndex][line[i]] = 0;
			for (unsigned j = options.firstSymbol(i); j <= options.lastSymbol(i); j++)
			{
				if (j == line[i])
					continue;
				Score[lineIndex][j] = varScoreOfRow(j, lineIndex, false);
			}
		}
	}
}

void CoveringArray::scoreInitialzeQ()
{
	const Options &options = specificationFile.getOptions();
	const unsigned strength = specificationFile.getStrenth();
	unsigned arraySize = array.size();
	Score.resize(arraySize);
	for (unsigned lineIndex = 0; lineIndex < arraySize; lineIndex++)
		Score[lineIndex].resize(valueCount, 0);
	unsigned tuplecount = coverage.tupleCount();
	for (unsigned i = 0; i < tuplecount; i++)
	{
		if (coverage.coverCount(i) != 1)
			continue;

		const std::vector<unsigned> &tuple = coverage.getTuple(i);
		const std::vector<unsigned> &columns = coverage.getColumns(i);

		for (unsigned lineIndex = 0; lineIndex < arraySize; lineIndex++)
		{
			std::vector<unsigned> &line = array[lineIndex];
			unsigned diffCount = 0;
			unsigned diffVar;
			bool notcover = false;
			for (unsigned i = 0; i < tuple.size() && !notcover; ++i)
			{
				if (line[columns[i]] != tuple[i])
					notcover = true;
			}
			if (notcover)
				continue;
			for (unsigned i = 0; i < tuple.size() && !notcover; ++i)
			{
				unsigned column = columns[i];
				for (unsigned j = options.firstSymbol(column); j <= options.lastSymbol(column); j++)
				{
					if (j == tuple[i])
						continue;
					Score[lineIndex][j] -= 1;
				}
			}
		}
	}
	printf("check\n");
	/*for (unsigned lineIndex = 0;lineIndex<array.size();lineIndex++)
	{
		for(unsigned v = 0; v < specificationFile.valueCounts();v++)
		{
			if(Score[lineIndex][v]!= varScoreOfRow(v,lineIndex,false))
				printf("wrong score!");
		}
	}*/
}

void CoveringArray::produceSatRow(std::vector<unsigned> &newLine, const unsigned encode)
{
	const unsigned strength = specificationFile.getStrenth();
	const Options &options = specificationFile.getOptions();
	const unsigned width = options.size();
	assert(width == newLine.size());

	InputKnown known;
	const std::vector<unsigned> &ranTuple = coverage.getTuple(encode);
	const std::vector<unsigned> &ranTupleColumns = coverage.getColumns(encode);
	for (unsigned i = 0; i < strength; ++i)
	{
		newLine[ranTupleColumns[i]] = ranTuple[i];
		if (option_constrained_indicator[ranTupleColumns[i]])
		{
			known.append(InputTerm(false, ranTuple[i]));
		}
	}
	std::vector<bool> columnStarted(width, false);
	std::vector<unsigned> columnBases(width);
	for (unsigned column = 0, passing = 0; column < width; ++column)
	{
		if (passing < strength && column == ranTupleColumns[passing])
		{
			passing++;
			continue;
		}
		columnBases[column] = mersenne.next(options.symbolCount(column));
		newLine[column] = options.firstSymbol(column) + columnBases[column] - 1;
	}
	for (long column = 0, passing = 0; column < width; ++column)
	{
		if (passing < strength && column == ranTupleColumns[passing])
		{
			passing++;
			continue;
		}
		const unsigned firstSymbol = options.firstSymbol(column);
		const unsigned lastSymbol = options.lastSymbol(column);
		while (true)
		{
			newLine[column]++;
			if (newLine[column] > lastSymbol)
			{
				newLine[column] = firstSymbol;
			}
			if (newLine[column] == firstSymbol + columnBases[column])
			{
				//backtrack
				if (columnStarted[column])
				{
					columnStarted[column] = false;
					//assign it the value before starting
					newLine[column]--;
					column--;
					while (passing > 0 && column == ranTupleColumns[passing - 1])
					{
						column--;
						passing--;
					}
					//undo column++ of the "for" loop
					column--;
					//the var of parent column is now unabled
					if (option_constrained_indicator[column])
					{
						known.undoAppend();
					}
					break;
				}
				else
				{
					columnStarted[column] = true;
				}
			}
			if (option_constrained_indicator[column])
			{
				known.append(InputTerm(false, newLine[column]));
				if (satSolver(known))
				{
					break;
				}
				known.undoAppend();
			}
			else
			{
				break;
			}
		}
	}
}

void CoveringArray::randomWalk(const unsigned randomSize)
{
	const Options &options = specificationFile.getOptions();
	unsigned strength = specificationFile.getStrenth();
	//printf("%d\n", uncoveredTuples.size());
	unsigned tupleEncode = uncoveredTuples.encode(mersenne.next(uncoveredTuples.size()));
	const std::vector<unsigned> &ranTuple = coverage.getTuple(tupleEncode);
	const std::vector<unsigned> &ranTupleColumns = coverage.getColumns(tupleEncode);
	std::vector<unsigned> validLine;
	std::vector<unsigned> changedVars;

	std::vector<unsigned> validtoChange;
	for (unsigned i = 0; i < array.size(); i++)
	{
		validtoChange.clear();
		unsigned optiontoChange = mersenne.next(options.size());

		unsigned oldVar = array[i][optiontoChange];
		for (unsigned j = options.firstSymbol(optiontoChange); j <= options.lastSymbol(optiontoChange); j++)
		{
			if (oldVar == j)
				continue;
			if (validater.valida_change(i, optiontoChange, oldVar, j))
				validtoChange.push_back(j);
		}
		if (validtoChange.size() != 0)
		{
			unsigned ran = mersenne.next(validtoChange.size());
			replace(validtoChange[ran], i);
		}
	}
	//std::cout<<"getone:"<< getone <<"  getnone:"<< getnone <<std::endl;
	entryTabu.initialize(Entry(array.size(), specificationFile.getOptions().size()));
}

void CoveringArray::removeUselessRows()
{
	//	printf("****************remove use lessrows**************\n");
	const Options &options = specificationFile.getOptions();
	const unsigned strength = specificationFile.getStrenth();
	std::vector<unsigned> tmpTuple(strength);

	for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex)
	{
		bool useless = true;
		for (unsigned j = 0; j < valueCount; j++)
		{
			if (Score[lineIndex][j] != 0)
				useless = false;
		}

		if (useless)
		{
			//printf("useless\n");
			const std::vector<unsigned> &line = array[lineIndex];
			for (std::vector<unsigned> columns = combinadic.begin(strength);
				 columns[strength - 1] < options.size();
				 combinadic.next(columns))
			{
				for (unsigned j = 0; j < strength; ++j)
				{
					tmpTuple[j] = line[columns[j]];
				}
				unsigned encode = coverage.encode(columns, tmpTuple);
				uncover(encode);
				std::vector<unsigned> changeoption;
				updateScoreofTuple(encode, false, lineIndex, changeoption);
			}
			std::swap(Score[lineIndex], Score[Score.size() - 1]);
			std::swap(array[lineIndex], array[array.size() - 1]);
			for (auto &entry : entryTabu)
			{
				if (entry.getRow() == array.size() - 1)
				{
					entry.setRow(lineIndex);
				}
				if (entry.getRow() == lineIndex)
				{
					entry.setRow(array.size() - 1);
				}
			}
			validater.exchange_row(lineIndex, array.size() - 1);
			validater.pop_back_row();
			array.pop_back();
			Score.pop_back();
			--lineIndex;
		}
	}
}

void CoveringArray::removeOneRow()
{
	//	printf("****************remove one row*****************\n");
	const Options &options = specificationFile.getOptions();
	const unsigned strength = specificationFile.getStrenth();
	unsigned minUncoverCount = std::numeric_limits<unsigned>::max();
	std::vector<unsigned> bestRowIndex;
	std::vector<unsigned> tmpTuple(strength);
	for (unsigned i = 0; i < array.size(); ++i)
	{
		unsigned uncoverCount = 0;
		const std::vector<unsigned> &line = array[i];
		for (std::vector<unsigned> columns = combinadic.begin(strength);
			 columns[strength - 1] < options.size();
			 combinadic.next(columns))
		{
			for (unsigned j = 0; j < strength; ++j)
			{
				tmpTuple[j] = line[columns[j]];
			}
			unsigned encode = coverage.encode(columns, tmpTuple);
			if (coverage.coverCount(encode) == 1)
			{
				uncoverCount++;
			}
		}
		if (uncoverCount < minUncoverCount)
		{
			minUncoverCount = uncoverCount;
			bestRowIndex.clear();
			bestRowIndex.push_back(i);
		}
		else if (uncoverCount == minUncoverCount)
		{
			bestRowIndex.push_back(i);
		}
	}

	unsigned rowToremoveIndex = bestRowIndex[mersenne.next(bestRowIndex.size())];
	for (std::vector<unsigned> columns = combinadic.begin(strength);
		 columns[strength - 1] < options.size();
		 combinadic.next(columns))
	{
		for (unsigned j = 0; j < strength; ++j)
		{
			tmpTuple[j] = array[rowToremoveIndex][columns[j]];
		}
		unsigned encode = coverage.encode(columns, tmpTuple);

		uncover(encode);
		std::vector<unsigned> changeoption = columns;
		//printf("encode:%u\n",encode);
		updateScoreofTuple(encode, false, rowToremoveIndex, changeoption);
	}

	std::swap(array[array.size() - 1], array[rowToremoveIndex]);
	std::swap(Score[Score.size() - 1], Score[rowToremoveIndex]);
	validater.exchange_row(rowToremoveIndex, array.size() - 1);
	validater.pop_back_row();
	for (auto &entry : entryTabu)
	{
		if (entry.getRow() == array.size() - 1)
		{
			entry.setRow(rowToremoveIndex);
		}
		if (entry.getRow() == rowToremoveIndex)
		{
			entry.setRow(array.size() - 1);
		}
	}
	array.pop_back();
	Score.pop_back();
}

void CoveringArray::removeOneRowRandom()
{
	//	printf("****************remove one row*****************\n");
	const Options &options = specificationFile.getOptions();
	const unsigned strength = specificationFile.getStrenth();
	std::vector<unsigned> tmpTuple(strength);

	unsigned rowToremoveIndex = mersenne.next(array.size());
	for (std::vector<unsigned> columns = combinadic.begin(strength);
		 columns[strength - 1] < options.size();
		 combinadic.next(columns))
	{
		for (unsigned j = 0; j < strength; ++j)
		{
			tmpTuple[j] = array[rowToremoveIndex][columns[j]];
		}
		unsigned encode = coverage.encode(columns, tmpTuple);

		uncover(encode);
		std::vector<unsigned> changeoption = columns;
		//printf("encode:%u\n",encode);
		updateScoreofTuple(encode, false, rowToremoveIndex, changeoption);
	}

	std::swap(array[array.size() - 1], array[rowToremoveIndex]);
	std::swap(Score[Score.size() - 1], Score[rowToremoveIndex]);
	validater.exchange_row(rowToremoveIndex, array.size() - 1);
	validater.pop_back_row();
	for (auto &entry : entryTabu)
	{
		if (entry.getRow() == array.size() - 1)
		{
			entry.setRow(rowToremoveIndex);
		}
		if (entry.getRow() == rowToremoveIndex)
		{
			entry.setRow(array.size() - 1);
		}
	}
	array.pop_back();
	Score.pop_back();
}

void CoveringArray::updateWeight()
{
	for (unsigned i = 0; i < uncoveredTuples.size(); i++)
	{
		unsigned encode = uncoveredTuples.encode(i);
		coverage.addWeight(encode, 1);
		const std::vector<unsigned> &tuple = coverage.getTuple(encode);
		const std::vector<unsigned> &columns = coverage.getColumns(encode);
		unsigned arraysize = array.size();
		for (unsigned lineIndex = 0; lineIndex < arraysize; lineIndex++)
		{
			std::vector<unsigned> &line = array[lineIndex];
			unsigned diffCount = 0;
			unsigned diffVar;
			for (unsigned j = 0; j < tuple.size(); j++)
			{
				if (line[columns[j]] != tuple[j])
				{
					//printf("bingo!");
					diffCount++;
					diffVar = tuple[j];
				}
				if (diffCount > 1)
					break;
			}
			if (diffCount == 1)
			{
				Score[lineIndex][diffVar] += 1;
			}
		}
	}
}

void CoveringArray::optimize() //propability random
{
	long long tabucount = 0;
	std::vector<std::vector<unsigned>> bestArray; // = array;
	long long bestScoreforNow = coverage.tupleCount() - uncoveredTuples.size();
	long long tempScore;
	const Options &options = specificationFile.getOptions();
	bool tabuSuc;
	printf("start\n");
	while (true)
	{
		if ((double)(clock() - clock_start_a) / CLOCKS_PER_SEC > maxTime)
		{
			break;
		}
		if (uncoveredTuples.size() == 0)
		{
			bestArray = array;
			tmpPrint();
			removeUselessRows();
			removeOneRowRandom();
			bestScoreforNow = coverage.tupleCount() - uncoveredTuples.size();
		}

		tabuSuc = tabuStep();
		tempScore = coverage.tupleCount() - uncoveredTuples.size();
		if (tempScore > bestScoreforNow)
		{
			bestScoreforNow = tempScore;
			stepIndex = 0;
		}
		stepIndex++;
		tabucount++;
		continue;
	}

	printf("end\n");
	//printf("tabucount:%lld\n", tabucount);
	if (uncoveredTuples.size() == 0)
	{

		removeUselessRows();
		bestArray = array;
		tmpPrint(tabucount);
	}

	if (!verify(bestArray))
	{
		std::cout << "wrong answer!!!!!" << std::endl;
		return;
	}

	// #ifndef NDEBUG
	std::cerr << "********Debuging CoveringArray::optimize*********" << std::endl;
	std::cerr << "printing bestArray..." << std::endl;
	for (unsigned i = 0; i < bestArray.size(); ++i)
	{
		std::cerr << i << "th  ";
		for (auto x : bestArray[i])
		{
			std::cerr << ' ' << x;
		}
		std::cerr << std::endl;
	}
	std::cerr << "total size : " << bestArray.size() << std::endl;
	std::cerr << "********End of Debuing CoveringArray::optimize********" << std::endl;
	// #endif
}

bool CoveringArray::tabuStep() //first best neighbor without worse acception with random
{
	if (stepIndex > 1000)
	{
		//printf("random\n");
		randomWalk(1);
		stepIndex = 0;
	}
	if (uncoveredTuples.size() == 0)
		return true;
	const Options &options = specificationFile.getOptions();
	unsigned uncoversize = uncoveredTuples.size();
	unsigned base = mersenne.next(uncoversize);

	bool candidateIsNull = true;
	for (unsigned t = 0; t < uncoversize; t++)
	{
		const unsigned tupleEncode = uncoveredTuples.encode((t + base) % uncoversize);
		std::vector<unsigned> bestRows;
		std::vector<unsigned> bestVars;
		long long bestScore = std::numeric_limits<long long>::min();
		const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
		const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);

		for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex)
		{
			std::vector<unsigned> &line = array[lineIndex];
			unsigned diffCount = 0;
			unsigned diffVar;
			for (unsigned i = 0; i < tuple.size(); ++i)
			{
				if (line[columns[i]] != tuple[i])
				{
					//printf("bingo!");
					diffCount++;
					diffVar = tuple[i];
				}
			}
			if (diffCount > 1)
			{
				continue;
			}			
			unsigned diffOption = options.option(diffVar);
			//Tabu
			if (entryTabu.isTabu(Entry(lineIndex, diffOption)))
			{
				continue;
			}
			
			if (!validater.valida_change(lineIndex, diffOption, line[diffOption], diffVar))
			{
				continue;
			}
			candidateIsNull = false;
			//long long realScore = varScoreOfRow(diffVar, lineIndex,false);
			long long tmpScore = Score[lineIndex][diffVar];
			//if(realScore == tmpScore) printf("rightscore!\n");
			if (tmpScore <= 0)
				continue;
			if (bestScore < tmpScore)
			{
				bestScore = tmpScore;
				bestRows.clear();
				bestRows.push_back(lineIndex);
				bestVars.clear();
				bestVars.push_back(diffVar);
			}
			else if (bestScore == tmpScore)
			{
				bestRows.push_back(lineIndex);
				bestVars.push_back(diffVar);
			}
		}
		if (bestRows.size() != 0 )
		{
			unsigned ran = mersenne.next(bestRows.size());
			replace(bestVars[ran], bestRows[ran]);
			return true;
		}
	}
	updateWeight();
	if (candidateIsNull)
	{
		for (unsigned t = 0; t < uncoversize; t++)
		{
			const unsigned tupleEncode = uncoveredTuples.encode((t + base) % uncoversize);
			std::vector<unsigned> bestRows;
			std::vector<unsigned> bestVars;
			long long bestScore = std::numeric_limits<long long>::min();
			const std::vector<unsigned> &tuple = coverage.getTuple(tupleEncode);
			const std::vector<unsigned> &columns = coverage.getColumns(tupleEncode);

			for (unsigned lineIndex = 0; lineIndex < array.size(); ++lineIndex)
			{
				std::vector<unsigned> &line = array[lineIndex];
				std::vector<unsigned> changeVars;
				bool isTabu = false;
				for (unsigned i = 0; i < tuple.size(); ++i)
				{
					if (line[columns[i]] != tuple[i])
					{
						//printf("bingo!");
						changeVars.push_back(tuple[i]);
						if (entryTabu.isTabu(Entry(lineIndex, columns[i])))
							isTabu = true;
					}
				}
				if (isTabu)
					continue;

				bool need_to_check = false;
				for (auto v : changeVars)
				{
					const Options &options = specificationFile.getOptions();
					if (option_constrained_indicator[options.option(v)])
					{
						need_to_check = true;
						break;
					}
				}

				// my check
				if (need_to_check && !validater.valida_row(array[lineIndex], changeVars))
				{
					continue;
				}

				for (auto v : changeVars)
				{
					replace(v, lineIndex);
				}
				return true;
			}
		}
	}

	return false;
}

long long CoveringArray::varScoreOfRow(const unsigned var, const unsigned lineIndex, bool flag)
{
	//printf("start cal");
	const Options &options = specificationFile.getOptions();
	const unsigned strength = specificationFile.getStrenth();
	std::vector<unsigned> &line = array[lineIndex];
	const unsigned varOption = options.option(var);
	if (line[varOption] == var)
	{
		return 0;
	}
	std::swap(line[line.size() - 1], line[varOption]);

	long long coverChangeCount = 0;
	std::vector<unsigned> tmpSortedColumns(strength);
	std::vector<unsigned> tmpSortedTupleToCover(strength);
	std::vector<unsigned> tmpSortedTupleToUncover(strength);
	for (std::vector<unsigned> columns = combinadic.begin(strength - 1);
		 columns[strength - 2] < line.size() - 1;
		 combinadic.next(columns))
	{
		for (unsigned i = 0; i < columns.size(); ++i)
		{
			tmpSortedTupleToUncover[i] = tmpSortedTupleToCover[i] = line[columns[i]];
		}
		tmpSortedTupleToCover[strength - 1] = var;
		tmpSortedTupleToUncover[strength - 1] = line[line.size() - 1];
		std::sort(tmpSortedTupleToCover.begin(), tmpSortedTupleToCover.end());
		std::sort(tmpSortedTupleToUncover.begin(), tmpSortedTupleToUncover.end());
		for (unsigned i = 0; i < tmpSortedTupleToCover.size(); ++i)
		{
			tmpSortedColumns[i] = options.option(tmpSortedTupleToCover[i]);
		}
		unsigned tmpTupleToCoverEncode = coverage.encode(tmpSortedColumns, tmpSortedTupleToCover);
		unsigned tmpTupleToUncoverEncode = coverage.encode(tmpSortedColumns, tmpSortedTupleToUncover);
		if (coverage.coverCount(tmpTupleToCoverEncode) == 0)
		{
			if (flag)
			{
				printf("+++:");
				printTuple(tmpTupleToCoverEncode);
				//printf("covercode:%u\n",tmpTupleToCoverEncode);
				//printf("the different stuple:%u\n",coverage.coverCount(tmpTupleToCoverEncode));
			}
			coverChangeCount += coverage.getWeight(tmpTupleToCoverEncode);
			;
		}
		if (coverage.coverCount(tmpTupleToUncoverEncode) == 1)
		{
			if (flag)
			{
				printf("---:");
				printTuple(tmpTupleToUncoverEncode);
			}
			coverChangeCount -= coverage.getWeight(tmpTupleToUncoverEncode);
			;
		}
	}
	std::swap(line[line.size() - 1], line[varOption]);

	return coverChangeCount;
}

void CoveringArray::printTuple(unsigned encode)
{
	const std::vector<unsigned> &tuple = coverage.getTuple(encode);
	for (unsigned i = 0; i < tuple.size(); i++)
	{
		printf("%u ", tuple[i]);
	}
	printf("\n");
}

//quite similar to varScoreOfRow function
void CoveringArray::replace(const unsigned var, const unsigned lineIndex)
{

	//std::cout << "line:" << lineIndex << " change to: " << var << std::endl;
	const Options &options = specificationFile.getOptions();
	const unsigned strength = specificationFile.getStrenth();
	std::vector<unsigned> &line = array[lineIndex];
	const unsigned varOption = options.option(var);
	const unsigned oldvar = line[varOption];

	entryTabu.insert(Entry(lineIndex, varOption));

	if (line[varOption] == var)
	{
		return;
	}
	std::vector<unsigned> changeoption;
	changeoption.push_back(varOption);

	validater.change_var(lineIndex, varOption, line[varOption], var);

	std::swap(line[line.size() - 1], line[varOption]);

	std::vector<unsigned> tmpSortedColumns(strength);
	std::vector<unsigned> tmpSortedTupleToCover(strength);
	std::vector<unsigned> tmpSortedTupleToUncover(strength);
	for (std::vector<unsigned> columns = combinadic.begin(strength - 1);
		 columns[strength - 2] < line.size() - 1;
		 combinadic.next(columns))
	{
		for (unsigned i = 0; i < columns.size(); ++i)
		{
			tmpSortedTupleToUncover[i] = tmpSortedTupleToCover[i] = line[columns[i]];
		}
		tmpSortedTupleToCover[strength - 1] = var;
		tmpSortedTupleToUncover[strength - 1] = line[line.size() - 1];
		std::sort(tmpSortedTupleToCover.begin(), tmpSortedTupleToCover.end());
		std::sort(tmpSortedTupleToUncover.begin(), tmpSortedTupleToUncover.end());
		for (unsigned i = 0; i < tmpSortedTupleToCover.size(); ++i)
		{
			tmpSortedColumns[i] = options.option(tmpSortedTupleToCover[i]);
		}
		unsigned tmpTupleToCoverEncode = coverage.encode(tmpSortedColumns, tmpSortedTupleToCover);
		unsigned tmpTupleToUncoverEncode = coverage.encode(tmpSortedColumns, tmpSortedTupleToUncover);

		//need not check coverCount, cover(encode) will do this

		//printf("do uncover\n");
		//printTuple(tmpTupleToUncoverEncode);
		uncover(tmpTupleToUncoverEncode);
		std::swap(line[line.size() - 1], line[varOption]);
		updateScoreofTuple(tmpTupleToUncoverEncode, false, lineIndex, changeoption);
		std::swap(line[line.size() - 1], line[varOption]);

		//printf("do cover\n");
		//printTuple(tmpTupleToCoverEncode);
		cover(tmpTupleToCoverEncode);
		//printTuple(tmpTupleToCoverEncode);
		std::swap(line[line.size() - 1], line[varOption]);
		updateScoreofTuple(tmpTupleToCoverEncode, true, lineIndex, changeoption);
		std::swap(line[line.size() - 1], line[varOption]);
	}

	std::swap(line[line.size() - 1], line[varOption]);
	updateScoreofOldRow(lineIndex, oldvar, oldvar);
	line[varOption] = var;
	updateScoreofNewRow(lineIndex, var, oldvar);
	//printf("Replace in line:%u, in option:%u, in var:%u ,oldvar:%u\n", lineIndex, varOption, var, oldvar);
}

void CoveringArray::updateScoreofTuple(unsigned encode, bool tobecover, unsigned changeline, std::vector<unsigned> &changeoption)
{

	int weight = coverage.getWeight(encode);
	int covercount = coverage.coverCount(encode);
	if (covercount > 2 || covercount < 0)
		return;
	unsigned arraySize = array.size();
	const Options &options = specificationFile.getOptions();

	if (tobecover)
	{
		if (covercount < 1 || covercount > 2)
		{
			return;
		}
		const std::vector<unsigned> &tuple = coverage.getTuple(encode);
		const std::vector<unsigned> &columns = coverage.getColumns(encode);
		if (covercount == 1)
		{
			//printf("cover:covercount = 1\n");
			for (unsigned lineIndex = 0; lineIndex < arraySize; lineIndex++)
			{
				std::vector<unsigned> &line = array[lineIndex];
				unsigned diffCount = 0;
				unsigned diffVar;
				//
				if (lineIndex == changeline)
				{
					for (unsigned i = 0; i < tuple.size(); ++i)
					{
						std::vector<unsigned>::iterator it = find(changeoption.begin(), changeoption.end(), columns[i]);
						for (unsigned j = options.firstSymbol(columns[i]); j <= options.lastSymbol(columns[i]); j++)
						{
							//printf("replace to be cover:%u\n", encode);
							if (it == changeoption.end() && j == line[columns[i]])
							{
								//printf("release j:%u\n", j);
								continue;
							}
							Score[lineIndex][j] -= weight;
						}
					}
				}
				else
				{
					for (unsigned i = 0; i < tuple.size(); ++i)
					{
						if (line[columns[i]] != tuple[i])
						{
							diffVar = tuple[i];
							diffCount++;
						}
					}
					if (diffCount != 1)
						continue;
					Score[lineIndex][diffVar] -= weight;

					//printf("lineIndex:%u,changeline:%u,diffCount:%u\n", lineIndex, changeline, diffCount);
				}
			}
		}

		if (covercount == 2)
		{
			//printf("cover:covercount = 2\n");
			for (unsigned lineIndex = 0; lineIndex < arraySize; lineIndex++)
			{
				if (lineIndex == changeline)
				{
					continue;
				}
				std::vector<unsigned> &line = array[lineIndex];
				unsigned diffCount = 0;
				unsigned diffVar;
				for (unsigned i = 0; i < tuple.size(); ++i)
				{
					if (line[columns[i]] != tuple[i])
					{
						diffVar = tuple[i];
						diffCount++;
					}
				}
				if (diffCount != 0)
					continue;

				//score to othervalue --
				for (unsigned i = 0; i < tuple.size(); ++i)
				{
					for (unsigned j = options.firstSymbol(columns[i]); j <= options.lastSymbol(columns[i]); j++)
					{
						if (j == tuple[i])
							continue;
						Score[lineIndex][j] += weight;
					}
				}
			}
		}
	}
	else
	{
		const std::vector<unsigned> &tuple = coverage.getTuple(encode);
		const std::vector<unsigned> &columns = coverage.getColumns(encode);
		if (covercount > 1)
		{
			return;
		}
		if (covercount == 0)
		{
			//printf("uncover:covercount = 0\n");
			for (unsigned lineIndex = 0; lineIndex < arraySize; lineIndex++)
			{
				std::vector<unsigned> &line = array[lineIndex];
				unsigned diffCount = 0;
				unsigned diffVar;
				if (lineIndex == changeline)
				{
					for (unsigned i = 0; i < tuple.size(); ++i)
					{
						std::vector<unsigned>::iterator it = find(changeoption.begin(), changeoption.end(), columns[i]);
						for (unsigned j = options.firstSymbol(columns[i]); j <= options.lastSymbol(columns[i]); j++)
						{
							//printf("option:%u, value: %u\n", columns[i], j);
							//printf("ischanged:%d\n", it == changeoption.end());
							//printf("changeoption is null:%d\n", changeoption.size() == 0);

							if (it == changeoption.end() && j == line[columns[i]])
							{
								continue;
							}
							//printf("option:%u, value: %u\n\n", columns[i], j);
							Score[lineIndex][j] += weight;
						}
					}
				}
				else
				{
					for (unsigned i = 0; i < tuple.size(); ++i)
					{
						if (line[columns[i]] != tuple[i])
						{
							diffVar = tuple[i];
							diffCount++;
						}
					}
					if (diffCount != 1)
						continue;
					//printf("diffCount == 1");
					Score[lineIndex][diffVar] += weight;
				}
			}
		}
		if (covercount == 1)
		{
			//printf("uncover:covercount = 1\n");
			for (unsigned lineIndex = 0; lineIndex < arraySize; lineIndex++)
			{
				std::vector<unsigned> &line = array[lineIndex];
				unsigned diffCount = 0;
				unsigned diffVar;
				if (lineIndex == changeline)
				{
					continue;
				}
				for (unsigned i = 0; i < tuple.size(); ++i)
				{
					if (line[columns[i]] != tuple[i])
					{
						diffVar = tuple[i];
						diffCount++;
					}
				}
				//printf("diffvar:%d\n",diffCount);
				if (diffCount != 0)
					continue;

				//printf("find the one\n");
				//score to othervalue --
				for (unsigned i = 0; i < tuple.size(); ++i)
				{
					std::vector<unsigned>::iterator it = find(changeoption.begin(), changeoption.end(), columns[i]);
					for (unsigned j = options.firstSymbol(columns[i]); j <= options.lastSymbol(columns[i]); j++)
					{
						if (j == tuple[i])
						{
							continue;
							//printf("release j:%u\n", j);
						}
						Score[lineIndex][j] -= weight;
					}
				}
			}
		}
	}
}

void CoveringArray::updateScoreofOldRow(unsigned lineIndex, unsigned var, unsigned oldvar)
{

	const Options &options = specificationFile.getOptions();
	const unsigned strength = specificationFile.getStrenth();
	std::vector<unsigned> &line = array[lineIndex];
	const unsigned varOption = options.option(var);
	for (unsigned linei = 0; linei < uncoveredTuples.size(); linei++)
	{
		unsigned encode = uncoveredTuples.encode(linei);
		int weight = coverage.getWeight(encode);
		const std::vector<unsigned> &tuple = coverage.getTuple(encode);
		const std::vector<unsigned> &columns = coverage.getColumns(encode);
		std::vector<unsigned> columnss;
		columnss.assign(columns.begin(), columns.end());
		unsigned diffCount = 0;
		unsigned diffVar, diffOption;
		for (unsigned i = 0; i < tuple.size(); ++i)
		{
			if (line[columns[i]] != tuple[i])
			{
				diffVar = tuple[i];
				diffOption = options.option(diffVar);
				diffCount++;
			}
		}
		if (diffCount == 1)
		{
			std::vector<unsigned>::iterator it = find(columnss.begin(), columnss.end(), varOption);
			if (it != columnss.end() && diffOption != varOption)
				Score[lineIndex][diffVar] -= weight;
		}
	}
}

void CoveringArray::updateScoreofNewRow(unsigned lineIndex, unsigned var, unsigned oldvar)
{
	const Options &options = specificationFile.getOptions();
	const unsigned strength = specificationFile.getStrenth();
	std::vector<unsigned> &line = array[lineIndex];
	const unsigned varOption = options.option(var);
	for (unsigned linei = 0; linei < uncoveredTuples.size(); linei++)
	{
		unsigned encode = uncoveredTuples.encode(linei);
		int weight = coverage.getWeight(encode);
		const std::vector<unsigned> &tuple = coverage.getTuple(encode);
		const std::vector<unsigned> &columns = coverage.getColumns(encode);
		std::vector<unsigned> columnss;
		columnss.assign(columns.begin(), columns.end());
		unsigned diffCount = 0;
		unsigned diffVar, diffOption;
		for (unsigned i = 0; i < tuple.size(); ++i)
		{
			if (line[columns[i]] != tuple[i])
			{
				diffVar = tuple[i];
				diffOption = options.option(diffVar);
				diffCount++;
			}
		}
		if (diffCount == 1)
		{
			std::vector<unsigned>::iterator it = find(columnss.begin(), columnss.end(), varOption);
			if (it != columnss.end() && diffOption != varOption)
				Score[lineIndex][diffVar] += weight;
		}
	}
}

void CoveringArray::cover(const unsigned encode)
{
	coverage.cover(encode);
	unsigned coverCount = coverage.coverCount(encode);
	if (coverCount == 1)
	{
		uncoveredTuples.pop(encode);
	}
}

void CoveringArray::uncover(const unsigned encode)
{
	coverage.uncover(encode);
	unsigned coverCount = coverage.coverCount(encode);
	if (coverCount == 0)
	{
		uncoveredTuples.push(encode);
	}
}

void CoveringArray::tmpPrint(long long tabucount)
{

	std::cout << (double)(clock() - clock_start_a) / CLOCKS_PER_SEC << '\t' << array.size() << '\t' << tabucount << std::endl;
}

void CoveringArray::tmpPrint()
{

	std::cout << (double)(clock() - clock_start_a) / CLOCKS_PER_SEC << '\t' << array.size() << std::endl;
}

void CoveringArray::printArray()
{
	for (unsigned i = 0; i < array.size(); i++)
	{
		for (unsigned j = 0; j < array[i].size(); j++)
			std::cout << array[i][j] << " ";
		std::cout << std::endl;
	}
}
void CoveringArray::printScore()
{
	std::cout<<"print score"<<std::endl;
	for (unsigned i = 0; i < Score.size(); i++)
	{
		for (unsigned j = 0; j < Score[i].size(); j++)
			std::cout<<Score[i][j]<<" ";
		std::cout<<std::endl;
	}
}
void CoveringArray::printUncovered()
{
	std::cout << "uncovered tuples:" << std::endl;
	for (unsigned i = 0; i < uncoveredTuples.size(); i++)
	{
		unsigned encode = uncoveredTuples.encode(i);
		printTuple(encode);
	}
}

bool CoveringArray::verify(const std::vector<std::vector<unsigned>> &resultArray)
{
	const unsigned strength = specificationFile.getStrenth();
	const Options &options = specificationFile.getOptions();
	Coverage tmpCoverage(specificationFile);
	tmpCoverage.initialize(satSolver);
	std::vector<unsigned> tuple(strength);
	unsigned lineIndex = 0;
	for (auto &line : resultArray)
	{
		for (unsigned column = 0; column < line.size(); ++column)
		{
			if (line[column] < options.firstSymbol(column) ||
				line[column] > options.lastSymbol(column))
			{
				std::cerr << "error line: " << lineIndex;
				std::cerr << " option: " << column << std::endl;
				std::cerr << "should be " << options.firstSymbol(column) << " <= var <= " << options.lastSymbol(column) << std::endl;
				for (auto var : line)
				{
					std::cerr << var << ' ';
				}
				std::cerr << std::endl;
				return false;
			}
		}
		for (std::vector<unsigned> columns = combinadic.begin(strength);
			 columns[strength - 1] < line.size();
			 combinadic.next(columns))
		{
			for (unsigned i = 0; i < strength; ++i)
			{
				tuple[i] = line[columns[i]];
			}
			unsigned encode = tmpCoverage.encode(columns, tuple);
			if (tmpCoverage.coverCount(encode) < 0)
			{
				std::cerr << "violate constraints" << std::endl;
				return false;
			}
			tmpCoverage.cover(encode);
		}
		++lineIndex;
	}
	return tmpCoverage.allIsCovered();
}
