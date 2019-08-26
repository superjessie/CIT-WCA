// =====================================================================================
//
//       Filename:  LocalSearch.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  10/27/2014 02:53:09 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Jinkun Lin, jkunlin@gmail.com
//   Organization:  School of EECS, Peking University
//
// =====================================================================================

#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H


#include "SpecificationFile.h"
#include "ConstraintFile.H"

void localSearch(const SpecificationFile &specificationFile,
		const ConstraintFile &constrFile, const unsigned long long maxTime, int seed);

#endif /* end of include guard: LOCALSEARCH_H */
