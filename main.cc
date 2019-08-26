// =====================================================================================
//
//       Filename:  main.cc
//
//    Description:  
//
//        Version:  1.0
//        Created:  10/27/2014 10:26:54 AM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Jinkun Lin, jkunlin@gmail.com
//   Organization:  School of EECS, Peking University
//
// =====================================================================================

#include <iostream>
#include <string>

#include "SpecificationFile.h"
#include "ConstraintFile.H"
#include "LocalSearch.h"


using namespace std;

int main(int argc, char const *argv[]) {
    cout<<"pre"<<endl;
	if (argc == 0) {
		return 1;
	}
	string modelFile(argv[1]);
	string constrFile;
	unsigned long long maxTime;
	int seed;
   //cout<<argc<<endl;
	if (argc == 5) {
		constrFile = argv[2];
		maxTime = atoi(argv[3]);
		seed = atoi(argv[4]);

	}
	else {
		maxTime = atoi(argv[2]);
		seed = atoi(argv[3]);

	}
   // cout<<"a"<<endl;
	SpecificationFile specificationFile(modelFile);
   // cout<<"b"<<endl;
	ConstraintFile constraintFile(constrFile);
   // cout<<"c"<<endl;
	localSearch(specificationFile, constrFile, maxTime, seed);
   // cout<<"d"<<endl;
	return 0;
}
