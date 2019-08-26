// Copyright 2008, 2009 Brady J. Garvin

// This file is part of Covering Arrays by Simulated Annealing (CASA).

// CASA is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// CASA is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with CASA.  If not, see <http://www.gnu.org/licenses/>.

#include "Combinadic.h"
#include <cassert>


std::vector<unsigned> Combinadic::begin(const unsigned size) const {
	std::vector<unsigned> vec(size);
	for (unsigned i = 0; i < size; ++i) {
		vec[i] = i;
	}
	return std::move(vec);
}

void Combinadic::next(std::vector<unsigned> &sortedSubset) const {
	assert(sortedSubset.size());
	unsigned limit = sortedSubset.size() - 1, ceiling = sortedSubset[0];
	for (unsigned i = 0; i < limit; ++i) {
		unsigned entry = ceiling + 1;
		ceiling = sortedSubset[i + 1];
		if (entry < ceiling) {
			sortedSubset[i] = entry;
			return;
		}
		sortedSubset[i] = i;
	}
	++sortedSubset[limit];
}

void Combinadic::previous(std::vector<unsigned> &sortedSubset) const {
	assert(sortedSubset.size());
	unsigned limit = sortedSubset.size();
	for(unsigned i = 0; i < limit; ++i) {
		unsigned entry = sortedSubset[i];
		if (entry > i) {
			do {
				sortedSubset[i] = --entry;
			} while (i-- > 0);
			return;
		}
	}
}

unsigned Combinadic::encode(const std::vector<unsigned> &sortedSubset) const {
	unsigned result = 0;
	for (unsigned i = 0; i < sortedSubset.size(); ++i) {
		result += pascalTriangle.nCr(sortedSubset[i], i+1);
	}
	return result;
}

Combinadic combinadic;
