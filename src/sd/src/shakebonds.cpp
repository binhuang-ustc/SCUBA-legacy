/*
 * shakebonds.cpp
 *
 *  Created on: 2017年11月2日
 *      Author: hyliu
 */
#include "sd/shakebonds.h"
#include "sd/forcefield.h"

#include <cmath>
using namespace NSPsd;

bool ShakeBonds::shake(const std::vector<double> & masses,
		const std::vector<double> &crdref, std::vector<double> &crd) const {

	if(!ison_) return true;

	int natoms = crd.size() / ndim_;
	std::vector<bool> moved(natoms, true);
	std::vector<bool> moved_next(natoms, false);
	double tolerance = 1.e-4;
	bool done = false;
	while (!done) {
		done = true;
		for (int ib = 0; ib < bonds.size(); ++ib) {
			int i = bonds[ib].i;
			int j = bonds[ib].j;
			if (!moved[i] && !moved[j])
				continue;
			std::vector<double> r = rij(crd, i, j);
			double b2 = norm2(r);
			double diff2 = bonds[ib].b02 - b2;
			if (fabs(diff2) <= tolerance * bonds[ib].b02)
				continue;
			std::vector<double> rref = rij(crdref, i, j);
			double sp = dot(r, rref);
			if (sp < 1.e-10)
				return false;
			double invmassi = 1.0 / (masses[i * ndim_]);
			double invmassj = 1.0 / (masses[j * ndim_]);
			double lamda = diff2 / (2.0 * sp * (invmassi + invmassj));
			for (int d = 0; d < ndim_; ++d) {
				crd[i * ndim_ + d] -= lamda * rref[d] * invmassi;
				crd[j * ndim_ + d] += lamda * rref[d] * invmassj;

			}
			moved_next[i] = true;
			moved_next[j] = true;
			done = false;
		}
		moved = moved_next;
		moved_next.assign(moved_next.size(), false);
	}
	return true;
}

ShakeBonds NSPsd::make_shakebonds(const ForceField &ff){
	ShakeBonds sbds;
	auto & bondterms=ff.bondterms();
	for(auto iter=bondterms.begin(); iter !=bondterms.end();++iter){
		sbds.addbond(iter->i,iter->j,iter->b0);
	}
	return sbds;
}

