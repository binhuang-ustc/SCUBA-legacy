/*
 * analytsasa.h
 *
 *  Created on: 2018年3月29日
 *      Author: hyliu
 */

#ifndef SD_ANLYTSASA_H_
#define SD_ANLYTSASA_H_
#include "sd/forcefield.h"
#include "sd/bsinchain.h"
#include "sd/sidechainff.h"
namespace NSPsd{
class Anlytsasa{
public:
	struct Result {
		std::vector<double> atomic_sasa;
		std::vector<std::vector<DvDxi>> dsadx;
	};
	void loopneighbors(const std::vector<NSPgeometry::XYZ> &xyz,
		const std::vector<std::vector<int>> &neighbors,
		double p,
		std::vector<std::vector<double>> &facti,
		std::vector<std::vector<NSPgeometry::XYZ>> &dfactidxi,
		std::vector<std::vector<DvDxi>> & dfactidxj);
	Result calc_anlytsasa(const std::vector<NSPgeometry::XYZ> &xyz,
		const NeighborList &nbl);
	Anlytsasa(const ForceField& ff);
private:
	double p12,p13,p14,pij,r_solv;
	std::map<std::string,double> r_atom;
	std::map<std::string,double> p0_atom;
	std::vector<std::vector<int>> nbor12;
	std::vector<std::vector<int>> nbor13;
	std::vector<std::vector<int>> nbor14;
	std::vector<double> ri;
	std::vector<double> pi;
	std::vector<double> si;
};
}



#endif /* SD_ANLYTSASA_H_ */
