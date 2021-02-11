/*
 * shadowterm.h
 *
 *  Created on: 2018年2月4日
 *      Author: hyliu
 */

#ifndef SD_SHADOWTERM_H_
#define SD_SHADOWTERM_H_

#include "sd/bsinchain.h"
#include "sd/sidechainff.h"
#include "sd/enefunc1d.h"
#include "pdbstatistics/virtualsidechain.h"
namespace NSPsd{
class ShadowTerm{
public:
	static std::vector<EneFunc1D *> enefuncs;
	static NullEne1D nullene1d;
	static std::map<std::string, ShiftedG> shiftgfuncs;
	static std::map<std::string,SineEne> sinefuncs;
	static void setupenefuncs();
//	static double atompairenergy(const std::string &a1, const std::string &a2, double r, double *dedr);
	static double sitepairenergy(const std::vector<NSPgeometry::XYZ> & bcrd1,
			const NSPpdbstatistics::VirtualSideChain & vsc1,
			const std::vector<NSPgeometry::XYZ> &bcrd2,const NSPpdbstatistics::VirtualSideChain &vsc2,
			std::vector<NSPgeometry::XYZ> *forces1, std::vector<NSPgeometry::XYZ> *forces2);
	static double sitepairenergy(const std::vector<NSPgeometry::XYZ> & bcrd1,
			const NSPpdbstatistics::VirtualSideChain & vsc1,
			const std::vector<NSPgeometry::XYZ> & crds, const BSInChain & bs2, const SCInChain & sc2,
			std::vector<NSPgeometry::XYZ> *forces1, std::vector<DvDxi> *dedx2);
};
}



#endif /* SD_SHADOWTERM_H_ */
