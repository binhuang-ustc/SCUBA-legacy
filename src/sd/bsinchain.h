/*
 * bsliteinchain.h
 *
 *  Created on: 2017年12月12日
 *      Author: hyliu
 */

#ifndef SD_BSINCHAIN_H_
#define SD_BSINCHAIN_H_
#include "geometry/xyz.h"
#include "backbone/backbonesite.h"
#include <vector>
#include <cassert>

namespace NSPsd{
struct BSInChain{
	BSInChain(){;}
	BSInChain(int n,int ca, int c, int o): nid(n),caid(ca),cid(c),oid(o){;}
	int nid{0},caid{1},cid{2},oid{3};
	std::vector<int> atomids() const{
		std::vector<int> res;
		res.push_back(nid);res.push_back(caid);res.push_back(cid);res.push_back(oid);
		return res;
	}
};
typedef std::pair<int, NSPgeometry::XYZ> DvDxi;
inline DvDxi operator*(double c, const DvDxi &dvdxi){
	return std::make_pair(dvdxi.first,c*(dvdxi.second));
}
inline DvDxi operator+(const DvDxi & d1, const DvDxi &d2){
	assert(d1.first==d2.first);
	return std::make_pair(d1.first, d1.second+d2.second);
}
inline NSPgeometry::XYZ getxyz(const std::vector<double> &crd, int posi){
	int idx=3*posi;
	return NSPgeometry::XYZ(crd[idx],crd[idx+1],crd[idx+2]);
}
struct PhiPsiCodes{
    static std::vector<double> gettriangcodes(double ang,const std::vector<DvDxi> &dadx,
    		std::vector<std::vector<DvDxi>>  *dcdx);
    PhiPsiCodes(){;}
    PhiPsiCodes(const std::vector<double> &crd, std::vector<BSInChain>::const_iterator bs,
    		bool Nterm=false,bool Cterm=false);
    double phi{360.0};
    double psi{360.0};
    std::vector<DvDxi> dphidx;
    std::vector<DvDxi> dpsidx;
	std::vector<double> phicodes;
	std::vector<double> psicodes;
	std::vector<std::vector<DvDxi>> dphicodesdx; //dimension: code_length X 4
	std::vector<std::vector<DvDxi>> dpsicodesdx;
};
struct SSCode{
	std::vector<double> p3;
	std::vector<std::vector<DvDxi>> dp3dx;
	int ssid;
};
std::vector<BSInChain> makeBSInChain(int nsites,int idx_offset=0);
std::vector<PhiPsiCodes> makephipsicodes(const std::vector<double> &crds,
		const std::vector<BSInChain> & bsinchain);
std::vector<SSCode> estimatess(const std::vector<PhiPsiCodes> &phipsicodes);
std::vector<std::vector<double>> estimatess(const std::vector<NSPproteinrep::BackBoneSite> & chain);
}




#endif /* SD_BSINCHAIN_H_ */
