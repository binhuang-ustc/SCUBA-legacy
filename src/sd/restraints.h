/*
 * restraints.h
 *
 *  Created on: 2017年11月10日
 *      Author: hyliu
 */

#ifndef SD_RESTRAINTS_H_
#define SD_RESTRAINTS_H_
#include "geometry/xyz.h"
#include <cassert>
#include <vector>
namespace NSPsd {
class StructRestraint{
public:
	enum {POSIMODE,TOTALMODE};
	StructRestraint(){;}
	StructRestraint(const std::vector<double> &crd, const std::vector<double> &weights,double kres,
			int mode=POSIMODE,double rmsdref=0.0);
	StructRestraint(const std::vector<double> &crd, const std::vector<std::pair<int,int>> &poslen,
			double kres,int mode=POSIMODE,double rmsdref=0.0);
	double posi_energy(const std::vector<NSPgeometry::XYZ> &crd,std::vector<NSPgeometry::XYZ> *forces) const;
	double rmsd2() const{return rmsd2_;}
	double totalrmsd_energy(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const;
	double energy(const std::vector<NSPgeometry::XYZ> &crd,std::vector<NSPgeometry::XYZ> *forces) const{
		if(mode_==POSIMODE) return posi_energy(crd,forces);
		else if(mode_==TOTALMODE) return totalrmsd_energy(crd,forces);
	}
private:
    std::vector<NSPgeometry::XYZ> crdref_;
	std::vector<double> weights_;
	double wtot_{0.0};
	double kres_{0.0};
	double rmsdref_{0.0};
	mutable double rmsd2_{0.0};
	int mode_{POSIMODE};
};
class RgRestraint{
public:
	RgRestraint(double b=0,double k=0):rgbound_(b),kres_(k){;}
	double energy(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces
			,const std::vector<double> & w=std::vector<double>()) const;
	double energy1(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces
			,const std::vector<double> & w=std::vector<double>()) const;
	double energy2(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces
			,const std::vector<double> & w=std::vector<double>()) const;
private:
	double rgbound_;
	double kres_;
};

struct DisRestraint{
	int a1,a2;
	double kres;
	double r0;
	double r1;
	DisRestraint(int i=0,int j=0, double rij0=0,double rij1=0, double k=0):
		a1(i),a2(j),r0(rij0),r1(rij1),kres(k){;}
	double energy(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const{
		if(r0>0) return energy_attract(crd,forces);
		else return energy_repulsion(crd,forces);
	}
	double energy_attract(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const;
	double energy_repulsion(const std::vector<NSPgeometry::XYZ> &crd, std::vector<NSPgeometry::XYZ> *forces) const;
};
class BSInChain;
std::vector<std::pair<int,int>> capairsinstrandpair(int length, int parallel,int s1s,
		int s2s, const std::vector<BSInChain> & bsc1, const std::vector<BSInChain> &bsc2);
StructRestraint makermsdrestraint(const std::string &filename,double rmsdref,double kres);
}



#endif /* SD_RESTRAINTS_H_ */
