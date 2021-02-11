/*
 * sidechainff.h
 *
 *  Created on: 2018年1月31日
 *      Author: hyliu
 */

#ifndef SD_SIDECHAINFF_H_
#define SD_SIDECHAINFF_H_
#include "sd/bsinchain.h"
#include "sd/enefunc1d.h"
#include <vector>
#include <string>
#include <map>
#include <set>

namespace NSPsd{
class PackingEneFuncs{
public:
	EneFunc1D & getenefunc(int atype,int btype) const {
		return *(functable_[eidx(atype,btype)]);
	}
	void setup(int ntypes){
		ntypes_=ntypes;
		functable_.assign(ntypes*ntypes,nullptr);
	}
	void addfunction(int a,int b,std::shared_ptr<EneFunc1D> func){
		assert(a<ntypes_ && b<ntypes_);
		functable_[eidx(a,b)]=func;
	}
private:
	std::vector<std::shared_ptr<EneFunc1D>> functable_;
	int ntypes_{0};
	int eidx(int a, int b) const {return a*ntypes_+b;}
};
struct PackingAtomType{
	double radius;
	int hbtype;
	int aromatic;
	int nconnections;
	static bool hbond(int hb1,int hb2){
		if(hb1 == 1 ){
			if(hb2==2 || hb2==3) return true;
		} else if(hb2==1){
			if(hb1==2 || hb1==3) return true;
		}
		return false;
	}
};
struct VSCType{
	static const VSCType & getVSCType(const std::string & resname, const std::string &filename=std::string());
	static std::map<std::string,VSCType> readVSCTypes(const std::string & filename);
	static std::vector<PackingAtomType>packingatomtypes;
	static PackingEneFuncs packingenefuncs;
	static PackingEneFuncs softpackingenefuncs;
	static std::map<std::string,int> stericatomtypes;
	static std::map<char,std::string> resnamefrom1letter;
	static std::map<std::string,std::set<std::string>> rotatablescatoms;
	static bool isrotatablescatom(const std::string & resname, const std::string & atomname) {
	    return rotatablescatoms.find(resname) != rotatablescatoms.end()
	            && rotatablescatoms.at(resname).find(atomname) != rotatablescatoms.at(resname).end();
	}
	static std::string resnameof( char letter){
		if(resnamefrom1letter.empty()) getVSCType("ALA");
		return resnamefrom1letter.at(letter);
	}
	static int getstericatomtype(const std::string &resname,const std::string &atomname){
		std::string fullname=resname+":"+atomname;
		auto it=stericatomtypes.find(fullname);
		if(it== stericatomtypes.end()){
			fullname=std::string("ANY:")+atomname;
			it=stericatomtypes.find(fullname);
		}
		if(it==stericatomtypes.end()) return -1;
		return  it->second;
	}
	int nscatoms;
	std::string resname;
	std::string pdbname;
	char oneletter{'\0'};
	std::vector<std::string> atomnames;
	std::vector<std::vector<std::pair<int,double>>> internalcrds;
	std::vector<int> rotameratoms;
	std::vector<std::pair<int,int>> newbonds;
	std::vector<double> b0;
	std::vector<double> kb0;
	std::vector<std::vector<int>> newangles;
	std::vector<double> a0;
	std::vector<double> ka0;
	std::vector<std::vector<int>> newimpdihs;
	std::vector<double> imp0;
	std::vector<double> kimp0;
	std::vector<std::vector<int>> newtorsions;
};
struct SCInChain {
	std::string restype;
	int poffset{-1}; //atom ID of the CB atom;
	int nscatoms{0};
	bool softpacking{false};
	std::vector<std::vector<int>> kaiatoms;
	SCInChain(){;}
	SCInChain(const std::string & name,
			const std::vector<std::vector<int>> &kais,int caposi,int nasc):restype(name),
					kaiatoms(kais),poffset(caposi),nscatoms(nasc){;}
};
struct ConformerCode {
	PhiPsiCodes *phipsicodes;
	std::vector<double> sidechaintorsions;
	std::vector<std::vector<double>> torsioncodes;
	std::vector<std::vector<std::vector<DvDxi>>> dtorsioncodesdx;
	std::string restype;
	static std::vector<double> gettorsioncodes(double ang,
			const std::vector<DvDxi> &dadx, std::vector<std::vector<DvDxi>> *dcdx);
};
std::vector<ConformerCode> makeconformercodes(const std::vector<double> &crds,
		const std::vector<SCInChain> &scinchains,
		std::vector<PhiPsiCodes> &phipsicodes);
double mcscpackingenergy(const std::vector<NSPgeometry::XYZ> &crds, const std::vector<int> &atomtypes,
		const BSInChain &mc,const SCInChain &sc,
		int sep,std::vector<DvDxi> *dedx,bool terminal=false);
double scscpackingenergy(const std::vector<NSPgeometry::XYZ> &crds,
		const std::vector<int> &atomtypes,
		const SCInChain &sc1, const SCInChain &sc2,
		int sep,std::vector<DvDxi> *dedx);
}



#endif /* SD_SIDECHAINFF_H_ */
