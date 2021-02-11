/*
 * forcefield.h
 *
 *  Created on: 2017年11月2日
 *      Author: hyliu
 */

#ifndef SD_FORCEFIELD_H_
#define SD_FORCEFIELD_H_

#include "dataio/parameters.h"
#include "dataio/controlfile.h"
#include "geometry/calculators.h"
#include "backbone/backbonesite.h"
#include "pdbstatistics/phipsidistr.h"
#include "sd/restraints.h"
#include "sd/ssrestraint.h"
#include "sd/bsinchain.h"
#include "sd/sidechainff.h"
#include "sd/activeselections.h"
#include <vector>
#include <functional>
#include <set>
#include <cassert>
#define A2NM 0.1
#define MAXEXCL  15
#define KBT 2.4942
#define KANG_FAC 3282.72
//#define CBPAIRENE
namespace NSPsd{

class GenChain;
//Derivative of a function of a set of atom coordinates with respect the
// coordinates of a specific atom, which is useful when the function
// depends only on a small subset of atom coordinates

struct  BondTerm
{
	double kb;
	double b0;
	int i;
	int j;
	bool off(const std::vector<bool> *forceoff) const{
		if(forceoff) return (forceoff->at(i) && forceoff->at(j));
		else return false;
	}
	BondTerm(int ii,int jj,double bb, double kk): i(ii),j(jj),b0(bb*A2NM),kb(kk/(A2NM*A2NM)){;}
	double energy(const std::vector<NSPgeometry::XYZ> &crd,
			std::vector<NSPgeometry::XYZ> *forces) const;
};

struct AngleTerm{
	double kt;
	double cost0;
	int i;
	int j;
	int k;
	bool off(const std::vector<bool> *forceoff) const{
		if(forceoff) return (forceoff->at(i) && forceoff->at(j) &&forceoff->at(k));
		else return false;
	}
	AngleTerm(int ii, int jj, int kk, double theta0, double kkt):
		i(ii),j(jj),k(kk),kt(kkt){cost0=cos(theta0*3.14159265/180.0);}
	double energy(const std::vector<NSPgeometry::XYZ> &crd,
			std::vector<NSPgeometry::XYZ> *forces) const;
};

struct ImpDihTerm {
	double kp;
	double p0;
	int i;
	int j;
	int k;
	int l;
	ImpDihTerm(int ii,int jj,int kk, int ll, double pp0, double kkp):
		i(ii),j(jj),k(kk),l(ll),kp(kkp) {p0=pp0*3.14159265/180.0;}
	double energy(const std::vector<NSPgeometry::XYZ> &crd,
			std::vector<NSPgeometry::XYZ> *forces) const;
	bool off(const std::vector<bool> *forceoff) const {
		if(forceoff) return (forceoff->at(i) && forceoff->at(j) &&forceoff->at(k)
				&& forceoff->at(l));
		else return false;
	}
};

struct DihTerm {
	double kp;
	double p0;
	int m;
	int i;
	int j;
	int k;
	int l;
	DihTerm(int ii,int jj,int kk, int ll, double pp0, double kkp,int mm):
		i(ii),j(jj),k(kk),l(ll),kp(kkp),m(mm) {p0=pp0*3.14159265/180.0;}
	double energy(const std::vector<NSPgeometry::XYZ> &crd,
			std::vector<NSPgeometry::XYZ> *forces) const;
};
struct PhiPsi {
	enum MODE {USUAL,NTERM,CTERM};
	const NSPpdbstatistics::PhiPsiDistr *distr;
//	std::vector<int> atoms;
	int mode;
	int siteid{-1};
	int chainid{0};
	PhiPsi(int sid,int cid,bool isnterm=false,bool iscterm=false):siteid(sid),chainid(cid){
		if(isnterm) mode=NTERM;
		else if(iscterm) mode=CTERM;
		else mode=USUAL;
//		distr=&(NSPpdbstatistics::PhiPsiDistr::mixcoildistr());
		distr=&(NSPpdbstatistics::PhiPsiDistr::coildistr());
	}
	PhiPsi(int sid,int cid,std::string resname, std::string resnext,bool isnterm=false,bool iscterm=false):siteid(sid),chainid(cid){
		if(isnterm) mode=NTERM;
		else if(iscterm) mode=CTERM;
		else mode=USUAL;
		distr=&(NSPpdbstatistics::PhiPsiDistr::phipsidistr(resname,resnext));
	}
/*	PhiPsi(int cpre,int n, int ca, int c,int nnext,const std::string & resname,const std::string &nextresname){
		atoms.push_back(cpre); atoms.push_back(n);atoms.push_back(ca);atoms.push_back(c);atoms.push_back(nnext);
		if(cpre<0) mode=NTERM;
		else if (nnext < 0 ) mode=CTERM;
		else mode=USUAL;
		distr=&(NSPpdbstatistics::PhiPsiDistr::phipsidistr(resname,nextresname));
	}
	PhiPsi(int cpre,int n, int ca, int c,int nnext){
		atoms.push_back(cpre); atoms.push_back(n);atoms.push_back(ca);atoms.push_back(c);atoms.push_back(nnext);
		if(cpre<0) mode=NTERM;
		else if (nnext < 0 ) mode=CTERM;
		else mode=USUAL;
		distr=&(NSPpdbstatistics::PhiPsiDistr::mixcoildistr());
	}*/

//	double energy(const std::vector<NSPgeometry::XYZ> &crd, double phipsiwght,
//			std::vector<NSPgeometry::XYZ> *forces) const;
	double energy(const std::vector<PhiPsiCodes> &phipsicodes,double phipsiwght,
			std::vector<NSPgeometry::XYZ> *forces) const;
};

struct StericAtom{
	double sigma;
	double eps;
	double sigmahb;
	static double DEDRMAX;
	bool hbdonor{false};
	bool hbacceptor{false};
	StericAtom(double ss,double ee):sigma(ss*A2NM),eps(ee),hbdonor(false),hbacceptor(false),sigmahb(ss*A2NM){;}
	StericAtom(double ss, double ee,bool donor,bool acceptor,double sshb):
		sigma(ss*A2NM),eps(ee),hbdonor(donor),hbacceptor(acceptor),sigmahb(sshb*A2NM){;}
	static double pairenergy(const NSPgeometry::XYZ & posi1, const NSPgeometry::XYZ &pos2,
				const StericAtom & sa1, const StericAtom &sa2,
				std::vector<NSPgeometry::XYZ> *forces);
};
class NeighborList;
struct EneComp{
	std::vector<double> ephipsi;
	std::vector<double> els;
	std::vector<double> elocalhb;
	std::vector<double> epacking;
	std::vector<double> escconf;
	std::vector<std::vector<double>> esitepair;
	int siteoffset{0}; //a working variable
	EneComp(){;}
	void init(int nsites){
		ephipsi.assign(nsites,0.0);
		els.assign(nsites,0.0);
		epacking.assign(nsites,0.0);
		elocalhb.assign(nsites,0.0);
		escconf.assign(nsites,0.0);
		esitepair.assign(nsites,std::vector<double>());
		for(int i=0;i<nsites;++i){
			esitepair[i].assign(nsites,0.0);
		}
	}
	void print(std::ostream &os);
};
class ForceField {
public:
	enum ENECOMP {ETOT,EBOND,EANG,EIMPDIH,EPHIPSI,ESCCONF,ESTERIC,ESCPACKING,ELOCALSTRUCTURE, ELOCALHB,
		ESITEPAIRS,ESHADOW,ESTRUCTREST,ERGRESTRAINT,ESSRESTRAINT,ETERMS};
	void init(int natoms) {
		natoms_=natoms;
		bondterms_.clear();
		angleterms_.clear();
		impdihterms_.clear();
		dihterms_.clear();
		stericatoms_.clear();
		stericatoms_.resize(natoms,StericAtom(0.0,0.0));
		excluded_.clear();
		excluded_.resize(natoms);
		list14_.clear();
		list14_.resize(natoms);
	}
	void addbond(int i, int j,double b0,double kb){
		assert(i<natoms_ && j<natoms_);
		bondterms_.push_back(BondTerm(i,j,b0,kb));
	}
	void addangle(int i,int j, int k, double theta0, double kt){
		assert(i<natoms_ &&j<natoms_ && k<natoms_);
		angleterms_.push_back(AngleTerm(i,j,k,theta0,kt));
	}
	void addimpdih(int i, int j, int k, int l, double p0, double kp){
		assert(i<natoms_ &&j<natoms_ &&k<natoms_ &&l<natoms_);
		impdihterms_.push_back(ImpDihTerm(i,j,k,l,p0,kp));
	}
	void adddih(int i, int j, int k, int l, double p0=0.0, double kp=0.0,int m=0){
		assert(i<natoms_ &&j<natoms_ &&k<natoms_ &&l<natoms_);
		dihterms_.push_back(DihTerm(i,j,k,l,p0,kp,m));
	}
	void addphipsi(int sid,int cid,const std::string &resname, const std::string &resnext,
			bool nterm,bool cterm){
		phipsis_.push_back(PhiPsi(sid,cid,resname,resnext,nterm,cterm));
	}
	void addphipsi(int sid,int cid, bool nterm,bool cterm){
		phipsis_.push_back(PhiPsi(sid,cid,nterm,cterm));
	}
	void addstructrestraint(const StructRestraint &restraint){
		structrestraints_.push_back(restraint);
	}
	void addssrestraint(int chainid,int sitebegin,int siteend,
			int ssid, double targetp,double forceconstant){
		ssrestraints_.push_back(SSRestraint(chainid,sitebegin,siteend,ssid,targetp,forceconstant));
	}
	void addbsinchain(int chainidx,int nid,int caid,int cid,int oid){
		if(chainidx+1>bsinchains_.size()) bsinchains_.resize(chainidx+1);
		bsinchains_[chainidx].push_back(BSInChain(nid,caid,cid,oid));
	}
	void setstericatom(int i,double sigma,double eps){
		assert(i<natoms_);
		stericatoms_.at(i)=StericAtom(sigma,eps,false,false,sigma);
	}
	void setstericatom(int i,double sigma,double eps,bool hbdonor,bool hbacceptor,double sigmahb){
		assert(i<natoms_);
		stericatoms_.at(i)=StericAtom(sigma,eps,hbdonor,hbacceptor,sigmahb);
	}
	void setexcln14();
	void setlstermon(bool val=true) {lstermon_=val;}
	void setsitepairtermon(bool val=true){sitepairtermon_=val;}
	void setanalysismode(bool val=true){eneanalysismode_=val;}
	void setrgrestraint(double rgbound,double k){
		rgrestraint_=RgRestraint(rgbound,k);
	}
	std::vector<double> forces(const std::vector<double> &crd,
			std::vector<double> *potenergies,const std::vector<bool> *forceoff=nullptr) const;
//	std::vector<double> forces(const std::vector<double> &crd, const NeighborList &nbl,
//			std::vector<double> *potenergies,const std::vector<bool> *forceoff=nullptr) const;
	std::vector<double> forces(const std::vector<double> &crd, const NeighborList &nbl,
				std::vector<double> *potenergies,ActiveSelections &acts) const;
	const std::vector<BondTerm> & bondterms() const {return bondterms_;}
	const std::vector<AngleTerm> & angleterms() const {return angleterms_;}
	const std::vector<std::set<int>> &list14() const {return list14_;}
	bool excluded(int i,int j) const {
		return excluded_[i].find(j) != excluded_[i].end();
	}
	bool is14(int i,int j) const{
	 	return list14_[i].find(j) != list14_[i].end();
	}
//	double lsenergy(const std::vector<PhiPsiCodes> & phipsicodes, const std::vector<double> & siteweight,
//			std::vector<NSPgeometry::XYZ> *xyzf) const;
	double lsenergy(const std::vector<PhiPsiCodes> & phipsicodes, const std::vector<double> & siteweight,
			std::vector<NSPgeometry::XYZ> *xyzf,
			const std::vector<bool> &lsactive) const;
#ifdef _OPENMP
//	std::vector<double> lsenergy(const std::vector<PhiPsiCodes> & phipsicodes, const std::vector<double> & siteweight,
//			std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread) const;
	std::vector<double> lsenergy(const std::vector<PhiPsiCodes> & phipsicodes, const std::vector<double> & siteweight,
			std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread,
			const std::vector<bool> &lsactive) const;
#endif
	double shadowenergy(const std::vector<NSPgeometry::XYZ> &xyz,std::vector<NSPgeometry::XYZ> *xyzf) const;
	template <typename CVLTERM>
	double covalentenergy(const std::vector<CVLTERM> &cvlterms,
			const std::vector<NSPgeometry::XYZ> &xyz, std::vector<NSPgeometry::XYZ> *xyzf,
			 const std::vector<bool> *forceoff=nullptr) const {
		double ene=0.0;
		for(auto &t:cvlterms){
			if(t.off(forceoff))continue;
			ene += t.energy(xyz,xyzf);
		}
		return ene;
	}
#ifdef _OPENMP
	template <typename CVLTERM>
	std::vector<double> covalentenergy(const std::vector<CVLTERM> &cvlterms,
			const std::vector<NSPgeometry::XYZ> &xyz,
			std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread,
			 const std::vector<bool> *forceoff=nullptr) const {
		int nthread=omp_get_max_threads();
		std::vector<double> ene_thread(nthread,0.0);
		int nterms=cvlterms.size();
#pragma omp parallel for schedule(dynamic,1)
		for(int i=0;i<nterms; ++i){
			int myid=omp_get_thread_num();
			double & ene=ene_thread[myid];
			std::vector<NSPgeometry::XYZ> *xyzf=&(xyzf_thread[myid]);
			const CVLTERM & t=cvlterms.at(i);
			if(t.off(forceoff))continue;
			ene += t.energy(xyz,xyzf);
		}
		return ene_thread;
	}
#endif
	 double stericenergy (const NeighborList &nbl,const std::vector<NSPgeometry::XYZ> &xyz,
			 std::vector<NSPgeometry::XYZ> *xyzf) const;
#ifdef _OPENMP
	 std::vector<double> stericenergy (const NeighborList &nbl,const std::vector<NSPgeometry::XYZ> &xyz,
			 std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread) const;
#endif
/*	double sitepairenergy(const std::vector<double> & crd,
			const std::vector<std::vector<PhiPsiCodes>> & phipsicodes,
			const std::vector<std::vector<SSCode>> &sscodes,std::vector<NSPgeometry::XYZ> *xyzf,
			double *etrap) const;*/
	double sitepairenergy(const std::vector<double> & crd,
			const std::vector<std::vector<PhiPsiCodes>> & phipsicodes,
			const std::vector<std::vector<SSCode>> &sscodes,std::vector<NSPgeometry::XYZ> *xyzf,
			double *etrap,const ActiveSelections &acts) const;
#ifdef _OPENMP
/*	std::vector<double> sitepairenergy(const std::vector<double> & crd,
			const std::vector<std::vector<PhiPsiCodes>> & phipsicodes,
			const std::vector<std::vector<SSCode>> &sscodes,
			std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread,
			std::vector<double> & etrap_thread) const;*/
	std::vector<double> sitepairenergy(const std::vector<double> & crd,
			const std::vector<std::vector<PhiPsiCodes>> & phipsicodes,
			const std::vector<std::vector<SSCode>> &sscodes,
			std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread,
			std::vector<double> & etrap_thread,
			const ActiveSelections &acts) const;
#endif
/*	double sideconfene(const std::vector<double> &crd,
		    std::vector<std::vector<PhiPsiCodes>> &phipsicodes,
			std::vector<NSPgeometry::XYZ> *xyzf) const;*/
	double sideconfene(const std::vector<double> &crd,
			std::vector<std::vector<PhiPsiCodes>> &phipsicodes,
			std::vector<NSPgeometry::XYZ> *xyzf,const ActiveSelections &acts) const;
#ifdef _OPENMP
/*	std::vector<double> sideconfene(const std::vector<double> &crd,
		    std::vector<std::vector<PhiPsiCodes>> &phipsicodes,
			std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread) const;*/
	std::vector<double> sideconfene(const std::vector<double> &crd,
			    const std::vector<std::vector<PhiPsiCodes>> &phipsicodes,
				std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread,
				const ActiveSelections & acts) const;
#endif
/*	double scpackingene(const std::vector<NSPgeometry::XYZ> &xyz,
			std::vector<NSPgeometry::XYZ> *xyzf) const;*/
	double scpackingene(const std::vector<NSPgeometry::XYZ> &xyz,
				std::vector<NSPgeometry::XYZ> *xyzf,const ActiveSelections &acts) const;
#ifdef _OPENMP
/*	std::vector<double> scpackingene(const std::vector<NSPgeometry::XYZ> &xyz,
			std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread) const;*/
	std::vector<double> scpackingene(const std::vector<NSPgeometry::XYZ> &xyz,
			std::vector<std::vector<NSPgeometry::XYZ>> &xyzf_thread,
			const ActiveSelections &acts) const;
#endif
	double localhbenergy(const std::vector<double> &crd,std::vector<NSPgeometry::XYZ> *xyzf) const;
	double localbbhbenergy(const std::vector<NSPgeometry::XYZ> &xyz, std::vector<NSPgeometry::XYZ> *xyzf) const;
	bool sitepairoff(const std::vector<SSCode> & sscodes,int posi1,int posi2) const;
	void usecontrols(std::string controlname){
		usecontrols(controlname,nullptr);
	}
	void usecontrols(std::string controlname, const GenChain *genchain);
	int siteid(int chainid,int posi) const {
		if(chainid==0) return posi;
		int offset=0;
		for(int i=0;i<chainid;++i){
			offset+=bsinchains_[i].size();
		}
		return offset+posi;
	}
	EneComp &enecomp() const {return enecomp_;}
	std::vector<int> & stericatomtypes(){return stericatomtypes_;}
	std::vector<std::vector<SCInChain>> &scinchains(){return scinchains_;}
	const std::vector<std::vector<SCInChain>> &scinchains() const {return scinchains_;}
	std::vector<std::vector<BSInChain>> &bsinchains(){return bsinchains_;}
	const std::vector<std::vector<BSInChain>> &bsinchains() const {return bsinchains_;}
	std::vector<int> mainchainatoms() const {
		std::vector<int> mca;
		for(auto &c:bsinchains_){
			for(auto &p:c) {
				std::vector<int> atms=p.atomids();
				for(int a:atms) mca.push_back(a);
			}
		}
		return mca;
	}
	std::vector<int> sidechainatoms() const {
		std::vector<int> sca;
		for(auto &c:scinchains_){
			for(auto &p:c){
				for(int i=0;i<p.nscatoms;++i) sca.push_back(p.poffset+i);
			}
		}
		return sca;
	}
	std::vector<int> atomsinresidue(int chainid,int posi) const {
		std::vector<int> ia=bsinchains_.at(chainid).at(posi).atomids();
		for(int i=0;i<scinchains_.at(chainid).at(posi).nscatoms;++i)
			ia.push_back(scinchains_.at(chainid).at(posi).poffset+i);
		return ia;
	}
	const bool & scmmsteric() const {return scmmsteric_;}
	bool & scmmsteric() {return scmmsteric_;}
	int natoms() const {return natoms_;}
        static double highescale;
private:
	std::vector<BondTerm> bondterms_;
	std::vector<AngleTerm> angleterms_;
	std::vector<ImpDihTerm> impdihterms_;
	std::vector<DihTerm> dihterms_;
	std::vector<StericAtom> stericatoms_;
	std::vector<int> stericatomtypes_;
	std::vector<std::set<int>> excluded_;
	std::vector<std::set<int>> list14_;
	std::vector<std::vector<BSInChain>> bsinchains_;
	std::vector<std::vector<SCInChain>> scinchains_;
	std::vector<PhiPsi> phipsis_;
	std::vector<StructRestraint> structrestraints_;
	std::vector<SSRestraint> ssrestraints_;
	std::vector<DisRestraint> disrestraints_;
	std::vector<double> coreweights_;
	std::vector<std::vector<double>> siteweights_;
//	std::vector<std::vector<bool>> shadowon_;
	RgRestraint rgrestraint_;
	bool lstermon_{false};
	bool sitepairtermon_{false};
//	mutable const std::vector<bool> *forceoff_{nullptr};
	double phipsiwght_{KBT};
	double scconfwght_{KBT};
	double stericwght_{1.0};
	double lsweight_{0.0};
	double shadowweight_{0.0};
	double scpackingweight_{0.0};
	double localhbweight_{0.0};
	double sitepairweight_{0.0};
	double eattract_{0.0};
	double rattractoff_{10000.0};
	double rattractswitch_{1000.0};
	int natoms_{0};
	double etrap(double r,double *dfdr) const;
	bool eneanalysismode_{false};
	bool scmmsteric_{false};
	mutable EneComp enecomp_;
};

struct NeighborList {
	std::vector<std::vector<int>> neighbors;
	NeighborList(const std::vector<double> &crd, const ForceField &ff,
			const std::vector<bool> *forceoff=nullptr);
};

typedef NSPdataio::TypedMultiInstanceControls<ForceField> ForceFieldControls;
inline NSPdataio::ParameterSet &getforcefieldcontrol(const std::string &name=std::string()){
	return ForceFieldControls::getparameterset(name);
}
void defineforcefieldcontrol(std::string name, const std::vector<std::string> &controllines);
std::vector<std::string> ssrestraincontrolfromchain(const std::vector<
		NSPproteinrep::BackBoneSite> &chain,int chainid,double p0,double kres);
std::string rgrestraincontrolfromchain(const std::vector<NSPproteinrep::BackBoneSite> &chain,
		double krgres);
//ForceField make_forcefield_backbone(int nsites);
//ForceField make_forcefield_backbone(const std::vector<NSPproteinrep::BackBoneSite> & chain);
ForceField make_forcefield_backbone(const std::vector<int> & nsites,
		std::vector<std::vector<int>> cissites=std::vector<std::vector<int>>());
ForceField make_forcefield_backbone(
		const std::vector<std::vector<NSPproteinrep::BackBoneSite>> & chain);
}




#endif /* SD_FORCEFIELD_H_ */
