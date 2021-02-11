/*
 * aasidechain.h
 *
 *  Created on: 2018年1月17日
 *      Author: hyliu
 */

#ifndef FULLSITE_FULLSITE_H_
#define FULLSITE_FULLSITE_H_
#include "backbone/backbonesite.h"
namespace NSPproteinrep {

struct AtomTopo{
	static const std::vector<AtomTopo> & definedatomtypes(const std::string &filename="");
	static int getatomtype(const std::string &residuename, const std::string &atomname);
	static const AtomTopo & getatomtopo(const std::string &residuename, const std::string &atomname){
		return definedatomtypes().at(getatomtype(residuename,atomname));
	}
	bool hbdonor{false};
	bool hbacceptor{false};
	bool aromatic{false};
	int netcharged{0};
	int elementid{-1};
	int atomtype{-1};
	std::vector<std::string> containing_residues;
	std::vector<std::string> atomnames;
	bool readtopo(const std::vector<std::string> & atomtypeline);
};

class SideChainTopo {
public:
	static const std::map<std::string,SideChainTopo> & definedsidechaintypes(
			const std::string & filename="");
	static bool residuetypedefined(const std::string &resname){
		return(definedsidechaintypes().find(resname) != definedsidechaintypes().end());
	}
	struct SideChainAtomTopo{
		int atomtype{-1};
		const std::vector<AtomTopo> *atomtopos{nullptr}; //for standard side chain atoms this is null
		std::string upatom;
		std::map<std::string,std::vector<double>> downatoms;
	};
	const AtomTopo & getatomtopo(const std::string &atomname) const {
		const SideChainAtomTopo & atom=atoms_.at(atomname);
		if(atom.atomtopos){
			return atom.atomtopos->at(atom.atomtype);
		} else {
			return AtomTopo::definedatomtypes().at(atom.atomtype);
		}
	}
	const SideChainAtomTopo & getsidechainatomtopo(const std::string &atomname) const {
		return atoms_.at(atomname);
	}
	std::string upperatom(const std::string &atomname) const {
		return atoms_.at(atomname).upatom;
	}
	std::map<std::string,std::vector<double>> downatoms(const std::string &atomname) const {
		return atoms_.at(atomname).downatoms;
	}
	const std::map<std::string, SideChainAtomTopo> atomtopos() const{
		return atoms_;
	}
	std::string atomname(int i) const {return atomnames_[i];}
	int natoms() const {return atomnames_.size();}
	const std::vector<std::string> & atomnames() const{
		return atomnames_;}
private:
	std::map<std::string,SideChainAtomTopo> atoms_;
	std::vector<std::string> atomnames_;
	std::string restype;
	bool readnewatom(const std::vector<std::string> &newatomline);
};


class FullSite{
public:
	static bool mainchainatom(const std::string &atomname) {
		return(atomname =="N" || atomname=="CA" ||atomname=="C"||atomname=="O");
	}
	void addatomcrd(const std::string &atomname,NSPgeometry::XYZ);
	void updateAtomCrd(const std::string & atomName, NSPgeometry::XYZ crd);
	void copybackbonecrd(const BackBoneSite &bs);
	bool sidechaincrdcomplete(){
		auto &sct=getsidechaintopo();
		for(int i=0;i<sct.natoms();++i){
			if(crds_.find(sct.atomname(i)) == crds_.end()) return false;
		}
		return true;
	}
	BackBoneSite getbackbonesite() const;
	void buildsidechain(const std::string &beginatom="CB",
			const std::map<std::string,double> &newtorsions=std::map<std::string,double>());
	NSPgeometry::XYZ getcrd(const std::string &atomname) const {
		return crds_.at(atomname);
	}
	double gettorsionat(const std::string &atomname) const ;
//	std::map<std::string,std::vector<double>> internalcrds()const;
	const SideChainTopo & getsidechaintopo() const {
		if(sidechaintopos_) return sidechaintopos_->at(restype_);
		else return SideChainTopo::definedsidechaintypes().at(restype_);
	}
	const SideChainTopo::SideChainAtomTopo & getsidechainatomtopo(const std::string &atomname)
	const {
		return getsidechaintopo().getsidechainatomtopo(atomname);
	}
	const std::string & resname() const {return restype_;}
	std::string &resname() {return restype_;}
	char & chainid() {return chainid_;}
	const char &chainid() const {return chainid_;}
	const std::string & pdbid() const {return pdbid_;}
	std::string &pdbid() {return pdbid_;}
	const int & resid() const {return resid_;}
	int &resid() {return resid_;}
	int &resseq() {return resseq_;}
	const int &resseq() const {return resseq_;}
	int natoms_sidechaintopo() const {
		return getsidechaintopo().natoms();
	}
	const std::vector<std::string> &atomnames_sidechaintopo() const {
		return getsidechaintopo().atomnames();
	}
	std::vector<std::string> atomnames_all() const;
	std::vector<std::string> atomnames_withcrd() const;
	std::string alltopdb(int atomoffset=1) const;
	std::string sidechaintopdb(int atomoffset=1) const;
	PdbRecord makepdbrecord(const std::string &atomname,int atomid) const;
	void readsite_cartesian(std::ifstream &ifs);
	void writesite_cartesian(std::ofstream &ofs) const;
	void writesite_internal(std::ofstream &ofs) const;

	std::vector<NSPgeometry::XYZ> extractcrd()const;
	bool hasatomcrd(const std::string & atomname) const{
		return crds_.find(atomname) != crds_.end();
	}
	std::vector<double> sidechaintorsions() const;
	std::vector<double> getatomsasa(int nsurfpoints=256)const;
	std::vector<std::vector<std::pair<int,double>>> internalcrds() const;
	int getcrds(std::vector<NSPgeometry::XYZ> &crds) const;
	int getradii(std::vector<double> &radii) const;
	int natomswithcrds() const {return crds_.size();}
private:
	const std::map<std::string,SideChainTopo> *sidechaintopos_{nullptr};
	std::map<std::string,NSPgeometry::XYZ>crds_;
	std::string restype_;
	char chainid_{'A'};
	std::string pdbid_{"NA"};
	int resid_{0};
	int resseq_{0};
};
std::vector<std::vector<FullSite>>  readfullsitesfrompdb(const std::string &pdbfilename,
		bool keepnonstandard=false);

std::vector<BackBoneSite> extractbackbonesegment(std::vector<FullSite>::const_iterator begin, int length);
inline std::vector<BackBoneSite> backbone(const std::vector<FullSite> & fullsites){
	return extractbackbonesegment(fullsites.begin(),fullsites.size());
}
std::vector<std::map<std::string,double>> calc_sasa(
		const std::vector<std::vector<FullSite>> &chains,int nsurfpoints=256);

FullSite make_fullsite(const BackBoneSite &bs,
		std::vector<std::vector<std::pair<int,double>>> internalcrds
		=std::vector<std::vector<std::pair<int,double>>>());
void writetopdb(const std::vector<std::vector<FullSite>> &sites,std::ostream &os);

}
#endif /* FULLSITE_FULLSITE_H_ */
