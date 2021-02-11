/*
 * fullsite.cpp
 *
 *  Created on: 2018年1月17日
 *      Author: hyliu
 */

#include "fullsite/fullsite.h"
#include "sd/sidechainff.h"
#include "dataio/datapaths.h"
#include "dataio/inputlines.h"
#include "proteinrep/pdbreader.h"
#include "geometry/atomsasa.h"
using namespace NSPproteinrep;
static std::string atomname(const NSPsd::VSCType &vsc,int i){
	if(i==0) return "N";
	else if(i==1) return "CA";
	else if (i-vsc.nscatoms==2) return "C";
	else if (i-vsc.nscatoms==3) return "O";
	else return vsc.atomnames[i-2];
}
std::vector<double> FullSite::sidechaintorsions() const {
	const NSPsd::VSCType & vsc=NSPsd::VSCType::getVSCType(restype_);
	std::vector<double> torsions;
	double rad=180.0/3.14159265;
	for(auto i:vsc.rotameratoms){
			std::string la=atomname(vsc,i+2);
			std::string ka=atomname(vsc,vsc.internalcrds[i][0].first);
			std::string ja=atomname(vsc,vsc.internalcrds[i][1].first);
			std::string ia=atomname(vsc,vsc.internalcrds[i][2].first);
			if(hasatomcrd(ia),hasatomcrd(ja) && hasatomcrd(ka) && hasatomcrd(la)){
				torsions.push_back(
						rad*NSPgeometry::torsion(getcrd(ia),getcrd(ja),getcrd(ka),getcrd(la)));
			} else {
				torsions.clear();
				break;
			}
	}
	return torsions;
}
int FullSite::getcrds(std::vector<NSPgeometry::XYZ> &crds) const{
	int nold=crds.size();
	if(hasatomcrd("N")) crds.push_back(getcrd("N"));
	if(hasatomcrd("CA")) crds.push_back(getcrd("CA"));
	const NSPsd::VSCType & vsc=NSPsd::VSCType::getVSCType(restype_);
	for(auto &a:vsc.atomnames){
		if(hasatomcrd(a)) crds.push_back(getcrd(a));
	}
	if(hasatomcrd("C")) crds.push_back(getcrd("C"));
	if(hasatomcrd("O")) crds.push_back(getcrd("O"));
	return crds.size()-nold;
}
int FullSite::getradii(std::vector<double> &radii) const{
	std::vector<std::string> atoms;
	int nold=radii.size();
	atoms.push_back("N");
	atoms.push_back("CA");
	const NSPsd::VSCType & vsc=NSPsd::VSCType::getVSCType(restype_);
	for(auto &a:vsc.atomnames) atoms.push_back(a);
	atoms.push_back("C");
	atoms.push_back("O");
	for(auto &a:atoms){
		if(!hasatomcrd(a)) continue;
		int at=NSPsd::VSCType::getstericatomtype(restype_,a);
		radii.push_back(NSPsd::VSCType::packingatomtypes[at].radius*0.91);
	}
	return radii.size()-nold;
}

std::vector<double> FullSite::getatomsasa(int nsurfpoints) const{
	std::vector<NSPgeometry::XYZ> crd;
	std::vector<double> radii;
	int nexpect=4+NSPsd::VSCType::getVSCType(restype_).nscatoms;
	if(getcrds(crd) != nexpect) return std::vector<double>();
	getradii(radii);
	return NSPgeometry::calc_sasa(crd,radii,1.4,nsurfpoints);
}
std::vector<std::vector<std::pair<int,double>>> FullSite::internalcrds() const{
	const NSPsd::VSCType & vsc=NSPsd::VSCType::getVSCType(restype_);
	std::vector<std::vector<std::pair<int,double>>> ics=vsc.internalcrds;
	double rad=180.0/3.14159265;
	for(int i=0;i<vsc.nscatoms;++i){
			std::string la=atomname(vsc,i+2);
			std::string ka=atomname(vsc,vsc.internalcrds[i][0].first);
			std::string ja=atomname(vsc,vsc.internalcrds[i][1].first);
			std::string ia=atomname(vsc,vsc.internalcrds[i][2].first);
			if(hasatomcrd(ia),hasatomcrd(ja) && hasatomcrd(ka) && hasatomcrd(la)){
				ics[i][0].second=NSPgeometry::distance(getcrd(la),getcrd(ka));
				ics[i][1].second=rad*NSPgeometry::angle(getcrd(la),getcrd(ka),getcrd(ja));
				ics[i][2].second=
						rad*NSPgeometry::torsion(getcrd(ia),getcrd(ja),getcrd(ka),getcrd(la));
				if(vsc.resname=="ARG" && la=="NH1"){
					double torsion=ics[i][2].second;
					if(torsion<-90.0) torsion +=180.0;
					else if(torsion>90.0) torsion-=180.0;
					ics[i][2].second=torsion;
				}
			} else {
				ics.clear();
				break;
			}
	}
	return ics;
}
const std::vector<AtomTopo> & AtomTopo::definedatomtypes(const std::string &filename){
	static std::vector<AtomTopo> atomtypes;
	if(atomtypes.empty()){
#ifdef _OPENMP
#pragma omp critical(atomtopoglobal)
		{
			if(atomtypes.empty()){
#endif
		std::string file=filename;
		if(file.empty()) file="proteinatomtypes.dat";
//		file=NSPdataio::datapath()+file;
		file=NSPdataio::datafilename(file);
		NSPdataio::InputLines inputlines;
		inputlines.init(file,'#');
		int natomtypes=std::stoi(inputlines[0][0]);
		atomtypes.resize(natomtypes);
		for(int i=1;i<inputlines.size();++i){
			AtomTopo at;
			if(!at.readtopo(inputlines[i])){
				std::cout<<"Error reading atomtypes from " << file <<std::endl;
				for(auto w:inputlines[i]) std::cout<<" "<<w<<std::endl;
				exit(1);
			}
			atomtypes[at.atomtype]=at;
		}
#ifdef _OPENMP
#pragma flush
		}
	}
#endif
	}
	return atomtypes;
}
bool AtomTopo::readtopo(const std::vector<std::string> &inputline){
	try{
		atomtype=std::stoi(inputline[0]);
		elementid=std::stoi(inputline[1]);
		hbdonor=(std::stoi(inputline[2])==1);
		hbacceptor=(std::stoi(inputline[3])==1);
		aromatic=(std::stoi(inputline[4])==1);
		netcharged=std::stoi(inputline[5]);
		int ncontainingres=std::stoi(inputline[6]);
		int i=7;
		while(i<7+2*ncontainingres) {
			containing_residues.push_back(inputline[i++]);
			atomnames.push_back(inputline[i++]);
		}
	}catch (std::exception &e){
		return false;
	}
	return true;
}
int AtomTopo::getatomtype(const std::string &residuename, const std::string &atomname){
	static std::map<std::string,int> atomtypemap;
	if(atomtypemap.empty()){
#ifdef _OPENMP
#pragma omp critical(atomtype_global)
		{
			if(atomtypemap.empty()){
#endif
		const std::vector<AtomTopo> & atomtypes=AtomTopo::definedatomtypes();
		for( const AtomTopo &at:atomtypes){
			for(int i=0;i<at.containing_residues.size();++i){
				std::string key=at.containing_residues[i]+at.atomnames[i];
				atomtypemap.insert(std::make_pair(key,at.atomtype));
			}
		}
#ifdef _OPENMP
#pragma omp flush
			}
		}
#endif
	}
	auto it=atomtypemap.find(residuename+atomname);
	if(it!=atomtypemap.end())return it->second;
	else{
		auto ita=atomtypemap.find("ANY"+atomname);
		if(ita!=atomtypemap.end()) return ita->second;
	}
	return -1;
}
const std::map<std::string,SideChainTopo> & SideChainTopo::definedsidechaintypes(
			const std::string & filename){
	static std::map<std::string,SideChainTopo> sidechaintopos;
	if(sidechaintopos.empty()){
#ifdef _OPENMP
#pragma omp critical(sidechainttopo_global)
		{
			if(sidechaintopos.empty()){
#endif
		std::string file=filename;
//		if(filename.empty()) file=NSPdataio::datapath()+"sidechaintopos.dat";
		if(filename.empty()) file=NSPdataio::datafilename("sidechaintopos.dat");
		NSPdataio::InputLines inputlines;
		inputlines.init(file,'#');
		int lidx=0;
		while (lidx<inputlines.size()){
			std::string resname=inputlines[lidx][0];
			int natoms=std::stoi(inputlines[lidx][1]);
			++lidx;
			sidechaintopos.insert(std::make_pair(resname,SideChainTopo()));
			SideChainTopo& st=sidechaintopos.at(resname);
			st.restype=resname;
			for(int i=0;i<natoms;++i){
				if(!st.readnewatom(inputlines[lidx++])){
					std::cout <<"Error reading atom line " <<resname<<" atom "<<i<<std::endl;
					exit(1);
				}
			}
		}
#ifdef _OPENMP
#pragma omp flush
			}
		}
#endif
	}
	return sidechaintopos;
}

bool SideChainTopo::readnewatom(const std::vector<std::string> & atomline){
	std::string atomname=atomline[0];
	int atomtype=AtomTopo::getatomtype(restype,atomname);
	if(atomtype < 0) {
		std::cout <<"Undefined atomtopo for " <<restype <<" "<<atomname<<std::endl;
		return false;
	}

	std::string upatom=atomline[1];
	if(atoms_.find(upatom) == atoms_.end() && atomname != "CB"){
		std::cout <<"Wrong order of atoms in residue topo" <<restype <<" "<<atomname<<std::endl;
		return false;
	}
	if(atomname !="CB"){
		std::vector<double> internalcrd;
		internalcrd.push_back(std::stod(atomline[2]));
		internalcrd.push_back(std::stod(atomline[3]));
		internalcrd.push_back(std::stod(atomline[4]));
		atoms_.at(upatom).downatoms.insert(std::make_pair(atomname,internalcrd));
	}
	SideChainAtomTopo sat;
	sat.atomtype=atomtype;
	sat.upatom=upatom;
	atoms_.insert(std::make_pair(atomname,sat));
	atomnames_.push_back(atomname);
	return true;
}
void FullSite::addatomcrd(const std::string &atomname,NSPgeometry::XYZ crd){
	crds_.insert(std::make_pair(atomname,crd));
}
void FullSite::updateAtomCrd(const std::string & atomName, NSPgeometry::XYZ crd) {
    crds_[atomName] = crd;
}
void FullSite::copybackbonecrd(const BackBoneSite &bs){
	crds_.erase("N");
	crds_.insert(std::make_pair("N",bs.ncrd()));
	crds_.erase("CA");
	crds_.insert(std::make_pair("CA",bs.cacrd()));
	crds_.erase("C");
	crds_.insert(std::make_pair("C",bs.ccrd()));
	crds_.erase("O");
	crds_.insert(std::make_pair("O",bs.ncrd()));
}

BackBoneSite FullSite::getbackbonesite() const{
	BackBoneSite bs;
	std::vector<NSPgeometry::XYZ> bcrd;
	bcrd.push_back(crds_.at("N"));
	bcrd.push_back(crds_.at("CA"));
	bcrd.push_back(crds_.at("C"));
	bcrd.push_back(crds_.at("O"));
	bs.changecrd(bcrd);
	bs.resname=restype_;
	bs.chainid=chainid_;
	bs.pdbid=pdbid_;
	bs.resseq=resseq_;
	bs.resid=resid_;
	bs.data_[BackBoneSite::PHI]=360.0;
	bs.data_[BackBoneSite::PSI]=360.0;
	bs.data_[BackBoneSite::OMIGA]=180.0;
	return bs;
}
double FullSite::gettorsionat(const std::string &atomname) const {
	std::string upatom1=getsidechainatomtopo(atomname).upatom;
	std::string upatom2=getsidechainatomtopo(upatom1).upatom;
	std::string upatom3=getsidechainatomtopo(upatom2).upatom;
	return NSPgeometry::torsion(crds_.at(atomname),crds_.at(upatom1),
			crds_.at(upatom2),crds_.at(upatom3))*180.0/3.14159265;
}

PdbRecord FullSite::makepdbrecord(const std::string & atomname,int atomid) const{
	PdbRecord  record;
	record.label = "ATOM";
	record.chainid = chainid_;
	record.atomname = atomname;
	record.namesymbol = record.atomname.substr(0, 1);
	record.elementname[1] = record.namesymbol[0];
	record.namemodifier = record.atomname.substr(1);
	record.residuename = restype_;
	record.atomid = atomid;
	record.residueid = resid_;
	NSPgeometry::XYZ crd = getcrd(atomname);
	record.x = crd.x_;
	record.y = crd.y_;
	record.z = crd.z_;
	return record;
}
std::string FullSite::alltopdb(int atomoffset) const{
		std::vector<std::string> allatoms=atomnames_all();
		std::vector<PdbRecord> records;
		for (auto a:allatoms) {
			if(crds_.find(a)==crds_.end()) continue;
			records.push_back(makepdbrecord(a,atomoffset++));
		}
		std::string result;
		for(auto & r:records){
			result = result + r.toString()+"\n";
		}
		return result;
}
std::vector<std::string> FullSite::atomnames_all()const {
	std::vector<std::string> allatoms;
	allatoms.push_back("N");
	allatoms.push_back("CA");
	if(SideChainTopo::residuetypedefined(restype_)) {
		for(auto &a:atomnames_sidechaintopo()){
			allatoms.push_back(a);
		}
	}
	allatoms.push_back("C");
	allatoms.push_back("O");
	return allatoms;
}
std::vector<std::string> FullSite::atomnames_withcrd()const {
	std::vector<std::string> atoms;
	if(hasatomcrd("N")) atoms.push_back("N");
	if(hasatomcrd("CA")) atoms.push_back("CA");
	auto &vsc=NSPsd::VSCType::getVSCType(restype_);
	for(auto &a:vsc.atomnames){
		if(hasatomcrd(a)) atoms.push_back(a);
	}
	if(hasatomcrd("C")) atoms.push_back("C");
	if(hasatomcrd("O")) atoms.push_back("O");
	return atoms;
}
std::vector<NSPgeometry::XYZ> FullSite::extractcrd() const{
		std::vector<std::string> allatoms=atomnames_all();
		std::vector<NSPgeometry::XYZ> crd;
		for (auto a:allatoms) {
			auto it=crds_.find(a);
			if(it==crds_.end()) continue;
			crd.push_back(it->second);
		}
		return crd;
}

std::string FullSite::sidechaintopdb(int atomoffset) const{
		std::vector<PdbRecord> records;
		for (auto &a:atomnames_sidechaintopo()) {
			if(crds_.find(a) == crds_.end()) continue;
			records.push_back(makepdbrecord(a,atomoffset++));
		}
		std::string result;
		for(auto & r:records){
			result = result + r.toString()+"\n";
		}
		return result;
}
std::vector<std::vector<FullSite>>  NSPproteinrep::readfullsitesfrompdb(
		const std::string &pdbfilename,bool keepnonstandard){
	PdbReader reader;
	reader.readpdb(pdbfilename);
	std::string chainids=reader.chainids();
	std::vector<FullSite>  sites;
	for(int i=0;i<chainids.size();++i){
		char chainid=chainids[i];
		std::vector<std::string> seq=reader.getaminoacidsequence(chainid);
		for(int resinumber=0;resinumber<seq.size();++resinumber){
			typename PdbReader::ResKeyType reskey=reader.mappdbkeyint()->pdbResKey(resinumber,i);
			std::vector<PdbRecord> &records=reader.records().at(chainid).at(reskey);
			if(!keepnonstandard && records[0].residuename=="MSE"){
				for(auto &r:records){
					r.residuename="MET";
					if(r.atomname=="SE") r.atomname="SD";
				}
			}
			if(!keepnonstandard)
				if(!SideChainTopo::residuetypedefined(records[0].residuename)) continue;
			FullSite fs;
			fs.resname()=records[0].residuename;
			fs.resid()=records[0].residueid;
			fs.chainid()=records[0].chainid;
			int  nmainchain=0;
			for(auto &r:records) {
				NSPgeometry::XYZ crd(r.x,r.y,r.z);
				if(FullSite::mainchainatom(r.atomname)) nmainchain++;
				fs.addatomcrd(r.atomname,crd);
			}
			if(nmainchain>=4) sites.push_back(fs);
		}
	}
	std::vector<std::vector<FullSite>> chains;
	chains.assign(1,std::vector<FullSite>());
	int chainnumber=0;
	int resseq=0;
	for(auto it=sites.begin();it !=sites.end()-1;++it){
		chains[chainnumber].push_back(*it);
		it->resseq()=resseq++;
		NSPgeometry::XYZ c0=it->getcrd("C");
		NSPgeometry::XYZ n1=(it+1)->getcrd("N");
		double bond = NSPgeometry::distance(c0, n1);
		if (bond > 2.0){
			chains.push_back(std::vector<FullSite>());
			++chainnumber;
			resseq=0;
		}
	}
	chains[chainnumber].push_back(sites.back());
	return chains;
}
std::vector<BackBoneSite> NSPproteinrep::extractbackbonesegment(std::vector<FullSite>::const_iterator begin, int length){
	std::vector<BackBoneSite> bsites;
	for(int i=0;i<length;++i){
		bsites.push_back(begin->getbackbonesite());
		++begin;
	}
	for(int i=0;i<length;++i){
		if(i>0) bsites[i].phi(bsites[i-1]);
		if(i<length-1) {
			bsites[i].psi(bsites[i+1]);
			bsites[i].omiga(bsites[i+1]);
		}
	}
	return bsites;
}
std::vector<std::map<std::string,double>> NSPproteinrep::calc_sasa(const std::vector<std::vector<FullSite>> &chains,int nsurfpoints) {
	std::vector<NSPgeometry::XYZ> crds;
	std::vector<double> radii;
	std::vector<std::vector<std::string>> atomnames;
	for(auto &c:chains){
		for(auto &site:c){
			site.getcrds(crds);
			site.getradii(radii);
			atomnames.push_back(site.atomnames_withcrd());
		}
	}
	std::vector<double> sasa=NSPgeometry::calc_sasa(crds,radii,1.4,nsurfpoints);
	int atomidx=0;
	int residx=0;
	std::vector<std::map<std::string,double>> result;
	for(auto &c:chains){
		for(auto &site:c){
			result.push_back(std::map<std::string,double>());
			auto &r=result.back();
			for(auto &a:atomnames.at(residx++)){
				r.insert(std::make_pair(a,sasa.at(atomidx++)));
			}
		}
	}
	return result;
}
#include "sd/sidechainff.h"
FullSite NSPproteinrep::make_fullsite(const BackBoneSite &bs,
		std::vector<std::vector<std::pair<int,double>>> internalcrds){
	FullSite result;

	result.resname()=bs.resname;
	result.resid()=bs.resid;
	result.chainid()=bs.chainid;
	result.addatomcrd("N",bs.ncrd());
	result.addatomcrd("CA",bs.cacrd());
	result.addatomcrd("C",bs.ccrd());
	result.addatomcrd("O",bs.ocrd());
	auto &vsc=NSPsd::VSCType::getVSCType(bs.resname);
	if(vsc.nscatoms==0) return result;
	std::vector<NSPgeometry::XYZ> r(4+vsc.nscatoms);
	r[0]=bs.ncrd();
	r[1]=bs.cacrd();
	r[2+vsc.nscatoms]=bs.ccrd();
	r[3+vsc.nscatoms]=bs.ocrd();
	if(internalcrds.empty()){
		internalcrds=vsc.internalcrds;
		double deg=3.14159265/180.0;
		for(int i=0;i<vsc.nscatoms;++i){
			auto &ic=internalcrds[i];
			r[i+2]=NSPgeometry::InternaltoXYZ(r[ic[0].first],r[ic[1].first],
												r[ic[2].first],
												ic[0].second, ic[1].second*deg,ic[2].second*deg);
			result.addatomcrd(vsc.atomnames[i],r[i+2]);
		}
	}
	return result;
}
#include "proteinrep/intatomkey.h"
void NSPproteinrep::writetopdb(const std::vector<std::vector<FullSite>> &sites, std::ostream &os){
	std::vector<PdbRecord> records;
	int atomseq=0;
	for(int c=0;c<sites.size();++c){
		for (int i = 0; i < sites[c].size(); ++i) {
			FullSite fs=sites[c][i];
			const NSPsd::VSCType &vsc=NSPsd::VSCType::getVSCType(fs.resname());
			std::string resname=vsc.pdbname;
			auto nkey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1, "N", c, resname,0);
			NSPgeometry::XYZ ncrd=fs.getcrd("N");
			records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(nkey,
					ncrd, ++atomseq));
			auto cakey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1,"CA", c, resname,0);
			NSPgeometry::XYZ cacrd=fs.getcrd("CA");
			records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(cakey,
					cacrd, ++atomseq));
			if (vsc.nscatoms>0) {
				for(int m=0;m<vsc.atomnames.size();++m){
					auto akey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1, vsc.atomnames[m], c,
					resname,0);
					NSPgeometry::XYZ r=fs.getcrd(vsc.atomnames[m]);
					records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(akey,
							r,++atomseq));
				}
			}
			auto ckey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1, "C", c, resname,0);
			NSPgeometry::XYZ ccrd=fs.getcrd("C");
			records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(ckey,
					ccrd, ++atomseq));
			auto okey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1, "O", c, resname,0);
			NSPgeometry::XYZ ocrd=fs.getcrd("O");
			records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(okey,
					ocrd, ++atomseq));
		}
	}
	for (auto &r : records) {
		os << r.toString()<<std::endl;
	}
}
//std::vector<BackBoneSite> backbone(const std::vector<FullSite> & fullsites);

