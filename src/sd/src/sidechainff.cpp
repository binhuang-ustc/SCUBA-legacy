/*
 * sidechainff.cpp
 *
 *  Created on: 2018年1月31日
 *      Author: hyliu
 */
#include "sd/sidechainff.h"
#include "sd/forcefield.h"
#include "dataio/inputlines.h"
#include "dataio/datapaths.h"
using namespace NSPsd;
std::map<std::string,int> VSCType::stericatomtypes;
std::map<char,std::string> VSCType::resnamefrom1letter;
std::vector<PackingAtomType> VSCType::packingatomtypes;
PackingEneFuncs VSCType::packingenefuncs;
PackingEneFuncs VSCType::softpackingenefuncs;
std::map<std::string,std::set<std::string>> VSCType::rotatablescatoms =
{
    {"ALA",{}}, {"CYS",{"SG"}}, {"ASP",{"CG","OD1"}}, {"GLU",{"CG","CD","OE1"}},
    {"PHE",{"CG","CD1"}}, {"GLY",{}}, {"HIS",{"CG","ND1"}},
    {"ILE",{"CG1","CD1"}}, {"LYS",{"CG","CD","CE","NZ"}}, {"LEU",{"CG","CD1"}},
    {"MET",{"CG","SD","CE"}}, {"ASN",{"CG","OD1"}}, {"PRO",{}},
    {"GLN", {"CG","CD","OE1"}}, {"ARG",{"CG","CD","NE","CZ"}}, {"SER",{"OG"}},
    {"THR",{"OG1"}}, {"VAL",{"CG1"}}, {"TRP",{"CG","CD1"}}, {"TYR",{"CG","CD1"}}
};
const VSCType & VSCType::getVSCType(const std::string &resname, const std::string &filename){
	static std::map<std::string,VSCType> vsctypes;
	static bool initialized{false};
	if(!initialized){
#ifdef _OPENMP
#pragma omp critical(vsctype_global)
		{
			if(!initialized){
#endif
		std::string file=filename;
		if(file.empty()) file=std::string("sidechainff.dat");
//		file =NSPdataio::datapath()+file;
		file=NSPdataio::datafilename(file);
		vsctypes=readVSCTypes(file);
		initialized=true;
#ifdef _OPENMP
#pragma omp flush
			}
		}
#endif
	}
	return vsctypes.at(resname);
}
std::map<std::string,VSCType> VSCType::readVSCTypes(const std::string & filename){
	std::map<std::string,VSCType> readtypes;
	NSPdataio::InputLines inputlines;
	inputlines.init(filename,'#');
	int lidx=0;
	std::vector<std::string> &words=inputlines[lidx++];
	int nstericatomtypes=std::stoi(words[0]);
	packingatomtypes.resize(2*nstericatomtypes);
	for(int i=0;i<nstericatomtypes;++i){
		words=inputlines[lidx++];
		int widx=0;
		int t=std::stoi(words[widx++]);
		assert(t<2*nstericatomtypes);
		packingatomtypes[t].radius=std::stod(words[widx++]);
		packingatomtypes[t].hbtype=std::stoi(words[widx++]);
		packingatomtypes[t].aromatic=std::stoi(words[widx++]);
		packingatomtypes[t].nconnections=std::stoi(words[widx++]);
		for(int w=widx;w<words.size();++w){
			std::string key=words[w];
			stericatomtypes.insert(std::make_pair(key,t));
		}
	}
	packingenefuncs.setup(nstericatomtypes);
	softpackingenefuncs.setup(nstericatomtypes);
//	std::shared_ptr<EneFunc1D> hbfunc(new HarmGauEne(0.3,0.0,0.01,0.3,6));
	std::shared_ptr<EneFunc1D> hbfunc(new LJGauEne(0.29,0.0,0.01));
	std::shared_ptr<EneFunc1D> softhbfunc(new SwitchGauEne(0.3,0.0,0.21,10.0,0.14));
	for(int i=0;i<nstericatomtypes;++i){
		auto & pi=packingatomtypes[i];
		bool inpolar=(pi.hbtype==0);
		for(int j=i;j<nstericatomtypes;++j){
			auto & pj=packingatomtypes[j];
			if(PackingAtomType::hbond(pi.hbtype,pj.hbtype)){
				packingenefuncs.addfunction(i,j,hbfunc);
				packingenefuncs.addfunction(j,i,hbfunc);
				softpackingenefuncs.addfunction(i,j,softhbfunc);
				softpackingenefuncs.addfunction(j,i,softhbfunc);
				continue;
			}
			bool jnpolar=(pj.hbtype==0);
			double rmin=0.1*(packingatomtypes[i].radius+packingatomtypes[j].radius);
			double emin;
			if(inpolar && jnpolar){
				if(pi.aromatic!=0 || pj.aromatic !=0) emin=1.0;
				else emin=1.0;
			} else if (inpolar || jnpolar){
				emin=0.2;
			} else {
				emin=0.0;
			}
			double ncon=pi.nconnections+pj.nconnections;
			if(ncon>2) emin=emin*1.5/(ncon-1.0);
//			std::shared_ptr<EneFunc1D> func(new HarmGauEne(rmin,emin,0.09,0.14,6));
			std::shared_ptr<EneFunc1D> func(new LJGauEne(rmin,emin,0.14));
			packingenefuncs.addfunction(i,j,func);
			packingenefuncs.addfunction(j,i,func);
			double rcore=rmin*0.7;
			double emax=10.0;
			double emin_s=-1.0;
			std::shared_ptr<EneFunc1D> softfunc(new SwitchGauEne(rmin,emin_s,rcore,emax,0.14));
			softpackingenefuncs.addfunction(i,j,softfunc);
			softpackingenefuncs.addfunction(j,i,softfunc);
		}
	}
	while (lidx<inputlines.size()){
		words=inputlines[lidx++];
		std::string resname=words[0];
		readtypes.insert(std::make_pair(resname,VSCType()));
		VSCType &vsc=readtypes.at(resname);
		vsc.resname=resname;
		vsc.pdbname=words[1];
		vsc.oneletter=words[2][0];
		resnamefrom1letter.insert(std::make_pair(vsc.oneletter,vsc.resname));
		int nscatoms=std::stoi(words[3]);
		vsc.nscatoms=nscatoms;
		if(nscatoms==0) continue;
		for(int i=0;i<nscatoms;++i){
			int widx=0;
			words=inputlines[lidx++];
			std::string atomname=words[widx++];
			int atype=getstericatomtype(resname,atomname);
			if(atype <0){
				std::cout <<"Steric atomtype for "<<resname<<":"<<atomname <<" is not defined." <<std::endl;
				exit(1);
			}
			double sigma=std::stod(words[widx++]);
			int hbtype=std::stoi(words[widx++]);
			assert(hbtype==packingatomtypes[atype].hbtype);
			int rotameratom=std::stoi(words[widx++]);
			int aj=std::stoi(words[widx++]);
			double b=std::stod(words[widx++]);
			double ak=std::stoi(words[widx++]);
			double a=std::stod(words[widx++]);
			double al=std::stoi(words[widx++]);
			double p=std::stod(words[widx++]);
			vsc.atomnames.push_back(atomname);
			if(rotameratom==1) vsc.rotameratoms.push_back(i);
			std::vector<std::pair<int,double>> intcrd;
			intcrd.push_back(std::make_pair(aj,b));
			intcrd.push_back(std::make_pair(ak,a));
			intcrd.push_back(std::make_pair(al,p));
			vsc.internalcrds.push_back(intcrd);
		}
		words=inputlines[lidx++];
		int widx=0;

		for(int i=0;i<words.size()/4;++i){
			int ai=std::stoi(words[widx++]);
			int aj=std::stoi(words[widx++]);
			double b=std::stod(words[widx++]);
			double b2=std::stod(words[widx++]);
			b2=2*KBT/(b2*b2);
			vsc.newbonds.push_back(std::make_pair(ai,aj));
			vsc.b0.push_back(b);
//			vsc.kb0.push_back(b2);
			vsc.kb0.push_back(1000.0);
		}
		words=inputlines[lidx++];
		widx=0;
//		std::cout <<resname <<std::endl;
		for(int i=0;i<words.size()/5;++i){
			std::vector<int> aijk;
			for(int m=0;m<3;++m) aijk.push_back(std::stoi(words[widx++]));
			vsc.newangles.push_back(aijk);
			double a=std::stod(words[widx++]);
			double a2=std::stod(words[widx++]);
			a2=2.*KBT*KANG_FAC/(9.0*a2*a2);
			vsc.a0.push_back(a);
			if(a2>2000.0) a2=2000.0;
			vsc.ka0.push_back(a2);
//			std::cout <<"ka "<< a2<<std::endl;
		}
		words=inputlines[lidx++];
		widx=0;
		for(int i=0;i<words.size()/6;++i){
			std::vector<int> aijkl;
			for(int m=0;m<4;++m) aijkl.push_back(std::stoi(words[widx++]));
			vsc.newimpdihs.push_back(aijkl);
			double p=std::stod(words[widx++]);
			double p2=std::stod(words[widx++]);
			p2=2.0*KBT*KANG_FAC/(9.0*p2*p2);
			vsc.imp0.push_back(p);
			vsc.kimp0.push_back(p2);
			if(p2>2000.0) p2=2000.0;
//			std::cout <<"kimp " <<p2 <<std::endl;
		}
		words=inputlines[lidx++];
		widx=0;
		for(int i=0;i<words.size()/4;++i){
			std::vector<int> aijkl;
			for(int m=0;m<4;++m) aijkl.push_back(std::stoi(words[widx++]));
			vsc.newtorsions.push_back(aijkl);
		}
	}
	return readtypes;
}

std::vector<double> ConformerCode::gettorsioncodes(double ang,
	const std::vector<DvDxi> &dadx, std::vector<std::vector<DvDxi>> *dcdx){
	const std::vector<double> nang{1.0,2.0,4.0};
	std::vector<double> res;
	dcdx->clear();
	for( double n:nang){
		double c=cos(n*ang);
		double s=sin(n*ang);
		res.push_back(c);
		dcdx->push_back(std::vector<DvDxi>());
		auto & dccdx=dcdx->back();
		for(auto &d:dadx)dccdx.push_back(-(s*n)*d);
		res.push_back(s);
		dcdx->push_back(std::vector<DvDxi>());
		auto & dcsdx=dcdx->back();
		for(auto &d:dadx)dcsdx.push_back((c*n)*d);
	}
	return res;
}
std::vector<ConformerCode> NSPsd::makeconformercodes(const std::vector<double> &crds,const std::vector<SCInChain> &scinchains,
		std::vector<PhiPsiCodes> &phipsicodes){
		int nscs=scinchains.size();
		std::vector<ConformerCode> result(nscs);
#pragma omp parallel for schedule(dynamic,1)
		for(int i=0;i<nscs;++i){
			if(scinchains[i].kaiatoms.empty()) continue;
			ConformerCode &cc=result[i];
			cc.restype=scinchains[i].restype;
			cc.phipsicodes=&(phipsicodes[i]);
			for(auto & ijkl:scinchains[i].kaiatoms){
				std::vector<NSPgeometry::XYZ> dtdx;
				double kai=NSPgeometry::torsion(getxyz(crds,ijkl[0]),
						getxyz(crds,ijkl[1]),getxyz(crds,ijkl[2]),getxyz(crds,ijkl[3]),&dtdx);
				cc.sidechaintorsions.push_back(kai);
				std::vector<DvDxi> dadx;
				dadx.push_back(std::make_pair(ijkl[0],dtdx[0]));
				dadx.push_back(std::make_pair(ijkl[1],dtdx[1]));
				dadx.push_back(std::make_pair(ijkl[2],dtdx[2]));
				dadx.push_back(std::make_pair(ijkl[3],dtdx[3]));
				std::vector<std::vector<DvDxi>> dcdx;
				cc.dtorsioncodesdx.push_back(std::vector<std::vector<DvDxi>>());
				cc.torsioncodes.push_back(
						ConformerCode::gettorsioncodes(kai,dadx,&(cc.dtorsioncodesdx.back())));
			}
		}
		return result;
}
double NSPsd::mcscpackingenergy(const std::vector<NSPgeometry::XYZ> &crds, const std::vector<int> &atomtypes,
		const BSInChain &mc,const SCInChain &sc,
		int sep,std::vector<DvDxi> *dedx,bool terminal){
	dedx->clear();
	double ene=0.0;
	int nscatoms=sc.nscatoms;
	if(nscatoms==0) return ene;
	std::vector<int> mcatoms;
	if(sep==1){
		mcatoms.push_back(mc.nid);
		if(terminal) mcatoms.push_back(mc.oid);
	} else if (sep==-1){
		if(terminal) mcatoms.push_back(mc.nid);
		mcatoms.push_back(mc.cid);
		mcatoms.push_back(mc.oid);
	} else {
		mcatoms=mc.atomids();
	}
	std::vector<NSPgeometry::XYZ> dvdmc(mcatoms.size(),NSPgeometry::XYZ(0,0,0));
	std::vector<NSPgeometry::XYZ> dvdsc(nscatoms,NSPgeometry::XYZ(0,0,0));
	bool soft=sc.softpacking;
	for(int i=0;i<mcatoms.size();++i){
		int im=mcatoms[i];
		for(int j=0;j<nscatoms;++j){
			int js=sc.poffset+j;
			std::vector<NSPgeometry::XYZ> drdx;
			double r=NSPgeometry::distance(crds[im],crds[js],&drdx);
			double  dedr;
			if(!soft){
				ene +=VSCType::packingenefuncs.getenefunc(atomtypes[im],atomtypes[js]).energy(r,&dedr);
			} else {
				ene +=VSCType::softpackingenefuncs.getenefunc(atomtypes[im],atomtypes[js]).energy(r,&dedr);
			}
			dvdmc[i] = dvdmc[i]+dedr*drdx[0];
			dvdsc[j] = dvdsc[j]+dedr*drdx[1];
		}
	}
	for(int i=0;i<mcatoms.size();++i) (*dedx).push_back(std::make_pair(mcatoms[i],dvdmc[i]));
	for(int j=0;j<nscatoms;++j) (*dedx).push_back(std::make_pair(sc.poffset+j,dvdsc[j]));
	return ene;
}
double NSPsd::scscpackingenergy(const std::vector<NSPgeometry::XYZ> &crds,
		const std::vector<int> &atomtypes,
		const SCInChain &sc1, const SCInChain &sc2,
		int sep,std::vector<DvDxi> *dedx){
	dedx->clear();
	double ene=0.0;
	int nsc1=sc1.nscatoms;
	int nsc2=sc2.nscatoms;
	if(nsc1==0 ||nsc2==0) return ene;
	std::vector<NSPgeometry::XYZ> dvdx1(nsc1,NSPgeometry::XYZ(0.0,0.0,0.0));
	std::vector<NSPgeometry::XYZ> dvdx2(nsc2,NSPgeometry::XYZ(0.0,0.0,0.0));
	bool soft=sc1.softpacking ||sc2.softpacking;
	for(int i=0;i<nsc1;++i){
		int ia=sc1.poffset+i;
		for(int j=0;j<nsc2;++j){
			int ja=sc2.poffset+j;
			std::vector<NSPgeometry::XYZ> drdx;
			double r=NSPgeometry::distance(crds[ia],crds[ja],&drdx);
			double dedr;
			if(!soft){
				if(sep==1 || sep==-1)
					ene += VSCType::packingenefuncs.getenefunc(atomtypes[ia],atomtypes[ja]).energyalt(r,&dedr);
				else
					ene += VSCType::packingenefuncs.getenefunc(atomtypes[ia],atomtypes[ja]).energy(r,&dedr);
			} else {
				if(sep==1 ||sep==-1) continue;
				ene += VSCType::softpackingenefuncs.getenefunc(atomtypes[ia],atomtypes[ja]).energy(r,&dedr);
			}
			dvdx1[i]=dvdx1[i]+dedr*drdx[0];
			dvdx2[j]=dvdx2[j]+dedr*drdx[1];
		}
	}
	for(int i=0;i<nsc1;++i) (*dedx).push_back(std::make_pair(sc1.poffset+i,dvdx1[i]));
	for(int j=0;j<nsc2;++j) (*dedx).push_back(std::make_pair(sc2.poffset+j,dvdx2[j]));
	return ene;
}
