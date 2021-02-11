/*
 * anlytsasa.cpp
 *
 *  Created on: 2018年3月23日
 *      Author: hyliu
 */
#include "sd/anlytsasa.h"
#include "dataio/datapaths.h"
using namespace NSPsd;
using namespace NSPgeometry;

static double bij_r(double ri,double rj, double solv,double rij,double *dbijdr){
	double rsum=ri+rj+2*solv;
	if(rij>=rsum) return 0.0;
	double res=3.14159265*(ri+solv)*(rsum-rij)*(1+(rj-ri)/rij);
	*dbijdr=3.14159265*(ri+solv)*(-1-(rj-ri)/rij-(rsum-rij)*(rj-ri)/(rij*rij));
	return res;
}
void Anlytsasa::loopneighbors(const std::vector<NSPgeometry::XYZ> &xyz,
		const std::vector<std::vector<int>> &neighbors,
		double p,
		std::vector<std::vector<double>> &facti,
		std::vector<std::vector<NSPgeometry::XYZ>> &dfactidxi,
		std::vector<std::vector<DvDxi>> & dfactidxj){
	int natoms=xyz.size();
	for(int i=0;i<natoms;++i){
		for(auto j:neighbors[i]){
			if(j<i) continue;
			std::vector<XYZ> drdx;
			double rij=distance(xyz[i],xyz[j],&drdx);
			double dbijdr;
			double bij=bij_r(ri[i],ri[j],r_solv,rij,&dbijdr);
			if(bij<=0.0) continue;
			double mt=pi[i]*p/si[i];
			facti[i].push_back(1-mt*bij);
			dfactidxi[i].push_back(-(mt*dbijdr)*drdx[0]);
			dfactidxj[i].push_back(std::make_pair(j,-mt*dbijdr*drdx[1]));
			double dbjidr;
			double bji=bij_r(ri[j],ri[i],r_solv,rij,&dbjidr);
			mt=pi[j]*p/si[j];
			facti[j].push_back(1-mt*bji);
			dfactidxi[j].push_back(-(mt*dbjidr)*drdx[1]);
			dfactidxj[j].push_back(std::make_pair(i,-mt*dbjidr*drdx[0]));
		}
	}
}
Anlytsasa::Result Anlytsasa::calc_anlytsasa(const std::vector<XYZ> &xyz,
		const NeighborList &nbl){
	Anlytsasa::Result result;
	result.atomic_sasa.resize(xyz.size());
	result.dsadx.assign(xyz.size(),std::vector<DvDxi>());
	int natoms=ri.size();
	std::vector<std::vector<double>> facti(natoms);
	std::vector<std::vector<DvDxi>> dfactidxj(natoms);
	std::vector<std::vector<XYZ>> dfactidxi(natoms);
	loopneighbors(xyz,nbor12,p12,facti,dfactidxi,dfactidxj);
	loopneighbors(xyz,nbor13,p13,facti,dfactidxi,dfactidxj);
	loopneighbors(xyz,nbor14,p14,facti,dfactidxi,dfactidxj);
	loopneighbors(xyz,nbl.neighbors,pij,facti,dfactidxi,dfactidxj);
	for(int i=0;i<natoms;++i){
		double fm=si[i];
		for(auto fac:facti[i]) fm*=fac;
		if(fm>0){
			XYZ dsadxi;
			for(int k=0;k<facti[i].size();++k){
				double fother=fm/facti[i][k];
				dsadxi = dsadxi+fother*dfactidxi[i][k];
				result.dsadx[i].push_back(fother*dfactidxj[i][k]);
			}
			result.dsadx[i].push_back(std::make_pair(i,dsadxi));
		} else{
			fm=0.0;
		}
		result.atomic_sasa[i]=fm;
	}
	return result;
}
Anlytsasa::Anlytsasa(const ForceField&ff){
//	std::string filename=NSPdataio::datapath()+"analytsasa.par";
	std::string filename = NSPdataio::datafilename("analytsasa.par");
	std::ifstream ifs;
	ifs.open(filename.c_str());
	ifs >>p12 >>p13 >>p14 >>pij >>r_solv;
	r_solv*=A2NM;
	int natomtypes;
	ifs >>natomtypes;
	std::map<std::string,double> r_atom;
	std::map<std::string,double> p0_atom;
	std::string resname;
	std::string atomname;
	for(int i=0;i<natomtypes;++i){
		double ra,p0a;
		ifs >>resname >>atomname>>ra>>p0a;
		ra*=A2NM;
		std::string key=resname+"_"+atomname;
		r_atom.insert(std::make_pair(key,ra));
		p0_atom.insert(std::make_pair(key,p0a));
	}
	nbor12.assign(ff.natoms(),std::vector<int>());
	for(auto &b:ff.bondterms()) {
		int i=b.i;
		int j=b.j;
		if(i>j){
			int tmp=i;
			i=j;
			j=tmp;
		}
		nbor12[i].push_back(j);
	}
	nbor13.assign(ff.natoms(),std::vector<int>());
	for(auto &a:ff.angleterms()){
		int i=a.i;
		int k=a.k;
		if(i>k){
			int tmp=i;
			i=k;
			k=tmp;
		}
		nbor13[i].push_back(k);
	}
	nbor14.assign(ff.natoms(),std::vector<int>());
	int aidx=0;
	for(auto &l14:ff.list14()){
		for(int i:l14){
			if(i>aidx) nbor14[aidx].push_back(i);
		}
		aidx++;
	}
	int nchains=ff.bsinchains().size();
	ri.assign(ff.natoms(),0.16);
	pi.assign(ff.natoms(),1.10);
	for(int c=0;c<nchains;++c){
		int nresidues=ff.bsinchains().at(c).size();
		for(int r=0;r<nresidues;++r){
			const BSInChain & bs=ff.bsinchains().at(c).at(r);
			const SCInChain & sc=ff.scinchains().at(c).at(r);
			std::string resname=sc.restype;
			std::string atom=resname+"_N";
			ri[bs.nid]=r_atom.at(atom);
			pi[bs.nid]=p0_atom.at(atom);
			atom=resname+"_CA";
			ri[bs.caid]=r_atom.at(atom);
			pi[bs.caid]=p0_atom.at(atom);
			atom=resname+"_C";
			ri[bs.cid]=r_atom.at(atom);
			pi[bs.cid]=p0_atom.at(atom);
			atom=resname+"_O";
			ri[bs.oid]=r_atom.at(atom);
			pi[bs.oid]=p0_atom.at(atom);
			const VSCType & vsc=VSCType::getVSCType(resname);
			for(int a=0;a<vsc.nscatoms;++a){
				atom=resname+"_"+vsc.atomnames[a];
				ri[sc.poffset+a]=r_atom.at(atom);
				pi[sc.poffset+a]=p0_atom.at(atom);
			}
		}
	}
	double fourpi=4*3.14159265;
	si.clear();
	for(auto r:ri){
		double ris=r+r_solv;
		si.push_back(fourpi*ris*ris);
	}
}
