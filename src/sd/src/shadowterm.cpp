/*
 * shawdowterm.cpp
 *
 *  Created on: 2018年2月4日
 *      Author: hyliu
 */

#include "sd/shadowterm.h"

using namespace NSPsd;
using namespace NSPgeometry;
using namespace NSPpdbstatistics;
std::vector<EneFunc1D *> ShadowTerm::enefuncs;
std::map<std::string, ShiftedG> ShadowTerm::shiftgfuncs;
std::map<std::string, SineEne> ShadowTerm::sinefuncs;
NullEne1D ShadowTerm::nullene1d;
//std::vector<std::string> allatoms { "N", "CA", "C", "O", "VB", "VG1", "VG2",
//		"VD1", "VD2" };
std::vector<std::string> mcatoms { "N", "CA", "C"};
std::vector<std::string> vscatoms { "VB", "VG1", "VG2", "VD1", "VD2" };
/*double ShadowTerm::atompairenergy(const std::string &a1, const std::string &a2,
		double r, double *dedr) {
	static bool initialized { false };
	if (!initialized) {
		setupenefuncs();
		initialized=true;
	}
	auto found=enefuncs.find(a1+a2);
	if(found ==enefuncs.end()){
		*dedr=0.0;
		return 0.0;
	}
	return found->second->energy(r,dedr);
}*/
double ShadowTerm::sitepairenergy(const std::vector<NSPgeometry::XYZ> & bcrd1,
		const NSPpdbstatistics::VirtualSideChain & vsc1,
		const std::vector<NSPgeometry::XYZ> &bcrd2,const NSPpdbstatistics::VirtualSideChain &vsc2,
		std::vector<NSPgeometry::XYZ> *forces1, std::vector<NSPgeometry::XYZ> *forces2){
	double rca2=distance(bcrd1[1],bcrd2[1]);
	forces1->assign(mcatoms.size(),XYZ());
	forces2->assign(mcatoms.size(),XYZ());
	if(rca2>1.69) return 0.0;
	double ene=0;
	int eidx=-1;
	for(int i=0;i<mcatoms.size();++i){
		for(int j=0;j<vscatoms.size();++j){
			++eidx;
			std::vector<XYZ> drdx;
			double r=distance(bcrd1[i],vsc2.getcrd(vscatoms[j]),&drdx);
			double dedr;
			double eap=enefuncs[eidx]->energy(r,&dedr);
			if(eap==0.0) continue;
			(*forces1)[i]=(*forces1)[i]-dedr*drdx[0];
			std::vector<XYZ> fdistr=vsc2.redistributederiv(vscatoms[j],-dedr*drdx[1]);
			(*forces2)[1]=(*forces2)[1]+fdistr[0]; //CA
			(*forces2)[0]=(*forces2)[0]+fdistr[1]; //N
			(*forces2)[2]=(*forces2)[2]+fdistr[2]; //C
			ene+=eap;
		}
	}
	eidx=-1;
	for(int i=0;i<mcatoms.size();++i){
		for(int j=0;j<vscatoms.size();++j){
			++eidx;
			std::vector<XYZ> drdx;
			double r=distance(bcrd2[i],vsc1.getcrd(vscatoms[j]),&drdx);
			double dedr;
			double eap=enefuncs[eidx]->energy(r,&dedr);
			if(eap==0.0) continue;
			(*forces2)[i]=(*forces2)[i]-dedr*drdx[0];
			std::vector<XYZ> fdistr=vsc1.redistributederiv(vscatoms[j],-dedr*drdx[1]);
			(*forces1)[1]=(*forces1)[1]+fdistr[0]; //CA
			(*forces1)[0]=(*forces1)[0]+fdistr[1]; //N
			(*forces1)[2]=(*forces1)[2]+fdistr[2]; //C
			ene+=eap;
		}
	}
	for(int i=0;i<vscatoms.size();++i){
		for(int j=0;j<vscatoms.size();++j){
			++eidx;
			std::vector<XYZ> drdx;
			double r=distance(vsc1.getcrd(vscatoms[i]),vsc2.getcrd(vscatoms[j]),&drdx);
			double dedr;
			double eap=enefuncs[eidx]->energy(r,&dedr);
			if(eap==0.0) continue;
			std::vector<XYZ> fdistr=vsc1.redistributederiv(vscatoms[i],-dedr*drdx[0]);
			(*forces1)[1]=(*forces1)[1]+fdistr[0]; //CA
			(*forces1)[0]=(*forces1)[0]+fdistr[1]; //N
			(*forces1)[2]=(*forces1)[2]+fdistr[2]; //C
			fdistr=vsc2.redistributederiv(vscatoms[j],-dedr*drdx[1]);
			(*forces2)[1]=(*forces2)[1]+fdistr[0]; //CA
			(*forces2)[0]=(*forces2)[0]+fdistr[1]; //N
			(*forces2)[2]=(*forces2)[2]+fdistr[2]; //C
			ene+=eap;
		}
	}
	return ene;
}
double ShadowTerm::sitepairenergy(const std::vector<NSPgeometry::XYZ> & bcrd1,
		const NSPpdbstatistics::VirtualSideChain & vsc1,
		const std::vector<NSPgeometry::XYZ> & crds, const BSInChain & bs2, const SCInChain & sc2,
		std::vector<NSPgeometry::XYZ> *forces1, std::vector<DvDxi> *dedx2){
	int eidx=-1;
	forces1->assign(mcatoms.size(),XYZ(0,0,0));
	dedx2->clear();
	double rca2=(bcrd1[1]-crds[bs2.caid]).squarednorm();
	if(rca2>1.69) return 0.0;
	std::vector<int> mc2=bs2.atomids();
	double ene=0.0;
	for(int i=0;i<mcatoms.size();++i){
		int ia=mc2[i];
		std::vector<NSPgeometry::XYZ> f2(3,XYZ(0,0,0));
		for(int j=0; j<vscatoms.size();++j){
//		for(int j=0;j<1;++j){
			++eidx;
			std::vector<XYZ> drdx;
			double r=distance(crds[ia],vsc1.getcrd(vscatoms[j]),&drdx);
			double dedr;
			double eap=enefuncs[eidx]->energy(r,&dedr);
			if(eap==0.0) continue;
			f2[i]=f2[i]+dedr*drdx[0];
			std::vector<XYZ> fdistr=vsc1.redistributederiv(vscatoms[j],-dedr*drdx[1]);
			(*forces1)[1]=(*forces1)[1]+fdistr[0]; //CA
			(*forces1)[0]=(*forces1)[0]+fdistr[1]; //N
			(*forces1)[2]=(*forces1)[2]+fdistr[2]; //C
			ene+=eap;
		}
		dedx2->push_back(std::make_pair(ia,f2[i]));
	}
	auto & vsc=VSCType::getVSCType(sc2.restype);
	for(int i=0;i<sc2.nscatoms;++i){
		int ia=i+sc2.poffset;
		char posi=vsc.atomnames[i][1];
		std::string asc;
		if(posi=='B') asc="VB";
		else if(posi=='G') asc="VD1";
		else asc="VG1";
		std::vector<NSPgeometry::XYZ> f2(sc2.nscatoms,XYZ(0,0,0));
		for(int j=0; j<vscatoms.size();++j){
//		for(int j=0;j<1;++j){
			std::vector<XYZ> drdx;
			double r=distance(crds[ia],vsc1.getcrd(vscatoms[j]),&drdx);
			double dedr;
			auto found=shiftgfuncs.find(asc+vscatoms[j]);
			if(found ==shiftgfuncs.end()) found=shiftgfuncs.find(vscatoms[j]+asc);
			double eap=found->second.energy(r,&dedr);
			if(eap==0.0) continue;
			f2[i]=f2[i]+dedr*drdx[0];
			std::vector<XYZ> fdistr=vsc1.redistributederiv(vscatoms[j],-dedr*drdx[1]);
			(*forces1)[1]=(*forces1)[1]+fdistr[0]; //CA
			(*forces1)[0]=(*forces1)[0]+fdistr[1]; //N
			(*forces1)[2]=(*forces1)[2]+fdistr[2]; //C
			ene+=eap;
		}
		dedx2->push_back(std::make_pair(ia,f2[i]));
	}
	return ene;
}
void ShadowTerm::setupenefuncs() {
#ifdef _OPENMP
	static bool setupdone{false};
#pragma omp critical(shadowene_global)
	{
		if(!setupdone){
#endif
	shiftgfuncs.insert(std::make_pair(std::string("NVB"), ShiftedG(0.29, 0.11, 2.8)));
	shiftgfuncs.insert(std::make_pair(std::string("NVG1"), ShiftedG(0.18, 0.18, 2.7)));
	shiftgfuncs.insert(std::make_pair(std::string("NVG2"), ShiftedG(0.18, 0.18, 2.7)));
	shiftgfuncs.insert(std::make_pair(std::string("NVD1"), ShiftedG(0.12, 0.23, 1.9)));
	shiftgfuncs.insert(std::make_pair(std::string("NVD2"), ShiftedG(0.12, 0.23, 1.9)));
	shiftgfuncs.insert(std::make_pair(std::string("CAVB"), ShiftedG(0.29, 0.11, 2.3)));
	shiftgfuncs.insert(std::make_pair(std::string("CAVG1"), ShiftedG(0.18, 0.18, 2.0)));
	shiftgfuncs.insert(std::make_pair(std::string("CAVG2"), ShiftedG(0.13,0.21, 1.8)));
	shiftgfuncs.insert(std::make_pair(std::string("CAVD1"), ShiftedG(0.10, 0.22, 1.8)));
	shiftgfuncs.insert(std::make_pair(std::string("CAVD2"), ShiftedG(0.10, 0.22, 1.8)));
	shiftgfuncs.insert(std::make_pair(std::string("CVB"), ShiftedG(0.29, 0.11, 2.3)));
	shiftgfuncs.insert(std::make_pair(std::string("CVG1"), ShiftedG(0.18, 0.18, 2.0)));
	shiftgfuncs.insert(std::make_pair(std::string("CVG2"), ShiftedG(0.13, 0.21, 1.8)));
	shiftgfuncs.insert(std::make_pair(std::string("CVD1"), ShiftedG(0.10, 0.22, 1.8)));
	shiftgfuncs.insert(std::make_pair(std::string("CVD2"), ShiftedG(0.10, 0.22, 1.8)));
	shiftgfuncs.insert(std::make_pair(std::string("VBVB"), ShiftedG(0.21, 0.13, 2.6)));
	shiftgfuncs.insert(std::make_pair(std::string("VBVG1"), ShiftedG(0.11, 0.15, 2.5)));
	shiftgfuncs.insert(std::make_pair(std::string("VBVG2"), ShiftedG(0.11, 0.15, 2.5)));
	shiftgfuncs.insert(std::make_pair(std::string("VBVD1"), ShiftedG(0.05, 0.18, 1.7)));
	shiftgfuncs.insert(std::make_pair(std::string("VBVD2"), ShiftedG(0.05, 0.18, 1.7)));
	shiftgfuncs.insert(std::make_pair(std::string("VG1VG1"), ShiftedG(0.12, 0.11, 1.5)));
	shiftgfuncs.insert(std::make_pair(std::string("SCVB"),ShiftedG(0.25,0.05,5.0)));
	shiftgfuncs.insert(
			std::make_pair(std::string("VG1VG2"), ShiftedG(0.12, 0.11, 1.5)));
	shiftgfuncs.insert(
			std::make_pair(std::string("VG1VD1"), ShiftedG(0.12, 0.1, 0.5)));
	shiftgfuncs.insert(
			std::make_pair(std::string("VG1VD2"), ShiftedG(0.12, 0.1, 0.5)));
	shiftgfuncs.insert(
			std::make_pair(std::string("VG2VG2"), ShiftedG(0.12, 0.11, 1.5)));
	shiftgfuncs.insert(
			std::make_pair(std::string("VG2VD1"), ShiftedG(0.12, 0.1, 0.5)));
	shiftgfuncs.insert(
			std::make_pair(std::string("VG2VD2"), ShiftedG(0.12, 0.1, 0.5)));
	sinefuncs.insert(
			std::make_pair(std::string("VD1VD1"), SineEne(-0.15, 1.0)));
	sinefuncs.insert(
			std::make_pair(std::string("VD1VD2"), SineEne(-0.15, 1.0)));
	sinefuncs.insert(
			std::make_pair(std::string("VD2VD2"), SineEne(-0.15, 1.0)));
	int enefuncsize=(mcatoms.size()+vscatoms.size())*vscatoms.size();
	enefuncs.assign(enefuncsize,nullptr);
	int eidx=-1;
	for (int i = 0; i < mcatoms.size(); ++i) {
		for (int j = 0; j < vscatoms.size(); ++j) {
			++eidx;
			std::string entry = mcatoms[i] + vscatoms[j];
//			if (enefuncs.find(entry) != enefuncs.end())
//				continue;
			std::string entrys = vscatoms[j] + mcatoms[i];
			auto found = shiftgfuncs.find(entry);
			if (found == shiftgfuncs.end())
				found = shiftgfuncs.find(entrys);
			if (found != shiftgfuncs.end()) {
				enefuncs[eidx]=&(found->second);
			} else {
				auto founds = sinefuncs.find(entry);
				if (founds == sinefuncs.end())
					founds = sinefuncs.find(entrys);
				if (founds != sinefuncs.end()) {
					enefuncs[eidx]=&(founds->second);
				}
			}
			if(enefuncs[eidx]==nullptr) enefuncs[eidx]=&nullene1d;
		}
	}
	for (int i = 0; i < vscatoms.size(); ++i) {
		for (int j = 0; j < vscatoms.size(); ++j) {
			++eidx;
			std::string entry = vscatoms[i] + vscatoms[j];
//			if (enefuncs.find(entry) != enefuncs.end())
//				continue;
			std::string entrys = vscatoms[j] + vscatoms[i];
			auto found = shiftgfuncs.find(entry);
			if (found == shiftgfuncs.end())
				found = shiftgfuncs.find(entrys);
			if (found != shiftgfuncs.end()) {
				enefuncs[eidx]=&(found->second);
			} else {
				auto founds = sinefuncs.find(entry);
				if (founds == sinefuncs.end())
					founds = sinefuncs.find(entrys);
				if (founds != sinefuncs.end()) {
					enefuncs[eidx]=&(founds->second);
				}
			}
			if(enefuncs[eidx]==nullptr) enefuncs[eidx]=&nullene1d;
		}
	}
#ifdef _OPENMP
		setupdone=true;
#pragma omp flush
		}
	}
#endif
}
