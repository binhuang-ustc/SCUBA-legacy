/*
 * bsinchain.cpp
 *
 *  Created on: 2017年12月13日
 *      Author: hyliu
 */

#include "sd/bsinchain.h"
#include "sd/nnterm.h"
#include "geometry/calculators.h"
#include <cmath>
using namespace NSPsd;
using namespace NSPgeometry;

std::vector<BSInChain> NSPsd::makeBSInChain(int nsites, int idx_offset){
	std::vector<BSInChain> result;
	for(int i=0;i<nsites;++i){
		int n=idx_offset+4*i;
		result.push_back(BSInChain(n,n+1,n+2,n+3));
	}
	return result;
}
std::vector<PhiPsiCodes> NSPsd::makephipsicodes(const std::vector<double> &crds,
		const std::vector<BSInChain> & bsinchain){
	int nsite=bsinchain.size();
	std::vector<PhiPsiCodes> result(nsite);
#pragma omp parallel for schedule(dynamic,1)
	for(int i=0;i<nsite;++i){
//	for(auto it=bsinchain.begin(); it!=bsinchain.end();++it){
		auto it=bsinchain.begin()+i;
		result[i]=PhiPsiCodes(crds,it,it==bsinchain.begin(), it ==
				bsinchain.begin()+bsinchain.size()-1);
	}
	return result;
}
std::vector<double> PhiPsiCodes::gettriangcodes(double ang,
	const std::vector<DvDxi> &dadx, std::vector<std::vector<DvDxi>> *dcdx){
	const std::vector<double> nang{1.0,2.0,4.0,8.0};
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
PhiPsiCodes::PhiPsiCodes(const std::vector<double> &crd, std::vector<BSInChain>::const_iterator bs,
		bool Nterm,bool Cterm){
	if(!Nterm){
		std::vector<XYZ> dphi;
		phi=torsion(getxyz(crd,(bs-1)->cid),
				getxyz(crd,bs->nid),getxyz(crd,bs->caid),
				getxyz(crd,bs->cid),&dphi);
		dphidx.push_back(std::make_pair((bs-1)->cid,dphi[0]));
		dphidx.push_back(std::make_pair(bs->nid,dphi[1]));
		dphidx.push_back(std::make_pair(bs->caid,dphi[2]));
		dphidx.push_back(std::make_pair(bs->cid,dphi[3]));
		phicodes=gettriangcodes(phi,dphidx,&dphicodesdx);
	} if(!Cterm){
		std::vector<XYZ> dpsi;
		psi=torsion(getxyz(crd,bs->nid),
				getxyz(crd,bs->caid),getxyz(crd,bs->cid),
				getxyz(crd,(bs+1)->nid),&dpsi);
		dpsidx.push_back(std::make_pair(bs->nid,dpsi[0]));
		dpsidx.push_back(std::make_pair(bs->caid,dpsi[1]));
		dpsidx.push_back(std::make_pair(bs->cid,dpsi[2]));
		dpsidx.push_back(std::make_pair((bs+1)->nid,dpsi[3]));
		psicodes=gettriangcodes(psi,dpsidx,&dpsicodesdx);
	}
}

std::vector<SSCode> NSPsd::estimatess(const std::vector<PhiPsiCodes> &phipsicodes){
	int nsite=phipsicodes.size();
	std::vector<SSCode> result(nsite);
#pragma omp parallel for schedule(dynamic,1)
	for(int i=0;i<nsite;++i) {
		auto & sc=result.at(i);
		sc.p3=NN_SSTerm().probabilities(phipsicodes,i,&(sc.dp3dx));
		sc.ssid=NN_SSTerm::sstype(sc.p3);
	}
	return result;
}
std::vector<std::vector<double>> NSPsd::estimatess(const std::vector<NSPproteinrep::BackBoneSite> & chain){
	std::vector<BSInChain> bsinchain=makeBSInChain(chain.size());
	std::vector<double> crd=NSPproteinrep::extractcrd(chain);
	std::vector<PhiPsiCodes> phipsicodes=makephipsicodes(crd,bsinchain);
	std::vector<SSCode> sscodes=estimatess(phipsicodes);
	std::vector<std::vector<double>> res;
	for(auto &ssc:sscodes) {
		res.push_back(ssc.p3);
	}
	return res;
}
