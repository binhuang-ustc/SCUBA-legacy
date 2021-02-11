/*
 * ssrestraint.cpp
 *
 *  Created on: 2017年12月20日
 *      Author: hyliu
 */

#include "sd/ssrestraint.h"
#include "sd/nnterm.h"
using namespace NSPsd;
using namespace NSPgeometry;

double SSRestraint::energy(const std::vector<std::vector<SSCode>> &sscodes,
		std::vector<XYZ> *forces) const {
	auto & codes=sscodes.at(chainid_);
	double psum=0;
	double nsites=(double)(siteend_-sitebegin_+1);
	for(auto it=codes.begin()+sitebegin_; it!=codes.begin()+siteend_+1; ++it){
		if(it->ssid != NN_SSTerm::TERMINUS)	psum+=it->p3[targetssid_];
	}
	psum /= nsites;
	if(psum >= targetssratio_) return 0.0;
	double dp=psum-targetssratio_;
	double ene=0.5*kres_*dp*dp;
	double dedp=-kres_*dp;
	if(dp<-0.2){
		ene=0.02*kres_-0.2*kres_*(dp+0.2);
		dedp=0.2*kres_;
	}
	dedp /= nsites;
	for(auto it=codes.begin()+sitebegin_; it!=codes.begin()+siteend_+1; ++it){
		if(it->ssid ==NN_SSTerm::TERMINUS) continue;
		auto & dpdxs=it->dp3dx[targetssid_];
		for(auto &d:dpdxs){
			(*forces)[d.first]= (*forces)[d.first] + dedp*d.second;
		}
	}
	return ene;
}

