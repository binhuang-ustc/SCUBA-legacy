/*
 * restraints.cpp
 *
 *  Created on: 2017年11月10日
 *      Author: hyliu
 */

#include "sd/restraints.h"
#include "sd/bsinchain.h"
#include "geometry/quatfit.h"
#include "sd/forcefield.h"
using namespace NSPsd;

StructRestraint::StructRestraint(const std::vector<double> &crd,
		const std::vector<double> &weights, double kres,int mode,double rmsdref) :
		weights_(weights), kres_(kres),mode_(mode),rmsdref_(rmsdref) {
	assert(crd.size() / 3 == weights_.size());
	wtot_=0.0;
	for (int i = 0; i < crd.size(); i += 3) {
		crdref_.push_back(
				NSPgeometry::XYZ(crd.at(i), crd.at(i + 1), crd.at(i + 2)));
		wtot_+=weights[i/3];
	}
}
StructRestraint::StructRestraint(const std::vector<double> &crd,
		const std::vector<std::pair<int,int>> &poslengths, double kres,int mode,double rmsdref) :
		 kres_(kres),mode_(mode),rmsdref_(rmsdref) {
	weights_.assign(crd.size()/3,0.0);
	wtot_=0.0;
	for(auto iter=poslengths.begin(); iter != poslengths.end(); ++iter){
		for(int i=iter->first; i<iter->first+iter->second;++i){
			weights_[i]=1.0;
			wtot_+=1.0;
		}
	}
	for (int i = 0; i < crd.size(); i += 3) {
		crdref_.push_back(
				NSPgeometry::XYZ(crd.at(i), crd.at(i + 1), crd.at(i + 2)));
	}
}
double StructRestraint::posi_energy(const std::vector<NSPgeometry::XYZ> &crd,std::vector<NSPgeometry::XYZ> *forces) const{
	NSPgeometry::QuatFit qfit;
	std::vector<NSPgeometry::XYZ> crdtmp=crdref_;
	double rmsd2=qfit.fitting(crd,crdtmp,weights_);//changing the refcrd!!!
	const double kresh=0.5*kres_;
	double etot=0.0;
//	double wtot=0.0;
//	double d2=0.0;
	for(int i=0;i<crd.size();++i) {
		if(weights_[i]<1.e-10) continue;
		std::vector<NSPgeometry::XYZ> deriv;
		double b = NSPgeometry::distance(crd[i], crdtmp[i], &deriv);
		if(b<1.e-10) continue;
//		d2+=b*b;
//		wtot+=weights_[i];
		etot += weights_[i]*kresh * b*b;
		double kbdb = -weights_[i]*kres_ * b;
//		if(kbdb >1.e-3 || kbdb<-1.e-3){
//			std::cout <<kbdb<<" "<<kres_<<" "<<b<<" "<<weights_[i]<<std::endl;
//			std::cout<<i <<" " <<crdref_.size()<<std::endl;
//			std::cout<<crd.at(i).toString() <<crdtmp.at(i).toString()<<std::endl;
		forces->at(i) = forces->at(i) + kbdb * deriv[0];
	}
//	double  rmsd2_recalc=d2/wtot;
//	std::cout <<"rmsd2: "<<rmsd2_<< " "<<rmsd2_recalc<<std::endl;
	rmsd2_=rmsd2;
	return etot;
}
double StructRestraint::totalrmsd_energy(const std::vector<NSPgeometry::XYZ> &crd,std::vector<NSPgeometry::XYZ> *forces) const{
	NSPgeometry::QuatFit qfit;
	std::vector<NSPgeometry::XYZ> crdtmp=crdref_;
	double rmsd2=qfit.fitting(crd,crdtmp,weights_);//changing the refcrd!!!
	if(rmsd2<=rmsdref_*rmsdref_) {
		rmsd2_=rmsd2;
		return 0.0;
	}
	double rmsd=sqrt(rmsd2);
	double drmsd=rmsd-rmsdref_;
	double dedrmsd=kres_*drmsd;
	double ene=0.5*dedrmsd*drmsd;
	dedrmsd /=rmsd;
	for(int i=0;i<crd.size();++i) {
		if(weights_[i]<1.e-10) continue;
		forces->at(i) = forces->at(i) - (dedrmsd *weights_[i]/wtot_)* (crd[i]-crdtmp[i]);
	}
//	double  rmsd2_recalc=d2/wtot;
//	std::cout <<"rmsd2: "<<rmsd2_<< " "<<rmsd2_recalc<<std::endl;
	rmsd2_=rmsd2;
	return ene;
}

double RgRestraint::energy(const std::vector<NSPgeometry::XYZ> &crd,
		std::vector<NSPgeometry::XYZ> *forces,const std::vector<double> &w) const {
	if(kres_>0) return energy1(crd,forces,w);
	else return energy2(crd,forces,w);
}
double RgRestraint::energy1(const std::vector<NSPgeometry::XYZ> &crd,
		std::vector<NSPgeometry::XYZ> *forces,const std::vector<double> &w) const{
		if(kres_==0.0) return 0.0;
		NSPgeometry::XYZ center=NSPgeometry::center(crd,w);
		double rtot2=0.0;
		std::vector<NSPgeometry::XYZ> drt2dx;
		std::vector<NSPgeometry::XYZ> deriv;
		int idx=0;
		double wtot=0.0;
		for(auto &c:crd){
			NSPgeometry::XYZ dc=c-center;
			double r2=dc.squarednorm();
			double wgt=1.0;
			if(!w.empty()){
				wgt=w[idx];
			}
			rtot2+=wgt*r2;
			drt2dx.push_back(2.0*wgt*dc);  //drtot2/dx
			wtot +=wgt;
			++idx;
		}
		double rg=sqrt(rtot2/wtot);
		if(rg<=rgbound_) return 0.0;
		double drg=rg-rgbound_;
		double ene=0.5*kres_*drg*drg;
		double dedrtot2=-kres_*drg/(2.0*rg*wtot); //kres*drg * (1/(2*rg)*1/crd.size()
		for(int i=0;i<crd.size();++i){
			(*forces)[i] = (*forces)[i] +dedrtot2*drt2dx[i];
		}
		return ene;
}

/*
 * ene=-kres_*log(rg^2)
 * -kres_: should be proportional to temperature and the number of degrees of freedom
 */
double RgRestraint::energy2(const std::vector<NSPgeometry::XYZ> &crd,
		std::vector<NSPgeometry::XYZ> *forces,const std::vector<double> & w) const{
		if(kres_==0.0) return 0.0;
		NSPgeometry::XYZ center=NSPgeometry::center(crd,w);
		double rtot2=0.0;
		std::vector<NSPgeometry::XYZ> drt2dx;
		std::vector<NSPgeometry::XYZ> deriv;
		int idx=0;
		double wtot=0.0;
		for(auto &c:crd){
			double wgt=1.0;
			if(!w.empty()){
				wgt=w[idx];
			}
			NSPgeometry::XYZ dc=c-center;
			double r2=dc.squarednorm();
			rtot2+=wgt*r2;
			wtot +=wgt;
			drt2dx.push_back(2.0*wgt*dc);  //drtot2/dx
			++idx;
		}
		double rg2=rtot2/wtot;
		if(rg2<=rgbound_*rgbound_) return 0.0;
		double ene=-kres_*log(rg2/(rgbound_*rgbound_));
		double dedrtot2=kres_/(rg2*wtot); //kres*drg * (1/(2*rg)*1/crd.size()
		for(int i=0;i<crd.size();++i){
			(*forces)[i] = (*forces)[i] +dedrtot2*drt2dx[i];
		}
		return ene;
}
double DisRestraint::energy_attract(const std::vector<NSPgeometry::XYZ> & crd, std::vector<NSPgeometry::XYZ> *forces) const{
	std::vector<NSPgeometry::XYZ> drdx;
	double r=NSPgeometry::distance(crd[a1],crd[a2],&drdx);
	if (r<=r0) return 0.0;
	double dr=r-r0;
	double ene,dedr;
	if(r<r1) {
		dedr=kres*dr;
		ene=0.5*dedr*dr;
	} else{
		dedr=kres*(r1-r0);
		ene=0.5*dedr*(r1-r0)+dedr*(r-r1);
	}
	(*forces)[a1]=(*forces)[a1] - dedr*(drdx[0]);
	(*forces)[a2]=(*forces)[a2] - dedr*(drdx[1]);
	return ene;
}
double DisRestraint::energy_repulsion(const std::vector<NSPgeometry::XYZ> & crd, std::vector<NSPgeometry::XYZ> *forces) const{
	std::vector<NSPgeometry::XYZ> drdx;
	double r=NSPgeometry::distance(crd[a1],crd[a2],&drdx);
	if (r>=-r0) return 0.0;
	double dr=-r0-r;
	double ene,dedr;
	if( dr>r1){
		dedr=kres*dr;
		ene=0.5*dedr*dr;
	}
	else {
		dedr=kres*(-r0-r1);
		ene=0.5*dedr*(-r0-r1)+dedr*(r1-r);
	}
	(*forces)[a1]=(*forces)[a1] + dedr*(drdx[0]);
	(*forces)[a2]=(*forces)[a2] + dedr*(drdx[1]);
	return ene;
}

std::vector<std::pair<int,int>> NSPsd::capairsinstrandpair(int length, int parallel,int s1s,
		int s2s, const std::vector<BSInChain> & bcs1,const std::vector<BSInChain> &bcs2){
		assert(parallel==1 || parallel == -1);
		std::vector<std::pair<int,int>> result;
		for(int i=0;i<length; i++){
			int ca1=bcs1[s1s+i].caid;
			int ca2=bcs2[s2s+parallel*i].caid;
			result.push_back(std::make_pair(ca1,ca2));
		}
		return result;
}
StructRestraint NSPsd::makermsdrestraint(const std::string &filename,double rmsdref,double kres){
	std::vector<std::vector<NSPproteinrep::BackBoneSite>> chains;
	NSPproteinrep::readbackbonefrompdb(filename,chains);
	assert(chains.size()==1);
	std::vector<double> crdref=NSPproteinrep::extractcrd(chains[0]);
	for(auto &c:crdref) c *=A2NM;
	std::vector<double> w(crdref.size()/3,1.0);
	return StructRestraint(crdref,w,kres,StructRestraint::TOTALMODE,rmsdref);
}

