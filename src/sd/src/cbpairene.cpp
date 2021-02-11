/*
 * cbpairene.cpp
 *
 *  Created on: 2018年1月24日
 *      Author: hyliu
 */
#include "sd/cbpairene.h"
#include "dataio/datapaths.h"
#include <cassert>
#include <iostream>
#include <fstream>
using namespace NSPsd;

double CBPairEne::cbforce(double rcb,double *dedrcb){
	static CBPairEne cbpairene;
	if(!cbpairene.tableread_){
#ifdef _OPENMP
#pragma omp critical(cbpairene_global)
		{
			if(!cbpairene.tableread_) {
#endif
//		std::string filename=NSPdataio::datapath()+"cbpairene.dat";
				std::string filename=NSPdataio::datafilename("cbpairene.dat");
				cbpairene.readetable(filename);
#ifdef _OPENMP
#pragma omp flush
			}
		}
#endif
		}
	return cbpairene.force(rcb,dedrcb);
}
void CBPairEne::readetable(const std::string &filename){
	etable_.resize(nsteps_+1);
	std::ifstream ifs;
	ifs.open(filename.c_str());
	double r;
	for(int i=0;i<nsteps_+1;++i){
		ifs >>r >>etable_[i];
	}
	assert(etable_[nsteps_]==0.0);
	tableread_=true;
}
double CBPairEne::force(double rcb,double *dedrcb) const{
	double hstep=0.5*rstep_;
	if(rcb<=hstep) {
		*dedrcb=(etable_[1]-etable_[0])/rstep_;
		return etable_[0];
	}
	int i1=(rcb-hstep)/rstep_;
	if(i1>=nsteps_) {
		*dedrcb=0.0;
		return 0.0;
	}
	double r1=rstep_*i1+hstep;
	double alpha=(rcb-r1)/rstep_;
	double e=etable_[i1]*(1-alpha)+etable_[i1+1]*alpha;
	*dedrcb=(etable_[i1+1]-etable_[i1])/rstep_;
	return e;
}

