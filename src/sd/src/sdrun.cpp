/*
 * sdrun.cpp
 *
 *  Created on: 2017年11月9日
 *      Author: hyliu
 */

#include "sd/sdrun.h"
using namespace NSPsd;
SDRun NSPsd::make_backbone_sdrun(const SDRun::SDRunIn &in,unsigned int seed) {
	int nsite = in.crd->size()/12;
	auto ff = std::shared_ptr < ForceField > (new ForceField);
	auto &pset=SDControls::getparameterset(in.sdcontrolname+"_sd");
	std::vector<int> nsites;
	nsites.push_back(nsite);
	std::vector<std::vector<int>> cissites;
	cissites.push_back(in.cissites);
	*ff = make_forcefield_backbone(nsites,cissites);
	ff->usecontrols(in.sdcontrolname+"_ff");
	auto shakebds = std::shared_ptr < ShakeBonds > (new ShakeBonds);
	int doshake;
	pset.getval("DoShake",&doshake);
	*shakebds = make_shakebonds(*ff);
	shakebds->seton(doshake!=0);
	auto sd = std::shared_ptr < StochasticDynamics > (new StochasticDynamics);
	*sd = make_sd_backbone(nsite,in.sdcontrolname+"_sd");
	SDRun sdrun(sd, ff, shakebds);
	sdrun.shakeon()=doshake!=0;
	int nblsteps;
	pset.getval("NeighborListSteps",&nblsteps);
	sdrun.nblsteps()=nblsteps;
	sdrun.initrandomengine(seed);
	sdrun.initstate(in);
	return sdrun;
}
void NSPsd::sdreadcontrols(const std::string &filename,std::string name){
	NSPdataio::ControlFile cf;
	cf.readfile(filename);
	std::vector<std::string> sdcontrolines=cf.getcontrolines("SD");
	std::vector<std::string> ffcontrolines=cf.getcontrolines("ForceField");
	definesdcontrol(name+"_sd",sdcontrolines);
	defineforcefieldcontrol(name+"_ff",ffcontrolines);
}
void NSPsd::sdprintcontrols(std::string name,std::ostream &ofs){
	ofs <<"START SD"<<std::endl;
	SDControls::getparameterset(name+"_sd").printparameters(ofs);
	ofs<<"END SD"<<std::endl;
	ofs<<"START ForceField"<<std::endl;
	ForceFieldControls::getparameterset(name+"_ff").printparameters(ofs);
	ofs<<"END ForceField"<<std::endl;
}
