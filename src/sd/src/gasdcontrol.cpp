/*
 * gasdcontrol.cpp
 *
 *  Created on: 2018年1月5日
 *      Author: hyliu
 */

#include "sd/gasdcontrol.h"
#include "sd/forcefield.h"
#include "sd/sdrun.h"
using namespace NSPsd;
GASDPar::GASDPar(const std::string &controlname){
	auto &pset=GASDControls::getparameterset(controlname+"_gasd");
	pset.getval("RMSD_SepConfig",&rmsd_sep);
	rmsd_sep *=A2NM;
	pset.getval("Temperature_OPT",&temperature_opt);
	pset.getval("Temperature_Explore",&temperature_explore);
	pset.getval("MaxGeneration",&maxgenerations);
	pset.getval("MaxStoredConfig",&maxstored_config);
	pset.getval("NChildrenPerConfig",&nchildren);
	pset.getval("MaxSDSteps_Opt",&steps_opt);
	pset.getval("MaxSDSteps_Expl",&steps_explore);
	pset.getval("Population_Size",&popsize);
	pset.getval("SubPopSize",&subsize);
	pset.getval("ShrinkGenerations",&shrink_generations);
}

void NSPsd::definegasdcontrol(const std::string & name,
		const std::vector<std::string> &controllines){
	std::map<std::string,double> doublepars{{"RMSD_SepConfig",4.0},{"Temperature_OPT",0.3},
			{"Temperature_Explore",3.0}};
	std::map<std::string,std::vector<std::string>> stringvecpars;
	std::map<std::string,std::vector<double>> doublevecpars;
	std::map<std::string,std::string> stringpars;
	std::map<std::string,std::vector<int>> intvecpars;
	std::map<std::string,int>intpars{{"MaxGeneration",0},{"MaxStoredConfig",0},
		{"NChildrenPerConfig",1},{"MaxSDSteps_Opt",0},{"MaxSDSteps_Expl",0},
		{"Population_Size",1},{"SubPopSize",1},{"ShrinkGenerations",5}};
	GASDControls::initdefaultkeyvals(name,doublepars,stringpars,intpars,doublevecpars,
			stringvecpars,intvecpars);
	int nsuccess=GASDControls::adjustvalues(name,controllines);
	if(nsuccess!= controllines.size()) {
		exit(1);
	}
}
void NSPsd::gasdreadcontrols(const std::string &filename,std::string name){
	NSPdataio::ControlFile cf;
	cf.readfile(filename);
	std::vector<std::string> sdcontrolines=cf.getcontrolines("SD");
	std::vector<std::string> ffcontrolines=cf.getcontrolines("ForceField");
	std::vector<std::string> gasdcontrolines=cf.getcontrolines("GASD");
	definesdcontrol(name+"_sd",sdcontrolines);
	defineforcefieldcontrol(name+"_ff",ffcontrolines);
	definegasdcontrol(name+"_gasd",gasdcontrolines);
}
void NSPsd::gasdprintcontrols(std::string name,std::ostream &ofs){
	ofs<<"START GASD" <<std::endl;
	GASDControls::getparameterset(name+"_gasd").printparameters(ofs);
	ofs <<"END GASD" <<std::endl;
	ofs <<"START SD"<<std::endl;
	SDControls::getparameterset(name+"_sd").printparameters(ofs);
	ofs<<"END SD"<<std::endl;
	ofs<<"START ForceField"<<std::endl;
	ForceFieldControls::getparameterset(name+"_ff").printparameters(ofs);
	ofs<<"END ForceField"<<std::endl;
}

void sdprintcontrols(std::string name,std::ostream &ofs=std::cout);
