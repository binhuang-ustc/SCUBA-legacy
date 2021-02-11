/*
 * testchainsd.cpp
 *
 *  Created on: 2018年1月25日
 *      Author: hyliu
 */

#include "sd/sdparam.h"
#include "sd/genchain.h"
#include "sd/sidechainff.h"
#include "fullsite/fullsite.h"
#include "fullsite/structplus.h"
#include "sd/sdrun.h"
#include "geometry/quatfit.h"
#include "sd/trajio.h"
#include "dstl/randomengine.h"
#include <iostream>
#include <fstream>
#include <ctime>

using namespace NSPproteinrep;
using namespace NSPsd;

struct TemperatureAnnealing{
	double t_high{0.0};
	double t_low{0.0};
	int nsteps_cycle{0},nsteps_high{0}, nsteps_drop{0};
	TemperatureAnnealing(){;}
	TemperatureAnnealing(double th,double tl, double cycle,double sh,double sd):
	t_high(th),t_low(tl),nsteps_cycle(cycle),nsteps_high(sh),nsteps_drop(sd){
		assert(sh+sd <= cycle);
	}
	double temperature(int step) const {
		step=(step-1)%nsteps_cycle;
		if(step<nsteps_high) return t_high;
		else if(step>=(nsteps_high+nsteps_drop)) return t_low;
		else return t_high - (t_high-t_low)*(double)(step-nsteps_high)/(double)(nsteps_drop);
	}
};

/**Called back in SDRun, this functor monitor quantities during SD.
 * It can also be used to effect conditional earlier ending of the
 * calling SDRun (by returning true)
 *
 */
struct SDStepCallBack {
	/**
	 * called in SDRun at every step, after the force calculation
	 */
	bool operator()(SDRun &sdrun) const {
		if(trajstep_>0 &&(sdrun.nstepsrun()+1)%trajstep_==0){
			trajio_.writeframe(sdrun.state().crd);
			trajio_.writeframe(sdrun.potenergies());
		}
		if(doannealing_){
			sdrun.changetemperature(ta_.temperature(sdrun.nstepsrun()+1),tagroup_);
		}
		//no earlier ending of the calling SDRun
		return false;
	}
	void setrefcrd(const std::vector<double> &ref){
		refcrd_=ref;
	}
	void setannealing(const std::string &controlname){
		auto &pset=SDControls::getparameterset(controlname);
		std::vector<double> annealingscheme;
		pset.getval("AnnealingScheme",&annealingscheme);
		if(annealingscheme.empty()) return;
		doannealing_=true;
		int widx=0;
		tagroup_=(int) (annealingscheme[widx++]);
		double th=annealingscheme[widx++];
		double tl=annealingscheme[widx++];
		int nstep_cycle= (int) annealingscheme[widx++];
		int step_h=(int) annealingscheme[widx++];
		int step_d=(int) annealingscheme[widx++];
		ta_=TemperatureAnnealing(th,tl,nstep_cycle,step_h,step_d);
	}
	void createtrajfile(const std::string &controlname){
		auto &pset=SDControls::getparameterset(controlname);
		std::string trajfile;
		pset.getval("TrajFile",&trajfile);
		pset.getval("NTrajSteps", &trajstep_);
		if(!trajfile.empty() && trajstep_>0) trajio_.createwriter(trajfile);
		else trajstep_=-1;
	}
	SDStepCallBack(const std::string & controlname){
		setannealing(controlname);
		createtrajfile(controlname);
	}
	std::vector<double> refcrd_;
	bool doannealing_{false};
	TemperatureAnnealing ta_;
	int tagroup_{0};
	TrajIO trajio_;
	int trajstep_{-1};
};

std::string getCurrentLocalTimeString() {
    time_t rawtime;
    struct tm * timeinfo;
    char tibuffer[80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(tibuffer,sizeof(tibuffer),"%Y-%m-%d %H:%M:%S", timeinfo);
    std::string timestring(tibuffer);
    return timestring;
}

std::string getVersion() {
    return "0.1";
}

void printUsage(const std::string & selfname)
{
    std::cout << "Usage:" << std::endl;
    std::cout << "$ " << selfname << " <ParameterFile>"
            << std::endl;
}

int main(int argc, char **argv) {
	std::string selfname(argv[0]);
	if (argc != 2) {
	    printUsage(selfname);
	    exit(1);
	}

	std::ostream & oslog = std::cout;

	oslog << "SCUBA SD version " << getVersion() << " started at "
	        << getCurrentLocalTimeString() << std::endl;

	std::string controlfile(argv[1]);
	std::string controlname="sdffcontrol";
	genchainreadcontrols(controlfile,controlname);
	std::string sdiocontrolname=controlname + "_sdinputoutput";
	SDInputOutput::readiocontrolset(controlfile, sdiocontrolname);

	oslog << "Echo of input control parameters:" << std::endl;
	oslog << std::endl;

	SDInputOutput::printiocontrol(controlname, oslog);
	genchainprintcontrols(controlname,oslog);
	oslog << std::endl;

	SDInputOutput sdio(sdiocontrolname);

	unsigned int seed=std::abs(sdio.sdseed());
	NSPdstl::RandomEngine<>::getinstance().reseed(seed);

	std::string gccontrolname = controlname + "_genchain";
    GenChain genchain(gccontrolname);
    std::shared_ptr<std::vector<double>> initcrd = genchain.initcrd();

	SDRun::SDRunIn sdrunin(*initcrd, controlname);
	SDRun sdrun = genchain.make_sdrun(sdrunin,seed);
	int totalsteps = sdio.sdtotalsteps();

	std::string sdcontrolname = controlname + "_sd";
	SDStepCallBack callback(sdcontrolname);
	std::shared_ptr<std::vector<double>> rmsdrefcrd = genchain.rmsdrefcrd();
    callback.setrefcrd(*rmsdrefcrd);

    std::string echostartpdbfile = sdio.echostartpdbfile();
    if (!echostartpdbfile.empty()) {
        // write actual starting configuration to a file in PDB format if configured
        std::ofstream ofsstart(echostartpdbfile);
        if (ofsstart.is_open()) {
            genchain.writepdb(*initcrd, ofsstart, 1.0/A2NM);
            ofsstart.close();
        }
    }

	std::string outputpdbfile = sdio.outputpdbfile();

	int stepinterval = 100;
	int numruns = totalsteps / stepinterval;
	std::vector<int> eneterms2print = {ForceField::ENECOMP::ETOT, ForceField::ENECOMP::EBOND,
	        ForceField::ENECOMP::EANG,ForceField::ENECOMP::EIMPDIH,ForceField::ENECOMP::EPHIPSI,
	        ForceField::ENECOMP::ESCCONF,ForceField::ENECOMP::ESTERIC,ForceField::ENECOMP::ESCPACKING,
	        ForceField::ENECOMP::ELOCALSTRUCTURE, ForceField::ENECOMP::ELOCALHB,ForceField::ENECOMP::ESITEPAIRS,
	        ForceField::ENECOMP::ESTRUCTREST,ForceField::ENECOMP::ERGRESTRAINT,ForceField::ENECOMP::ESSRESTRAINT
	};

	oslog << "Summarizing results will be printed every " << stepinterval
	            << " SD steps." << std::endl;
	oslog << "The order of terms in the printed potential energy lists:" << std::endl;
	oslog << "TOTAL, bond, angle, improper_dihedral, ramachandran, rotamer, "
	      << "steric, sidechain_packing, mainchain_local_correlation, "
	      << "mainchain_local_hydrogenbond, mainchain_nonlocal_packing, "
	      << "other_structure_restraint, radius_of_gyration_restraint, "
	      << "secondary_structure_restraint" << std::endl;
	oslog << std::endl;
    oslog << "Potential energies at starting point: " << std::endl;
    for(int idx : eneterms2print) {
        oslog << sdrun.potenergies()[idx] << "  ";
    }
    oslog << std::endl;
    oslog << std::endl;
	for(int i = 0; i < numruns; ++i) {
		if(!(sdrun.runsteps(stepinterval, callback))){
			oslog<<"Shake failure occurred."<<std::endl;
			exit(1);
		}

		double temp=sdrun.temperature();
		oslog << "run " << (i+1) << " : nstepsrun: " << sdrun.nstepsrun()
		            << " , temperature: " << temp << std::endl;
		oslog << "energies:  ";
		for(int idx : eneterms2print) {
		    oslog << "  " << sdrun.potenergies()[idx];
		}
		oslog << std::endl;

        // print Rg value to log
        double rg = NSPgeometry::radiusgyr(sdrun.state().crd)/A2NM;
        oslog << "radius of gyration: " << rg << " ";

        // print main-chain and all-atom RMSDs to log
        std::vector<double> w(rmsdrefcrd->size()/3, 0.0);
        std::vector<int> mcatoms=sdrun.ff()->mainchainatoms();
        for(auto i:mcatoms) w[i]=1.0;
        NSPgeometry::QuatFit qfmc,qfall;
        double rmsd2mc=qfmc.setup(sdrun.state().crd, *rmsdrefcrd, w);
        double rmsd2=qfall.setup(sdrun.state().crd, *rmsdrefcrd);
        oslog << ", RMSD from reference, main chain: " << sqrt(rmsd2mc)/A2NM
                << " all atom: " << sqrt(rmsd2)/A2NM << std::endl;
        oslog << std::endl;

		if (!outputpdbfile.empty()) {
		    // write the real-time state in PDB format if configured
            std::ofstream ofsoutpdb(outputpdbfile);
            if (ofsoutpdb.is_open()) {
                genchain.writepdb(sdrun.state().crd, ofsoutpdb, 1.0/A2NM);
                ofsoutpdb.close();
            }
		}
	}

	std::string outputcrdfile = sdio.outputcrdfile();
	if (!outputcrdfile.empty()) {
	    // write the final coordinates in raw text format if configured
        std::ofstream ofsoutcrd(outputcrdfile);
        if (ofsoutcrd.is_open()) {
            sdio.writecrd(sdrun.state().crd,ofsoutcrd);
            ofsoutcrd.close();
        }
	}

	oslog << "SCUBA SD version " << getVersion() << " successfully finished at "
	        << getCurrentLocalTimeString() << std::endl;
}
