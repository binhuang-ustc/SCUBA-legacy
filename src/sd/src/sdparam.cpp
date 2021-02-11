/*
 * sdparam.cpp
 *
 *  Created on: 2019/2/15
 *      Author: hbfrank
 */

#include "sd/sdparam.h"

using namespace NSPsd;
using namespace NSPproteinrep;

SDInputOutput::SDInputOutput(const std::string &controlname)
    :controlname_(controlname) {
    auto & pset = SDInputOutputControls::getparameterset(controlname);
    pset.getval("SDSeed", &sdseed_);
    pset.getval("SDTotalSteps", &sdtotalsteps_);
    pset.getval("EchoStartPDBFile", &echostartpdbfile_);
    pset.getval("OutputPDBFile", &outputpdbfile_);
    pset.getval("OutputCoordFile", &outputcrdfile_);
}

void SDInputOutput::readiocontrolset(const std::string controlfile,
                                            const std::string controlname) {
    NSPdataio::ControlFile cf;
    cf.readfile(controlfile);
    std::vector<std::string> iocontrolines = cf.getcontrolines("SDInputOutput");
    SDInputOutput::defineiocontrol(controlname, iocontrolines);
}

void SDInputOutput::defineiocontrol(std::string name,
                                    const std::vector<std::string> &controllines) {
    std::map<std::string,std::string> stringpars{
        {"InputCoordFile",""}, {"EchoStartPDBFile",""},
        {"OutputPDBFile",""}, {"OutputCoordFile",""}};
    std::map<std::string,int> intpars{{"SDSeed",0}, {"SDTotalSteps",0}};
    std::map<std::string,std::vector<int>> intvecpars;
    std::map<std::string,double> doublepars;
    std::map<std::string,std::vector<double>> doublevecpars;
    std::map<std::string,std::vector<std::string>> stringvecpars;
    SDInputOutputControls::initdefaultkeyvals(name,doublepars,stringpars,intpars,doublevecpars,
            stringvecpars,intvecpars);
    int nsuccess=SDInputOutputControls::adjustvalues(name,controllines);
    if(nsuccess!= controllines.size()) {
        std::cout << "Execution aborted due to SDInputOutput control lines read error." << std::endl;
        exit(1);
    }
    auto & pset=SDInputOutputControls::getparameterset(name);
}

void SDInputOutput::printiocontrol(std::string name,std::ostream &ofs) {
    ofs<<"START SDInputOutput"<<std::endl;
    SDInputOutputControls::getparameterset(name+"_sdinputoutput").printparameters(ofs);
    ofs<<"END SDInputOutput"<<std::endl;
}

void SDInputOutput::writecrd(const std::vector<double> &crd,std::ostream & os){
    os<<crd.size();
    for(auto &c:crd) os<<" "<<c;
}

bool SDInputOutput::readcrd(std::istream &is,std::vector<double> &crd){
    int nrecord;
    is>>nrecord;
    crd.resize(nrecord);
    for(int i=0;i<nrecord;++i) is >>crd[i];
    return(!is.fail());
}

std::vector<double> SDInputOutput::extractcrd(
    const std::vector<std::vector<NSPproteinrep::FullSite>> &sites) {
    std::vector<NSPgeometry::XYZ> res;
    double deg=3.14159265/180.0;
    for (int c=0;c<sites.size();++c){
        for(int s=0;s<sites[c].size();++s){
            const FullSite &fs=sites[c][s];
            const VSCType &vsc=VSCType::getVSCType(fs.resname());
            std::vector<NSPgeometry::XYZ> r(vsc.nscatoms+4);
            r[0]=fs.getcrd("N");
            r[1]=fs.getcrd("CA");
            r[vsc.nscatoms+2]=fs.getcrd("C");
            r[vsc.nscatoms+3]=fs.getcrd("O");
            bool sccomplete{true};
            for(int a=0;a<vsc.nscatoms;++a){
                if(!fs.hasatomcrd(vsc.atomnames[a])) {
                    sccomplete=false;
                    break;
                }
            }
            if(sccomplete){
                for(int a=0;a<vsc.nscatoms;++a){
                    r[a+2]=fs.getcrd(vsc.atomnames[a]);
                }
            } else {
                auto & ics=vsc.internalcrds;
                assert(ics.size()==vsc.nscatoms);
                int aidx=2;
                for(auto &ic:ics){
                    r[aidx++]=NSPgeometry::InternaltoXYZ(r[ic[0].first],r[ic[1].first],
                                r[ic[2].first],
                                ic[0].second, ic[1].second*deg,ic[2].second*deg);
                        }
            }
            for(auto &x:r) res.push_back(x);
        }
    }
    std::vector<double> crd;
    for (auto &r : res) {
        for (int m = 0; m < 3; ++m)
            crd.push_back(r[m]);
    }
    return crd;
}


LoopParam::LoopParam(const std::string & controlname)
    :controlname_(controlname) {
    auto & pset = LoopParamControls::getparameterset(controlname);
    pset.getval("NumLoops", &numloops_);
    pset.getval("LoopLocation", &looplocation_);
    pset.getval("MaxLinkTryTimes", &maxlinktrytimes_);
    pset.getval("NumRunsThreshold", &numrunsthreshold_);
    pset.getval("FlucEneCutoff", &flucenecutoff_);
    pset.getval("MaxNumRuns", &maxnumruns_);
    pset.getval("NumRunsEneCheckPoint", &numrunsenecheckpoint_);
    pset.getval("TotEnergyCheckPoint", &totenergycheckpoint_);
    pset.getval("LowEnergyCutoff", &lowenergycutoff_);
    pset.getval("RedundantRMSDCutoff", &redundantrmsdcutoff_);
    pset.getval("MaxSampleTimes", &maxsampletimes_);
    pset.getval("FilterInterval", &filterinterval_);
    pset.getval("ConvergentCutoff", &convergentcutoff_);
    pset.getval("AllSampledLoopsFile", &allsampledloopsfile_);
    pset.getval("NonredundantLoopsFile", &nonredundantloopsfile_);
}

void LoopParam::readcontrolset(const std::string & controlfile,
                               const std::string & controlname) {
    NSPdataio::ControlFile cf;
    cf.readfile(controlfile);
    std::vector<std::string> iocontrolines = cf.getcontrolines("LoopParam");
    LoopParam::definecontrol(controlname, iocontrolines);
}

void LoopParam::definecontrol(const std::string & name,
                              const std::vector<std::string> & controllines) {
    std::map<std::string,std::string> stringpars {
        {"AllSampledLoopsFile", ""}, {"NonredundantLoopsFile", ""}
    };
    std::map<std::string,int> intpars{
        {"NumLoops", 0}, {"MaxLinkTryTimes", 0}, {"NumRunsThreshold", 0},
        {"MaxNumRuns", 0}, {"NumRunsEneCheckPoint", 0}, {"MaxSampleTimes", 0},
        {"FilterInterval", 1}, {"ConvergentCutoff", 1}
    };
    std::map<std::string,std::vector<int>> intvecpars{
        {"LoopLocation", {-1,-1,-1,-1}}
    };
    std::map<std::string,double> doublepars {
        {"FlucEneCutoff", 1000.0}, {"TotEnergyCheckPoint", 1000.0},
        {"LowEnergyCutoff", 0.0}, {"RedundantRMSDCutoff", 1.0},
    };
    std::map<std::string,std::vector<double>> doublevecpars;
    std::map<std::string,std::vector<std::string>> stringvecpars;
    LoopParamControls::initdefaultkeyvals(name,doublepars,stringpars,intpars,doublevecpars,
            stringvecpars,intvecpars);
    int nsuccess=LoopParamControls::adjustvalues(name,controllines);
    if(nsuccess!= controllines.size()) {
        std::cout << "Execution aborted due to LoopParam control lines read error." << std::endl;
        exit(1);
    }
    auto & pset=LoopParamControls::getparameterset(name);
}

void LoopParam::printcontrol(const std::string & name, std::ostream & ofs) {
    ofs<<"START LoopParam"<<std::endl;
    LoopParamControls::getparameterset(name+"_loopparam").printparameters(ofs);
    ofs<<"END LoopParam"<<std::endl;
}

LoopOptimizeParam::LoopOptimizeParam(const std::string & controlname)
    :controlname_(controlname) {
    auto & pset = LoopOptimizeParamControls::getparameterset(controlname);
    pset.getval("LoopLocations", &looplocations_);
    pset.getval("MaxLinkTryTimes", &maxlinktrytimes_);
    pset.getval("NumRunsThreshold", &numrunsthreshold_);
    pset.getval("FlucEneCutoff", &flucenecutoff_);
    pset.getval("MaxNumRuns", &maxnumruns_);
    pset.getval("NumRunsEneCheckPoint", &numrunsenecheckpoint_);
    pset.getval("TotEnergyCheckPoint", &totenergycheckpoint_);
    pset.getval("AllKeptLoopsFile", &allkeptloopsfile_);
    pset.getval("LowestEnergyPDBFile", &lowestEnergyPDBFile_);
    pset.getval("NoneViolatedLoopsFile", &noneViolatedLoopsFile_);
    pset.getval("MaxReconstructLoopsOnce", &maxReconstructLoopsOnce_);
    pset.getval("McNumThreads", &mcnumthreads_);
    pset.getval("McMaxNumSteps", &mcmaxnumsteps_);
    pset.getval("McContDiscardStepsCutoff", &mccontdiscardstepscutoff_);
    pset.getval("McTemperature", &mctemperature_);
}

void LoopOptimizeParam::readcontrolset(const std::string & controlfile,
                               const std::string & controlname) {
    NSPdataio::ControlFile cf;
    cf.readfile(controlfile);
    std::vector<std::string> iocontrolines = cf.getcontrolines("LoopOptimizeParam");
    LoopOptimizeParam::definecontrol(controlname, iocontrolines);
}

void LoopOptimizeParam::definecontrol(const std::string & name,
                              const std::vector<std::string> & controllines) {
    std::map<std::string,std::string> stringpars {
        {"AllKeptLoopsFile", ""}, {"LowestEnergyPDBFile", ""},
        {"NoneViolatedLoopsFile", ""}
    };
    std::map<std::string,int> intpars{
        {"MaxLinkTryTimes", 0}, {"NumRunsThreshold", 0},
        {"MaxNumRuns", 0}, {"NumRunsEneCheckPoint", 0},
        {"McNumThreads", 1}, {"McMaxNumSteps", 1000},
        {"McContDiscardStepsCutoff", 50},
        {"MaxReconstructLoopsOnce", 1}
    };
    std::map<std::string,std::vector<int>> intvecpars{
        {"LoopLocations", {-1,-1,-1,-1}}
    };
    std::map<std::string,double> doublepars {
        {"FlucEneCutoff", 1000.0}, {"McTemperature", 1.0},
        {"TotEnergyCheckPoint", 1000.0}
    };
    std::map<std::string,std::vector<double>> doublevecpars;
    std::map<std::string,std::vector<std::string>> stringvecpars;
    LoopOptimizeParamControls::initdefaultkeyvals(name,doublepars,stringpars,intpars,doublevecpars,
            stringvecpars,intvecpars);
    int nsuccess=LoopOptimizeParamControls::adjustvalues(name,controllines);
    if(nsuccess!= controllines.size()) {
        std::cout << "Execution aborted due to LoopOptimizeParam control lines read error." << std::endl;
        exit(1);
    }
    auto & pset=LoopOptimizeParamControls::getparameterset(name);
}

void LoopOptimizeParam::printcontrol(const std::string & name, std::ostream & ofs) {
    ofs<<"START LoopOptimizeParam"<<std::endl;
    LoopOptimizeParamControls::getparameterset(name+"_loopexploreparam").printparameters(ofs);
    ofs<<"END LoopOptimizeParam"<<std::endl;
}

LoopExploreParam::LoopExploreParam(const std::string & controlname)
    :controlname_(controlname) {
    auto & pset = LoopExploreParamControls::getparameterset(controlname);
    pset.getval("LoopDescriptors", &loopDescriptors_);
    pset.getval("MaxLinkTryTimes", &maxLinkTryTimes_);
    pset.getval("NumRunsThreshold", &numRunsThreshold_);
    pset.getval("EneFlucCutoff", &eneFlucCutoff_);
    pset.getval("MaxNumRuns", &maxNumRuns_);
    pset.getval("NumRunsEneCheckPoint", &numRunsEneCheckPoint_);
    pset.getval("TotEnergyCheckPoint", &totEnergyCheckPoint_);
    pset.getval("AllKeptLoopsFile", &allKeptLoopsFile_);
    pset.getval("NoneViolatedLoopsFile", &noneViolatedLoopsFile_);
    pset.getval("MaxExtendLen", &maxExtendLen_);
    pset.getval("MaxShiftLen", &maxShiftLen_);
    pset.getval("MinLoopLen", &minLoopLen_);
    pset.getval("MaxLoopLen", &maxLoopLen_);
    pset.getval("MaxNumOptmLoopsEach", &maxNumOptmLoopsEach_);
    pset.getval("MaxSDTimesOfOneLoop", &maxSDTimesOfOneLoop_);
}

void LoopExploreParam::readcontrolset(const std::string & controlfile,
                               const std::string & controlname) {
    NSPdataio::ControlFile cf;
    cf.readfile(controlfile);
    std::vector<std::string> iocontrolines = cf.getcontrolines("LoopExploreParam");
    LoopExploreParam::definecontrol(controlname, iocontrolines);
}

void LoopExploreParam::definecontrol(const std::string & name,
                              const std::vector<std::string> & controllines) {
    std::map<std::string,std::string> stringpars {
        {"AllKeptLoopsFile", ""},
        {"NoneViolatedLoopsFile", ""}
    };
    std::map<std::string,int> intpars{
        {"MaxLinkTryTimes", 0}, {"NumRunsThreshold", 0},
        {"MaxNumRuns", 0}, {"NumRunsEneCheckPoint", 0},
        {"MaxNumOptmLoopsEach", 1}, {"MaxSDTimesOfOneLoop", 1},
        {"MaxExtendLen", 0}, {"MaxShiftLen", 1},
        {"MinLoopLen", -1}, {"MaxLoopLen", -1}
    };
    std::map<std::string,std::vector<int>> intvecpars{
        {"LoopDescriptors", {-1,-1,-1,-1}}
    };
    std::map<std::string,double> doublepars {
        {"EneFlucCutoff", 1000.0},
        {"TotEnergyCheckPoint", 1000.0}
    };
    std::map<std::string,std::vector<double>> doublevecpars;
    std::map<std::string,std::vector<std::string>> stringvecpars;
    LoopExploreParamControls::initdefaultkeyvals(name,doublepars,stringpars,intpars,doublevecpars,
            stringvecpars,intvecpars);
    int nsuccess=LoopExploreParamControls::adjustvalues(name,controllines);
    if(nsuccess!= controllines.size()) {
        std::cout << "Execution aborted due to LoopExploreParam control lines read error." << std::endl;
        exit(1);
    }
    auto & pset=LoopExploreParamControls::getparameterset(name);
}

void LoopExploreParam::printcontrol(const std::string & name, std::ostream & ofs) {
    ofs<<"START LoopExploreParam"<<std::endl;
    LoopExploreParamControls::getparameterset(name+"_loopexploreparam").printparameters(ofs);
    ofs<<"END LoopExploreParam"<<std::endl;
}
