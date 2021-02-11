/*
 * sdparam.h
 *
 *  Created on: 2019/2/15
 *      Author: hbfrank
 */

#ifndef SD_SDPARAM_H_
#define SD_SDPARAM_H_

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <memory>
#include "dataio/parameters.h"
#include "fullsite/fullsite.h"
#include "sd/genchain.h"

namespace NSPsd{
class SDInputOutput {
public:
    SDInputOutput(const std::string &controlname);
    typedef NSPdataio::TypedMultiInstanceControls<SDInputOutput> SDInputOutputControls;
    static void readiocontrolset(const std::string controlfile, const std::string controlname);
    static void defineiocontrol(std::string name, const std::vector<std::string> &controllines);
    static void printiocontrol(std::string name,std::ostream &ofs);
    int sdseed() const { return sdseed_; }
    int sdtotalsteps() const { return sdtotalsteps_; }
    std::string echostartpdbfile() const { return echostartpdbfile_; }
    std::string outputpdbfile() const { return outputpdbfile_; }
    std::string outputcrdfile() const { return outputcrdfile_; }
    static void writecrd(const std::vector<double> &crd,std::ostream & os);
    static bool readcrd(std::istream &is,std::vector<double> &crd);
    static std::vector<double> extractcrd(const std::vector<std::vector<NSPproteinrep::FullSite>> &sites);
private:
    std::string controlname_;
    int sdseed_;
    int sdtotalsteps_;
    std::string echostartpdbfile_;
    std::string outputpdbfile_;
    std::string outputcrdfile_;
};

class LoopParam {
public :
    LoopParam(const std::string & controlname);
    typedef NSPdataio::TypedMultiInstanceControls<LoopParam> LoopParamControls;
    static void readcontrolset(const std::string & controlfile, const std::string & controlname);
    static void definecontrol(const std::string & name, const std::vector<std::string> & controllines);
    static void printcontrol(const std::string & name, std::ostream &ofs);
    int numloops() const { return numloops_; }
    std::vector<int> loopLocation() const { return looplocation_; }
    int maxLinkTryTimes() const { return maxlinktrytimes_; }
    int numRunsThreshold() const { return numrunsthreshold_; }
    double flucEneCutoff() const { return flucenecutoff_; }
    double maxNumRuns() const { return maxnumruns_; }
    int numRunsEneCheckPoint() const { return numrunsenecheckpoint_; }
    double totEnergyCheckPoint() const { return totenergycheckpoint_; }
    double lowEnergyCutoff() const { return lowenergycutoff_; }
    double redundantRMSDCutoff() const { return redundantrmsdcutoff_; }
    int maxSampleTimes() const { return maxsampletimes_; }
    int filterInterval() const { return filterinterval_; }
    int convergentCutoff() const { return convergentcutoff_; }
    std::string allSampledLoopsFile() const { return allsampledloopsfile_; }
    std::string nonredundantLoopsFile() const { return nonredundantloopsfile_; }
private:
    std::string controlname_;
    int numloops_; // max number of valid (low energy) loops
    std::vector<int> looplocation_;
    int maxlinktrytimes_;
    int numrunsthreshold_;
    double flucenecutoff_;
    int maxnumruns_;
    int numrunsenecheckpoint_;
    double totenergycheckpoint_;
    double lowenergycutoff_; // low energy cutoff
    double redundantrmsdcutoff_; // RMSD cutoff for definition of two redundant loops
    int maxsampletimes_;    // max total sample times (including failed ones)
    int filterinterval_;    // interval between two filtering
    int convergentcutoff_; // number of filter intervals that filtered loops count stays unchanged
    std::string allsampledloopsfile_;
    std::string nonredundantloopsfile_;
};

class LoopOptimizeParam {
public:
    LoopOptimizeParam(const std::string & controlname);
    typedef NSPdataio::TypedMultiInstanceControls<LoopOptimizeParam> LoopOptimizeParamControls;
    static void readcontrolset(const std::string & controlfile, const std::string & controlname);
    static void definecontrol(const std::string & name, const std::vector<std::string> & controllines);
    static void printcontrol(const std::string & name, std::ostream &ofs);
    std::vector<int> loopLocations() const { return looplocations_; }
    int maxLinkTryTimes() const { return maxlinktrytimes_; }
    int maxReconstructLoopsOnce() const { return maxReconstructLoopsOnce_; };
    double maxNumRuns() const { return maxnumruns_; }
    int numRunsThreshold() const { return numrunsthreshold_; }
    double flucEneCutoff() const { return flucenecutoff_; }
    int numRunsEneCheckPoint() const { return numrunsenecheckpoint_; }
    double totEnergyCheckPoint() const { return totenergycheckpoint_; }
    int mcNumThreads() const { return mcnumthreads_; }
    int mcMaxNumSteps() const { return mcmaxnumsteps_; }
    int mcContDiscardStepsCutoff() const { return mccontdiscardstepscutoff_; }
    double mcTemperature() const { return mctemperature_; }
    std::string allKeptLoopsFile() const { return allkeptloopsfile_; }
    std::string lowestEnergyPDBFile() const { return lowestEnergyPDBFile_; }
    std::string noneViolatedLoopsFile() const { return noneViolatedLoopsFile_; }
private:
    std::string controlname_;
    std::vector<int> looplocations_;
    int maxlinktrytimes_;
    int maxReconstructLoopsOnce_;
    int numrunsthreshold_;
    double flucenecutoff_;
    int maxnumruns_;
    int numrunsenecheckpoint_;
    double totenergycheckpoint_;
    int mcnumthreads_;
    int mcmaxnumsteps_;
    int mccontdiscardstepscutoff_;
    double mctemperature_;
    std::string allkeptloopsfile_;
    std::string lowestEnergyPDBFile_;
    std::string noneViolatedLoopsFile_;
};

class LoopExploreParam {
public:
    LoopExploreParam(const std::string & controlname);
    typedef NSPdataio::TypedMultiInstanceControls<LoopExploreParam> LoopExploreParamControls;
    static void readcontrolset(const std::string & controlfile, const std::string & controlname);
    static void definecontrol(const std::string & name, const std::vector<std::string> & controllines);
    static void printcontrol(const std::string & name, std::ostream &ofs);
    std::vector<int> loopDescriptors() const { return loopDescriptors_; }
    int maxLinkTryTimes() const { return maxLinkTryTimes_; }
    int maxExtendLen() const { return maxExtendLen_; }
    int maxShiftLen() const { return maxShiftLen_; }
    int minLoopLen() const { return minLoopLen_; }
    int maxLoopLen() const { return maxLoopLen_; }
    double maxNumRuns() const { return maxNumRuns_; }
    int numRunsThreshold() const { return numRunsThreshold_; }
    double eneFlucCutoff() const { return eneFlucCutoff_; }
    int numRunsEneCheckPoint() const { return numRunsEneCheckPoint_; }
    double totEnergyCheckPoint() const { return totEnergyCheckPoint_; }
    int maxNumOptmLoopsEach() const { return maxNumOptmLoopsEach_; }
    int maxSDTimesOfOneLoop() const { return maxSDTimesOfOneLoop_; }
    std::string allKeptLoopsFile() const { return allKeptLoopsFile_; }
    std::string noneViolatedLoopsFile() const { return noneViolatedLoopsFile_; }
private:
    std::string controlname_;
    std::vector<int> loopDescriptors_;
    int maxLinkTryTimes_;
    int maxExtendLen_;
    int maxShiftLen_;
    int minLoopLen_;
    int maxLoopLen_;
    int numRunsThreshold_;
    double eneFlucCutoff_;
    int maxNumRuns_;
    int numRunsEneCheckPoint_;
    double totEnergyCheckPoint_;
    int maxNumOptmLoopsEach_; // max number of optimized loops for each loop descriptor
    int maxSDTimesOfOneLoop_;
    std::string allKeptLoopsFile_;
    std::string noneViolatedLoopsFile_;
};
}
#endif /* SD_SDPARAM_H_ */
