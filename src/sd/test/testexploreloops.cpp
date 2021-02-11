/*
 * testexploreloops.cpp
 *
 *  Created on: 2020.6.26
 *      Author: hbfrank
 */

#include "sd/sdparam.h"
#include "sd/loopexplorationio.h"
#include "sd/genchain.h"
#include "sd/sidechainff.h"
#include "fullsite/fullsite.h"
#include "fullsite/structplus.h"
#include "sd/sdrun.h"
#include "geometry/quatfit.h"
#include "sd/trajio.h"
#include "dstl/randomengine.h"
#include "backbone/backbonebuilder.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <time.h>

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

typedef std::vector<FullSite> Peptide;
typedef std::vector<std::vector<FullSite>> Protein;
typedef std::vector<BackBoneSite> BackboneChain;
typedef std::pair<BackBoneSite, BackBoneSite> AnchorPair;



std::vector<BackBoneSite> buildLinkerBetween(
        const std::pair<BackBoneSite, BackBoneSite> & anchorSites,
        const int length, const int maxTryTimes)
{
    std::vector<std::shared_ptr<std::vector<BackBoneSite>>> linker;
    int trytimes = 0;
    while (trytimes++ < maxTryTimes) {
        linker = BackBoneBuilder::buildlinkers(
                                  length, anchorSites.first, anchorSites.second,
                                  std::vector<std::pair<int,int>>(),
                                  std::vector<std::pair<int,int>>(),
                                  std::set<int>());
        if (!linker.empty()) {
            break;
        }
    }
    if (linker.empty()) {
        return std::vector<BackBoneSite>();
    }

    return *linker.at(0);
}

AnchorPair extractAnchorSites(const Protein & chains, const LoopLocation & loopLoc)
{
    const Peptide & chain = chains.at(loopLoc.idxChain);
    FullSite fLeft = chain.at(loopLoc.idxNTerminal-1);
    FullSite fLeftNext = chain.at(loopLoc.idxNTerminal);
    FullSite fRight = chain.at(loopLoc.idxCTerminal+1);
    FullSite fRightPrev = chain.at(loopLoc.idxCTerminal);
    BackBoneSite bLeft = fLeft.getbackbonesite();
    bLeft.psi(fLeftNext.getbackbonesite());
    BackBoneSite bRight = fRight.getbackbonesite();
    bRight.phi(fRightPrev.getbackbonesite());
    AnchorPair anchors(bLeft, bRight);
    return anchors;
}

int constructLoops(const Protein & oriChains, Protein * newChains,
                    const std::vector<LoopDescriptor> & oriLoopDess,
                    std::vector<LoopLocation> * newLoopLocs,
                    const int maxTryTimes) {
    /// sort loop locations by order of chain_index->residue_index
    std::vector<LoopDescriptor> sortedLoopDess;
    for (int c = 0; c < oriChains.size(); c++) {
        const auto chain = oriChains.at(c);
        int s = 0;
        while (s < chain.size()) {
            for (auto des : oriLoopDess) {
                if (des.loopLoc.idxChain == c && des.loopLoc.idxNTerminal == s) {
                    sortedLoopDess.push_back(des);
                    s = des.loopLoc.idxCTerminal + 1;
                    continue;
                }
            }
            s++;
        }
    }
    /// get loop anchor sites
    std::vector<AnchorPair> anchorPairs;
    for (auto des : sortedLoopDess) {
        anchorPairs.push_back(extractAnchorSites(oriChains, des.loopLoc));
    }
    /// build loops defined by sortedLoopDess and generate newLoopDess
    std::vector<BackboneChain> loops;
    newLoopLocs->clear();
    int curChainIdx = 0;
    int curOffset = 0;
    for (int i = 0; i < anchorPairs.size(); i++) {
        auto & anchor = anchorPairs.at(i);
        auto & des = sortedLoopDess.at(i);
        BackboneChain loop;
        loop = buildLinkerBetween(anchor, des.newLen, maxTryTimes);
        if (loop.empty()) {
            return -1;
        }
        loops.push_back(loop);
        if (curChainIdx != des.loopLoc.idxChain) {
            curOffset = 0;
        }
        LoopLocation newLoc;
        newLoc.idxChain = des.loopLoc.idxChain;
        newLoc.idxNTerminal = des.loopLoc.idxNTerminal + curOffset;
        newLoc.idxCTerminal = des.loopLoc.idxNTerminal + des.newLen - 1;
        newLoopLocs->push_back(newLoc);
        curOffset += des.newLen - (des.loopLoc.idxCTerminal - des.loopLoc.idxNTerminal + 1);
    }

    /// assemble full chains
    newChains->clear();
    curChainIdx = 0;
    int curSiteIdx = 0;
    for (int i = 0; i < sortedLoopDess.size(); i++) {
        auto & des = sortedLoopDess.at(i);
        if (curChainIdx != des.loopLoc.idxChain) {
            for (int c = curChainIdx; c < des.loopLoc.idxChain; c++) {
                newChains->push_back(oriChains.at(c));
            }
            curChainIdx = des.loopLoc.idxChain;
            curSiteIdx = 0;
        }
        if (newChains->size() < curChainIdx + 1) {
            std::vector<FullSite> curChain;
            newChains->push_back(curChain);
        }
        for (int s = curSiteIdx; s < des.loopLoc.idxNTerminal; s++) {
            newChains->at(curChainIdx).push_back(oriChains.at(curChainIdx).at(s));
        }
        curSiteIdx = des.loopLoc.idxCTerminal + 1;
        for (auto & bs : loops.at(i)) {
            bs.resname = "GLY";
            newChains->at(curChainIdx).push_back(make_fullsite(bs));
        }
    }
    // add the rest residues of the current chain
    for (int s = curSiteIdx; s < oriChains.at(curChainIdx).size(); s++) {
        newChains->at(curChainIdx).push_back(oriChains.at(curChainIdx).at(s));
    }
    curChainIdx++;
    for (int c = curChainIdx; c < oriChains.size(); c++) {
        newChains->push_back(oriChains.at(c));
    }
    return 0;
}


int runSD(std::ostream & oslog, std::string & controlname,
          GenChain & genchain, LoopExploreParam & lpexplparam, SDInputOutput & sdio,
          const std::vector<int> & eneterms2print,
          std::vector<double> * sdendcrds, std::vector<double> * sdendenes,
          EneComp & ec) {
    std::shared_ptr<std::vector<double>> initcrd = genchain.initcrd();
    SDRun::SDRunIn sdrunin(*initcrd, controlname);
    int sdioseed = sdio.sdseed();
    int currsdseed = NSPdstl::RandomEngine<>::getinstance()
                                    .intrng(sdioseed, 3*sdioseed)();
    //upscale high energy terms
    ForceField::highescale=4.0;
    SDRun sdrun = genchain.make_sdrun(sdrunin, currsdseed);
    std::string sdcontrolname = controlname + "_sd";
    SDStepCallBack callback(sdcontrolname);
    std::shared_ptr<std::vector<double>> rmsdrefcrd = genchain.rmsdrefcrd();
    callback.setrefcrd(*rmsdrefcrd);

    oslog << "[Info][RunSD] Seed Number: " << currsdseed << std::endl;
    oslog << "[Info][RunSD] Potential energies at starting point: " << std::endl;
    for(int idx : eneterms2print) {
        oslog << sdrun.potenergies()[idx] << "  ";
    }
    oslog << std::endl;

    int PeriodRunsNumThreshold = lpexplparam.numRunsThreshold();
    double PeriodEneFlucCutoff = lpexplparam.eneFlucCutoff();
    int NumRunsEneCheckPoint = lpexplparam.numRunsEneCheckPoint();
    double TotEnergyCheckPoint = lpexplparam.totEnergyCheckPoint();
    const int MaxNumRuns = lpexplparam.maxNumRuns();
    int stepinterval = 100;
    int numruns = 0;
    std::vector<double> periodEnes = {};
    bool keeprun = true;
    while (keeprun) {
        if(!(sdrun.runsteps(stepinterval, callback))){
            oslog<<"[Info][RunSD] Shake failure occurred."<<std::endl;
            return 0;
        }
        double temp=sdrun.temperature();
        oslog << "run " << (numruns+1) << " : nstepsrun: " << sdrun.nstepsrun()
                    << " , temperature: " << temp << std::endl;
        oslog << "energies:  ";
        for(int idx : eneterms2print) {
            oslog << "  " << sdrun.potenergies()[idx];
        }
        oslog << std::endl;

        double currTotEne = sdrun.potenergies()[ForceField::ENECOMP::ETOT];
        periodEnes.insert(periodEnes.begin(), currTotEne);
        if (periodEnes.size() > PeriodRunsNumThreshold) {
            periodEnes.pop_back();
        }
        if (periodEnes.size() == PeriodRunsNumThreshold) {
            double minEne = DBL_MAX;
            double maxEne = -DBL_MAX; // NOTE: NOT DBL_MIN, which is a very small POSITIVE value
            for (double ene : periodEnes) {
                if (ene > maxEne) maxEne = ene;
                if (ene < minEne) minEne = ene;
            }
            double flucTotEne = maxEne - minEne;
            if (flucTotEne <= PeriodEneFlucCutoff) {
                keeprun = false;
            }
        }

        numruns++;
        if (numruns > NumRunsEneCheckPoint && currTotEne > TotEnergyCheckPoint) {
            oslog << "[Info][RunSD] total energy exceeded the limit at check point, abort." << std::endl;
            return 0;
        }
        if (numruns > MaxNumRuns) {
            oslog << "[Info][RunSD] reached MaxNumRuns limit, abort." << std::endl;
            return 0;
        }
    }
    sdendcrds->clear();
    for (double d : sdrun.state().crd) {
        sdendcrds->push_back(d/A2NM);
    }
    sdendenes->clear();
    for (int i : eneterms2print) {
        sdendenes->push_back(sdrun.potenergies()[i]);
    }

    ec = sdrun.eneanalyze();

    return sdrun.nstepsrun();
}


bool filter(const EneComp & ec, int & nppViolation, int & nlsViolation,
            const std::vector<LoopLocation> & loopLocs) {
    std::set<int> tmpViolatedLoopsSet;
    auto & ephipsi = ec.ephipsi;
    const int NON_LOOP_SITE = -1;
    std::vector<int> siteTypes(ephipsi.size(), NON_LOOP_SITE);
    for(int idx = 0; idx < loopLocs.size(); idx++) {
        auto & loc = loopLocs.at(idx);
        for (int i = loc.idxNTerminal; i <= loc.idxCTerminal; i++) {
            siteTypes.at(i) = idx;
        }
    }
    nppViolation = 0;
    bool result = true;
    for (int i = 0; i < ephipsi.size(); i++) {
        double e = ephipsi.at(i);
        if (e > 0.001) {
            nppViolation++;
            if (siteTypes.at(i) != NON_LOOP_SITE) {
                tmpViolatedLoopsSet.insert(siteTypes.at(i));
            }
            result = false;
        }
    }
    auto & elocal = ec.els;
    nlsViolation = 0;
    for (int i = 0; i < elocal.size(); i++) {
        double e = elocal.at(i);
        if (e > 0.5) {
            nlsViolation++;
            if (siteTypes.at(i) != NON_LOOP_SITE) {
                tmpViolatedLoopsSet.insert(siteTypes.at(i));
            }
            result = false;
        }
        if (i > 0 && i < elocal.size()-1) {
            double e1 = elocal.at(i-1);
            double e2 = elocal.at(i+1);
            double ave = (e1 + e + e2) / 3;
            if (ave > 0.001) {
                nlsViolation++;
                if (siteTypes.at(i) != NON_LOOP_SITE) {
                    tmpViolatedLoopsSet.insert(siteTypes.at(i));
                }
                result = false;
            }
        }
    }

    return result;
}


bool acceptLoop(double deltaE, double mcT) {
    if (deltaE <= 0) {
        return true;
    }
    double pAc = exp(-1.0 * deltaE / mcT);
    double r = 1.0 * rand() / RAND_MAX;
    return r < pAc;
}

std::string genFixedSitesString(const Protein & sites, const std::vector<LoopLocation> & locs) {
    int curChainIdx = 0;
    int curSiteIdx = 0;
    std::vector<LoopLocation> fixedRegions;
    for (int i = 0; i < locs.size(); i++) {
        auto & loc = locs.at(i);
        if (curChainIdx != loc.idxChain) {
            for (int c = curChainIdx; c < loc.idxChain; c++) {
                LoopLocation locFixed;
                locFixed.idxChain = c;
                locFixed.idxNTerminal = 0;
                locFixed.idxCTerminal = sites.at(c).size() - 1;
                fixedRegions.push_back(locFixed);
            }
            curChainIdx = loc.idxChain;
            curSiteIdx = 0;
        }
        LoopLocation locFixed;
        locFixed.idxChain = curChainIdx;
        locFixed.idxNTerminal = curSiteIdx;
        locFixed.idxCTerminal = loc.idxNTerminal - 1;
        fixedRegions.push_back(locFixed);
        curSiteIdx = loc.idxCTerminal + 1;
    }
    if (curSiteIdx < sites.at(curChainIdx).size() - 1) {
        LoopLocation locFixed;
        locFixed.idxChain = curChainIdx;
        locFixed.idxNTerminal = curSiteIdx;
        locFixed.idxCTerminal = sites.at(curChainIdx).size() - 1;
        fixedRegions.push_back(locFixed);
    }
    curChainIdx++;
    for (int c = curChainIdx; c < sites.size(); c++) {
        LoopLocation locFixed;
        locFixed.idxChain = c;
        locFixed.idxNTerminal = 0;
        locFixed.idxCTerminal = sites.at(c).size() - 1;
        fixedRegions.push_back(locFixed);
    }
    std::string fixedsitesstr = "";
    for (const auto & loc : fixedRegions) {
        fixedsitesstr.append("chain");
        fixedsitesstr.append(std::to_string(loc.idxChain));
        fixedsitesstr.append(" ");
        fixedsitesstr.append(std::to_string(loc.idxNTerminal));
        fixedsitesstr.append("-");
        fixedsitesstr.append(std::to_string(loc.idxCTerminal));
        fixedsitesstr.append(", ");
    }
    fixedsitesstr = fixedsitesstr.substr(0, fixedsitesstr.size()-2);
    return fixedsitesstr;
}

void exploreLoops(const std::string & controlfile) {
    std::ostream & oslog = std::cout;

    ///
    std::string controlname = "sdffcontrol";
    genchainreadcontrols(controlfile, controlname);
    std::string sdiocontrolname = controlname + "_sdinputoutput";
    SDInputOutput::readiocontrolset(controlfile, sdiocontrolname);
    std::string lpexplcontrolname = controlname + "_loopexploreparam";
    LoopExploreParam::readcontrolset(controlfile, lpexplcontrolname);

    SDInputOutput sdio(sdiocontrolname);
    LoopExploreParam lpexplparam(lpexplcontrolname);

    unsigned int seed=std::abs(sdio.sdseed());
    NSPdstl::RandomEngine<>::getinstance().reseed(seed);

    std::string gccontrolname = controlname + "_genchain";
    GenChain genchain(gccontrolname);

    std::vector<int> eneterms2print = {ForceField::ENECOMP::ETOT, ForceField::ENECOMP::EBOND,
            ForceField::ENECOMP::EANG,ForceField::ENECOMP::EIMPDIH,ForceField::ENECOMP::EPHIPSI,
            ForceField::ENECOMP::ESCCONF,ForceField::ENECOMP::ESTERIC,ForceField::ENECOMP::ESCPACKING,
            ForceField::ENECOMP::ELOCALSTRUCTURE, ForceField::ENECOMP::ELOCALHB,ForceField::ENECOMP::ESITEPAIRS,
            ForceField::ENECOMP::ESTRUCTREST,ForceField::ENECOMP::ERGRESTRAINT,ForceField::ENECOMP::ESSRESTRAINT
    };

    std::vector<double> sdendcrds;
    std::vector<double> sdendenes;
    std::ofstream ofsAllKeptLoops(lpexplparam.allKeptLoopsFile());
    std::ofstream ofsNoneViolatedLoops(lpexplparam.noneViolatedLoopsFile());

    /// get loop locations
    const std::vector<int> & oriLoopDesSeq = lpexplparam.loopDescriptors();
    std::vector<LoopDescriptor> oriLoopDess;
    for (int i = 0; i < oriLoopDesSeq.size(); i += 4) {
        LoopDescriptor des;
        des.loopLoc.idxChain = oriLoopDesSeq[i];
        des.loopLoc.idxNTerminal = oriLoopDesSeq[i+1];
        des.loopLoc.idxCTerminal = oriLoopDesSeq[i+2];
        des.newLen = oriLoopDesSeq[i+3];
        oriLoopDess.push_back(des);
    }

    /// write reference PDB block
    ofsAllKeptLoops << "START_PDBCONF" << std::endl;
    genchain.writepdb(*genchain.initcrd(), ofsAllKeptLoops, 1.0/A2NM);
    ofsAllKeptLoops << "END_PDBCONF" << std::endl;
    ofsNoneViolatedLoops << "START_PDBCONF" << std::endl;
    genchain.writepdb(*genchain.initcrd(), ofsNoneViolatedLoops, 1.0/A2NM);
    ofsNoneViolatedLoops << "END_PDBCONF" << std::endl;

    /// echo contorls, only print once
    oslog << "Echo of input control parameters:" << std::endl;
    oslog << std::endl;
    SDInputOutput::printiocontrol(controlname, oslog);
    genchainprintcontrols(controlname, oslog);
    LoopOptimizeParam::printcontrol(controlname, oslog);
    oslog << std::endl;
    oslog << "The order of terms in the printed potential energy lists:" << std::endl;
    oslog << "TOTAL, bond, angle, improper_dihedral, ramachandran, rotamer, "
          << "steric, sidechain_packing, mainchain_local_correlation, "
          << "mainchain_local_hydrogenbond, mainchain_nonlocal_packing, "
          << "other_structure_restraint, radius_of_gyration_restraint, "
          << "secondary_structure_restraint" << std::endl;
    oslog << std::endl;
    ///

    /// start explore
    oslog << "[Info] Start Exploring" << std::endl;
    std::map<std::vector<double>, std::vector<double>> allKeptLoops;
    EneComp ec;
    int minLoopLen = lpexplparam.minLoopLen();
    if ( minLoopLen < 0 /* not set */
         || minLoopLen < 3   ) {
        minLoopLen = 3;
    }

    /// enumerate all single-loop descriptors
    std::vector<LoopDescriptor> newLoopDess;
    for (auto & des : oriLoopDess) {
        int oriLen = des.loopLoc.idxCTerminal - des.loopLoc.idxNTerminal + 1;
        int ntEdgeBack = des.loopLoc.idxNTerminal - lpexplparam.maxShiftLen();
        int ntEdgeForward = des.loopLoc.idxNTerminal;
        int ctEdgeBack = des.loopLoc.idxCTerminal;
        int ctEdgeForward = des.loopLoc.idxCTerminal + lpexplparam.maxShiftLen();
        for (int n = ntEdgeBack; n <= ntEdgeForward; n++) {
            int ntShift = des.loopLoc.idxNTerminal - n;
            for (int c = ctEdgeBack; c <= ctEdgeForward; c++) {
                int ctShift = c - des.loopLoc.idxCTerminal;
                int maxLen = oriLen + ntShift + ctShift + lpexplparam.maxExtendLen();
                if (lpexplparam.maxLoopLen() > 0
                        && maxLen > lpexplparam.maxLoopLen()) {
                    maxLen = lpexplparam.maxLoopLen();
                }
                for (int len = minLoopLen; len <= maxLen; len++) {
                    LoopDescriptor d;
                    d.loopLoc.idxChain = des.loopLoc.idxChain;
                    d.loopLoc.idxNTerminal = n;
                    d.loopLoc.idxCTerminal = c;
                    d.newLen = len;
                    newLoopDess.push_back(d);
                    oslog << "[Info][EnumerateLoopDescriptors] " << d.toString() << std::endl;
                }
            }
        }
    }

    Protein refconf = *genchain.initpdb();

    /// explore and locally optimize loops for each possible loop descriptor
    for (auto & des : newLoopDess) {
        Protein newSites;
        std::vector<LoopDescriptor> curLoopDess;
        curLoopDess.push_back(des);
        std::vector<LoopLocation> newLoopLocs;
        int constructResult = constructLoops(refconf, &newSites, curLoopDess, &newLoopLocs,
                   lpexplparam.maxLinkTryTimes());
        if (constructResult != 0) {
            oslog << "[Warn][DiscardInvalidLoopDescriptor] " << des.toString() << std::endl;
            continue;
        }
        genchain.resetStartConf(newSites);

        /// setup fixed sites
        std::string fixsitesstr = genFixedSitesString(newSites, newLoopLocs);
        auto & pset=GenChainControls::getparameterset(gccontrolname);
        pset.StringPar.replace("AllFixedSites", fixsitesstr);

        int numContDiscarded = 0;
        int numAccepted = 0;
        bool filterEverPassed = false;
        // store energies and loop coordinates
        std::vector<int> loopAtomIndices;
        for (int s = 0; s < lpexplparam.maxNumOptmLoopsEach(); s++) {
            std::vector<double> loopCrds;
            int loopNum = s+1;
            oslog << "[Info][OptimizeLoop] optimize loop " << loopNum << std::endl;
            Protein startconf;
            int sdTryTimes = 0;
            bool sdConvergent = false;
            while (!sdConvergent) {
                std::vector<LoopLocation> locs;
                constructResult = constructLoops(refconf, &startconf, curLoopDess, &locs,
                                   lpexplparam.maxLinkTryTimes());
                if (constructResult != 0) {
                    oslog << "[Error][ReconstructLoops] failed to construct loop(s), continue ..." << std::endl;
                    continue;
                }
                genchain.resetStartConf(startconf);
                /// run SD quench
                int sdsteps = runSD(oslog, controlname, genchain, lpexplparam, sdio, eneterms2print,
                      &sdendcrds, &sdendenes, ec);
                sdTryTimes++;
                if (sdsteps > 0) {
                    sdConvergent = true;
                    break;
                } else {
                    if (sdTryTimes >= lpexplparam.maxSDTimesOfOneLoop()) {
                        break;
                    }
                }
            }

            if (!sdConvergent) {
                oslog << "[Error][SDOptimizeLoop][Misconvergence] failed to optimize loop {"
                        << loopNum << "}" << std::endl;
                continue;
            }

            if (loopAtomIndices.empty()) {
                ActiveSelections* acts = genchain.getactiveselections();
                std::vector<bool> atomsFixed = acts->atomfixed();
                int natoms = atomsFixed.size();
                for (int i = 0; i < natoms; i++) {
                    if (!atomsFixed[i]) {
                        loopAtomIndices.push_back(i);
                    }
                }
            }

            for (int idx : loopAtomIndices) {
                loopCrds.push_back(sdendcrds[idx*3]);
                loopCrds.push_back(sdendcrds[idx*3+1]);
                loopCrds.push_back(sdendcrds[idx*3+2]);
            }

            double energy = sdendenes.at(ForceField::ENECOMP::ETOT);
            int nppViolation = 0;
            int nlsViolation = 0;
            bool filterPassed = filter(ec, nppViolation, nlsViolation,
                                       newLoopLocs);

            ofsAllKeptLoops << "LOOPOPTM " << des.toString() << " " << loopNum << std::endl;
            ofsAllKeptLoops << "COORDINATES " << loopCrds.size() / 3 << std::endl;
            for (int i = 0; i < loopCrds.size(); i += 3) {
                ofsAllKeptLoops << loopCrds.at(i) << " " << loopCrds.at(i+1) << " "
                        << loopCrds.at(i+2) << std::endl;
            }
            ofsAllKeptLoops << "ENERGIES";
            for (double d : sdendenes) {
                ofsAllKeptLoops << " " << d;
            }
            ofsAllKeptLoops << std::endl;

            if (filterPassed) {
                ofsNoneViolatedLoops << "LOOPOPTM " << des.toString() << " " << loopNum << std::endl;
                ofsNoneViolatedLoops << "COORDINATES " << loopCrds.size() / 3 << std::endl;
                for (int i = 0; i < loopCrds.size(); i += 3) {
                    ofsNoneViolatedLoops << loopCrds.at(i) << " " << loopCrds.at(i+1) << " "
                            << loopCrds.at(i+2) << std::endl;
                }
                ofsNoneViolatedLoops << "ENERGIES";
                for (double d : sdendenes) {
                    ofsNoneViolatedLoops << " " << d;
                }
                ofsNoneViolatedLoops << std::endl;
            }
        } // end exploration
    }
    ofsAllKeptLoops.close();
    ofsNoneViolatedLoops.close();
}

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
        return 1;
    }

    std::string controlfile(argv[1]);

    std::cout << "[Info][Main] loop exploration started at " << getCurrentLocalTimeString()
                << std::endl;

    exploreLoops(controlfile);

    std::cout << "[Info][Main] loop exploration finished at " << getCurrentLocalTimeString()
            << std::endl;

    return 0;
}


