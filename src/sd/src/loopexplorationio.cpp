/*
 * loopexplorationio.cpp
 *
 *  Created on: 2020.6.28
 *      Author: hbfrank
 */

#include "sd/loopexplorationio.h"
#include "proteinrep/pdbreader.h"
#include "backbone/backbonesite.h"
#include <float.h>
#include <time.h>

using namespace NSPsd;
using namespace NSPproteinrep;

bool LoopLocation::operator < (const LoopLocation &rValue) const {
    if (idxChain < rValue.idxChain) {
        return true;
    }
    if (idxChain == rValue.idxChain) {
        return idxCTerminal < rValue.idxNTerminal;
    }
    return false;
}

bool LoopLocation::operator > (const LoopLocation &rValue) const {
    if (idxChain > rValue.idxChain) {
        return true;
    }
    if (idxChain == rValue.idxChain) {
        return idxNTerminal > rValue.idxCTerminal;
    }
    return false;
}

bool LoopLocation::operator == (const LoopLocation &loc) const {
    return loc.idxChain == idxChain
            && loc.idxNTerminal == idxNTerminal
            && loc.idxCTerminal == idxCTerminal;
}

LoopLocation & LoopLocation::operator = (const LoopLocation &rValue) {
    idxChain = rValue.idxChain;
    idxNTerminal = rValue.idxNTerminal;
    idxCTerminal = rValue.idxCTerminal;
    return *this;
}

bool LoopLocation::overlap(const LoopLocation &loc) const {
    if (idxChain != loc.idxChain) {
        return false;
    }
    int lenTot = (idxCTerminal - idxNTerminal + 1)
                + (loc.idxCTerminal - loc.idxNTerminal + 1);
    return loc.idxCTerminal > idxCTerminal ? (loc.idxCTerminal - idxNTerminal <= lenTot - 1)
            : (idxCTerminal - loc.idxNTerminal <= lenTot - 1);
}

bool LoopLocation::cover(int chain, int site) const {
    if (chain != idxChain) {
        return false;
    }
    return (idxCTerminal >= site) && (idxNTerminal <= site);
}

bool LoopLocation::inside(const LoopLocation & loc) const {
    if (loc.idxChain != idxChain) {
        return false;
    }
    return loc.idxNTerminal < idxNTerminal
            && loc.idxCTerminal > idxCTerminal;
}

bool LoopLocation::include(const LoopLocation & loc) const {
    if (loc.idxChain != idxChain) {
        return false;
    }
    return loc.idxNTerminal > idxNTerminal
            && loc.idxCTerminal < idxCTerminal;
}

std::string LoopLocation::toString() const {
    return std::to_string(idxChain) + " " + std::to_string(idxNTerminal) + " "
            + std::to_string(idxCTerminal);
}

bool LoopDescriptor::operator == (const LoopDescriptor &des) const {
    return des.loopLoc == loopLoc
            && des.newLen == newLen;
}

LoopDescriptor & LoopDescriptor::operator = (const LoopDescriptor &des) {
    loopLoc = des.loopLoc;
    newLen = des.newLen;
    return *this;
}

bool LoopDescriptor::operator < (const LoopDescriptor &des) const {
    return toString() < des.toString();
}

std::string LoopDescriptor::toString() const {
    return std::to_string(loopLoc.idxChain) + " "
            + std::to_string(loopLoc.idxNTerminal) + " "
            + std::to_string(loopLoc.idxCTerminal) + " "
            + std::to_string(newLen);
}

bool LoopOptm::operator < (const LoopOptm &lo) const {
    return toString() < lo.toString();
}

double LoopOptm::rmsdTo(const LoopOptm & other) const {
    const auto & crds1 = coords;
    const auto & crds2 = other.coords;
    assert(crds1.size() == crds2.size());
    int natoms = crds1.size();
    double rmsd2 = 0.0;
    for (int i = 0; i < natoms; i++) {
        const XYZ & crd1 = crds1.at(i);
        const XYZ & crd2 = crds2.at(i);
        rmsd2 += (crd1.x_-crd2.x_)*(crd1.x_-crd2.x_)
                + (crd1.y_-crd2.y_)*(crd1.y_-crd2.y_)
                + (crd1.z_-crd2.z_)*(crd1.z_-crd2.z_);
    }
    rmsd2 /= natoms;
    return std::sqrt(rmsd2);
}

std::string LoopOptm::toString() const {
    return loopDesp.toString() + " "
            + std::to_string(number);
}

void LoopExploration::addLoopOptm(LoopOptm loop) {
    std::vector<LoopOptm> & loops = optmLoops_[loop.loopDesp];
    loops.push_back(loop);
    EnergySortedLoopOptms & sortedLoops = sortedOptmLoops_[loop.loopDesp];
    sortedLoops.insert(std::make_pair(loop.energies, loop));
}

void LoopExploration::updateLoopLocationIndex() {
    loopLocIndex_.clear();
    auto & firstDes = sortedOptmLoops_.begin()->first;
    loopLocIndex_.push_back(std::make_pair(firstDes.loopLoc, std::vector<LoopDescriptor>()));
    loopLocIndex_.back().second.push_back(firstDes);
    for (auto it = std::next(sortedOptmLoops_.begin(), 1); it != sortedOptmLoops_.end(); it++) {
        auto & des = it->first;
        auto & locQ = des.loopLoc;
        bool handled = false;
        for (int i = 0; i < loopLocIndex_.size(); i++) {
            auto & pair = loopLocIndex_.at(i);
            auto & locV = pair.first;
            if (locQ < locV) {
                loopLocIndex_.insert(loopLocIndex_.begin()+i,
                             std::make_pair(locQ, std::vector<LoopDescriptor>{des}));
                handled = true;
                break;
            }
            if (locQ.overlap(locV)) {
                pair.second.push_back(des);
                if (locQ.inside(locV)) {
                    locV = locQ;    // use the most inside location as key
                }
                handled = true;
                break;
            }
        }
        if (!handled) {
            loopLocIndex_.push_back(std::make_pair(locQ,
                                           std::vector<LoopDescriptor>{des}));
        }
    }
}

void LoopExploration::statLoopAreas() {
    if (loopLocIndex_.empty()) {
        updateLoopLocationIndex();
    }
    loopAreas_.clear(); // do a clean process
    for (auto & li : loopLocIndex_) {
        const auto & ll = li.first;
        loopAreas_.push_back(std::make_pair(ll, LoopArea()));
        LoopArea & la = loopAreas_.back().second;
        const auto & lds = li.second;
        for (const auto & ld : lds) {
            const auto & loops = sortedOptmLoops_.at(ld);
            int loopLen = ld.newLen;
            la.insert(std::make_pair(loopLen, std::multimap<double, LoopOptm>()));
            auto & lenCluster = la.at(loopLen);
            for (const auto & lo : loops) {
                double uniEtotal = lo.first[0] / loopLen;   // Etotal/LoopLength
                lenCluster.insert(std::make_pair(uniEtotal, lo.second));
            }
        }
    }
}

void LoopExploration::refineLoopAreas(const int maxLoopsOfLowestE,
                                      const int shorterLenExtend,
                                      const double resRMSDCutoff) {
    refinedLoopAreas_.clear();
    if (loopAreas_.empty()) {
        statLoopAreas();
    }
    for (const auto & la : loopAreas_) {
        const auto & area = la.second;
        refinedLoopAreas_.push_back(std::make_pair(la.first, LoopArea()));
        LoopArea & curArea = refinedLoopAreas_.back().second;
        /// find the lowest uniEnergy loop length
        int lenOfLowestE = -1;
        double lowestE = DBL_MAX;
        for (const auto & le : area) {
            if (le.second.begin()->first < lowestE) {
                lenOfLowestE = le.first;
                lowestE = le.second.begin()->first;
            }
        }
        ///---
        /// use maxLoopsOfLowestE loops of the loop length above
        curArea.insert(std::make_pair(lenOfLowestE, std::multimap<double, LoopOptm>()));
        auto & loops = curArea.at(lenOfLowestE);
        const auto & oriLoops = area.at(lenOfLowestE);
        double upperLimitE = -DBL_MAX;
        int n = 0;
        for (const auto & el : oriLoops) {
            bool isRedundant = false;
            for (const auto & exist : loops) {
                if (el.second.loopDesp == exist.second.loopDesp) {
                    // filter redundant loops by RMSD per residue
                    if (el.second.rmsdTo(exist.second)/lenOfLowestE < resRMSDCutoff) {
                        isRedundant = true;
                        break;
                    }
                }
            }
            if (isRedundant) {
                continue;
            }
            loops.insert(std::make_pair(el.first, el.second));
            n++;
            upperLimitE = el.first; // assign value without comparison here because the items in oriLoops are already sorted by item.first
            if (n >= maxLoopsOfLowestE) {
                break;
            }
        }///---
        /// cover loops of shorter length with lower energy than the highest energy of loop above
        for (int s = 1; s <= shorterLenExtend; s++) {
            int lenExt = lenOfLowestE - s;
            if (area.find(lenExt) == area.end()) {
                continue;
            }
            const auto & extLoops = area.at(lenExt);
            for (const auto & el : extLoops) {
                if (el.first <= upperLimitE) {
                    bool isRedundant = false;
                    for (const auto & exist : loops) {
                        if (el.second.loopDesp == exist.second.loopDesp) {
                            // filter redundant loops by RMSD per residue
                            if (el.second.rmsdTo(exist.second)/lenOfLowestE < resRMSDCutoff) {
                                isRedundant = true;
                                break;
                            }
                        }
                    }
                    if (isRedundant) {
                        continue;
                    }
                    loops.insert(std::make_pair(el.first, el.second));
                } else {
                    break;
                }
            }
        }
        ///---
    }
}

EnergySortedProteins LoopExploration::assembleLoops(const int numConfs) {
    EnergySortedProteins confs;
    if (refinedLoopAreas_.empty()) {
        return confs;
    }
    std::vector<std::vector<LoopOptm>> refinedMixLoops;
    for (const auto & area : refinedLoopAreas_) {
        refinedMixLoops.push_back(std::vector<LoopOptm>());
        auto & loops = refinedMixLoops.back();
        for (const auto & lloops : area.second) {
            for (const auto & eloop : lloops.second) {
                loops.push_back(eloop.second);
            }
        }
    }
    int actualNumConfs = numConfs;
    long numUpperLimit = 1;
    for (const auto & loops : refinedMixLoops) {
        numUpperLimit *= loops.empty() ? 1 : loops.size();
        if (numUpperLimit > actualNumConfs) {
            break;
        }
    }
    if (numUpperLimit < actualNumConfs) {
        actualNumConfs = (int)numUpperLimit;
    }
    std::set<std::string> visitedCombinations;
    int n = 0;
    srand((unsigned int)(time(NULL)));
    while (n < actualNumConfs) {
        std::string loopcomb = "";
        std::vector<LoopOptm> selectedLoops;
        for (const auto & area : refinedMixLoops) {
            int selindex = rand() % area.size();
            const auto & loop = area.at(selindex);
            selectedLoops.push_back(loop);
            loopcomb += loop.toString() + " ";
        }
        if (visitedCombinations.find(loopcomb) != visitedCombinations.end()) {
            continue;
        } else {
            visitedCombinations.insert(loopcomb);
        }
        double totalLoopE = 0.0;
        int totalLoopLen = 0;
        std::vector<double> pseudoE;
        for (const auto & loop : selectedLoops) {
            totalLoopE += loop.energies.at(0);
            totalLoopLen += loop.loopDesp.newLen;
            pseudoE.push_back(loop.loopDesp.loopLoc.idxChain);
            pseudoE.push_back(loop.loopDesp.loopLoc.idxNTerminal);
            pseudoE.push_back(loop.loopDesp.loopLoc.idxCTerminal);
            pseudoE.push_back(loop.loopDesp.newLen);
        }
        double uniTotalLoopE = totalLoopE / totalLoopLen;
        pseudoE.insert(pseudoE.begin(), uniTotalLoopE);
        Protein conf = restorePDBFromLoopOptms(selectedLoops);
        confs.insert(std::make_pair(pseudoE, conf));
        n++;
    }
    return confs;
}

Protein LoopExploration::restorePDBFromLoopOptm(const LoopOptm & loop) const {
    Protein pr;
    int idxChain = loop.loopDesp.loopLoc.idxChain;
    int idxNT = loop.loopDesp.loopLoc.idxNTerminal;
    int idxCT = loop.loopDesp.loopLoc.idxCTerminal;
    int loopLen = loop.loopDesp.newLen;
    for (int c = 0; c < idxChain; c++) {
        pr.push_back(refPDB_.at(c));
    }
    const Peptide & refchain = refPDB_.at(idxChain);
    Peptide newchain;
    for (int i = 0; i < idxNT; i++) {
        newchain.push_back(refchain.at(i));
    }
    std::vector<XYZ> bbCoords;
    for (int i = 0; i < loopLen; i++) {
        bbCoords.clear();
        bbCoords.push_back(loop.coords[i*4+0]);
        bbCoords.push_back(loop.coords[i*4+1]);
        bbCoords.push_back(loop.coords[i*4+2]);
        bbCoords.push_back(loop.coords[i*4+3]);
        BackBoneSite bs;
        bs.changecrd(bbCoords);
        bs.resname = "GLY";
        newchain.push_back(make_fullsite(bs));
    }
    for (int i = idxCT+1; i < refchain.size(); i++) {
        newchain.push_back(refchain.at(i));
    }
    pr.push_back(newchain);
    for (int c = idxChain+1; c < refPDB_.size(); c++) {
        pr.push_back(refPDB_.at(c));
    }
    return pr;
}

Protein LoopExploration::restorePDBFromLoopOptms(const std::vector<LoopOptm> & loops) const {
    Protein pr;
    int lastChain = -1;
    int lastRes = -1;
    for (const auto & loop : loops) {
        int idxChain = loop.loopDesp.loopLoc.idxChain;
        if (idxChain != lastChain) {
            if (lastChain >= 0 && lastChain < refPDB_.size()) {
                for (int r = lastRes + 1; r < refPDB_.at(lastChain).size(); r++) {
                    pr.back().push_back(refPDB_.at(lastChain).at(r));
                }
                lastRes = -1;
            }
            for (int c = lastChain+1; c < idxChain; c++) {
                pr.push_back(refPDB_.at(c));
            }
            lastChain = idxChain;
            pr.push_back(Peptide());
        }
        int idxNT = loop.loopDesp.loopLoc.idxNTerminal;
        int idxCT = loop.loopDesp.loopLoc.idxCTerminal;
        int loopLen = loop.loopDesp.newLen;
        const Peptide & refchain = refPDB_.at(idxChain);
        Peptide & curChain = pr.back();
        for (int i = lastRes+1; i < idxNT; i++) {
            curChain.push_back(refchain.at(i));
        }
        std::vector<XYZ> bbCoords;
        for (int i = 0; i < loopLen; i++) {
            bbCoords.clear();
            bbCoords.push_back(loop.coords[i*4+0]);
            bbCoords.push_back(loop.coords[i*4+1]);
            bbCoords.push_back(loop.coords[i*4+2]);
            bbCoords.push_back(loop.coords[i*4+3]);
            BackBoneSite bs;
            bs.changecrd(bbCoords);
            bs.resname = "GLY";
            curChain.push_back(make_fullsite(bs));
        }
        lastRes = idxCT;
    }
    if (lastChain >= 0 && lastChain < refPDB_.size()) {
        for (int r = lastRes + 1; r < refPDB_.at(lastChain).size(); r++) {
            pr.back().push_back(refPDB_.at(lastChain).at(r));
        }
    }
    for (int c = lastChain+1; c < refPDB_.size(); c++) {
        pr.push_back(refPDB_.at(c));
    }
    return pr;
}

bool LoopExplorationIO::read(const std::string & file) {
    std::ifstream ifs(file);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open file: " << file << std::endl;
        return false;
    }
    return read(ifs);
}

/**
 * TEMP use
 * copied from fullsite/src/fullsite.cpp: NSPproteinrep::readfullsitesfrompdb
 */
Protein fullsitesFromReader(PdbReader & reader) {
    Protein protein;
    std::string chainids=reader.chainids();
    std::vector<FullSite>  sites;
    for(int i=0;i<chainids.size();++i){
        char chainid=chainids[i];
        std::vector<std::string> seq=reader.getaminoacidsequence(chainid);
        for(int resinumber=0;resinumber<seq.size();++resinumber){
            typename PdbReader::ResKeyType reskey=reader.mappdbkeyint()->pdbResKey(resinumber,i);
            std::vector<PdbRecord> &records=reader.records().at(chainid).at(reskey);
            if(records[0].residuename=="MSE"){
                for(auto &r:records){
                    r.residuename="MET";
                    if(r.atomname=="SE") r.atomname="SD";
                }
            }
            if(!SideChainTopo::residuetypedefined(records[0].residuename)) continue;
            FullSite fs;
            fs.resname()=records[0].residuename;
            fs.resid()=records[0].residueid;
            fs.chainid()=records[0].chainid;
            int  nmainchain=0;
            for(auto &r:records) {
                NSPgeometry::XYZ crd(r.x,r.y,r.z);
                if(FullSite::mainchainatom(r.atomname)) nmainchain++;
                fs.addatomcrd(r.atomname,crd);
            }
            if(nmainchain>=4) sites.push_back(fs);
        }
    }
    protein.assign(1,std::vector<FullSite>());
    int chainnumber=0;
    int resseq=0;
    for(auto it=sites.begin();it !=sites.end()-1;++it){
        protein[chainnumber].push_back(*it);
        it->resseq()=resseq++;
        NSPgeometry::XYZ c0=it->getcrd("C");
        NSPgeometry::XYZ n1=(it+1)->getcrd("N");
        double bond = NSPgeometry::distance(c0, n1);
        if (bond > 2.0){
            protein.push_back(std::vector<FullSite>());
            ++chainnumber;
            resseq=0;
        }
    }
    protein[chainnumber].push_back(sites.back());
    return protein;
}

bool lineStartMatch(const std::string & line, const std::string & target) {
    int len = target.size();
    return line.substr(0, len) == target;
}

bool LoopExplorationIO::read(std::ifstream & ifs) {
    std::string line;
    std::istringstream iss;
    while (std::getline(ifs, line)) {
        if (lineStartMatch(line, "START_PDBCONF")) {
            PdbReader pdbReader;
            std::vector<std::string> pdbLines;
            while (std::getline(ifs, line)) {
                if (lineStartMatch(line, "END_PDBCONF")) {
                    break;
                }
                pdbLines.push_back(line);
            }
            pdbReader.readpdb(pdbLines);
            lpExploration_.setRefPDB(fullsitesFromReader(pdbReader));
            continue;
        }
        if (lineStartMatch(line, "LOOPOPTM")) {
            LoopOptm lo;
            iss = std::istringstream(line.substr(9));
            int idxChain, idxNT, idxCT, newLen, loopNum;
            iss >> idxChain >> idxNT >> idxCT >> newLen >> loopNum;
            lo.loopDesp.loopLoc.idxChain = idxChain;
            lo.loopDesp.loopLoc.idxNTerminal = idxNT;
            lo.loopDesp.loopLoc.idxCTerminal = idxCT;
            lo.loopDesp.newLen = newLen;
            lo.number = loopNum;
            std::getline(ifs, line);
            int nline = std::stod(line.substr(12)); //"COORDINATES"
            for (int i = 0; i < nline; i++) {
                std::getline(ifs, line);
                iss = std::istringstream(line);
                double x, y, z;
                iss >> x >> y >> z;
                XYZ xyz(x, y, z);
                lo.coords.push_back(xyz);
            }
            std::getline(ifs, line);
            iss = std::istringstream(line.substr(9)); //"ENERGIES"
            double e;
            while (iss >> e) {
                lo.energies.push_back(e);
            }
            lpExploration_.addLoopOptm(lo);
            continue;
        }
    }
    lpExploration_.updateLoopLocationIndex();
    return true;
}

