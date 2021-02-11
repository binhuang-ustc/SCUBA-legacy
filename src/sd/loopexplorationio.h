/*
 * loopexplorationio.h
 *
 *  Created on: 2020.6.28
 *      Author: hbfrank
 */

#ifndef SD_LOOPEXPLORATIONIO_H_
#define SD_LOOPEXPLORATIONIO_H_

#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include "fullsite/fullsite.h"
#include "geometry/xyz.h"
#include "sd/forcefield.h"

namespace NSPsd{
using namespace NSPproteinrep;
using namespace NSPgeometry;

typedef std::vector<FullSite> Peptide;
typedef std::vector<std::vector<FullSite>> Protein;
typedef std::vector<BackBoneSite> BackboneChain;
typedef std::pair<BackBoneSite, BackBoneSite> AnchorPair;

class LoopLocation {
public:
    int loopNumber;
    int idxChain;
    int idxNTerminal;
    int idxCTerminal;
    LoopLocation() : loopNumber(-1), idxChain(-1),
            idxNTerminal(-1), idxCTerminal(-1) {;}
    bool operator == (const LoopLocation &loc) const;
    bool operator < (const LoopLocation &rValue) const;
    bool operator > (const LoopLocation &rValue) const;
    LoopLocation & operator = (const LoopLocation &rValue);
    bool overlap(const LoopLocation &loc) const;
    bool cover(int idxChain, int idxSite) const;
    bool inside(const LoopLocation &loc) const;
    bool include(const LoopLocation &loc) const;
    std::string toString() const;
};

class LoopDescriptor {
public:
    LoopLocation loopLoc;
    int newLen;
    LoopDescriptor(): newLen(-1) {;}
    bool operator == (const LoopDescriptor &des) const;
    LoopDescriptor & operator = (const LoopDescriptor &des);
    bool operator < (const LoopDescriptor &des) const;
    std::string toString() const;
private:
};

class LoopOptm {
public:
    int number;
    LoopDescriptor loopDesp;
    std::vector<XYZ> coords = std::vector<XYZ>();
    std::vector<double> energies = std::vector<double>();
    LoopOptm(): loopDesp(), number(-1) {;}
    bool operator < (const LoopOptm &lo) const;
    double rmsdTo(const LoopOptm & other) const;
    std::string toString() const;
private:
};

class EnergyLessThan : std::vector<double> {
public:
    bool operator() (const std::vector<double> & ene1,
                    const std::vector<double> & ene2) const {
        if (ene1.size() == 0 || ene2.size() == 0
                || ene1.size() != ene2.size()) {
            return false; // trick for exceptions
        }
        if (ene1[ForceField::ENECOMP::ETOT] < ene2[ForceField::ENECOMP::ETOT]) {
            return true;
        }
        return false;
    }
};
typedef std::multimap<std::vector<double>, LoopOptm,
        EnergyLessThan> EnergySortedLoopOptms;
typedef std::multimap<std::vector<double>, Protein,
        EnergyLessThan> EnergySortedProteins;
typedef std::map<int, std::multimap<double, LoopOptm>> LoopArea;

class LoopExploration {
public:
    LoopExploration() {;}
    void clear() { optmLoops_.clear(); }
    Protein restorePDBFromLoopOptm(const LoopOptm & loop) const;
    Protein restorePDBFromLoopOptms(const std::vector<LoopOptm> & loops) const;
    void setRefPDB(Protein pdb) { refPDB_ = pdb; }
    void addLoopOptm(LoopOptm loop);
    void updateLoopLocationIndex();
    void statLoopAreas();
    void refineLoopAreas(const int maxLoopsOfLowestE, const int shorterLenExtend = 0,
                         const double resRMSDCutoff = 0.1);
    EnergySortedProteins assembleLoops(const int numConfs = 1);
    const std::map<LoopDescriptor, std::vector<LoopOptm>> & loops() const { return optmLoops_; }
    const std::map<LoopDescriptor, EnergySortedLoopOptms> & sortedLoops() const { return sortedOptmLoops_; }
    const std::vector<std::pair<LoopLocation, std::vector<LoopDescriptor>>> & loopLocationIndex() const { return loopLocIndex_; }
    const std::vector<std::pair<LoopLocation, LoopArea>> & loopAreas() const { return loopAreas_; }
    const std::vector<std::pair<LoopLocation, LoopArea>> & refinedLoopAreas() const { return refinedLoopAreas_; }
private:
    Protein refPDB_;
    std::map<LoopDescriptor, std::vector<LoopOptm>> optmLoops_;
    std::map<LoopDescriptor, EnergySortedLoopOptms> sortedOptmLoops_;
    std::vector<std::pair<LoopLocation, std::vector<LoopDescriptor>>> loopLocIndex_;   // to maintain representative loop locations in order
    std::vector<std::pair<LoopLocation, LoopArea>> loopAreas_;
    std::vector<std::pair<LoopLocation, LoopArea>> refinedLoopAreas_;
    std::map<int, std::vector<LoopOptm>> loopCandidates_;    // key: serial number of loop; value: loop candidates belongs to the loop of the key
    std::map<int, std::map<int, double>> loopUnitEneDistr_; // key: serial number of loop; value: loop_length=>unit_energy pairs
};

class LoopExplorationIO {
public:
    LoopExplorationIO(){;}
    bool read(const std::string & file);
    bool read(std::ifstream & ifs);
    LoopExploration & exploration() { return lpExploration_; }
    bool write(const std::string & file) const;
    bool write(std::ofstream ofs) const;
private:
    LoopExploration lpExploration_;
};

}



#endif /* SD_LOOPEXPLORATIONIO_H_ */
