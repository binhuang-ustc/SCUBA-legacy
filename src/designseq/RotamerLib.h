/*
 * RotamerLib.h
 *
 *  Created on: 2017��10��23��
 *      Author: notxp
 */

#ifndef DESIGNSEQ_ROTAMERLIB_H_
#define DESIGNSEQ_ROTAMERLIB_H_

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include "designseq/Rotamer.h"
#include "designseq/RotamerGroup.h"
#include "designseq/RotamerEnergyBBDep.h"
#include "designseq/ResName.h"
#include "dataio/datapaths.h"
#include "designseq/ProteinRep.h"


namespace NSPdesignseq {

using namespace std;

class RotamerLib {
public:
    vector<RotamerGroup*> aaRotGroups;

    RotamerGroup allRots;
    map<string,RotamerEnergyBBDep*> rotEnergyMap;
    map<string, int> rotNameToIndex;
    PhipsiLib ppLib;
    ResName rn;

    RotamerLib();
    RotamerLib(const string& rotLibType);

    RotamerGroup* getAAGroup(string& triName);
    float getRotamerEnergy(string& rotName, Phipsi* pp);
    float getRotamerEnergy(string& rotName, int ppType);

    float getNearestRotamerEnergy(Rotamer* rot, int ppType);

    virtual ~RotamerLib();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_ROTAMERLIB_H_ */
