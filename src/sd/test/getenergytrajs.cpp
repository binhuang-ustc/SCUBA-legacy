/*
 * getenergytrajs.cpp
 *
 *  Created on: 2019/3/14
 *      Author: hbfrank
 */

#include "sd/trajio.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace NSPsd;

void printUsage(const std::string & selfname)
{
    std::cout << "Usage:" << std::endl;
    std::cout << "$ " << selfname << " <TrajFile> [Interval]"
            << std::endl;
}

int main(int argc, char **argv) {
    std::string selfname = argv[0];
    if (argc < 2) {
        std::cerr << "Wrong number of arguments." << std::endl;
        printUsage(selfname);
        return 1;
    }

    std::string trajfile(argv[1]);

    int interval = 1;
    if (argc == 3) {
        interval = std::stoi(argv[2]);
    }

    TrajIO trajio;
    trajio.createreader(trajfile);
    int nframe = 0;
    std::ofstream ofs;
    ofs.open("EnergyTraj_" + trajfile);
    /*{ForceField::ENECOMP::ETOT, ForceField::ENECOMP::EBOND,
    ForceField::ENECOMP::EANG,ForceField::ENECOMP::EIMPDIH,ForceField::ENECOMP::EPHIPSI,
    ForceField::ENECOMP::ESCCONF,ForceField::ENECOMP::ESTERIC,ForceField::ENECOMP::ESCPACKING,
    ForceField::ENECOMP::ELOCALSTRUCTURE, ForceField::ENECOMP::ELOCALHB,ForceField::ENECOMP::ESITEPAIRS,
    ForceField::ENECOMP::ESTRUCTREST,ForceField::ENECOMP::ERGRESTRAINT,ForceField::ENECOMP::ESSRESTRAINT
    }*/
    std::vector<int> eneterms2print = {0,1,2,3,4,5,6,7,8,9,10,12,13,14}; // except ESHADOW
    while(true) {
        std::vector<double> crdframe;
        std::vector<double> eneframe;
        bool success = trajio.readframe(&crdframe);
        if (success) success = trajio.readframe(&eneframe);
        if (!success) break;
        nframe++;
        if ((nframe % interval) == 0) {
            for (int idx : eneterms2print) {
                ofs << eneframe[idx] << "  ";
            }
            ofs << std::endl;
        }
    }
    ofs.close();
}
