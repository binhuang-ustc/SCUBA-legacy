/*
 * traj2pdbs.cpp
 *
 *  Created on: 2018/4/2
 *      Author: hbfrank
 */

#include "sd/genchain.h"
#include "sd/trajio.h"
#include <iostream>
#include <fstream>
using namespace NSPproteinrep;
using namespace NSPsd;
using namespace NSPgeometry;

void printUsage(const std::string & selfname)
{
    std::cout << "Usage:" << std::endl;
    std::cout << "$ " << selfname << " <RunSD_Parameter_File> <TrajFile> [Interval]"
            << std::endl;
}

int main(int argc, char **argv) {
    std::string selfname = argv[0];
    if (argc < 3) {
        std::cerr << "Wrong number of arguments." << std::endl;
        printUsage(selfname);
        return 1;
    }

    std::string sdparamfile(argv[1]);
    std::string controlname="sdffcontrol";
    genchainreadcontrols(sdparamfile, controlname);
    GenChain genchain(controlname + "_genchain");

    std::string trajfile(argv[2]);

    int interval = 1;
    if (argc == 4) {
        interval = std::stoi(argv[3]);
    }

    TrajIO trajio;
    trajio.createreader(trajfile);
    int nframe = 0;
    while(true) {
        std::vector<double> crdframe;
        std::vector<double> eneframe;
        bool success = trajio.readframe(&crdframe);
        if (success) success = trajio.readframe(&eneframe);
        if (!success) break;
        nframe++;
        if ((nframe % interval) == 0) {
            std::ofstream ofs;
            ofs.open(trajfile + "_" + std::to_string(nframe) + ".pdb");
            genchain.writepdb(crdframe, ofs, 10.0);
            ofs.close();
        }
    }
}
