/*
 * processloopexplorationresult.cpp
 *
 *  Created on: 2020.6.28
 *      Author: hbfrank
 */

#include "proteinrep/pdbrecord.h"
#include "sd/forcefield.h"
#include "sd/loopexplorationio.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <assert.h>
#include <regex>
#include <cstdio>
#include <float.h>

using namespace NSPproteinrep;
using namespace NSPsd;

std::vector<std::string> split(const std::string & line) {
    std::regex ws_re("\\s+");
    std::vector<std::string> v(
            std::sregex_token_iterator(line.begin(), line.end(), ws_re, -1),
            std::sregex_token_iterator());
    return v;
}

double rmsdBetweenLoops(std::vector<XYZ> & crds1, std::vector<XYZ> & crds2) {
    assert(!crds1.empty());
    assert(crds1.size() == crds2.size());
    int natoms = crds1.size();
    double rmsd2 = 0.0;
    for (int i = 0; i < natoms; i++) {
        XYZ & crd1 = crds1.at(i);
        XYZ & crd2 = crds2.at(i);
        rmsd2 += crd1.squaredDistance(crd2);
    }
    rmsd2 /= natoms;
    return std::sqrt(rmsd2);
}

LoopExploration readExploration(const std::string & file) {
    LoopExplorationIO io;
    io.read(file);
    return io.exploration();
}

void stat(LoopExploration & expl) {
    if (expl.loopAreas().empty()) {
        expl.statLoopAreas();
    }
    for (int li = 0; li < expl.loopLocationIndex().size(); li++) {
        auto & area = expl.loopAreas()[li].second;
        for (const auto & lc : area) {
            const auto & loops = lc.second;
            for (const auto & lo : loops) {
                std::cout << "LoopArea " << (li+1) << " :" << lo.second.toString() << " " << lo.first << std::endl;
            }
        }
        /// find the lowest uniEnergy loop length
        int lenOfLowestE = -1;
        double lowestE = DBL_MAX;
        for (const auto & le : area) {
            if (le.second.begin()->first < lowestE) {
                lenOfLowestE = le.first;
                lowestE = le.second.begin()->first;
            }
        }
        std::cout << "LoopArea " << (li+1) << " :LenOfLowestE: " << lenOfLowestE << std::endl;
        ///---
    }
}

int serve(LoopExploration & expl) {
    std::cout << "> ";
    std::map<LoopDescriptor, Protein> prTemplates;
    std::string cmdline;
    while(std::getline(std::cin, cmdline)) {
        bool sortEnergy = true;
        std::vector<std::string> args = split(cmdline);
        if (args.empty()) {
            // do nothing
        } else if ("help" == args[0]) {
            // print help
            std::cout << "Supported commands:" << std::endl;
            std::cout << "quit/exit : exit program" << std::endl;
            std::cout << "sort on/off : use sorted loops (on) or serial loops (off)"
                        << std::endl;
            std::cout << "list [all] : print summary of the dataset" << std::endl;
            std::cout << "list LOOPDESCRIPTOR : print loops with TotalEnergy "
                        "belong to the specified LOOPDESCRIPTOR" << std::endl;
            std::cout << "show LOOPOPTM : print coordinates and energies of the "
                        "specified LOOPOPTM" << std::endl;
            std::cout << "show eneterms : print names of energy terms" << std::endl;
            std::cout << "save LOOPOPTM PDBfilename : restore total conformation "
                        "based on the specified LOOPOPTM and save as PDB format "
                        "to file PDBfilename" << std::endl;
            std::cout << "rmsd LOOPDESCRIPTOR loopNum1 loopNum2 : print RMSD "
                        "between loopNum1 and loopNum2 belong to the specified "
                        "LOOPDESCRIPTOR" << std::endl;
            std::cout << "merge datafile1 [datafile2 datafile3 ...] : merge data "
                        "in datafile1 (and datafile2 datafile3 ... if more file "
                        "specified) into the existing dataset" << std::endl;
            std::cout << "stat [all]/loopArea1 [loopArea2 loopArea3 ...] : print "
                        "LoopOptm ID with its unified Etotal line by line of "
                        "specified loop area(s)" << std::endl;
            std::cout << "refine [maxLoopsEachArea shorterLenExtend resRMSDCutoff] : "
                        "refine loops in each loop area(s) according to the optional "
                        "params" << std::endl;
            std::cout << "assemble numConfs outputDir [pdbNamePrefix] : pick loops "
                        "from each loop area and assemble to a complete conformation "
                        "and write to PDB file" << std::endl;
            std::cout << "" << std::endl;
        } else if ("quit" == args[0] || "exit" == args[0]) {
            std::cout << "Bye" << std::endl;
            return 0;
        } else if ("sort" == args[0] && args.size() == 2) {
            if ("on" == args[1]) {
                sortEnergy = true;
            } else if ("off" == args[1]) {
                sortEnergy = false;
            } else {
                std::cerr << "Invalid sort operation: " << args[1] << std::endl;
            }
        } else if ("list" == args[0]) {
            if (args.size() == 1 || "all" == args[1]) {
                std::cout << "This dataset contains " << expl.loopLocationIndex().size()
                        << " loop regions: " << std::endl;
                for (const auto & li : expl.loopLocationIndex()) {
                    std::cout << "LoopLocation: " << li.first.toString() << std::endl;
                }
                std::cout << "And there are " << expl.loops().size() << " loop descriptors totally: " << std::endl;
                for (const auto & dl : expl.loops()) {
                    std::cout << "LOOPDESCRIPTOR " << dl.first.toString()
                            << " has " << dl.second.size() << " loops" << std::endl;
                }
            } else {
                LoopDescriptor ldes;
                ldes.loopLoc.idxChain = std::stoi(args[1]);
                ldes.loopLoc.idxNTerminal = std::stoi(args[2]);
                ldes.loopLoc.idxCTerminal = std::stoi(args[3]);
                ldes.newLen = std::stoi(args[4]);
                if (sortEnergy) {
                    if (expl.sortedLoops().find(ldes) == expl.sortedLoops().end()) {
                        std::cerr << "LOOPDESCRIPTOR record not found" << std::endl;
                    } else {
                        auto sloops = expl.sortedLoops().at(ldes);
                        for (auto & slp : sloops) {
                            auto & lp = slp.second;
                            std::cout << "LOOPOPTM " << lp.toString() << " TotEne:"
                                    << lp.energies.at(0) << std::endl;
                        }
                    }
                } else {
                    if (expl.loops().find(ldes) == expl.loops().end()) {
                        std::cerr << "LOOPDESCRIPTOR record not found" << std::endl;
                    } else {
                        auto loops = expl.loops().at(ldes);
                        for (auto & lp : loops) {
                            std::cout << "LOOPOPTM " << lp.toString() << " TotEne:"
                                    << lp.energies.at(0) << std::endl;
                        }
                    }
                }
            }
        } else if ("show" == args[0] && args.size() == 6) {
            LoopDescriptor ldes;
            ldes.loopLoc.idxChain = std::stoi(args[1]);
            ldes.loopLoc.idxNTerminal = std::stoi(args[2]);
            ldes.loopLoc.idxCTerminal = std::stoi(args[3]);
            ldes.newLen = std::stoi(args[4]);
            int loopNum = std::stoi(args[5]);
            if (sortEnergy){
                if (expl.sortedLoops().find(ldes) == expl.sortedLoops().end()) {
                    std::cerr << "LOOPDESCRIPTOR record not found" << std::endl;
                } else {
                    auto sloops = expl.sortedLoops().at(ldes);
                    for (auto & slp : sloops) {
                        auto & lp = slp.second;
                        if (lp.number == loopNum) {
                            std::cout << "LOOPOPTM " << lp.toString() << std::endl;
                            std::cout << "COORDINATES " << lp.coords.size() << std::endl;
                            for (const XYZ & xyz : lp.coords) {
                                std::cout << xyz.x_ << " " << xyz.y_ << " " << xyz.z_ << std::endl;
                            }
                            std::cout << "ENERGIES";
                            for (double d : lp.energies) {
                                std::cout << " " << d;
                            }
                            std::cout << std::endl;
                            break;
                        }
                    }
                }
            } else {
                if (expl.loops().find(ldes) == expl.loops().end()) {
                    std::cerr << "LOOPDESCRIPTOR record not found" << std::endl;
                } else {
                    auto loops = expl.loops().at(ldes);
                    for (auto & lp : loops) {
                        if (lp.number == loopNum) {
                            std::cout << "LOOPOPTM " << lp.toString() << std::endl;
                            std::cout << "COORDINATES " << lp.coords.size() << std::endl;
                            for (const XYZ & xyz : lp.coords) {
                                std::cout << xyz.x_ << " " << xyz.y_ << " " << xyz.z_ << std::endl;
                            }
                            std::cout << "ENERGIES";
                            for (double d : lp.energies) {
                                std::cout << " " << d;
                            }
                            std::cout << std::endl;
                            break;
                        }
                    }
                }
            }
        } else if ("show" == args[0] && args.size() >1 && "eneterms" == args[1]) {
            std::cout << "Etot Ebond Eang Eimpdih Ephipsi Escconf Esteric "
                    "Escpacking Elocalstructure Elocalhb Esitepairs Estructrest "
                    "Ergrestraint Essrestraint" << std::endl;
        } else if ("save" == args[0]) {
            if (args.size() == 7) {
                LoopDescriptor ldes;
                ldes.loopLoc.idxChain = std::stoi(args[1]);
                ldes.loopLoc.idxNTerminal = std::stoi(args[2]);
                ldes.loopLoc.idxCTerminal = std::stoi(args[3]);
                ldes.newLen = std::stoi(args[4]);
                int loopNum = std::stoi(args[5]);
                if (expl.loops().find(ldes) == expl.loops().end()) {
                    std::cerr << "LOOPDESCRIPTOR record not found" << std::endl;
                } else {
                    auto loops = expl.loops().at(ldes);
                    std::string & pdbfile = args[6];
                    for (auto & lp : loops) {
                        if (lp.number == loopNum) {
                            Protein pr = expl.restorePDBFromLoopOptm(lp);
                            std::ofstream ofs(pdbfile);
                            writetopdb(pr, ofs);
                            std::cout << "LOOPOPTM " << lp.toString()
                                    << " saved to " << pdbfile << std::endl;
                            break;
                        }
                    }
                }
            }
        } else if ("rmsd" == args[0]) {
            if (args.size() == 7) {
                LoopDescriptor ldes;
                ldes.loopLoc.idxChain = std::stoi(args[1]);
                ldes.loopLoc.idxNTerminal = std::stoi(args[2]);
                ldes.loopLoc.idxCTerminal = std::stoi(args[3]);
                ldes.newLen = std::stoi(args[4]);
                int loopNum1 = std::stoi(args[5]);
                int loopNum2 = std::stoi(args[6]);
                LoopOptm* loop1;
                LoopOptm* loop2;
                bool foundLoop1{false};
                bool foundLoop2{false};
                if (expl.loops().find(ldes) == expl.loops().end()) {
                    std::cerr << "LOOPDESCRIPTOR record not found" << std::endl;
                } else {
                    auto loops = expl.loops().at(ldes);
                    for (auto & lp : loops) {
                        if (lp.number == loopNum1) {
                            loop1 = &lp;
                            foundLoop1 = true;
                        } else if (lp.number == loopNum2) {
                            loop2 = &lp;
                            foundLoop2 = true;
                        }
                        if (foundLoop1 && foundLoop2) {
                            break;
                        }
                    }
                    if (foundLoop1 && foundLoop2) {
                        std::cout << "RMSD: " << rmsdBetweenLoops(loop1->coords, loop2->coords) << std::endl;
                    } else {
                        std::cerr << "Loop not found" << std::endl;
                    }
                }
            }
        } else if ("merge" == args[0]) {
            if (args.size() < 2) {
                std::cerr << "Missing another data to merge" << std::endl;
            } else {
                for (int i = 1; i < args.size(); i++) {
                    const std::string & mfile = args[i];
                    LoopExplorationIO io;
                    bool rs = io.read(mfile);
                    if (!rs) {
                        continue;
                    }
                    LoopExploration tmpLE = io.exploration();
                    int numLD = 0; // number of LoopDescriptors
                    int numLoops = 0; // number of loops
                    for (auto & dloops : tmpLE.loops()) {
                        numLD++;
                        int numLoopsOne = 0;
                        for (auto & loop : dloops.second) {
                            expl.addLoopOptm(loop);
                            numLoops++;
                            numLoopsOne++;
                        }
                        std::cout << "... LoopDescriptor " << dloops.first.toString()
                                << " has " << numLoopsOne << " loops" << std::endl;
                    }
                    std::cout << "Merged " << numLD << " LoopDescriptors containing "
                            << numLoops << " loops from file " << mfile << std::endl;
                }
                expl.updateLoopLocationIndex();
            }
        } else if ("stat" == args[0]) {
            std::vector<int> loopIndices;
            if (args.size() < 2 || "all" == args[1]) {
                for (int i = 0; i < expl.loopLocationIndex().size(); i++) {
                    loopIndices.push_back(i);
                }
            } else {
                for (int i = 1; i < args.size(); i++) {
                    loopIndices.push_back(std::stoi(args[i]));
                }
            }
            auto & loopAreas = expl.loopAreas();
            if (loopAreas.empty()) {
                expl.statLoopAreas();
            }
            for (int li : loopIndices) {
                std::cout << "LoopArea " << (li+1) << ":" << std::endl;
                auto & area = expl.loopAreas()[li].second;
                for (const auto & lc : area) {
                    const auto & loops = lc.second;
                    for (const auto & lo : loops) {
                        std::cout << lo.second.toString() << " " << lo.first << std::endl;
                    }
                }
                std::cout << "--------------------" << std::endl;
            }
        } else if ("refine" == args[0]) {
            int maxLoopsEachArea = 10;
            int shorterLenExtend = 1;
            double resRMSDCutoff = 0.1;
            if (args.size() == 2) {
                maxLoopsEachArea = std::stoi(args[1]);
            } else if (args.size() == 3) {
                shorterLenExtend = std::stoi(args[2]);
            } else if (args.size() == 4) {
                resRMSDCutoff = std::stod(args[3]);
            }
            expl.refineLoopAreas(maxLoopsEachArea, shorterLenExtend, resRMSDCutoff);
            int numAreas = expl.refinedLoopAreas().size();
            for (int li = 0; li < numAreas; li++) {
                std::cout << "LoopArea " << (li+1) << ":" << std::endl;
                auto & area = expl.refinedLoopAreas()[li].second;
                for (const auto & lc : area) {
                    const auto & loops = lc.second;
                    for (const auto & lo : loops) {
                        std::cout << lo.second.toString() << " " << lo.first << std::endl;
                    }
                }
                std::cout << "--------------------" << std::endl;
            }
        } else if ("assemble" == args[0] && args.size() > 2) {
            if (expl.refinedLoopAreas().empty()) {
                std::cerr << "Need to <refine> loops first" << std::endl;
                std::cout << "> ";
                continue;
            }
            int numConfs = std::stoi(args[1]);
            std::string outputDir = args[2];
            if (outputDir.at(outputDir.size()-1) != '/') {
                outputDir.push_back('/');
            }
            std::string pdbPrefix = "Loop-";
            if (args.size() > 3) {
                pdbPrefix = args[3];
            }
            std::ofstream ofs;
            /// Test availability of output directory
            std::string testfile = outputDir + "test.tmp";
            ofs.open(testfile);
            if (!ofs.is_open()) {
                std::cerr << "Cannot write to output directory" << std::endl;
                std::cout << "> ";
                continue;
            }
            ofs.close();
            std::remove(testfile.c_str());
            ///---
            auto eConfs = expl.assembleLoops(numConfs);
            int tmpNum = numConfs;
            int bitNum = 0;
            while (tmpNum > 0) {
                tmpNum = tmpNum / 10;
                bitNum++;
            }
            char buf[bitNum+1];
            std::string format = "%0" + std::to_string(bitNum) + "d";
            int loopNum = 0;
            for (const auto & eConf : eConfs) {
                loopNum++;
                const auto & pseudoE = eConf.first;
                const auto & conf = eConf.second;
                snprintf(buf, bitNum+1, format.c_str(), loopNum);
                std::string pdbNum(buf);
                std::string pdbName = pdbPrefix + pdbNum;
                ofs.open(outputDir + pdbName + ".pdb");
                std::string pdbHeader = "HEADER ";
                pdbHeader += std::to_string(pseudoE.at(0)) + " ";
                for (int i = 1; i < pseudoE.size(); i++) {
                    pdbHeader += std::to_string((int)pseudoE.at(i)) + " ";
                }
                pdbHeader += pdbName;
                ofs << pdbHeader << std::endl;
                writetopdb(conf, ofs);
                ofs.close();
                std::cout << "Saved conf " << pdbName << " to PDB file" << std::endl;
            }
        }
        else {
            std::cerr << "Unkown arguments" << std::endl;
        }
        std::cout << "> ";
    }
    return 0;
}

void assemble(int maxLoopsEachArea, int shorterLenExtend, double resRMSDCutoff,
              int numConfs, std::string & pdbPrefix, std::string & outputDir,
              std::vector<std::string> & loopOptmFiles) {
    LoopExploration expl = readExploration(loopOptmFiles[0]);
    for (int i = 1; i< loopOptmFiles.size(); i++) {
        const std::string & mfile = loopOptmFiles[i];
        LoopExplorationIO io;
        bool rs = io.read(mfile);
        if (!rs) {
            continue;
        }
        LoopExploration tmpLE = io.exploration();
        int numLD = 0; // number of LoopDescriptors
        int numLoops = 0; // number of loops
        for (auto & dloops : tmpLE.loops()) {
            numLD++;
            int numLoopsOne = 0;
            for (auto & loop : dloops.second) {
                expl.addLoopOptm(loop);
                numLoops++;
                numLoopsOne++;
            }
        }
        expl.updateLoopLocationIndex();
    }
    expl.refineLoopAreas(maxLoopsEachArea, shorterLenExtend, resRMSDCutoff);

    if (outputDir.at(outputDir.size()-1) != '/') {
        outputDir.push_back('/');
    }
    std::ofstream ofs;
    /// Test availability of output directory
    std::string testfile = outputDir + "test.tmp";
    ofs.open(testfile);
    if (!ofs.is_open()) {
        std::cerr << "Cannot write to output directory" << std::endl;
        return;
    }
    ofs.close();
    std::remove(testfile.c_str());
    ///---
    auto eConfs = expl.assembleLoops(numConfs);
    int tmpNum = numConfs;
    int bitNum = 0;
    while (tmpNum > 0) {
        tmpNum = tmpNum / 10;
        bitNum++;
    }
    char buf[bitNum+1];
    std::string format = "%0" + std::to_string(bitNum) + "d";
    int loopNum = 0;
    for (const auto & eConf : eConfs) {
        loopNum++;
        const auto & pseudoE = eConf.first;
        const auto & conf = eConf.second;
        snprintf(buf, bitNum+1, format.c_str(), loopNum);
        std::string pdbNum(buf);
        std::string pdbName = pdbPrefix + pdbNum;
        ofs.open(outputDir + pdbName + ".pdb");
        std::string pdbHeader = "HEADER ";
        pdbHeader += std::to_string(pseudoE.at(0)) + " ";
        for (int i = 1; i < pseudoE.size(); i++) {
            pdbHeader += std::to_string((int)pseudoE.at(i)) + " ";
        }
        pdbHeader += pdbName;
        ofs << pdbHeader << std::endl;
        writetopdb(conf, ofs);
        ofs.close();
        std::cout << "Saved conf " << pdbName << " to PDB file" << std::endl;
    }

}

void exportLoopConfs(int loopAreaNumber, int numConfs, double resRMSDCutoff,
                     std::string pdbPrefix, std::string outputDir,
                     std::string loopOptmFile) {
    LoopExploration expl = readExploration(loopOptmFile);
    expl.refineLoopAreas(numConfs, 0, resRMSDCutoff);
    auto loopArea = expl.refinedLoopAreas().at(loopAreaNumber-1).second;
    auto eLoops = loopArea.begin()->second;
    int actualNumConfs = numConfs > eLoops.size() ? eLoops.size() : numConfs;

    if (outputDir.at(outputDir.size()-1) != '/') {
        outputDir.push_back('/');
    }
    std::ofstream ofs;
    /// Test availability of output directory
    std::string testfile = outputDir + "test.tmp";
    ofs.open(testfile);
    if (!ofs.is_open()) {
        std::cerr << "Cannot write to output directory" << std::endl;
        return;
    }
    ofs.close();
    std::remove(testfile.c_str());
    ///---
    int tmpNum = numConfs;
    int bitNum = 0;
    while (tmpNum > 0) {
        tmpNum = tmpNum / 10;
        bitNum++;
    }
    char buf[bitNum+1];
    std::string format = "%0" + std::to_string(bitNum) + "d";
    int nc = 0;
    for (const auto el : eLoops) {
        nc++;
        double uniE = el.first;
        const auto lp = el.second;
        int loopLen =lp.loopDesp.newLen;
        std::printf("%.3f", uniE);
        std::cout << " " << loopLen;
        Protein pr = expl.restorePDBFromLoopOptm(lp);
        snprintf(buf, bitNum+1, format.c_str(), nc);
        std::string pdbNum(buf);
        std::string pdbName = pdbPrefix + pdbNum;
        std::cout << " " << pdbName << std::endl;
        ofs.open(outputDir + pdbName + ".pdb");
        ofs << "HEADER " << uniE << " "
            << lp.loopDesp.loopLoc.idxChain << " "
            << lp.loopDesp.loopLoc.idxNTerminal << " "
            << lp.loopDesp.loopLoc.idxNTerminal + loopLen - 1 << " "
            << loopLen << " "
            << pdbName << std::endl;
        writetopdb(pr, ofs);
        ofs.close();

        if (nc == actualNumConfs) {
            break;
        }
    }
}

void printUsage(const std::string & selfname)
{
    std::cout << "Usage:" << std::endl;
    std::cout << "$ " << selfname << " <LoopExplorationFile>"
            << std::endl;
}

int main(int argc, char **argv) {
    std::string selfname = argv[0];
    if (argc < 2) {
        std::cerr << "Wrong number of arguments." << std::endl;
        printUsage(selfname);
        return 1;
    }

    if (argc == 2) {
        std::string combofile(argv[1]);
        LoopExploration expl = readExploration(combofile);
        serve(expl);
    } else if (argc >= 3) {
        if (std::string(argv[1]) == "stat") {
            std::string combofile(argv[2]);
            LoopExploration expl = readExploration(combofile);
            stat(expl);
        } else if (std::string(argv[1]) == "assemble") {
            int maxLoopsEachArea = std::stoi(std::string(argv[2]));
            int shorterLenExtend = std::stoi(std::string(argv[3]));
            double resRMSDCutoff = std::stod(std::string(argv[4]));
            int numConfs = std::stoi(std::string(argv[5]));
            std::string pdbPrefix(argv[6]);
            std::string outputDir(argv[7]);
            std::vector<std::string> loopOptmFiles;
            for (int i = 8; i < argc; i++) {
                loopOptmFiles.push_back(std::string(argv[i]));
            }
            assemble(maxLoopsEachArea, shorterLenExtend, resRMSDCutoff, numConfs,
                     pdbPrefix, outputDir, loopOptmFiles);
        } else if (std::string(argv[1]) == "export") {
            int loopAreaNumber = std::stoi(std::string(argv[2]));
            int numConfs = std::stoi(std::string(argv[3]));
            double resRMSDCutoff = std::stod(std::string(argv[4]));
            std::string pdbPrefix(argv[5]);
            std::string outputDir(argv[6]);
            std::string loopOptmFile(argv[7]);
            exportLoopConfs(loopAreaNumber, numConfs, resRMSDCutoff, pdbPrefix, outputDir, loopOptmFile);
        }
    }

    return 0;
}

