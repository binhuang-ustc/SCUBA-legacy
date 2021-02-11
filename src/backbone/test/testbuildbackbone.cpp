/*
 * testbuildbackbone.cpp
 *
 *  Created on: 2018/10/18
 *      Author: zhaozian & hbfrank
 */

#include <string>
#include <vector> 
#include <iostream>
#include <fstream>
#include <sstream>

#include <set>
#include <utility>
#include <memory>

#include "backbone/backbonebuilder.h"
#include "dstl/randomengine.h"

using namespace NSPproteinrep;
using namespace NSPgeometry;

using namespace std;

typedef vector< vector<double> > matrix;


//skip the lines started with '#' in data file
void skipline(ifstream & f)
{
    string stringline;
    streampos position;
    do
    {
        position = f.tellg();
        getline(f, stringline);
    }while (stringline[0] == '#'|| stringline.size() == 0);
    f.seekg(position);
}

//read a matrix form file f to variable v
//The matrix reader will stop when it meets '#' or '\n' as the first char of the line
void readmatrix(ifstream & f, matrix & v)
{
    string stringline;
    getline(f, stringline);
    while(stringline[0] != '#' && stringline.size() != 0)
    {
        vector<double> doubleline;
        istringstream iss (stringline);
        double di;

        while (iss >> di)
            doubleline.push_back(di);
        v.push_back(doubleline);
        getline(f, stringline);
    }
}


int buildchain(matrix& StructInfo, matrix& Vdirection, string& PDBFileName)
{
    //Set backbones
    vector<double> structinfoi;
    int N = StructInfo.size();
    vector< vector<BackBoneSite> >BackBones(N);
    int i = 0;

    for (; i < N; ++i)
    {
        structinfoi = StructInfo[i];
        XYZ startpoint (structinfoi);
        int length = structinfoi[3];
        XYZ direction (Vdirection[structinfoi[4]-1]) ;
        bool forward = (bool)structinfoi[5];
        if (structinfoi[6] == 1)//strand ?
            BackBones[i]=BackBoneBuilder::buildstrandat(length, startpoint, direction, forward);
        else
            BackBones[i]=BackBoneBuilder::buildhelixat(length, startpoint, direction, forward);
    }
    cout << "build backbones for " << i << " times" << '\n';

    //Link backbones
    std::vector<BackBoneSite> chain;
    for(auto &s:BackBones[0]) chain.push_back(s);

    for (int i = 0; i < N-1; ++i)
    {
        int j = i+1;
        std::vector<std::shared_ptr<std::vector<BackBoneSite>>> linki;
        int ntry = 0;
        int length = (int)StructInfo[i][7];
        std::cout << "Try linking backbone " << i << " and " << j <<std::endl;
        do
        {
            linki = BackBoneBuilder::buildlinkers(length,BackBones[i].back(),BackBones[j][0],
                    std::vector<std::pair<int,int>> (),
                    std::vector<std::pair<int,int>>(),std::set<int>());
            ++ntry;
        }   while(linki.empty() && ntry < 2000);
        if (linki.empty())
        {
            std::cout << "Failed to link backbone " << i << " and " << j <<std::endl;
            return 1;
        }
        //save to chain
        for(auto &s:*(linki[0])) {
            s.resname = "GLY";
            chain.push_back(s);
        }
        for(auto &s:BackBones[j]) {
            s.resname = "GLY";
            chain.push_back(s);
        }
    }


    //Save pdb file
    std::ofstream ofs;
    ofs.open(PDBFileName.c_str());
    writeSitesToPDB(ofs,chain);
    ofs.close();

    return 0;
}

int main(int argc,char **argv)
{
    if (argc != 4) {
        std::cerr << "Usage: " << std::endl << "testbuildbackbone <sketchdeffile> <outputpdbfile> <randseed>" << std::endl;
        return 1;
    }

    ifstream sketchdef (std::string(argv[1]).c_str());
    string PDBFileName(argv[2]);
    int seed=std::stoi(std::string(argv[3]));

    NSPdstl::RandomEngine<>::getinstance().reseed(seed);

    skipline(sketchdef);
    skipline(sketchdef);
    matrix Vdirection;
    readmatrix(sketchdef, Vdirection);

    skipline(sketchdef);
    matrix StructInfo;
    readmatrix(sketchdef, StructInfo);
    sketchdef.close();

    int count = 0;
    while(true)
    {
        if  (buildchain(StructInfo,Vdirection,PDBFileName) == 0) {
            break;
        }
        if (count++ > 1000) {
            break;
        }
    }
}
