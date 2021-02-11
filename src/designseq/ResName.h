/*
 * ResName.h
 *
 *  Created on: 2017Äê10ÔÂ17ÈÕ
 *      Author: notxp
 */

#ifndef DESIGNSEQ_RESNAME_H_
#define DESIGNSEQ_RESNAME_H_

#include <string>
#include <vector>
#include <map>

namespace NSPdesignseq {

using namespace std;
class ResName {
private:
    string aaSeq;
    vector<string> triNameList;
    map<string,int> triToIntMap;
    map<char,int> sinToIntMap;

public:
    ResName();
    char intToSin(int i) const;
    int sinToInt(char sin) const;
    char triToSin(const string& tri) const;
    string sinToTri(char sin) const;
    int triToInt(const string& tri) const;
    string intToTri(int i) const;
    bool isStandardAminoAcid(const string& tri);
    bool isAminoAcid(const string& tri);
    int chiNum(const string& tri) const;
    virtual ~ResName();
};

} /* namespace NSPdesignseq */

#endif /* DESIGNSEQ_RESNAME_H_ */
