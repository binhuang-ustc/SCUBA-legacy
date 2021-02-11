/*
 * activeselections.h
 *
 *  Created on: 2018年7月9日
 *      Author: hyliu
 */

#ifndef ACTIVESELECTIONS_H_
#define ACTIVESELECTIONS_H_
#include "sd/bsinchain.h"
#include <map>
#include <string>
namespace NSPsd{
class ForceField;
class ConformerCode;
//parse a string to obtain selected positions in chains
//for example, for an input string "chain0 5, 7-10, chain2 3,5,7,9-10",
// the resulting map will contain:
//{{0,{5,7,8,9,10}},{2,{3,5,7,9,10}}}
std::map<int,std::vector<int>> chainpositionselections(const std::string & str);
class ActiveSelections{
public:
	ActiveSelections(){;}
	ActiveSelections(const ForceField *ff,const std::string & allfixed,
			const std::string &mcfixed){
		init(ff,allfixed,mcfixed);
	}
	void init(const ForceField *ff,const std::string &allfixed,
			const std::string &mcfixed);
	void calcinitcodes(const std::vector<double> &crd);
	void calccodes(const std::vector<double> &crd);
	const std::vector<std::vector<PhiPsiCodes>> & phipsicodes() const {
		return phipsicodes_;
	}
	const std::vector<std::vector<SSCode>> & sscodes()const {
		return sscodes_;
	}
	const std::vector<std::vector<ConformerCode>> &conformercodes()const {
		return conformercodes_;
	}
	struct PosiStatus{
		enum {ACTIVE,SCACTIVE,NOACTIVE};
		int state{ACTIVE};
	};

	bool phipsiactive(int chain, int posi) const{
		return phipsiactive_[chain][posi];
	}
	bool localstructactive(int chain,int posi) const {
		return localstructactive_[chain][posi];
	}
	const std::vector<bool> & lsactive(int chain) const {
		return localstructactive_[chain];
	}
	bool rotameractive(int chain,int posi) const {
		return rotameractive_[chain][posi];
	}
	bool atomactive(int i) const {
		return atomactive_[i];
	}
	bool sitepairactive(int c1,int p1,int c2,int p2) const{
		return (posistatus_[c1][p1].state ==PosiStatus::ACTIVE) ||
				(posistatus_[c2][p2].state==PosiStatus::ACTIVE);
	}
	bool sidechainpairactive(int c1,int p1,int c2,int p2) const{
		return (posistatus_[c1][p1].state <= PosiStatus::SCACTIVE)||
				(posistatus_[c2][p2].state <= PosiStatus::SCACTIVE);
	}
	std::vector<bool> atomfixed() const{
		std::vector<bool> res;
		for(auto b:atomactive_) res.push_back(!b);
		return res;
	}
private:
	std::vector<std::vector<PosiStatus>> posistatus_;
	std::vector<std::vector<bool>> phipsiactive_;
	std::vector<std::vector<bool>> ssactive_;
	std::vector<std::vector<bool>> localstructactive_;
	std::vector<std::vector<bool>> rotameractive_;
	std::vector<bool> atomactive_;
	std::vector<std::vector<PhiPsiCodes>> phipsicodes_;
	std::vector<std::vector<SSCode>> sscodes_;
	std::vector<std::vector<ConformerCode>> conformercodes_;
	const ForceField *ff_{nullptr};
	bool doinitcalc_{true};
};

}



#endif /* ACTIVESELECTIONS_H_ */
