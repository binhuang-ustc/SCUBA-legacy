/*
 * proteincrdplus.h
 *
 *  Created on: 2018年1月19日
 *      Author: hyliu
 */

#ifndef FULLSITE_STRUCTPLUS_H_
#define FULLSITE_STRUCTPLUS_H_
#include "sd/bsinchain.h"
#include "fullsite/fullsite.h"
#include <set>
namespace NSPproteinrep{

class StructPlus{
public:
	struct Position : public std::pair<int,int> {
		int & chainid() {return first;}
		const int & chainid() const{return first;}
		int & posi(){return second;}
		const int &posi() const {return second;}
		Position (int cid=-1, int posi=-1): std::pair<int,int>(cid,posi){;}
		bool valid() const {return first>=0 && second >=0;}
	};
	struct Element{
		Position posi;
		int length;
	};
	void init(const std::vector<std::vector<FullSite>> &sites,bool doextra);
	StructPlus(const std::vector<std::vector<FullSite>> &sites,bool doextra=false){
		init(sites,doextra);
	}
	StructPlus(const std::string &pdbfilename,bool doextra=false){
		mysites_= std::shared_ptr<std::vector<std::vector<FullSite>>>(
				new std::vector<std::vector<FullSite>>(readfullsitesfrompdb(pdbfilename)));
//		*mysites_=;
		init(*mysites_,doextra);
	}
	void computeextra() {
		determinesstypes();
		computehbs();
		setssregions();
		determinestrandneighbors();
	}
	bool insamesse(const Position &p1,const Position &p2) const {
		int e1=elementid(p1);
		int e2=elementid(p2);
		if(e1 >=0 &&e1==e2) return true;
		return false;
	}
	int elementid(const Position & p) const {
		int id=-1;
		for(auto &e:elements_){
			++id;
			if(p.first != e.posi.first) continue;
			if(p.second >=e.posi.second && p.second <e.posi.second+e.length)return id;
		}
		return -1;
	}
	char sstype(const Position &p) const {
		int eleid=elementid(p);
		if(eleid>=0){
			return elementtypes_[eleid];
		} else {
			return 'C';
		}
	}
	bool hbonded(const Position &posi) const {
		return  hbdonors_.at(posi.chainid()).at(posi.posi()).valid() ||
				hbacceptors_.at(posi.chainid()).at(posi.posi()).valid();
	}
	bool inneighborstrands(const Position &p1, const  Position &p2) const {
		int e1=elementid(p1);
		int e2=elementid(p2);
		if(e1>e2){
			int tmp=e1;
			e1=e2;
			e2=tmp;
		}
		if(neighborstrands_.find(std::make_pair(e1,e2)) != neighborstrands_.end()) return true;
		return false;
	}
	/**
	 * returns 0if no HB found
	 * 1 p1 N-H..O=C p2
	 * 2 p1 C=O..H-N P2
	 * 3 both
	 */
	int mainchainhbond(const Position &p1,const Position &p2);
	void printhbonds() const{
		std::cout<<"Backbone hbonds: "<<std::endl;
		for(int cid=0;cid<hbdonors_.size();++cid){
			for(int p=0;p<hbdonors_[cid].size();++p){
				Position p1(cid,p);
				std::cout <<posistring(p1)<<":";
				if(hbdonors_[cid][p].valid()) std::cout<<posistring(hbdonors_[cid][p])<<":";
				if(hbacceptors_[cid][p].valid()) std::cout<<posistring(hbacceptors_[cid][p]);
				std::cout<<std::endl;
			}
		}
	}
	void printelements()const {
		std::cout<<"SS elements(TYPE:posi:length)"<<std::endl;
		for(int e=0;e<elements_.size();++e){
			std::cout<<e<<" "<<elementtypes_[e]<<":"<<posistring(elements_[e].posi)
					<<":"<<elements_[e].length<<std::endl;
		}
		std::cout<<"Strand pairs: elementid1:elementid2";
		for(auto &p:neighborstrands_){
			std::cout<<" "<<p.first<<":"<<p.second;
		}
		std::cout<<std::endl;
	}
	std::string posistring(const Position &posi) const {
		char chainid=sites_->at(posi.chainid()).at(posi.posi()).chainid();
		int resid=sites_->at(posi.chainid()).at(posi.posi()).resid();
		return std::to_string(resid)+std::string(1,chainid);
	}
	int nchains() const {return bsinchains_.size();}
	int length(int chainid) const {return bsinchains_.at(chainid).size();}
private:
	std::shared_ptr<std::vector<std::vector<FullSite>>> mysites_;
	const std::vector<std::vector<FullSite>> *sites_;
	std::vector<std::vector<NSPsd::BSInChain>> bsinchains_;
	std::vector<double> crd_;
	std::vector<std::vector<std::vector<double>>> sstypes_;
	std::vector<std::vector<Position>> hbdonors_;
	std::vector<std::vector<Position>> hbacceptors_;
	std::vector<Element> elements_;
	std::vector<char> elementtypes_;
	std::set<std::pair<int,int>> neighborstrands_;
	void determinesstypes();
	void setssregions();
	void determinestrandneighbors();
	void computehbs();
};
void ssregions(const std::vector<std::vector<double>> &sstypes,
		const std::vector<bool> & strandhbonded,
		std::vector<std::pair<int, int>> &helixregions,
		std::vector<std::pair<int, int>>& strandregions);
}
#endif /* FULLSITE_STRUCTPLUS_H_ */
