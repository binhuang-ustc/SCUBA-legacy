/*
 * structplus.cpp
 *
 *  Created on: 2018年1月19日
 *      Author: hyliu
 */

#include "fullsite/structplus.h"
#define A2NM 0.1
using namespace NSPproteinrep;
using namespace NSPsd;
void StructPlus::init(const std::vector<std::vector<FullSite>> &sites,bool doextra){
	int atomid=0;
	sites_=&sites;
	for(auto &chain:sites){
		bsinchains_.push_back(std::vector<BSInChain>());
		auto &bicvec=bsinchains_.back();
		for(auto &s:chain){
			bicvec.push_back(BSInChain());
			auto &bic=bicvec.back();
			std::vector<NSPgeometry::XYZ> scrd=s.extractcrd();
			for( auto &r:scrd){
				crd_.push_back(r.x_*A2NM);
				crd_.push_back(r.y_*A2NM);
				crd_.push_back(r.z_*A2NM);

			}
			std::vector<std::string> atomnames=s.atomnames_all();
			for(auto &a:atomnames){
				if(a=="N") bic.nid=atomid;
				else if(a=="CA")bic.caid=atomid;
				else if(a=="C") bic.cid=atomid;
				else if(a=="O") bic.oid=atomid;
				if(s.hasatomcrd(a))++atomid;
			}
		}
	}
	if(doextra)computeextra();
}
void StructPlus::determinesstypes(){
	sstypes_.clear();
	for(auto &bsvec:bsinchains_){
		sstypes_.push_back(std::vector<std::vector<double>>());
		std::vector<NSPsd::PhiPsiCodes> phipsicodes=makephipsicodes(crd_,bsvec);
		std::vector<NSPsd::SSCode> sscodes=estimatess(phipsicodes);
		std::vector<std::vector<double>> &res=sstypes_.back();
		for(auto &ssc:sscodes) {
			res.push_back(ssc.p3);
		}
	}
}
void StructPlus::computehbs(){
	hbdonors_.clear();
	hbacceptors_.clear();
	for(int cid=0;cid<bsinchains_.size();++cid){
		hbdonors_.push_back(std::vector<Position>());
		hbacceptors_.push_back(std::vector<Position>());
		auto &  hbd=hbdonors_.back();
		auto & hba=hbacceptors_.back();
		for(int pid=0;pid<bsinchains_[cid].size();++pid){
			Position p1(cid,pid);
			Position phbd(-1,-1);
			Position phba(-1,-1);
			for( int cid2=0;cid2<bsinchains_.size();++cid2){
				for(int pid2=0;pid2<bsinchains_[cid2].size();++pid2){
					Position p2(cid2,pid2);
					if(p1==p2) continue;
					int hb=mainchainhbond(p1,p2);
					if(hb==2 ||hb==3) phbd=p2;
					if(hb==1|| hb==3) phba=p2;
					if(phbd.valid()&&phba.valid()) break;
				}
				if(phbd.valid()&&phba.valid()) break;
			}
			hbd.push_back(phbd);
			hba.push_back(phba);
		}
	}
}
void StructPlus::setssregions(){
	for(int cid=0;cid<bsinchains_.size();++cid){
		std::vector<std::pair<int, int>> helixregions;
		std::vector<std::pair<int, int>> strandregions;
		std::vector<bool> strandhbonded;
		for(int p=0;p<bsinchains_.at(cid).size();++p){
			if(sstypes_[cid][p].size()==3){
				if(sstypes_[cid][p][1]<0.2){
					strandhbonded.push_back(false);
					continue;
				}
			}
			Position phbd=hbdonors_[cid][p];
			bool shb=false;
			if(phbd.valid()){
				if(sstypes_[phbd.chainid()][phbd.posi()].size()!=3) shb=true;
				else shb=sstypes_[phbd.chainid()][phbd.posi()][1]>0.2;
			}
			Position phba=hbacceptors_[cid][p];
			if(phba.valid()) {
				if(sstypes_[phba.chainid()][phba.posi()].size()!=3) shb=true;
				else shb=shb || (sstypes_[phba.chainid()][phba.posi()][1]>0.2);
			}
			strandhbonded.push_back(shb);
		}
		ssregions(sstypes_[cid],strandhbonded,helixregions,strandregions);
		for(auto &h:helixregions){
			elements_.push_back(Element());
			elementtypes_.push_back('H');
			Element & ele=elements_.back();
			ele.posi.chainid()=cid;
			ele.posi.posi()=h.first;
			ele.length=h.second-h.first;
		}
		for(auto &s:strandregions){
			elements_.push_back(Element());
			elementtypes_.push_back('E');
			Element & ele=elements_.back();
			ele.posi.chainid()=cid;
			ele.posi.posi()=s.first;
			ele.length=s.second-s.first;
		}
	}
}
void StructPlus::determinestrandneighbors(){
	if(elements_.empty()) return;
	for(int eid=0;eid<elements_.size()-1;++eid){
		if(elementtypes_[eid] !='E') continue;
		for (int eid2=eid+1;eid2<elements_.size();++eid2){
			if(elementtypes_[eid] !='E') continue;
			int nhb=0;
			for(int i=0;i<elements_[eid].length;++i){
				Position p1;
				p1.chainid()=elements_[eid].posi.chainid();
				p1.posi()=elements_[eid].posi.posi()+i;
				for(int j=0;j<elements_[eid2].length;++j){
					Position p2;
					p2.chainid()=elements_[eid2].posi.chainid();
					p2.posi()=elements_[eid2].posi.posi()+j;
					int hb=mainchainhbond(p1,p2);
					nhb += (int) ((hb+1)/2);
				}
			}
			if(nhb>=2) neighborstrands_.insert(std::make_pair(eid,eid2));
		}
	}
}
int StructPlus::mainchainhbond(const Position &posi1,const Position &posi2){
	NSPgeometry::XYZ ncrd1=getxyz(crd_,bsinchains_[posi1.chainid()][posi1.posi()].nid);
	NSPgeometry::XYZ ccrd1=getxyz(crd_,bsinchains_[posi1.chainid()][posi1.posi()].cid);
	NSPgeometry::XYZ ocrd1=getxyz(crd_,bsinchains_[posi1.chainid()][posi1.posi()].oid);
	NSPgeometry::XYZ ncrd2=getxyz(crd_,bsinchains_[posi2.chainid()][posi2.posi()].nid);
	NSPgeometry::XYZ ccrd2=getxyz(crd_,bsinchains_[posi2.chainid()][posi2.posi()].cid);
	NSPgeometry::XYZ ocrd2=getxyz(crd_,bsinchains_[posi2.chainid()][posi2.posi()].oid);
	double rcut2=0.35*0.35;
	double acut=110.0*3.14159265/180.0;
	double r2_n1o2=(ncrd1-ocrd2).squarednorm();
	bool h12=false;
	if(r2_n1o2<=rcut2){
		double a=NSPgeometry::angle(ncrd1,ocrd2,ccrd2);
		if(a>acut)h12=true;;
	}
	bool h21=false;
	double r2_o1n2=(ocrd1-ncrd2).squarednorm();
	if(r2_o1n2<=rcut2){
		double  a=NSPgeometry::angle(ccrd1,ocrd1,ncrd2);
		if(a>acut) h21=true;
	}
	if(h12&&h21)return 3;
	else if(h12) return 1;
	else if(h21) return 2;
	else return 0;
}
void NSPproteinrep::ssregions(const std::vector<std::vector<double>> &sstypes,
		const std::vector<bool> &strandhbonded,
		std::vector<std::pair<int, int>> &helixregions,
		std::vector<std::pair<int, int>>& strandregions) {
	int hstart = -1;
	int hlen = 0;
	int hc=0;
	for (int i = 0; i < sstypes.size(); ++i) {
		if(sstypes[i].size() != 3) continue;
		if (sstypes[i][0] > 0.5) {
			if (hstart == -1) {
				hstart = i;
				hlen = 0;
			}
			hlen++;
			if(sstypes[i][0]>0.75) hc++;
		} else {
			if (hstart != -1) {
				if (hc >= 4) {
					int hend = hstart + hlen;
					helixregions.push_back(std::make_pair(hstart, hend));
				}
			}
			hstart = -1;
			hc=0;
		}
	}
	if(hstart !=-1 && hc>=4){
		helixregions.push_back(std::make_pair(hstart, hstart+hlen));
	}
	int estart = -1;
	int elen = 0;
	int ec=0;
	for (int i = 0; i < sstypes.size(); ++i) {
//		if(sstypes[i].size() !=3) continue;
		if (strandhbonded[i] ) {
			if (estart == -1) {
				estart = i;
				elen = 0;
			}
			elen++;
			ec++;
		} else {
			if (estart != -1) {
					if(!strandhbonded[i-1]){
//				if (ec >= 3) {
						if(elen>=3 && ec>=2){
							int eend = estart + elen;
							strandregions.push_back(std::make_pair(estart, eend));
						}
					} else {
						elen++;
						continue;
					}
			}
			estart = -1;
			ec=0;
		}
	}
	if(estart !=-1 && elen>=3 && ec>2){
		strandregions.push_back(std::make_pair(estart,estart+elen));
	}
/*	std::cout << " Helix Regions: ";
	for (auto &h : helixregions)
		std::cout << h.first << ":" << h.second - 1 << " ";
	std::cout << std::endl;
	std::cout << " Strand Regions: ";
	for (auto &e : strandregions)
		std::cout << e.first << ":" << e.second - 1 << " ";
	std::cout << std::endl;*/
	return;
}

