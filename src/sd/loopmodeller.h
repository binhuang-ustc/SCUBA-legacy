/*
 * loopmodeller.h
 *
 *  Created on: 2018年7月10日
 *      Author: hyliu
 */

#ifndef LOOPMODELLER_H_
#define LOOPMODELLER_H_
#include "fullsite/fullsite.h"
#include <memory>

namespace NSPsd{
// a position in possibly multi-chain protein
struct ChainPosi{
	int cid{0}; //chain id
	int posi{-1}; //position id
	ChainPosi(int p):posi(p){;}
	ChainPosi(int c,int p):cid(c),posi(p){;}
	//make it comparable
	friend bool operator<(const ChainPosi &cp1, const ChainPosi &cp2){
		if(cp1.cid<cp2.cid) return true;
		else if(cp1.cid==cp2.cid) return cp1.posi<cp2.posi;
		else return false;
	}
};

// Specifies a part of the basechains_, which is to be substituted by the rebuilt loop
class TargetLoop{
public:
	TargetLoop(){;}
	TargetLoop(const std::vector<std::vector<NSPproteinrep::FullSite>>&basec,
			int pos1,int pos2): basechains_(&basec),
		startcp_(pos1),endcp_(pos2){;}
	TargetLoop(const std::vector<std::vector<NSPproteinrep::FullSite>>&basec,
			int cid1,int pos1,int cid2,int pos2):basechains_(&basec),
		startcp_(cid1,pos1),endcp_(cid2,pos2){;}
	friend bool operator<(const TargetLoop &tl1, const TargetLoop &tl2){
			return tl1.startcp_<tl2.startcp_;
		}
	bool isstartposi(int c,int p){
		return (c==startcp_.cid) &&(p==startcp_.posi);
	}
	bool isendposi(int c,int p){
		return (c==endcp_.cid) &&(p==endcp_.posi);
	}
	bool contains(int c, int p){
		if(isstartposi(c,p)) return true;
		ChainPosi cp(c,p);
		return (startcp_<cp && cp<endcp_);
	}
	//extract the structures in the basechain
	std::vector<NSPproteinrep::FullSite> extractbaseconfiguration() const;//TODO
	//retrieve a generated loop, remove it from the reservior of new conformers_
	std::shared_ptr<std::vector<NSPproteinrep::FullSite>> popconformer(const
	std::vector<std::string> &newsequence=std::vector<std::string>()){
		std::shared_ptr<std::vector<NSPproteinrep::FullSite>> res=getconformer(newsequence);
		if(res) newconformers_.resize(newconformers_.size()-1);
		return res;
	}
	//retrieves a generated loop
	std::shared_ptr<std::vector<NSPproteinrep::FullSite>> getconformer(
			const std::vector<std::string> &newsequence=std::vector<std::string>()){
			bool sequencechanged{false};
			if(!newsequence.empty()){
				if(newsctypes_!=newsequence) sequencechanged=true;
				newlength_=newsequence.size();
				newsctypes_=newsequence;
			}
			assert(newlength_>0);
			if(sequencechanged || newconformers_.empty()){
				if(!buildnewconformers()) return nullptr;
			}
			return newconformers_.back();
		}
private:
	//* starting and ending positions of the target loops in the basechains
	ChainPosi startcp_,endcp_;
	//* base chains
	const std::vector<std::vector<NSPproteinrep::FullSite>> *basechains_{nullptr};
	//* intended length of the new loop
	int newlength_{-1};
	//* intended amino acid sequences of the new loop
	std::vector<std::string> newsctypes_;
	//* stores loop conformations generated form backbone builder
	std::vector<std::shared_ptr
		<std::vector<NSPproteinrep::FullSite>>> newconformers_;
	//calls backbonebuilder to generate random but properly closed loop structures
	//the generated structures are stored in newconformers_
	bool buildnewconformers(); //TODO
};

//stores the previously optimized conformations of all the target loops.
class ArchivedLoops{
public:
	ArchivedLoops(double rcut):rmsdcut_(rcut){;}
	//try to add another optimized conformer to the archive
	bool addmodel(std::shared_ptr<NSPproteinrep::FullSite> model,
			const std::vector<double> &energy);//TODO
private:
	double rmsdcut_;
	std::vector<std::shared_ptr<std::vector<NSPproteinrep::FullSite>>> conformers_;
	std::vector<std::vector<double>> energies_;
};

class LoopModeller{
public:
	LoopModeller(){;}
	LoopModeller(const NSPdataio::ParameterSet &pset);//TODO

	std::vector<std::vector<NSPproteinrep::FullSite>>  & basechains() {
		return basechains_;
	}
	const std::vector<std::vector<NSPproteinrep::FullSite>>  & basechains() const {
		return basechains_;
	}
	void run();//TODO
/*	std::vector<std::vector<NSPproteinrep::FullSite>> newinitchains();
	std::vector<std::string> genchaincontrolines();
	std::vector<std::string> quenchsdcontrolines();
	std::shared_ptr<std::vector<std::vector<NSPproteinrep::FullSite>>>
		extractmodelledloops(const std::vector<double> &crd);*/

private:
 std::vector<std::vector<NSPproteinrep::FullSite>> basechains_;
 std::shared_ptr<GenChain> genchain_;
 ArchivedLoops archivedloops_;
 std::vector<TargetLoop> loops_;
 NSPdataio::ParameterSet pset_;
 //call from constructor
void addtargetloop(int p1,int p2){
		loops_.push_back(TargetLoop(basechains_,p1,p2));
	}
void addtargetloop(int c1, int p1,int c2,int p2){
		loops_.push_back(TargetLoop(basechains_,c1,p1,c2,p2));
	}

};

}



#endif /* LOOPMODELLER_H_ */
