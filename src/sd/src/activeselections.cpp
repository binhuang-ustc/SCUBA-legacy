/*
 * activeselections.cpp
 *
 *  Created on: 2018年7月9日
 *      Author: hyliu
 */
#include "sd/activeselections.h"
#include "sd/forcefield.h"
#include "sd/nnterm.h"
#include "sd/sidechainff.h"
#include "dataio/inputlines.h"
using namespace NSPsd;
void ActiveSelections::calcinitcodes(const std::vector<double> & crd) {
	const std::vector<std::vector<BSInChain>> &bsinchains = ff_->bsinchains();
	int nchains = bsinchains.size();
	phipsicodes_.clear();
	for (int i = 0; i < nchains; ++i) {
		phipsicodes_.push_back(makephipsicodes(crd, bsinchains[i]));
	}
	sscodes_.clear();
	for (int i = 0; i < nchains; ++i) {
		sscodes_.push_back(estimatess(phipsicodes_[i]));
	}
	conformercodes_.clear();
	const std::vector<std::vector<SCInChain>> &scinchains = ff_->scinchains();
    for (int c = 0; c < nchains; ++c) {
            conformercodes_.push_back(
                            makeconformercodes(crd, scinchains[c], phipsicodes_[c]));
    }

}
void ActiveSelections::calccodes(const std::vector<double> &crd) {
	if (doinitcalc_) {
		calcinitcodes(crd);
		doinitcalc_ = false;
		return;
	} else {
		const std::vector<std::vector<BSInChain>> &bsinchains = ff_->bsinchains();
		int nchains = bsinchains.size();
		for (int ic = 0; ic < nchains; ++ic) {
			int nsite = bsinchains[ic].size();
#pragma omp parallel for schedule(dynamic,1)
			for (int i = 0; i < nsite; ++i) {
				if (!phipsiactive_[ic][i])
					continue;
				auto it = bsinchains[ic].begin() + i;
				phipsicodes_[ic][i] = PhiPsiCodes(crd, it,
						it == bsinchains[ic].begin(),
						it
								== bsinchains[ic].begin()
										+ bsinchains[ic].size() - 1);
			}
#pragma omp parallel for schedule(dynamic,1)
			for (int i = 0; i < nsite; ++i) {
				if (!ssactive_[ic][i])
					continue;
				auto & sc = sscodes_[ic].at(i);
				sc.p3 = NN_SSTerm().probabilities(phipsicodes_[ic], i,
						&(sc.dp3dx));
				sc.ssid = NN_SSTerm::sstype(sc.p3);
			}
			const std::vector<SCInChain> & scinchains=ff_->scinchains().at(ic);
#pragma omp parallel for schedule(dynamic,1)
			for(int i=0;i<nsite;++i){
				if(!rotameractive_[ic][i]) continue;
				if(scinchains[i].kaiatoms.empty()) continue;
				ConformerCode &cc=conformercodes_[ic][i];
				cc.restype=scinchains[i].restype;
				cc.phipsicodes=&(phipsicodes_[ic][i]);
				cc.sidechaintorsions.clear();
				cc.dtorsioncodesdx.clear();
				cc.torsioncodes.clear();
				for(auto & ijkl:scinchains[i].kaiatoms){
					std::vector<NSPgeometry::XYZ> dtdx;
					double kai=NSPgeometry::torsion(getxyz(crd,ijkl[0]),
							getxyz(crd,ijkl[1]),getxyz(crd,ijkl[2]),getxyz(crd,ijkl[3]),&dtdx);
					cc.sidechaintorsions.push_back(kai);
					std::vector<DvDxi> dadx;
					dadx.push_back(std::make_pair(ijkl[0],dtdx[0]));
					dadx.push_back(std::make_pair(ijkl[1],dtdx[1]));
					dadx.push_back(std::make_pair(ijkl[2],dtdx[2]));
					dadx.push_back(std::make_pair(ijkl[3],dtdx[3]));
					std::vector<std::vector<DvDxi>> dcdx;
					cc.dtorsioncodesdx.push_back(std::vector<std::vector<DvDxi>>());
					cc.torsioncodes.push_back(
							ConformerCode::gettorsioncodes(kai,dadx,&(cc.dtorsioncodesdx.back())));
				}
			}
		}
	}
}
void ActiveSelections::init(const ForceField *ff, const std::string & allfixed,
		const std::string &mcfixed) {
	ff_ = ff;
	const std::vector<std::vector<BSInChain>> &bsinchains = ff_->bsinchains();
	int nchains = bsinchains.size();
	std::map<int, std::vector<int>> allfixedcp = chainpositionselections(
			allfixed);
	std::map<int, std::vector<int>> mcfixedcp = chainpositionselections(
			mcfixed);
	posistatus_.assign(nchains, std::vector<PosiStatus>());
	//initially, all positions in all chains active
	for (int ic = 0; ic < nchains; ++ic) {
		int cl = bsinchains[ic].size();
		posistatus_[ic].assign(cl, PosiStatus());
	}
	//set noactive for fixed positions
	for (auto & af : allfixedcp) {
		int c = af.first;
		for (auto p : af.second)
			posistatus_[c][p].state = PosiStatus::NOACTIVE;
	}
	//set scactive for mainchain-fixed positions
	for (auto &mcf : mcfixedcp) {
		int c = mcf.first;
		for (auto p : mcf.second)
			if (posistatus_[c][p].state == PosiStatus::ACTIVE)
				posistatus_[c][p].state = PosiStatus::SCACTIVE;
	}
	//phipsi ai position i active if any of residues i, i-1 and i+1 is main chain active
	phipsiactive_.assign(nchains, std::vector<bool>());
	for (int ic = 0; ic < nchains; ++ic) {
		for (int i = 0; i < bsinchains[ic].size(); ++i) {
			bool a = false;
			if (posistatus_[ic][i].state == PosiStatus::ACTIVE)
				a = true;
			if (i > 0)
				if (posistatus_[ic][i - 1].state == PosiStatus::ACTIVE)
					a = true;
			if (i < bsinchains[ic].size() - 1)
				if (posistatus_[ic][i + 1].state == PosiStatus::ACTIVE)
					a = true;
			phipsiactive_[ic].push_back(a);
		}
	}
	//localstruct active if any of residues from i-2 to i+2 active
	localstructactive_.assign(nchains, std::vector<bool>());
	for (int ic = 0; ic < nchains; ++ic) {
		for (int i = 0; i < bsinchains[ic].size(); ++i) {
			bool a = false;
			int left = (i - 2) < 0 ? 0 : (i - 2);
			int right =
					(i + 3) > bsinchains[ic].size() ?
							bsinchains[ic].size() : (i + 3);
			for (int p = left; p < right; ++p)
				if (posistatus_[ic][p].state == PosiStatus::ACTIVE)
					a = true;
			localstructactive_[ic].push_back(a);
		}
	}
	rotameractive_.assign(nchains, std::vector<bool>());
	for (int ic = 0; ic < nchains; ++ic) {
		for (int i = 0; i < bsinchains[ic].size(); ++i) {
			bool a = false;
			if (phipsiactive_[ic][i]
					|| posistatus_[ic][i].state <= PosiStatus::SCACTIVE)
				a = true;
			rotameractive_[ic].push_back(a);
		}
	}
	//ssative if any of residues i-3 to i+3 active
	ssactive_.assign(nchains, std::vector<bool>());
	for (int ic = 0; ic < nchains; ++ic) {
		for (int i = 0; i < bsinchains[ic].size(); ++i) {
			bool a = false;
			int left = (i - 3) < 0 ? 0 : (i - 3);
			int right =
					(i + 4) > bsinchains[ic].size() ?
							bsinchains[ic].size() : (i + 4);
			for (int p = left; p < right; ++p)
				if (posistatus_[ic][p].state == PosiStatus::ACTIVE)
					a = true;
			ssactive_[ic].push_back(a);
		}
	}
	const std::vector<std::vector<SCInChain>> &scinchains = ff_->scinchains();
	atomactive_.assign(ff_->natoms(), true);
	for (int ic = 0; ic < nchains; ++ic) {
		for (int p = 0; p < bsinchains[ic].size(); ++p) {
			if (PosiStatus::ACTIVE == posistatus_[ic][p].state)
				continue;
			//inactivate mainchain atoms
			if (PosiStatus::SCACTIVE <= posistatus_[ic][p].state) {
				for (auto &a : bsinchains[ic][p].atomids())
					atomactive_[a] = false;
			}
			//in addition, inactivate sidechain atoms
			if (PosiStatus::NOACTIVE <= posistatus_[ic][p].state) {
				for (int a = 0; a < scinchains[ic][p].nscatoms; ++a) {
					atomactive_[scinchains[ic][p].poffset + a] = false;
				}
			}
		}
	}
}
std::map<int, std::vector<int>> NSPsd::chainpositionselections(
		const std::string& str) {
	std::vector<std::string> words = NSPdataio::parseline(str,
			std::vector<int>());
	std::map<int, std::vector<int>> result;
	std::map<int, std::string> chainselections;
	int cid = -1;
	try {
		for (auto & w : words) {
			if (w.substr(0, 5) == "chain") {
				cid = std::stoi(w.substr(5));
				if (chainselections.find(cid) == chainselections.end()) {
					chainselections[cid] = std::string();
				}
			} else {
				assert(cid >= 0);
				chainselections.at(cid) += w + ",";
			}
		}
		for (auto &cs : chainselections) {
			result[cs.first] = NSPdataio::stringtoselection(cs.second);
		}
	} catch (std::exception &e) {
		std::cout << "Error procession selection string: " << str << std::endl;
		exit(1);
	}
	return result;
}

