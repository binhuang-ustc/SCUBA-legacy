/*
 * forcefield.cpp
 *
 *  Created on: 2017年11月2日
 *      Author: hyliu
 */

#include "sd/forcefield.h"
#include "sd/genchain.h"
#include "sd/activeselections.h"
#ifdef CBPAIRENE
#include "sd/cbpairene.h"
#define DOCB
#endif
#include "sd/nnterm.h"
#include "sd/localhbterm.h"
#include "sd/backboneff.h"
#include "fullsite/structplus.h"
#include "sd/shadowterm.h"
#ifdef USE_CB
#define DOCB
#endif
#include <cassert>
#include <iomanip>
#include <sstream>
using namespace NSPsd;
using namespace NSPgeometry;
double StericAtom::DEDRMAX{500.0};
double ForceField::highescale{1.0};
void NSPsd::defineforcefieldcontrol(std::string name,
		const std::vector<std::string> &controllines) {
	std::map<std::string, double> doublepars { { "StericWeight", 1.0 }, {
			"PhiPsiWeight", 0.0 }, { "SideConfWeight", 0.0 },
			{ "LSWeight", 0.0 }, { "PairPackingWeight", 0.0 }, {"SCPackingWeight",0.0},
			{ "EAttraction",0.0 },
			{ "RAttractionOff", 10000.0 }, { "RAttractionSwitch",
					1000.0 }, { "LocalHBWeight", 0.0 }, { "ShadowWeight", 0.0 },
			{ "KresStrandPair", 50.0 } };
	std::map<std::string, std::vector<std::string>> stringvecpars { {
			"SubStructureRestrains", { } }, { "TotalRMSDRestrain", { } } };
	std::map<std::string, std::vector<double>> doublevecpars { { "RgRestrain", {
			100.0, 0.0 } }, { "HelixRestrains", { } },
			{ "StrandRestrains", { } }, { "DisRestrains", { } }, {
					"LoopRegionWeights", { } } };
	std::map<std::string, std::string> stringpars { { "CoreRegionBasePDB", "" } };
	std::map<std::string, std::vector<int>> intvecpars { {
			"StrandPairRestrains", { } }, { "ExtendRestrains", { } }, {
			"CoreRegion", { } },{"RMSDRestrainedRegion",{} }};
	std::map<std::string, int> intpars{{"SideChainMMSteric",0}};
	ForceFieldControls::initdefaultkeyvals(name, doublepars, stringpars,
			intpars, doublevecpars, stringvecpars, intvecpars);
	int nsuccess = ForceFieldControls::adjustvalues(name, controllines);
	if (nsuccess != controllines.size()) {
		exit(1);
	}
}
static std::string makessreststr(int chainid, int begin, int end, double p0,
		double kres) {
	std::ostringstream oss;
	oss << std::setw(2) << chainid << std::setw(5) << begin << std::setw(5)
			<< end << std::setw(5) << std::setiosflags(std::ios::right)
			<< std::setiosflags(std::ios::fixed) << std::setprecision(2) << p0
			<< std::setw(10) << std::setiosflags(std::ios::right)
			<< std::setiosflags(std::ios::fixed) << std::setprecision(2)
			<< kres;
	return oss.str();
}
static void ssregions(const std::vector<std::vector<double>> &sstypes,
		std::vector<std::pair<int, int>> &helixregions,
		std::vector<std::pair<int, int>>& strandregions) {
	int hstart = -1;
	int hlen = 0;
	int hc = 0;
	for (int i = 0; i < sstypes.size(); ++i) {
		if (sstypes[i].size() != 3)
			continue;
		if (sstypes[i][0] > 0.5) {
			if (hstart == -1) {
				hstart = i;
				hlen = 0;
			}
			hlen++;
			if (sstypes[i][0] > 0.75)
				hc++;
		} else {
			if (hstart != -1) {
				if (hc >= 4) {
					int hend = hstart + hlen;
					helixregions.push_back(std::make_pair(hstart, hend));
				}
			}
			hstart = -1;
			hc = 0;
		}
	}
	int estart = -1;
	int elen = 0;
	int ec = 0;
	for (int i = 0; i < sstypes.size(); ++i) {
		if (sstypes[i].size() != 3)
			continue;
		if (sstypes[i][1] > 0.5) {
			if (estart == -1) {
				estart = i;
				elen = 0;
			}
			elen++;
			if (sstypes[i][1] > 0.75)
				ec++;
		} else {
			if (estart != -1) {
				if (ec >= 3) {
					int eend = estart + elen;
					strandregions.push_back(std::make_pair(estart, eend));
				}
			}
			estart = -1;
			ec = 0;
		}
	}
	std::cout << " Helix Regions: ";
	for (auto &h : helixregions)
		std::cout << h.first << ":" << h.second - 1 << " ";
	std::cout << std::endl;
	std::cout << " Strand Regions: ";
	for (auto &e : strandregions)
		std::cout << e.first << ":" << e.second - 1 << " ";
	std::cout << std::endl;
	return;
}
std::vector<std::string> NSPsd::ssrestraincontrolfromchain(
		const std::vector<NSPproteinrep::BackBoneSite> &chain, int chainid,
		double p0, double kres) {
	std::vector<std::vector<double>> sstypes = estimatess(chain);
	std::vector<std::pair<int, int>> helixregions;
	std::vector<std::pair<int, int>> strandregions;
	ssregions(sstypes, helixregions, strandregions);
	std::vector<std::string> extraffcontrolls;
	if (!helixregions.empty()) {
		std::string line = "HelixRestrains= ";
		for (auto &h : helixregions) {
			line = line + makessreststr(chainid, h.first, h.second, p0, kres);
		}
		extraffcontrolls.push_back(line);
	}
	if (!strandregions.empty()) {
		std::string line = "StrandRestrains= ";
		for (auto &s : strandregions) {
			line = line + makessreststr(chainid, s.first, s.second, p0, kres);
		}
		extraffcontrolls.push_back(line);
	}
	return extraffcontrolls;
}

static double estimatergbound(
		const std::vector<NSPproteinrep::BackBoneSite> & chain) {
	std::vector<NSPgeometry::XYZ> crd;
	for (auto &c : chain)
		c.getcrd(crd);
	double rg = radiusgyr(crd);
	std::cout << "Estimated RG bound: 1.5 *" << rg << " = " << 1.5 * rg
			<< std::endl;
	return 1.5 * rg;
}

std::string NSPsd::rgrestraincontrolfromchain(
		const std::vector<NSPproteinrep::BackBoneSite> &chain, double krgres) {
	double rgbound = estimatergbound(chain);
	std::ostringstream oss;
	oss << "RgRestrain= " << rgbound << " " << krgres;
	return oss.str();
}
void EneComp::print(std::ostream &os) {
	os << "Energies of individual backbone sites: " << std::endl;
	for (int i = 0; i < ephipsi.size(); ++i) {
		os << i << " " << ephipsi[i] << " " << els[i] << " " << epacking[i]
				<< " " << elocalhb[i] << std::endl;
	}
	os << "Non-zero Pair-Pair packing energies:" << std::endl;
	for (int i = 0; i < esitepair.size(); ++i) {
		os << i;
		for (int j = 0; j < esitepair.size(); ++j) {
			if (esitepair[i][j] != 0)
				os << " (" << j << "," << esitepair[i][j] << ")";
		}
		std::cout << std::endl;
	}
}
void ForceField::usecontrols(std::string controlname, const GenChain *genchain) {
	auto &pset = ForceFieldControls::getparameterset(controlname);
	pset.getval("StericWeight", &stericwght_);
	pset.getval("PhiPsiWeight", &phipsiwght_);
	phipsiwght_ *= KBT;
	pset.getval("SideConfWeight", &scconfwght_);
	scconfwght_ *= KBT;
	pset.getval("LSWeight", &lsweight_);
	lsweight_ *= KBT;
	if (lsweight_ > 0.00000001)
		setlstermon();
	pset.getval("PairPackingWeight", &sitepairweight_);
	pset.getval("SCPackingWeight",&scpackingweight_);
	sitepairweight_ *= KBT;
	if (sitepairweight_ > 0.00000001)
		setsitepairtermon();
	pset.getval("LocalHBWeight", &localhbweight_);
	localhbweight_ *= KBT;
	pset.getval("ShadowWeight", &shadowweight_);
	if (shadowweight_ != 0.0) {
		shadowweight_ *= KBT;
		ShadowTerm::setupenefuncs();
		/*		shadowon_.assign(bsinchains_.size(),std::vector<bool>());
		 std::string shadowregions;
		 pset.getval("ShadowRegions",&shadowregions);
		 if(shadowregions.empty()){
		 shadowon_[0].assign(bsinchains_[0].size(),true);
		 } else	{
		 if( bsinchains_.size()>1){
		 std::cout <<"ShadowRegions will only be applied to the first chain."<<std::endl;
		 }
		 std::vector<int> iregions=NSPdataio::stringtoselection(shadowregions);
		 shadowon_[0].assign(bsinchains_[0].size(),false);
		 for(auto s:iregions) shadowon_[0][s]=true;
		 }
		 for(int cidx=1;cidx<bsinchains_.size();++cidx)
		 shadowon_[cidx].assign(bsinchains_[cidx].size(),true);*/
	}
	int scmm;
	pset.getval("SideChainMMSteric",&scmm);
	if(scmm !=0) scmmsteric_=true;
	pset.getval("EAttraction", &eattract_);
	pset.getval("RAttractionOff", &rattractoff_);
	rattractoff_ *= A2NM;
	pset.getval("RAttractionSwitch", &rattractswitch_);
	rattractswitch_ *= A2NM;
	std::vector<int> coreregion;
	pset.getval("CoreRegion", &coreregion);
	std::string coreregionpdb;
	pset.getval("CoreRegionBasePDB", &coreregionpdb);
	if (!coreregion.empty()) {
		coreweights_.assign(natoms_, 0.0);
		for (int i = 0; i < coreregion.size() / 3; ++i) {
		    // fields (index starts from 0): chain_idx, start_idx, end_idx
			int chain = coreregion[3 * i];
			int start = coreregion[3 * i + 1];
			int end = coreregion[3 * i + 2];
			for (int j = start; j <= end; ++j) {
				BSInChain &bs = bsinchains_.at(chain).at(j);
				for (int a : bs.atomids())
					coreweights_[a] = 1.0;
			}
		}
	} else if (!coreregionpdb.empty()) {
		NSPproteinrep::StructPlus sp(coreregionpdb, true);
		assert(sp.nchains() == bsinchains_.size());
		coreweights_.assign(natoms_, 0.0);
		for (int cid = 0; cid < sp.nchains(); ++cid) {
			assert(sp.length(cid) == bsinchains_[cid].size());
			for (int p = 0; p < sp.length(cid); ++p) {
				char sstype = sp.sstype(
						NSPproteinrep::StructPlus::Position(cid, p));
				if (sstype == 'H' || sstype == 'E') {
					for (int a : bsinchains_[cid][p].atomids())
						coreweights_[a] = 1.0;
				}
			}
		}
	} else {
		coreweights_.clear();
	}
	std::vector<std::string> rmsdrst;
	pset.getval("TotalRMSDRestrain", &rmsdrst);
	if (!rmsdrst.empty()) {
		std::string restraintype = rmsdrst[0];
		assert(genchain != nullptr);
		NSPproteinrep::StructPlus refsp(*(genchain->rmsdrefpdb()), true);
		double rmsdref = std::stod(rmsdrst[1]) * A2NM;
		double krmsdres = std::stod(rmsdrst[2]) / (A2NM * A2NM);
		bool backboneonly=false;
		bool sseonly=false;
		if(rmsdrst[3]=="backboneonly") backboneonly=true;
		else if(rmsdrst[3]=="allatoms") backboneonly=false;
		else if(rmsdrst[3]=="sseonly") {
			backboneonly=true;
			sseonly=true;
		}
		std::vector<int> rregion;
		pset.getval("RMSDRestrainedRegion", &rregion);
		std::vector<double> rmsdweights(natoms_,0.0);
		if (!rregion.empty()) {
			for (int i = 0; i < rregion.size() / 3; ++i) {
			    // fields (index starts from 0): chain_idx, start_idx, end_idx
				int chain = rregion[3 * i];
				int start = rregion[3 * i + 1];
				int end = rregion[3 * i + 2];
				for (int j = start; j <= end; ++j) {
					if(sseonly){
							if(refsp.sstype(NSPproteinrep::StructPlus::Position(chain,j))
									=='C') continue;
						}
					BSInChain &bs = bsinchains_.at(chain).at(j);
					for (int a : bs.atomids())
						rmsdweights[a] = 1.0;
					if(!backboneonly){
						SCInChain &sc = scinchains_.at(chain).at(j);
						for(int a=0;a<sc.nscatoms;++a){
							rmsdweights[sc.poffset+a]=1.0;
						}
					}
				}
			}
		} else {
			for(int c=0;c<bsinchains_.size();++c){
				for(int j=0;j<bsinchains_.at(c).size();++j){
					if(sseonly){
						if(refsp.sstype(NSPproteinrep::StructPlus::Position(c,j))
								=='C') continue;
					}
					BSInChain &bs=bsinchains_.at(c).at(j);
					for(int a:bs.atomids()) rmsdweights[a]=1.0;
					if(!backboneonly){
						SCInChain &sc=scinchains_.at(c).at(j);
						for(int a=0;a<sc.nscatoms;++a){
							rmsdweights[sc.poffset+a]=1.0;
						}
					}
				}
			}
		}
		std::vector<double> refcrd(*(genchain->rmsdrefcrd()));

		if (restraintype == "posi_mode") {
		    // posi_mode
            structrestraints_.push_back(
                    StructRestraint(refcrd,rmsdweights,krmsdres,
                            StructRestraint::POSIMODE,rmsdref));
		} else {
		    // total_mode
            structrestraints_.push_back(
                    StructRestraint(refcrd,rmsdweights,krmsdres,
                            StructRestraint::TOTALMODE,rmsdref));
		}
	}
	std::vector<double> rgrest;
	pset.getval("RgRestrain", &rgrest);
	double krg = rgrest[1];
	if (krg > 0.00001 || krg < -0.00001) {
		if (krg > 0)
			krg /= (A2NM * A2NM);
		setrgrestraint(rgrest[0] * A2NM, krg);
	}
	std::vector<double> helixrst;
	pset.getval("HelixRestrains", &helixrst);
	assert(helixrst.size() % 5 == 0);
	if (!helixrst.empty()) {
		for (int i = 0; i < helixrst.size() / 5; ++i) {
		    // fields (index starts from 0): chain_idx, start_idx, end_idx, keep_prob, force_const
			int c = (int) helixrst[5 * i];
			int s = (int) helixrst[5 * i + 1];
			int e = (int) helixrst[5 * i + 2];
			assert(c<bsinchains_.size());
			assert(s<bsinchains_.at(c).size());
			assert(e<bsinchains_.at(c).size());
			addssrestraint(c, s, e, 0, helixrst[5 * i + 3],
					helixrst[5 * i + 4]);
		}
	}

	std::vector<int> extendrst;
	pset.getval("ExtendRestrains", &extendrst);
	if (!extendrst.empty()) {
		for (int i = 0; i < extendrst.size() / 3; ++i) {
		    // fields (index starts from 0): chain_idx, start_idx, end_idx
			int chain = extendrst[3 * i];
			int start = extendrst[3 * i + 1];
			int end = extendrst[3 * i + 2];
			const int RstCaSeqDist = 3; // the sequence distance of the restrained two CA atoms
			assert((end-start) >= RstCaSeqDist);
			for(int s1 = start; s1 <= end-RstCaSeqDist; ++s1) {
				int a1 = bsinchains_.at(chain).at(s1).caid;
				int a2 = bsinchains_.at(chain).at(s1 + RstCaSeqDist).caid;
				double r0 = 10.0 * A2NM;
				double kres = 50000.0;
				disrestraints_.push_back(DisRestraint(a1, a2, -r0, 0, kres));
			}
		}
	}
	std::vector<double> strandrst;
	pset.getval("StrandRestrains", &strandrst);
	if (!strandrst.empty()) {
		for (int i = 0; i < strandrst.size() / 5; ++i) {
		    // fields (index starts from 0): chain_idx, start_idx, end_idx
			int c = (int) strandrst[5 * i];
			int s = (int) strandrst[5 * i + 1];
			int e = (int) strandrst[5 * i + 2];
			assert(c<bsinchains_.size());
			assert(s<bsinchains_[c].size());
			assert(e<bsinchains_[c].size());
			addssrestraint(c, s, e, 1, strandrst[5 * i + 3],
					strandrst[5 * i + 4]);
		}
	}
	std::vector<double> disrst;
	pset.getval("DisRestrains", &disrst);
	if (!disrst.empty()) {
		for (int i = 0; i < disrst.size() / 5; ++i) {
			int a1 = (int) disrst[5 * i];
			int a2 = (int) disrst[5 * i + 1];
			assert(a1<natoms_);
			assert(a2<natoms_);
			double r0 = disrst[5 * i + 2] * A2NM;
			double r1 = disrst[5 * i + 3] * A2NM;
			double kres = disrst[5 * i + 4] / (A2NM * A2NM);
			disrestraints_.push_back(DisRestraint(a1, a2, r0, r1, kres));
		}
	}
	std::vector<double> lwt;
	pset.getval("LoopRegionWeights", &lwt);
	assert(lwt.size() % 4 == 0);
	siteweights_.assign(bsinchains_.size(), std::vector<double>());
	for (int i = 0; i < bsinchains_.size(); ++i) {
		siteweights_[i].assign(bsinchains_[i].size(), 1.0);
	}
	if (!lwt.empty()) {
		for (int i = 0; i < lwt.size() / 4; ++i) {
		    // fields (index starts from 0): chain_idx, start_idx, end_idx
			int c = (int) lwt[4 * i];
			int s = (int) lwt[4 * i + 1];
			int e = (int) lwt[4 * i + 2];
			double w = lwt[4 * i + 3];
			for (int j = s; j <= e; ++j)
				siteweights_.at(c).at(j) = w;
		}
	}
	std::vector<int> strandpairrst;
	pset.getval("StrandPairRestrains", &strandpairrst);
	if (!strandpairrst.empty()) {
		double kres;
		pset.getval("KresStrandPair", &kres);
		kres /= (A2NM * A2NM);
		for (int i = 0; i < strandpairrst.size() / 6; ++i) {
			int idx0 = 6 * i;
			int length = strandpairrst[idx0];
			int para = strandpairrst[idx0 + 1];
			int chain1 = strandpairrst[idx0 + 2];
			int s1 = strandpairrst[idx0 + 3];
			int chain2 = strandpairrst[idx0 + 4];
			int s2 = strandpairrst[idx0 + 5];
			std::vector<std::pair<int, int>> capairs = capairsinstrandpair(
					length, para, s1, s2, bsinchains_[chain1],
					bsinchains_[chain2]);
			for (auto &cp : capairs) {
				disrestraints_.push_back(
						DisRestraint(cp.first, cp.second, 0.55, 0.75, kres));
			}
		}
	}
}
#ifndef _OPENMP
NeighborList::NeighborList(const std::vector<double> &crd, const ForceField &ff,
		const std::vector<bool> *forceoff) {
	std::vector<XYZ> xyz;
	for (int i = 0; i < crd.size(); i += 3) {
		xyz.push_back(XYZ(crd[i], crd[i + 1], crd[i + 2]));
	}
	int natoms = crd.size() / 3;
	neighbors.assign(natoms, std::vector<int>());
	double rcut2 = 1.00; //nm^2
	if(ff.scmmsteric()){
		for (int i = 0; i < natoms; ++i) {
			std::vector<int> & nbrsi = neighbors[i];
			if (forceoff) {
				if (forceoff->at(i))
					continue;
			}
			for (int j = 0; j < natoms; ++j) {
				if (i == j)
					continue;
				if (forceoff) {
					if (!forceoff->at(j) && (j < i))
						continue;
				} else if (j < i)
					continue;
				if (fabs(j - i) < MAXEXCL) {
					if (ff.excluded(i, j))
						continue;
					if (ff.is14(i, j))
						continue;
				}
				double diff2 = (xyz[i] - xyz[j]).squarednorm();
				if (diff2 < rcut2) {
					nbrsi.push_back(j);
					neighbors[j].push_back(i);
				}//rcut2
			} // j
		} //i
	} else {
		for(int m=0;m<ff.bsinchains().size();++m){
			const std::vector<BSInChain> & bsm=ff.bsinchains()[m];
			for(int i=0;i<bsm.size();++i){
				std::vector<int> ia=bsm[i].atomids();
				for(int n=m;n<ff.bsinchains().size();++n){
					const std::vector<BSInChain> &bsn=ff.bsinchains()[n];
					int jstart=0;
					int jend=bsn.size();
					if(n==m) {
						jstart=i+1;
						if(jstart >=jend) continue;
					}
					for(int j=jstart;j<jend;++j){
						std::vector<int> ja=bsn[j].atomids();
						for(auto ii:ia){
							for(auto jj:ja){
								if(forceoff)
									if(forceoff->at(ii) && forceoff->at(jj))continue;
								if(fabs(jj-ii) <MAXEXCL){
									if(ff.excluded(ii,jj)) continue;
									if(ff.is14(ii,jj)) continue;
								}
								double diff2=(xyz[ii]-xyz[jj]).squarednorm();
								if(diff2<rcut2){
									neighbors[ii].push_back(jj);
									neighbors[jj].push_back(ii);
								} //rcut2
							} //ja
						}//ia
					} //j
				} //n
			} //i
		} //m
	} //scmmsteric
}
#endif
double BondTerm::energy(const std::vector<XYZ> & crd,
		std::vector<XYZ> *forces) const {
	std::vector<XYZ> deriv;
	double b = distance(crd[i], crd[j], &deriv);
	double db = b - b0;
	double e = 0.5 * kb * db * db;
	double kbdb = -kb * db;
	forces->at(i) = forces->at(i) + kbdb * deriv[0];
	forces->at(j) = forces->at(j) + kbdb * deriv[1];
	return e;
}
double AngleTerm::energy(const std::vector<XYZ> & crd,
		std::vector<XYZ> *forces) const {
	std::vector<XYZ> deriv;
	double ct = cos_angle(crd[i], crd[j], crd[k], &deriv);
	double dct = ct - cost0;
	double e = 0.5 * kt * dct * dct;
	double ktdct = -kt * dct;
	forces->at(i) = forces->at(i) + ktdct * deriv[0];
	forces->at(j) = forces->at(j) + ktdct * deriv[1];
	forces->at(k) = forces->at(k) + ktdct * deriv[2];
	return e;
}
double ImpDihTerm::energy(const std::vector<XYZ> & crd,
		std::vector<XYZ> *forces) const {
	std::vector<XYZ> deriv;
	double p = torsion(crd[i], crd[j], crd[k], crd[l], &deriv);
	double dp = p - p0;
	if (dp > 3.14159265)
		dp -= 2 * 3.14159265;
	else if (dp < -3.14159265)
		dp += 2 * 3.14159265;
	double e = 0.5 * kp * dp * dp;
	double kf = -kp * dp;
	forces->at(i) = forces->at(i) + kf * deriv[0];
	forces->at(j) = forces->at(j) + kf * deriv[1];
	forces->at(k) = forces->at(k) + kf * deriv[2];
	forces->at(l) = forces->at(l) + kf * deriv[3];
	return e;
}
bool ForceField::sitepairoff(const std::vector<SSCode> &sscodes, int posi1,
		int posi2) const {
	static const int minsep = 5;
	int p1 = posi1 < posi2 ? posi1 : posi2;
	int p2 = posi2 > posi1 ? posi2 : posi1;
	if (p2 - p1 <= minsep)
		return true;
	int ssid1 = sscodes[p1].ssid;
	int ssid2 = sscodes[p2].ssid;
	if (ssid1 == NN_SSTerm::COIL || ssid2 == NN_SSTerm::COIL)
		return false;
	int ssid_p = ssid1;
	if (ssid_p == NN_SSTerm::TERMINUS)
		ssid_p = sscodes[p1 + 1].ssid;
	for (int p = p1 + 1; p < p2; ++p) {
		if (sscodes[p].ssid != ssid_p)
			return false;
	}
	if (ssid2 != ssid_p && ssid2 != NN_SSTerm::TERMINUS)
		return false;
	return true;
}
double ForceField::etrap(double r, double *dfdr) const {
	double ron = rattractoff_;
	double rflat = rattractswitch_;
	double e;
	double dr = ron - rflat;
	if (r >= ron) {
		*dfdr = 0.0;
		e = 1.0;
	} else if (r <= rflat) {
		*dfdr = 0.0;
		e = 0.0;
	} else {
		double x = (r - rflat) / dr;
		if (x <= 0.5) {
			e = 2.0 * x * x;
			*dfdr = 4.0 * x / dr;
		} else {
			e = 1.0 - 2.0 * (1 - x) * (1 - x);
			*dfdr = 4.0 * (1 - x) / dr;
		}
	}
	return e;
}
double ForceField::sideconfene(const std::vector<double> &crd,
		std::vector<std::vector<PhiPsiCodes>> &phipsicodes,
		std::vector<NSPgeometry::XYZ> *xyzf,const ActiveSelections &acts) const {
	const std::vector<std::vector<ConformerCode>> &sccodes=acts.conformercodes();
	int nchains = phipsicodes.size();
	double ene = 0.0;
	int pidx = -1;
	for (int c = 0; c < nchains; ++c) {
		for (int p = 0; p < sccodes[c].size(); ++p) {
			++pidx;
			if(!acts.rotameractive(c,p)) continue;
			if (scinchains_[c][p].kaiatoms.empty())
				continue;
//			assert(p != 0 && p != sccodes[c].size() - 1); //current code does not yet handle terminals properly
			std::vector<DvDxi> dedx;
			double e;
			if(p==0 || p==sccodes[c].size()-1){
				NN_KaiTerm_T kaiterm(sccodes[c][p]);
				e=kaiterm.outvalue(&dedx);
			} else {
				NN_KaiTerm kaiterm(sccodes[c][p]);
				e = kaiterm.outvalue(&dedx);
			}
			ene += e;
			for (auto &d : dedx) {
				xyzf->at(d.first) = xyzf->at(d.first) - scconfwght_ * d.second;
			}
			if (eneanalysismode_) {
				enecomp_.escconf[pidx] = e * scconfwght_;
			}
		}
	}
	return ene * scconfwght_;
}
double ForceField::localhbenergy(const std::vector<double> &crd,
		std::vector<NSPgeometry::XYZ> *xyzf) const {
	int minsep = 5;
	double ene = 0.0;
	int chainid = -1;
	for (auto &bss : bsinchains_) {
		chainid++;
		for (int s1 = 0; s1 < bss.size() - minsep - 4; s1++) {
			for (int s2 = s1 + 4; s2 <= s1 + minsep; ++s2) {
				double w =
						siteweights_[chainid][s1] < siteweights_[chainid][s2] ?
								siteweights_[chainid][s1] :
								siteweights_[chainid][s2];
				if (w == 0)
					continue;
				std::vector<int> atoms;
				atoms.push_back(bss[s1].cid);
				atoms.push_back(bss[s1].oid);
				atoms.push_back(bss[s2].nid);
				atoms.push_back(bss[s2].caid);
				atoms.push_back(bss[s2 - 1].cid);
				std::vector<NSPgeometry::XYZ> rconcac;
				for (auto a : atoms) {
					rconcac.push_back(getxyz(crd, a));
				}
				std::vector<NSPgeometry::XYZ> dedr;
				double ehb = LocalHBTerm::getinstance().forces(rconcac, dedr);
				if (ehb == 0)
					continue;
				ehb *= w;
				w *= localhbweight_;
				if (ehb != 0.0) {
					int idx = 0;
					for (auto a : atoms)
						xyzf->at(a) = xyzf->at(a) - w * dedr[idx++];
				}
				if (eneanalysismode_) {
					enecomp_.elocalhb[s1] += 0.5 * ehb * localhbweight_;
					enecomp_.elocalhb[s2] += 0.5 * ehb * localhbweight_;
				}
				ene += ehb;
			}
		}
		for (int s1 = 1; s1 < bss.size() - minsep - 3; s1++) {
			for (int s2 = s1 + 3; s2 <= s1 + minsep; ++s2) {
				double w =
						siteweights_[chainid][s1] < siteweights_[chainid][s2] ?
								siteweights_[chainid][s1] :
								siteweights_[chainid][s2];
				if (w == 0)
					continue;
				std::vector<int> atoms;
				atoms.push_back(bss[s2].cid);
				atoms.push_back(bss[s2].oid);
				atoms.push_back(bss[s1].nid);
				atoms.push_back(bss[s1].caid);
				atoms.push_back(bss[s1 - 1].cid);
				std::vector<NSPgeometry::XYZ> rconcac;
				for (auto a : atoms) {
					rconcac.push_back(getxyz(crd, a));
				}
				std::vector<NSPgeometry::XYZ> dedr;
				double ehb = LocalHBTerm::getinstance().forces(rconcac, dedr);
				if (ehb == 0)
					continue;
				ehb *= w;
				w *= localhbweight_;
				if (ehb != 0.0) {
					int idx = 0;
					for (auto a : atoms)
						xyzf->at(a) = xyzf->at(a) - w * dedr[idx++];
				}
				if (eneanalysismode_) {
					enecomp_.elocalhb[s2] += 0.5 * ehb * localhbweight_;
					enecomp_.elocalhb[s1] += 0.5 * ehb * localhbweight_;
				}
				ene += ehb;
			}
		}
	}
	return ene * localhbweight_;
}
double ForceField::localbbhbenergy(const std::vector<XYZ> &xyz,
        std::vector<NSPgeometry::XYZ> *xyzf) const {
    int minsep = 5;
    double ene = 0.0;
    int chainid = -1;
    for (auto &bss : bsinchains_) {
        chainid++;
        for (int s1 = 1; s1 < bss.size() -3; s1++) {
            for (int s2 = s1 + 3; s2 <= s1 + minsep; ++s2) {
                if (s2 >= bss.size()) continue;
                double w =
                        siteweights_[chainid][s1] < siteweights_[chainid][s2] ?
                                siteweights_[chainid][s1] :
                                siteweights_[chainid][s2];
                if (w == 0)
                    continue;

                std::vector<int> atoms;
                atoms.push_back(bss[s1-1].cid);
                atoms.push_back(bss[s1-1].oid);
                atoms.push_back(bss[s1].nid);
                atoms.push_back(bss[s2-1].cid);
                atoms.push_back(bss[s2-1].oid);
                atoms.push_back(bss[s2].nid);
                std::vector<DvDxi> dvdx;
                LocalBbHBGeoNNTerm nnmodel;
                nnmodel.setup(xyz, atoms);
                double ehb = nnmodel.outvalue(&dvdx);
                //std::cout << s1 << ' ' << s2 << ' ' << ehb << std::endl;
                if (ehb == 0)
                    continue;
                w *= localhbweight_;
                ehb *= w;
                if (ehb != 0.0) {
                    int idx = 0;
                    for (auto & d : dvdx)
                        xyzf->at(d.first) = xyzf->at(d.first) - w * d.second;
                }
                if (eneanalysismode_) {
                    enecomp_.elocalhb[s1] += 0.5 * ehb;
                    enecomp_.elocalhb[s2] += 0.5 * ehb;
                }
                ene += ehb;
            }
        }
    }
    return ene;
}
double ForceField::shadowenergy(const std::vector<NSPgeometry::XYZ> &xyz,
		std::vector<NSPgeometry::XYZ> *xyzf) const {
	int nchains = bsinchains_.size();
	std::vector<std::vector<NSPpdbstatistics::VirtualSideChain>> vscs;
	std::vector<std::vector<std::vector<NSPgeometry::XYZ>>>bcrds(nchains,std::vector<std::vector<NSPgeometry::XYZ>>());
	vscs.assign(nchains, std::vector<NSPpdbstatistics::VirtualSideChain>());
	for (int n = 0; n < nchains; ++n) {
		const std::vector<BSInChain> &bsn = bsinchains_[n];
		std::vector<NSPpdbstatistics::VirtualSideChain> &vn = vscs[n];
		for (int m = 0; m < bsn.size(); ++m) {
			vn.push_back(
					NSPpdbstatistics::VirtualSideChain(xyz.at(bsn[m].caid),
							xyz.at(bsn[m].nid), xyz.at(bsn[m].cid), 0.1));
			bcrds[n].push_back(std::vector<NSPgeometry::XYZ>());
			bcrds[n].back().push_back(xyz[bsn[m].nid]);
			bcrds[n].back().push_back(xyz[bsn[m].caid]);
			bcrds[n].back().push_back(xyz[bsn[m].cid]);
			bcrds[n].back().push_back(xyz[bsn[m].oid]);
		}
	}
	int minsep = 5;
	double esum = 0.0;
	for (int n = 0; n < nchains; ++n) {
		const std::vector<BSInChain> &bsn = bsinchains_[n];
		for (int m = n; m < nchains; ++m) {
			const std::vector<BSInChain> &bsm = bsinchains_[m];
			for (int i = 0; i < bsn.size(); ++i) {
				int jstart = 0;
				int jend = bsm.size();
				if (m == n) {
					jstart = i + minsep;
				}
				bool ishadow = (scinchains_[n][i].restype == "ZMC");
				if (jstart > jend)
					continue;
				for (int j = jstart; j < jend; ++j) {
					bool jshadow = (scinchains_[m][j].restype == "ZMC");
					if (!(ishadow || jshadow))
						continue;
					double sw =
							siteweights_[n][i] < siteweights_[m][j] ?
									siteweights_[n][i] : siteweights_[m][j];
					if (sw == 0.0)
						continue;
					sw *= shadowweight_;
					double esd;
					if (ishadow && jshadow) {
						std::vector<NSPgeometry::XYZ> f1, f2;
						esd = ShadowTerm::sitepairenergy(bcrds[n][i],
								vscs[n][i], bcrds[m][j], vscs[m][j], &f1, &f2);
						if (esd == 0.0)
							continue;
						(*xyzf)[bsn[i].nid] = (*xyzf)[bsn[i].nid] + sw * f1[0];
						(*xyzf)[bsn[i].caid] = (*xyzf)[bsn[i].caid]
								+ sw * f1[1];
						(*xyzf)[bsn[i].cid] = (*xyzf)[bsn[i].cid] + sw * f1[2];
						(*xyzf)[bsm[j].nid] = (*xyzf)[bsm[j].nid] + sw * f2[0];
						(*xyzf)[bsm[j].caid] = (*xyzf)[bsm[j].caid]
								+ sw * f2[1];
						(*xyzf)[bsm[j].cid] = (*xyzf)[bsm[j].cid] + sw * f2[2];
					} else {
						std::vector<NSPgeometry::XYZ> fs;
						std::vector<DvDxi> dedx;
						if (ishadow) {
							esd = ShadowTerm::sitepairenergy(bcrds[n][i],
									vscs[n][i], xyz, bsinchains_[m][j],
									scinchains_[m][j], &fs, &dedx);
							if (esd == 0.0)
								continue;
							(*xyzf)[bsn[i].nid] = (*xyzf)[bsn[i].nid]
									+ sw * fs[0];
							(*xyzf)[bsn[i].caid] = (*xyzf)[bsn[i].caid]
									+ sw * fs[1];
							(*xyzf)[bsn[i].cid] = (*xyzf)[bsn[i].cid]
									+ sw * fs[2];
						} else {
							esd = ShadowTerm::sitepairenergy(bcrds[m][j],
									vscs[m][j], xyz, bsinchains_[n][i],
									scinchains_[n][i], &fs, &dedx);
							if (esd == 0.0)
								continue;
							(*xyzf)[bsm[j].nid] = (*xyzf)[bsm[j].nid]
									+ sw * fs[0];
							(*xyzf)[bsm[j].caid] = (*xyzf)[bsm[j].caid]
									+ sw * fs[1];
							(*xyzf)[bsm[j].cid] = (*xyzf)[bsm[j].cid]
									+ sw * fs[2];
						}
						for (auto &d : dedx)
							(*xyzf)[d.first] = (*xyzf)[d.first] - sw * d.second;
					}
					esum += esd * sw;
				}
			}
		}
	}
	return esum;
}
double ForceField::scpackingene(const std::vector<NSPgeometry::XYZ> &xyz,
		std::vector<NSPgeometry::XYZ> *xyzf,const ActiveSelections &acts) const {
	int nchains = bsinchains_.size();
	int minsep=1;
	double esum = 0.0;
	double rcut2=1.21;
	for (int n = 0; n < nchains; ++n) {
		const std::vector<BSInChain> &bsn = bsinchains_[n];
		for (int m = n; m < nchains; ++m) {
			const std::vector<BSInChain> &bsm = bsinchains_[m];
			for (int i = 0; i < bsn.size(); ++i) {
				int jstart = 0;
				int jend = bsm.size();
				if (m == n) {
					jstart = i + minsep;
				}
				bool isc = (scinchains_[n][i].nscatoms >0);
				if (jstart > jend)
					continue;
				for (int j = jstart; j < jend; ++j) {
					if(!acts.sidechainpairactive(n,i,m,j)) continue;
					bool jsc = (scinchains_[m][j].nscatoms > 0);
					double rca2=(xyz[bsn[i].caid]-xyz[bsm[j].caid]).squarednorm();
					if(rca2>rcut2) continue;
					if (!(isc || jsc))
						continue;
					double sw =
							siteweights_[n][i] < siteweights_[m][j] ?
									siteweights_[n][i] : siteweights_[m][j];
					if (sw == 0.0)
						continue;
					sw*=scpackingweight_;
					int sep=1000;
					if(m==n){
						sep=j-i;
					}
					double esc=0.0;
					if(jsc){
						std::vector<DvDxi> dedx;
						bool terminal=(j==jend);
						esc +=mcscpackingenergy(xyz, stericatomtypes_,
								bsn[i],scinchains_[m][j],
								sep,&dedx,terminal);
						for(auto &d:dedx)
							(*xyzf)[d.first]=(*xyzf)[d.first]-sw*d.second;
					}
					if(isc){
						std::vector<DvDxi> dedx;
						bool terminal=(i==0);
						esc +=mcscpackingenergy(xyz, stericatomtypes_,
								bsm[j],scinchains_[n][i],
								-sep,&dedx,terminal);
						for(auto &d:dedx)
							(*xyzf)[d.first]=(*xyzf)[d.first]-sw*d.second;
					}
					if (isc && jsc) {
						std::vector<DvDxi> dedx;
						esc +=scscpackingenergy(xyz,stericatomtypes_,
								scinchains_[n][i],scinchains_[m][j],sep,&dedx);
						for (auto &d : dedx)
							(*xyzf)[d.first] = (*xyzf)[d.first] - sw * d.second;
					}
					esum += esc * sw;
				}
			}
		}
	}
	return esum;
}
double ForceField::sitepairenergy(const std::vector<double> & crd,
		const std::vector<std::vector<PhiPsiCodes>> & phipsicodes,
		const std::vector<std::vector<SSCode>> &sscodes,
		std::vector<NSPgeometry::XYZ> *xyzf, double *eattr,
		const ActiveSelections &acts) const {
	double ene = 0.0;
	*eattr = 0.0;
	int nchains = bsinchains_.size();
#ifdef DOCB
	std::vector<std::vector<CBData>> cbdata;
	for(int n=0;n<nchains;++n) {
		const std::vector<BSInChain> &bsn=bsinchains_[n];
		cbdata.push_back(std::vector<CBData>());
		std::vector<CBData> &ccb=cbdata.back();
		for (int m=0;m<bsn.size();++m) {
			std::vector<NSPgeometry::XYZ> crdncac;
			std::vector<int> index = bsn[m].atomids();
			for (auto id : index)
			crdncac.push_back(getxyz(crd, id));
//test
			/*			CBData cbt(crdncac);
			 for(int ii=0;ii<3;++ii){
			 std::vector<NSPgeometry::XYZ> dcbdi(3);
			 for(int jj=0;jj<3;++jj){
			 crdncac[ii][jj] +=0.0001;
			 NSPgeometry::XYZ cb1=CBData(crdncac).cbcrd();
			 crdncac[ii][jj] -=0.0002;
			 NSPgeometry::XYZ cb2=CBData(crdncac).cbcrd();
			 dcbdi[jj]=(cb1-cb2)/0.0002;
			 crdncac[ii][jj] +=0.0001;
			 std::cout <<dcbdi[jj].toString() <<" "<< cbt.drcb()[ii][jj].toString() <<std::endl;
			 }
			 }
			 exit(1);*/
			ccb.push_back(CBData(crdncac));
		}
	}
#endif
	for (int n = 0; n < nchains; ++n) {
		const std::vector<BSInChain> &bsn = bsinchains_[n];
		for (int m = n; m < nchains; ++m) {
			const std::vector<BSInChain> &bsm = bsinchains_[m];
			for (int i = 1; i < bsn.size() - 1; ++i) {
				int jstart = 1;
				int jend = bsm.size() - 1;
				if (m == n) {
					jstart = i + 1;
				}
				for (int j = jstart; j < jend; ++j) {
					if(!acts.sitepairactive(n,i,m,j)) continue;
					if (m == n) {
						if (sitepairoff(sscodes[n], i, j))
							continue;
					}
					SitePairNNTerm spterm;
#ifdef USE_CB
					spterm.setupwithcbdata(crd,bsn,cbdata[n][i],phipsicodes[n],i,bsm,
							cbdata[m][j],phipsicodes[m],j);
#else
					spterm.setup(crd, bsn, phipsicodes[n], i, bsm,
							phipsicodes[m], j);
#endif
					double rca = spterm.rca();
					const std::vector<DvDxi> &drcadx = spterm.drcadx();
					double detrapdr;
					double wij = 1.0;
					if (!coreweights_.empty()) {
						wij = coreweights_[bsn[i].caid]
								* coreweights_[bsm[j].caid];
					}
					*eattr += wij * eattract_ * etrap(rca, &detrapdr);
					detrapdr *= wij * eattract_;
					for (auto &d : drcadx) {
						xyzf->at(d.first) = xyzf->at(d.first)
								- detrapdr * d.second;
					}
#ifdef CBPAIRENE
					std::vector<NSPgeometry::XYZ> drdxcb;
					double rcb=NSPgeometry::distance(cbdata[n][i].cbcrd(),
							cbdata[m][j].cbcrd(),&drdxcb);
					double dedrcb=0.0;
//					double ecbp=CBPairEne::cbforce(rcb+0.0001,&dedrcb);
//					double ecbm=CBPairEne::cbforce(rcb-0.0001,&dedrcb);
//					std::cout <<rcb<<" "<<(ecbp-ecbm)/0.0002 << " "<< dedrcb<<std::endl;
					double ecb=KBT*CBPairEne::cbforce(rcb,&dedrcb);
					/*					double ecb;
					 if (rcb>0.5){
					 ecb=KBT*0.5*100*(rcb-0.5)*(rcb-0.5);
					 dedrcb=100*(rcb-0.5);
					 }*/
//					if((ecbp-ecbm)/0.0002 !=dedrcb){ //for debug
//						ecb=0.0;
//						dedrcb=0.0;
//					}
					dedrcb *=KBT;
					if(dedrcb !=0.0 ) {
						std::vector<NSPgeometry::XYZ> f1=cbdata[n][i].
						distributederiv(dedrcb*drdxcb[0]);
						std::vector<NSPgeometry::XYZ> f2=cbdata[m][j].
						distributederiv(dedrcb*drdxcb[1]);
						std::vector<int> index1=bsn[i].atomids();
						std::vector<int> index2=bsm[j].atomids();
						for(int ii=0;ii<3;++ii) {
							xyzf->at(index1[ii]) =xyzf->at(index1[ii])-f1[ii];
							xyzf->at(index2[ii]) =xyzf->at(index2[ii])-f2[ii];
						}
						if(eneanalysismode_) {
							std::cout<<n <<":"<<i<<"-"<<m<<":"<<j<<" rcb "<<rcb<<" ecb "<<ecb<<std::endl;
						}
					}
					ene +=ecb/sitepairweight_;
#endif
					double sw =
							siteweights_[n][i] < siteweights_[m][j] ?
									siteweights_[n][i] : siteweights_[m][j];
					std::vector<DvDxi> dedx;
					double e = spterm.outvalue(&dedx);

					ene += e * sw;
					if (eneanalysismode_) {
						e *= 0.5 * sitepairweight_;
						int ni = siteid(n, i);
						int mj = siteid(m, j);
						enecomp_.esitepair[ni][mj] = e;
						enecomp_.esitepair[mj][ni] = e;
						enecomp_.epacking[ni] += e;
						enecomp_.epacking[mj] += e;
					}
					sw *= sitepairweight_;
					for (auto &d : dedx) {
						xyzf->at(d.first) = xyzf->at(d.first) - sw * d.second;
					}
				}
			}
		}
	}
	return ene * sitepairweight_;
}
double ForceField::stericenergy(const NeighborList &nbl,
		const std::vector<XYZ> &xyz, std::vector<XYZ> *xyzf) const {
	std::vector<XYZ> pairf(2);
	double ene = 0.0;
	for (int i = 0; i < natoms_; ++i) {
		for (auto j : nbl.neighbors[i]) {
			if (j <= i)
				continue;
			ene += StericAtom::pairenergy(xyz[i], xyz[j], stericatoms_[i],
					stericatoms_[j], &pairf);
			(*xyzf)[i] = (*xyzf)[i] + stericwght_ * pairf[0];
			(*xyzf)[j] = (*xyzf)[j] + stericwght_ * pairf[1];
		}
	}
	return ene * stericwght_;
}

double PhiPsi::energy(const std::vector<PhiPsiCodes> &phipsicodes,
		double phipsiwght, std::vector<XYZ> *forces) const {
	static const double rad = 180.0 / 3.14159265;
	std::vector<XYZ> dphidr;
	std::vector<XYZ> dpsidr;
	auto & codes = phipsicodes[siteid];
	double dedpsi = 0.0;
	double dedphi = 0.0;
	double ene;
        double scale=1.0;
	if (mode == NTERM) {
		ene = distr->intplene_psi(codes.psi * rad, &dedpsi);
                if(ene>0.0) scale=ForceField::highescale;
                ene*=scale;
		dedpsi *= scale*phipsiwght * rad;
		for (auto &d : codes.dpsidx)
			(*forces)[d.first] = (*forces)[d.first] - dedpsi * d.second;
	} else if (mode == CTERM) {
		ene = distr->intplene_phi(codes.phi * rad, &dedphi);
                if(ene>0.0) scale=ForceField::highescale;
                ene*=scale;
		dedphi *= scale*phipsiwght * rad;
		for (auto &d : codes.dpsidx)
			(*forces)[d.first] = (*forces)[d.first] - dedpsi * d.second;
	} else {
		PhiPsiNNTerm<PhiPsiNNTerm<>::MIXCOIL> phipsinnterm;
		std::vector<DvDxi> dedxi;
		ene = phipsinnterm.phipsiene(phipsicodes, siteid, &dedxi);
                if(ene>0.0) scale=ForceField::highescale;
                ene*=scale;
		for (auto &d : dedxi) {
			(*forces)[d.first] = (*forces)[d.first] - scale*phipsiwght * d.second;
		}
	}
	return ene * phipsiwght;
}
double StericAtom::pairenergy(const NSPgeometry::XYZ & pos1,
		const NSPgeometry::XYZ &pos2, const StericAtom & sa1,
		const StericAtom &sa2, std::vector<NSPgeometry::XYZ> *forces) {
	static const double fac = pow(2, 1 / 6.0);
	static const double rcut2 = 0.25; // in nm**2
	forces->assign(2, (0.0, 0.0, 0.0));
	XYZ dr = pos2 - pos1;
	if (dr.squarednorm() > rcut2) {
		return 0.0;
	}
	double sigma;
	if ((sa1.hbdonor && sa2.hbacceptor) || (sa1.hbacceptor && sa2.hbdonor)) {
		sigma = 0.5 * (sa1.sigmahb + sa2.sigmahb);
	} else if ((sa1.hbdonor && sa2.hbdonor)
			|| (sa1.hbacceptor && sa2.hbacceptor)) {
		sigma = 0.55 * (sa1.sigma + sa2.sigma);
	} else {
		sigma = 0.5 * (sa1.sigma + sa2.sigma);
	}
	double eps = sqrt(sa1.eps * sa2.eps);
	double sl = fac * sigma;
	if (dr.squarednorm() > sl * sl) {
		return 0.0;
	}
	double r = distance(pos1, pos2, forces);
	double r2 = r * r;
	double r3 = r2 * r;
	double r5 = r2 * r3;
	double r6 = r3 * r3;
	double r7 = r * r6;
	double r12 = r6 * r6;
	double r13 = r * r12;
	double sigma3 = sigma * sigma * sigma;
	double sigma6 = sigma3 * sigma3;
	double sigma12 = sigma6 * sigma6;
	double ene = 4.0 * eps * (sigma12 / r12 - sigma6 / r6) + eps;
	double dedr = -4.0 * eps * (-12.0 * sigma12 / r13 + 6.0 * sigma6 / r7);
	if (dedr > DEDRMAX)
		dedr = DEDRMAX;
	(*forces)[0] = dedr * (*forces)[0];
	(*forces)[1] = dedr * (*forces)[1];
	return ene;
}
void ForceField::setexcln14() {
	for (auto it = bondterms_.begin(); it != bondterms_.end(); ++it) {
		excluded_[it->i].insert(it->j);
		excluded_[it->j].insert(it->i);
	}
	for (auto it = angleterms_.begin(); it != angleterms_.end(); ++it) {
		excluded_[it->i].insert(it->j);
		excluded_[it->i].insert(it->k);
		excluded_[it->j].insert(it->i);
		excluded_[it->j].insert(it->k);
		excluded_[it->k].insert(it->i);
		excluded_[it->k].insert(it->j);
	}
	for (auto it = impdihterms_.begin(); it != impdihterms_.end(); ++it) {
		excluded_[it->i].insert(it->j);
		excluded_[it->i].insert(it->k);
		excluded_[it->i].insert(it->l);
		excluded_[it->j].insert(it->i);
		excluded_[it->j].insert(it->k);
		excluded_[it->j].insert(it->l);
		excluded_[it->k].insert(it->i);
		excluded_[it->k].insert(it->j);
		excluded_[it->k].insert(it->l);
		excluded_[it->l].insert(it->i);
		excluded_[it->l].insert(it->j);
		excluded_[it->l].insert(it->k);
	}
	for (auto it = dihterms_.begin(); it != dihterms_.end(); ++it) {
		list14_[it->i].insert(it->l);
		list14_[it->l].insert(it->i);
	}
}
double ForceField::lsenergy(const std::vector<PhiPsiCodes> & phipsicodes,
		const std::vector<double> &siteweight,
		std::vector<NSPgeometry::XYZ> *xyzf,
		const std::vector<bool> &lsactive) const {
	double els = 0.0;
	for (int i = 2; i < phipsicodes.size() - 2; ++i) {
		if(!lsactive[i]) continue;
		LSNNTerm lsnnterm;
		std::vector<DvDxi> dedx;
		lsnnterm.setup(phipsicodes, i);
		double w = siteweight[i] * lsweight_;
		double e = lsnnterm.outvalue(&dedx);
		if (e > 0.0)
			w = w * ForceField::highescale;
		els += e * w;
		if (eneanalysismode_) {
			enecomp_.els[enecomp_.siteoffset + i] = e * w;
		}
		for (auto &d : dedx) {
			(*xyzf)[d.first] = (*xyzf)[d.first] - w * d.second;
		}
	}
	return els;
}
#ifndef _OPENMP
/*std::vector<double> ForceField::forces(const std::vector<double> &crd,
		const NeighborList &nbl, std::vector<double> *potenergies,
		const std::vector<bool> *forceoff) const {*/
std::vector<double> ForceField::forces(const std::vector<double> &crd,
		const NeighborList &nbl, std::vector<double> *potenergies,
		const ActiveSelections &acts) const {
	if (eneanalysismode_) {
		int nsites = 0;
		for (auto &bs : bsinchains_) {
			nsites += bs.size();
		}
		enecomp_.init(nsites);
	}
	std::vector<XYZ> xyz;
	assert(natoms_ == crd.size() / 3);
	for (int i = 0; i < crd.size(); i += 3) {
		xyz.push_back(XYZ(crd[i], crd[i + 1], crd[i + 2]));
	}
	std::vector<XYZ> xyzf(natoms_, XYZ(0.0, 0.0, 0.0));
	potenergies->assign(ETERMS, 0.0);
	std::vector<bool> forceoff=acts.atomfixed();
	potenergies->at(EBOND) += covalentenergy(bondterms_, xyz, &xyzf, &forceoff);
	potenergies->at(EANG) += covalentenergy(angleterms_, xyz, &xyzf, &forceoff);
	potenergies->at(EIMPDIH) += covalentenergy(impdihterms_, xyz, &xyzf,
			&forceoff);
	int nchains = bsinchains_.size();
	/*
	std::vector<std::vector<PhiPsiCodes>> phipsicodes;
	for (int i = 0; i < nchains; ++i) {
		phipsicodes.push_back(makephipsicodes(crd, bsinchains_[i]));
	}*/
	/*	std::vector<std::vector<SSCode>> sscodes;
	for (int i = 0; i < nchains; ++i) {
		sscodes.push_back(estimatess(phipsicodes[i]));
	}*/
	acts.calccodes(crd);
	const std::vector<std::vector<PhiPsiCodes>> &phipsicodes=acts.phipsicodes();
	const std::vector<std::vector<SSCode>> & sscodes=acts.sscodes();
	int siteidx = 0;
	for (auto it = phipsis_.begin(); it != phipsis_.end(); ++it) {
		if(!acts.phipsiactive(it->chainid,it->siteid)) {
			siteidx++;
			continue;
		}
		double w = phipsiwght_ * siteweights_[it->chainid][it->siteid];
		double e = it->energy(phipsicodes[it->chainid], w, &xyzf);
		potenergies->at(EPHIPSI) += e;
		if (eneanalysismode_) {
			enecomp_.ephipsi[siteidx++] = e;
		}
	}
	if (lstermon_) {
		enecomp_.siteoffset = 0;
		for (int i = 0; i < nchains; ++i) {
			potenergies->at(ELOCALSTRUCTURE) += lsenergy(phipsicodes[i],
					siteweights_[i], &xyzf,acts.lsactive(i));
			enecomp_.siteoffset += bsinchains_.size();
		}
	}
	if (sitepairtermon_) {
		double etrap;
		potenergies->at(ESITEPAIRS) = sitepairenergy(crd, phipsicodes, sscodes,
				&xyzf, &etrap,acts);
		potenergies->at(ESTRUCTREST) += etrap;
	}
	if (scconfwght_ != 0.0) {
		potenergies->at(ESCCONF) = sideconfene(crd, phipsicodes, &xyzf,acts);
	}
	if (shadowweight_ != 0.0) {
		potenergies->at(ESHADOW) = shadowenergy(xyz, &xyzf);
	}
	for (auto it = structrestraints_.begin(); it != structrestraints_.end();
			++it) {
		potenergies->at(ESTRUCTREST) += it->energy(xyz, &xyzf);
	}
	for (auto it = disrestraints_.begin(); it != disrestraints_.end(); ++it) {
		potenergies->at(ESTRUCTREST) += it->energy(xyz, &xyzf);
	}
	potenergies->at(ERGRESTRAINT) += rgrestraint_.energy(xyz, &xyzf,
			coreweights_);
	for (auto it = ssrestraints_.begin(); it != ssrestraints_.end(); ++it) {
		potenergies->at(ESSRESTRAINT) += it->energy(sscodes, &xyzf);
	}
	if (stericwght_ > 0.0)
		potenergies->at(ESTERIC) = stericenergy(nbl, xyz, &xyzf);
	if (scpackingweight_>0.0)
		potenergies->at(ESCPACKING) = scpackingene(xyz,&xyzf,acts);
	if (localhbweight_ > 0.0)
		potenergies->at(ELOCALHB) = localbbhbenergy(xyz, &xyzf);
	potenergies->at(ETOT) = 0.0;
	for (int i = 1; i < ETERMS; ++i) {
		potenergies->at(ETOT) += potenergies->at(i);
	}
	std::vector<double> forces;
	for (auto it = xyzf.begin(); it != xyzf.end(); ++it) {
		forces.push_back(it->x_);
		forces.push_back(it->y_);
		forces.push_back(it->z_);
	}
	return forces;
}
#endif
BackBoneFF backboneff;
static void make_ff_backbone_nophipsi(const std::vector<int> & nsites,
		ForceField &ff, const std::vector<std::vector<int>> &cissites) {
	int offset = 0;
	int chainid = 0;
	for (int ns : nsites) {
		for (int i = 0; i < ns; ++i) {
			if (i > 0) {
				ff.addbond(offset - 2, offset, backboneff.b0_cn,
						2*KBT*backboneff.kb_cn);
				ff.addangle(offset - 3, offset - 2, offset, backboneff.t0_cacn,
						2*KBT*KANG_FAC*backboneff.kt_cacn);
				ff.addangle(offset - 1, offset - 2, offset, backboneff.t0_ocn,
						2*KBT*KANG_FAC*backboneff.kt_ocn);
				ff.addangle(offset - 2, offset, offset + 1, backboneff.t0_cnca,
						2*KBT*KANG_FAC*backboneff.kt_cnca);
				ff.addimpdih(offset - 2, offset, offset - 3, offset - 1,
						backboneff.p_cncao, 2*KBT*KANG_FAC*backboneff.kp_cncao);
				double p0 = backboneff.p_cacnca;
				double p1 = backboneff.p_ocnca;
				for (auto s : cissites[chainid]) {
					if (i == s) {
						p0 = 0.0;
						p1 = 180.0;
						break;
					}
				}
				ff.addimpdih(offset - 3, offset - 2, offset, offset + 1, p0,
						2*KBT*KANG_FAC*backboneff.kp_cacnca);
				ff.addimpdih(offset - 1, offset - 2, offset, offset + 1, p1,
						2*KBT*KANG_FAC*backboneff.kp_ocnca);
				ff.adddih(offset - 4, offset - 3, offset - 2, offset);
				ff.adddih(offset - 2, offset, offset + 1, offset + 2);
			}
			ff.addbsinchain(chainid, offset, offset + 1, offset + 2,
					offset + 3);
			ff.addbond(offset, offset + 1, backboneff.b0_nca,
					2*KBT*backboneff.kb_nca);
			ff.addbond(offset + 1, offset + 2, backboneff.b0_cac,
					2*KBT*backboneff.kb_cac);
			ff.addbond(offset + 2, offset + 3, backboneff.b0_co,
					2*KBT*backboneff.kb_co);
			ff.addangle(offset, offset + 1, offset + 2, backboneff.t0_ncac,
					2*KBT*KANG_FAC*backboneff.kt_ncac);
			ff.addangle(offset + 1, offset + 2, offset + 3, backboneff.t0_caco,
					2*KBT*KANG_FAC*backboneff.kt_caco);
			ff.adddih(offset, offset + 1, offset + 2, offset + 3);
			ff.setstericatom(offset, backboneff.sigma_n, backboneff.eps, true,
					false, backboneff.sigma_nhb);
			ff.setstericatom(offset + 1, backboneff.sigma_ca, backboneff.eps);
			ff.setstericatom(offset + 2, backboneff.sigma_c, backboneff.eps);
			ff.setstericatom(offset + 3, backboneff.sigma_o, backboneff.eps,
					false, true, backboneff.sigma_ohb);
			offset += backboneff.natomspersites;
		}
		++chainid;
	}
}
ForceField NSPsd::make_forcefield_backbone(const std::vector<int> & nsites,
		std::vector<std::vector<int>> cissites) {
	ForceField ff;
	int natoms = 0;
	for (auto ns : nsites)
		natoms += backboneff.natomspersites * ns;
	ff.init(natoms);
	if (cissites.empty())
		cissites.assign(nsites.size(), std::vector<int>());
	make_ff_backbone_nophipsi(nsites, ff, cissites);
	int chainid = 0;
	for (int ns : nsites) {
		for (int i = 0; i < ns; ++i) {
			if (i == 0) {
				ff.addphipsi(i, chainid, true, false);
			} else if (i == ns - 1) {
				ff.addphipsi(i, chainid, false, true);
			} else {
				ff.addphipsi(i, chainid, false, false);
			}
		}
		++chainid;
	}
	ff.setexcln14();
	return ff;
}
ForceField NSPsd::make_forcefield_backbone(
		const std::vector<std::vector<NSPproteinrep::BackBoneSite>> &chains) {
	ForceField ff;
	std::vector<int> nsites;
	for (auto &c : chains)
		nsites.push_back(c.size());
	int natoms = 0;
	for (auto ns : nsites)
		natoms += backboneff.natomspersites * ns;
	ff.init(natoms);
	std::vector<std::vector<int>> cissites;
	for (auto &c : chains)
		cissites.push_back(NSPproteinrep::findcissites(c));
	make_ff_backbone_nophipsi(nsites, ff, cissites);
	for (int chainid = 0; chainid < chains.size(); ++chainid) {
		auto &chain = chains.at(chainid);
		for (int i = 0; i < nsites[chainid]; ++i) {
			std::string resname = chain[i].resname;
			std::string nextresname("ALA");
			if (i < nsites[chainid] - 1)
				nextresname = chain[i + 1].resname;
			if (i == 0) {
				ff.addphipsi(i, chainid, resname, nextresname, true, false);
			} else if (i == nsites[chainid] - 1) {
				ff.addphipsi(i, chainid, resname, nextresname, false, true);
			} else {
				ff.addphipsi(i, chainid, resname, nextresname, false, false);
			}
		}
	}
	ff.setexcln14();
	return ff;
}

