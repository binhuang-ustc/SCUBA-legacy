/*
 * genchain.cpp
 *
 *  Created on: 2018年1月25日
 *      Author: hyliu
 */
#include "sd/genchain.h"
#include "proteinrep/pdbrecord.h"
#include "sd/backboneff.h"
#include "proteinrep/intatomkey.h"
#include "sd/sidechainff.h"
#include "fullsite/structplus.h"
#include "dataio/splitstring.h"
#include "sd/sdparam.h"

using namespace NSPsd;
using namespace NSPproteinrep;
void NSPsd::definegenchaincontrol(std::string name,const std::vector<std::string> &controllines){
	std::map<std::string,double> doublepars{};
	std::map<std::string,std::vector<std::string>> stringvecpars;
	std::map<std::string,std::vector<double>> doublevecpars;
	std::map<std::string,std::string> stringpars{
	    {"StartPDBFile",""}, {"SCTypeFrom","refpdb"}, {"RefPDBFile",""},
	    {"SCCoordFrom",""}, {"RMSDRefPDBFile",""}, {"InputCoordFile", ""},
	    {"TypeInCoil","UMC"}, {"TypeInHelix","ZMC"},{"TypeInStrand","ZMC"},
	    {"SeqString",""}, {"AllFixedSites",""}, {"MCFixedSites",""} };
	std::map<std::string,std::vector<int>> intvecpars;
	std::map<std::string,int>intpars;
	GenChainControls::initdefaultkeyvals(name,doublepars,stringpars,intpars,doublevecpars,
			stringvecpars,intvecpars);
	int nsuccess=GenChainControls::adjustvalues(name,controllines);
	if(nsuccess!= controllines.size()) {
		exit(1);
	}
}
GenChain::GenChain(const std::string &controlname):
	controlname_(controlname){
	auto & pset=GenChainControls::getparameterset(controlname);
	std::string startpdbfile, sctypefrom, refpdbfile, seqstring, sccoordfrom,
	            rmsdrefpdbfile, inputcoordfile;
	pset.getval("StartPDBFile", &startpdbfile);
	pset.getval("SCTypeFrom",&sctypefrom);
	pset.getval("RefPDBFile",&refpdbfile);
	pset.getval("SeqString", &seqstring);
	pset.getval("SCCoordFrom", &sccoordfrom);
	pset.getval("RMSDRefPDBFile", &rmsdrefpdbfile);
	pset.getval("InputCoordFile", &inputcoordfile);

	bool scTypeChanged = false; // whether the type of sidechains in startpdb is changed

	if (startpdbfile.empty()) {
	    std::cout <<"No StartPDBFile defined. Execution aborted."<<std::endl;
        exit(1);
	}
    std::vector<std::vector<FullSite>> startpdb =
            NSPproteinrep::readfullsitesfrompdb(startpdbfile);

    initpdb_ = std::shared_ptr<std::vector<std::vector<FullSite>>>(
                new std::vector<std::vector<FullSite>>); // REQUIRED: initialize the pointer

    refpdb_ = std::shared_ptr<std::vector<std::vector<FullSite>>>(
            new std::vector<std::vector<FullSite>>); // REQUIRED: initialize the pointer

	if(refpdbfile.empty()){
	    *refpdb_ = startpdb;
	} else {
        *refpdb_ = NSPproteinrep::readfullsitesfrompdb(refpdbfile);
        sprefpdb_ = std::shared_ptr<StructPlus>(new
                StructPlus(*refpdb_,true));
	}
	if(sctypefrom=="refpdb"){
		for(auto &c:*refpdb_){
			sctypes_.push_back(std::vector<std::string>());
			softsc_.push_back(std::vector<bool>());
			for(auto &aa:c){
				std::string upper;
				for(char a:aa.resname()) upper.push_back(toupper(a));
				if(upper==aa.resname()) softsc_.back().push_back(false);
				else softsc_.back().push_back(true);
				try {
					VSCType::getVSCType(upper);
				} catch (std::exception &e){
					std::cout <<"Undefined residue type in PDB:" <<aa.resname()<<std::endl;
					exit(1);
				}
				sctypes_.back().push_back(upper);
			}
		}
		scTypeChanged = true;
	}
	if(sctypefrom=="sstypeinrefpdb"){
		std::string htype,stype,ctype;
		pset.getval("TypeInHelix",&htype);
		pset.getval("TypeInStrand",&stype);
		pset.getval("TypeInCoil",&ctype);
		StructPlus & sp=*sprefpdb_;
		for(int cid=0;cid<refpdb_->size();++cid){
			sctypes_.push_back(std::vector<std::string>());
			softsc_.push_back(std::vector<bool>());
			for(int s=0;s<(*refpdb_)[cid].size();++s){
				char ss=sp.sstype(StructPlus::Position(cid,s));
				std::string sc;
				if(ss=='H') sc=htype;
				else if(ss=='E') sc=stype;
				else if(ss=='C') sc=ctype;
				std::string upper;
				for(char a:sc) upper.push_back(toupper(a));
				if(upper==sc) softsc_.back().push_back(false);
				else softsc_.back().push_back(true);
				try{
					VSCType::getVSCType(upper);
				} catch (std::exception &e){
					std::cout<<"Undefined sidechain type " <<sc<<std::endl;
					exit(1);
				}
				sctypes_.back().push_back(upper);
			}
		}
		scTypeChanged = true;
	}
	if(sctypefrom=="seqstring"){
		if(seqstring.empty()){
			std::cout<<"No SeqString defined."<<std::endl;
			exit(1);
		}
		std::vector<std::string> subseqs=NSPdataio::wordsInString(seqstring,",;:");
		sctypes_.assign(subseqs.size(),std::vector<std::string>());
		for(int c=0;c<subseqs.size();++c){
			softsc_.push_back(std::vector<bool>());
			for(int s=0;s<subseqs[c].size();++s){
				char aa=subseqs[c][s];
				char upper=toupper(aa);
				if(upper==aa) softsc_.back().push_back(false);
				else softsc_.back().push_back(true);
				sctypes_[c].push_back(VSCType::resnameof(upper));
			}
		}
		scTypeChanged = true;
	}
	std::vector<std::vector<int>> scCoordGenMode;
	for (auto & c : startpdb) {
	    scCoordGenMode.push_back(std::vector<int>());
	    for (auto & s : c) {
	        scCoordGenMode.back().push_back(1); // 1 - default mode initially, use coordinates from startpdb
	    }
	}
    if (sccoordfrom == "startpdb_only") {
        std::string firstchain;
        for (auto & fs : startpdb.at(0)) {
            firstchain.append(fs.resname());
        }
        for (char & a : firstchain) {
            a = toupper(a);
        }
        std::string aastringin;
        for (std::string & name : sctypes_.at(0)) {
            aastringin.append(name);
        }
        for (char & a : aastringin) {
            a = toupper(a);
        }
        if (firstchain != aastringin) {
            std::cout << "AA Sequence is not consistent with the reference PDB, "
                    << "while the sidechain coordinates are derived from startpdb_only"
                    << std::endl;
            exit(1);
        }
        scTypeChanged = false;
    } else if (sccoordfrom == "startpdb_plus_random") {
        std::string firstchain;
        for (auto & fs : startpdb.at(0)) {
            firstchain.append(fs.resname());
        }
        for (char & a : firstchain) {
            a = toupper(a);
        }
        std::string aastringin;
        for (std::string & name : sctypes_.at(0)) {
            aastringin.append(name);
        }
        for (char & a : aastringin) {
            a = toupper(a);
        }
        int numChains = sctypes_.size();
        for (int c = 0; c < numChains; c++) {
            const auto & seq = sctypes_.at(c);
            int numSites = seq.size();
            const auto & chain = startpdb.at(c);
            auto & modeseq = scCoordGenMode.at(c);
            for (int s = 0; s < numSites; s++) {
                const std::string & resRef = seq.at(s);
                const std::string & resIn = chain.at(s).resname();
                if (strncasecmp(resRef.c_str(), resIn.c_str(), resRef.size()) != 0) {
                    modeseq.at(s) = 0;
                }
            }
        }
        scTypeChanged = true;
    } else if (sccoordfrom == "all_random") {
        for (auto & c : scCoordGenMode) {
            for (int & s : c) {
                s = 0; // 0 - use random kai angle for rotatable side chain atoms
            }
        }
        scTypeChanged = true;
    } else {
        std::cout <<"Execution aborted due to invalid GenChain parameter: "
                << sccoordfrom <<std::endl;
        exit(1);
    }

    if (scTypeChanged) {
        std::vector<std::vector<FullSite>> initpdb
                    = buildInitPDB(startpdb, scCoordGenMode);
        *initpdb_ = initpdb;
    } else {
        *initpdb_ = startpdb;
    }

    if(!inputcoordfile.empty()){
        std::ifstream ifs;
        inputcrd_ = std::shared_ptr<std::vector<double>>(new std::vector<double>());
        ifs.open(inputcoordfile.c_str());
        SDInputOutput::readcrd(ifs, *inputcrd_);
        ifs.close();
        initcrd_ = inputcrd_; // user specified input coordinates has higher priority
    } else {
        // otherwise, generate initial coordinates according to the related parameters
        std::vector<double> initcrd = SDInputOutput::extractcrd(*initpdb_);
        for (auto &c : initcrd) c *= A2NM; // convert unit from Angstrom to Nanometer
        initcrd_ = std::shared_ptr<std::vector<double>>(
                        new std::vector<double>(initcrd));
    }

    pset.getval("RMSDRefPDBFile", &rmsdrefpdbfile);
    if(rmsdrefpdbfile.empty()) {
        rmsdrefpdb_ = initpdb_;
        rmsdrefcrd_ = initcrd_;
    } else {
        rmsdrefpdb_ = std::shared_ptr<std::vector<std::vector<FullSite>>>(
                    new std::vector<std::vector<FullSite>>); // REQUIRED: initialize the pointer
        auto sites = NSPproteinrep::readfullsitesfrompdb(rmsdrefpdbfile);
        *rmsdrefpdb_=sites;
        std::vector<double> rmsdrefcrd = SDInputOutput::extractcrd(sites);
        for(auto & c : rmsdrefcrd) c *= A2NM; // convert unit from Angstrom to Nanometer
        rmsdrefcrd_ = std::shared_ptr<std::vector<double>>(
                            new std::vector<double>(rmsdrefcrd));
        assert(rmsdrefcrd_->size() == initcrd_->size());
    }

}
std::string GenChain::seqstring() const {
	std::string seq;
	for(int c=0;c<sctypes_.size();++c){
		if(c>0) seq.push_back(';');
		for(auto & sc:sctypes_[c]){
			auto &vsc=VSCType::getVSCType(sc);
			seq.push_back(vsc.oneletter);
		}
	}
	return seq;
}
std::vector<double> GenChain::getcrd(const std::vector<std::vector<BackBoneSite>> &bssites,
		bool userefpdb) const {
	std::vector<NSPgeometry::XYZ> res;
	double deg = 3.14159265 / 180.0;
	assert(bssites.size()==sctypes_.size());
	for (int c=0;c<bssites.size();++c){
		assert(bssites[c].size()==sctypes_[c].size());
		for(int s=0;s<bssites[c].size();++s){
			const VSCType &vsc=VSCType::getVSCType(sctypes_[c][s]);
			const BackBoneSite &bs=bssites[c][s];
			std::vector<NSPgeometry::XYZ> r(vsc.nscatoms+4);
			r[0]=bs.ncrd();
			r[1]=bs.cacrd();
			r[vsc.nscatoms+2]=bs.ccrd();
			r[vsc.nscatoms+3]=bs.ocrd();
			std::vector<std::vector<std::pair<int,double>>> ics;
			if(userefpdb && refpdb_){
				if(sctypes_[c][s]==(*refpdb_)[c][s].resname()){
					ics=(*refpdb_)[c][s].internalcrds();
				}
			}
			if(ics.empty()){
				ics=vsc.internalcrds;
			}
			assert(ics.size()==vsc.nscatoms);
			int aidx=2;
			for(auto &ic:ics){
				r[aidx++]=NSPgeometry::InternaltoXYZ(r[ic[0].first],r[ic[1].first],
									r[ic[2].first],
									ic[0].second, ic[1].second*deg,ic[2].second*deg);
			}
			for(auto &x:r) res.push_back(x);
		}
	}
	std::vector<double> crd;
	for (auto &r : res) {
		for (int m = 0; m < 3; ++m)
			crd.push_back(r[m]);
	}
	return crd;
}
std::vector<double> GenChain::getcrd(const std::string &pdbfile) const {
	std::vector<std::vector<FullSite>> sites=readfullsitesfrompdb(pdbfile);
	return getcrd(sites);
}
std::vector<double> GenChain::getcrd(const
	std::vector<std::vector<FullSite>> &sites) const{
	std::vector<NSPgeometry::XYZ> res;
	double deg=3.14159265/180.0;
	for (int c=0;c<sites.size();++c){
		assert(sites[c].size()==sctypes_[c].size());
		for(int s=0;s<sites[c].size();++s){
			const VSCType &vsc=VSCType::getVSCType(sctypes_[c][s]);
			std::vector<NSPgeometry::XYZ> r(vsc.nscatoms+4);
			const FullSite &fs=sites[c][s];
			r[0]=fs.getcrd("N");
			r[1]=fs.getcrd("CA");
			r[vsc.nscatoms+2]=fs.getcrd("C");
			r[vsc.nscatoms+3]=fs.getcrd("O");
			bool sccomplete{true};
			for(int a=0;a<vsc.nscatoms;++a){
				if(!fs.hasatomcrd(vsc.atomnames[a])) {
					sccomplete=false;
					break;
				}
			}
			if( sccomplete){
				for(int a=0;a<vsc.nscatoms;++a){
					r[a+2]=fs.getcrd(vsc.atomnames[a]);
				}
			} else {
				auto & ics=vsc.internalcrds;
				assert(ics.size()==vsc.nscatoms);
				int aidx=2;
				for(auto &ic:ics){
					r[aidx++]=NSPgeometry::InternaltoXYZ(r[ic[0].first],r[ic[1].first],
								r[ic[2].first],
								ic[0].second, ic[1].second*deg,ic[2].second*deg);
						}
			}
			for(auto &x:r) res.push_back(x);
		}
	}
	std::vector<double> crd;
	for (auto &r : res) {
		for (int m = 0; m < 3; ++m)
			crd.push_back(r[m]);
	}
	return crd;
}

void GenChain::resetStartConf(const std::vector<std::vector<FullSite>> & newStartConf) {
    *initpdb_ = newStartConf;
    refpdb_ = initpdb_;
    std::vector<double> initcrd = SDInputOutput::extractcrd(*initpdb_);
    for (auto &c : initcrd) c *= A2NM; // convert unit from Angstrom to Nanometer
    initcrd_ = std::shared_ptr<std::vector<double>>(
                            new std::vector<double>(initcrd));
    rmsdrefcrd_ = initcrd_;
    sctypes_.clear();
    for(auto & c : *initpdb_) {
        sctypes_.push_back(std::vector<std::string>());
        for (auto & s : c) {
            sctypes_.back().push_back(s.resname());
        }
    }
}

const std::vector<std::vector<FullSite>> GenChain::buildInitPDB(
        const std::vector<std::vector<FullSite>> & srcChains,
        const std::vector<std::vector<int>> & scCoordGenMode)
{
    std::vector<std::vector<FullSite>> dstChains;

    double deg=3.14159265/180.0;
    NSPdstl::RandomEngine<> &re = NSPdstl::RandomEngine<>::getinstance();
    re.setrealrng(-1.0, 1.0);
    for (int c = 0; c < srcChains.size(); c++) {
        assert(srcChains[c].size() == sctypes_[c].size());
        assert(srcChains[c].size() == scCoordGenMode[c].size());
        dstChains.push_back(std::vector<FullSite>());
        for(int s = 0; s < srcChains[c].size(); s++){
            const FullSite &fs=srcChains[c][s];
            dstChains.back().push_back(FullSite());
            FullSite & dstSite = dstChains.back().back();
            dstSite.pdbid() = fs.pdbid();
            dstSite.chainid() = fs.chainid();
            dstSite.resid() = fs.resid();
            dstSite.resseq() = fs.resseq();
            dstSite.resname() = fs.resname();
            const int mode = scCoordGenMode[c][s];
            const VSCType &vsc=VSCType::getVSCType(sctypes_[c][s]);
            std::vector<NSPgeometry::XYZ> r(vsc.nscatoms+4);
            r[0]=fs.getcrd("N");
            dstSite.addatomcrd("N", r[0]);
            r[1]=fs.getcrd("CA");
            dstSite.addatomcrd("CA", r[1]);
            r[vsc.nscatoms+2]=fs.getcrd("C");
            dstSite.addatomcrd("C", r[vsc.nscatoms+2]);
            r[vsc.nscatoms+3]=fs.getcrd("O");
            dstSite.addatomcrd("O", r[vsc.nscatoms+3]);
            if (mode == 0) {
                // use random value for rotatable kai angles to build side chain
                auto & ics=vsc.internalcrds;
                assert(ics.size()==vsc.nscatoms);
                int aidx=2;
                int scaidx = 0;
                const std::string & resname = vsc.resname;
                //std::cout << "<DEBUG>1dstSite.resname:" << dstSite.resname() << std::endl;
                dstSite.resname() = resname; // override
                //std::cout << "<DEBUG>2dstSite.resname:" << dstSite.resname() << std::endl;
                for (auto &ic : ics) {
                    const std::string & atomname = vsc.atomnames[scaidx++];
/*                    double rand = 1.0;
                    if (VSCType::isrotatablescatom(resname, atomname)) {
                        rand = re.realrng()();
                    }
                    double phi = ic[2].second * deg * rand;*/
                    int tmpaidx = aidx++;
                    r[tmpaidx] = NSPgeometry::InternaltoXYZ(
                                    r[ic[0].first], r[ic[1].first], r[ic[2].first],
                                    ic[0].second, ic[1].second*deg, ic[2].second*deg);
//                                    ic[0].second, ic[1].second*deg, phi);
                    dstSite.addatomcrd(atomname, r[tmpaidx]);
                }
            } else {
                // try to use coordinates from source structure, if side chain incomplete, use ideal VSCType instead as default
                bool sccomplete{true};
                for(int a=0;a<vsc.nscatoms;++a){
                    if(!fs.hasatomcrd(vsc.atomnames[a])) {
                        sccomplete=false;
                        break;
                    }
                }
                if(sccomplete){
                    for(int a=0;a<vsc.nscatoms;++a){
                        r[a+2]=fs.getcrd(vsc.atomnames[a]);
                        dstSite.addatomcrd(vsc.atomnames[a], r[a+2]);
                    }
                } else {
                    auto & ics=vsc.internalcrds;
                    assert(ics.size()==vsc.nscatoms);
                    int aidx=2;
                    int scaidx = 0;
                    for(auto &ic:ics){
                        int tmpaidx = aidx++;
                        r[tmpaidx]=NSPgeometry::InternaltoXYZ(r[ic[0].first],r[ic[1].first],
                                    r[ic[2].first],
                                    ic[0].second, ic[1].second*deg,ic[2].second*deg);
                        dstSite.addatomcrd(vsc.atomnames[scaidx++], r[tmpaidx]);
                    }
                }
            }
        }
    }

    return dstChains;
}

void GenChain::writepdb(const std::vector<double> & crd, std::ostream &os,double crdtoangstrom) const{
	int offset = 0;
	std::vector<PdbRecord> records;
	for(int c=0;c<sctypes_.size();++c){
		for (int i = 0; i < sctypes_[c].size(); ++i) {
			const VSCType &vsc=VSCType::getVSCType(sctypes_[c][i]);
			std::string resname=vsc.pdbname;
			auto nkey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1, "N", c, resname,0);
			NSPgeometry::XYZ ncrd(crd[offset], crd[offset + 1], crd[offset + 2]);
			records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(nkey, ncrd*crdtoangstrom, offset / 3 + 1));
			offset += 3;
			auto cakey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1,"CA", c, resname,0);
			NSPgeometry::XYZ cacrd(crd[offset], crd[offset + 1], crd[offset + 2]);
			records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(cakey,
					cacrd*crdtoangstrom, offset / 3 + 1));
			offset += 3;
			if (vsc.nscatoms>0) {
				for(int m=0;m<vsc.atomnames.size();++m){
					auto akey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1, vsc.atomnames[m], c,
					resname,0);
					NSPgeometry::XYZ r(crd[offset], crd[offset + 1],
							crd[offset + 2]);
					records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(akey,
							r*crdtoangstrom, offset / 3 + 1));
					offset += 3;
				}
			}
			auto ckey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1, "C", c, resname,0);
			NSPgeometry::XYZ ccrd(crd[offset], crd[offset + 1], crd[offset + 2]);
			records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(ckey,
					ccrd*crdtoangstrom, offset / 3 + 1));
			offset += 3;
			auto okey = NSPproteinrep::AtomKeyTypeA::genKey(i + 1, "O", c, resname,0);
			NSPgeometry::XYZ ocrd(crd[offset], crd[offset + 1], crd[offset + 2]);
			records.push_back(make_pdbrecord<AtomKeyTypeA,NSPgeometry::XYZ>(okey,
					ocrd*crdtoangstrom, offset / 3 + 1));
			offset += 3;
		}
	}
	for (auto &r : records) {
		os << r.toString()<<std::endl;
	}
}

std::vector<std::vector<FullSite>> GenChain::crds2fullsites(const std::vector<double> & crds,
                                                 double crdtoangstrom) const {
    int offset = 0;
    std::vector<std::vector<FullSite>> newSites(*initpdb_);
    for(int c = 0; c < newSites.size(); c++){
        for (int i = 0; i < newSites[c].size(); i++) {
            FullSite & site = newSites[c][i];
            const VSCType & vsc=VSCType::getVSCType(sctypes_[c][i]);
            NSPgeometry::XYZ ncrd(crds[offset], crds[offset + 1], crds[offset + 2]);
            offset += 3;
            site.updateAtomCrd("N", ncrd * crdtoangstrom);
            NSPgeometry::XYZ cacrd(crds[offset], crds[offset + 1], crds[offset + 2]);
            offset += 3;
            site.updateAtomCrd("CA", cacrd * crdtoangstrom);
            if (vsc.nscatoms>0) {
                for(int m=0;m<vsc.atomnames.size();++m){
                    NSPgeometry::XYZ r(crds[offset], crds[offset + 1],
                            crds[offset + 2]);
                    site.updateAtomCrd(vsc.atomnames[m], r * crdtoangstrom);
                    offset += 3;
                }
            }
            NSPgeometry::XYZ ccrd(crds[offset], crds[offset + 1], crds[offset + 2]);
            offset += 3;
            site.updateAtomCrd("C", ccrd * crdtoangstrom);
            NSPgeometry::XYZ ocrd(crds[offset], crds[offset + 1], crds[offset + 2]);
            offset += 3;
            site.updateAtomCrd("O", ocrd * crdtoangstrom);
        }
    }
    return newSites;
}

static void make_ff_chain_nophipsi(const std::vector<int> & nsites,
		ForceField &ff, const std::vector<std::vector<std::string>> &sctypes,
		const std::vector<std::vector<int>> &cissites) {
	int offset = 0;
	int chainid = 0;
	int nscatoms = 0;
	std::vector<int> & stericatomtypes=ff.stericatomtypes();
	stericatomtypes.clear();
	for (int ns : nsites) {
		for (int i = 0; i < ns; ++i) {
			int nscatoms_p=nscatoms;
			const VSCType & vsc=VSCType::getVSCType(sctypes[chainid][i]);
			nscatoms = vsc.nscatoms;
			if (i > 0) {
				ff.addbond(offset - 2, offset, backboneff.b0_cn,
						2*KBT*backboneff.kb_cn);
				ff.addangle(offset - 3 - nscatoms_p, offset - 2, offset,
						backboneff.t0_cacn,2*KBT*KANG_FAC*backboneff.kt_cacn);
				ff.addangle(offset - 1, offset - 2, offset, backboneff.t0_ocn,
						2*KBT*KANG_FAC*backboneff.kt_ocn);
				ff.addangle(offset - 2, offset, offset + 1, backboneff.t0_cnca,
						2*KBT*KANG_FAC*backboneff.kt_cnca);
				ff.addimpdih(offset - 2, offset, offset - 3-nscatoms_p, offset - 1,
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
				ff.addimpdih(offset - 3 - nscatoms_p, offset - 2, offset,
						offset + 1, p0, 2*KBT*KANG_FAC*backboneff.kp_cacnca);
				ff.addimpdih(offset - 1, offset - 2, offset, offset + 1, p1,
						2*KBT*KANG_FAC*backboneff.kp_ocnca);
				ff.adddih(offset - 4 - nscatoms_p, offset - 3 - nscatoms_p,
						offset - 2, offset);
				ff.adddih(offset - 2, offset, offset + 1,
						offset + 2 + nscatoms);
			}
			ff.addbsinchain(chainid, offset, offset + 1, offset + 2 + nscatoms,
					offset + 3 + nscatoms);
			std::vector<std::vector<int>> kaiatoms;
			for(int rsc:vsc.rotameratoms){
				kaiatoms.push_back(std::vector<int>());
				kaiatoms.back().push_back(rsc+offset+2);
				for(int i=0;i<3;++i){
					kaiatoms.back().push_back(vsc.internalcrds[rsc][i].first+offset);
				}
			}
/*			if(!kaiatoms.empty() &&(i==0 || i==ns-1)){
				std::cout <<"Current code does not support residues with side chain torsions at chain terminus"<<std::endl;
				exit(1);
			}*/
			if(ff.scinchains().size()<chainid+1) ff.scinchains().resize(chainid+1);
			ff.scinchains()[chainid].push_back(SCInChain(vsc.resname,kaiatoms,offset+2,nscatoms));
			ff.addbond(offset, offset + 1, backboneff.b0_nca,
					2*KBT*backboneff.kb_nca);
			ff.addbond(offset + 1, offset + 2 + nscatoms, backboneff.b0_cac,
					2*KBT*backboneff.kb_cac);
			ff.addbond(offset + 2+nscatoms, offset + 3 + nscatoms, backboneff.b0_co,
					2*KBT*backboneff.kb_co);
			ff.addangle(offset, offset + 1, offset + 2 + nscatoms,
					backboneff.t0_ncac, 2*KBT*KANG_FAC*backboneff.kt_ncac);
			ff.addangle(offset + 1, offset + 2 + nscatoms,
					offset + 3 + nscatoms, backboneff.t0_caco,
					2*KBT*KANG_FAC*backboneff.kt_caco);
			ff.adddih(offset, offset + 1, offset + 2 + nscatoms,
					offset + 3 + nscatoms);
			ff.setstericatom(offset, backboneff.sigma_n, backboneff.eps, true,
					false, backboneff.sigma_nhb);
			ff.setstericatom(offset + 1, backboneff.sigma_ca, backboneff.eps);
			stericatomtypes.push_back(VSCType::getstericatomtype(vsc.resname,"N"));
			stericatomtypes.push_back(VSCType::getstericatomtype(vsc.resname,"CA"));
			if (nscatoms>0) {
				int aidx=2;
				for(auto & a:vsc.atomnames){
					int atype=VSCType::getstericatomtype(vsc.resname,a);
					double sigma=VSCType::packingatomtypes[atype].radius;
					int hbtype =VSCType::packingatomtypes[atype].hbtype;
					bool hbdonor=(hbtype==1 ||hbtype==3);
					bool hbacceptor=(hbtype==2||hbtype==3);
					double sigmahb=2.8;
					ff.setstericatom(offset + aidx++, sigma, 0.2,hbdonor,hbacceptor,sigmahb);
					stericatomtypes.push_back(atype);
				}
				for(int b=0;b<vsc.newbonds.size();++b){
					ff.addbond(offset+vsc.newbonds[b].first,offset+vsc.newbonds[b].second,
							vsc.b0[b],vsc.kb0[b]);
				}
				for(int a=0;a<vsc.newangles.size();++a){
						ff.addangle(offset+vsc.newangles[a][0],
								offset+vsc.newangles[a][1], offset+vsc.newangles[a][2],
								vsc.a0[a],vsc.ka0[a]);
					}

				for(int a=0;a<vsc.newimpdihs.size();++a){
						ff.addimpdih(offset+vsc.newimpdihs[a][0],
								offset+vsc.newimpdihs[a][1], offset+vsc.newimpdihs[a][2],
								offset+vsc.newimpdihs[a][3],
								vsc.imp0[a],vsc.kimp0[a]);
					}
				for(int a=0;a<vsc.newtorsions.size();++a){
					bool keep=true;
					for(auto k:vsc.newtorsions[a]){
						if(k<0 && i==0) keep=false;
						if(k>=nscatoms+4 && i==ns-1) keep=false;
					}
					if(!keep) continue;
					ff.adddih(offset+vsc.newtorsions[a][0],
								offset+vsc.newtorsions[a][1], offset+vsc.newtorsions[a][2],
								offset+vsc.newtorsions[a][3]);
				}
			}
			ff.setstericatom(offset + 2 + nscatoms, backboneff.sigma_c,
					backboneff.eps);
			ff.setstericatom(offset + 3 + nscatoms, backboneff.sigma_o,
					backboneff.eps, false, true, backboneff.sigma_ohb);
			stericatomtypes.push_back(VSCType::getstericatomtype(vsc.resname,"C"));
			stericatomtypes.push_back(VSCType::getstericatomtype(vsc.resname,"O"));
			offset += backboneff.natomspersites + nscatoms;
		}
		++chainid;
	}
}
ForceField GenChain::make_forcefield(const std::string &controlname) const{
	ForceField ff;
	std::vector<std::vector<int>> cissites;
	if(refpdb_){
		for(auto &chain:*refpdb_){
			std::vector<BackBoneSite> bs=backbone(chain);
			cissites.push_back(findcissites(bs));
		}
	}
	if (cissites.empty())
		cissites.assign(sctypes_.size(), std::vector<int>());
	int natoms = 0;
	std::vector<int>nsites;
	for (int i = 0; i < sctypes_.size(); ++i) {
		nsites.push_back(sctypes_[i].size());
		for (int j = 0; j < sctypes_[i].size(); ++j) {
			natoms += 4+VSCType::getVSCType(sctypes_[i][j]).nscatoms;
		}
	}
	ff.init(natoms);
	make_ff_chain_nophipsi(nsites, ff, sctypes_, cissites);
	auto & scs=ff.scinchains();
	for(int c=0;c<scs.size();++c){
		for(int p=0;p< scs[c].size();++p){
			scs[c][p].softpacking=softsc_[c][p];
		}
	}
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
	ff.usecontrols(controlname,this);
	return ff;
}
StochasticDynamics GenChain::make_sd(const std::string &controlname,
		const std::vector<std::vector<int>> &atomgroups, const std::vector<double> & temperatures) const{
	std::vector<double> masses;
	std::vector<double> gammas;
	double gamma;
	NSPdataio::ParameterSet &pset=SDControls::getparameterset(controlname);
	pset.getval("FrictionCoeff",&gamma);
	for(int c=0;c<sctypes_.size();++c){
		for (int i = 0; i < sctypes_[c].size(); ++i) {
			int natoms=4+VSCType::getVSCType(sctypes_[c][i]).nscatoms;
			for(int a=0;a<natoms;++a)
				masses.push_back(14.0);
			for (int d = 0; d < natoms; ++d)
				gammas.push_back(gamma);
		}
	}
	StochasticDynamics sd;
	double timestep;
	pset.getval("TimeStep",&timestep);
	sd.init(masses, gammas, timestep, atomgroups, temperatures,3);
	return sd;
}
ActiveSelections * GenChain::setactiveselections(const ForceField *ff){
	std::string allfixed;
	std::string mcfixed;
	auto & pset=GenChainControls::getparameterset(controlname_);
	pset.getval("AllFixedSites",&allfixed);
	pset.getval("MCFixedSites",&mcfixed);
	acts_=std::shared_ptr<ActiveSelections>(new ActiveSelections(ff,allfixed,mcfixed));
	return acts_.get();
}
SDRun GenChain::make_sdrun(const SDRun::SDRunIn &in, unsigned int seed){
	auto ff = std::shared_ptr < ForceField > (new ForceField);
	auto &pset = SDControls::getparameterset(in.sdcontrolname + "_sd");
	*ff=make_forcefield(in.sdcontrolname + "_ff");
	auto shakebds = std::shared_ptr < ShakeBonds > (new ShakeBonds);
	int doshake;
	pset.getval("DoShake", &doshake);
	*shakebds = make_shakebonds(*ff);
	shakebds->seton(doshake != 0);
	std::vector<std::string> tgroups;
	pset.getval("TemperatureGroups",&tgroups);
	auto agrps=std::shared_ptr<std::vector<std::vector<int>>>(
			new std::vector<std::vector<int>>(tgroups.size()));
	int na=0;
	for(int i=0;i<tgroups.size();++i){
		if(tgroups[i]=="all") {
			for(int a=0;a<ff->natoms();++a) (*agrps)[i].push_back(a);
		} else	if(tgroups[i]=="mainchain"){
			(*agrps)[i]=ff->mainchainatoms();
		} else if(tgroups[i] =="sidechain"){
			(*agrps)[i]=ff->sidechainatoms();
		}  else if(tgroups[i] =="ssregions"){
			assert(sprefpdb_ !=nullptr);
			for(int c=0;c<sctypes_.size();++c){
				for(int p=0;p<sctypes_[c].size();++p){
					if(sprefpdb_->sstype(StructPlus::Position(c,p)) =='C') continue;
					std::vector<int> ia=ff->atomsinresidue(c,p);
					for(int a:ia) (*agrps)[i].push_back(a);
				}
			}
		} else if(tgroups[i]=="coilregions"){
			for(int c=0;c<sctypes_.size();++c){
				for(int p=0;p<sctypes_[c].size();++p){
					if(sprefpdb_->sstype(StructPlus::Position(c,p)) !='C') continue;
					std::vector<int> ia=ff->atomsinresidue(c,p);
					for(int a:ia) (*agrps)[i].push_back(a);
				}
			}
		} else {
			std::cout<<"Undefined TemperatureGroup Type: " <<tgroups[i]<<std::endl;
			exit(1);
		}
		na +=(*agrps)[i].size();
	}
	std::vector<double> temperatures;
	pset.getval("Temperatures",&temperatures);
	for(auto &t:temperatures) t*=KBT;
	assert (na==ff->natoms());
	auto sd = std::shared_ptr < StochasticDynamics > (new StochasticDynamics);
	*sd = make_sd(in.sdcontrolname + "_sd", *agrps,temperatures);
	SDRun sdrun(sd, ff, shakebds);
	sdrun.temperaturegroups()=agrps;
	sdrun.bathtemperatures()=temperatures;
	sdrun.shakeon() = doshake != 0;
	int nblsteps;
	pset.getval("NeighborListSteps", &nblsteps);
	sdrun.nblsteps() = nblsteps;
	sdrun.initrandomengine(seed);
	std::vector<bool> fixatoms(ff->natoms(),false);
	int  fixmainchain;
	pset.getval("FixMainChain",&fixmainchain);
	if(fixmainchain!=0){
		std::vector<int> mcatoms=ff->mainchainatoms();
		for(int a:mcatoms) fixatoms[a]=true;
	}
	SDRun::SDRunIn newin(*(in.crd),fixatoms);
	this->setactiveselections(ff.get());
	sdrun.setactiveselections(acts_.get());
	sdrun.initstate(newin);
	return sdrun;
}
void NSPsd::genchainreadcontrols(const std::string &filename,std::string name){
	NSPdataio::ControlFile cf;
	cf.readfile(filename);
	std::vector<std::string> sdcontrolines=cf.getcontrolines("SD");
	std::vector<std::string> ffcontrolines=cf.getcontrolines("ForceField");
	std::vector<std::string> genchaincontrollines=cf.getcontrolines("GenChain");
	definesdcontrol(name+"_sd",sdcontrolines);
	defineforcefieldcontrol(name+"_ff",ffcontrolines);
	definegenchaincontrol(name+"_genchain",genchaincontrollines);
}
void NSPsd::genchainprintcontrols(std::string name,std::ostream &ofs){
	ofs<<"START GenChain"<<std::endl;
	GenChainControls::getparameterset(name+"_genchain").printparameters(ofs);
	ofs<<"END GenChain"<<std::endl;
	ofs<<"START SD"<<std::endl;
	SDControls::getparameterset(name+"_sd").printparameters(ofs);
	ofs<<"END SD"<<std::endl;
	ofs<<"START ForceField"<<std::endl;
	ForceFieldControls::getparameterset(name+"_ff").printparameters(ofs);
	ofs<<"END ForceField"<<std::endl;
}


