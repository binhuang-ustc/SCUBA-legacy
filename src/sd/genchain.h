/*
 * genchain.h
 *
 *  Created on: 2018年1月25日
 *      Author: hyliu
 */

#ifndef SD_GENCHAIN_H_
#define SD_GENCHAIN_H_
#include "sd/forcefield.h"
#include "sd/sdrun.h"
#include "backbone/backbonesite.h"
#include "fullsite/fullsite.h"
#include "fullsite/structplus.h"
#include <memory>
namespace NSPsd{
using namespace NSPproteinrep;
class GenChain{
public:
	GenChain(const std::string &controlname);
	std::vector<std::vector<std::string>> &sctypes(){return sctypes_;}
	const std::vector<std::vector<std::string>> &sctypes() const {return sctypes_;}
	std::string seqstring() const;
	void writepdb(const std::vector<double> & crd, std::ostream &os,double crdtoangstrom) const;
	std::vector<std::vector<FullSite>> crds2fullsites(const std::vector<double> & crds, double crdtoangstrom) const;
	std::vector<double> getcrd(const std::vector<std::vector<BackBoneSite>> &sites,
			bool userefpdb=true) const;
	std::vector<double> getcrd(const std::string &pdbfile) const;
	std::vector<double> getcrd(const std::vector<std::vector<FullSite>> & ) const;
	ForceField make_forcefield(const std::string &controlname) const;
	StochasticDynamics make_sd(const std::string &controlname,
			const std::vector<std::vector<int>> &atomgroups,
			const std::vector<double> &temperatures) const;
	SDRun make_sdrun(const SDRun::SDRunIn &in, unsigned int seed);
	ActiveSelections *setactiveselections(const ForceField *ff);
	ActiveSelections *getactiveselections() const {
		return acts_.get();
	}
	std::shared_ptr<std::vector<std::vector<FullSite>>> initpdb()
	        const { return initpdb_; }
	std::shared_ptr<std::vector<double>> initcrd() const { return initcrd_; }
	std::shared_ptr<std::vector<std::vector<FullSite>>> rmsdrefpdb()
	        const { return rmsdrefpdb_; }
	std::shared_ptr<std::vector<double>> rmsdrefcrd() const { return rmsdrefcrd_; }
	void resetStartConf(const std::vector<std::vector<FullSite>> & newStartConf);
private:
	const std::vector<std::vector<FullSite>> buildInitPDB(const std::vector<std::vector<FullSite>> & srcChains,
	                  const std::vector<std::vector<int>> & scCoordGenMode);
private:
	std::vector<std::vector<std::string>> sctypes_;
	std::vector<std::vector<bool>> softsc_;
	std::string controlname_;
	std::shared_ptr<std::vector<std::vector<NSPproteinrep::FullSite>>> refpdb_{nullptr};
	std::shared_ptr<std::vector<std::vector<NSPproteinrep::FullSite>>> initpdb_{nullptr};
    std::shared_ptr<std::vector<double>> initcrd_{nullptr};
	std::shared_ptr<NSPproteinrep::StructPlus> sprefpdb_{nullptr};
	std::shared_ptr<ActiveSelections> acts_;
	std::shared_ptr<std::vector<double>> inputcrd_{nullptr};
	std::shared_ptr<std::vector<std::vector<NSPproteinrep::FullSite>>> rmsdrefpdb_{nullptr};
	std::shared_ptr<std::vector<double>> rmsdrefcrd_{nullptr};
};

typedef NSPdataio::TypedMultiInstanceControls<GenChain> GenChainControls;
void definegenchaincontrol(std::string name,const std::vector<std::string> &controllines);
void genchainreadcontrols(const std::string &filename,std::string name);
void genchainprintcontrols(std::string name,std::ostream &ofs);
}




#endif /* SD_GENCHAIN_H_ */
