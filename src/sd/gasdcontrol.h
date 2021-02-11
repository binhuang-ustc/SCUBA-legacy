/*
 * gasdcontrols.h
 *
 *  Created on: 2018年1月5日
 *      Author: hyliu
 */

#ifndef SD_GASDCONTROL_H_
#define SD_GASDCONTROL_H_
#include "geometry/quatfit.h"
#include "dataio/controlfile.h"
#include "dataio/parameters.h"
#include <memory>

namespace NSPsd{
struct GASDPar {
	int maxgenerations;
	int maxstored_config;
	int shrink_generations;
	int nchildren;
	int steps_opt;
	int steps_explore;
	double temperature_opt;
	double temperature_explore;
	int popsize;
	int subsize;
	double rmsd_sep;
	GASDPar(const std::string &controlname);
};
typedef NSPdataio::TypedMultiInstanceControls<GASDPar> GASDControls;
class RmsdScoreForGA {
public:
	RmsdScoreForGA(const GASDPar & gasdpar){
		rmsdcut2_=gasdpar.rmsd_sep*gasdpar.rmsd_sep;
		scoremin_=scale_*rmsdcut2_;
	}
	double operator()(std::shared_ptr<std::vector<double>> crd1,
			std::shared_ptr<std::vector<double>> crd2) const {
		if (crd1 == crd2)
			return scoremin_;
		double rmsd2 = NSPgeometry::QuatFit().setup(*crd1, *crd2);
		if (rmsd2 < rmsdcut2_)
			return scale_ * rmsd2; //larger rmsd,lower score, better fitness
		return scoremin_;
	}

private:
	double rmsdcut2_;
	double scale_{-1.e8};
	double scoremin_;
};
void definegasdcontrol(const std::string & name, const std::vector<std::string> &controllines);
void gasdreadcontrols(const std::string &filename,std::string name="");
void gasdprintcontrols(std::string name,std::ostream &ofs=std::cout);
}

#endif /* SD_GASDCONTROL_H_ */
