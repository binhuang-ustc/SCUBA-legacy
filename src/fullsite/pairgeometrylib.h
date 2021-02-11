/*
 * pairgeometry.h
 *
 *  Created on: 2018年1月19日
 *      Author: hyliu
 */

#ifndef FULLSITE_PAIRGEOMETRYLIB_H_
#define FULLSITE_PAIRGEOMETRYLIB_H_
#include "dataio/inputlines.h"
#include "geometry/xyz.h"
#include "dstl/vectortree.h"
#include "fullsite/fullsite.h"
#include <memory>
#include <set>

namespace NSPproteinrep{

struct PairGeometry{
	std::vector<double> torsions; //backbone torsions of two segments in one vector
	std::vector<NSPgeometry::XYZ> crd; //backbone atom coordinates of two residues
	std::vector<double>shiftedtorsions()const;
	std::vector<double>distancematrix() const;
	double rca2() const{
		return (crd[1]-crd[5]).squarednorm();
	}
	PairGeometry(){;}
	PairGeometry(std::vector<FullSite>::const_iterator iters1,
			std::vector<FullSite>::const_iterator iters2);
};
std::vector<PairGeometry> readpairgeometries(std::istream & is);
double pairpairrmsd(const PairGeometry &p1, const PairGeometry &p2);

class PairGeometryLib {
public:
	static std::set<std::string> SpecificPairs;
	static PairGeometryLib & getpglib(const std::string &specificpair);
	struct PairDeviation{
		std::vector<double> devtorsions;
		std::vector<double> devdm;
		double rmsd{100.0};
		PairDeviation(){;}
		PairDeviation(const std::vector<double> &ttor, const std::vector<double> & tdm,
				const std::vector<double> &qtor,const std::vector<double> &qdm);
		bool similar(double tcutrms, double dcutrms, double tcutmax,double dcutmax) const;
		double devtrms()const {
			double devtot2=0.0;
			for(auto t:devtorsions) devtot2 += t*t;
			return sqrt(devtot2/(double) devtorsions.size());
		}
		double devdrms()const {
			double devtot2=0.0;
			for(auto t:devdm) devtot2 += t*t;
			return sqrt(devtot2/(double) devdm.size());
		}
		double devtmax() const{
			double dmax=-1.0;
			for(auto t:devtorsions) if(fabs(t)>dmax) dmax=fabs(t);
			return dmax;
		}
		double devdmax() const{
			double dmax=-1.0;
			for(auto t:devdm) if(fabs(t)>dmax) dmax=fabs(t);
			return dmax;
		}
	};
	void buildtrees(const std::vector<PairGeometry> &pgs);
	std::vector<int> findneighbors(const PairGeometry &query, double rmsdcut,double tcutrms, double dcutrms,
			double tcutmax, double dcutmax,
			std::vector<PairDeviation> &devs);
	PairGeometryLib(const std::string &filename);
	PairGeometryLib(){;}
	double rca2max() const {return rca2max_;}
private:
	domaintree::VectorTree dmtree_;
	domaintree::VectorTree torsiontree_;
	std::shared_ptr<std::vector<std::vector<double>>> dms_;
	std::shared_ptr<std::vector<std::vector<double>>> torsions_;
	std::shared_ptr<std::vector<PairGeometry>> pgs_;
	int npairs_{0};
	double rca2max_{-1.0};
};

}



#endif /* FULLSITE_PAIRGEOMETRYLIB_H_ */
