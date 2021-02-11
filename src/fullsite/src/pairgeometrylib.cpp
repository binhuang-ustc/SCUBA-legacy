/*
 * pairgeometrylib.cpp
 *
 *  Created on: 2018年1月19日
 *      Author: hyliu
 */
#include "fullsite/pairgeometrylib.h"
#include "geometry/quatfit.h"
#include "dataio/datapaths.h"
#include <set>
#include <fstream>
using namespace NSPproteinrep;
std::set<std::string> PairGeometryLib::SpecificPairs{"ASPARGbridges",
	"ASPLYSbridges","ASPHISbridges","GLUARGbridges","GLULYSbridges","GLUHISbridges"
	};
PairGeometry::PairGeometry(std::vector<FullSite>::const_iterator iters1,
		std::vector<FullSite>::const_iterator iters2){
	std::vector<BackBoneSite> bsa = extractbackbonesegment(
							iters1 - 1, 3);
	std::vector<BackBoneSite> bsb = extractbackbonesegment(
							iters2 - 1, 3);
	torsions.resize(8);
	torsions[0]=bsa[0].psi();
	torsions[1]=bsa[1].phi();
	torsions[2]=bsa[1].psi();
	torsions[3]=bsa[2].phi();
	torsions[4]=bsb[0].psi();
	torsions[5]=bsb[1].phi();
	torsions[6]=bsb[1].psi();
	torsions[7]=bsb[2].phi();
	crd.clear();
	bsa[1].getcrd(crd);
	crd.push_back(bsa[1].cbcrd());
	bsb[1].getcrd(crd);
	crd.push_back(bsb[1].cbcrd());
}
PairGeometryLib::PairGeometryLib(const std::string & filename){
	std::ifstream ifs;
	ifs.open(filename.c_str());
	pgs_=std::shared_ptr<std::vector<PairGeometry>>(new std::vector<PairGeometry>
	(readpairgeometries(ifs)));
	ifs.close();
	if(pgs_->empty()){
		std::cout <<"Failed to read pairgeometrylib from " << filename<<std::endl;
		exit(1);
	}
	buildtrees(*pgs_);
}

PairGeometryLib & PairGeometryLib::getpglib(const std::string &specificpair){
	static std::map<std::string,std::shared_ptr<PairGeometryLib>> pglibs;
	auto res=pglibs.find(specificpair);
	if(res !=pglibs.end()) return *(res->second);
//	std::string file=NSPdataio::datapath()+specificpair+".dat";
	std::string file=NSPdataio::datafilename(specificpair+".dat");
	auto pglib=std::shared_ptr<PairGeometryLib>(new PairGeometryLib(file));
	pglibs.insert(std::make_pair(specificpair,pglib));
	return *pglib;
}
std::vector<double> PairGeometry::shiftedtorsions() const {
	std::vector<double> res;
	for (auto t : torsions) {
		while (t < 0)
			t += 360.0;
		while (t > 360.0)
			t -= 360.0;
		res.push_back(t);
	}
	return res;
}
std::vector<double> PairGeometry::distancematrix() const {
	std::vector<double> res(25);
	int idx = 0;
	for (int i = 0; i < 5; ++i) {
		for (int j = 5; j < 10; ++j) {
			res[idx++] = sqrt((crd[i] - crd[j]).squarednorm());
		}
	}
	return res;
}
double NSPproteinrep::pairpairrmsd(const PairGeometry &p1,const PairGeometry &p2){
	NSPgeometry::QuatFit qfit;
	return sqrt(qfit.setup(p1.crd,p2.crd));
}
std::vector<PairGeometry> NSPproteinrep::readpairgeometries(std::istream & is) {
	char buffer[120];
	int nsection = 0;
	std::vector<PairGeometry> pgs;
	while (is.getline(buffer, 120)) {
		if(buffer[0]=='#') continue;
		std::vector<std::string> pairlines;
		pairlines.push_back(std::string(buffer));
		for (int i = 0; i < 11; i++) {
			is.getline(buffer, 120);
			pairlines.push_back(std::string(buffer));
		}
		pgs.push_back(PairGeometry());
		PairGeometry &pg = pgs.back();
		pg.torsions.resize(8);
		pg.crd.resize(10);
		int idx = 0;
		int cidx = 0;
		for (int i = 0; i < 12; ++i) {
			std::stringstream sstr(pairlines[i]);
			if (i == 0 || i == 6) {
				for (int j = 0; j < 4; ++j)
					sstr >> pg.torsions[idx++];
			} else {
				double x, y, z;
				sstr >> x >> y >> z;
				pg.crd[cidx++] = NSPgeometry::XYZ(x, y, z);
			}
		}
	}
	return pgs;
}

void PairGeometryLib::buildtrees(const std::vector<PairGeometry> &pgs) {
	dms_ = std::shared_ptr < std::vector<std::vector<double>>>(new std::vector<
			std::vector<double>>());
	torsions_ =
			std::shared_ptr < std::vector<std::vector<double>>>(new std::vector<
					std::vector<double>>());
	npairs_ = 0;
	for (auto &pg : pgs) {
		std::vector<double> tor = pg.shiftedtorsions();
		torsions_->push_back(tor);
		std::vector<double> dm = pg.distancematrix();
		dms_->push_back(dm);
		if (npairs_ == 0) {
//			torsiontree_.init(torsions_, 1.0, 0.0, 360.0);
			dmtree_.init(dms_, 0.1, 0.0, 20.0);
		}
//		torsiontree_.insertpoint(npairs_);
		dmtree_.insertpoint(npairs_);
		double rca2=pg.rca2();
		if(rca2max_<rca2)rca2max_=rca2;
		++npairs_;
	}
}
using namespace domaintree;
std::vector<int> PairGeometryLib::findneighbors(const PairGeometry &query, double rmsdcut,
		double tcutrms, double dcutrms, double tcutmax, double dcutmax,
		std::vector<PairDeviation> &devs) {
	std::vector<double> qtor = query.shiftedtorsions();
	std::vector<double> qdm = query.distancematrix();
	std::vector<int> neighborindices;
	double bnd2 = 25 * dcutrms * dcutrms;
	D2Leaf<long, std::vector<std::vector<double>>, UsualCrd> d2leaf(dms_.get(),
			1000000, bnd2);
	dmtree_.gettree().findneighbor(qdm, d2leaf, bnd2);
	std::vector<std::pair<long, double>> &neighbors =
			d2leaf.nnearest().neighbors();
	for (auto &n : neighbors) {
		if (n.second < 0.000001)
			continue;  //ignore self

		const std::vector<double> & ttor = torsions_->at(n.first);
		const std::vector<double> & tdm = dms_->at(n.first);
		PairDeviation dev(ttor, tdm, qtor, qdm);
		if (dev.similar(tcutrms, dcutrms, tcutmax, dcutmax)) {
			double rmsd=pairpairrmsd(query, pgs_->at(n.first));
			if(rmsd <=rmsdcut){
				dev.rmsd=rmsd;
				devs.push_back(dev);
				neighborindices.push_back(n.first);
			}
		}
	}
	return neighborindices;
}

PairGeometryLib::PairDeviation::PairDeviation(const std::vector<double> &ttor,
		const std::vector<double> & tdm, const std::vector<double> &qtor,
		const std::vector<double> &qdm) {
	devtorsions.clear();
	devdm.clear();
	for (int i = 0; i < qtor.size(); ++i) {
		double dtor = qtor[i] - ttor[i];
		if (dtor < -180.0)
			dtor += 360.0;
		if (dtor > 180.0)
			dtor -= 360.0;
		devtorsions.push_back(dtor);
	}
	for (int i = 0; i < qdm.size(); ++i) {
		devdm.push_back(qdm[i] - tdm[i]);
	}
}
bool PairGeometryLib::PairDeviation::similar(double tcutrms, double dcutrms,
		double tcutmax, double dcutmax) const {
	return (devtrms() <= tcutrms) && (devtmax() <= tcutmax)
			&& (devdrms() <= dcutrms) && (devdmax() <= dcutmax);
}
