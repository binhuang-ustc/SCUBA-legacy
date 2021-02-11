/*
 * scoresideconf.cpp
 *
 *  Created on: 2018年2月8日
 *      Author: hyliu
 */

#include "fullsite/scoresideconf.h"
#include "dstl/randomengine.h"
using namespace NSPproteinrep;
using namespace domaintree;

void ScoreSideConf::buildtree(std::shared_ptr<std::vector<std::vector<double>>> torsions){
	tree_.init(torsions,2.0);
	torsions_=torsions;
	for (long i=0;i<(*torsions).size();++i) {
		tree_.insertpoint(i);
	}
}
void ScoreSideConf::buildreftree(int ntimes){
	reftorsions_=std::shared_ptr<std::vector<std::vector<double>>> (new std::vector<std::vector<double>>());
	auto & realrng=NSPdstl::RandomEngine<>::getinstance().realrng(0.0,1.0);
	int ii=0;
	for (long i=0;i<torsions_->size();++i) {
		std::vector<double> reft=(*torsions_)[i];
		for(int k=0;k<ntimes;++k){
			for(int m=2;m<reft.size();++m){
				reft[m]= 359.9999999*(-0.5+realrng());
			}
			(*reftorsions_).push_back(reft);
			if(ii==0)reftree_.init(reftorsions_,2.0);
			reftree_.insertpoint(ii);
			++ii;
		}
	}
}
const double myd20=64.0;
const double myd2max=900.0;
static double neighborcut2(int n){
	return (double) (n)*myd2max;
}
double ScoreSideConf::pairscore(const std::vector<double> &conf, const std::vector<double> &tmpl){
	double s=1.0;
	for (int i=0;i<tmpl.size();++i) {
		double si;
		double diff2=conf[i]-tmpl[i];
		while(diff2 <-180.0) diff2+=360.0;
		while(diff2 >180.0) diff2 -=360.0;
		diff2=diff2*diff2;
		if(diff2 <= myd20) si=1.0;
		else if(diff2 >=myd2max) si=0.0;
		else {
//			double x=(diff2-myd20)/(myd2max-myd20);
//			si=1-x*x;
			double x=(diff2-myd2max)/(myd20-myd2max);
			si=x*x;
		}
		s*=si;
	}
	return s;
}
double ScoreSideConf::neighborsum(const std::vector<double> &conf,
		TorsionVectorTree &tree){
	double bound2=neighborcut2(conf.size());
	std::vector<std::vector<double>> *tvdata=&(tree.gettorsionvectors());
	domaintree::D2Leaf<long,std::vector<std::vector<double>>,AngleCrd>
		d2leaf(tvdata,1000000,bound2);
	tree.gettree().findneighbor(conf,d2leaf,bound2);
	std::vector<std::pair<long,double>> &neighbors=d2leaf.nnearest().neighbors();
	double s=0.0;
//	double s=NMIN;
	for(auto &n:neighbors){
		if(n.second < 0.00001) continue;  //ignore self
		const std::vector<double> & temp=(*tvdata)[n.first];
		s+=pairscore(conf,temp);
	}
	return s;
}

double ScoreSideConf::score(const std::vector<double> &conf){
	double mins=0.1;
	double s=(neighborsum(conf,tree_)+mins)/(double) tree_.size();
	double sref=(neighborsum(conf,reftree_)+mins)/(double) reftree_.size();
	return s/sref;
}

