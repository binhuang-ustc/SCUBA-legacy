/*
 * scoresideconf.h
 *
 *  Created on: 2018年2月8日
 *      Author: hyliu
 */

#ifndef FULLSITE_SCORESIDECONF_H_
#define FULLSITE_SCORESIDECONF_H_

#include "backbone/torsionvectortree.h"
#include <memory>
namespace NSPproteinrep {
class ScoreSideConf {
public:
	void buildtree(std::shared_ptr<std::vector<std::vector<double>>> torsions);
	void buildreftree(int ntimes);
	double score(const std::vector<double> & conf);
	double neighborsum(const std::vector<double> & conf,TorsionVectorTree &tree);
private:
	TorsionVectorTree tree_;
	TorsionVectorTree reftree_;
	std::shared_ptr<std::vector<std::vector<double>>> torsions_;
	std::shared_ptr<std::vector<std::vector<double>>> reftorsions_;
	static double pairscore(const std::vector<double> &conf, const std::vector<double> & tmpl);
};
}



#endif /* FULLSITE_SCORESIDECONF_H_ */
