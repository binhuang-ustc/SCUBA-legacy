/*
 * cbpairene.h
 *
 *  Created on: 2018年1月24日
 *      Author: hyliu
 */

#ifndef SD_CBPAIRENE_H_
#define SD_CBPAIRENE_H_
#include <vector>
#include <string>
namespace NSPsd{

class CBPairEne{
public:
	static double cbforce(double rcb,double *dedrcb);
private:
	void readetable(const std::string &filename);
	double force(double rcb,double *dedrcb) const;
	std::vector<double> etable_;
	double rstep_{0.01};
	int nsteps_{80};
	bool tableread_{false};
};

}



#endif /* SD_CBPAIRENE_H_ */
