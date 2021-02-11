/*
 * localhbterm.h
 *
 *  Created on: 2018年1月26日
 *      Author: hyliu
 */

#ifndef SD_LOCALHBTERM_H_
#define SD_LOCALHBTERM_H_
#include "geometry/calculators.h"
inline double switch1d(double x0, double x1, double x, double *deriv) {
	*deriv = 0.0;
	if (x <= x0)
		return 1.0;
	else if (x >= x1)
		return 0.0;
	if (x <= x0 + 0.5 * (x1 - x0)) {
		double xs = (x - x0) / (x1 - x0);
		double fx = 1.0 - 2.0 * xs * xs;
		*deriv = -4.0 * xs/(x1-x0);
		return fx;
	} else {
		double xs = (x1 - x) / (x1 - x0);
		double fx = 2.0 * xs * xs;
		*deriv = -4.0 * xs/(x1-x0);
		return fx;
	}
}

class LocalHBTerm {
public:
	static LocalHBTerm & getinstance(){
		static LocalHBTerm localhbterm;
		return localhbterm;
	}
	LocalHBTerm() :
			rcut1_(0.34), rcut2_(0.45),e0_(-1.0),kres_(1600.0),ron0_(0.30) {
		double rad = 3.14159265 / 180.0;
		costcut2_ = cos(100.0 * rad);
		costcut1_ = cos(120.0 * rad);
		pcut1_ = 20.0 * rad;
		pcut2_ = 30.0 * rad;
	}
	double forces(const std::vector<NSPgeometry::XYZ>& rconcac,
			std::vector<NSPgeometry::XYZ> &dedr) {
		std::vector<NSPgeometry::XYZ> drdon;
		dedr.assign(5, NSPgeometry::XYZ(0, 0, 0));
		double ron = NSPgeometry::distance(rconcac[1], rconcac[2], &drdon);
		double dfdr, dfdcos;
		double fron = switch1d(rcut1_, rcut2_, ron, &dfdr);
		if (fron == 0.0)
			return 0.0;
		std::vector<NSPgeometry::XYZ> dadcon;
		double costcon = NSPgeometry::cos_angle(rconcac[0], rconcac[1],
				rconcac[2], &dadcon);
		double fcost = switch1d(costcut1_, costcut2_, costcon, &dfdcos);
		if (fcost == 0.0)
			return 0.0;
		std::vector<NSPgeometry::XYZ> dpdncaco;
		double p = NSPgeometry::torsion(rconcac[2], rconcac[3], rconcac[4],
				rconcac[1], &dpdncaco);
		double fp, dfpdp;
		if (p > 0.0) {
			fp = switch1d(pcut1_, pcut2_, p, &dfpdp);
		} else {
			fp = switch1d(pcut1_, pcut2_, -p, &dfpdp);
			dfpdp = -dfpdp;
		}
		if (fp == 0.0)
			return 0.0;
		double eres=e0_;
		double deresdr=0.0;
		if(ron>=ron0_){
			eres +=0.5*kres_*(ron-ron0_)*(ron-ron0_);
			deresdr=kres_*(ron-ron0_);
		}
		double ene = eres * fron * fcost * fp;
		if(ene==e0_) return ene;
		dedr[0] = (eres * fron * fp * dfdcos) * dadcon[0];
		dedr[1] = eres
				* ((fron * fp * dfdcos) * dadcon[1]
						+ (fcost * fp * dfdr) * drdon[0]
						+ (fcost * fron * dfpdp) * dpdncaco[3]) +  (fron*fcost*fp*deresdr)*drdon[0];
		dedr[2] = eres
				* ((fron * fp * dfdcos) * dadcon[2]
						+ (fcost * fp * dfdr) * drdon[1]
						+ (fcost * fron * dfpdp) * dpdncaco[0])+(fron*fcost*fp*deresdr)*drdon[1];
		dedr[3] = (eres * fron * fcost * dfpdp) * dpdncaco[1];
		dedr[4] = (eres * fron * fcost * dfpdp) * dpdncaco[2];
		return ene;
	}
private:
	double rcut1_, rcut2_;
	double costcut1_, costcut2_;
	double pcut1_, pcut2_;
	double e0_;
	double kres_;
	double ron0_;
};

#endif /* SD_LOCALHBTERM_H_ */
