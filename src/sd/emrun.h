/*
 * emrun.h
 *
 *  Created on: 2018年7月5日
 *      Author: hyliu
 */

#ifndef EMRUN_H_
#define EMRUN_H_
#include "sd/forcefield.h"
#include "gsl/gsl_multimin.h"
#include "Eigen/Core"
#include "LBFGS/LBFGS.h"

namespace NSPsd {
class EMRun;
void enededx(const gsl_vector *x,void *emrun,double *ene,gsl_vector *dedx);
double ene(const gsl_vector *x, void *emrun);
void dedx(const gsl_vector *x,void *emrun,gsl_vector *dedx);

class EneFunctor {
public:
	EneFunctor(int n, EMRun *emrun):n_(n),emrun_(emrun){;}
	double operator()(const Eigen::VectorXd &x, Eigen::VectorXd& grad);
private:
	EMRun *emrun_;
	int n_;
};
class GenChain;
class EMRun {
public:
	struct Param{
		double gsl_initstepsize{0.1};
		double gsl_linesearchtol{0.1};
		double lbfsg_epsilon{1e-6};
	};
	Param param;
	enum {BFGS2,BFGS,CONJUGATE_PR,CONJUGATE_FR,SDESCENT,LBFGS};
	EMRun(const GenChain *genchain,const ForceField *ff,
			std::vector<bool> fixedatom=std::vector<bool>());
	bool run(int algorithm,
			int maxsteps,std::vector<double> *outcrd,double *ene,
				std::vector<double> *gradient=nullptr){
		if(algorithm==LBFGS)
			return run_lbfgs(maxsteps,outcrd,ene,gradient);
		else
			return run_gsl(algorithm,maxsteps,outcrd,ene,gradient);
	}
	void writepdb(std::ostream &ofs);
	bool run_gsl(int algorithm,int maxsteps,std::vector<double> *outcrd,double *ene,
			std::vector<double> *gradient);
	bool run_lbfgs(int maxsteps,std::vector<double> *outcrd,double *ene,
			std::vector<double> *gradient=nullptr);
	const ForceField *forcefield() const {return ff_;}
	std::vector<double> &crdbck() {return crdbck_;}
	std::vector<bool> &atomfixed(){return atomfixed_;}
	std::vector<double> &potenergies() {return potenergies_;}
	std::vector<double> &forces() {return forces_;}
	ActiveSelections &acts(){ return *acts_;}
protected:
	const ForceField* ff_;
	const GenChain *genchain_;
	std::vector<bool> atomfixed_;
	std::vector<double> crdbck_;
	std::vector<double> potenergies_;
	std::vector<double> forces_;
	ActiveSelections *acts_;
	int ndof_;
	void set_gsl_vector(const std::vector<double> &source, gsl_vector *dest);
	void copy_gsl_vector(const gsl_vector *source, std::vector<double> &dest);
};

}



#endif /* EMRUN_H_ */
