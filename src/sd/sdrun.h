/*
 * SDRun.h
 *
 *  Created on: 2017年11月9日
 *      Author: hyliu
 */

#ifndef SD_SDRUN_H_
#define SD_SDRUN_H_
#include "sd/stochasticdynamics.h"
#include "sd/forcefield.h"
#include "sd/shakebonds.h"
#include "dstl/randomengine.h"
namespace NSPsd {
class SDRun;
class EmptySDCallBack{
public:
	bool operator() (SDRun &run) const { return false;}
};
class SDRun {
public:
	struct SDRunIn {
		const std::vector<double> *crd { nullptr };
		const std::vector<bool> *fixatom { nullptr };
		bool shakeinitcrd { true };
		SDRunIn(const std::vector<double> & initcrd,std::string controlname="") :
				crd(&initcrd),sdcontrolname(controlname) {
			;
		}
		SDRunIn(const std::vector<double> & initcrd,
				const std::vector<bool> &fixatoms) :
				crd(&initcrd), fixatom(&fixatoms) {
			;
		}
		std::vector<int> cissites;
		std::string sdcontrolname { "" };
		std::vector<std::string> sctypes;
	};
	SDRun() {
		;
	}
	SDRun(std::shared_ptr<StochasticDynamics> sd,
			std::shared_ptr<ForceField> ff,
			std::shared_ptr<ShakeBonds> shakebds) :
			sd_(sd), ff_(ff), shakebds_(shakebds) {
		;
	}
	/*	SDRun(const SDRun &run)  {
	 sd_=run.sd_;
	 ff_=run.ff_;
	 shakebds_=run.shakebds_;
	 }*/
	template<typename SEED>
	void initrandomengine(SEED seed) {
		rng_ = std::shared_ptr<NSPdstl::RandomEngine<>>(
				new NSPdstl::RandomEngine<>);
		rng_->init(seed);
		//std::cout << "First random number:" << rng_->realrng()() << std::endl;
	}
	std::shared_ptr<StochasticDynamics> sd() {
		return sd_;
	}
	std::shared_ptr<ForceField> ff() {
		return ff_;
	}
	EneComp eneanalyze() {
        std::shared_ptr<NeighborList> nbl;
        ff_->setanalysismode(true);
        nbl = std::shared_ptr < NeighborList
                > (new NeighborList(state_->crd, *ff_, &forceoff_));
            std::vector<double> forces = ff_->forces(state_->crd, *nbl,
                    &potenergies_, *acts_);
        ff_->setanalysismode(false);
        return ff_->enecomp();
	}
	std::shared_ptr<ShakeBonds> shakebds() {
		return shakebds_;
	}
	bool initstate(const SDRunIn & in) {
		auto &pset=SDControls::getparameterset(in.sdcontrolname+"_sd");
			vmasses_ = sd_->masses();
		int natoms = in.crd->size() / 3;
//		forceoff_.assign(natoms, false);
		forceoff_=acts_->atomfixed();
		nfixedcrd_=0;
		if(in.fixatom){
			for(int i=0;i<natoms;++i) forceoff_[i]= forceoff_[i] || in.fixatom->at(i);
		}
		for (int i = 0; i < natoms; ++i) {
				if (forceoff_[i]) {
					for (int idx = 3 * i; idx < 3 * i + 3; ++idx) {
						vmasses_[idx] = MASS_MAX;
						++nfixedcrd_;
					}
				}
			if(shakeon_) nfixedcrd_-=(nfixedcrd_/3-1);
		}
		state_ = sd_->make_initstate(*(in.crd), *rng_, &vmasses_);
		buffstate_ = std::shared_ptr < StochasticDynamics::State
				> (new StochasticDynamics::State);
		*buffstate_ = *state_;
		nstepsrun_ = 0;
		// calculate potential energies at starting point
		bool shakesuccess=sd_->shakestate(*state_, shakebds_->shakefunction(),
				in.shakeinitcrd, &vmasses_);
		std::shared_ptr<NeighborList> nbl = std::shared_ptr<NeighborList>(
		                        new NeighborList(state_->crd, *ff_, &forceoff_));
		std::vector<double> forces = ff_->forces(state_->crd, *nbl,
                                                 &potenergies_, *acts_);
		return shakesuccess;
	}

	template<typename CALLBACK=EmptySDCallBack>
	bool runsteps(int nsteps, const CALLBACK &stepcallback=CALLBACK()) {
		std::shared_ptr<NeighborList> nbl;
		for (int i = 0; i < nsteps; ++i) {
			if (i % nblsteps_ == 0)
				nbl = std::shared_ptr < NeighborList
						> (new NeighborList(state_->crd, *ff_, &forceoff_));
			std::vector<double> forces = ff_->forces(state_->crd, *nbl,
					&potenergies_, *acts_);
			bool done = stepcallback(*this);
			if (done)
				return true;
			if (!sd_->leapstep(*state_, *buffstate_, *rng_, forces,
					shakebds_->shakefunction(), &vmasses_))
				return false;
			auto temp = state_;
			state_ = buffstate_;
			buffstate_ = temp;
			++nstepsrun_;
		}
		return true;
	}
	StochasticDynamics::State & state() {
		return *state_;
	}
	const StochasticDynamics::State &state() const {
		return *state_;
	}
	int nstepsrun() const {
		return nstepsrun_;
	}
	double temperature() const {
		double nfreedof;
		if (shakeon_)
			nfreedof = state_->crd.size() -nfixedcrd_- shakebds_->nbonds();
		else
			nfreedof = state_->crd.size()-nfixedcrd_;
		return 2.0 * (sd_->ekin(*state_)) / (nfreedof * KBT);
	}
	const std::vector<double> &potenergies() const {
		return potenergies_;
	}
	bool & shakeon() {
		return shakeon_;
	}
	const bool & shakeon() const {
		return shakeon_;
	}
	int &nblsteps() {return nblsteps_;}
	const int & nblsteps() const {return nblsteps_;}
	void changetemperature(double kbt_new){
		double scale=kbt_new*KBT/sd_->kbT()[0];
		state_->scaletemp(scale);
		sd_->scaletemperatures(scale);
	}
	void changetemperature(double kbt_new,const std::vector<int> &atoms){
		double scale=kbt_new*KBT/sd_->kbT()[3*atoms[0]];
		state_->scaletemp(scale,atoms);
		sd_->scaletemperatures(scale,atoms,3);
	}
	void changetemperature(double kbt_new,int tgrp){
		changetemperature(kbt_new,temperaturegroups_->at(tgrp));
	}
	std::shared_ptr<std::vector<std::vector<int>>> & temperaturegroups(){
		return temperaturegroups_;}
	const std::shared_ptr<std::vector<std::vector<int>>> & temperaturegroups() const {
		return temperaturegroups_;}
	std::vector<double> &bathtemperatures() {return bathtemperatures_;}
	void setactiveselections(ActiveSelections *acts){
		acts_=acts;
	}
private:
	std::shared_ptr<StochasticDynamics> sd_ { nullptr };
	std::shared_ptr<ForceField> ff_ { nullptr };
	std::shared_ptr<ShakeBonds> shakebds_ { nullptr };
	std::shared_ptr<StochasticDynamics::State> state_ { nullptr };
	std::shared_ptr<StochasticDynamics::State> buffstate_ { nullptr };
	std::shared_ptr<NSPdstl::RandomEngine<>> rng_ { nullptr };
	std::shared_ptr<std::vector<std::vector<int>>> temperaturegroups_{nullptr};
	ActiveSelections* acts_;
	std::vector<double> bathtemperatures_;
	std::vector<double> potenergies_;
	std::vector<double> vmasses_;
	std::vector<bool> forceoff_;
	int nfixedcrd_{0};
	int nstepsrun_ { 0 };
	int nblsteps_{50};
	bool shakeon_ { true };
};

SDRun make_backbone_sdrun(const SDRun::SDRunIn &in, unsigned int randomseed);
void sdreadcontrols(const std::string &filename="sdffcontrols.par",std::string name="");
void sdprintcontrols(std::string name,std::ostream &ofs=std::cout);
}

#endif /* SD_SDRUN_H_ */
