/*
 * contacts.cpp

 *
 *  Created on: 2018年1月18日
 *      Author: hyliu
 */
#include "fullsite/contacts.h"
#include <set>
using namespace NSPproteinrep;
const typename AtomContact::ContactFilter AtomContact::hbondfilter{2,0,-1,-1,-1};
const typename AtomContact::ContactFilter AtomContact::nonpolarfilter{0,3,-1,4,-1};
const typename AtomContact::ContactFilter AtomContact::sidesidefilter{-1,-1,0,-1,-1};
const typename AtomContact::ContactFilter AtomContact::sidemainfilter{-1,-1,1,-1,-1};
const typename AtomContact::ContactFilter AtomContact::mainmainfilter{-1,-1,2,-1,-1};
const typename AtomContact::ContactFilter AtomContact::saltbridgefilter{2,-1,-1,0,-1};
const typename AtomContact::ContactFilter AtomContact::aromaromfilter{0,-1,-1,-1,0};
const typename AtomContact::ContactFilter AtomContact::cgaromfilter{1,-1,-1,3,1};
std::vector<int> AtomContact::findcontacttype(const std::string & residue1,
		const std::string &atom1,
		const std::string &residue2, const std::string & atom2){
		static const std::set<std::string> mcatoms{"N","CA","C","O"};
		const AtomTopo &at1=AtomTopo::getatomtopo(residue1,atom1);
		const AtomTopo &at2=AtomTopo::getatomtopo(residue2,atom2);
		std::vector<int> ct(5,-1);
		bool p1=at1.hbacceptor ||at1.hbdonor ||at1.netcharged !=0;
		bool p2=at2.hbacceptor ||at2.hbdonor ||at2.netcharged !=0;
		if(p1 &&p2) ct[0]=PP;
		else if(p1 || p2) ct[0]=PNP;
		else ct[0]=NPNP;
		if((at1.hbacceptor && at2.hbdonor)||(at1.hbdonor && at2.hbacceptor)){
			ct[1]=DA;
		} else if (at1.hbacceptor &&at2.hbacceptor) {ct[1]=AA;}
		else if(at1.hbdonor &&at2.hbdonor) {ct[1]=DD;}
		else ct[1]=NOHB;
		bool mc1=mcatoms.find(atom1) !=mcatoms.end();
		bool mc2=mcatoms.find(atom2) !=mcatoms.end();
		if(mc1 &&mc2) ct[2]=MM;
		else if(mc1 || mc2) ct[2]=SM;
		else ct[2]=SS;
		if(at1.netcharged*at2.netcharged<0)ct[3]=0;
		else if(at1.netcharged>0 && at2.netcharged>0) ct[3]=1;
		else if(at1.netcharged<0 &&at2.netcharged<0) ct[3]=2;
		else if(at1.netcharged*at2.netcharged==0 && at1.netcharged+at2.netcharged!=0) ct[3]=3;
		else ct[3]=4;
		if(at1.aromatic &&at2.aromatic) ct[4]=0;
		else if((at1.aromatic &&at2.netcharged!=0) ||
				(at1.netcharged!=0 &&at2.aromatic)) ct[4]=1;
		else if((at1.aromatic &&p2) ||(p1 &&at2.aromatic)) ct[4]=2;
		else if(at1.aromatic ||at2.aromatic) ct[4]=3;
		else ct[4]=4;
		return ct;

}
std::vector<AtomContact> AtomContact::applyfilter(const std::vector<AtomContact> &contacts,
			const ContactFilter &filter){
		std::vector<AtomContact> result;
		for(auto &c:contacts){
			bool keep=true;
			for(int i=0;i<NTYPESCHEMES;++i){
				if(filter[i]>=0){
					if(filter[i] != c.contacttype[i]) {
						keep=false;
						break;
					}
				}
			}
			if(keep) result.push_back(c);
		}
		return result;
}
InterSiteContacts::InterSiteContacts(const FullSite site1, const FullSite site2){
	site1_=&site1;
	site2_=&site2;
	const std::vector<std::string> mcatomnames{"N","CA","C","O"};
	std::vector<std::string> atomnames1,atomnames2;
	for(auto a:mcatomnames){
		atomnames1.push_back(a);
		atomnames2.push_back(a);
	}
	for(auto a:site1_->atomnames_sidechaintopo()) atomnames1.push_back(a);
	for(auto a:site2_->atomnames_sidechaintopo()) atomnames2.push_back(a);
	const double rcut2=4.5*4.5;
	for(auto a1:atomnames1){
		NSPgeometry::XYZ crd1=site1_->getcrd(a1);
		for(auto a2:atomnames2){
			NSPgeometry::XYZ crd2=site2_->getcrd(a2);
			double dist2=(crd1-crd2).squarednorm();
			if(dist2>=rcut2) continue;
			std::vector<int> ct=AtomContact::findcontacttype(site1_->resname(),a1,
					site2_->resname(),a2);
			if(ct[AtomContact::HBTYPE]==AtomContact::DA && dist2 >3.5*3.5) continue;
			if(ct[AtomContact::ARTYPE]==AtomContact::NOAR &&dist2 >4.0*4.0) continue;
			atomcontacts_.push_back(AtomContact(a1,a2,ct,sqrt(dist2)));
		}
	}
}
