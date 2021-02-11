/*
 * contacts.h
 *
 *  Created on: 2018年1月18日
 *      Author: hyliu
 */

#ifndef FULLSITE_CONTACTS_H_
#define FULLSITE_CONTACTS_H_

#include "fullsite/fullsite.h"

namespace NSPproteinrep{
struct AtomContact{

	enum contacttypes{POLARTYPE,HBTYPE,SMTYPE,CGTYPE,ARTYPE,NTYPESCHEMES};
	enum polartypes{NPNP,PNP,PP};  //polar, nonpolar
	enum hbtypes{DA,DD,AA,NOHB};   //hydrogen bonding
	enum smtypes{SS,SM,MM}; //side-chain main chain
	enum cgtypes{POSINEG, POSIPOSI,NEGNEG,CGNEUTRAL,NEUTRALNEUTRAL}; //charged group
	enum artypes{ARAR,CGAR,PAR,NPAR,NOAR};
	typedef std::vector<int> ContactFilter;
	static const ContactFilter hbondfilter;
	static const ContactFilter nonpolarfilter;
	static const ContactFilter sidesidefilter;
	static const ContactFilter sidemainfilter;
	static const ContactFilter mainmainfilter;
	static const ContactFilter saltbridgefilter;
	static const ContactFilter aromaromfilter;
	static const ContactFilter cgaromfilter;
	static std::vector<AtomContact> applyfilter(const std::vector<AtomContact> &contacts,
			const ContactFilter &filter);
	static std::vector<AtomContact> applyfilters(const std::vector<AtomContact> &contacts,
			const std::vector<ContactFilter> &filters){
		std::vector<AtomContact> result=contacts;
		for(auto & f:filters){
			result=applyfilter(result,f);
		}
		return result;
	}
	static std::vector<int> findcontacttype(const std::string & residue1, const std::string &atom1,
			const std::string &residue2, const std::string & atom2);
	AtomContact(){;}
	AtomContact(const std::string &a1,const std::string &a2,const std::vector<int> &ct,
			double d):atomname1(a1),atomname2(a2),contacttype(ct),distance(d){;}
	std::string atomname1;
	std::string atomname2;
	std::vector<int> contacttype;
	double distance{-1.0};
};

class InterSiteContacts{
public:
	InterSiteContacts(){;}
	InterSiteContacts(const FullSite site1, const FullSite site2);
	const std::vector<AtomContact> & atomcontacts()const {return atomcontacts_;}
private:
	const FullSite * site1_{nullptr};
	const FullSite * site2_{nullptr};
	std::vector<AtomContact> atomcontacts_;
};
}



#endif /* FULLSITE_CONTACTS_H_ */
