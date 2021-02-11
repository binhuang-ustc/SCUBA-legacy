/*
 * nrdtopn.h
 *
 *  Created on: 2019/5/15
 *      Author: hyliu
 */

#ifndef DSTL_NRDTOPN_H_
#define DSTL_NRDTOPN_H_

#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <cmath>

/*
 * Store score-sorted, non-redundant items.
 * It makes use of the automated key sorting property of multimap.
 *
 * The functor IsRedundant is required to meet the following requirements:
 *  IsRedundant(double cut) constructs the object, the cutoff parameter defining redundancy may be passed here
 *  bool operator()(const Item1,const Item2, double *diff) determine if two items are redundant
 *  and *diff returns the "distance" between them.
 *
 */
template <typename Item,class IsRedundant>
class NrdTopN:public std::multimap<double,Item>{
public:
    NrdTopN(double neighborcut): isredundant_(neighborcut){;}
    //the new item is compared with existing items, sequentially in ascending score(energy)
    // orders. If it is redundant in comparison with an existing item of lower energy, the new item will not
    // be added. Otherwise the new item replaces the existing item.
    // If there are multiple redundant existing items of higher energies, the one of the smallest difference
    // relative to the new item will be replaces.
    bool additem(double score,Item item) {
        auto toreplace_it=this->end();
        double mindiff=DBL_MAX;
        bool keep=true;
        for(auto iter=this->begin();iter!=this->end();++iter){
            double diff;
            if(isredundant_(item, iter->second,&diff)){
                if(score >= iter->first) {
                    //an existing item is of lower energy than the new item
                    keep=false;
                    break;
                }
                if(diff < mindiff){
                    mindiff=diff;
                    toreplace_it=iter;
                }
            }
        }
        if(!keep) return false;
        if(toreplace_it !=this->end()) this->erase(toreplace_it);
        this->insert(std::make_pair(score,item));
        return true;
    }
private:
    IsRedundant isredundant_;
};

/*
 * This is an example functor, which can be used to define NrdTopN
 * Here the Item type is std::shared_ptr<std::vector<double>>
 */
class IsRedundant {
public:
    IsRedundant(double cut):neighborcut_(cut){;}

//This is an example of how to define the functor
    bool operator()
    (const std::shared_ptr<std::vector<double>> v1,
            const std::shared_ptr<std::vector<double>> v2,double *diff) const {
        double diff2=0.0;
        for(int i=0;i<v1->size();++i){
            double d=v1->at(i)-v2->at(i);
            diff2 +=d*d;
        }
        *diff=sqrt(diff2);
        if(*diff <neighborcut_)return true;
        return false;
    }

private:
    double neighborcut_;
};

#endif /* DSTL_NRDTOPN_H_ */
