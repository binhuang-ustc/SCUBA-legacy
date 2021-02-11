/*
 * trjprocessor.h
 *
 *  Created on: 2018年3月28日
 *      Author: hyliu
 */

#ifndef SD_TRAJIO_H_
#define SD_TRAJIO_H_
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
namespace NSPsd{

class TrajIO {
public:
	TrajIO(){;}
	void createwriter(const std::string & filename){
		ofs_=std::shared_ptr<std::ofstream>(new std::ofstream());
		ofs_->open(filename.c_str());
	}
	void writeframe(const std::vector<double> &data) const {
		int nrecord=data.size();
		(*ofs_)<<nrecord<<std::endl;
		for(int i=0;i<nrecord;++i){
			(*ofs_)<<data[i]<< std::endl;
		}
	}
	void createreader(const std::string &filename){
		ifs_=std::shared_ptr<std::ifstream>(new std::ifstream());
		ifs_->open(filename.c_str());
	}
	bool readframe(std::vector<double> *data) const {
		try{
			int nrecord;
			(*ifs_)>>nrecord;
			if(ifs_->fail()) return false;
//			std::cout<<"number of records in trjframe: "<<nrecord<<std::endl;
			data->resize(nrecord);
			for(int i=0;i<nrecord;++i){
				(*ifs_) >>(*data)[i];
			}
			if(ifs_->fail()) return false;
		} catch (std::exception &e){
			return false;
		}
		return true;
	}
private:
	std::shared_ptr<std::ifstream> ifs_;
	std::shared_ptr<std::ofstream> ofs_;
};

}




#endif /* SD_TRAJIO_H_ */
