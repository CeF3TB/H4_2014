#include "../interface/EnergyCalibration.h"

#include <iostream>
#include <fstream>


EnergyCalibration::EnergyCalibration( const std::string& fileName ) {

  setCalibrationFile(fileName);

}


EnergyCalibration::~EnergyCalibration() {};



void EnergyCalibration::setCalibrationFile( const std::string& fileName ) {

  calibConstants_.clear();

  calibFile_ = fileName;

  std::ifstream ifs(fileName.c_str());

  while( ifs.good() ) {

    float calib;
    ifs >> calib;

    calibConstants_.push_back(calib);

  }

  calibConstants_.pop_back(); // remove last additional line

}




void EnergyCalibration::applyCalibration( std::vector<float>& raw ) const {


  if( raw.size() != calibConstants_.size() ) {
    std::cout << "[EnergyCalibration] :: ERROR! The number of channels is different than expected!! Will not calibrate!" << std::endl;
    return;
  }


  for( unsigned i=0; i<calibConstants_.size(); ++i )
    raw[i] *= calibConstants_[i];


}
