#include "../interface/CalibrationUtility.h"

#include <iostream>
#include <cstdlib>



CalibrationUtility::CalibrationUtility( const std::string& tag ) {

  setTag(tag);

}





CalibrationUtility::~CalibrationUtility() {

}



void CalibrationUtility::setTag( const std::string& tag ) {

  tag_ = tag;

  if( tag=="V00" ) {

    tag_cef3_ ="V0";
    tag_bgo_ ="V0";

  }else if( tag=="V01" ) {

    tag_cef3_ ="V1";
    tag_bgo_ ="V0";

  }else if( tag=="dev" ) {

    tag_cef3_ ="dev";
    tag_bgo_ ="V1";

  } else {
    std::cout << "[CalibrationUtility] :: Tag " << tag << " does not exist. Exiting" << std::endl;
    exit(12345);

  }
  
}




std::string CalibrationUtility::getCeF3FileName() const {

  std::string fileName = "CeF3Calibration/constants_" + tag_cef3_ + ".txt";

  return fileName;

}


std::string CalibrationUtility::getBGOFileName() const {

  std::string fileName = "BGOCalibration/constants_" + tag_bgo_ + ".txt";

  return fileName;

}
