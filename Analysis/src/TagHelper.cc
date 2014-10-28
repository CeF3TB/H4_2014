#include "../interface/TagHelper.h"

#include <iostream>
#include <cstdlib>



TagHelper::TagHelper( const std::string& tag ) {

  setTag(tag);

}





TagHelper::~TagHelper() {

}



void TagHelper::setTag( const std::string& tag ) {

  tag_ = tag;

  if( tag=="V00" ) {

    tag_cef3_ ="V0";
    tag_bgo_ ="V0";
    tag_align_ ="V0";

  }else if( tag=="V01" ) {

    tag_cef3_ ="V1";
    tag_bgo_ ="V0";
    tag_align_ ="V1";

  }else if( tag=="V02" ) {

    tag_cef3_ ="V1";
    tag_bgo_ ="V0";
    tag_align_ ="V2";

  }else if( tag=="dev" ) {

    tag_cef3_ ="V1";
    tag_bgo_ ="V0";
    tag_align_ ="dev";

  } else {
    std::cout << "[TagHelper] :: Tag " << tag << " does not exist. Exiting" << std::endl;
    exit(12345);

  }
  
}




std::string TagHelper::getCeF3FileName() const {

  std::string fileName = "CeF3Calibration/constants_" + tag_cef3_ + ".txt";

  return fileName;

}


std::string TagHelper::getBGOFileName() const {

  std::string fileName = "BGOCalibration/constants_" + tag_bgo_ + ".txt";

  return fileName;

}


std::string TagHelper::getAlignmentFileName() const {

  std::string fileName = "Alignment/offsets_" + tag_align_ + ".txt";

  return fileName;

}
