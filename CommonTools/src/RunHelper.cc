#include "../interface/RunHelper.h"

#include <iostream>
#include "TString.h"


void RunHelper::getBeamPosition( const std::string& runName, float& beamX, float& beamY ) {

  TString runName_tstr(runName.c_str());

  float xTable=0.;
  float yTable=0.;

  beamX = xTable;
  beamY = yTable;

}


