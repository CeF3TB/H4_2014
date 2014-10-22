#include "../interface/HodoCluster.h"
#include <iostream>




void HodoCluster::addFibre( int i ) {

  float xmin = 0.5*(hodoTotFibres_-1)*fibreWidth_;
  float thisPos = fibreWidth_*((float)i) - xmin;

  if( size_==0 ) {

    pos_ += thisPos;
    size_+=1;

  } else {

    pos_ *= (float)size_;
    pos_ += thisPos;
    size_+=1;
    pos_ /= size_;

  }

}

