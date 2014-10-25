#include "../interface/AlignmentOfficer.h"

#include <fstream>
#include <iostream>



AlignmentOfficer::AlignmentOfficer( const std::string& fileName ) {

  this->setOffsetsFile(fileName);

}



void AlignmentOfficer::setOffsetsFile( const std::string& fileName ) {

  offsets_.clear();

  alignmentFile_ = fileName;

  std::ifstream ifs(fileName.c_str());

  while( ifs.good() ) {

    std::string name;
    float offset;
    ifs >> name >> offset;

    offsets_[name] = offset;

  }

  std::cout << "-> Loaded alignment offsets from file: " << fileName << std::endl;

}



float AlignmentOfficer::getOffset( const std::string& name ) {

 return offsets_[name];

}

 

void AlignmentOfficer::fix( const std::string& name, int n, float* x ) {

  for( unsigned i=0; i<n; ++i ) 
    x[i] += offsets_[name];

}
