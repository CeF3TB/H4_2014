#ifndef AlignmentOfficer_h
#define AlignmentOfficer_h

#include <string>
#include <vector>
#include <map>



class AlignmentOfficer {

 public:

  AlignmentOfficer( const std::string& fileName );
  ~AlignmentOfficer() {};

  void setOffsetsFile( const std::string& fileName );

  void fix( const std::string& name, int n, float *x );

  float getOffset( const std::string& name );



 private:

  std::map<std::string, float>  offsets_;
  std::string alignmentFile_;

};



#endif
