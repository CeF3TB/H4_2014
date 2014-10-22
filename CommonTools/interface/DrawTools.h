#ifndef DrawTools_h
#define DrawTools_h


#include "TStyle.h"
#include "TPaveText.h"



class DrawTools {

 public:

  static TStyle* setStyle();

  static TPaveText* getLabelTop( const std::string& text="H4 Test Beam 2014" );
  static TPaveText* getLabelRun( const std::string& runName, bool top=true );


};


#endif
