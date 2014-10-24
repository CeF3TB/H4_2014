#ifndef CalibrationUtility_h
#define CalibrationUtility_h


#include <string>



class CalibrationUtility {


 public:

  CalibrationUtility( const std::string& tag );
  ~CalibrationUtility();

  void setTag( const std::string& tag );

  std::string getCeF3FileName() const;
  std::string getBGOFileName() const;

  std::string getTag()     const { return tag_; };
  std::string getCeF3Tag() const { return tag_cef3_; };
  std::string getBGOTag()  const { return tag_bgo_; };



 private:

  std::string tag_;
  std::string tag_cef3_;
  std::string tag_bgo_;


};



#endif
