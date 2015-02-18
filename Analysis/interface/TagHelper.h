#ifndef TagHelper_h
#define TagHelper_h


#include <string>



class TagHelper {


 public:

  TagHelper( const std::string& tag, const std::string& Energy );
  ~TagHelper();

  void setTag( const std::string& tag, const std::string& Energy );

  std::string getCeF3FileName() const;
  std::string getBGOFileName() const;
  std::string getAlignmentFileName() const;

  std::string getTag()     const { return tag_; };
  std::string getCeF3Tag() const { return tag_cef3_; };
  std::string getBGOTag()  const { return tag_bgo_; };
  std::string getAlignmentTag()  const { return tag_align_; };



 private:

  std::string tag_;
  std::string tag_cef3_;
  std::string tag_bgo_;
  std::string tag_align_;


};



#endif
