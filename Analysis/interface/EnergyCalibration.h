#ifndef EnergyCalibration_h
#define EnergyCalibration_h


#include <vector>
#include <string>


class EnergyCalibration {

 public:

  EnergyCalibration( const std::string& fileName );
  ~EnergyCalibration();

  void setCalibrationFile( const std::string& fileName );

  void applyCalibration( std::vector<float>& raw ) const;



 private:

  std::string calibFile_;
  std::vector<float> calibConstants_;


};



#endif
