#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

#include "channelInfo.h"

#include "CommonTools/interface/RunHelper.h"



int main( int argc, char* argv[] ) {


   if( argc<2 ) {
     std::cout << "ERROR. You need to specify the name of the run you want to process." << std::endl;  
     exit(1);
   }

   std::string runName = "";
   std::string tag = "V00";

   if( argc>1 ) {

     std::string runName_str(argv[1]);
     runName = runName_str;

     if( argc>2 ) {
       std::string tag_str(argv[2]);
       tag = tag_str;
     }

   } else {

     std::cout << "Usage:" << std::endl;
     std::cout << "./makeAnalysisTree [runName]" << std::endl;
     exit(12345);

   }

   std::string fileName = "data/run_" + runName + ".root";
   TFile* file = TFile::Open(fileName.c_str());
   if( file==0 ) {
     std::cout << "ERROR! Din't find file " << fileName << std::endl;
     std::cout << "Exiting." << std::endl;
     exit(11);
   }
   TTree* tree = (TTree*)file->Get("outputTree");




   // Declaration of leaf types
   std::vector<float>   *digi_gr0_ch0_X;
   std::vector<float>   *digi_gr0_ch0_Y;
   std::vector<float>   *digi_gr0_ch0_max_amplitude;
   std::vector<float>   *digi_gr0_ch0_pedestal;
   std::vector<float>   *digi_gr0_ch0_pedestal_rms;
   std::vector<float>   *digi_gr0_ch0_time_at_frac30;
   std::vector<float>   *digi_gr0_ch0_time_at_frac50;
   std::vector<float>   *digi_gr0_ch0_time_at_max;
   std::vector<float>   *digi_gr0_ch1_X;
   std::vector<float>   *digi_gr0_ch1_Y;
   std::vector<float>   *digi_gr0_ch1_max_amplitude;
   std::vector<float>   *digi_gr0_ch1_pedestal;
   std::vector<float>   *digi_gr0_ch1_pedestal_rms;
   std::vector<float>   *digi_gr0_ch1_time_at_frac30;
   std::vector<float>   *digi_gr0_ch1_time_at_frac50;
   std::vector<float>   *digi_gr0_ch1_time_at_max;
   std::vector<float>   *digi_gr0_ch2_X;
   std::vector<float>   *digi_gr0_ch2_Y;
   std::vector<float>   *digi_gr0_ch2_max_amplitude;
   std::vector<float>   *digi_gr0_ch2_pedestal;
   std::vector<float>   *digi_gr0_ch2_pedestal_rms;
   std::vector<float>   *digi_gr0_ch2_time_at_frac30;
   std::vector<float>   *digi_gr0_ch2_time_at_frac50;
   std::vector<float>   *digi_gr0_ch2_time_at_max;
   std::vector<float>   *digi_gr0_ch3_X;
   std::vector<float>   *digi_gr0_ch3_Y;
   std::vector<float>   *digi_gr0_ch3_max_amplitude;
   std::vector<float>   *digi_gr0_ch3_pedestal;
   std::vector<float>   *digi_gr0_ch3_pedestal_rms;
   std::vector<float>   *digi_gr0_ch3_time_at_frac30;
   std::vector<float>   *digi_gr0_ch3_time_at_frac50;
   std::vector<float>   *digi_gr0_ch3_time_at_max;
   std::vector<float>   *digi_gr0_ch4_X;
   std::vector<float>   *digi_gr0_ch4_Y;
   std::vector<float>   *digi_gr0_ch4_max_amplitude;
   std::vector<float>   *digi_gr0_ch4_pedestal;
   std::vector<float>   *digi_gr0_ch4_pedestal_rms;
   std::vector<float>   *digi_gr0_ch4_time_at_frac30;
   std::vector<float>   *digi_gr0_ch4_time_at_frac50;
   std::vector<float>   *digi_gr0_ch4_time_at_max;
   std::vector<float>   *digi_gr0_ch5_X;
   std::vector<float>   *digi_gr0_ch5_Y;
   std::vector<float>   *digi_gr0_ch5_max_amplitude;
   std::vector<float>   *digi_gr0_ch5_pedestal;
   std::vector<float>   *digi_gr0_ch5_pedestal_rms;
   std::vector<float>   *digi_gr0_ch5_time_at_frac30;
   std::vector<float>   *digi_gr0_ch5_time_at_frac50;
   std::vector<float>   *digi_gr0_ch5_time_at_max;
   std::vector<float>   *digi_gr0_ch6_X;
   std::vector<float>   *digi_gr0_ch6_Y;
   std::vector<float>   *digi_gr0_ch6_max_amplitude;
   std::vector<float>   *digi_gr0_ch6_pedestal;
   std::vector<float>   *digi_gr0_ch6_pedestal_rms;
   std::vector<float>   *digi_gr0_ch6_time_at_frac30;
   std::vector<float>   *digi_gr0_ch6_time_at_frac50;
   std::vector<float>   *digi_gr0_ch6_time_at_max;
   std::vector<float>   *digi_gr0_ch7_X;
   std::vector<float>   *digi_gr0_ch7_Y;
   std::vector<float>   *digi_gr0_ch7_max_amplitude;
   std::vector<float>   *digi_gr0_ch7_pedestal;
   std::vector<float>   *digi_gr0_ch7_pedestal_rms;
   std::vector<float>   *digi_gr0_ch7_time_at_frac30;
   std::vector<float>   *digi_gr0_ch7_time_at_frac50;
   std::vector<float>   *digi_gr0_ch7_time_at_max;
   std::vector<float>   *digi_gr0_ch8_X;
   std::vector<float>   *digi_gr0_ch8_Y;
   std::vector<float>   *digi_gr0_ch8_max_amplitude;
   std::vector<float>   *digi_gr0_ch8_pedestal;
   std::vector<float>   *digi_gr0_ch8_pedestal_rms;
   std::vector<float>   *digi_gr0_ch8_time_at_frac30;
   std::vector<float>   *digi_gr0_ch8_time_at_frac50;
   std::vector<float>   *digi_gr0_ch8_time_at_max;
   std::vector<float>   *ADC_board_11301_0;
   std::vector<float>   *ADC_board_11301_1;
   std::vector<float>   *ADC_board_11301_2;
   std::vector<float>   *ADC_board_11301_3;
   std::vector<float>   *ADC_board_11301_4;
   std::vector<float>   *ADC_board_11301_5;
   std::vector<float>   *ADC_board_11301_6;
   std::vector<float>   *ADC_board_11301_7;
   std::vector<float>   *ADC_board_6301_0;
   std::vector<float>   *ADC_board_6301_1;
   std::vector<float>   *ADC_board_6301_10;
   std::vector<float>   *ADC_board_6301_11;
   std::vector<float>   *ADC_board_6301_12;
   std::vector<float>   *ADC_board_6301_13;
   std::vector<float>   *ADC_board_6301_14;
   std::vector<float>   *ADC_board_6301_15;
   std::vector<float>   *ADC_board_6301_16;
   std::vector<float>   *ADC_board_6301_17;
   std::vector<float>   *ADC_board_6301_18;
   std::vector<float>   *ADC_board_6301_19;
   std::vector<float>   *ADC_board_6301_2;
   std::vector<float>   *ADC_board_6301_20;
   std::vector<float>   *ADC_board_6301_21;
   std::vector<float>   *ADC_board_6301_22;
   std::vector<float>   *ADC_board_6301_23;
   std::vector<float>   *ADC_board_6301_24;
   std::vector<float>   *ADC_board_6301_25;
   std::vector<float>   *ADC_board_6301_26;
   std::vector<float>   *ADC_board_6301_27;
   std::vector<float>   *ADC_board_6301_28;
   std::vector<float>   *ADC_board_6301_29;
   std::vector<float>   *ADC_board_6301_3;
   std::vector<float>   *ADC_board_6301_30;
   std::vector<float>   *ADC_board_6301_31;
   std::vector<float>   *ADC_board_6301_4;
   std::vector<float>   *ADC_board_6301_5;
   std::vector<float>   *ADC_board_6301_6;
   std::vector<float>   *ADC_board_6301_7;
   std::vector<float>   *ADC_board_6301_8;
   std::vector<float>   *ADC_board_6301_9;
   std::vector<float>   *MULTILINE_time0;
   std::vector<float>   *MULTILINE_time1;
   std::vector<float>   *MULTILINE_time2;
   std::vector<float>   *TDCinputTime1;
   std::vector<float>   *TDCinputTime2;
   std::vector<float>   *TDCinputTime3;
   std::vector<float>   *TDCinputTime4;
   std::vector<float>   *TDCrecoPos_X;
   std::vector<float>   *TDCrecoX;
   std::vector<float>   *TDCrecoY;
   std::vector<float>   *beamPositionSmallX;
   std::vector<float>   *beamPositionSmallY;
   std::vector<float>   *beamPositionX;
   std::vector<float>   *beamPositionX1;
   std::vector<float>   *beamPositionX2;
   std::vector<float>   *beamPositionY;
   std::vector<float>   *beamPositionY1;
   std::vector<float>   *beamPositionY2;
   std::vector<float>   *beamProfileSmallX;
   std::vector<float>   *beamProfileSmallY;
   std::vector<float>   *beamProfileX1;
   std::vector<float>   *beamProfileX2;
   std::vector<float>   *beamProfileY1;
   std::vector<float>   *beamProfileY2;
   std::vector<float>   *deltaTime10;
   std::vector<float>   *deltaTime20;
   std::vector<float>   *deltaTime21;
   std::vector<float>   *fractionTakenTrigHisto;
   std::vector<float>   *nFibersOnSmallX;
   std::vector<float>   *nFibersOnSmallY;
   std::vector<float>   *nFibersOnX1;
   std::vector<float>   *nFibersOnX2;
   std::vector<float>   *nFibersOnY1;
   std::vector<float>   *nFibersOnY2;
   std::vector<float>   *nTotalEvts;

   // List of branches
   TBranch        *b_digi_gr0_ch0_X;   //!
   TBranch        *b_digi_gr0_ch0_Y;   //!
   TBranch        *b_digi_gr0_ch0_max_amplitude;   //!
   TBranch        *b_digi_gr0_ch0_pedestal;   //!
   TBranch        *b_digi_gr0_ch0_pedestal_rms;   //!
   TBranch        *b_digi_gr0_ch0_time_at_frac30;   //!
   TBranch        *b_digi_gr0_ch0_time_at_frac50;   //!
   TBranch        *b_digi_gr0_ch0_time_at_max;   //!
   TBranch        *b_digi_gr0_ch1_X;   //!
   TBranch        *b_digi_gr0_ch1_Y;   //!
   TBranch        *b_digi_gr0_ch1_max_amplitude;   //!
   TBranch        *b_digi_gr0_ch1_pedestal;   //!
   TBranch        *b_digi_gr0_ch1_pedestal_rms;   //!
   TBranch        *b_digi_gr0_ch1_time_at_frac30;   //!
   TBranch        *b_digi_gr0_ch1_time_at_frac50;   //!
   TBranch        *b_digi_gr0_ch1_time_at_max;   //!
   TBranch        *b_digi_gr0_ch2_X;   //!
   TBranch        *b_digi_gr0_ch2_Y;   //!
   TBranch        *b_digi_gr0_ch2_max_amplitude;   //!
   TBranch        *b_digi_gr0_ch2_pedestal;   //!
   TBranch        *b_digi_gr0_ch2_pedestal_rms;   //!
   TBranch        *b_digi_gr0_ch2_time_at_frac30;   //!
   TBranch        *b_digi_gr0_ch2_time_at_frac50;   //!
   TBranch        *b_digi_gr0_ch2_time_at_max;   //!
   TBranch        *b_digi_gr0_ch3_X;   //!
   TBranch        *b_digi_gr0_ch3_Y;   //!
   TBranch        *b_digi_gr0_ch3_max_amplitude;   //!
   TBranch        *b_digi_gr0_ch3_pedestal;   //!
   TBranch        *b_digi_gr0_ch3_pedestal_rms;   //!
   TBranch        *b_digi_gr0_ch3_time_at_frac30;   //!
   TBranch        *b_digi_gr0_ch3_time_at_frac50;   //!
   TBranch        *b_digi_gr0_ch3_time_at_max;   //!
   TBranch        *b_digi_gr0_ch4_X;   //!
   TBranch        *b_digi_gr0_ch4_Y;   //!
   TBranch        *b_digi_gr0_ch4_max_amplitude;   //!
   TBranch        *b_digi_gr0_ch4_pedestal;   //!
   TBranch        *b_digi_gr0_ch4_pedestal_rms;   //!
   TBranch        *b_digi_gr0_ch4_time_at_frac30;   //!
   TBranch        *b_digi_gr0_ch4_time_at_frac50;   //!
   TBranch        *b_digi_gr0_ch4_time_at_max;   //!
   TBranch        *b_digi_gr0_ch5_X;   //!
   TBranch        *b_digi_gr0_ch5_Y;   //!
   TBranch        *b_digi_gr0_ch5_max_amplitude;   //!
   TBranch        *b_digi_gr0_ch5_pedestal;   //!
   TBranch        *b_digi_gr0_ch5_pedestal_rms;   //!
   TBranch        *b_digi_gr0_ch5_time_at_frac30;   //!
   TBranch        *b_digi_gr0_ch5_time_at_frac50;   //!
   TBranch        *b_digi_gr0_ch5_time_at_max;   //!
   TBranch        *b_digi_gr0_ch6_X;   //!
   TBranch        *b_digi_gr0_ch6_Y;   //!
   TBranch        *b_digi_gr0_ch6_max_amplitude;   //!
   TBranch        *b_digi_gr0_ch6_pedestal;   //!
   TBranch        *b_digi_gr0_ch6_pedestal_rms;   //!
   TBranch        *b_digi_gr0_ch6_time_at_frac30;   //!
   TBranch        *b_digi_gr0_ch6_time_at_frac50;   //!
   TBranch        *b_digi_gr0_ch6_time_at_max;   //!
   TBranch        *b_digi_gr0_ch7_X;   //!
   TBranch        *b_digi_gr0_ch7_Y;   //!
   TBranch        *b_digi_gr0_ch7_max_amplitude;   //!
   TBranch        *b_digi_gr0_ch7_pedestal;   //!
   TBranch        *b_digi_gr0_ch7_pedestal_rms;   //!
   TBranch        *b_digi_gr0_ch7_time_at_frac30;   //!
   TBranch        *b_digi_gr0_ch7_time_at_frac50;   //!
   TBranch        *b_digi_gr0_ch7_time_at_max;   //!
   TBranch        *b_digi_gr0_ch8_X;   //!
   TBranch        *b_digi_gr0_ch8_Y;   //!
   TBranch        *b_digi_gr0_ch8_max_amplitude;   //!
   TBranch        *b_digi_gr0_ch8_pedestal;   //!
   TBranch        *b_digi_gr0_ch8_pedestal_rms;   //!
   TBranch        *b_digi_gr0_ch8_time_at_frac30;   //!
   TBranch        *b_digi_gr0_ch8_time_at_frac50;   //!
   TBranch        *b_digi_gr0_ch8_time_at_max;   //!
   TBranch        *b_ADC_board_11301_0;   //!
   TBranch        *b_ADC_board_11301_1;   //!
   TBranch        *b_ADC_board_11301_2;   //!
   TBranch        *b_ADC_board_11301_3;   //!
   TBranch        *b_ADC_board_11301_4;   //!
   TBranch        *b_ADC_board_11301_5;   //!
   TBranch        *b_ADC_board_11301_6;   //!
   TBranch        *b_ADC_board_11301_7;   //!
   TBranch        *b_ADC_board_6301_0;   //!
   TBranch        *b_ADC_board_6301_1;   //!
   TBranch        *b_ADC_board_6301_10;   //!
   TBranch        *b_ADC_board_6301_11;   //!
   TBranch        *b_ADC_board_6301_12;   //!
   TBranch        *b_ADC_board_6301_13;   //!
   TBranch        *b_ADC_board_6301_14;   //!
   TBranch        *b_ADC_board_6301_15;   //!
   TBranch        *b_ADC_board_6301_16;   //!
   TBranch        *b_ADC_board_6301_17;   //!
   TBranch        *b_ADC_board_6301_18;   //!
   TBranch        *b_ADC_board_6301_19;   //!
   TBranch        *b_ADC_board_6301_2;   //!
   TBranch        *b_ADC_board_6301_20;   //!
   TBranch        *b_ADC_board_6301_21;   //!
   TBranch        *b_ADC_board_6301_22;   //!
   TBranch        *b_ADC_board_6301_23;   //!
   TBranch        *b_ADC_board_6301_24;   //!
   TBranch        *b_ADC_board_6301_25;   //!
   TBranch        *b_ADC_board_6301_26;   //!
   TBranch        *b_ADC_board_6301_27;   //!
   TBranch        *b_ADC_board_6301_28;   //!
   TBranch        *b_ADC_board_6301_29;   //!
   TBranch        *b_ADC_board_6301_3;   //!
   TBranch        *b_ADC_board_6301_30;   //!
   TBranch        *b_ADC_board_6301_31;   //!
   TBranch        *b_ADC_board_6301_4;   //!
   TBranch        *b_ADC_board_6301_5;   //!
   TBranch        *b_ADC_board_6301_6;   //!
   TBranch        *b_ADC_board_6301_7;   //!
   TBranch        *b_ADC_board_6301_8;   //!
   TBranch        *b_ADC_board_6301_9;   //!
   TBranch        *b_MULTILINE_time0;   //!
   TBranch        *b_MULTILINE_time1;   //!
   TBranch        *b_MULTILINE_time2;   //!
   TBranch        *b_TDCinputTime1;   //!
   TBranch        *b_TDCinputTime2;   //!
   TBranch        *b_TDCinputTime3;   //!
   TBranch        *b_TDCinputTime4;   //!
   TBranch        *b_TDCrecoPos_X;   //!
   TBranch        *b_TDCrecoX;   //!
   TBranch        *b_TDCrecoY;   //!
   TBranch        *b_beamPositionSmallX;   //!
   TBranch        *b_beamPositionSmallY;   //!
   TBranch        *b_beamPositionX;   //!
   TBranch        *b_beamPositionX1;   //!
   TBranch        *b_beamPositionX2;   //!
   TBranch        *b_beamPositionY;   //!
   TBranch        *b_beamPositionY1;   //!
   TBranch        *b_beamPositionY2;   //!
   TBranch        *b_beamProfileSmallX;   //!
   TBranch        *b_beamProfileSmallY;   //!
   TBranch        *b_beamProfileX1;   //!
   TBranch        *b_beamProfileX2;   //!
   TBranch        *b_beamProfileY1;   //!
   TBranch        *b_beamProfileY2;   //!
   TBranch        *b_deltaTime10;   //!
   TBranch        *b_deltaTime20;   //!
   TBranch        *b_deltaTime21;   //!
   TBranch        *b_fractionTakenTrigHisto;   //!
   TBranch        *b_nFibersOnSmallX;   //!
   TBranch        *b_nFibersOnSmallY;   //!
   TBranch        *b_nFibersOnX1;   //!
   TBranch        *b_nFibersOnX2;   //!
   TBranch        *b_nFibersOnY1;   //!
   TBranch        *b_nFibersOnY2;   //!
   TBranch        *b_nTotalEvts;   //!

   // Set object pointer
   digi_gr0_ch0_X = 0;
   digi_gr0_ch0_Y = 0;
   digi_gr0_ch0_max_amplitude = 0;
   digi_gr0_ch0_pedestal = 0;
   digi_gr0_ch0_pedestal_rms = 0;
   digi_gr0_ch0_time_at_frac30 = 0;
   digi_gr0_ch0_time_at_frac50 = 0;
   digi_gr0_ch0_time_at_max = 0;
   digi_gr0_ch1_X = 0;
   digi_gr0_ch1_Y = 0;
   digi_gr0_ch1_max_amplitude = 0;
   digi_gr0_ch1_pedestal = 0;
   digi_gr0_ch1_pedestal_rms = 0;
   digi_gr0_ch1_time_at_frac30 = 0;
   digi_gr0_ch1_time_at_frac50 = 0;
   digi_gr0_ch1_time_at_max = 0;
   digi_gr0_ch2_X = 0;
   digi_gr0_ch2_Y = 0;
   digi_gr0_ch2_max_amplitude = 0;
   digi_gr0_ch2_pedestal = 0;
   digi_gr0_ch2_pedestal_rms = 0;
   digi_gr0_ch2_time_at_frac30 = 0;
   digi_gr0_ch2_time_at_frac50 = 0;
   digi_gr0_ch2_time_at_max = 0;
   digi_gr0_ch3_X = 0;
   digi_gr0_ch3_Y = 0;
   digi_gr0_ch3_max_amplitude = 0;
   digi_gr0_ch3_pedestal = 0;
   digi_gr0_ch3_pedestal_rms = 0;
   digi_gr0_ch3_time_at_frac30 = 0;
   digi_gr0_ch3_time_at_frac50 = 0;
   digi_gr0_ch3_time_at_max = 0;
   digi_gr0_ch4_X = 0;
   digi_gr0_ch4_Y = 0;
   digi_gr0_ch4_max_amplitude = 0;
   digi_gr0_ch4_pedestal = 0;
   digi_gr0_ch4_pedestal_rms = 0;
   digi_gr0_ch4_time_at_frac30 = 0;
   digi_gr0_ch4_time_at_frac50 = 0;
   digi_gr0_ch4_time_at_max = 0;
   digi_gr0_ch5_X = 0;
   digi_gr0_ch5_Y = 0;
   digi_gr0_ch5_max_amplitude = 0;
   digi_gr0_ch5_pedestal = 0;
   digi_gr0_ch5_pedestal_rms = 0;
   digi_gr0_ch5_time_at_frac30 = 0;
   digi_gr0_ch5_time_at_frac50 = 0;
   digi_gr0_ch5_time_at_max = 0;
   digi_gr0_ch6_X = 0;
   digi_gr0_ch6_Y = 0;
   digi_gr0_ch6_max_amplitude = 0;
   digi_gr0_ch6_pedestal = 0;
   digi_gr0_ch6_pedestal_rms = 0;
   digi_gr0_ch6_time_at_frac30 = 0;
   digi_gr0_ch6_time_at_frac50 = 0;
   digi_gr0_ch6_time_at_max = 0;
   digi_gr0_ch7_X = 0;
   digi_gr0_ch7_Y = 0;
   digi_gr0_ch7_max_amplitude = 0;
   digi_gr0_ch7_pedestal = 0;
   digi_gr0_ch7_pedestal_rms = 0;
   digi_gr0_ch7_time_at_frac30 = 0;
   digi_gr0_ch7_time_at_frac50 = 0;
   digi_gr0_ch7_time_at_max = 0;
   digi_gr0_ch8_X = 0;
   digi_gr0_ch8_Y = 0;
   digi_gr0_ch8_max_amplitude = 0;
   digi_gr0_ch8_pedestal = 0;
   digi_gr0_ch8_pedestal_rms = 0;
   digi_gr0_ch8_time_at_frac30 = 0;
   digi_gr0_ch8_time_at_frac50 = 0;
   digi_gr0_ch8_time_at_max = 0;
   ADC_board_11301_0 = 0;
   ADC_board_11301_1 = 0;
   ADC_board_11301_2 = 0;
   ADC_board_11301_3 = 0;
   ADC_board_11301_4 = 0;
   ADC_board_11301_5 = 0;
   ADC_board_11301_6 = 0;
   ADC_board_11301_7 = 0;
   ADC_board_6301_0 = 0;
   ADC_board_6301_1 = 0;
   ADC_board_6301_10 = 0;
   ADC_board_6301_11 = 0;
   ADC_board_6301_12 = 0;
   ADC_board_6301_13 = 0;
   ADC_board_6301_14 = 0;
   ADC_board_6301_15 = 0;
   ADC_board_6301_16 = 0;
   ADC_board_6301_17 = 0;
   ADC_board_6301_18 = 0;
   ADC_board_6301_19 = 0;
   ADC_board_6301_2 = 0;
   ADC_board_6301_20 = 0;
   ADC_board_6301_21 = 0;
   ADC_board_6301_22 = 0;
   ADC_board_6301_23 = 0;
   ADC_board_6301_24 = 0;
   ADC_board_6301_25 = 0;
   ADC_board_6301_26 = 0;
   ADC_board_6301_27 = 0;
   ADC_board_6301_28 = 0;
   ADC_board_6301_29 = 0;
   ADC_board_6301_3 = 0;
   ADC_board_6301_30 = 0;
   ADC_board_6301_31 = 0;
   ADC_board_6301_4 = 0;
   ADC_board_6301_5 = 0;
   ADC_board_6301_6 = 0;
   ADC_board_6301_7 = 0;
   ADC_board_6301_8 = 0;
   ADC_board_6301_9 = 0;
   MULTILINE_time0 = 0;
   MULTILINE_time1 = 0;
   MULTILINE_time2 = 0;
   TDCinputTime1 = 0;
   TDCinputTime2 = 0;
   TDCinputTime3 = 0;
   TDCinputTime4 = 0;
   TDCrecoPos_X = 0;
   TDCrecoX = 0;
   TDCrecoY = 0;
   beamPositionSmallX = 0;
   beamPositionSmallY = 0;
   beamPositionX = 0;
   beamPositionX1 = 0;
   beamPositionX2 = 0;
   beamPositionY = 0;
   beamPositionY1 = 0;
   beamPositionY2 = 0;
   beamProfileSmallX = 0;
   beamProfileSmallY = 0;
   beamProfileX1 = 0;
   beamProfileX2 = 0;
   beamProfileY1 = 0;
   beamProfileY2 = 0;
   deltaTime10 = 0;
   deltaTime20 = 0;
   deltaTime21 = 0;
   fractionTakenTrigHisto = 0;
   nFibersOnSmallX = 0;
   nFibersOnSmallY = 0;
   nFibersOnX1 = 0;
   nFibersOnX2 = 0;
   nFibersOnY1 = 0;
   nFibersOnY2 = 0;
   nTotalEvts = 0;

   tree->SetBranchAddress("digi_gr0_ch0_X", &digi_gr0_ch0_X, &b_digi_gr0_ch0_X);
   tree->SetBranchAddress("digi_gr0_ch0_Y", &digi_gr0_ch0_Y, &b_digi_gr0_ch0_Y);
   tree->SetBranchAddress("digi_gr0_ch0_max_amplitude", &digi_gr0_ch0_max_amplitude, &b_digi_gr0_ch0_max_amplitude);
   tree->SetBranchAddress("digi_gr0_ch0_pedestal", &digi_gr0_ch0_pedestal, &b_digi_gr0_ch0_pedestal);
   tree->SetBranchAddress("digi_gr0_ch0_pedestal_rms", &digi_gr0_ch0_pedestal_rms, &b_digi_gr0_ch0_pedestal_rms);
   tree->SetBranchAddress("digi_gr0_ch0_time_at_frac30", &digi_gr0_ch0_time_at_frac30, &b_digi_gr0_ch0_time_at_frac30);
   tree->SetBranchAddress("digi_gr0_ch0_time_at_frac50", &digi_gr0_ch0_time_at_frac50, &b_digi_gr0_ch0_time_at_frac50);
   tree->SetBranchAddress("digi_gr0_ch0_time_at_max", &digi_gr0_ch0_time_at_max, &b_digi_gr0_ch0_time_at_max);
   tree->SetBranchAddress("digi_gr0_ch1_X", &digi_gr0_ch1_X, &b_digi_gr0_ch1_X);
   tree->SetBranchAddress("digi_gr0_ch1_Y", &digi_gr0_ch1_Y, &b_digi_gr0_ch1_Y);
   tree->SetBranchAddress("digi_gr0_ch1_max_amplitude", &digi_gr0_ch1_max_amplitude, &b_digi_gr0_ch1_max_amplitude);
   tree->SetBranchAddress("digi_gr0_ch1_pedestal", &digi_gr0_ch1_pedestal, &b_digi_gr0_ch1_pedestal);
   tree->SetBranchAddress("digi_gr0_ch1_pedestal_rms", &digi_gr0_ch1_pedestal_rms, &b_digi_gr0_ch1_pedestal_rms);
   tree->SetBranchAddress("digi_gr0_ch1_time_at_frac30", &digi_gr0_ch1_time_at_frac30, &b_digi_gr0_ch1_time_at_frac30);
   tree->SetBranchAddress("digi_gr0_ch1_time_at_frac50", &digi_gr0_ch1_time_at_frac50, &b_digi_gr0_ch1_time_at_frac50);
   tree->SetBranchAddress("digi_gr0_ch1_time_at_max", &digi_gr0_ch1_time_at_max, &b_digi_gr0_ch1_time_at_max);
   tree->SetBranchAddress("digi_gr0_ch2_X", &digi_gr0_ch2_X, &b_digi_gr0_ch2_X);
   tree->SetBranchAddress("digi_gr0_ch2_Y", &digi_gr0_ch2_Y, &b_digi_gr0_ch2_Y);
   tree->SetBranchAddress("digi_gr0_ch2_max_amplitude", &digi_gr0_ch2_max_amplitude, &b_digi_gr0_ch2_max_amplitude);
   tree->SetBranchAddress("digi_gr0_ch2_pedestal", &digi_gr0_ch2_pedestal, &b_digi_gr0_ch2_pedestal);
   tree->SetBranchAddress("digi_gr0_ch2_pedestal_rms", &digi_gr0_ch2_pedestal_rms, &b_digi_gr0_ch2_pedestal_rms);
   tree->SetBranchAddress("digi_gr0_ch2_time_at_frac30", &digi_gr0_ch2_time_at_frac30, &b_digi_gr0_ch2_time_at_frac30);
   tree->SetBranchAddress("digi_gr0_ch2_time_at_frac50", &digi_gr0_ch2_time_at_frac50, &b_digi_gr0_ch2_time_at_frac50);
   tree->SetBranchAddress("digi_gr0_ch2_time_at_max", &digi_gr0_ch2_time_at_max, &b_digi_gr0_ch2_time_at_max);
   tree->SetBranchAddress("digi_gr0_ch3_X", &digi_gr0_ch3_X, &b_digi_gr0_ch3_X);
   tree->SetBranchAddress("digi_gr0_ch3_Y", &digi_gr0_ch3_Y, &b_digi_gr0_ch3_Y);
   tree->SetBranchAddress("digi_gr0_ch3_max_amplitude", &digi_gr0_ch3_max_amplitude, &b_digi_gr0_ch3_max_amplitude);
   tree->SetBranchAddress("digi_gr0_ch3_pedestal", &digi_gr0_ch3_pedestal, &b_digi_gr0_ch3_pedestal);
   tree->SetBranchAddress("digi_gr0_ch3_pedestal_rms", &digi_gr0_ch3_pedestal_rms, &b_digi_gr0_ch3_pedestal_rms);
   tree->SetBranchAddress("digi_gr0_ch3_time_at_frac30", &digi_gr0_ch3_time_at_frac30, &b_digi_gr0_ch3_time_at_frac30);
   tree->SetBranchAddress("digi_gr0_ch3_time_at_frac50", &digi_gr0_ch3_time_at_frac50, &b_digi_gr0_ch3_time_at_frac50);
   tree->SetBranchAddress("digi_gr0_ch3_time_at_max", &digi_gr0_ch3_time_at_max, &b_digi_gr0_ch3_time_at_max);
   tree->SetBranchAddress("digi_gr0_ch4_X", &digi_gr0_ch4_X, &b_digi_gr0_ch4_X);
   tree->SetBranchAddress("digi_gr0_ch4_Y", &digi_gr0_ch4_Y, &b_digi_gr0_ch4_Y);
   tree->SetBranchAddress("digi_gr0_ch4_max_amplitude", &digi_gr0_ch4_max_amplitude, &b_digi_gr0_ch4_max_amplitude);
   tree->SetBranchAddress("digi_gr0_ch4_pedestal", &digi_gr0_ch4_pedestal, &b_digi_gr0_ch4_pedestal);
   tree->SetBranchAddress("digi_gr0_ch4_pedestal_rms", &digi_gr0_ch4_pedestal_rms, &b_digi_gr0_ch4_pedestal_rms);
   tree->SetBranchAddress("digi_gr0_ch4_time_at_frac30", &digi_gr0_ch4_time_at_frac30, &b_digi_gr0_ch4_time_at_frac30);
   tree->SetBranchAddress("digi_gr0_ch4_time_at_frac50", &digi_gr0_ch4_time_at_frac50, &b_digi_gr0_ch4_time_at_frac50);
   tree->SetBranchAddress("digi_gr0_ch4_time_at_max", &digi_gr0_ch4_time_at_max, &b_digi_gr0_ch4_time_at_max);
   tree->SetBranchAddress("digi_gr0_ch5_X", &digi_gr0_ch5_X, &b_digi_gr0_ch5_X);
   tree->SetBranchAddress("digi_gr0_ch5_Y", &digi_gr0_ch5_Y, &b_digi_gr0_ch5_Y);
   tree->SetBranchAddress("digi_gr0_ch5_max_amplitude", &digi_gr0_ch5_max_amplitude, &b_digi_gr0_ch5_max_amplitude);
   tree->SetBranchAddress("digi_gr0_ch5_pedestal", &digi_gr0_ch5_pedestal, &b_digi_gr0_ch5_pedestal);
   tree->SetBranchAddress("digi_gr0_ch5_pedestal_rms", &digi_gr0_ch5_pedestal_rms, &b_digi_gr0_ch5_pedestal_rms);
   tree->SetBranchAddress("digi_gr0_ch5_time_at_frac30", &digi_gr0_ch5_time_at_frac30, &b_digi_gr0_ch5_time_at_frac30);
   tree->SetBranchAddress("digi_gr0_ch5_time_at_frac50", &digi_gr0_ch5_time_at_frac50, &b_digi_gr0_ch5_time_at_frac50);
   tree->SetBranchAddress("digi_gr0_ch5_time_at_max", &digi_gr0_ch5_time_at_max, &b_digi_gr0_ch5_time_at_max);
   tree->SetBranchAddress("digi_gr0_ch6_X", &digi_gr0_ch6_X, &b_digi_gr0_ch6_X);
   tree->SetBranchAddress("digi_gr0_ch6_Y", &digi_gr0_ch6_Y, &b_digi_gr0_ch6_Y);
   tree->SetBranchAddress("digi_gr0_ch6_max_amplitude", &digi_gr0_ch6_max_amplitude, &b_digi_gr0_ch6_max_amplitude);
   tree->SetBranchAddress("digi_gr0_ch6_pedestal", &digi_gr0_ch6_pedestal, &b_digi_gr0_ch6_pedestal);
   tree->SetBranchAddress("digi_gr0_ch6_pedestal_rms", &digi_gr0_ch6_pedestal_rms, &b_digi_gr0_ch6_pedestal_rms);
   tree->SetBranchAddress("digi_gr0_ch6_time_at_frac30", &digi_gr0_ch6_time_at_frac30, &b_digi_gr0_ch6_time_at_frac30);
   tree->SetBranchAddress("digi_gr0_ch6_time_at_frac50", &digi_gr0_ch6_time_at_frac50, &b_digi_gr0_ch6_time_at_frac50);
   tree->SetBranchAddress("digi_gr0_ch6_time_at_max", &digi_gr0_ch6_time_at_max, &b_digi_gr0_ch6_time_at_max);
   tree->SetBranchAddress("digi_gr0_ch7_X", &digi_gr0_ch7_X, &b_digi_gr0_ch7_X);
   tree->SetBranchAddress("digi_gr0_ch7_Y", &digi_gr0_ch7_Y, &b_digi_gr0_ch7_Y);
   tree->SetBranchAddress("digi_gr0_ch7_max_amplitude", &digi_gr0_ch7_max_amplitude, &b_digi_gr0_ch7_max_amplitude);
   tree->SetBranchAddress("digi_gr0_ch7_pedestal", &digi_gr0_ch7_pedestal, &b_digi_gr0_ch7_pedestal);
   tree->SetBranchAddress("digi_gr0_ch7_pedestal_rms", &digi_gr0_ch7_pedestal_rms, &b_digi_gr0_ch7_pedestal_rms);
   tree->SetBranchAddress("digi_gr0_ch7_time_at_frac30", &digi_gr0_ch7_time_at_frac30, &b_digi_gr0_ch7_time_at_frac30);
   tree->SetBranchAddress("digi_gr0_ch7_time_at_frac50", &digi_gr0_ch7_time_at_frac50, &b_digi_gr0_ch7_time_at_frac50);
   tree->SetBranchAddress("digi_gr0_ch7_time_at_max", &digi_gr0_ch7_time_at_max, &b_digi_gr0_ch7_time_at_max);
   tree->SetBranchAddress("digi_gr0_ch8_X", &digi_gr0_ch8_X, &b_digi_gr0_ch8_X);
   tree->SetBranchAddress("digi_gr0_ch8_Y", &digi_gr0_ch8_Y, &b_digi_gr0_ch8_Y);
   tree->SetBranchAddress("digi_gr0_ch8_max_amplitude", &digi_gr0_ch8_max_amplitude, &b_digi_gr0_ch8_max_amplitude);
   tree->SetBranchAddress("digi_gr0_ch8_pedestal", &digi_gr0_ch8_pedestal, &b_digi_gr0_ch8_pedestal);
   tree->SetBranchAddress("digi_gr0_ch8_pedestal_rms", &digi_gr0_ch8_pedestal_rms, &b_digi_gr0_ch8_pedestal_rms);
   tree->SetBranchAddress("digi_gr0_ch8_time_at_frac30", &digi_gr0_ch8_time_at_frac30, &b_digi_gr0_ch8_time_at_frac30);
   tree->SetBranchAddress("digi_gr0_ch8_time_at_frac50", &digi_gr0_ch8_time_at_frac50, &b_digi_gr0_ch8_time_at_frac50);
   tree->SetBranchAddress("digi_gr0_ch8_time_at_max", &digi_gr0_ch8_time_at_max, &b_digi_gr0_ch8_time_at_max);
   tree->SetBranchAddress("ADC_board_11301_0", &ADC_board_11301_0, &b_ADC_board_11301_0);
   tree->SetBranchAddress("ADC_board_11301_1", &ADC_board_11301_1, &b_ADC_board_11301_1);
   tree->SetBranchAddress("ADC_board_11301_2", &ADC_board_11301_2, &b_ADC_board_11301_2);
   tree->SetBranchAddress("ADC_board_11301_3", &ADC_board_11301_3, &b_ADC_board_11301_3);
   tree->SetBranchAddress("ADC_board_11301_4", &ADC_board_11301_4, &b_ADC_board_11301_4);
   tree->SetBranchAddress("ADC_board_11301_5", &ADC_board_11301_5, &b_ADC_board_11301_5);
   tree->SetBranchAddress("ADC_board_11301_6", &ADC_board_11301_6, &b_ADC_board_11301_6);
   tree->SetBranchAddress("ADC_board_11301_7", &ADC_board_11301_7, &b_ADC_board_11301_7);
   tree->SetBranchAddress("ADC_board_6301_0", &ADC_board_6301_0, &b_ADC_board_6301_0);
   tree->SetBranchAddress("ADC_board_6301_1", &ADC_board_6301_1, &b_ADC_board_6301_1);
   tree->SetBranchAddress("ADC_board_6301_10", &ADC_board_6301_10, &b_ADC_board_6301_10);
   tree->SetBranchAddress("ADC_board_6301_11", &ADC_board_6301_11, &b_ADC_board_6301_11);
   tree->SetBranchAddress("ADC_board_6301_12", &ADC_board_6301_12, &b_ADC_board_6301_12);
   tree->SetBranchAddress("ADC_board_6301_13", &ADC_board_6301_13, &b_ADC_board_6301_13);
   tree->SetBranchAddress("ADC_board_6301_14", &ADC_board_6301_14, &b_ADC_board_6301_14);
   tree->SetBranchAddress("ADC_board_6301_15", &ADC_board_6301_15, &b_ADC_board_6301_15);
   tree->SetBranchAddress("ADC_board_6301_16", &ADC_board_6301_16, &b_ADC_board_6301_16);
   tree->SetBranchAddress("ADC_board_6301_17", &ADC_board_6301_17, &b_ADC_board_6301_17);
   tree->SetBranchAddress("ADC_board_6301_18", &ADC_board_6301_18, &b_ADC_board_6301_18);
   tree->SetBranchAddress("ADC_board_6301_19", &ADC_board_6301_19, &b_ADC_board_6301_19);
   tree->SetBranchAddress("ADC_board_6301_2", &ADC_board_6301_2, &b_ADC_board_6301_2);
   tree->SetBranchAddress("ADC_board_6301_20", &ADC_board_6301_20, &b_ADC_board_6301_20);
   tree->SetBranchAddress("ADC_board_6301_21", &ADC_board_6301_21, &b_ADC_board_6301_21);
   tree->SetBranchAddress("ADC_board_6301_22", &ADC_board_6301_22, &b_ADC_board_6301_22);
   tree->SetBranchAddress("ADC_board_6301_23", &ADC_board_6301_23, &b_ADC_board_6301_23);
   tree->SetBranchAddress("ADC_board_6301_24", &ADC_board_6301_24, &b_ADC_board_6301_24);
   tree->SetBranchAddress("ADC_board_6301_25", &ADC_board_6301_25, &b_ADC_board_6301_25);
   tree->SetBranchAddress("ADC_board_6301_26", &ADC_board_6301_26, &b_ADC_board_6301_26);
   tree->SetBranchAddress("ADC_board_6301_27", &ADC_board_6301_27, &b_ADC_board_6301_27);
   tree->SetBranchAddress("ADC_board_6301_28", &ADC_board_6301_28, &b_ADC_board_6301_28);
   tree->SetBranchAddress("ADC_board_6301_29", &ADC_board_6301_29, &b_ADC_board_6301_29);
   tree->SetBranchAddress("ADC_board_6301_3", &ADC_board_6301_3, &b_ADC_board_6301_3);
   tree->SetBranchAddress("ADC_board_6301_30", &ADC_board_6301_30, &b_ADC_board_6301_30);
   tree->SetBranchAddress("ADC_board_6301_31", &ADC_board_6301_31, &b_ADC_board_6301_31);
   tree->SetBranchAddress("ADC_board_6301_4", &ADC_board_6301_4, &b_ADC_board_6301_4);
   tree->SetBranchAddress("ADC_board_6301_5", &ADC_board_6301_5, &b_ADC_board_6301_5);
   tree->SetBranchAddress("ADC_board_6301_6", &ADC_board_6301_6, &b_ADC_board_6301_6);
   tree->SetBranchAddress("ADC_board_6301_7", &ADC_board_6301_7, &b_ADC_board_6301_7);
   tree->SetBranchAddress("ADC_board_6301_8", &ADC_board_6301_8, &b_ADC_board_6301_8);
   tree->SetBranchAddress("ADC_board_6301_9", &ADC_board_6301_9, &b_ADC_board_6301_9);
   tree->SetBranchAddress("MULTILINE_time0", &MULTILINE_time0, &b_MULTILINE_time0);
   tree->SetBranchAddress("MULTILINE_time1", &MULTILINE_time1, &b_MULTILINE_time1);
   tree->SetBranchAddress("MULTILINE_time2", &MULTILINE_time2, &b_MULTILINE_time2);
   tree->SetBranchAddress("TDCinputTime1", &TDCinputTime1, &b_TDCinputTime1);
   tree->SetBranchAddress("TDCinputTime2", &TDCinputTime2, &b_TDCinputTime2);
   tree->SetBranchAddress("TDCinputTime3", &TDCinputTime3, &b_TDCinputTime3);
   tree->SetBranchAddress("TDCinputTime4", &TDCinputTime4, &b_TDCinputTime4);
   tree->SetBranchAddress("TDCrecoPos_X", &TDCrecoPos_X, &b_TDCrecoPos_X);
   tree->SetBranchAddress("TDCrecoX", &TDCrecoX, &b_TDCrecoX);
   tree->SetBranchAddress("TDCrecoY", &TDCrecoY, &b_TDCrecoY);
   tree->SetBranchAddress("beamPositionSmallX", &beamPositionSmallX, &b_beamPositionSmallX);
   tree->SetBranchAddress("beamPositionSmallY", &beamPositionSmallY, &b_beamPositionSmallY);
   tree->SetBranchAddress("beamPositionX", &beamPositionX, &b_beamPositionX);
   tree->SetBranchAddress("beamPositionX1", &beamPositionX1, &b_beamPositionX1);
   tree->SetBranchAddress("beamPositionX2", &beamPositionX2, &b_beamPositionX2);
   tree->SetBranchAddress("beamPositionY", &beamPositionY, &b_beamPositionY);
   tree->SetBranchAddress("beamPositionY1", &beamPositionY1, &b_beamPositionY1);
   tree->SetBranchAddress("beamPositionY2", &beamPositionY2, &b_beamPositionY2);
   tree->SetBranchAddress("beamProfileSmallX", &beamProfileSmallX, &b_beamProfileSmallX);
   tree->SetBranchAddress("beamProfileSmallY", &beamProfileSmallY, &b_beamProfileSmallY);
   tree->SetBranchAddress("beamProfileX1", &beamProfileX1, &b_beamProfileX1);
   tree->SetBranchAddress("beamProfileX2", &beamProfileX2, &b_beamProfileX2);
   tree->SetBranchAddress("beamProfileY1", &beamProfileY1, &b_beamProfileY1);
   tree->SetBranchAddress("beamProfileY2", &beamProfileY2, &b_beamProfileY2);
   tree->SetBranchAddress("deltaTime10", &deltaTime10, &b_deltaTime10);
   tree->SetBranchAddress("deltaTime20", &deltaTime20, &b_deltaTime20);
   tree->SetBranchAddress("deltaTime21", &deltaTime21, &b_deltaTime21);
   tree->SetBranchAddress("fractionTakenTrigHisto", &fractionTakenTrigHisto, &b_fractionTakenTrigHisto);
   tree->SetBranchAddress("nFibersOnSmallX", &nFibersOnSmallX, &b_nFibersOnSmallX);
   tree->SetBranchAddress("nFibersOnSmallY", &nFibersOnSmallY, &b_nFibersOnSmallY);
   tree->SetBranchAddress("nFibersOnX1", &nFibersOnX1, &b_nFibersOnX1);
   tree->SetBranchAddress("nFibersOnX2", &nFibersOnX2, &b_nFibersOnX2);
   tree->SetBranchAddress("nFibersOnY1", &nFibersOnY1, &b_nFibersOnY1);
   tree->SetBranchAddress("nFibersOnY2", &nFibersOnY2, &b_nFibersOnY2);
   tree->SetBranchAddress("nTotalEvts", &nTotalEvts, &b_nTotalEvts);




   std::string outdir = "analysisTrees_" + tag;
   system( Form("mkdir -p %s", outdir.c_str()) );

   std::string outfileName = outdir + "/Reco_" + runName + ".root";
   TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );
   TTree* outTree = new TTree("recoTree", "recoTree");



   int run;
   outTree->Branch( "run", &run, "run/i" );
   int event;
   outTree->Branch( "event", &event, "event/i" );

   float s1;
   outTree->Branch( "s1", &s1, "s1/F" );
   float s3;
   outTree->Branch( "s3", &s3, "s3/F" );
   float s4;
   outTree->Branch( "s4", &s4, "s4/F" );
   float s6;
   outTree->Branch( "s6", &s6, "s6/F" );

   int cef3_chan(CEF3_CHANNELS);
   outTree->Branch( "cef3_chan", &cef3_chan, "cef3_chan/I" );
   float cef3[CEF3_CHANNELS];
   outTree->Branch( "cef3", cef3, "cef3[cef3_chan]/F" );

   int bgo_chan(BGO_CHANNELS);
   outTree->Branch( "bgo_chan", &bgo_chan, "bgo_chan/I" );
   float bgo[BGO_CHANNELS];
   outTree->Branch( "bgo", bgo, "bgo[bgo_chan]/F" );

   float xBeam;
   outTree->Branch( "xBeam", &xBeam, "xBeam/F");
   float yBeam;
   outTree->Branch( "yBeam", &yBeam, "yBeam/F");


   int hodoSmallX_chan(HODOSMALLX_CHANNELS);
   outTree->Branch( "hodoSmallX_chan", &hodoSmallX_chan, "hodoSmallX_chan/I");
   float hodoSmallX[HODOSMALLX_CHANNELS];
   outTree->Branch( "hodoSmallX", hodoSmallX, "hodoSmallX[hodoSmallX_chan]/F" );
   int hodoSmallY_chan(HODOSMALLY_CHANNELS);
   outTree->Branch( "hodoSmallY_chan", &hodoSmallY_chan, "hodoSmallY_chan/I");
   float hodoSmallY[HODOSMALLY_CHANNELS];
   outTree->Branch( "hodoSmallY", hodoSmallY, "hodoSmallY[hodoSmallY_chan]/F" );

   int hodoX1_chan(HODOX1_CHANNELS);
   outTree->Branch( "hodoX1_chan", &hodoX1_chan, "hodoX1_chan/I");
   float hodoX1[HODOX1_CHANNELS];
   outTree->Branch( "hodoX1", hodoX1, "hodoX1[hodoX1_chan]/F" );
   int hodoY1_chan(HODOY1_CHANNELS);
   outTree->Branch( "hodoY1_chan", &hodoY1_chan, "hodoY1_chan/I");
   float hodoY1[HODOY1_CHANNELS];
   outTree->Branch( "hodoY1", hodoY1, "hodoY1[hodoY1_chan]/F" );

   int hodoX2_chan(HODOX2_CHANNELS);
   outTree->Branch( "hodoX2_chan", &hodoX2_chan, "hodoX2_chan/I");
   float hodoX2[HODOX2_CHANNELS];
   outTree->Branch( "hodoX2", hodoX2, "hodoX2[hodoX2_chan]/F" );
   int hodoY2_chan(HODOY2_CHANNELS);
   outTree->Branch( "hodoY2_chan", &hodoY2_chan, "hodoY2_chan/I");
   float hodoY2[HODOY2_CHANNELS];
   outTree->Branch( "hodoY2", hodoY2, "hodoY2[hodoY2_chan]/F" );



   int nentries = tree->GetEntries();
   RunHelper::getBeamPosition( runName, xBeam, yBeam );

   for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

     tree->GetEntry( iEntry );

     if( iEntry %  5000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

     outTree->Fill();

   } // for entries


   
   outfile->cd();
   outTree->Write();
   outfile->Close();

   std::cout << "-> Analysis Tree saved in: " << outfile->GetName() << std::endl;
   return 0;

}
  
