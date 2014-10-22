#include <iostream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>



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
   tree = (TChain*)file->Get("outputTree");



   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<float>   *ADC_board_11301_0;
   vector<float>   *ADC_board_11301_1;
   vector<float>   *ADC_board_11301_2;
   vector<float>   *ADC_board_11301_3;
   vector<float>   *ADC_board_11301_4;
   vector<float>   *ADC_board_11301_5;
   vector<float>   *ADC_board_11301_6;
   vector<float>   *ADC_board_11301_7;
   vector<float>   *ADC_board_6301_0;
   vector<float>   *ADC_board_6301_1;
   vector<float>   *ADC_board_6301_10;
   vector<float>   *ADC_board_6301_11;
   vector<float>   *ADC_board_6301_12;
   vector<float>   *ADC_board_6301_13;
   vector<float>   *ADC_board_6301_14;
   vector<float>   *ADC_board_6301_15;
   vector<float>   *ADC_board_6301_16;
   vector<float>   *ADC_board_6301_17;
   vector<float>   *ADC_board_6301_18;
   vector<float>   *ADC_board_6301_19;
   vector<float>   *ADC_board_6301_2;
   vector<float>   *ADC_board_6301_20;
   vector<float>   *ADC_board_6301_21;
   vector<float>   *ADC_board_6301_22;
   vector<float>   *ADC_board_6301_23;
   vector<float>   *ADC_board_6301_24;
   vector<float>   *ADC_board_6301_25;
   vector<float>   *ADC_board_6301_26;
   vector<float>   *ADC_board_6301_27;
   vector<float>   *ADC_board_6301_28;
   vector<float>   *ADC_board_6301_29;
   vector<float>   *ADC_board_6301_3;
   vector<float>   *ADC_board_6301_30;
   vector<float>   *ADC_board_6301_31;
   vector<float>   *ADC_board_6301_4;
   vector<float>   *ADC_board_6301_5;
   vector<float>   *ADC_board_6301_6;
   vector<float>   *ADC_board_6301_7;
   vector<float>   *ADC_board_6301_8;
   vector<float>   *ADC_board_6301_9;
   vector<float>   *MULTILINE_time0;
   vector<float>   *MULTILINE_time1;
   vector<float>   *MULTILINE_time2;
   vector<float>   *TDCinputTime1;
   vector<float>   *TDCinputTime2;
   vector<float>   *TDCinputTime3;
   vector<float>   *TDCinputTime4;
   vector<float>   *TDCrecoPos_X;
   vector<float>   *TDCrecoPos_X;
   vector<float>   *TDCrecoX;
   vector<float>   *TDCrecoY;
   vector<float>   *beamPositionSmallX;
   vector<float>   *beamPositionSmallY;
   vector<float>   *beamPositionX;
   vector<float>   *beamPositionX1;
   vector<float>   *beamPositionX2;
   vector<float>   *beamPositionY;
   vector<float>   *beamPositionY1;
   vector<float>   *beamPositionY2;
   vector<float>   *beamProfileSmallX;
   vector<float>   *beamProfileSmallY;
   vector<float>   *beamProfileX1;
   vector<float>   *beamProfileX2;
   vector<float>   *beamProfileY1;
   vector<float>   *beamProfileY2;
   vector<float>   *deltaTime10;
   vector<float>   *deltaTime20;
   vector<float>   *deltaTime21;
   vector<float>   *fractionTakenTrigHisto;
   vector<float>   *nFibersOnSmallX;
   vector<float>   *nFibersOnSmallY;
   vector<float>   *nFibersOnX1;
   vector<float>   *nFibersOnX2;
   vector<float>   *nFibersOnY1;
   vector<float>   *nFibersOnY2;
   vector<float>   *nTotalEvts;

   // List of branches
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
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ADC_board_11301_0", &ADC_board_11301_0, &b_ADC_board_11301_0);
   fChain->SetBranchAddress("ADC_board_11301_1", &ADC_board_11301_1, &b_ADC_board_11301_1);
   fChain->SetBranchAddress("ADC_board_11301_2", &ADC_board_11301_2, &b_ADC_board_11301_2);
   fChain->SetBranchAddress("ADC_board_11301_3", &ADC_board_11301_3, &b_ADC_board_11301_3);
   fChain->SetBranchAddress("ADC_board_11301_4", &ADC_board_11301_4, &b_ADC_board_11301_4);
   fChain->SetBranchAddress("ADC_board_11301_5", &ADC_board_11301_5, &b_ADC_board_11301_5);
   fChain->SetBranchAddress("ADC_board_11301_6", &ADC_board_11301_6, &b_ADC_board_11301_6);
   fChain->SetBranchAddress("ADC_board_11301_7", &ADC_board_11301_7, &b_ADC_board_11301_7);
   fChain->SetBranchAddress("ADC_board_6301_0", &ADC_board_6301_0, &b_ADC_board_6301_0);
   fChain->SetBranchAddress("ADC_board_6301_1", &ADC_board_6301_1, &b_ADC_board_6301_1);
   fChain->SetBranchAddress("ADC_board_6301_10", &ADC_board_6301_10, &b_ADC_board_6301_10);
   fChain->SetBranchAddress("ADC_board_6301_11", &ADC_board_6301_11, &b_ADC_board_6301_11);
   fChain->SetBranchAddress("ADC_board_6301_12", &ADC_board_6301_12, &b_ADC_board_6301_12);
   fChain->SetBranchAddress("ADC_board_6301_13", &ADC_board_6301_13, &b_ADC_board_6301_13);
   fChain->SetBranchAddress("ADC_board_6301_14", &ADC_board_6301_14, &b_ADC_board_6301_14);
   fChain->SetBranchAddress("ADC_board_6301_15", &ADC_board_6301_15, &b_ADC_board_6301_15);
   fChain->SetBranchAddress("ADC_board_6301_16", &ADC_board_6301_16, &b_ADC_board_6301_16);
   fChain->SetBranchAddress("ADC_board_6301_17", &ADC_board_6301_17, &b_ADC_board_6301_17);
   fChain->SetBranchAddress("ADC_board_6301_18", &ADC_board_6301_18, &b_ADC_board_6301_18);
   fChain->SetBranchAddress("ADC_board_6301_19", &ADC_board_6301_19, &b_ADC_board_6301_19);
   fChain->SetBranchAddress("ADC_board_6301_2", &ADC_board_6301_2, &b_ADC_board_6301_2);
   fChain->SetBranchAddress("ADC_board_6301_20", &ADC_board_6301_20, &b_ADC_board_6301_20);
   fChain->SetBranchAddress("ADC_board_6301_21", &ADC_board_6301_21, &b_ADC_board_6301_21);
   fChain->SetBranchAddress("ADC_board_6301_22", &ADC_board_6301_22, &b_ADC_board_6301_22);
   fChain->SetBranchAddress("ADC_board_6301_23", &ADC_board_6301_23, &b_ADC_board_6301_23);
   fChain->SetBranchAddress("ADC_board_6301_24", &ADC_board_6301_24, &b_ADC_board_6301_24);
   fChain->SetBranchAddress("ADC_board_6301_25", &ADC_board_6301_25, &b_ADC_board_6301_25);
   fChain->SetBranchAddress("ADC_board_6301_26", &ADC_board_6301_26, &b_ADC_board_6301_26);
   fChain->SetBranchAddress("ADC_board_6301_27", &ADC_board_6301_27, &b_ADC_board_6301_27);
   fChain->SetBranchAddress("ADC_board_6301_28", &ADC_board_6301_28, &b_ADC_board_6301_28);
   fChain->SetBranchAddress("ADC_board_6301_29", &ADC_board_6301_29, &b_ADC_board_6301_29);
   fChain->SetBranchAddress("ADC_board_6301_3", &ADC_board_6301_3, &b_ADC_board_6301_3);
   fChain->SetBranchAddress("ADC_board_6301_30", &ADC_board_6301_30, &b_ADC_board_6301_30);
   fChain->SetBranchAddress("ADC_board_6301_31", &ADC_board_6301_31, &b_ADC_board_6301_31);
   fChain->SetBranchAddress("ADC_board_6301_4", &ADC_board_6301_4, &b_ADC_board_6301_4);
   fChain->SetBranchAddress("ADC_board_6301_5", &ADC_board_6301_5, &b_ADC_board_6301_5);
   fChain->SetBranchAddress("ADC_board_6301_6", &ADC_board_6301_6, &b_ADC_board_6301_6);
   fChain->SetBranchAddress("ADC_board_6301_7", &ADC_board_6301_7, &b_ADC_board_6301_7);
   fChain->SetBranchAddress("ADC_board_6301_8", &ADC_board_6301_8, &b_ADC_board_6301_8);
   fChain->SetBranchAddress("ADC_board_6301_9", &ADC_board_6301_9, &b_ADC_board_6301_9);
   fChain->SetBranchAddress("MULTILINE_time0", &MULTILINE_time0, &b_MULTILINE_time0);
   fChain->SetBranchAddress("MULTILINE_time1", &MULTILINE_time1, &b_MULTILINE_time1);
   fChain->SetBranchAddress("MULTILINE_time2", &MULTILINE_time2, &b_MULTILINE_time2);
   fChain->SetBranchAddress("TDCinputTime1", &TDCinputTime1, &b_TDCinputTime1);
   fChain->SetBranchAddress("TDCinputTime2", &TDCinputTime2, &b_TDCinputTime2);
   fChain->SetBranchAddress("TDCinputTime3", &TDCinputTime3, &b_TDCinputTime3);
   fChain->SetBranchAddress("TDCinputTime4", &TDCinputTime4, &b_TDCinputTime4);
   fChain->SetBranchAddress("TDCrecoPos_X", &TDCrecoPos_X, &b_TDCrecoPos_X);
//    fChain->SetBranchAddress("TDCrecoPos_X", &TDCrecoPos_X, &b_TDCrecoPos_X);
   fChain->SetBranchAddress("TDCrecoX", &TDCrecoX, &b_TDCrecoX);
   fChain->SetBranchAddress("TDCrecoY", &TDCrecoY, &b_TDCrecoY);
   fChain->SetBranchAddress("beamPositionSmallX", &beamPositionSmallX, &b_beamPositionSmallX);
   fChain->SetBranchAddress("beamPositionSmallY", &beamPositionSmallY, &b_beamPositionSmallY);
   fChain->SetBranchAddress("beamPositionX", &beamPositionX, &b_beamPositionX);
   fChain->SetBranchAddress("beamPositionX1", &beamPositionX1, &b_beamPositionX1);
   fChain->SetBranchAddress("beamPositionX2", &beamPositionX2, &b_beamPositionX2);
   fChain->SetBranchAddress("beamPositionY", &beamPositionY, &b_beamPositionY);
   fChain->SetBranchAddress("beamPositionY1", &beamPositionY1, &b_beamPositionY1);
   fChain->SetBranchAddress("beamPositionY2", &beamPositionY2, &b_beamPositionY2);
   fChain->SetBranchAddress("beamProfileSmallX", &beamProfileSmallX, &b_beamProfileSmallX);
   fChain->SetBranchAddress("beamProfileSmallY", &beamProfileSmallY, &b_beamProfileSmallY);
   fChain->SetBranchAddress("beamProfileX1", &beamProfileX1, &b_beamProfileX1);
   fChain->SetBranchAddress("beamProfileX2", &beamProfileX2, &b_beamProfileX2);
   fChain->SetBranchAddress("beamProfileY1", &beamProfileY1, &b_beamProfileY1);
   fChain->SetBranchAddress("beamProfileY2", &beamProfileY2, &b_beamProfileY2);
   fChain->SetBranchAddress("deltaTime10", &deltaTime10, &b_deltaTime10);
   fChain->SetBranchAddress("deltaTime20", &deltaTime20, &b_deltaTime20);
   fChain->SetBranchAddress("deltaTime21", &deltaTime21, &b_deltaTime21);
   fChain->SetBranchAddress("fractionTakenTrigHisto", &fractionTakenTrigHisto, &b_fractionTakenTrigHisto);
   fChain->SetBranchAddress("nFibersOnSmallX", &nFibersOnSmallX, &b_nFibersOnSmallX);
   fChain->SetBranchAddress("nFibersOnSmallY", &nFibersOnSmallY, &b_nFibersOnSmallY);
   fChain->SetBranchAddress("nFibersOnX1", &nFibersOnX1, &b_nFibersOnX1);
   fChain->SetBranchAddress("nFibersOnX2", &nFibersOnX2, &b_nFibersOnX2);
   fChain->SetBranchAddress("nFibersOnY1", &nFibersOnY1, &b_nFibersOnY1);
   fChain->SetBranchAddress("nFibersOnY2", &nFibersOnY2, &b_nFibersOnY2);
   fChain->SetBranchAddress("nTotalEvts", &nTotalEvts, &b_nTotalEvts);




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
   float s2;
   outTree->Branch( "s2", &s2, "s2/F" );
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
   RunHelper::getBeamPosition( runName, xBeam_, yBeam_ );

   for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {
