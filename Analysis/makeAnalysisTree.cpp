#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TChain.h"


#include "channelInfo.h"
#include "interface/HodoCluster.h"
#include "interface/TagHelper.h"
#include "interface/EnergyCalibration.h"
#include "interface/AlignmentOfficer.h"

#include "CommonTools/interface/RunHelper.h"


void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos );
void assignValuesBool( std::vector<bool> &target, std::vector<bool> source, unsigned int startPos );
void doHodoReconstructionBool( std::vector<bool> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut );
std::vector<HodoCluster*> getHodoClustersBool( std::vector<bool> hodo, float fibreWidth, int nClusterMax, float Cut );

void doHodoReconstruction( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut );
std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax, float Cut );
void copyArray( int n, float *source, float *target );






int main( int argc, char* argv[] ) {


   if( argc<2 ) {
     std::cout << "ERROR. You need to specify the name of the run you want to process." << std::endl;  
     exit(1);
   }

   std::string runName = "";
   std::string tag = "V01";

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

   std::string fileName = "data/Corr04_12/run_" + runName + ".root";
   TFile* file = TFile::Open(fileName.c_str());
   if( file==0 ) {
     std::cout << "ERROR! Din't find file " << fileName << std::endl;
     std::cout << "Exiting." << std::endl;
     exit(11);
   }
   TTree* tree = (TTree*)file->Get("outputTree");

   //   TChain * chain = new TChain("tree","");
   //   tree = chain;
   Int_t           fCurrent;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
 

   // Declaration of leaf types
   UInt_t          runNumber;
   UInt_t          spillNumber;
   UInt_t          evtNumber;
   std::vector<float>   *BGOvalues;
   std::vector<float>   *SCINTvalues;
   std::vector<float>   *TDCreco;
   std::vector<float>   *digi_charge_integrated;
   std::vector<float>   *digi_max_amplitude;
   std::vector<float>   *digi_pedestal;
   std::vector<float>   *digi_pedestal_rms;
   std::vector<float>   *digi_time_at_frac30;
   std::vector<float>   *digi_time_at_frac50;
   std::vector<float>   *digi_time_at_max;
   std::vector<bool>    *HODOX1;
   std::vector<bool>    *HODOX2;
   std::vector<bool>    *HODOY1;
   std::vector<bool>    *HODOY2;
   Float_t         TableX;
   Float_t         TableY;
   Float_t         CeF3HV;
   Float_t         BGOHV;
   Float_t         BeamEnergy;
   Float_t         BeamTilt;
   Int_t           IsPhysics;
   std::vector<int>     *nTdcHits;
   std::vector<float>   *digi_charge_integrated_sub;
   std::vector<float>   *digi_max_amplitude_sub;
   std::vector<float>   *digi_pedestal_sub;
   std::vector<float>   *digi_pedestal_rms_sub;
   std::vector<float>   *digi_charge_integrated_corr1;
   std::vector<float>   *digi_max_amplitude_corr1;
   std::vector<float>   *digi_charge_integrated_corr2;
   std::vector<float>   *digi_max_amplitude_corr2;

   std::vector<float>   *digi_max_amplitude_bare;
   std::vector<float>   *digi_charge_integrated_bare;
   std::vector<float>   *digi_charge_integrated_frac10;
   std::vector<float>   *digi_charge_integrated_frac30;
   std::vector<float>   *digi_charge_integrated_frac50;


   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_spillNumber;   //!
   TBranch        *b_evtNumber;   //!
   TBranch        *b_BGOvalues;   //!
   TBranch        *b_SCINTvalues;   //!
   TBranch        *b_TDCreco;   //!
   TBranch        *b_digi_charge_integrated;   //!
   TBranch        *b_digi_max_amplitude;   //!
   TBranch        *b_digi_pedestal;   //!
   TBranch        *b_digi_pedestal_rms;   //!
   TBranch        *b_digi_time_at_frac30;   //!
   TBranch        *b_digi_time_at_frac50;   //!
   TBranch        *b_digi_time_at_max;   //!
   TBranch        *b_HODOX1;   //!
   TBranch        *b_HODOX2;   //!
   TBranch        *b_HODOY1;   //!
   TBranch        *b_HODOY2;   //!
   TBranch        *b_TableX;   //!
   TBranch        *b_TableY;   //!
   TBranch        *b_CeF3HV;   //!
   TBranch        *b_BGOHV;   //!
   TBranch        *b_BeamEnergy;   //!
   TBranch        *b_BeamTilt;   //!
   TBranch        *b_IsPhysics;   //!
   TBranch        *b_nTdcHits;   //!
   TBranch        *b_digi_charge_integrated_sub;   //!
   TBranch        *b_digi_max_amplitude_sub;   //!
   TBranch        *b_digi_pedestal_sub;   //!
   TBranch        *b_digi_pedestal_rms_sub;   //!
   TBranch        *b_digi_charge_integrated_corr1;   //!
   TBranch        *b_digi_max_amplitude_corr1;   //!
   TBranch        *b_digi_charge_integrated_corr2;   //!
   TBranch        *b_digi_max_amplitude_corr2;   //!

   TBranch *b_digi_max_amplitude_bare;
   TBranch *b_digi_charge_integrated_bare;
   TBranch *b_digi_charge_integrated_frac10;
   TBranch *b_digi_charge_integrated_frac30;
   TBranch *b_digi_charge_integrated_frac50;

   // Set object pointer
   BGOvalues = 0;
   SCINTvalues = 0;
   TDCreco = 0;
   digi_charge_integrated = 0;
   digi_max_amplitude = 0;
   digi_pedestal = 0;
   digi_pedestal_rms = 0;
   digi_time_at_frac30 = 0;
   digi_time_at_frac50 = 0;
   digi_time_at_max = 0;
   HODOX1 = 0;
   HODOX2 = 0;
   HODOY1 = 0;
   HODOY2 = 0;
   nTdcHits = 0;
   digi_charge_integrated_sub = 0;
   digi_max_amplitude_sub = 0;
   digi_pedestal_sub = 0;
   digi_pedestal_rms_sub = 0;
   digi_charge_integrated_corr1 = 0;
   digi_max_amplitude_corr1 = 0;
   digi_charge_integrated_corr2 = 0;
   digi_max_amplitude_corr2 = 0;

   digi_max_amplitude_bare = 0;
   digi_charge_integrated_bare = 0;
   digi_charge_integrated_frac10 = 0;
   digi_charge_integrated_frac30 = 0;
   digi_charge_integrated_frac50 = 0;
   

   // Set branch addresses and branch pointers
 fChain = tree;

   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("spillNumber", &spillNumber, &b_spillNumber);
   fChain->SetBranchAddress("evtNumber", &evtNumber, &b_evtNumber);
   fChain->SetBranchAddress("BGOvalues", &BGOvalues, &b_BGOvalues);
   fChain->SetBranchAddress("SCINTvalues", &SCINTvalues, &b_SCINTvalues);
   fChain->SetBranchAddress("TDCreco", &TDCreco, &b_TDCreco);
   fChain->SetBranchAddress("digi_charge_integrated", &digi_charge_integrated, &b_digi_charge_integrated);
   fChain->SetBranchAddress("digi_max_amplitude", &digi_max_amplitude, &b_digi_max_amplitude);
   fChain->SetBranchAddress("digi_pedestal", &digi_pedestal, &b_digi_pedestal);
   fChain->SetBranchAddress("digi_pedestal_rms", &digi_pedestal_rms, &b_digi_pedestal_rms);
   fChain->SetBranchAddress("digi_time_at_frac30", &digi_time_at_frac30, &b_digi_time_at_frac30);
   fChain->SetBranchAddress("digi_time_at_frac50", &digi_time_at_frac50, &b_digi_time_at_frac50);
   fChain->SetBranchAddress("digi_time_at_max", &digi_time_at_max, &b_digi_time_at_max);
   fChain->SetBranchAddress("HODOX1", &HODOX1, &b_HODOX1);
   fChain->SetBranchAddress("HODOX2", &HODOX2, &b_HODOX2);
   fChain->SetBranchAddress("HODOY1", &HODOY1, &b_HODOY1);
   fChain->SetBranchAddress("HODOY2", &HODOY2, &b_HODOY2);
   fChain->SetBranchAddress("TableX", &TableX, &b_TableX);
   fChain->SetBranchAddress("TableY", &TableY, &b_TableY);
   fChain->SetBranchAddress("CeF3HV", &CeF3HV, &b_CeF3HV);
   fChain->SetBranchAddress("BGOHV", &BGOHV, &b_BGOHV);
   fChain->SetBranchAddress("BeamEnergy", &BeamEnergy, &b_BeamEnergy);
   fChain->SetBranchAddress("BeamTilt", &BeamTilt, &b_BeamTilt);
   fChain->SetBranchAddress("IsPhysics", &IsPhysics, &b_IsPhysics);
   fChain->SetBranchAddress("nTdcHits", &nTdcHits, &b_nTdcHits);
   fChain->SetBranchAddress("digi_charge_integrated_sub", &digi_charge_integrated_sub, &b_digi_charge_integrated_sub);
   fChain->SetBranchAddress("digi_max_amplitude_sub", &digi_max_amplitude_sub, &b_digi_max_amplitude_sub);
   fChain->SetBranchAddress("digi_pedestal_sub", &digi_pedestal_sub, &b_digi_pedestal_sub);
   fChain->SetBranchAddress("digi_pedestal_rms_sub", &digi_pedestal_rms_sub, &b_digi_pedestal_rms_sub);
   fChain->SetBranchAddress("digi_charge_integrated_corr1", &digi_charge_integrated_corr1, &b_digi_charge_integrated_corr1);
   fChain->SetBranchAddress("digi_max_amplitude_corr1", &digi_max_amplitude_corr1, &b_digi_max_amplitude_corr1);
   fChain->SetBranchAddress("digi_charge_integrated_corr2", &digi_charge_integrated_corr2, &b_digi_charge_integrated_corr2);
   fChain->SetBranchAddress("digi_max_amplitude_corr2", &digi_max_amplitude_corr2, &b_digi_max_amplitude_corr2);

   fChain->SetBranchAddress("digi_max_amplitude_bare", &digi_max_amplitude_bare, &b_digi_max_amplitude_bare);
   fChain->SetBranchAddress("digi_charge_integrated_bare", &digi_charge_integrated_bare, &b_digi_charge_integrated_bare);
   fChain->SetBranchAddress("digi_charge_integrated_frac10", &digi_charge_integrated_frac10, &b_digi_charge_integrated_frac10);
   fChain->SetBranchAddress("digi_charge_integrated_frac30", &digi_charge_integrated_frac30, &b_digi_charge_integrated_frac30);
   fChain->SetBranchAddress("digi_charge_integrated_frac50", &digi_charge_integrated_frac50, &b_digi_charge_integrated_frac50);



   /*
   tree->GetEntry( 42 );
   
   std::string theBeamEnergy = Form("%.0f",BeamEnergy);
   if( runNumber > 272 && runNumber < 298){
     theBeamEnergy = "273"; //For the long position scan prior to tdc adjustment
   }
   std::cout << "The used constant file has label = "<< theBeamEnergy  << std::endl;
   

   //set the tag for calibration
   TagHelper tagHelper(tag,theBeamEnergy);
   EnergyCalibration cef3Calib(tagHelper.getCeF3FileName());
   EnergyCalibration bgoCalib(tagHelper.getBGOFileName());
   AlignmentOfficer alignOfficer(tagHelper.getAlignmentFileName());
   */






  ///AND FINALLY THE RECO TREE/////////  
   std::string outdir = "analysisTrees_" + tag;
   system( Form("mkdir -p %s", outdir.c_str()) );

   std::string outfileName = outdir + "/Reco_" + runName + ".root";
   TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );
   TTree* outTree = new TTree("recoTree", "recoTree");


   //and now naming the stuff for the recoTree//
   outTree->Branch( "run", &runNumber, "run/i" );
   outTree->Branch( "spill", &spillNumber, "spill/i" );
   outTree->Branch( "event", &evtNumber, "event/i" );

   float s1;
   outTree->Branch( "s1", &s1, "s1/F" );
   float s3;
   outTree->Branch( "s3", &s3, "s3/F" );
   float s4;
   outTree->Branch( "s4", &s4, "s4/F" );
   float s6;
   outTree->Branch( "s6", &s6, "s6/F" );

   //Original
   std::vector<float> cef3( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3", &cef3 ); //for backwards compatibility
   std::vector<float> cef3_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_corr", &cef3_corr );


   std::vector<float> cef3_maxAmpl( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl", &cef3_maxAmpl );
   std::vector<float> cef3_chaInt( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt", &cef3_chaInt );

   std::vector<float> bgo( BGO_CHANNELS, -1. );
   outTree->Branch( "bgo", &bgo );

   //Original Intercalibrated
   std::vector<float> cef3_maxAmpl_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_corr", &cef3_maxAmpl_corr );
   std::vector<float> cef3_chaInt_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_corr", &cef3_chaInt_corr );

   std::vector<float> bgo_corr( BGO_CHANNELS, -1. );
   outTree->Branch( "bgo_corr", &bgo_corr );

   /*
   //Pedestal Subtracted
   std::vector<float> cef3_sub( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_sub", &cef3_sub );
   //Pedestal Subtracted Intercalibrated
   std::vector<float> cef3_sub_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_sub_corr", &cef3_sub_corr );
   */

   //With Corr1
   std::vector<float> cef3_maxAmpl_corr1( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_corr1", &cef3_maxAmpl_corr1 );
   std::vector<float> cef3_chaInt_corr1( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_corr1", &cef3_chaInt_corr1 );

  //With Corr1 Intercalibrated
   std::vector<float> cef3_maxAmpl_corr1_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_corr1_corr", &cef3_maxAmpl_corr1_corr );
   std::vector<float> cef3_chaInt_corr1_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_corr1_corr", &cef3_chaInt_corr1_corr );


   std::vector<float> cef3_maxAmpl_bare( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_bare", &cef3_maxAmpl_bare );
   std::vector<float> cef3_chaInt_bare( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_bare", &cef3_chaInt_bare );
   std::vector<float> cef3_maxAmpl_bare_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_maxAmpl_bare_corr", &cef3_maxAmpl_bare_corr );
   std::vector<float> cef3_chaInt_bare_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_bare_corr", &cef3_chaInt_bare_corr );


   std::vector<float> cef3_chaInt_frac10( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_frac10", &cef3_chaInt_frac10 );
   std::vector<float> cef3_chaInt_frac30( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_frac30", &cef3_chaInt_frac30 );
   std::vector<float> cef3_chaInt_frac50( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_frac50", &cef3_chaInt_frac50 );

   std::vector<float> cef3_chaInt_frac10_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_frac10_corr", &cef3_chaInt_frac10_corr );
   std::vector<float> cef3_chaInt_frac30_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_frac30_corr", &cef3_chaInt_frac30_corr );
   std::vector<float> cef3_chaInt_frac50_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_chaInt_frac50_corr", &cef3_chaInt_frac50_corr );



   float xTable;
   outTree->Branch( "xTable", &xTable, "xTable/F");
   float yTable;
   outTree->Branch( "yTable", &yTable, "yTable/F");

   float beamEnergy;
   outTree->Branch( "beamEnergy", &beamEnergy, "xbeamEnergy/F");

   float angle; //aka BeamTilt..
   outTree->Branch("angle", &angle, "angle/F");

   float HVCeF3;
   outTree->Branch( "HVCeF3", &HVCeF3, "HVCeF3/F");
   float HVBGO;
   outTree->Branch( "HVBGO", &HVBGO, "HVBGO/F");

   float xBeam;
   outTree->Branch( "xBeam", &xBeam, "xBeam/F");
   float yBeam;
   outTree->Branch( "yBeam", &yBeam, "yBeam/F");

   int nClusters_hodoX1;
   outTree->Branch( "nClusters_hodoX1", &nClusters_hodoX1, "nClusters_hodoX1/I" );
   int nFibres_hodoX1[HODOX1_CHANNELS];
   outTree->Branch( "nFibres_hodoX1", nFibres_hodoX1, "nFibres_hodoX1[nClusters_hodoX1]/I" );
   float pos_hodoX1[HODOX1_CHANNELS];
   outTree->Branch( "pos_hodoX1", pos_hodoX1, "pos_hodoX1[nClusters_hodoX1]/F" );
   float pos_corr_hodoX1[HODOX1_CHANNELS];
   outTree->Branch( "pos_corr_hodoX1", pos_corr_hodoX1, "pos_corr_hodoX1[nClusters_hodoX1]/F" );

   int nClusters_hodoY1;
   outTree->Branch( "nClusters_hodoY1", &nClusters_hodoY1, "nClusters_hodoY1/I" );
   int nFibres_hodoY1[HODOY1_CHANNELS];
   outTree->Branch( "nFibres_hodoY1", nFibres_hodoY1, "nFibres_hodoY1[nClusters_hodoY1]/I" );
   float pos_hodoY1[HODOY1_CHANNELS];
   outTree->Branch( "pos_hodoY1", pos_hodoY1, "pos_hodoY1[nClusters_hodoY1]/F" );
   float pos_corr_hodoY1[HODOY1_CHANNELS];
   outTree->Branch( "pos_corr_hodoY1", pos_corr_hodoY1, "pos_corr_hodoY1[nClusters_hodoY1]/F" );

   int nClusters_hodoX2;
   outTree->Branch( "nClusters_hodoX2", &nClusters_hodoX2, "nClusters_hodoX2/I" );
   int nFibres_hodoX2[HODOX2_CHANNELS];
   outTree->Branch( "nFibres_hodoX2", nFibres_hodoX2, "nFibres_hodoX2[nClusters_hodoX2]/I" );
   float pos_hodoX2[HODOX2_CHANNELS];
   outTree->Branch( "pos_hodoX2", pos_hodoX2, "pos_hodoX2[nClusters_hodoX2]/F" );
   float pos_corr_hodoX2[HODOX2_CHANNELS];
   outTree->Branch( "pos_corr_hodoX2", pos_corr_hodoX2, "pos_corr_hodoX2[nClusters_hodoX2]/F" );

   int nClusters_hodoY2;
   outTree->Branch( "nClusters_hodoY2", &nClusters_hodoY2, "nClusters_hodoY2/I" );
   int nFibres_hodoY2[HODOY2_CHANNELS];
   outTree->Branch( "nFibres_hodoY2", nFibres_hodoY2, "nFibres_hodoY2[nClusters_hodoY2]/I" );
   float pos_hodoY2[HODOY2_CHANNELS];
   outTree->Branch( "pos_hodoY2", pos_hodoY2, "pos_hodoY2[nClusters_hodoY2]/F" );
   float pos_corr_hodoY2[HODOY2_CHANNELS];
   outTree->Branch( "pos_corr_hodoY2", pos_corr_hodoY2, "pos_corr_hodoY2[nClusters_hodoY2]/F" );

   int nClusters_hodoSmallX;
   outTree->Branch( "nClusters_hodoSmallX", &nClusters_hodoSmallX, "nClusters_hodoSmallX/I" );
   int nFibres_hodoSmallX[HODOSMALLX_CHANNELS];
   outTree->Branch( "nFibres_hodoSmallX", nFibres_hodoSmallX, "nFibres_hodoSmallX[nClusters_hodoSmallX]/I" );
   float pos_hodoSmallX[HODOSMALLX_CHANNELS];
   outTree->Branch( "pos_hodoSmallX", pos_hodoSmallX, "pos_hodoSmallX[nClusters_hodoSmallX]/F" );
   float pos_corr_hodoSmallX[HODOSMALLX_CHANNELS];
   outTree->Branch( "pos_corr_hodoSmallX", pos_corr_hodoSmallX, "pos_corr_hodoSmallX[nClusters_hodoSmallX]/F" );

   int nClusters_hodoSmallY;
   outTree->Branch( "nClusters_hodoSmallY", &nClusters_hodoSmallY, "nClusters_hodoSmallY/I" );
   int nFibres_hodoSmallY[HODOSMALLY_CHANNELS];
   outTree->Branch( "nFibres_hodoSmallY", nFibres_hodoSmallY, "nFibres_hodoSmallY[nClusters_hodoSmallY]/I" );
   float pos_hodoSmallY[HODOSMALLY_CHANNELS];
   outTree->Branch( "pos_hodoSmallY", pos_hodoSmallY, "pos_hodoSmallY[nClusters_hodoSmallY]/F" );
   float pos_corr_hodoSmallY[HODOSMALLY_CHANNELS];
   outTree->Branch( "pos_corr_hodoSmallY", pos_corr_hodoSmallY, "pos_corr_hodoSmallY[nClusters_hodoSmallY]/F" );



   
   std::vector<int> nTDCHits( 4, -1. );
   outTree->Branch( "nTDCHits", &nTDCHits );




   float pos_2FibClust_hodoX1;
   outTree->Branch( "pos_2FibClust_hodoX1", &pos_2FibClust_hodoX1, "pos_2FibClust_hodoX1/F" );
   float pos_2FibClust_corr_hodoX1;
   outTree->Branch( "pos_2FibClust_corr_hodoX1", &pos_2FibClust_corr_hodoX1, "pos_2FibClust_corr_hodoX1/F" );

   float pos_2FibClust_hodoY1;
   outTree->Branch( "pos_2FibClust_hodoY1", &pos_2FibClust_hodoY1, "pos_2FibClust_hodoY1/F" );
   float pos_2FibClust_corr_hodoY1;
   outTree->Branch( "pos_2FibClust_corr_hodoY1", &pos_2FibClust_corr_hodoY1, "pos_2FibClust_corr_hodoY1/F" );

   float pos_2FibClust_hodoX2;
   outTree->Branch( "pos_2FibClust_hodoX2", &pos_2FibClust_hodoX2, "pos_2FibClust_hodoX2/F" );
   float pos_2FibClust_corr_hodoX2;
   outTree->Branch( "pos_2FibClust_corr_hodoX2", &pos_2FibClust_corr_hodoX2, "pos_2FibClust_corr_hodoX2/F" );

   float pos_2FibClust_hodoY2;
   outTree->Branch( "pos_2FibClust_hodoY2", &pos_2FibClust_hodoY2, "pos_2FibClust_hodoY2/F" );
   float pos_2FibClust_corr_hodoY2;
   outTree->Branch( "pos_2FibClust_corr_hodoY2", &pos_2FibClust_corr_hodoY2, "pos_2FibClust_corr_hodoY2/F" );


   float cluster_pos_hodoX1;
   outTree->Branch( "cluster_pos_hodoX1", &cluster_pos_hodoX1, "cluster_pos_hodoX1/F" );
   float cluster_pos_corr_hodoX1;
   outTree->Branch( "cluster_pos_corr_hodoX1", &cluster_pos_corr_hodoX1, "cluster_pos_corr_hodoX1/F" );
   float cluster_pos_hodoX2;
   outTree->Branch( "cluster_pos_hodoX2", &cluster_pos_hodoX2, "cluster_pos_hodoX2/F" );
   float cluster_pos_corr_hodoX2;
   outTree->Branch( "cluster_pos_corr_hodoX2", &cluster_pos_corr_hodoX2, "cluster_pos_corr_hodoX2/F" );

   float cluster_pos_hodoY1;
   outTree->Branch( "cluster_pos_hodoY1", &cluster_pos_hodoY1, "cluster_pos_hodoY1/F" );
   float cluster_pos_corr_hodoY1;
   outTree->Branch( "cluster_pos_corr_hodoY1", &cluster_pos_corr_hodoY1, "cluster_pos_corr_hodoY1/F" );
   float cluster_pos_hodoY2;
   outTree->Branch( "cluster_pos_hodoY2", &cluster_pos_hodoY2, "cluster_pos_hodoY2/F" );
   float cluster_pos_corr_hodoY2;
   outTree->Branch( "cluster_pos_corr_hodoY2", &cluster_pos_corr_hodoY2, "cluster_pos_corr_hodoY2/F" );



   float wc_x;
   outTree->Branch( "wc_x", &wc_x, "wc_x/F");
   float wc_y;
   outTree->Branch( "wc_y", &wc_y, "wc_y/F");
   float wc_x_corr;
   outTree->Branch( "wc_x_corr", &wc_x_corr, "wc_x_corr/F");
   float wc_y_corr;
   outTree->Branch( "wc_y_corr", &wc_y_corr, "wc_y_corr/F");


   int nentries = tree->GetEntries();


   RunHelper::getBeamPosition( runName, xBeam, yBeam );

   std::cout << nentries << std::endl;
 
 
   for(int  iEntry=0; iEntry<nentries; ++iEntry ) {

     tree->GetEntry( iEntry );     
     if( iEntry %  10000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;


   
   std::string theBeamEnergy = Form("%.0f",BeamEnergy);
   if( runNumber > 272 && runNumber < 298){
     theBeamEnergy = "273"; //For the long position scan prior to tdc adjustment
   }else if( runNumber  < 273){
     theBeamEnergy = "259"; //For the first runs, before adjustments
   }
   //   std::cout << "The used constant file has label = "<< theBeamEnergy  << std::endl;
   

   //set the tag for calibration
   TagHelper tagHelper(tag,theBeamEnergy);
   EnergyCalibration cef3Calib(tagHelper.getCeF3FileName());
   EnergyCalibration bgoCalib(tagHelper.getBGOFileName());
   AlignmentOfficer alignOfficer(tagHelper.getAlignmentFileName());







     //BACKWARDS COMPATIBILITY///
     assignValues( cef3, *digi_charge_integrated_sub, CEF3_START_CHANNEL);      

     assignValues( cef3_corr, *digi_charge_integrated_sub, CEF3_START_CHANNEL );  
     cef3Calib.applyCalibration(cef3_corr);
     //

     assignValues( cef3_maxAmpl, *digi_max_amplitude_sub, CEF3_START_CHANNEL);
     assignValues( cef3_chaInt, *digi_charge_integrated_sub, CEF3_START_CHANNEL);

     assignValues( cef3_maxAmpl_corr1, *digi_max_amplitude_corr1, CEF3_START_CHANNEL);
     assignValues( cef3_chaInt_corr1, *digi_charge_integrated_corr1, CEF3_START_CHANNEL);
 


     //   assignValues( cef3_maxAmpl_corr2, *digi_max_amplitude_corr2, CEF3_START_CHANNEL);
     //   assignValues( cef3_chaInt_corr2, *digi_charge_integrated_corr2, CEF3_START_CHANNEL);


     assignValues( cef3_maxAmpl_bare, *digi_max_amplitude_bare, CEF3_START_CHANNEL);
     assignValues( cef3_chaInt_bare, *digi_charge_integrated_bare, CEF3_START_CHANNEL);

     assignValues( cef3_chaInt_frac10, *digi_charge_integrated_frac10, CEF3_START_CHANNEL);
     assignValues( cef3_chaInt_frac30, *digi_charge_integrated_frac30, CEF3_START_CHANNEL);
     assignValues( cef3_chaInt_frac50, *digi_charge_integrated_frac50, CEF3_START_CHANNEL);


     assignValues( cef3_maxAmpl_bare_corr, *digi_max_amplitude_bare, CEF3_START_CHANNEL);
     assignValues( cef3_chaInt_bare_corr, *digi_charge_integrated_bare, CEF3_START_CHANNEL);
     cef3Calib.applyCalibration(cef3_maxAmpl_bare_corr);
     cef3Calib.applyCalibration(cef3_chaInt_bare_corr);

     assignValues( cef3_chaInt_frac10_corr, *digi_charge_integrated_frac10, CEF3_START_CHANNEL);
     cef3Calib.applyCalibration(cef3_chaInt_frac10_corr);
     assignValues( cef3_chaInt_frac30_corr, *digi_charge_integrated_frac30, CEF3_START_CHANNEL);
    cef3Calib.applyCalibration(cef3_chaInt_frac30_corr);
     assignValues( cef3_chaInt_frac50_corr, *digi_charge_integrated_frac50, CEF3_START_CHANNEL);   
 cef3Calib.applyCalibration(cef3_chaInt_frac50_corr);

 
 

     assignValues( bgo, *BGOvalues, 0);



  
     //"Intercalibrated" values
     assignValues( cef3_maxAmpl_corr, *digi_max_amplitude_sub, CEF3_START_CHANNEL );  
     cef3Calib.applyCalibration(cef3_maxAmpl_corr);

     assignValues( cef3_chaInt_corr, *digi_charge_integrated_sub, CEF3_START_CHANNEL );  
     cef3Calib.applyCalibration(cef3_chaInt_corr);

     assignValues( cef3_maxAmpl_corr1_corr, *digi_max_amplitude_corr1, CEF3_START_CHANNEL );  
     cef3Calib.applyCalibration(cef3_maxAmpl_corr1_corr);

     assignValues( cef3_chaInt_corr1_corr, *digi_charge_integrated_corr1, CEF3_START_CHANNEL );  
     cef3Calib.applyCalibration(cef3_chaInt_corr1_corr);

     //     assignValues( cef3_maxAmpl_corr2_corr, *digi_max_amplitude_corr2, CEF3_START_CHANNEL );  
     //     cef3Calib.applyCalibration(cef3_maxAmpl_corr2_corr);

     //    assignValues( cef3_chaInt_corr2_corr, *digi_charge_integrated_corr2, CEF3_START_CHANNEL );  
     //    cef3Calib.applyCalibration(cef3_chaInt_corr2_corr);


     assignValues( bgo_corr, *BGOvalues, 0 );
     bgoCalib.applyCalibration(bgo_corr);

     std::vector<bool> hodoX1_values(HODOX1_CHANNELS, -1.);
     std::vector<bool> hodoY1_values(HODOY1_CHANNELS, -1.);
     assignValuesBool( hodoX1_values, *HODOX1, 0. );
     assignValuesBool( hodoY1_values, *HODOY1, 0. );


     std::vector<bool> hodoX2_values(HODOX2_CHANNELS, -1.);
     std::vector<bool> hodoY2_values(HODOY2_CHANNELS, -1.);
     assignValuesBool( hodoX2_values, *HODOX2, 0 );
     assignValuesBool( hodoY2_values, *HODOY2, 0 );


     std::vector<float> hodoSmallX_values(HODOSMALLX_CHANNELS, -1.);
     std::vector<float> hodoSmallY_values(HODOSMALLY_CHANNELS, -1.);
     assignValues( hodoSmallX_values, *digi_max_amplitude, HODOSMALLX_ADC_START_CHANNEL );
     assignValues( hodoSmallY_values, *digi_max_amplitude, HODOSMALLY_ADC_START_CHANNEL );

     // hodo cluster reconstruction
     int clusterMaxFibres = 4;
     doHodoReconstructionBool( hodoX1_values    , nClusters_hodoX1    , nFibres_hodoX1    , pos_hodoX1    , 0.5, clusterMaxFibres, 0. );
     doHodoReconstructionBool( hodoY1_values    , nClusters_hodoY1    , nFibres_hodoY1    , pos_hodoY1    , 0.5, clusterMaxFibres, 0. );
     doHodoReconstructionBool( hodoX2_values    , nClusters_hodoX2    , nFibres_hodoX2    , pos_hodoX2    , 0.5, clusterMaxFibres , 0.);
     doHodoReconstructionBool( hodoY2_values    , nClusters_hodoY2    , nFibres_hodoY2    , pos_hodoY2    , 0.5, clusterMaxFibres, 0. );
     doHodoReconstruction( hodoSmallX_values, nClusters_hodoSmallX, nFibres_hodoSmallX, pos_hodoSmallX, 1.0, 1, 60.);
     doHodoReconstruction( hodoSmallY_values, nClusters_hodoSmallY, nFibres_hodoSmallY, pos_hodoSmallY, 1.0, 1, 60.);

     copyArray( nClusters_hodoX1, pos_hodoX1, pos_corr_hodoX1 );
     copyArray( nClusters_hodoY1, pos_hodoY1, pos_corr_hodoY1 );
     copyArray( nClusters_hodoX2, pos_hodoX2, pos_corr_hodoX2 );
     copyArray( nClusters_hodoY2, pos_hodoY2, pos_corr_hodoY2 );
     copyArray( nClusters_hodoSmallX, pos_hodoSmallX, pos_corr_hodoSmallX );
     copyArray( nClusters_hodoSmallY, pos_hodoSmallY, pos_corr_hodoSmallY );

     alignOfficer.fix("hodoX1", nClusters_hodoX1, pos_corr_hodoX1);
     alignOfficer.fix("hodoY1", nClusters_hodoY1, pos_corr_hodoY1);
     alignOfficer.fix("hodoX2", nClusters_hodoX2, pos_corr_hodoX2);
     alignOfficer.fix("hodoY2", nClusters_hodoY2, pos_corr_hodoY2);

     s1 = SCINTvalues->at(0);
     s3 = SCINTvalues->at(1);
     s4 = SCINTvalues->at(2);
     s6 = SCINTvalues->at(3);

     xTable = TableX;
     yTable = TableY;
     HVCeF3 = CeF3HV;
     HVBGO = BGOHV;
     beamEnergy = BeamEnergy;
     angle = BeamTilt;

     nTDCHits[0] = nTdcHits->at(0);
     nTDCHits[1] = nTdcHits->at(1);
     nTDCHits[2] = nTdcHits->at(2);
     nTDCHits[3] = nTdcHits->at(3);

     wc_x = TDCreco->at(0);
     wc_y = TDCreco->at(1);

     if( runNumber>=170 ) wc_y = -wc_y;

     wc_x_corr = wc_x + alignOfficer.getOffset("wc_x");
     wc_y_corr = wc_y + alignOfficer.getOffset("wc_y");

 



     int posOf2FibClustX1=0;
     int nrOf2FibreClustersX1 = 0;
     for( int i=0; i<nClusters_hodoX1; ++i ) {
       if( nFibres_hodoX1[i]==2  )  {  
	 ++nrOf2FibreClustersX1;
	 posOf2FibClustX1 = i ; } 
     }
     if( nrOf2FibreClustersX1 == 1){
     pos_2FibClust_hodoX1 =   pos_hodoX1[ posOf2FibClustX1] ;
     }else{
       pos_2FibClust_hodoX1 =  -999 ;}

     if(nClusters_hodoX1==1){
       cluster_pos_hodoX1 = pos_hodoX1[0];
     }else if(nrOf2FibreClustersX1==1){
       cluster_pos_hodoX1 = pos_hodoX1[ posOf2FibClustX1];
     }else{ cluster_pos_hodoX1 = -999;}


     int posOf2FibClustX2=0;
     int nrOf2FibreClustersX2 = 0;
     for( int i=0; i<nClusters_hodoX2; ++i ) {
       if( nFibres_hodoX2[i]==2  )  {  
	 ++nrOf2FibreClustersX2;
	 posOf2FibClustX2 = i ; } 
     }
     if( nrOf2FibreClustersX2 == 1){
     pos_2FibClust_hodoX2 =   pos_hodoX2[ posOf2FibClustX2] ;
     }else{
       pos_2FibClust_hodoX2 =  -999 ;}

     if(nClusters_hodoX2==1){
       cluster_pos_hodoX2 = pos_hodoX2[0];
     }else if(nrOf2FibreClustersX2==1){
       cluster_pos_hodoX2 = pos_hodoX2[ posOf2FibClustX2];
     }else{ cluster_pos_hodoX2 = -999;}

     int posOf2FibClustY1=0;
     int nrOf2FibreClustersY1 = 0;
     for( int i=0; i<nClusters_hodoY1; ++i ) {
       if( nFibres_hodoY1[i]==2  )  {  
	 ++nrOf2FibreClustersY1;
	 posOf2FibClustY1 = i ; } 
     }
     if( nrOf2FibreClustersY1 == 1){
     pos_2FibClust_hodoY1 =   pos_hodoY1[ posOf2FibClustY1] ;
     }else{
       pos_2FibClust_hodoY1 =  -999 ;}
     
     if(nClusters_hodoY1==1){
       cluster_pos_hodoY1 = pos_hodoY1[0];
     }else if(nrOf2FibreClustersY1==1){
       cluster_pos_hodoY1 = pos_hodoY1[ posOf2FibClustY1];
     }else{ cluster_pos_hodoY1 = -999;}
   

     int posOf2FibClustY2=0;
     int nrOf2FibreClustersY2 = 0;
     for( int i=0; i<nClusters_hodoY2; ++i ) {
       if( nFibres_hodoY2[i]==2  )  {  
	 ++nrOf2FibreClustersY2;
	 posOf2FibClustY2 = i ; } 
     }
     if( nrOf2FibreClustersY2 == 1){
     pos_2FibClust_hodoY2 =   pos_hodoY2[ posOf2FibClustY2] ;
     }else{
       pos_2FibClust_hodoY2 =  -999 ;}

     if(nClusters_hodoY2==1){
       cluster_pos_hodoY2 = pos_hodoY2[0];
     }else if(nrOf2FibreClustersY2==1){
       cluster_pos_hodoY2 = pos_hodoY2[ posOf2FibClustY2];
     }else{ cluster_pos_hodoY2 = -999;}


     pos_2FibClust_corr_hodoX1 = pos_2FibClust_hodoX1 + alignOfficer.getOffset("hodoX1");
     pos_2FibClust_corr_hodoY1 = pos_2FibClust_hodoY1 + alignOfficer.getOffset("hodoY1");   
     pos_2FibClust_corr_hodoX2 = pos_2FibClust_hodoX2 + alignOfficer.getOffset("hodoX2");   
     pos_2FibClust_corr_hodoY2 = pos_2FibClust_hodoY2 + alignOfficer.getOffset("hodoY2");


     cluster_pos_corr_hodoX1 = cluster_pos_hodoX1 + alignOfficer.getOffset("hodoX1");
     cluster_pos_corr_hodoY1 = cluster_pos_hodoY1 + alignOfficer.getOffset("hodoY1");   
     cluster_pos_corr_hodoX2 = cluster_pos_hodoX2 + alignOfficer.getOffset("hodoX2");   
     cluster_pos_corr_hodoY2 = cluster_pos_hodoY2 + alignOfficer.getOffset("hodoY2");





     outTree->Fill();
    
   
   } // for entries

 
   outfile->cd();
   outTree->Write();
   outfile->Close();

   std::cout << "-> Analysis Tree saved in: " << outfile->GetName() << std::endl;
   return 0;

}
  


void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos ) {

  for( unsigned i=0; i<target.size(); ++i ) 
    target[i] = source[startPos+i];

}


void assignValuesBool( std::vector<bool> &target, std::vector<bool> source, unsigned int startPos ) {

  for( unsigned i=0; i<target.size(); ++i ) 
    target[i] = source[startPos+i];

}




std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax, float Cut ) {

  std::vector<HodoCluster*> clusters;

  HodoCluster* currentCluster = new HodoCluster( hodo.size(), fibreWidth );

  for( unsigned i=0; i<hodo.size(); ++i ) {

    if( hodo[i] > Cut) { // hit

      if( currentCluster->getSize() < nClusterMax ) {

        currentCluster->addFibre( i );

      } else {

        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one
        currentCluster->addFibre( i );        // get that fibre!

      }

    } else { // as soon as you find a hole
      
      if( currentCluster->getSize() > 0 ) {
     
        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one

      }

    }


  } // for fibres


  if( currentCluster->getSize()>0 )
    clusters.push_back( currentCluster ); // store last cluster


  return clusters;

}




void doHodoReconstruction( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut ) {

  std::vector<HodoCluster*> clusters = getHodoClusters( values, fibreWidth, clusterMaxFibres, Cut );

  nClusters = clusters.size();
  for( unsigned i=0; i<clusters.size(); ++i ) {
    nFibres[i] = clusters[i]->getSize();
    pos[i] = clusters[i]->getPosition();
  }

}


void copyArray( int n, float *source, float *target ) {

  for( unsigned i=0; i<n; ++i ) 
    target[i] = source[i];

}













std::vector<HodoCluster*> getHodoClustersBool( std::vector<bool> hodo, float fibreWidth, int nClusterMax, float Cut ) {

  std::vector<HodoCluster*> clusters;

  HodoCluster* currentCluster = new HodoCluster( hodo.size(), fibreWidth );

  for( unsigned i=0; i<hodo.size(); ++i ) {

    if( hodo[i] > Cut) { // hit

      if( currentCluster->getSize() < nClusterMax ) {

        currentCluster->addFibre( i );

      } else {

        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one
        currentCluster->addFibre( i );        // get that fibre!

      }

    } else { // as soon as you find a hole
      
      if( currentCluster->getSize() > 0 ) {
     
        clusters.push_back( currentCluster ); // store old one
        currentCluster = new HodoCluster( hodo.size(), fibreWidth );   // create a new one

      }

    }


  } // for fibres


  if( currentCluster->getSize()>0 )
    clusters.push_back( currentCluster ); // store last cluster


  return clusters;

}


void doHodoReconstructionBool( std::vector<bool> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut ) {

  std::vector<HodoCluster*> clusters = getHodoClustersBool( values, fibreWidth, clusterMaxFibres, Cut );

  nClusters = clusters.size();
  for( unsigned i=0; i<clusters.size(); ++i ) {
    nFibres[i] = clusters[i]->getSize();
    pos[i] = clusters[i]->getPosition();
  }

}
