#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH2D.h"
#include "TGaxis.h"
#include "TLegend.h"

#include "DrawTools.h"
#include "HodoCluster.h"
#include "AlignmentOfficer.h"
#include "TagHelper.h"
#include "channelInfo.h"





void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos );


void doHodoReconstruction( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut  );
std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax, float Cut  );

void assignValuesBool( std::vector<bool> &target, std::vector<bool> source, unsigned int startPos );
void doHodoReconstructionBool( std::vector<bool> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres, float Cut );
std::vector<HodoCluster*> getHodoClustersBool( std::vector<bool> hodo, float fibreWidth, int nClusterMax, float Cut );


float fitAndDraw( const std::string& outputdir, TH1F* h1 );

void drawFibres( const std::string& outputdir,TH1F* h2, TH1F* h2BG );
void drawClusters(const std::string& outputdir, TH1F* h2, TH1F* h2BG );


float fitBG( const std::string& outputdir, TH1F* h1, const float offset );

void draw3Hodos( const std::string& outputdir, TH1F* h1, TH1F* h2, TH1F* h3 );
void draw2Hodos( const std::string& outputdir, TH1F* h1, TH1F* h2, const std::string nFibs );

int main( int argc, char* argv[] ) {


  DrawTools::setStyle();
  TGaxis::SetMaxDigits(3);

  std::string tag="V0";
  if( argc<2 ) {
    std::cout << "USAGE: ./alignTracking [startTag]" << std::endl;
    exit(11);
  } else {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }
  

  // this is the dir in which the offsets will be saved:
  std::string constDirName = "Alignment";
  system(Form("mkdir -p %s", constDirName.c_str()));

  std::string constFileName = constDirName+"/offsets_"+tag+".txt";
  AlignmentOfficer alignOfficer(constFileName);



  //  TFile* file = TFile::Open("data/run_central.root");
  //  TFile* file = TFile::Open("data/run_corr2.root");
  //  TFile* file = TFile::Open("data/run_481.root");
  //  TFile* file = TFile::Open("data/run_487.root");
  //TFile* file = TFile::Open("data/run_273.root");
  TFile* file = TFile::Open("data/run_398.root");
  // TFile* file = TFile::Open("data/run_corr3.root");
  TTree* tree = (TTree*)file->Get("outputTree");


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
   std::vector<float>   *digi_charge_integrated_corr2;
   std::vector<float>   *digi_max_amplitude_corr2;

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
   TBranch        *b_digi_charge_integrated_corr2;   //!
   TBranch        *b_digi_max_amplitude_corr2;   //!

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
   digi_charge_integrated_corr2 = 0;
   digi_max_amplitude_corr2 = 0;


   //Set Branch Adresses and branch pointers 
   tree->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   tree->SetBranchAddress("spillNumber", &spillNumber, &b_spillNumber);
   tree->SetBranchAddress("evtNumber", &evtNumber, &b_evtNumber);
   tree->SetBranchAddress("BGOvalues", &BGOvalues, &b_BGOvalues);
   tree->SetBranchAddress("SCINTvalues", &SCINTvalues, &b_SCINTvalues);
   tree->SetBranchAddress("TDCreco", &TDCreco, &b_TDCreco);
   tree->SetBranchAddress("digi_charge_integrated", &digi_charge_integrated, &b_digi_charge_integrated);
   tree->SetBranchAddress("digi_max_amplitude", &digi_max_amplitude, &b_digi_max_amplitude);
   tree->SetBranchAddress("digi_pedestal", &digi_pedestal, &b_digi_pedestal);
   tree->SetBranchAddress("digi_pedestal_rms", &digi_pedestal_rms, &b_digi_pedestal_rms);
   tree->SetBranchAddress("digi_time_at_frac30", &digi_time_at_frac30, &b_digi_time_at_frac30);
   tree->SetBranchAddress("digi_time_at_frac50", &digi_time_at_frac50, &b_digi_time_at_frac50);
   tree->SetBranchAddress("digi_time_at_max", &digi_time_at_max, &b_digi_time_at_max);
   tree->SetBranchAddress("HODOX1", &HODOX1, &b_HODOX1);
   tree->SetBranchAddress("HODOX2", &HODOX2, &b_HODOX2);
   tree->SetBranchAddress("HODOY1", &HODOY1, &b_HODOY1);
   tree->SetBranchAddress("HODOY2", &HODOY2, &b_HODOY2);
   tree->SetBranchAddress("TableX", &TableX, &b_TableX);
   tree->SetBranchAddress("TableY", &TableY, &b_TableY);
   tree->SetBranchAddress("CeF3HV", &CeF3HV, &b_CeF3HV);
   tree->SetBranchAddress("BGOHV", &BGOHV, &b_BGOHV);
   tree->SetBranchAddress("BeamEnergy", &BeamEnergy, &b_BeamEnergy);
   tree->SetBranchAddress("BeamTilt", &BeamTilt, &b_BeamTilt);
   tree->SetBranchAddress("IsPhysics", &IsPhysics, &b_IsPhysics);
   tree->SetBranchAddress("nTdcHits", &nTdcHits, &b_nTdcHits);
   tree->SetBranchAddress("digi_charge_integrated_sub", &digi_charge_integrated_sub, &b_digi_charge_integrated_sub);
   tree->SetBranchAddress("digi_max_amplitude_sub", &digi_max_amplitude_sub, &b_digi_max_amplitude_sub);
   tree->SetBranchAddress("digi_pedestal_sub", &digi_pedestal_sub, &b_digi_pedestal_sub);
   tree->SetBranchAddress("digi_pedestal_rms_sub", &digi_pedestal_rms_sub, &b_digi_pedestal_rms_sub);
   tree->SetBranchAddress("digi_charge_integrated_corr2", &digi_charge_integrated_corr2, &b_digi_charge_integrated_corr2);
   tree->SetBranchAddress("digi_max_amplitude_corr2", &digi_max_amplitude_corr2, &b_digi_max_amplitude_corr2);


  int nClusters_hodoX1;
  int nFibres_hodoX1[HODOX1_CHANNELS];
  float pos_hodoX1[HODOX1_CHANNELS];
  int nClusters_hodoY1;
  int nFibres_hodoY1[HODOY1_CHANNELS];
  float pos_hodoY1[HODOY1_CHANNELS];

  int nClusters_hodoX2;
  int nFibres_hodoX2[HODOX2_CHANNELS];
  float pos_hodoX2[HODOX2_CHANNELS];
  int nClusters_hodoY2;
  int nFibres_hodoY2[HODOY2_CHANNELS];
  float pos_hodoY2[HODOY2_CHANNELS];

  float wc_x;
  float wc_y;

  int nBins = 80;
  int nBinsWC = 80;
  float xMin = -20.;
  float xMax =  20.;

  TH1F* h1_wc_y_low = new TH1F("wc_y_low", "", nBinsWC, xMin, xMax);
  TH1F* h1_wc_y_hi  = new TH1F("wc_y_hi" , "", nBinsWC, xMin, xMax);
  TH1F* h1_wc_x_low = new TH1F("wc_x_low", "", nBinsWC, xMin, xMax);
  TH1F* h1_wc_x_hi  = new TH1F("wc_x_hi" , "", nBinsWC, xMin, xMax);

  TH1F* h1_hodoY1_low = new TH1F("hodoY1_low", "", nBins, xMin, xMax);
  TH1F* h1_hodoY1_hi  = new TH1F("hodoY1_hi" , "", nBins, xMin, xMax);
  TH1F* h1_hodoX1_low = new TH1F("hodoX1_low", "", nBins, xMin, xMax);
  TH1F* h1_hodoX1_hi  = new TH1F("hodoX1_hi" , "", nBins, xMin, xMax);

  TH1F* h1_hodoY2_low = new TH1F("hodoY2_low", "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_hi  = new TH1F("hodoY2_hi" , "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_low = new TH1F("hodoX2_low", "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_hi  = new TH1F("hodoX2_hi" , "", nBins, xMin, xMax);

  TH1F* h1_hodoY1_hi_n_low = new TH1F("hodoY1_hi_n_low", "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_hi_n_low = new TH1F("hodoY2_hi_n_low" , "", nBins, xMin, xMax);
  TH1F* h1_hodoX1_hi_n_low = new TH1F("hodoX1_hi_n_low", "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_hi_n_low = new TH1F("hodoX2_hi_n_low" , "", nBins, xMin, xMax);


  TH1F* h1_hodoY1_singleClust_1Fibre = new TH1F("hodoY1_singleClust_1Fibre", "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_singleClust_1Fibre = new TH1F("hodoY2_singleClust_1Fibre" , "", nBins, xMin, xMax);
  TH1F* h1_hodoX1_singleClust_1Fibre = new TH1F("hodoX1_singleClust_1Fibre", "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_singleClust_1Fibre = new TH1F("hodoX2_singleClust_1Fibre" , "", nBins, xMin, xMax);

  TH1F* h1_hodoY1_singleClust_2Fibre = new TH1F("hodoY1_singleClust_2Fibre", "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_singleClust_2Fibre = new TH1F("hodoY2_singleClust_2Fibre" , "", nBins, xMin, xMax);
  TH1F* h1_hodoX1_singleClust_2Fibre = new TH1F("hodoX1_singleClust_2Fibre", "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_singleClust_2Fibre = new TH1F("hodoX2_singleClust_2Fibre" , "", nBins, xMin, xMax);

  TH1F* h1_hodoY1_singleClust_34Fibre = new TH1F("hodoY1_singleClust_34Fibre", "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_singleClust_34Fibre = new TH1F("hodoY2_singleClust_34Fibre" , "", nBins, xMin, xMax);
  TH1F* h1_hodoX1_singleClust_34Fibre = new TH1F("hodoX1_singleClust_34Fibre", "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_singleClust_34Fibre = new TH1F("hodoX2_singleClust_34Fibre" , "", nBins, xMin, xMax);


  TH1F* h1_hodoX1_SingleEvents_1Fib = new TH1F("hodoX1_SingleEvents_1Fib" , "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_SingleEvents_1Fib = new TH1F("hodoX2_SingleEvents_1Fib" , "", nBins, xMin, xMax);
  TH1F* h1_hodoY1_SingleEvents_1Fib = new TH1F("hodoY1_SingleEvents_1Fib" , "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_SingleEvents_1Fib = new TH1F("hodoY2_SingleEvents_1Fib" , "", nBins, xMin, xMax);
 
  TH1F* h1_hodoX1_SingleEvents_2Fib = new TH1F("hodoX1_SingleEvents_2Fib" , "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_SingleEvents_2Fib = new TH1F("hodoX2_SingleEvents_2Fib" , "", nBins, xMin, xMax);
  TH1F* h1_hodoY1_SingleEvents_2Fib = new TH1F("hodoY1_SingleEvents_2Fib" , "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_SingleEvents_2Fib = new TH1F("hodoY2_SingleEvents_2Fib" , "", nBins, xMin, xMax);
   
  TH1F* h1_hodoX1_SingleEvents_34Fib = new TH1F("hodoX1_SingleEvents_34Fib" , "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_SingleEvents_34Fib = new TH1F("hodoX2_SingleEvents_34Fib" , "", nBins, xMin, xMax);
  TH1F* h1_hodoY1_SingleEvents_34Fib = new TH1F("hodoY1_SingleEvents_34Fib" , "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_SingleEvents_34Fib = new TH1F("hodoY2_SingleEvents_34Fib" , "", nBins, xMin, xMax);


  int counter_wo_single_clust_cut_X1 = 0;
  int counter_wo_single_clust_cut_X2 = 0;
  int counter_wo_single_clust_cut_Y1 = 0;
  int counter_wo_single_clust_cut_Y2 = 0;


  int nentries = tree->GetEntries();

  //THE FIRST BIIIIG LOOOOP ~(^.^)~
  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

     if( iEntry %  10000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

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
     doHodoReconstructionBool( hodoX1_values    , nClusters_hodoX1    , nFibres_hodoX1    , pos_hodoX1    , 0.5, clusterMaxFibres, 0.  );
     doHodoReconstructionBool( hodoY1_values    , nClusters_hodoY1    , nFibres_hodoY1    , pos_hodoY1    , 0.5, clusterMaxFibres, 0.  );
     doHodoReconstructionBool( hodoX2_values    , nClusters_hodoX2    , nFibres_hodoX2    , pos_hodoX2    , 0.5, clusterMaxFibres, 0.  );
     doHodoReconstructionBool( hodoY2_values    , nClusters_hodoY2    , nFibres_hodoY2    , pos_hodoY2    , 0.5, clusterMaxFibres, 0.  );

     alignOfficer.fix("hodoX1", nClusters_hodoX1, pos_hodoX1);
     alignOfficer.fix("hodoY1", nClusters_hodoY1, pos_hodoY1);
     alignOfficer.fix("hodoX2", nClusters_hodoX2, pos_hodoX2);
     alignOfficer.fix("hodoY2", nClusters_hodoY2, pos_hodoY2);

     wc_x = TDCreco->at(0);
     wc_y = TDCreco->at(1);
     if( runNumber>=170 ) wc_y = -wc_y; // temporary fix

     wc_x += alignOfficer.getOffset("wc_x");
     wc_y += alignOfficer.getOffset("wc_y");

     float hodoSmallY_low = digi_max_amplitude->at(6);
     float hodoSmallY_hi  = digi_max_amplitude->at(7);
     //float hodoSmallX_low = digi_max_amplitude->at(5);
     //float hodoSmallX_hi  = digi_max_amplitude->at(4);
     float hodoSmallX_low = digi_max_amplitude->at(4);
     float hodoSmallX_hi  = digi_max_amplitude->at(5);

     // y low
     if( hodoSmallY_low>60. ) {
       h1_wc_y_low->Fill( wc_y );
       for( int i=0; i<nClusters_hodoY1; ++i ) {
         if( true ) {  // cut away some noise
           h1_hodoY1_low->Fill( pos_hodoY1[i] );
	   h1_hodoY1_hi_n_low->Fill(pos_hodoY1[i]);
	 }
       } // for Y1
       for( int i=0; i<nClusters_hodoY2; ++i ) {
         if( true ) {  // cut away some noise
           h1_hodoY2_low->Fill( pos_hodoY2[i] );
	   h1_hodoY2_hi_n_low->Fill(pos_hodoY2[i]);
         }
       } // for Y2
     }

     // y hi
     if( hodoSmallY_hi>60. ) {
       h1_wc_y_hi->Fill( wc_y );
       for( int i=0; i<nClusters_hodoY1; ++i ) {
         if( true ) {  // cut away some noise
           h1_hodoY1_hi->Fill( pos_hodoY1[i] );
	   h1_hodoY1_hi_n_low->Fill(pos_hodoY1[i]);
         }
       } // for Y1
       for( int i=0; i<nClusters_hodoY2; ++i ) {
         if( true) {  // cut away some noise
           h1_hodoY2_hi->Fill( pos_hodoY2[i] );
	   h1_hodoY2_hi_n_low->Fill(pos_hodoY2[i]);
         }
       } // for Y2
     }

     // x low
     if( hodoSmallX_low>60. ) {
       h1_wc_x_low->Fill( wc_x );
       for( int i=0; i<nClusters_hodoX1; ++i ) {
         if( true ) {  // cut away some noise
           h1_hodoX1_low->Fill( pos_hodoX1[i] );
	   h1_hodoX1_hi_n_low->Fill(pos_hodoX1[i]);
         }
       } // for X1
       for( int i=0; i<nClusters_hodoX2; ++i ) {
         if( true ) {  // cut away some noise
           h1_hodoX2_low->Fill( pos_hodoX2[i] );
	   h1_hodoX2_hi_n_low->Fill(pos_hodoX2[i]);
         }
       } // for X2
     }

     // x hi
     if( hodoSmallX_hi>60. ) {
       h1_wc_x_hi->Fill( wc_x );
       for( int i=0; i<nClusters_hodoX1; ++i ) {
         if( true ) {  // cut away some noise
           h1_hodoX1_hi->Fill( pos_hodoX1[i] ); 
	   h1_hodoX1_hi_n_low->Fill(pos_hodoX1[i]);
         }
       } // for X1
       for( int  i=0; i<nClusters_hodoX2; ++i ) {
         if( true ) {  // cut away some noise
           h1_hodoX2_hi->Fill( pos_hodoX2[i] );
	   h1_hodoX2_hi_n_low->Fill(pos_hodoX2[i]);
         }
       } // for X2
     }




     if( hodoSmallX_hi>60. || hodoSmallX_low>60. ) { //for HODO X1
       for( int  j=0; j<nClusters_hodoX1; ++j ) {
       if((nFibres_hodoX1[j]==1 ) ){
	 counter_wo_single_clust_cut_X1++;
	 h1_hodoX1_singleClust_1Fibre->Fill(pos_hodoX1[j]);}
       if((nFibres_hodoX1[j]==2 ) ){
	 counter_wo_single_clust_cut_X1++;
	 h1_hodoX1_singleClust_2Fibre->Fill(pos_hodoX1[j]);}
       if( (nFibres_hodoX1[j]==3 ) ){
	 counter_wo_single_clust_cut_X1++;
	 h1_hodoX1_singleClust_34Fibre->Fill(pos_hodoX1[j]);}
       if( (nFibres_hodoX1[j]==4 ) ){
	 counter_wo_single_clust_cut_X1++;
	 h1_hodoX1_singleClust_34Fibre->Fill(pos_hodoX1[j]);}
       }
       if(nClusters_hodoX1==1 && nFibres_hodoX1[0]==1) h1_hodoX1_SingleEvents_1Fib->Fill(pos_hodoX1[0]);
      
       if(nClusters_hodoX1==1 && nFibres_hodoX1[0]==2) h1_hodoX1_SingleEvents_2Fib->Fill(pos_hodoX1[0]);
       
       if(nClusters_hodoX1==1 &&( nFibres_hodoX1[0]==3 || nFibres_hodoX1[0]==4) ) h1_hodoX1_SingleEvents_34Fib->Fill(pos_hodoX1[0]);


      for( int  j=0; j<nClusters_hodoX2; ++j ) {
       if( (nFibres_hodoX2[j]==1 )){
	 counter_wo_single_clust_cut_X2++;
	 h1_hodoX2_singleClust_1Fibre->Fill(pos_hodoX2[j]);} 
       if( (nFibres_hodoX2[j]==2 )){
	 counter_wo_single_clust_cut_X2++;
	 h1_hodoX2_singleClust_2Fibre->Fill(pos_hodoX2[j]);} 
       if((nFibres_hodoX2[j]==3 )){
	 counter_wo_single_clust_cut_X2++;
	 h1_hodoX2_singleClust_34Fibre->Fill(pos_hodoX2[j]);} 
       if( (nFibres_hodoX2[j]==4 )){
	 counter_wo_single_clust_cut_X2++;
	 h1_hodoX2_singleClust_34Fibre->Fill(pos_hodoX2[j]);}   
     }
       if(nClusters_hodoX2==1 && nFibres_hodoX2[0]==1) h1_hodoX2_SingleEvents_1Fib->Fill(pos_hodoX2[0]);
       if(nClusters_hodoX2==1 &&( nFibres_hodoX2[0]==3 || nFibres_hodoX2[0]==4) ) h1_hodoX2_SingleEvents_34Fib->Fill(pos_hodoX2[0]);
       if(nClusters_hodoX2==1 && nFibres_hodoX2[0]==2) h1_hodoX2_SingleEvents_2Fib->Fill(pos_hodoX2[0]);
       

     }


     if( hodoSmallY_hi>60. || hodoSmallY_low>60. ) { //for HODO X1
       for( int  j=0; j<nClusters_hodoY1; ++j ) {
       if((nFibres_hodoY1[j]==1 ) ){
	 counter_wo_single_clust_cut_Y1++;
	 h1_hodoY1_singleClust_1Fibre->Fill(pos_hodoY1[j]);}
       if((nFibres_hodoY1[j]==2 ) ){
	 counter_wo_single_clust_cut_Y1++;
	 h1_hodoY1_singleClust_2Fibre->Fill(pos_hodoY1[j]);}
       if( (nFibres_hodoY1[j]==3 ) ){
	 counter_wo_single_clust_cut_Y1++;
	 h1_hodoY1_singleClust_34Fibre->Fill(pos_hodoY1[j]);}
       if( (nFibres_hodoY1[j]==4 ) ){
	 counter_wo_single_clust_cut_Y1++;
	 h1_hodoY1_singleClust_34Fibre->Fill(pos_hodoY1[j]);}
       }
       if(nClusters_hodoY1==1 && nFibres_hodoY1[0]==1) h1_hodoY1_SingleEvents_1Fib->Fill(pos_hodoY1[0]);
       if(nClusters_hodoY1==1 &&( nFibres_hodoY1[0]==3 || nFibres_hodoY1[0]==4) ) h1_hodoY1_SingleEvents_34Fib->Fill(pos_hodoY1[0]);
       if(nClusters_hodoY1==1 && nFibres_hodoY1[0]==2) h1_hodoY1_SingleEvents_2Fib->Fill(pos_hodoY1[0]);
       

      for( int  j=0; j<nClusters_hodoY2; ++j ) {
       if( (nFibres_hodoY2[j]==1 )){
	 counter_wo_single_clust_cut_Y2++;
	 h1_hodoY2_singleClust_1Fibre->Fill(pos_hodoY2[j]);} 
       if( (nFibres_hodoY2[j]==2 )){
	 counter_wo_single_clust_cut_Y2++;
	 h1_hodoY2_singleClust_2Fibre->Fill(pos_hodoY2[j]);} 
       if((nFibres_hodoY2[j]==3 )){
	 counter_wo_single_clust_cut_Y2++;
	 h1_hodoY2_singleClust_34Fibre->Fill(pos_hodoY2[j]);} 
       if( (nFibres_hodoY2[j]==4 )){
	 counter_wo_single_clust_cut_Y2++;
	 h1_hodoY2_singleClust_34Fibre->Fill(pos_hodoY2[j]);}   
     }
       if(nClusters_hodoY2==1 && nFibres_hodoY2[0]==1) h1_hodoY2_SingleEvents_1Fib->Fill(pos_hodoY2[0]);
       if(nClusters_hodoY2==1 &&( nFibres_hodoY2[0]==3 || nFibres_hodoY2[0]==4) ) h1_hodoY2_SingleEvents_34Fib->Fill(pos_hodoY2[0]);
       if(nClusters_hodoY2==1 && nFibres_hodoY2[0]==2) h1_hodoY2_SingleEvents_2Fib->Fill(pos_hodoY2[0]);
       
     }
     /*
     // x
     if( hodoSmallX_hi>60. || hodoSmallX_low>60. ) { //for HODO X1
       if( nClusters_hodoX1==1 &&(nFibres_hodoX1[0]==1 ) ){
	 counter_wo_single_clust_cut_X1++;
	 h1_hodoX1_singleClust_1Fibre->Fill(pos_hodoX1[0]);}
       if( nClusters_hodoX1==1 &&(nFibres_hodoX1[0]==2 ) ){
	 counter_wo_single_clust_cut_X1++;
	 h1_hodoX1_singleClust_2Fibre->Fill(pos_hodoX1[0]);}
       if( nClusters_hodoX1==1 &&(nFibres_hodoX1[0]==3 ) ){
	 counter_wo_single_clust_cut_X1++;
	 h1_hodoX1_singleClust_3Fibre->Fill(pos_hodoX1[0]);}
       if( nClusters_hodoX1==1 &&(nFibres_hodoX1[0]==4 ) ){
	 counter_wo_single_clust_cut_X1++;
	 h1_hodoX1_singleClust_4Fibre->Fill(pos_hodoX1[0]);}

       if( nClusters_hodoX2==1 &&(nFibres_hodoX2[0]==1 )){
	 counter_wo_single_clust_cut_X2++;
	 h1_hodoX2_singleClust_1Fibre->Fill(pos_hodoX2[0]);} 
       if( nClusters_hodoX2==1 &&(nFibres_hodoX2[0]==2 )){
	 counter_wo_single_clust_cut_X2++;
	 h1_hodoX2_singleClust_2Fibre->Fill(pos_hodoX2[0]);} 
       if( nClusters_hodoX2==1 &&(nFibres_hodoX2[0]==3 )){
	 counter_wo_single_clust_cut_X2++;
	 h1_hodoX2_singleClust_3Fibre->Fill(pos_hodoX2[0]);} 
       if( nClusters_hodoX2==1 &&(nFibres_hodoX2[0]==4 )){
	 counter_wo_single_clust_cut_X2++;
	 h1_hodoX2_singleClust_4Fibre->Fill(pos_hodoX2[0]);}   
     }
     //Y
    if( hodoSmallY_hi>60. || hodoSmallY_low>60. ) { //for HODO Y1
      if( nClusters_hodoY1==1 &&(nFibres_hodoY1[0]==1 )){
	 counter_wo_single_clust_cut_Y1++;
	 h1_hodoY1_singleClust_1Fibre->Fill(pos_hodoY1[0]);}
      if( nClusters_hodoY1==1 &&(nFibres_hodoY1[0]==2 )){
	 counter_wo_single_clust_cut_Y1++;
	 h1_hodoY1_singleClust_2Fibre->Fill(pos_hodoY1[0]);}
      if( nClusters_hodoY1==1 &&(nFibres_hodoY1[0]==3 )){
	 counter_wo_single_clust_cut_Y1++;
	 h1_hodoY1_singleClust_3Fibre->Fill(pos_hodoY1[0]);}
      if( nClusters_hodoY1==1 &&(nFibres_hodoY1[0]==4 )){
	 counter_wo_single_clust_cut_Y1++;
	 h1_hodoY1_singleClust_4Fibre->Fill(pos_hodoY1[0]);}

      if( nClusters_hodoY2==1 &&(nFibres_hodoY2[0]==1 )){
	 counter_wo_single_clust_cut_Y2++;
	 h1_hodoY2_singleClust_1Fibre->Fill(pos_hodoY2[0]);}
      if( nClusters_hodoY2==1 &&(nFibres_hodoY2[0]==2 )){
	 counter_wo_single_clust_cut_Y2++;
	 h1_hodoY2_singleClust_2Fibre->Fill(pos_hodoY2[0]);}
      if( nClusters_hodoY2==1 &&(nFibres_hodoY2[0]==3 )){
	 counter_wo_single_clust_cut_Y2++;
	 h1_hodoY2_singleClust_3Fibre->Fill(pos_hodoY2[0]);}
      if( nClusters_hodoY2==1 &&(nFibres_hodoY2[0]==4 )){
	 counter_wo_single_clust_cut_Y2++;
	 h1_hodoY2_singleClust_4Fibre->Fill(pos_hodoY2[0]);}
    }
     */

  } // for entries of FIRST LOOOOOP


  std::string outputdir = "AlignmentStudiesPlots_"+tag;
  system( Form("mkdir -p %s", outputdir.c_str()));

  draw3Hodos(outputdir,h1_hodoX1_singleClust_1Fibre, h1_hodoX1_singleClust_2Fibre, h1_hodoX1_singleClust_34Fibre  );
  draw3Hodos(outputdir,h1_hodoX2_singleClust_1Fibre, h1_hodoX2_singleClust_2Fibre, h1_hodoX2_singleClust_34Fibre  );
  draw3Hodos(outputdir,h1_hodoY1_singleClust_1Fibre, h1_hodoY1_singleClust_2Fibre, h1_hodoY1_singleClust_34Fibre  );
  draw3Hodos(outputdir,h1_hodoY2_singleClust_1Fibre, h1_hodoY2_singleClust_2Fibre, h1_hodoY2_singleClust_34Fibre  );

  draw2Hodos( outputdir,h1_hodoX1_SingleEvents_1Fib,h1_hodoX1_singleClust_1Fibre,"1");
  draw2Hodos( outputdir,h1_hodoY2_SingleEvents_1Fib,h1_hodoY2_singleClust_1Fibre,"1");
  draw2Hodos( outputdir,h1_hodoX2_SingleEvents_1Fib,h1_hodoX2_singleClust_1Fibre,"1");
  draw2Hodos( outputdir,h1_hodoY1_SingleEvents_1Fib,h1_hodoY1_singleClust_1Fibre,"1");

  draw2Hodos( outputdir,h1_hodoX1_SingleEvents_2Fib,h1_hodoX1_singleClust_2Fibre,"2");
  draw2Hodos( outputdir,h1_hodoY2_SingleEvents_2Fib,h1_hodoY2_singleClust_2Fibre,"2");
  draw2Hodos( outputdir,h1_hodoX2_SingleEvents_2Fib,h1_hodoX2_singleClust_2Fibre,"2");
  draw2Hodos( outputdir,h1_hodoY1_SingleEvents_2Fib,h1_hodoY1_singleClust_2Fibre,"2");


  draw2Hodos( outputdir,h1_hodoX1_SingleEvents_34Fib,h1_hodoX1_singleClust_34Fibre,"3/4");
  draw2Hodos( outputdir,h1_hodoX2_SingleEvents_34Fib,h1_hodoX2_singleClust_34Fibre,"3/4");
  draw2Hodos( outputdir,h1_hodoY1_SingleEvents_34Fib,h1_hodoY1_singleClust_34Fibre,"3/4");
  draw2Hodos( outputdir,h1_hodoY2_SingleEvents_34Fib,h1_hodoY2_singleClust_34Fibre,"3/4");
  
  float offset_wc_y_low = fitAndDraw( outputdir, h1_wc_y_low );
  float offset_wc_y_hi  = fitAndDraw( outputdir, h1_wc_y_hi );
  float offset_wc_x_low = fitAndDraw( outputdir, h1_wc_x_low );
  float offset_wc_x_hi  = fitAndDraw( outputdir, h1_wc_x_hi );
  float offset_hodoY1_low = fitAndDraw( outputdir, h1_hodoY1_low );
  float offset_hodoY1_hi  = fitAndDraw( outputdir, h1_hodoY1_hi );
  float offset_hodoX1_low = fitAndDraw( outputdir, h1_hodoX1_low );
  float offset_hodoX1_hi  = fitAndDraw( outputdir, h1_hodoX1_hi );
  float offset_hodoY2_low = fitAndDraw( outputdir, h1_hodoY2_low );
  float offset_hodoY2_hi  = fitAndDraw( outputdir, h1_hodoY2_hi );
  float offset_hodoX2_low = fitAndDraw( outputdir, h1_hodoX2_low );
  float offset_hodoX2_hi  = fitAndDraw( outputdir, h1_hodoX2_hi );

  float offset_Y1_hi_n_low =  fitAndDraw( outputdir, h1_hodoY1_hi_n_low );
  float offset_X1_hi_n_low =  fitAndDraw( outputdir, h1_hodoX1_hi_n_low );
  float offset_Y2_hi_n_low =  fitAndDraw( outputdir, h1_hodoY2_hi_n_low );
  float offset_X2_hi_n_low =  fitAndDraw( outputdir, h1_hodoX2_hi_n_low );

  float BGconstY1_hi_n_low =  fitBG( outputdir, h1_hodoY1_hi_n_low, offset_Y1_hi_n_low );
  float BGconstX1_hi_n_low =  fitBG( outputdir, h1_hodoX1_hi_n_low, offset_X1_hi_n_low );
  float BGconstY2_hi_n_low =  fitBG( outputdir, h1_hodoY2_hi_n_low, offset_Y2_hi_n_low );
  float BGconstX2_hi_n_low =  fitBG( outputdir, h1_hodoX2_hi_n_low, offset_X2_hi_n_low );
  /*
  std::cout << "BG hodo X1 tot = 2*16*"<< BGconstX1_hi_n_low << " = " << 2.*16*BGconstX1_hi_n_low << std::endl;
  std::cout << "BG hodo X2 tot = 2*16*"<< BGconstX2_hi_n_low << " = " << 2.*16*BGconstX2_hi_n_low << std::endl;
  std::cout << "BG hodo Y1 tot = 2*16*"<< BGconstY1_hi_n_low << " = " << 2.*16*BGconstY1_hi_n_low << std::endl;
  std::cout << "BG hodo Y2 tot = 2*16*"<< BGconstY2_hi_n_low << " = " << 2.*16*BGconstY2_hi_n_low << std::endl;

  std::cout << "X1 Signal = " << h1_hodoX1_hi_n_low->GetEntries() - BGconstX1_hi_n_low * 2. * 16. << std::endl;
  std::cout << "X2 Signal / BG = " << (h1_hodoX1_hi_n_low->GetEntries() - BGconstX1_hi_n_low * 2. * 16.) / (BGconstX1_hi_n_low * 2. * 16.)  << std::endl;
  std::cout << "X2 Signal = " << h1_hodoX2_hi_n_low->GetEntries() - BGconstX2_hi_n_low * 2. * 16. << std::endl;
  std::cout << "X2 Signal / BG = " << (h1_hodoX2_hi_n_low->GetEntries() -BGconstX2_hi_n_low * 2. * 16.)/(BGconstX2_hi_n_low * 2. * 16.) << std::endl;
  std::cout << "Y1 Signal = " << h1_hodoY1_hi_n_low->GetEntries() - BGconstY1_hi_n_low * 2. * 16. << std::endl;
  std::cout << "Y1 Signal / BG = " << (h1_hodoY1_hi_n_low->GetEntries() -BGconstY1_hi_n_low * 2. * 16.) /( BGconstY1_hi_n_low * 2. * 16.) << std::endl;
  std::cout << "Y2 Signal = " << h1_hodoY2_hi_n_low->GetEntries() - BGconstY2_hi_n_low * 2. * 16. << std::endl;
  std::cout << "Y2 Signal / BG = " << (h1_hodoY2_hi_n_low->GetEntries() - BGconstY2_hi_n_low * 2. * 16. )/( BGconstY2_hi_n_low * 2. * 16. )  << std::endl;
  */

  float offset_wc_y   = 0.5*(offset_wc_y_low  +offset_wc_y_hi);
  float offset_wc_x   = 0.5*(offset_wc_x_low  +offset_wc_x_hi);
  float offset_hodoY1 = 0.5*(offset_hodoY1_low+offset_hodoY1_hi);
  float offset_hodoX1 = 0.5*(offset_hodoX1_low+offset_hodoX1_hi);
  float offset_hodoY2 = 0.5*(offset_hodoY2_low+offset_hodoY2_hi);
  float offset_hodoX2 = 0.5*(offset_hodoX2_low+offset_hodoX2_hi);



  float BGconstY1_hi =  fitBG( outputdir, h1_hodoY1_hi, offset_hodoY1_hi );
  float BGconstX1_hi =  fitBG( outputdir, h1_hodoX1_hi, offset_hodoX1_hi );
  float BGconstY2_hi =  fitBG( outputdir, h1_hodoY2_hi, offset_hodoY2_hi );
  float BGconstX2_hi =  fitBG( outputdir, h1_hodoX2_hi, offset_hodoX2_hi );
  float BGconst_wc_y_hi =  fitBG( outputdir, h1_wc_y_hi, offset_wc_y_hi );
  float BGconst_wc_x_hi =  fitBG( outputdir, h1_wc_x_hi, offset_wc_x_hi );

  float BGconstY1_low =  fitBG( outputdir, h1_hodoY1_low, offset_hodoY1_low );
  float BGconstX1_low =  fitBG( outputdir, h1_hodoX1_low, offset_hodoX1_low );
  float BGconstY2_low =  fitBG( outputdir, h1_hodoY2_low, offset_hodoY2_low );
  float BGconstX2_low =  fitBG( outputdir, h1_hodoX2_low, offset_hodoX2_low );
  float BGconst_wc_y_low =  fitBG( outputdir, h1_wc_y_low, offset_wc_y_low );
  float BGconst_wc_x_low =  fitBG( outputdir, h1_wc_x_low, offset_wc_x_low );

  float BG_hodoY1 = 0.5*(BGconstY1_low+BGconstY1_hi);
  float BG_hodoX1 = 0.5*(BGconstX1_low+BGconstX1_hi);
  float BG_hodoY2 = 0.5*(BGconstY2_low+BGconstY2_hi);
  float BG_hodoX2 = 0.5*(BGconstX2_low+BGconstX2_hi);
  float BGwc_y = 0.5*(BGconst_wc_y_low+BGconst_wc_y_hi);
  float BGwc_x = 0.5*(BGconst_wc_x_low+BGconst_wc_x_hi);



  std::cout << "BG hodo X1 = 2*16*"<< BG_hodoX1 << " = " << 2.*16*BG_hodoX1 << std::endl;
  std::cout << "BG hodo X2 = 2*16*"<< BG_hodoX2 << " = " << 2.*16*BG_hodoX2 << std::endl;
  std::cout << "BG hodo Y1 = 2*16*"<< BG_hodoY1 << " = " << 2.*16*BG_hodoY1 << std::endl;
  std::cout << "BG hodo Y2 = 2*16*"<< BG_hodoY2 << " = " << 2.*16*BG_hodoY2 << std::endl;
  /*
  std::cout << "offset_wc_y: "   <<  offset_wc_y << std::endl;
  std::cout << "offset_wc_x: "   <<  offset_wc_x << std::endl;
  std::cout << "offset_hodoY1: " <<  offset_hodoY1 << std::endl;
  std::cout << "offset_hodoX1: " <<  offset_hodoX1 << std::endl;
  std::cout << "offset_hodoY2: " <<  offset_hodoY2 << std::endl;
  std::cout << "offset_hodoX2: " <<  offset_hodoX2 << std::endl;
  
  std::string ofsName = constDirName + "/offsets_new.txt";
  ofstream ofs(ofsName.c_str());
  ofs << "wc_y "   << -offset_wc_y << std::endl;
  ofs << "wc_x "   << -offset_wc_x << std::endl;
  ofs << "hodoY1 " << -offset_hodoY1 << std::endl;
  ofs << "hodoX1 " << -offset_hodoX1 << std::endl;
  ofs << "hodoY2 " << -offset_hodoY2 << std::endl;
  ofs << "hodoX2 " << -offset_hodoX2 << std::endl;
  ofs.close();

  std::cout << "-> Saved offsets in: " << ofsName << std::endl;
  */




  TH1F* h1_hodoY1_tot = new TH1F("hodoY1_tot", "", nBins, xMin, xMax);
  TH1F* h1_hodoX1_tot = new TH1F("hodoX1_tot", "", nBins, xMin, xMax);

  TH1F* h1_hodoY2_tot = new TH1F("hodoY2_tot", "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_tot = new TH1F("hodoX2_tot", "", nBins, xMin, xMax);

  //Fibres
  TH1F* h1_hodoX1_Fibres = new TH1F("hodoX1_Fibres", "", 6, 0, 6);
  TH1F* h1_hodoX1_FibresBG = new TH1F("hodoX1_FibresBG", "", 6, 0, 6);
  TH1F* h1_hodoX2_Fibres = new TH1F("hodoX2_Fibres", "", 6, 0, 6);
  TH1F* h1_hodoX2_FibresBG = new TH1F("hodoX2_FibresBG", "", 6, 0, 6);
  TH1F* h1_hodoY1_Fibres = new TH1F("hodoY1_Fibres", "", 6, 0, 6);
  TH1F* h1_hodoY1_FibresBG = new TH1F("hodoY1_FibresBG", "", 6, 0, 6);
  TH1F* h1_hodoY2_Fibres = new TH1F("hodoY2_Fibres", "", 6, 0, 6);
  TH1F* h1_hodoY2_FibresBG = new TH1F("hodoY2_FibresBG", "", 6, 0, 6);
  //Positions of Clusters for nFibres == 1
  TH1F* h1_hodoX1_1Fibre = new TH1F("hodoX1_1Fibre", "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_1Fibre = new TH1F("hodoX2_1Fibre", "", nBins, xMin, xMax);
  TH1F* h1_hodoY1_1Fibre = new TH1F("hodoY1_1Fibre", "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_1Fibre = new TH1F("hodoY2_1Fibre", "", nBins, xMin, xMax);
 
  //Positions of Clusters for nFibres < 1
  TH1F* h1_hodoX1_MoreThan1Fib = new TH1F("hodoX1_MoreThan1Fib", "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_MoreThan1Fib = new TH1F("hodoX2_MoreThan1Fib", "", nBins, xMin, xMax);
  TH1F* h1_hodoY1_MoreThan1Fib = new TH1F("hodoY1_MoreThan1Fib", "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_MoreThan1Fib = new TH1F("hodoY2_MoreThan1Fib", "", nBins, xMin, xMax);
 
  //Positions of Clusters for nFibres == 2
  TH1F* h1_hodoX1_2Fibres = new TH1F("hodoX1_2Fibres", "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_2Fibres = new TH1F("hodoX2_2Fibres", "", nBins, xMin, xMax);
  TH1F* h1_hodoY1_2Fibres = new TH1F("hodoY1_2Fibres", "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_2Fibres = new TH1F("hodoY2_2Fibres", "", nBins, xMin, xMax);
 

  //Cluster
  TH1F* h1_hodoX1_Clusters = new TH1F("hodoX1_Clusters", "", 25, 0, 25);
  TH1F* h1_hodoX2_Clusters = new TH1F("hodoX2_Clusters", "", 25, 0, 25);
  TH1F* h1_hodoY1_Clusters = new TH1F("hodoY1_Clusters", "", 25, 0, 25);
  TH1F* h1_hodoY2_Clusters = new TH1F("hodoY2_Clusters", "", 25, 0, 25);

  TH1F* h1_hodoX1_2FibClusts = new TH1F("hodoX1_2FibClusts", "", 25, 0, 25);
  TH1F* h1_hodoX2_2FibClusts = new TH1F("hodoX2_2FibClusts", "", 25, 0, 25);
  TH1F* h1_hodoY1_2FibClusts = new TH1F("hodoY1_2FibClusts", "", 25, 0, 25);
  TH1F* h1_hodoY2_2FibClusts = new TH1F("hodoY2_2FibClusts", "", 25, 0, 25);

  TH1F* h1_hodoY1_hi_n_low_old_cut = new TH1F("hodoY1_hi_n_low_old_cut", "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_hi_n_low_old_cut = new TH1F("hodoY2_hi_n_low_old_cut" , "", nBins, xMin, xMax);  
  TH1F* h1_hodoX1_hi_n_low_old_cut = new TH1F("hodoX1_hi_n_low_old_cut", "", nBins, xMin, xMax);  
  TH1F* h1_hodoX2_hi_n_low_old_cut = new TH1F("hodoX2_hi_n_low_old_cut" , "", nBins, xMin, xMax);  
  TH1F* h1_wc_x_hi_n_low_old_cut = new TH1F("wc_x_hi_n_low_old_cut", "", nBins, xMin, xMax);
  TH1F* h1_wc_y_hi_n_low_old_cut = new TH1F("wc_y_hi_n_low_old_cut" , "", nBins, xMin, xMax);

  TH1F* h1_hodoY1_hi_n_low_cut = new TH1F("hodoY1_hi_n_low_cut", "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_hi_n_low_cut = new TH1F("hodoY2_hi_n_low_cut" , "", nBins, xMin, xMax);
  TH1F* h1_hodoX1_hi_n_low_cut = new TH1F("hodoX1_hi_n_low_cut", "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_hi_n_low_cut = new TH1F("hodoX2_hi_n_low_cut" , "", nBins, xMin, xMax);
  TH1F* h1_wc_x_hi_n_low_cut = new TH1F("wc_x_hi_n_low_cut", "", nBins, xMin, xMax);
  TH1F* h1_wc_y_hi_n_low_cut = new TH1F("wc_y_hi_n_low_cut" , "", nBins, xMin, xMax);

  //Signal+BG counters
  int SnB_hodoX1 = 0;
  int SnB_hodoX2 = 0;
  int SnB_hodoY1 = 0;
  int SnB_hodoY2 = 0;
  int SnB_wc_x = 0;
  int SnB_wc_y = 0;

  //Signal+BG counters for old cuuuts
  int SnB_old_cut_hodoX1 = 0;
  int SnB_old_cut_hodoX2 = 0;
  int SnB_old_cut_hodoY1 = 0;
  int SnB_old_cut_hodoY2 = 0;
  int SnB_old_cut_wc_x = 0;
  int SnB_old_cut_wc_y = 0;
  
  // nClusters_hodo ==1 || n2FibersClusters ==1  counters
  int counter_hodoX1 = 0;
  int counter_hodoX2 = 0;
  int counter_hodoY1 = 0;
  int counter_hodoY2 = 0;
  int counter_wc_x = 0;
  int counter_wc_y = 0;

  int counter_HodoALL = 0;

  int counter_ALL = 0;


  // nFibres!=1 i.e. old cut  counters
  int counter_old_hodoX1 = 0;
  int counter_old_hodoX2 = 0;
  int counter_old_hodoY1 = 0;
  int counter_old_hodoY2 = 0;
  int counter_old_wc_x = 0;
  int counter_old_wc_y = 0;

  int counter_old_HodoALL = 0;


  TH1F* h1_hodoY1_not2FibreSingleClust = new TH1F("hodoY1_not2FibreSingleClust", "", nBins, xMin, xMax);
  TH1F* h1_hodoY2_not2FibreSingleClust = new TH1F("hodoY2_not2FibreSingleClust" , "", nBins, xMin, xMax);
  TH1F* h1_hodoX1_not2FibreSingleClust = new TH1F("hodoX1_not2FibreSingleClust", "", nBins, xMin, xMax);
  TH1F* h1_hodoX2_not2FibreSingleClust = new TH1F("hodoX2_not2FibreSingleClust" , "", nBins, xMin, xMax);



  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

     if( iEntry %  10000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;
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
     doHodoReconstructionBool( hodoX1_values    , nClusters_hodoX1    , nFibres_hodoX1    , pos_hodoX1    , 0.5, clusterMaxFibres, 0.  );
     doHodoReconstructionBool( hodoY1_values    , nClusters_hodoY1    , nFibres_hodoY1    , pos_hodoY1    , 0.5, clusterMaxFibres, 0.  );
     doHodoReconstructionBool( hodoX2_values    , nClusters_hodoX2    , nFibres_hodoX2    , pos_hodoX2    , 0.5, clusterMaxFibres, 0.  );
     doHodoReconstructionBool( hodoY2_values    , nClusters_hodoY2    , nFibres_hodoY2    , pos_hodoY2    , 0.5, clusterMaxFibres, 0.  );

     alignOfficer.fix("hodoX1", nClusters_hodoX1, pos_hodoX1);
     alignOfficer.fix("hodoY1", nClusters_hodoY1, pos_hodoY1);
     alignOfficer.fix("hodoX2", nClusters_hodoX2, pos_hodoX2);
     alignOfficer.fix("hodoY2", nClusters_hodoY2, pos_hodoY2);

     float hodoSmallY_low = digi_max_amplitude->at(6);
     float hodoSmallY_hi  = digi_max_amplitude->at(7);
     //float hodoSmallX_low = digi_max_amplitude->at(5);
     //float hodoSmallX_hi  = digi_max_amplitude->at(4);
     float hodoSmallX_low = digi_max_amplitude->at(4);
     float hodoSmallX_hi  = digi_max_amplitude->at(5);

     wc_x = TDCreco->at(0);
     wc_y = TDCreco->at(1);
     if( runNumber>=170 ) wc_y = -wc_y; // temporary fix

     wc_x += alignOfficer.getOffset("wc_x");
     wc_y += alignOfficer.getOffset("wc_y");

     if(abs(wc_x) <200){ counter_old_wc_x++;}
     if(abs(wc_y) <200){ counter_old_wc_y++;}


     if(abs(wc_x) <200 &&nTdcHits->at(0)==1 && nTdcHits->at(1)==1 && nTdcHits->at(2)==1 && nTdcHits->at(3)==1 ){ counter_wc_x++;}
     if(abs(wc_y) <200 && nTdcHits->at(0)==1 && nTdcHits->at(1)==1 && nTdcHits->at(2)==1 && nTdcHits->at(3)==1){ counter_wc_y++;}

     int nrOf1FibreClustersX1 = 0;
     for( int i=0; i<nClusters_hodoX1; ++i ) {
         if( nFibres_hodoX1[i]==1  )  {  
	   ++nrOf1FibreClustersX1;
	 } 
      }
          h1_hodoX1_Clusters->Fill(nrOf1FibreClustersX1);

      if( nrOf1FibreClustersX1==1 && nClusters_hodoX1==1 ) {  // cut away some noise
	 counter_old_hodoX1++;
       }

     int nrOf2FibreClustersX1 = 0;
       for( int i=0; i<nClusters_hodoX1; ++i ) {
         if( nFibres_hodoX1[i]==2  )  {  
	   ++nrOf2FibreClustersX1;	 } 
       }
       h1_hodoX1_2FibClusts->Fill(nrOf2FibreClustersX1);
	
       //Counting good events
       if( (nClusters_hodoX1==1 || nrOf2FibreClustersX1==1 )){
	    counter_hodoX1++;}

     int nrOf1FibreClustersY1 = 0;
     for( int i=0; i<nClusters_hodoY1; ++i ) {
         if( nFibres_hodoY1[i]==1  )  {  
	   ++nrOf1FibreClustersY1; }   
     }
      if( nrOf1FibreClustersY1==1 && nClusters_hodoY1==1 ) {  // cut away some noise
	 counter_old_hodoY1++;
       }
          h1_hodoY1_Clusters->Fill(nrOf1FibreClustersY1);

     int nrOf2FibreClustersY1 = 0;
       for( int i=0; i<nClusters_hodoY1; ++i ) {
         if( nFibres_hodoY1[i]==2  )  {  
	   ++nrOf2FibreClustersY1;	 } 
       }
       h1_hodoY1_2FibClusts->Fill(nrOf2FibreClustersY1);

       //Counting good events
       if( (nClusters_hodoY1==1 || nrOf2FibreClustersY1==1 )){
	    counter_hodoY1++;       }

     int nrOf1FibreClustersX2 = 0;
     for( int i=0; i<nClusters_hodoX2; ++i ) {
         if( nFibres_hodoX2[i]==1  )  {  
	   ++nrOf1FibreClustersX2;	 } 
        }          
      if( nrOf1FibreClustersX2==1 && nClusters_hodoX2==1 ) {  // cut away some noise
	 counter_old_hodoX2++;
      }
      h1_hodoX2_Clusters->Fill(nrOf1FibreClustersX2);

     int nrOf2FibreClustersX2 = 0;
       for( int i=0; i<nClusters_hodoX2; ++i ) {
         if( nFibres_hodoX2[i]==2  )  {  
	   ++nrOf2FibreClustersX2;	 } 
       }
       h1_hodoX2_2FibClusts->Fill(nrOf2FibreClustersX2);

       //Counting good events
       if( (nClusters_hodoX2==1 || nrOf2FibreClustersX2==1 )){
	    counter_hodoX2++;       }

     int nrOf1FibreClustersY2 = 0;
     for( int i=0; i<nClusters_hodoY2; ++i ) {
         if( nFibres_hodoY2[i]==1  )  {  
	   ++nrOf1FibreClustersY2;	 }
       }
      if( nrOf1FibreClustersY2==1 && nClusters_hodoY2==1 ) {  // cut away some noise
	 counter_old_hodoY2++;
      }
          h1_hodoY2_Clusters->Fill(nrOf1FibreClustersY2);

     int nrOf2FibreClustersY2 = 0;
       for( int i=0; i<nClusters_hodoY2; ++i ) {
         if( nFibres_hodoY2[i]==2  )  {  
	   ++nrOf2FibreClustersY2;	 } 
       }
       h1_hodoY2_2FibClusts->Fill(nrOf2FibreClustersY2);

       //Counting good events
       if( (nClusters_hodoY2==1 || nrOf2FibreClustersY2==1 )){
	    counter_hodoY2++;       }

       if( ( (nClusters_hodoY2==1 || nrOf2FibreClustersY2==1) &&       (nClusters_hodoX2==1 || nrOf2FibreClustersX2==1) &&  (nClusters_hodoY1==1 || nrOf2FibreClustersY1==1) &&     (nClusters_hodoX1==1 || nrOf2FibreClustersX1==1)     ) ){
	    counter_HodoALL++;
       }
       




       //   if(abs(wc_y) <200 && abs(wc_x)<200 && nTdcHits->at(0)>0 && nTdcHits->at(1)>0 && nTdcHits->at(2)>0 && nTdcHits->at(3)>0 &&  nTdcHits->at(0)<3 && nTdcHits->at(1)<3 && nTdcHits->at(2)<3 && nTdcHits->at(3)<3 && ( nTdcHits->at(0) + nTdcHits->at(1) + nTdcHits->at(2) + nTdcHits->at(3))<7 && ( (nClusters_hodoY2==1 || nrOf2FibreClustersY2==1) &&       (nClusters_hodoX2==1 || nrOf2FibreClustersX2==1) &&  (nClusters_hodoY1==1 || nrOf2FibreClustersY1==1) &&     (nClusters_hodoX1==1 || nrOf2FibreClustersX1==1)     )       ){ counter_ALL++;}
       if(abs(wc_y) <200 && abs(wc_x)<200 && nTdcHits->at(0)>0 && nTdcHits->at(1)>0 && nTdcHits->at(2)>0 && nTdcHits->at(3)>0 &&  nTdcHits->at(0)<3 && nTdcHits->at(1)<3 && nTdcHits->at(2)<3 && nTdcHits->at(3)<3 && ( nTdcHits->at(0) + nTdcHits->at(1) + nTdcHits->at(2) + nTdcHits->at(3))<7 && ( (nClusters_hodoY2==1 || nrOf2FibreClustersY2==1) &&       (nClusters_hodoX2==1 || nrOf2FibreClustersX2==1) &&  (nClusters_hodoY1==1 || nrOf2FibreClustersY1==1) &&     (nClusters_hodoX1==1 || nrOf2FibreClustersX1==1)     )       ){ counter_ALL++;}



  

   // x
     if( hodoSmallX_hi>60. || hodoSmallX_low>60. ) { //for HODO X1
       for( int i=0; i<nClusters_hodoX1; ++i ) {
         if( ((pos_hodoX1[i]< (offset_hodoX1_low +2)) && (pos_hodoX1[i]>( offset_hodoX1_low-2.)))  ||  (pos_hodoX1[i] < ( offset_hodoX1_hi +2) && (pos_hodoX1[i]>( offset_hodoX1_hi-2.)))  )  {  
           h1_hodoX1_tot->Fill( pos_hodoX1[i] );
	   h1_hodoX1_Fibres->Fill(nFibres_hodoX1[i]);
	   } else if( pos_hodoX1[i] < offset_hodoX1_low-2 || pos_hodoX1[i] > offset_hodoX1_hi +2.   ) { 
	   h1_hodoX1_FibresBG->Fill(nFibres_hodoX1[i]); 
	 } 
	 if(nFibres_hodoX1[i] == 1){
	   h1_hodoX1_1Fibre->Fill(pos_hodoX1[i]);
	 }else if(nFibres_hodoX1[i] > 1){
	   h1_hodoX1_MoreThan1Fib->Fill(pos_hodoX1[i]);
	 }
	 if(nFibres_hodoX1[i] == 2){
	   h1_hodoX1_2Fibres->Fill(pos_hodoX1[i]);
	 }
       }

       for( int i=0; i<nClusters_hodoX2; ++i ) { // FOR HODO X2
         if( ((pos_hodoX2[i]< (offset_hodoX2_low +2)) && (pos_hodoX2[i]>( offset_hodoX2_low-2.)))  ||  (pos_hodoX2[i] < ( offset_hodoX2_hi +2) && (pos_hodoX2[i]>( offset_hodoX2_hi-2.)))  )  {  
           h1_hodoX2_tot->Fill( pos_hodoX2[i] );
	   h1_hodoX2_Fibres->Fill(nFibres_hodoX2[i]);
	   } else if( pos_hodoX2[i] < offset_hodoX2_low-2 || pos_hodoX2[i] > offset_hodoX2_hi +2.   ) { 
	   h1_hodoX2_FibresBG->Fill(nFibres_hodoX2[i]); 
	 }
	 if(nFibres_hodoX2[i] == 1){
	   h1_hodoX2_1Fibre->Fill(pos_hodoX2[i]);
	 }else if(nFibres_hodoX2[i]>1){
	   h1_hodoX2_MoreThan1Fib->Fill(pos_hodoX2[i]);
	 }
	 if(nFibres_hodoX2[i] == 2){
	   h1_hodoX2_2Fibres->Fill(pos_hodoX2[i]);
	 }
       }
   

       if( nClusters_hodoX1==1 &&(nFibres_hodoX1[0]==1 || nFibres_hodoX1[0]==3 || nFibres_hodoX1[0]==4) ){
	 h1_hodoX1_not2FibreSingleClust->Fill(pos_hodoX1[0]);}
      if( nClusters_hodoX2==1 &&(nFibres_hodoX2[0]==1 || nFibres_hodoX2[0]==3 || nFibres_hodoX2[0]==4)){
	 h1_hodoX2_not2FibreSingleClust->Fill(pos_hodoX2[0]);}    

 }

     //Y
    if( hodoSmallY_hi>60. || hodoSmallY_low>60. ) { //for HODO Y1
       for( int i=0; i<nClusters_hodoY1; ++i ) {
         if( ((pos_hodoY1[i]< (offset_hodoY1_low +2)) && (pos_hodoY1[i]>( offset_hodoY1_low-2.)))  ||  (pos_hodoY1[i] < ( offset_hodoY1_hi +2) && (pos_hodoY1[i]>( offset_hodoY1_hi-2.)))  ) 
	   {  
           h1_hodoY1_tot->Fill( pos_hodoY1[i] );
	   h1_hodoY1_Fibres->Fill(nFibres_hodoY1[i]);
	   } 

	 else if( pos_hodoY1[i] < offset_hodoY1_low-2 || pos_hodoY1[i] > offset_hodoY1_hi +2.   ) { 
	   h1_hodoY1_FibresBG->Fill(nFibres_hodoY1[i]); 
	 }
	 if(nFibres_hodoY1[i] == 1){
	   h1_hodoY1_1Fibre->Fill(pos_hodoY1[i]);
	 }else if(nFibres_hodoY1[i] > 1){
	   h1_hodoY1_MoreThan1Fib->Fill(pos_hodoY1[i]);
	 }
	 if(nFibres_hodoY1[i] == 2){
	   h1_hodoY1_2Fibres->Fill(pos_hodoY1[i]);
	 }

 }

	 for( int i=0; i<nClusters_hodoY2; ++i ) { //for Hodo Y2
         if( ((pos_hodoY2[i]< (offset_hodoY2_low +2)) && (pos_hodoY2[i]>( offset_hodoY2_low-2.)))  ||  (pos_hodoY2[i] < ( offset_hodoY2_hi +2) && (pos_hodoY2[i]>( offset_hodoY2_hi-2.)))  ) 
	   {  
           h1_hodoY2_tot->Fill( pos_hodoY2[i] );
	   h1_hodoY2_Fibres->Fill(nFibres_hodoY2[i]);
	   } 

	 else if( pos_hodoY2[i] < offset_hodoY2_low-2 || pos_hodoY2[i] > offset_hodoY2_hi +2.   ) { 
	   h1_hodoY2_FibresBG->Fill(nFibres_hodoY2[i]); 
	 }
	 if(nFibres_hodoY2[i] == 1){
	   h1_hodoY2_1Fibre->Fill(pos_hodoY2[i]);
	 }else if(nFibres_hodoY2[i] > 1){
	   h1_hodoY2_MoreThan1Fib->Fill(pos_hodoY2[i]);
	 }
	 if(nFibres_hodoY2[i] == 2){
	   h1_hodoY2_2Fibres->Fill(pos_hodoY2[i]);
	 }
       }
    if( nClusters_hodoY1==1 &&(nFibres_hodoY1[0]==1 || nFibres_hodoY1[0]==3 || nFibres_hodoY1[0]==4)){
	 h1_hodoY1_not2FibreSingleClust->Fill(pos_hodoY1[0]);}
    if( nClusters_hodoY2==1 &&(nFibres_hodoY2[0]==1 || nFibres_hodoY2[0]==3 || nFibres_hodoY2[0]==4)){
	 h1_hodoY2_not2FibreSingleClust->Fill(pos_hodoY2[0]);}
    }





         if( hodoSmallX_hi>60. || hodoSmallX_low>60. ) {
   
	   if( nTdcHits->at(0)==1 && nTdcHits->at(1)==1 && nTdcHits->at(2)==1 && nTdcHits->at(3)==1 ){
	   h1_wc_x_hi_n_low_cut->Fill(wc_x);
	
     if( ((wc_x < (offset_wc_x_low +2)) && (wc_x >( offset_wc_x_low-2.)))  ||  (wc_x  < ( offset_wc_x_hi +2) && (wc_x >( offset_wc_x_hi-2.)))  )   SnB_wc_x++;
	 }
	 int posOf2FibClustX1=0;
	 int NrOf2FibreClustersX1 = 0;
	 for( int i=0; i<nClusters_hodoX1; ++i ) {
     if( ((pos_hodoX1[i]< (offset_hodoX1_low +2)) && (pos_hodoX1[i]>( offset_hodoX1_low-2.)))  ||  (pos_hodoX1[i] < ( offset_hodoX1_hi +2) && (pos_hodoX1[i]>( offset_hodoX1_hi-2.)))  ){
	   if( nFibres_hodoX1[i]==2  )  {  
	     ++NrOf2FibreClustersX1;
	     posOf2FibClustX1 = i ; } 
     }	 }
	 if(nClusters_hodoX1==1){
	   h1_hodoX1_hi_n_low_cut->Fill( pos_hodoX1[0] );
	   SnB_hodoX1++;
}
	 else if( NrOf2FibreClustersX1 == 1){
	   h1_hodoX1_hi_n_low_cut ->Fill(  pos_hodoX1[ posOf2FibClustX1]) ;
	   SnB_hodoX1++;
	 }


   int posOf2FibClustX2=0;
   int NrOf2FibreClustersX2 = 0;
   for( int i=0; i<nClusters_hodoX2; ++i ) {
     if( ((pos_hodoX2[i]< (offset_hodoX2_low +2)) && (pos_hodoX2[i]>( offset_hodoX2_low-2.)))  ||  (pos_hodoX2[i] < ( offset_hodoX2_hi +2) && (pos_hodoX2[i]>( offset_hodoX2_hi-2.)))  ){
     if( nFibres_hodoX2[i]==2  )  {  
       ++NrOf2FibreClustersX2;
       posOf2FibClustX2 = i ; } 
   }
   }
   if( NrOf2FibreClustersX2 == 1){
     h1_hodoX2_hi_n_low_cut->Fill(   pos_hodoX2[ posOf2FibClustX2]);
     SnB_hodoX2++;
   }else if(nClusters_hodoX2==1){
     h1_hodoX2_hi_n_low_cut->Fill(   pos_hodoX2[0]) ;
     SnB_hodoX2++;
   }
	 }

   
        if( hodoSmallY_hi>60. || hodoSmallY_low>60. ) {
	  if( nTdcHits->at(0)==1 && nTdcHits->at(1)==1 && nTdcHits->at(2)==1 && nTdcHits->at(3)==1 ){
	    h1_wc_y_hi_n_low_cut->Fill(wc_y);
	    if( ((wc_y < (offset_wc_y_low +2)) && (wc_y >( offset_wc_y_low-2.)))  ||  (wc_y  < ( offset_wc_y_hi +2) && (wc_y >( offset_wc_y_hi-2.)))  )  { SnB_wc_y++;}
	 }

	  int posOf2FibClustY1=0;
	  int NrOf2FibreClustersY1 = 0;
	  for( int i=0; i<nClusters_hodoY1; ++i ) {
	    if( ((pos_hodoY1[i]< (offset_hodoY1_low +2)) && (pos_hodoY1[i]>( offset_hodoY1_low-2.)))  ||  (pos_hodoY1[i] < ( offset_hodoY1_hi +2) && (pos_hodoY1[i]>( offset_hodoY1_hi-2.)))  ){
	      if( nFibres_hodoY1[i]==2  )  {  
		++NrOf2FibreClustersY1;
		posOf2FibClustY1 = i ; } 
	    }    }
	  if( NrOf2FibreClustersY1 == 1){
	    h1_hodoY1_hi_n_low_cut->Fill(    pos_hodoY1[ posOf2FibClustY1]) ;
	    SnB_hodoY1++;
	  }else if(nClusters_hodoY1==1){
	    h1_hodoY1_hi_n_low_cut->Fill(  pos_hodoY1[0]) ;
	    SnB_hodoY1++;
	  }
     

   int posOf2FibClustY2=0;
   int NrOf2FibreClustersY2 = 0;
   for( int i=0; i<nClusters_hodoY2; ++i ) {
     if( ((pos_hodoY2[i]< (offset_hodoY2_low +2)) && (pos_hodoY2[i]>( offset_hodoY2_low-2.)))  ||  (pos_hodoY2[i] < ( offset_hodoY2_hi +2) && (pos_hodoY2[i]>( offset_hodoY2_hi-2.)))  ){
     if( nFibres_hodoY2[i]==2  )  {  
       ++NrOf2FibreClustersY2;
       posOf2FibClustY2 = i ; } 
   }  }
   if( NrOf2FibreClustersY2 == 1){
     h1_hodoY2_hi_n_low_cut->Fill(   pos_hodoY2[ posOf2FibClustY2]) ;
     SnB_hodoY2++;
   }else if(nClusters_hodoY2==1){
     h1_hodoY2_hi_n_low_cut->Fill(  pos_hodoY2[0]) ;
     SnB_hodoY2++;
   }
   
	}



	//FOR THE OLD CUUUUNTS

      if( hodoSmallX_hi>60. || hodoSmallX_low>60. ) {
    h1_wc_x_hi_n_low_old_cut->Fill(wc_x);
	
    if( ((wc_x < (offset_wc_x_low +2)) && (wc_x >( offset_wc_x_low-2.)))  ||  (wc_x  < ( offset_wc_x_hi +2) && (wc_x >( offset_wc_x_hi-2.)))  ){
   SnB_old_cut_wc_x++;
	 }

	 for( int i=0; i<nClusters_hodoX1; ++i ) {
	   if( nFibres_hodoX1[i]!=1  )  {  
	     	   h1_hodoX1_hi_n_low_old_cut->Fill( pos_hodoX1[i] ); } 

     if( ((pos_hodoX1[i]< (offset_hodoX1_low +2)) && (pos_hodoX1[i]>( offset_hodoX1_low-2.)))  ||  (pos_hodoX1[i] < ( offset_hodoX1_hi +2) && (pos_hodoX1[i]>( offset_hodoX1_hi-2.)))  ){
	   SnB_old_cut_hodoX1++;
     }	 }
	 for( int i=0; i<nClusters_hodoX2; ++i ) {
	   if( nFibres_hodoX2[i]!=1  )  {  
	     	   h1_hodoX2_hi_n_low_old_cut->Fill( pos_hodoX2[i] ); } 

     if( ((pos_hodoX2[i]< (offset_hodoX2_low +2)) && (pos_hodoX2[i]>( offset_hodoX2_low-2.)))  ||  (pos_hodoX2[i] < ( offset_hodoX2_hi +2) && (pos_hodoX2[i]>( offset_hodoX2_hi-2.)))  ){
	   SnB_old_cut_hodoX2++;
     }	 }

	 }

    if( hodoSmallY_hi>60. || hodoSmallY_low>60. ) {
    h1_wc_y_hi_n_low_old_cut->Fill(wc_y);
	
    if( ((wc_y < (offset_wc_y_low +2)) && (wc_y >( offset_wc_y_low-2.)))  ||  (wc_y  < ( offset_wc_y_hi +2) && (wc_y >( offset_wc_y_hi-2.)))  ){
   SnB_old_cut_wc_y++;
	 }

	 for( int i=0; i<nClusters_hodoY1; ++i ) {
	   if( nFibres_hodoY1[i]!=1  )  {  
	     	   h1_hodoY1_hi_n_low_old_cut->Fill( pos_hodoY1[i] ); } 

     if( ((pos_hodoY1[i]< (offset_hodoY1_low +2)) && (pos_hodoY1[i]>( offset_hodoY1_low-2.)))  ||  (pos_hodoY1[i] < ( offset_hodoY1_hi +2) && (pos_hodoY1[i]>( offset_hodoY1_hi-2.)))  ){
	   SnB_old_cut_hodoY1++;
     }	 }
	 for( int i=0; i<nClusters_hodoY2; ++i ) {
	   if( nFibres_hodoY2[i]!=1  )  {  
	     	   h1_hodoY2_hi_n_low_old_cut->Fill( pos_hodoY2[i] ); } 

     if( ((pos_hodoY2[i]< (offset_hodoY2_low +2)) && (pos_hodoY2[i]>( offset_hodoY2_low-2.)))  ||  (pos_hodoY2[i] < ( offset_hodoY2_hi +2) && (pos_hodoY2[i]>( offset_hodoY2_hi-2.)))  ){
	   SnB_old_cut_hodoY2++;
     }	 }

	 }

   




 

  }//end of loop over entries

  float offset_hodoX1_tot = fitAndDraw( outputdir, h1_hodoX1_tot );
  float offset_hodoX2_tot = fitAndDraw( outputdir, h1_hodoX2_tot );

  float offset_hodoY1_tot = fitAndDraw( outputdir, h1_hodoY1_tot );
  float offset_hodoY2_tot = fitAndDraw( outputdir, h1_hodoY2_tot );

  float offset_wo_singleClust_X1 = fitAndDraw( outputdir,h1_hodoX1_not2FibreSingleClust );
  float offset_wo_singleClust_X2 = fitAndDraw( outputdir,h1_hodoX2_not2FibreSingleClust );
  float offset_wo_singleClust_Y1 = fitAndDraw( outputdir,h1_hodoY1_not2FibreSingleClust );
  float offset_wo_singleClust_Y2 = fitAndDraw( outputdir,h1_hodoY2_not2FibreSingleClust );


  drawFibres(outputdir, h1_hodoX1_Fibres, h1_hodoX1_FibresBG);
  drawFibres(outputdir, h1_hodoX2_Fibres, h1_hodoX2_FibresBG);
  drawFibres(outputdir, h1_hodoY1_Fibres, h1_hodoY1_FibresBG);
  drawFibres(outputdir, h1_hodoY2_Fibres, h1_hodoY2_FibresBG);

  float offset_hodoX1_1Fibre = fitAndDraw(outputdir, h1_hodoX1_1Fibre);
  float offset_hodoY1_1Fibre = fitAndDraw(outputdir, h1_hodoY1_1Fibre);
  float offset_hodoY2_1Fibre = fitAndDraw(outputdir, h1_hodoY2_1Fibre);
  float offset_hodoX2_1Fibre = fitAndDraw(outputdir, h1_hodoX2_1Fibre);
  
  float offset_hodoX1_MoreThan1Fib = fitAndDraw(outputdir, h1_hodoX1_MoreThan1Fib);
  float offset_hodoY1_MoreThan1Fib = fitAndDraw(outputdir, h1_hodoY1_MoreThan1Fib);
  float offset_hodoY2_MoreThan1Fib = fitAndDraw(outputdir, h1_hodoY2_MoreThan1Fib);
  float offset_hodoX2_MoreThan1Fib = fitAndDraw(outputdir, h1_hodoX2_MoreThan1Fib);
  
  float offset_hodoX1_2Fibres = fitAndDraw(outputdir, h1_hodoX1_2Fibres);
  float offset_hodoY1_2Fibres = fitAndDraw(outputdir, h1_hodoY1_2Fibres);
  float offset_hodoY2_2Fibres = fitAndDraw(outputdir, h1_hodoY2_2Fibres);
  float offset_hodoX2_2Fibres = fitAndDraw(outputdir, h1_hodoX2_2Fibres);
 

  drawClusters(outputdir, h1_hodoX1_2FibClusts, h1_hodoX1_Clusters);
  drawClusters(outputdir, h1_hodoX2_2FibClusts, h1_hodoX2_Clusters);
  drawClusters(outputdir, h1_hodoY1_2FibClusts, h1_hodoY1_Clusters);
  drawClusters(outputdir, h1_hodoY2_2FibClusts, h1_hodoY2_Clusters);


  float offset_Y1_hi_n_low_cut =  fitAndDraw( outputdir, h1_hodoY1_hi_n_low_cut );
  float offset_X1_hi_n_low_cut =  fitAndDraw( outputdir, h1_hodoX1_hi_n_low_cut );
  float offset_Y2_hi_n_low_cut =  fitAndDraw( outputdir, h1_hodoY2_hi_n_low_cut );
  float offset_X2_hi_n_low_cut =  fitAndDraw( outputdir, h1_hodoX2_hi_n_low_cut );

  float offset_wc_x_hi_n_low_cut =  fitAndDraw( outputdir, h1_wc_x_hi_n_low_cut );
  float offset_wc_y_hi_n_low_cut =  fitAndDraw( outputdir, h1_wc_y_hi_n_low_cut );

  
  float BGconstY1_hi_n_low_cut =  fitBG( outputdir, h1_hodoY1_hi_n_low_cut, offset_Y1_hi_n_low_cut );
  float BGconstX1_hi_n_low_cut =  fitBG( outputdir, h1_hodoX1_hi_n_low_cut, offset_X1_hi_n_low_cut );
  float BGconstY2_hi_n_low_cut =  fitBG( outputdir, h1_hodoY2_hi_n_low_cut, offset_Y2_hi_n_low_cut );
  float BGconstX2_hi_n_low_cut =  fitBG( outputdir, h1_hodoX2_hi_n_low_cut, offset_X2_hi_n_low_cut );
 float BGconstwc_x_hi_n_low_cut =  fitBG( outputdir, h1_wc_x_hi_n_low_cut, offset_wc_x_hi_n_low_cut );
 float BGconstwc_y_hi_n_low_cut =  fitBG( outputdir, h1_wc_y_hi_n_low_cut, offset_wc_y_hi_n_low_cut );



 float offset_X1_hi_n_low_old_cut =  fitAndDraw( outputdir, h1_hodoX1_hi_n_low_old_cut );
 float offset_X2_hi_n_low_old_cut =  fitAndDraw( outputdir, h1_hodoX2_hi_n_low_old_cut );
 float offset_Y1_hi_n_low_old_cut =  fitAndDraw( outputdir, h1_hodoY1_hi_n_low_old_cut );
 float offset_Y2_hi_n_low_old_cut =  fitAndDraw( outputdir, h1_hodoY2_hi_n_low_old_cut );
 float offset_wc_x_hi_n_low_old_cut =  fitAndDraw( outputdir, h1_wc_x_hi_n_low_old_cut );
 float offset_wc_y_hi_n_low_old_cut =  fitAndDraw( outputdir, h1_wc_y_hi_n_low_old_cut );
  

  float BGconstY1_hi_n_low_old_cut =  fitBG( outputdir, h1_hodoY1_hi_n_low_old_cut, offset_Y1_hi_n_low_old_cut );
  float BGconstX1_hi_n_low_old_cut =  fitBG( outputdir, h1_hodoX1_hi_n_low_old_cut, offset_X1_hi_n_low_old_cut );
  float BGconstY2_hi_n_low_old_cut =  fitBG( outputdir, h1_hodoY2_hi_n_low_old_cut, offset_Y2_hi_n_low_old_cut );
  float BGconstX2_hi_n_low_old_cut =  fitBG( outputdir, h1_hodoX2_hi_n_low_old_cut, offset_X2_hi_n_low_old_cut );
  float BGconstwc_x_hi_n_low_old_cut =  fitBG( outputdir, h1_wc_x_hi_n_low_old_cut, offset_wc_x_hi_n_low_old_cut );
  float BGconstwc_y_hi_n_low_old_cut =  fitBG( outputdir, h1_wc_y_hi_n_low_old_cut, offset_wc_y_hi_n_low_old_cut );


  /*
  std::cout << "BG hodo X1 tot = 2*16*"<< BGconstX1_hi_n_low_cut << " = " << 2.*16*BGconstX1_hi_n_low_cut << std::endl;
  std::cout << "BG hodo X2 tot = 2*16*"<< BGconstX2_hi_n_low_cut << " = " << 2.*16*BGconstX2_hi_n_low_cut << std::endl;
  std::cout << "BG hodo Y1 tot = 2*16*"<< BGconstY1_hi_n_low_cut << " = " << 2.*16*BGconstY1_hi_n_low_cut << std::endl;
  std::cout << "BG hodo Y2 tot = 2*16*"<< BGconstY2_hi_n_low_cut << " = " << 2.*16*BGconstY2_hi_n_low_cut << std::endl;

  std::cout << "X1 Signal = " << h1_hodoX1_hi_n_low_cut->GetEntries() - BGconstX1_hi_n_low_cut * 2. * 16. << std::endl;
  std::cout << "X1 Signal / BG = " <<( h1_hodoX1_hi_n_low_cut->GetEntries()- BGconstX1_hi_n_low_cut * 2. * 16.) /( BGconstX1_hi_n_low_cut * 2. * 16.)  << std::endl;
  std::cout << "X2 Signal = " << h1_hodoX2_hi_n_low_cut->GetEntries() - BGconstX2_hi_n_low_cut * 2. * 16. << std::endl;
  std::cout << "X2 Signal / BG = " << (h1_hodoX2_hi_n_low_cut->GetEntries() - BGconstX2_hi_n_low_cut * 2. * 16.) / (BGconstX2_hi_n_low_cut * 2. * 16.)  << std::endl;
  std::cout << "Y1 Signal = " << h1_hodoY1_hi_n_low_cut->GetEntries() - BGconstY1_hi_n_low_cut * 2. * 16. << std::endl;
std::cout << "Y2 Signal / BG = " << (h1_hodoY1_hi_n_low_cut->GetEntries()- BGconstY1_hi_n_low_cut * 2. * 16.) / (BGconstY1_hi_n_low_cut * 2. * 16.)   << std::endl;
  std::cout << "Y2 Signal = " << h1_hodoY2_hi_n_low_cut->GetEntries() - BGconstY2_hi_n_low_cut * 2. * 16. << std::endl;
  std::cout << "Y2 Signal / BG = " << (h1_hodoY2_hi_n_low_cut->GetEntries()- BGconstY2_hi_n_low_cut * 2. * 16.) / (BGconstY2_hi_n_low_cut * 2. * 16.)   << std::endl;


  float offset_wc_x   = 0.5*(offset_wc_x_low  +offset_wc_x_hi);
  float offset_hodoY1 = 0.5*(offset_hodoY1_low+offset_hodoY1_hi);
  */
  std::cout << "( offset_hodoY1_hi- offset_hodoY1_low) = " << ( offset_hodoY1_hi- offset_hodoY1_low) << std::endl;
  std::cout << "Signal X1 = " <<  SnB_hodoX1 - (4. +( offset_hodoX1_hi- offset_hodoX1_low))*BGconstX1_hi_n_low_cut << std::endl;
  std::cout << "BG X1 = " <<  (4. +( offset_hodoX1_hi- offset_hodoX1_low)) *BGconstX1_hi_n_low_cut  << std::endl;
  std::cout << "S/B X1 = " << (  SnB_hodoX1 - (4 +( offset_hodoX1_hi- offset_hodoX1_low))*BGconstX1_hi_n_low_cut )/( (4 +( offset_hodoX1_hi- offset_hodoX1_low))*BGconstX1_hi_n_low_cut) << std::endl;

  std::cout << "Signal X2 = " <<  SnB_hodoX2 - (4. +( offset_hodoX2_hi- offset_hodoX2_low))*BGconstX2_hi_n_low_cut << std::endl;
  std::cout << "BG X2 = " <<  (4. +( offset_hodoX2_hi- offset_hodoX2_low))*BGconstX2_hi_n_low_cut  << std::endl;
  std::cout << "S/B X2 = " << (  SnB_hodoX2 - (4 +( offset_hodoX2_hi- offset_hodoX2_low))*BGconstX2_hi_n_low_cut )/( (4 +( offset_hodoX2_hi- offset_hodoX2_low))*BGconstX2_hi_n_low_cut) << std::endl;
 
  std::cout << "Signal Y1 = " <<  SnB_hodoY1 - (4. +( offset_hodoY1_hi- offset_hodoY1_low))*BGconstY1_hi_n_low_cut << std::endl;
  std::cout << "BG Y1 = " <<  (4. +( offset_hodoY1_hi- offset_hodoY1_low))*BGconstY1_hi_n_low_cut  << std::endl;
  std::cout << "S/B Y1 = " << (  SnB_hodoY1 - (4. +( offset_hodoY1_hi- offset_hodoY1_low))*BGconstY1_hi_n_low_cut )/( (4. +( offset_hodoY1_hi- offset_hodoY1_low))*BGconstY1_hi_n_low_cut) << std::endl;
 
 std::cout << "Signal Y2 = " <<  SnB_hodoY2 - (4. +( offset_hodoY2_hi- offset_hodoY2_low))*BGconstY2_hi_n_low_cut << std::endl;
  std::cout << "BG Y2 = " <<  (4. +( offset_hodoY2_hi- offset_hodoY2_low))*BGconstY2_hi_n_low_cut  << std::endl;
  std::cout << "S/B Y2 = " << (  SnB_hodoY2 - (4. +( offset_hodoY2_hi- offset_hodoY2_low))*BGconstY2_hi_n_low_cut )/( (4. +( offset_hodoY2_hi- offset_hodoY2_low))* BGconstY2_hi_n_low_cut) << std::endl;
  
std::cout << "Signal wc_x = " <<  SnB_wc_x - (4. +( offset_wc_x_hi- offset_wc_x_low))*BGconstwc_x_hi_n_low_cut << std::endl;
  std::cout << "BG wc_x = " <<  (4. +( offset_wc_x_hi- offset_wc_x_low))*BGconstwc_x_hi_n_low_cut  << std::endl;
  std::cout << "S/B wc_x = " << (  SnB_wc_x - (4. +( offset_wc_x_hi- offset_wc_x_low))*BGconstwc_x_hi_n_low_cut )/( (4. +( offset_wc_x_hi- offset_wc_x_low))*BGconstwc_x_hi_n_low_cut) << std::endl; 
 
  std::cout << "Signal wc_y = " <<  SnB_wc_y - (4. +( offset_wc_y_hi- offset_wc_y_low))*BGconstwc_y_hi_n_low_cut << std::endl;
  std::cout << "BG wc_y = " <<  (4. +( offset_wc_y_hi- offset_wc_y_low)) *BGconstwc_y_hi_n_low_cut  << std::endl;
  std::cout << "S/B wc_y = " << (  SnB_wc_y - (4. +( offset_wc_y_hi- offset_wc_y_low))*BGconstwc_y_hi_n_low_cut )/( (4. +( offset_wc_y_hi- offset_wc_y_low))*BGconstwc_y_hi_n_low_cut) << std::endl;


  //For the old cuts
  std::cout << "OOOOOLD CUUUUNTS" << std::endl;

  std::cout << "Signal X1 = " <<  SnB_hodoX1 - (4. +( offset_hodoX1_hi- offset_hodoX1_low))*BGconstX1_hi_n_low_old_cut << std::endl;
  std::cout << "BG X1 = " <<  (4. +( offset_hodoX1_hi- offset_hodoX1_low)) *BGconstX1_hi_n_low_old_cut  << std::endl;
  std::cout << "S/B X1 = " << (  SnB_hodoX1 - (4 +( offset_hodoX1_hi- offset_hodoX1_low))*BGconstX1_hi_n_low_old_cut )/( (4 +( offset_hodoX1_hi- offset_hodoX1_low))*BGconstX1_hi_n_low_old_cut) << std::endl;
 
 std::cout << "Signal X2 = " <<  SnB_hodoX2 - (4. +( offset_hodoX2_hi- offset_hodoX2_low))*BGconstX2_hi_n_low_old_cut << std::endl;
  std::cout << "BG X2 = " <<  (4. +( offset_hodoX2_hi- offset_hodoX2_low))*BGconstX2_hi_n_low_old_cut  << std::endl;
  std::cout << "S/B X2 = " << (  SnB_hodoX2 - (4 +( offset_hodoX2_hi- offset_hodoX2_low))*BGconstX2_hi_n_low_old_cut )/( (4 +( offset_hodoX2_hi- offset_hodoX2_low))*BGconstX2_hi_n_low_old_cut) << std::endl;

  std::cout << "Signal Y1 = " <<  SnB_hodoY1 - (4. +( offset_hodoY1_hi- offset_hodoY1_low))*BGconstY1_hi_n_low_old_cut << std::endl;
  std::cout << "BG Y1 = " <<  (4. +( offset_hodoY1_hi- offset_hodoY1_low))*BGconstY1_hi_n_low_old_cut  << std::endl;
  std::cout << "S/B Y1 = " << (  SnB_hodoY1 - (4. +( offset_hodoY1_hi- offset_hodoY1_low))*BGconstY1_hi_n_low_old_cut )/( (4. +( offset_hodoY1_hi- offset_hodoY1_low))*BGconstY1_hi_n_low_old_cut) << std::endl;

  std::cout << "Signal Y2 = " <<  SnB_hodoY2 - (4. +( offset_hodoY2_hi- offset_hodoY2_low))*BGconstY2_hi_n_low_old_cut << std::endl;
  std::cout << "BG Y2 = " <<  (4. +( offset_hodoY2_hi- offset_hodoY2_low))*BGconstY2_hi_n_low_old_cut  << std::endl;
  std::cout << "S/B Y2 = " << (  SnB_hodoY2 - (4. +( offset_hodoY2_hi- offset_hodoY2_low))*BGconstY2_hi_n_low_old_cut )/( (4. +( offset_hodoY2_hi- offset_hodoY2_low))* BGconstY2_hi_n_low_old_cut) << std::endl;
  std::cout << "Signal wc_x = " <<  SnB_wc_x - (4. +( offset_wc_x_hi- offset_wc_x_low))*BGconstwc_x_hi_n_low_old_cut << std::endl;
  std::cout << "BG wc_x = " <<  (4. +( offset_wc_x_hi- offset_wc_x_low))*BGconstwc_x_hi_n_low_old_cut  << std::endl;
  std::cout << "S/B wc_x = " << (  SnB_wc_x - (4. +( offset_wc_x_hi- offset_wc_x_low))*BGconstwc_x_hi_n_low_old_cut )/( (4. +( offset_wc_x_hi- offset_wc_x_low))*BGconstwc_x_hi_n_low_old_cut) << std::endl; 
 std::cout << "Signal wc_y = " <<  SnB_wc_y - (4. +( offset_wc_y_hi- offset_wc_y_low))*BGconstwc_y_hi_n_low_old_cut << std::endl;
  std::cout << "BG wc_y = " <<  (4. +( offset_wc_y_hi- offset_wc_y_low)) *BGconstwc_y_hi_n_low_old_cut  << std::endl;
  std::cout << "S/B wc_y = " << (  SnB_wc_y - (4. +( offset_wc_y_hi- offset_wc_y_low))*BGconstwc_y_hi_n_low_old_cut )/( (4. +( offset_wc_y_hi- offset_wc_y_low))*BGconstwc_y_hi_n_low_old_cut) << std::endl;



   std::cout << "Number of events with NClust==1 || n2FibresClusts ==1  = " << std::endl;

   std::cout << "New  Efficiencies " << std::endl;

   std::cout << "Hodo X1  = "<< float(counter_hodoX1)/nentries << std::endl;
   std::cout << "Hodo X2  = "<< float(counter_hodoX2)/nentries  << std::endl;
   std::cout << "Hodo Y1  = "<< float(counter_hodoY1)/nentries  << std::endl;
   std::cout << "Hodo Y2  = "<< float(counter_hodoY2)/nentries  << std::endl;
   std::cout << "WC X = "<< float(counter_wc_x)/nentries  << std::endl;
   std::cout << "WC Y = "<< float(counter_wc_y)/nentries  << std::endl;

   std::cout << "Number of events with NClust==1 || n2FibresClusts ==1 for all " << std::endl;

   std::cout << "Hodos  = "<< counter_HodoALL << std::endl;
   std::cout << "Hodos = "<< float(counter_HodoALL)/nentries << std::endl;

   std::cout << "ALLS = "<< float(counter_ALL)/nentries << std::endl;


    std::cout << "OLD the Efficiencies " << std::endl;

   std::cout << "Hodo X1  = "<< 1.-float(counter_old_hodoX1)/nentries << std::endl;
   std::cout << "Hodo X2  = "<< 1.-float(counter_old_hodoX2)/nentries  << std::endl;
   std::cout << "Hodo Y1  = "<< 1.-float(counter_old_hodoY1)/nentries  << std::endl;
   std::cout << "Hodo Y2  = "<< 1.-float(counter_old_hodoY2)/nentries  << std::endl;
   std::cout << "WC X = "<< float(counter_old_wc_x)/nentries  << std::endl;
   std::cout << "WC Y = "<< float(counter_old_wc_y)/nentries  << std::endl;



   std::cout << std::endl;
   std::cout << "Without single clust requriement Efficiencies " << std::endl;
   std::cout << "HodoX1 = " << float(counter_wo_single_clust_cut_X1)/nentries << std::endl;
   std::cout << "HodoX2 = " << float(counter_wo_single_clust_cut_X2)/nentries << std::endl;
   std::cout << "HodoY1 = " << float(counter_wo_single_clust_cut_Y1)/nentries << std::endl;
   std::cout << "HodoY2 = " << float(counter_wo_single_clust_cut_Y2)/nentries << std::endl;

  return 0;

}





void drawClusters(const std::string& outputdir, TH1F* h2, TH1F* h2BG ) {

  TCanvas* c2 = new TCanvas("c2", "", 600, 600);
  c2->cd();

  gPad->SetLogy();

  // h2BG ->Scale(1./h2BG->Integral());

  //  h2 ->Scale(1./h2->Integral());

 float yMaxBG = 5.1*h2BG->GetMaximum();
 float yMaxSignal = 5.1*h2->GetMaximum();

 float yMax2 = yMaxBG;
 if (yMaxSignal>yMaxBG) yMax2=yMaxSignal;

  TH2D* h2_axes2 = new TH2D( "axes2", "", 15, 0,15 , 10, 1, yMax2 );
  h2_axes2->SetXTitle( h2->GetName() );
  h2_axes2->SetYTitle( "Entries" );
  h2_axes2->Draw("");

  h2BG->SetFillColor(kGray);
  h2BG->SetFillStyle(3244);
  h2BG->SetLineColor(kGray);
  h2BG->SetLineWidth(2);
  h2BG->Draw("histo same");

  //  h2->SetFillColor(kOrange);
  h2->SetLineColor(kOrange);
  h2->SetLineWidth(2);

  h2->Draw("histo same");

  gPad->RedrawAxis();
  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");

  TLegend *leggo = new TLegend(0.7, 0.9, 0.9, 0.9-0.06*2);
  leggo->SetTextSize(0.035);
  leggo->SetFillColor(0);
  leggo->AddEntry( h2, "2 Fibers", "L" );
  leggo->AddEntry( h2BG, "1 Fiber", "L" );
  leggo->Draw("same");

  c2->SaveAs( Form("%s/%s.eps", outputdir.c_str() , h2->GetName()) );
  c2->SaveAs( Form("%s/%s.png", outputdir.c_str(), h2->GetName() ) );

  delete c2;
  delete h2_axes2;
  delete h2;
  delete h2BG;

}














void drawFibres(const std::string& outputdir, TH1F* h2, TH1F* h2BG ) {

  TCanvas* c2 = new TCanvas("c2", "", 600, 600);
  c2->cd();

  h2BG ->Scale(1./h2BG->Integral());

  h2 ->Scale(1./h2->Integral());

 float yMaxBG = 1.1*h2BG->GetMaximum();
 float yMaxSignal = 1.1*h2->GetMaximum();

 float yMax2 = yMaxBG;
 if (yMaxSignal>yMaxBG) yMax2=yMaxSignal;

  TH2D* h2_axes2 = new TH2D( "axes2", "", 6, 0,6, 10, 0., yMax2 );
  h2_axes2->SetXTitle( h2->GetName() );
  h2_axes2->SetYTitle( "Entries" );
  h2_axes2->Draw("");

  h2BG->SetFillColor(kGray);
  h2BG->SetFillStyle(3244);
  h2BG->SetLineColor(kGray);
  h2BG->SetLineWidth(2);
  h2BG->Draw("histo same");

  //  h2->SetFillColor(kOrange);
  h2->SetLineColor(kOrange);
  h2->SetLineWidth(2);

  h2->Draw("histo same");

  gPad->RedrawAxis();
  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");

  TLegend *leggo = new TLegend(0.7, 0.9, 0.9, 0.9-0.06*2);
  leggo->SetTextSize(0.035);
  leggo->SetFillColor(0);
  leggo->AddEntry( h2, "Signal", "L" );
  leggo->AddEntry( h2BG, "BG", "L" );
  leggo->Draw("same");

  c2->SaveAs( Form("%s/Fibers.eps", outputdir.c_str() ) );
  c2->SaveAs( Form("%s/%s.png", outputdir.c_str(), h2->GetName() ) );

  delete c2;
  delete h2_axes2;
  delete h2;
  delete h2BG;



}





float fitAndDraw( const std::string& outputdir, TH1F* h1 ) {


  float returnConst=0.;

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  float xMin = h1->GetXaxis()->GetXmin();
  float xMax = h1->GetXaxis()->GetXmax();
  float yMax = 1.1*h1->GetMaximum();


  TH2D* axes = new TH2D( "axes", "", 10, xMin, xMax, 10, 0., yMax );
  axes->SetXTitle(h1->GetName());
  axes->SetYTitle("Events");
  axes->Draw();

  h1->SetFillColor(29);
  h1->Draw("same");

  int opt=3;

  if( opt==1 ) { // gaussian fit

    float maxPos = h1->GetBinCenter(h1->GetMaximumBin());

    TF1* f1 = new TF1(Form("f1_%s", h1->GetName()), "gaus", -10., 10.);
    
    f1->SetParameter( 1, maxPos );
    h1->Fit(f1, "QRN");
    
    for( unsigned i=0; i<4; ++i ) {

      float m = f1->GetParameter(1);
      float s = f1->GetParameter(2);
      float nSigmas = 1.2;
      f1->SetRange( m-nSigmas*s, m+nSigmas*s );
      if( i==3 ) 
        h1->Fit(f1, "RQ");
      else
        h1->Fit(f1, "RNQ");

    }


    f1->SetLineColor(kRed);
    returnConst = f1->GetParameter(1);

  } else if( opt==2 ) { // average of high bins

    float maximum = h1->GetMaximum();
    float thresh = 0.25*maximum;

    float total=0.;
    float denom = 0.;
    
    int nBins = h1->GetNbinsX();
    TH1F* h1_usedBins = new TH1F( Form("usedBins_%s", h1->GetName()), "", nBins, xMin, xMax );

    for( unsigned ibin=1; ibin<nBins; ++ibin ) {

      if( h1->GetBinContent(ibin)<thresh ) continue;

      h1_usedBins->SetBinContent( ibin, h1->GetBinContent(ibin) );
  
      // average:
      total += h1->GetBinCenter(ibin)*h1->GetBinContent(ibin);
      denom += h1->GetBinContent(ibin);

      //// weighted average:
      //total += h1->GetBinCenter(ibin)*h1->GetBinContent(ibin);
      //denom += h1->GetBinContent(ibin);

    }

    TLine* lineThresh = new TLine( xMin, thresh, xMax, thresh );
    lineThresh->SetLineStyle(2);
    lineThresh->Draw("same");

    returnConst = total/denom;
    h1_usedBins->SetFillColor(kOrange);
    h1_usedBins->Draw("same");


  } else if( opt==3 ) { // maximum

    float maxBinCenter = h1->GetBinCenter(h1->GetMaximumBin());

    returnConst = maxBinCenter;

  }


  TLine* lineOffset = new TLine( returnConst, 0., returnConst, yMax );
  lineOffset->SetLineColor(kRed);
  lineOffset->SetLineWidth(3);
  lineOffset->Draw("same");

  TPaveText* offsetText = new TPaveText( 0.6, 0.7, 0.9, 0.9, "brNDC");
  offsetText->SetFillColor(0);
  offsetText->SetTextSize(0.035);
  if( returnConst>0. )
    offsetText->AddText( Form("offset = +%.1f mm", returnConst) );
  else
    offsetText->AddText( Form("offset = %.1f mm", returnConst) );
  offsetText->Draw("same");
  
  
  gPad->RedrawAxis();

  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");

  c1->SaveAs(Form("%s/%s.eps", outputdir.c_str(), h1->GetName()));
  c1->SaveAs(Form("%s/%s.png", outputdir.c_str(), h1->GetName()));

  delete c1;
  delete axes;

  return returnConst;

}




void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos ) {

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





void assignValuesBool( std::vector<bool> &target, std::vector<bool> source, unsigned int startPos ) {

  for( unsigned i=0; i<target.size(); ++i ) 
    target[i] = source[startPos+i];

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














float fitBG( const std::string& outputdir, TH1F* h1, const float offset ) {

  float returnConst=0.;

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  float xMin = h1->GetXaxis()->GetXmin();
  float xMax = h1->GetXaxis()->GetXmax();
  float yMax = 1.1*h1->GetMaximum();

  TH2D* axes = new TH2D( "axes", "", 10, xMin, xMax, 10, 0., yMax );
  axes->SetXTitle(h1->GetName());
  axes->SetYTitle("Events");
  axes->Draw();

  h1->SetFillColor(29);
  h1->Draw("same");

  TF1* f1_low = new TF1(Form("f1_low_%s", h1->GetName()), "[0]", -20., offset-3.);
  TF1* f1_up = new TF1(Form("f1_up_%s", h1->GetName()), "[0]", offset+3, 20.);
  
  //  f1->SetParameter( 1, maxPos );
  h1->Fit(f1_up, "QRN");
  h1->Fit(f1_low, "QRN");
  
  returnConst = 0.5* (f1_low->GetParameter(0) + f1_up->GetParameter(0));
  
  TLine* lineOffset = new TLine( -20. , returnConst, 20., returnConst  );
  lineOffset->SetLineColor(kRed);
  lineOffset->SetLineWidth(3);
  lineOffset->Draw("same");
  
  TPaveText* offsetText = new TPaveText( 0.6, 0.7, 0.9, 0.9, "brNDC");
  offsetText->SetFillColor(0);
  offsetText->SetTextSize(0.035);
  offsetText->AddText( Form("BG = %.1f ", returnConst) );
  offsetText->Draw("same");  
  
  gPad->RedrawAxis();
  
  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");
  
  c1->SaveAs(Form("%s/BGFit%s.eps", outputdir.c_str(), h1->GetName()));
  c1->SaveAs(Form("%s/BGFit%s.png", outputdir.c_str(), h1->GetName()));
  
  delete c1;
  delete axes;
  return returnConst;  
}
















void draw3Hodos( const std::string& outputdir, TH1F* h1, TH1F* h2, TH1F* h3 ) {

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  //gPad->SetLogy();
  float xMin = h1->GetXaxis()->GetXmin();
  float xMax = h1->GetXaxis()->GetXmax();
  float  yMax = 1.1 * h2->GetMaximum();

  TH2D* axes = new TH2D( "axes", "", 10, xMin, xMax, 10, 0.1, yMax );
  axes->SetXTitle("Cluster Position [mm]");
  axes->SetYTitle("Events");
  axes->Draw();

  h1->SetFillColor(kGray+2);
  h1->SetLineColor(kGray+2);
  h1->SetFillStyle(3001);
  h2->SetFillColor(kOrange);
  h2->SetLineColor(kOrange);
  h3->SetFillColor(kBlue);
  h3->SetLineColor(kBlue);

  h2->Draw("same");
  h1->Draw("same");
  h3->Draw("same");
 
  gPad->RedrawAxis();

  TLegend* legend = new TLegend( 0.2, 0.9-4.*0.06, 0.5, 0.9 );
  legend->SetTextSize(0.035);
  legend->SetFillColor(0);
  legend->AddEntry( h1, "1 Fibre", "L" );
  legend->AddEntry( h2, "2 Fibre", "L" );
  legend->AddEntry( h3, "3+4 Fibre", "L" );
  legend->Draw("same"); 

  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");

  c1->SaveAs(Form("%s/SingleClusters%s.eps", outputdir.c_str(), h1->GetName()));
  c1->SaveAs(Form("%s/SingleClusters%s.png", outputdir.c_str(), h1->GetName()));
  c1->SaveAs(Form("%s/SingleClusters%s.pdf", outputdir.c_str(), h1->GetName()));

  delete c1;
  delete axes;
}



void draw2Hodos( const std::string& outputdir, TH1F* h1, TH1F* h2, const std::string nFibs ) {

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  //gPad->SetLogy();
  float xMin = h1->GetXaxis()->GetXmin();
  float xMax = h1->GetXaxis()->GetXmax();
  float  yMax = 1.1 * h2->GetMaximum();

  TH2D* axes = new TH2D( "axes", "", 10, xMin, xMax, 10, 0.1, yMax );
  axes->SetXTitle("Cluster Position [mm]");
  axes->SetYTitle("Events");
  axes->Draw();

  h2->SetFillColor(kGray+2);
  h2->SetLineColor(kGray+2);
  h2->SetFillStyle(3001);
  h1->SetFillColor(kOrange);
  h1->SetLineColor(kOrange);

  h2->Draw("same");
  h1->Draw("same");

 
  gPad->RedrawAxis();

  TLegend* legend = new TLegend( 0.2, 0.9-2.*0.06, 0.5, 0.9 );
  legend->SetTextSize(0.035);
  legend->SetFillColor(0);
  legend->AddEntry( h2, Form("All Clusters with nFibres=%s",nFibs.c_str() ), "L" );
  legend->AddEntry( h1, Form("Single Clusters with nFibres=%s",nFibs.c_str() ), "L" );

  legend->Draw("same"); 

  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");

  c1->SaveAs(Form("%s/CSC_%s.eps", outputdir.c_str(), h1->GetName()));
  c1->SaveAs(Form("%s/CSC_%s.png", outputdir.c_str(), h1->GetName()));
  c1->SaveAs(Form("%s/CSC_%s.pdf", outputdir.c_str(), h1->GetName()));

  delete c1;
  delete axes;
}


