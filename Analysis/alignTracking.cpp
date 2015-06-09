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


  //  TFile* file = TFile::Open("data/run_487.root");
  // TFile* file = TFile::Open("data/run_428.root"); //200 GeV
  //  TFile* file = TFile::Open("data/run_407.root"); //150 GeV
  //  TFile* file = TFile::Open("data/run_418.root"); //100 GeV  
//TFile* file = TFile::Open("data/run_398.root"); //50 GeV
  // TFile* file = TFile::Open("data/run_455.root"); //20 GeV
  // TFile* file = TFile::Open("data/run_508.root"); //20 GeV seems to be better than run 455
  //  TFile* file = TFile::Open("data/run_457.root"); //15 GeV
  //TFile* file = TFile::Open("data/run_431.root"); //10 GeV
  TFile* file = TFile::Open("data/run_273.root");
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


  //int nBins = 80*2*2*2*2;
  int nBins = 80;
  //   int nBinsWC = 40*5;
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


  int nentries = tree->GetEntries();

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

     /* 
     // y low
     if( hodoSmallY_low>60. ) {
       // if( (( nTdcHits->at(0)==1 && nTdcHits->at(1)==1 ) ||   (nTdcHits->at(0)==1 && nTdcHits->at(2)==1 ) ||  (nTdcHits->at(0)==1 && nTdcHits->at(3)==1) || ( nTdcHits->at(1)==1 && nTdcHits->at(2)==1 ) || ( nTdcHits->at(1)==1 && nTdcHits->at(3)==1 ) || ( nTdcHits->at(2)==1 && nTdcHits->at(3)==1 )) && (  nTdcHits->at(0)<5 && nTdcHits->at(1)<5 &&  nTdcHits->at(2)<5 && nTdcHits->at(3)<5 )  ){       h1_wc_y_low->Fill( wc_y );}
       if( nTdcHits->at(0)==1 && nTdcHits->at(1)==1 && nTdcHits->at(2)==1 && nTdcHits->at(3)==1 ){       h1_wc_y_low->Fill( wc_y );}
       int posOf2FibClustY1Low = 0;
       int posOf2FibClustY2Low=0;

	 int nrOf2FibreClustersY1Low = 0;
	 for( int i=0; i<nClusters_hodoY1; ++i ) {
	   if( nFibres_hodoY1[i]==2  )  {  
	     ++nrOf2FibreClustersY1Low;
	     posOf2FibClustY1Low= i ; } 
	 }
	 if(nrOf2FibreClustersY1Low==1|| nClusters_hodoY1==1  )
	   h1_hodoY1_low->Fill( pos_hodoY1[ posOf2FibClustY1Low] );

	  //And for Y2 LOW:
     int nrOf2FibreClustersY2Low = 0;
       for( int i=0; i<nClusters_hodoY2; ++i ) {
         if( nFibres_hodoY2[i]==2  )  {  
	   ++nrOf2FibreClustersY2Low;
	   posOf2FibClustY2Low=i;  } 
       }
       if(nrOf2FibreClustersY2Low==1||nClusters_hodoY2==1)
	 h1_hodoY2_low->Fill( pos_hodoY2[posOf2FibClustY2Low] );
     }
       
 
    
    // y hi
     if( hodoSmallY_hi>60. ) {
       if( nTdcHits->at(0)==1 && nTdcHits->at(1)==1 && nTdcHits->at(2)==1 && nTdcHits->at(3)==1 ){       h1_wc_y_hi->Fill( wc_y );}
     int posOf2FibClustY1Hi = 0;
    int   posOf2FibClustY2Hi=0;

     int nrOf2FibreClustersY1Hi = 0;
       for( int i=0; i<nClusters_hodoY1; ++i ) {
         if( nFibres_hodoY1[i]==2  )  {  
	   ++nrOf2FibreClustersY1Hi;
	   posOf2FibClustY1Hi= i ;} 
       }
       if(nrOf2FibreClustersY1Hi==1 || nClusters_hodoY1==1)
	 h1_hodoY1_hi->Fill( pos_hodoY1[ posOf2FibClustY1Hi] ); 
	 
	  //And for Y2 HI:
     int nrOf2FibreClustersY2Hi = 0;
       for( int i=0; i<nClusters_hodoY2; ++i ) {
         if( nFibres_hodoY2[i]==2  )  {  
	   ++nrOf2FibreClustersY2Hi;
	   posOf2FibClustY2Hi=i;  } 
       }
       if(nrOf2FibreClustersY2Hi==1||nClusters_hodoY2==1)
	 h1_hodoY2_hi->Fill( pos_hodoY2[posOf2FibClustY2Hi] );     
     }
     */

  // y low
     if( hodoSmallY_low>60. ) {

     int posOf2FibClustY1Low = 0;
    int   posOf2FibClustY2Low=0;

     int nrOf2FibreClustersY1Low = 0;
       for( int i=0; i<nClusters_hodoY1; ++i ) {
         if( nFibres_hodoY1[i]==2  )  {  
	   ++nrOf2FibreClustersY1Low;
	   posOf2FibClustY1Low= i ;} 
       }
 
	  //And for Y2 LOW:
     int nrOf2FibreClustersY2Low = 0;
       for( int i=0; i<nClusters_hodoY2; ++i ) {
         if( nFibres_hodoY2[i]==2  )  {  
	   ++nrOf2FibreClustersY2Low;
	   posOf2FibClustY2Low=i;  } 
       }

       if((nrOf2FibreClustersY2Low==1||nClusters_hodoY2==1)&& (nrOf2FibreClustersY1Low==1 || nClusters_hodoY1==1) && abs(wc_y) <200 && nTdcHits->at(2)>0 && nTdcHits->at(3)>0 &&  nTdcHits->at(2)<3 && nTdcHits->at(3)<3 && ( nTdcHits->at(2) + nTdcHits->at(3))<4){
	 h1_wc_y_low->Fill( wc_y );
	 h1_hodoY1_low->Fill( pos_hodoY1[posOf2FibClustY1Low] );
	 h1_hodoY2_low->Fill( pos_hodoY2[posOf2FibClustY2Low] );     
     }
     }

   // y hi
     if( hodoSmallY_hi>60. ) {

     int posOf2FibClustY1Hi = 0;
    int   posOf2FibClustY2Hi=0;

     int nrOf2FibreClustersY1Hi = 0;
       for( int i=0; i<nClusters_hodoY1; ++i ) {
         if( nFibres_hodoY1[i]==2  )  {  
	   ++nrOf2FibreClustersY1Hi;
	   posOf2FibClustY1Hi= i ;} 
       }
 
	  //And for Y2 HI:
     int nrOf2FibreClustersY2Hi = 0;
       for( int i=0; i<nClusters_hodoY2; ++i ) {
         if( nFibres_hodoY2[i]==2  )  {  
	   ++nrOf2FibreClustersY2Hi;
	   posOf2FibClustY2Hi=i;  } 
       }

       if((nrOf2FibreClustersY2Hi==1||nClusters_hodoY2==1)&& (nrOf2FibreClustersY1Hi==1 || nClusters_hodoY1==1) && abs(wc_y) <200 && nTdcHits->at(2)>0 && nTdcHits->at(3)>0 &&  nTdcHits->at(2)<3 && nTdcHits->at(3)<3 && ( nTdcHits->at(2) + nTdcHits->at(3))<4){
	 h1_wc_y_hi->Fill( wc_y );
	 h1_hodoY1_hi->Fill( pos_hodoY1[ posOf2FibClustY1Hi] );
	 h1_hodoY2_hi->Fill( pos_hodoY2[posOf2FibClustY2Hi] );     
     }
     }


  // x low
     if( hodoSmallX_low>60. ) {

     int posOf2FibClustX1Low = 0;
    int   posOf2FibClustX2Low=0;

     int nrOf2FibreClustersX1Low = 0;
       for( int i=0; i<nClusters_hodoX1; ++i ) {
         if( nFibres_hodoX1[i]==2  )  {  
	   ++nrOf2FibreClustersX1Low;
	   posOf2FibClustX1Low= i ;} 
       }
 
	  //And for X2 LOW:
     int nrOf2FibreClustersX2Low = 0;
       for( int i=0; i<nClusters_hodoX2; ++i ) {
         if( nFibres_hodoX2[i]==2  )  {  
	   ++nrOf2FibreClustersX2Low;
	   posOf2FibClustX2Low=i;  } 
       }

       if((nrOf2FibreClustersX2Low==1||nClusters_hodoX2==1)&& (nrOf2FibreClustersX1Low==1 || nClusters_hodoX1==1) &&  abs(wc_x)<200 && nTdcHits->at(0)>0 && nTdcHits->at(1)>0 && nTdcHits->at(0)<3 && nTdcHits->at(1)<3 &&  ( nTdcHits->at(0) + nTdcHits->at(1) )<4){
	 h1_wc_x_low->Fill( wc_x );
	 h1_hodoX1_low->Fill( pos_hodoX1[ posOf2FibClustX1Low] );
	 h1_hodoX2_low->Fill( pos_hodoX2[posOf2FibClustX2Low] );     
     }
     }

   // x hi
     if( hodoSmallX_hi>60. ) {

     int posOf2FibClustX1Hi = 0;
    int   posOf2FibClustX2Hi=0;

     int nrOf2FibreClustersX1Hi = 0;
       for( int i=0; i<nClusters_hodoX1; ++i ) {
         if( nFibres_hodoX1[i]==2  )  {  
	   ++nrOf2FibreClustersX1Hi;
	   posOf2FibClustX1Hi= i ;} 
       }
 
	  //And for X2 HI:
     int nrOf2FibreClustersX2Hi = 0;
       for( int i=0; i<nClusters_hodoX2; ++i ) {
         if( nFibres_hodoX2[i]==2  )  {  
	   ++nrOf2FibreClustersX2Hi;
	   posOf2FibClustX2Hi=i;  } 
       }

       if((nrOf2FibreClustersX2Hi==1||nClusters_hodoX2==1)&& (nrOf2FibreClustersX1Hi==1 || nClusters_hodoX1==1) && abs(wc_x) <200 && nTdcHits->at(0)>0 && nTdcHits->at(1)>0 &&  nTdcHits->at(0)<3 && nTdcHits->at(1)<3 && ( nTdcHits->at(0) + nTdcHits->at(1) )<4){
	 h1_wc_x_hi->Fill( wc_x );
	 h1_hodoX1_hi->Fill( pos_hodoX1[ posOf2FibClustX1Hi] );
	 h1_hodoX2_hi->Fill( pos_hodoX2[posOf2FibClustX2Hi] );     
     }
     }

     /*
     // X low
     if( hodoSmallX_low>60. ) {
       if( nTdcHits->at(0)==1 && nTdcHits->at(1)==1 && nTdcHits->at(2)==1 && nTdcHits->at(3)==1 ) 
      h1_wc_x_low->Fill( wc_x );
       int posOf2FibClustX1Low = 0;
       int posOf2FibClustX2Low=0;

	 int nrOf2FibreClustersX1Low = 0;
	 for( int i=0; i<nClusters_hodoX1; ++i ) {
	   if( nFibres_hodoX1[i]==2  )  {  
	     ++nrOf2FibreClustersX1Low;
	     posOf2FibClustX1Low= i ; } 
	 }
	 if(nrOf2FibreClustersX1Low==1|| nClusters_hodoX1==1)
	   h1_hodoX1_low->Fill( pos_hodoX1[ posOf2FibClustX1Low] );

	  //And for X2 LOW:
     int nrOf2FibreClustersX2Low = 0;
       for( int i=0; i<nClusters_hodoX2; ++i ) {
         if( nFibres_hodoX2[i]==2  )  {  
	   ++nrOf2FibreClustersX2Low;
	   posOf2FibClustX2Low=i;  } 
       }
       if(nrOf2FibreClustersX2Low==1||nClusters_hodoX2==1)
	 h1_hodoX2_low->Fill( pos_hodoX2[posOf2FibClustX2Low] );
     }
       
     // X hi
     if( hodoSmallX_hi>60. ) {
       if( nTdcHits->at(0)==1 && nTdcHits->at(1)==1 && nTdcHits->at(2)==1 && nTdcHits->at(3)==1 ){       h1_wc_x_hi->Fill( wc_x );}
     int posOf2FibClustX1Hi = 0;
    int   posOf2FibClustX2Hi=0;

     int nrOf2FibreClustersX1Hi = 0;
       for( int i=0; i<nClusters_hodoX1; ++i ) {
         if( nFibres_hodoX1[i]==2  )  {  
	   ++nrOf2FibreClustersX1Hi;
	   posOf2FibClustX1Hi= i ;} 
       }
       if(nrOf2FibreClustersX1Hi==1 || nClusters_hodoX1==1)
	 h1_hodoX1_hi->Fill( pos_hodoX1[ posOf2FibClustX1Hi] ); 
	 
	  //And for X2 HI:
     int nrOf2FibreClustersX2Hi = 0;
       for( int i=0; i<nClusters_hodoX2; ++i ) {
         if( nFibres_hodoX2[i]==2  )  {  
	   ++nrOf2FibreClustersX2Hi;
	   posOf2FibClustX2Hi=i;  } 
       }
       if(nrOf2FibreClustersX2Hi==1||nClusters_hodoX2==1)
	 h1_hodoX2_hi->Fill( pos_hodoX2[posOf2FibClustX2Hi] );     
     }
     */

     /*
     // y low
     if( hodoSmallY_low>60. ) {
       h1_wc_y_low->Fill( wc_y );
       for( unsigned i=0; i<nClusters_hodoY1; ++i ) {
         if( nFibres_hodoY1[i]!=1 ) {  // cut away some noise
           h1_hodoY1_low->Fill( pos_hodoY1[i] );
         }
       } // for Y1
       for( int i=0; i<nClusters_hodoY2; ++i ) {
         if( nFibres_hodoY2[i]!=1 ) {  // cut away some noise
           h1_hodoY2_low->Fill( pos_hodoY2[i] );
         }
       } // for Y2
     }
      
     // x low
     if( hodoSmallX_low>60. ) {
       h1_wc_x_low->Fill( wc_x );
       for( int i=0; i<nClusters_hodoX1; ++i ) {
         if( nFibres_hodoX1[i]!=1 ) {  // cut away some noise
           h1_hodoX1_low->Fill( pos_hodoX1[i] );
         }
       } // for X1
       for( int i=0; i<nClusters_hodoX2; ++i ) {
         if( nFibres_hodoX2[i]!=1 ) {  // cut away some noise
           h1_hodoX2_low->Fill( pos_hodoX2[i] );
         }
       } // for X2
     }
   
     // x hi
     if( hodoSmallX_hi>60. ) {
       h1_wc_x_hi->Fill( wc_x );
       for( int i=0; i<nClusters_hodoX1; ++i ) {
         if( nFibres_hodoX1[i]!=1 ) {  // cut away some noise
           h1_hodoX1_hi->Fill( pos_hodoX1[i] );
         }
       } // for X1
       for( int i=0; i<nClusters_hodoX2; ++i ) {
         if( nFibres_hodoX2[i]!=1 ) {  // cut away some noise
           h1_hodoX2_hi->Fill( pos_hodoX2[i] );
         }
       } // for X2
     }
*/



  } // for entries


  std::string outputdir = "AlignmentPlots_"+tag;
  system( Form("mkdir -p %s", outputdir.c_str()));

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

  float offset_wc_y   = 0.5*(offset_wc_y_low  +offset_wc_y_hi);
  float offset_wc_x   = 0.5*(offset_wc_x_low  +offset_wc_x_hi);
  float offset_hodoY1 = 0.5*(offset_hodoY1_low+offset_hodoY1_hi);
  float offset_hodoX1 = 0.5*(offset_hodoX1_low+offset_hodoX1_hi);
  float offset_hodoY2 = 0.5*(offset_hodoY2_low+offset_hodoY2_hi);
  float offset_hodoX2 = 0.5*(offset_hodoX2_low+offset_hodoX2_hi);

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

  return 0;

}


float fitAndDraw( const std::string& outputdir, TH1F* h1 ) {


  float returnConst=0.;

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  float xMin = h1->GetXaxis()->GetXmin();
  float xMax = h1->GetXaxis()->GetXmax();
  float yMax = 1.1*h1->GetMaximum();


  TH2D* axes = new TH2D( "axes", "", 10, xMin, xMax, 10, 0., yMax );
  axes->SetXTitle("Cluster Position [mm]");
  //  axes->SetXTitle(h1->GetName());
  axes->SetYTitle("Events");
  axes->Draw();

  h1->SetFillColor(29);
  h1->Draw("same");

  int opt=2;

  if( opt==1 ) { // gaussian fit

    float maxPos = h1->GetBinCenter(h1->GetMaximumBin());

    TF1* f1 = new TF1(Form("f1_%s", h1->GetName()), "gaus", -10., 10.);
    
    f1->SetParameter( 1, maxPos );
    h1->Fit(f1, "QRN");
    
    for( int i=0; i<4; ++i ) {

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

    for( int ibin=1; ibin<nBins; ++ibin ) {

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
  c1->SaveAs(Form("%s/%s.pdf", outputdir.c_str(), h1->GetName()));

  delete c1;
  delete axes;

  return returnConst;

}




void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos ) {

  for( int i=0; i<target.size(); ++i ) 
    target[i] = source[startPos+i];


}



std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax, float Cut ) {

  std::vector<HodoCluster*> clusters;

  HodoCluster* currentCluster = new HodoCluster( hodo.size(), fibreWidth );

  for( int  i=0; i<hodo.size(); ++i ) {

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
  for( int i=0; i<clusters.size(); ++i ) {
    nFibres[i] = clusters[i]->getSize();
    pos[i] = clusters[i]->getPosition();
  }

}





void assignValuesBool( std::vector<bool> &target, std::vector<bool> source, unsigned int startPos ) {

  for( int  i=0; i<target.size(); ++i ) 
    target[i] = source[startPos+i];

}








std::vector<HodoCluster*> getHodoClustersBool( std::vector<bool> hodo, float fibreWidth, int nClusterMax, float Cut ) {

  std::vector<HodoCluster*> clusters;

  HodoCluster* currentCluster = new HodoCluster( hodo.size(), fibreWidth );

  for( int i=0; i<hodo.size(); ++i ) {

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
  for( int i=0; i<clusters.size(); ++i ) {
    nFibres[i] = clusters[i]->getSize();
    pos[i] = clusters[i]->getPosition();
  }

}











