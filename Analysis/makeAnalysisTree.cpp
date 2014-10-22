#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

#include "channelInfo.h"

#include "CommonTools/interface/RunHelper.h"



void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos );



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
   std::vector<float>   *ADCvalues;
   std::vector<float>   *digi_max_amplitude;
   std::vector<float>   *digi_pedestal;
   std::vector<float>   *digi_pedestal_rms;
   std::vector<float>   *digi_time_at_frac30;
   std::vector<float>   *digi_time_at_frac50;
   std::vector<float>   *digi_time_at_max;

   // List of branches
   TBranch        *b_ADCvalues;   //!
   TBranch        *b_digi_max_amplitude;   //!
   TBranch        *b_digi_pedestal;   //!
   TBranch        *b_digi_pedestal_rms;   //!
   TBranch        *b_digi_time_at_frac30;   //!
   TBranch        *b_digi_time_at_frac50;   //!
   TBranch        *b_digi_time_at_max;   //!

   // Set object pointer
   ADCvalues = 0;
   digi_max_amplitude = 0;
   digi_pedestal = 0;
   digi_pedestal_rms = 0;
   digi_time_at_frac30 = 0;
   digi_time_at_frac50 = 0;
   digi_time_at_max = 0;

   tree->SetBranchAddress("ADCvalues", &ADCvalues, &b_ADCvalues);
   tree->SetBranchAddress("digi_max_amplitude", &digi_max_amplitude, &b_digi_max_amplitude);
   tree->SetBranchAddress("digi_pedestal", &digi_pedestal, &b_digi_pedestal);
   tree->SetBranchAddress("digi_pedestal_rms", &digi_pedestal_rms, &b_digi_pedestal_rms);
   tree->SetBranchAddress("digi_time_at_frac30", &digi_time_at_frac30, &b_digi_time_at_frac30);
   tree->SetBranchAddress("digi_time_at_frac50", &digi_time_at_frac50, &b_digi_time_at_frac50);
   tree->SetBranchAddress("digi_time_at_max", &digi_time_at_max, &b_digi_time_at_max);




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

   std::vector<float> cef3( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3", &cef3 );

   std::vector<float> bgo( BGO_CHANNELS, -1. );
   outTree->Branch( "bgo", &bgo );

   float xBeam;
   outTree->Branch( "xBeam", &xBeam, "xBeam/F");
   float yBeam;
   outTree->Branch( "yBeam", &yBeam, "yBeam/F");


   std::vector<float> hodoSmallX( HODOSMALLX_CHANNELS, -1. );
   outTree->Branch( "hodoSmallX", &hodoSmallX );
   std::vector<float> hodoSmallY( HODOSMALLY_CHANNELS, -1. );
   outTree->Branch( "hodoSmallY", &hodoSmallY );

   std::vector<float> hodoX1( HODOX1_CHANNELS, -1. );
   outTree->Branch( "hodoX1", &hodoX1 );
   std::vector<float> hodoY1( HODOY1_CHANNELS, -1. );
   outTree->Branch( "hodoY1", &hodoY1 );

   std::vector<float> hodoX2( HODOX2_CHANNELS, -1. );
   outTree->Branch( "hodoX2", &hodoX2 );
   std::vector<float> hodoY2( HODOY2_CHANNELS, -1. );
   outTree->Branch( "hodoY2", &hodoY2 );


   float wc_x;
   outTree->Branch( "wc_x", &wc_x, "wc_x");
   float wc_y;
   outTree->Branch( "wc_y", &wc_y, "wc_y");



   int nentries = tree->GetEntries();
   RunHelper::getBeamPosition( runName, xBeam, yBeam );

   for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

     tree->GetEntry( iEntry );

     if( iEntry %  5000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

     assignValues( cef3, *digi_max_amplitude, CEF3_START_CHANNEL );

     assignValues( bgo, *ADCvalues, BGO_ADC_START_CHANNEL );

     assignValues( hodoX1, *ADCvalues, HODOX1_ADC_START_CHANNEL );
     assignValues( hodoY1, *ADCvalues, HODOY1_ADC_START_CHANNEL );

     assignValues( hodoX2, *ADCvalues, HODOX2_ADC_START_CHANNEL );
     assignValues( hodoY2, *ADCvalues, HODOY2_ADC_START_CHANNEL );

     assignValues( hodoSmallX, *ADCvalues, HODOSMALLX_ADC_START_CHANNEL );
     assignValues( hodoSmallY, *ADCvalues, HODOSMALLY_ADC_START_CHANNEL );

     s1 = ADCvalues->at(S1_ADC_START_CHANNEL);
     s3 = ADCvalues->at(S3_ADC_START_CHANNEL);
     s4 = ADCvalues->at(S4_ADC_START_CHANNEL);
     s6 = ADCvalues->at(S6_ADC_START_CHANNEL);

     wc_x = ADCvalues->at(WC_X_ADC_START_CHANNEL);
     wc_y = ADCvalues->at(WC_Y_ADC_START_CHANNEL);

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
