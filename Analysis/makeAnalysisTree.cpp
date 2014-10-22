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
  
