#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

#include "channelInfo.h"
#include "interface/HodoCluster.h"
#include "interface/CalibrationUtility.h"
#include "interface/EnergyCalibration.h"
#include "CommonTools/interface/RunHelper.h"

#include "TApplication.h"


void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos );

void doHodoReconstruction( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres );
std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax );







int main( int argc, char* argv[] ) {

  TApplication* a = new TApplication("a", 0, 0);


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

   std::string fileName = "data/run_" + runName + ".root";
   TFile* file = TFile::Open(fileName.c_str());
   if( file==0 ) {
     std::cout << "ERROR! Din't find file " << fileName << std::endl;
     std::cout << "Exiting." << std::endl;
     exit(11);
   }
   TTree* tree = (TTree*)file->Get("outputTree");


  //set the tag for calibration
  CalibrationUtility calibUtil(tag);
  EnergyCalibration cef3Calib(calibUtil.getCeF3FileName());
  EnergyCalibration bgoCalib(calibUtil.getBGOFileName());



   // Declaration of leaf types
   UInt_t          runNumber;
   UInt_t          spillNumber;
   UInt_t          evtNumber;

   std::vector<float>   *ADCvalues;
   std::vector<float>   *digi_max_amplitude;
   std::vector<float>   *digi_charge_integrated;
   std::vector<float>   *digi_pedestal;
   std::vector<float>   *digi_pedestal_rms;
   std::vector<float>   *digi_time_at_frac30;
   std::vector<float>   *digi_time_at_frac50;
   std::vector<float>   *digi_time_at_max;

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_spillNumber;   //!
   TBranch        *b_evtNumber;   //!
   TBranch        *b_ADCvalues;   //!
   TBranch        *b_digi_max_amplitude;   //!
   TBranch        *b_digi_charge_integrated;   //!
   TBranch        *b_digi_pedestal;   //!
   TBranch        *b_digi_pedestal_rms;   //!
   TBranch        *b_digi_time_at_frac30;   //!
   TBranch        *b_digi_time_at_frac50;   //!
   TBranch        *b_digi_time_at_max;   //!

   // Set object pointer
   ADCvalues = 0;
   digi_max_amplitude = 0;
   digi_charge_integrated = 0;
   digi_pedestal = 0;
   digi_pedestal_rms = 0;
   digi_time_at_frac30 = 0;
   digi_time_at_frac50 = 0;
   digi_time_at_max = 0;



   tree->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   tree->SetBranchAddress("spillNumber", &spillNumber, &b_spillNumber);
   tree->SetBranchAddress("evtNumber", &evtNumber, &b_evtNumber);
   tree->SetBranchAddress("ADCvalues", &ADCvalues, &b_ADCvalues);
   tree->SetBranchAddress("digi_max_amplitude", &digi_max_amplitude, &b_digi_max_amplitude); 
   tree->SetBranchAddress("digi_charge_integrated", &digi_charge_integrated, &b_digi_charge_integrated);
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

   std::vector<float> cef3( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3", &cef3 );

   std::vector<float> bgo( BGO_CHANNELS, -1. );
   outTree->Branch( "bgo", &bgo );


   std::vector<float> cef3_corr( CEF3_CHANNELS, -1. );
   outTree->Branch( "cef3_corr", &cef3_corr );

   std::vector<float> bgo_corr( BGO_CHANNELS, -1. );
   outTree->Branch( "bgo_corr", &bgo_corr );

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
   int nClusters_hodoY1;
   outTree->Branch( "nClusters_hodoY1", &nClusters_hodoY1, "nClusters_hodoY1/I" );
   int nFibres_hodoY1[HODOY1_CHANNELS];
   outTree->Branch( "nFibres_hodoY1", nFibres_hodoY1, "nFibres_hodoY1[nClusters_hodoY1]/I" );
   float pos_hodoY1[HODOY1_CHANNELS];
   outTree->Branch( "pos_hodoY1", pos_hodoY1, "pos_hodoY1[nClusters_hodoY1]/F" );

   int nClusters_hodoX2;
   outTree->Branch( "nClusters_hodoX2", &nClusters_hodoX2, "nClusters_hodoX2/I" );
   int nFibres_hodoX2[HODOX2_CHANNELS];
   outTree->Branch( "nFibres_hodoX2", nFibres_hodoX2, "nFibres_hodoX2[nClusters_hodoX2]/I" );
   float pos_hodoX2[HODOX2_CHANNELS];
   outTree->Branch( "pos_hodoX2", pos_hodoX2, "pos_hodoX2[nClusters_hodoX2]/F" );
   int nClusters_hodoY2;
   outTree->Branch( "nClusters_hodoY2", &nClusters_hodoY2, "nClusters_hodoY2/I" );
   int nFibres_hodoY2[HODOY2_CHANNELS];
   outTree->Branch( "nFibres_hodoY2", nFibres_hodoY2, "nFibres_hodoY2[nClusters_hodoY2]/I" );
   float pos_hodoY2[HODOY2_CHANNELS];
   outTree->Branch( "pos_hodoY2", pos_hodoY2, "pos_hodoY2[nClusters_hodoY2]/F" );

   int nClusters_hodoSmallX;
   outTree->Branch( "nClusters_hodoSmallX", &nClusters_hodoSmallX, "nClusters_hodoSmallX/I" );
   int nFibres_hodoSmallX[HODOSMALLX_CHANNELS];
   outTree->Branch( "nFibres_hodoSmallX", nFibres_hodoSmallX, "nFibres_hodoSmallX[nClusters_hodoSmallX]/I" );
   float pos_hodoSmallX[HODOSMALLX_CHANNELS];
   outTree->Branch( "pos_hodoSmallX", pos_hodoSmallX, "pos_hodoSmallX[nClusters_hodoSmallX]/F" );
   int nClusters_hodoSmallY;
   outTree->Branch( "nClusters_hodoSmallY", &nClusters_hodoSmallY, "nClusters_hodoSmallY/I" );
   int nFibres_hodoSmallY[HODOSMALLY_CHANNELS];
   outTree->Branch( "nFibres_hodoSmallY", nFibres_hodoSmallY, "nFibres_hodoSmallY[nClusters_hodoSmallY]/I" );
   float pos_hodoSmallY[HODOSMALLY_CHANNELS];
   outTree->Branch( "pos_hodoSmallY", pos_hodoSmallY, "pos_hodoSmallY[nClusters_hodoSmallY]/F" );


   float wc_x;
   outTree->Branch( "wc_x", &wc_x, "wc_x");
   float wc_y;
   outTree->Branch( "wc_y", &wc_y, "wc_y");



   int nentries = tree->GetEntries();
   RunHelper::getBeamPosition( runName, xBeam, yBeam );

   for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

     tree->GetEntry( iEntry );

     if( iEntry %  10000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;

     //     assignValues( cef3, *digi_max_amplitude, CEF3_START_CHANNEL );
     assignValues( cef3, *digi_charge_integrated, CEF3_START_CHANNEL );  

     assignValues( bgo, *ADCvalues, BGO_ADC_START_CHANNEL );

     assignValues( cef3_corr, *digi_charge_integrated, CEF3_START_CHANNEL );  

     assignValues( bgo_corr, *ADCvalues, BGO_ADC_START_CHANNEL );

     bgoCalib.applyCalibration(bgo_corr);
     cef3Calib.applyCalibration(cef3_corr);

     std::vector<float> hodoX1_values(HODOX1_CHANNELS, -1.);
     std::vector<float> hodoY1_values(HODOY1_CHANNELS, -1.);
     assignValues( hodoX1_values, *ADCvalues, HODOX1_ADC_START_CHANNEL );
     assignValues( hodoY1_values, *ADCvalues, HODOY1_ADC_START_CHANNEL );

     std::vector<float> hodoX2_values(HODOX2_CHANNELS, -1.);
     std::vector<float> hodoY2_values(HODOY2_CHANNELS, -1.);
     assignValues( hodoX2_values, *ADCvalues, HODOX2_ADC_START_CHANNEL );
     assignValues( hodoY2_values, *ADCvalues, HODOY2_ADC_START_CHANNEL );

     std::vector<float> hodoSmallX_values(HODOSMALLX_CHANNELS, -1.);
     std::vector<float> hodoSmallY_values(HODOSMALLY_CHANNELS, -1.);
     assignValues( hodoSmallX_values, *ADCvalues, HODOSMALLX_ADC_START_CHANNEL );
     assignValues( hodoSmallY_values, *ADCvalues, HODOSMALLY_ADC_START_CHANNEL );

     // hodo cluster reconstruction
     int clusterMaxFibres = 4;
     doHodoReconstruction( hodoX1_values    , nClusters_hodoX1    , nFibres_hodoX1    , pos_hodoX1    , 0.5, clusterMaxFibres );
     doHodoReconstruction( hodoY1_values    , nClusters_hodoY1    , nFibres_hodoY1    , pos_hodoY1    , 0.5, clusterMaxFibres );
     doHodoReconstruction( hodoX2_values    , nClusters_hodoX2    , nFibres_hodoX2    , pos_hodoX2    , 0.5, clusterMaxFibres );
     doHodoReconstruction( hodoY2_values    , nClusters_hodoY2    , nFibres_hodoY2    , pos_hodoY2    , 0.5, clusterMaxFibres );
     doHodoReconstruction( hodoSmallX_values, nClusters_hodoSmallX, nFibres_hodoSmallX, pos_hodoSmallX, 1.0, clusterMaxFibres );
     doHodoReconstruction( hodoSmallY_values, nClusters_hodoSmallY, nFibres_hodoSmallY, pos_hodoSmallY, 1.0, clusterMaxFibres );


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





std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax ) {

  std::vector<HodoCluster*> clusters;

  HodoCluster* currentCluster = new HodoCluster( hodo.size(), fibreWidth );

  for( unsigned i=0; i<hodo.size(); ++i ) {

    if( hodo[i] > 0.) { // hit

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




void doHodoReconstruction( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres ) {

  std::vector<HodoCluster*> clusters = getHodoClusters( values, fibreWidth, clusterMaxFibres );

  nClusters = clusters.size();
  for( unsigned i=0; i<clusters.size(); ++i ) {
    nFibres[i] = clusters[i]->getSize();
    pos[i] = clusters[i]->getPosition();
  }

}
