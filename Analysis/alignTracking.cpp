#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"

#include "HodoCluster.h"
#include "channelInfo.h"



void assignValues( std::vector<float> &target, std::vector<float> source, unsigned int startPos );

void doHodoReconstruction( std::vector<float> values, int &nClusters, int *nFibres, float *pos, float fibreWidth, int clusterMaxFibres );
std::vector<HodoCluster*> getHodoClusters( std::vector<float> hodo, float fibreWidth, int nClusterMax );

float fitAndDraw( const std::string& outputdir, TH1F* h1 );


int main() {

  // this is the dir in which the constants will be saved:
  std::string constDirName = "AlignmentConstants";
  system(Form("mkdir -p %s", constDirName.c_str()));


  TFile* file = TFile::Open("data/run_273.root");
  TTree* tree = (TTree*)file->Get("outputTree");


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


  int nBins = 50;
  float xMin = -20.;
  float xMax =  20.;

  TH1F* h1_wc_y_low = new TH1F("wc_y_low", "", nBins, xMin, xMax);
  TH1F* h1_wc_y_hi  = new TH1F("wc_y_hi" , "", nBins, xMin, xMax);
  TH1F* h1_wc_x_low = new TH1F("wc_x_low", "", nBins, xMin, xMax);
  TH1F* h1_wc_x_hi  = new TH1F("wc_x_hi" , "", nBins, xMin, xMax);

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
     assignValues( hodoSmallX_values, *digi_max_amplitude, HODOSMALLX_ADC_START_CHANNEL );
     assignValues( hodoSmallY_values, *digi_max_amplitude, HODOSMALLY_ADC_START_CHANNEL );

     // hodo cluster reconstruction
     int clusterMaxFibres = 4;
     doHodoReconstruction( hodoX1_values    , nClusters_hodoX1    , nFibres_hodoX1    , pos_hodoX1    , 0.5, clusterMaxFibres );
     doHodoReconstruction( hodoY1_values    , nClusters_hodoY1    , nFibres_hodoY1    , pos_hodoY1    , 0.5, clusterMaxFibres );
     doHodoReconstruction( hodoX2_values    , nClusters_hodoX2    , nFibres_hodoX2    , pos_hodoX2    , 0.5, clusterMaxFibres );
     doHodoReconstruction( hodoY2_values    , nClusters_hodoY2    , nFibres_hodoY2    , pos_hodoY2    , 0.5, clusterMaxFibres );


     wc_x = ADCvalues->at(WC_X_ADC_START_CHANNEL);
     wc_y = ADCvalues->at(WC_Y_ADC_START_CHANNEL);
     if( runNumber>=170 ) wc_y = -wc_y; // temporary fix

     float hodoSmallY_low = digi_max_amplitude->at(6);
     float hodoSmallY_hi  = digi_max_amplitude->at(7);
     float hodoSmallX_low = digi_max_amplitude->at(4);
     float hodoSmallX_hi  = digi_max_amplitude->at(5);

     // y low
     if( hodoSmallY_low>60. ) {
       h1_wc_y_low->Fill( wc_y );
       for( unsigned i=0; i<nClusters_hodoY1; ++i ) {
         if( nFibres_hodoY1[i]!=1 ) {  // cut away some noise
           h1_hodoY1_low->Fill( pos_hodoY1[i] );
         }
       } // for Y1
       for( unsigned i=0; i<nClusters_hodoY2; ++i ) {
         if( nFibres_hodoY2[i]!=1 ) {  // cut away some noise
           h1_hodoY2_low->Fill( pos_hodoY2[i] );
         }
       } // for Y2
     }

     // y hi
     if( hodoSmallY_hi>60. ) {
       h1_wc_y_hi->Fill( wc_y );
       for( unsigned i=0; i<nClusters_hodoY1; ++i ) {
         if( nFibres_hodoY1[i]!=1 ) {  // cut away some noise
           h1_hodoY1_hi->Fill( pos_hodoY1[i] );
         }
       } // for Y1
       for( unsigned i=0; i<nClusters_hodoY2; ++i ) {
         if( nFibres_hodoY2[i]!=1 ) {  // cut away some noise
           h1_hodoY2_hi->Fill( pos_hodoY2[i] );
         }
       } // for Y2
     }

     // x low
     if( hodoSmallX_low>60. ) {
       h1_wc_x_low->Fill( wc_x );
       for( unsigned i=0; i<nClusters_hodoX1; ++i ) {
         if( nFibres_hodoX1[i]!=1 ) {  // cut away some noise
           h1_hodoX1_low->Fill( pos_hodoX1[i] );
         }
       } // for X1
       for( unsigned i=0; i<nClusters_hodoX2; ++i ) {
         if( nFibres_hodoX2[i]!=1 ) {  // cut away some noise
           h1_hodoX2_low->Fill( pos_hodoX2[i] );
         }
       } // for X2
     }

     // x hi
     if( hodoSmallX_hi>60. ) {
       h1_wc_x_hi->Fill( wc_x );
       for( unsigned i=0; i<nClusters_hodoX1; ++i ) {
         if( nFibres_hodoX1[i]!=1 ) {  // cut away some noise
           h1_hodoX1_hi->Fill( pos_hodoX1[i] );
         }
       } // for X1
       for( unsigned i=0; i<nClusters_hodoX2; ++i ) {
         if( nFibres_hodoX2[i]!=1 ) {  // cut away some noise
           h1_hodoX2_hi->Fill( pos_hodoX2[i] );
         }
       } // for X2
     }

  } // for entries


  std::string outputdir = "AlignmentPlots";
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

  float offset_wc_y  = 0.5*(offset_wc_y_low+offset_wc_y_hi);
  float offset_wc_x  = 0.5*(offset_wc_x_low+offset_wc_x_hi);
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

  
  std::string ofsName = constDirName + "/constants_new.txt";
  ofstream ofs(ofsName.c_str());
  ofs << "wc_y "   << -offset_wc_y << std::endl;
  ofs << "wc_x "   << -offset_wc_x << std::endl;
  ofs << "hodoY1 " << -offset_hodoY1 << std::endl;
  ofs << "hodoX1 " << -offset_hodoX1 << std::endl;
  ofs << "hodoY2 " << -offset_hodoY2 << std::endl;
  ofs << "hodoX2 " << -offset_hodoX2 << std::endl;
  ofs.close();

  std::cout << "-> Saved constants in: " << ofsName << std::endl;

  return 0;

}


float fitAndDraw( const std::string& outputdir, TH1F* h1 ) {

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  float maxPos = h1->GetBinCenter(h1->GetMaximumBin());

  TF1* f1 = new TF1(Form("f1_%s", h1->GetName()), "gaus", -10., 10.);
  
  f1->SetParameter( 1, maxPos );
  h1->Fit(f1, "QRN");

  for( unsigned i=0; i<4; ++i ) {

    float m = f1->GetParameter(1);
    float s = f1->GetParameter(2);
    float nSigmas = 1.5;
    f1->SetRange( m-nSigmas*s, m+nSigmas*s );
    if( i==3 ) 
      h1->Fit(f1, "RQ");
    else
      h1->Fit(f1, "RNQ");

  }

  f1->SetLineColor(kRed);

  h1->Draw();

  c1->SaveAs(Form("%s/%s.eps", outputdir.c_str(), h1->GetName()));
  c1->SaveAs(Form("%s/%s.png", outputdir.c_str(), h1->GetName()));

  delete c1;

  return f1->GetParameter(1);

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