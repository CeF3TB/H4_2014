#include <iostream>

#include "DrawTools.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"




void drawSingleVariable( const std::string& outputdir, TTree* tree, const std::string& savename, const std::string& treeVar, const std::string& selection, int nbins, float xMin, float xMax, std::string axisName="" );


int main( int argc, char* argv[] ) {


  DrawTools::setStyle();

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

  std::string fileName = "analysisTrees_" + tag + "/Reco_" + runName + ".root";
  TFile* file = TFile::Open(fileName.c_str());
  if( file==0 ) {
    std::cout << "ERROR! Din't find file " << fileName << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  TTree* tree = (TTree*)file->Get("recoTree");


  std::string outputdir = "PositionPlots_" + tag + "/" + runName;
  system( Form("mkdir -p %s", outputdir.c_str()) );


  // uncorrected:

  drawSingleVariable( outputdir, tree, "wc_x", "wc_x", "", 50, -20., 20., "X( WC ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_y", "wc_y", "", 50, -20., 20., "Y( WC ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoX1", "cluster_pos_hodoX1", "  cluster_pos_hodoX1>-300  && nFibres_hodoSmallX >0  ", 50, -20., 20., "X( hodoX1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoY1", "cluster_pos_hodoY1", " nFibres_hodoSmallY >0", 50, -20., 20., "Y( hodoY1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoX2", "cluster_pos_hodoX2", " nFibres_hodoSmallX >0", 50, -20., 20., "X( hodoX2 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoY2", "cluster_pos_hodoY2", "nFibres_hodoSmallY >0", 50, -20., 20., "Y( hodoY2 ) [mm]"  );

  drawSingleVariable( outputdir, tree, "wc_x_vs_hodoX1", "wc_x - cluster_pos_hodoX1", " nFibres_hodoSmallX >0  && wc_x > -300 && nTDCHits[0]>0 && nTDCHits[1]>0 && nTDCHits[0]<3 && nTDCHits[1]<3 && (nTDCHits[0]+nTDCHits[1])<4 ", 50, -10., 10., "X( WC - HodoX1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_y_vs_hodoY1", "wc_y - cluster_pos_hodoY1", "nFibres_hodoSmallY >0", 50, -10., 10., "Y( WC - HodoY1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_x_vs_hodoX2", "wc_x - cluster_pos_hodoX2", "nFibres_hodoSmallX >0  && wc_x > -300 && nTDCHits[0]>0 && nTDCHits[1]>0 && nTDCHits[0]<3 && nTDCHits[1]<3 && (nTDCHits[0]+nTDCHits[1])<4", 50, -10., 10., "X( WC - HodoX2 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_y_vs_hodoY2", "wc_y - cluster_pos_hodoY2", "nFibres_hodoSmallY >0", 50, -10., 10., "Y( WC - HodoY2 ) [mm]"  );

  drawSingleVariable( outputdir, tree, "hodoX1_vs_hodoX2", "cluster_pos_hodoX1 - cluster_pos_hodoX2", "nFibres_hodoSmallX >0 && cluster_pos_hodoX1>-300 && cluster_pos_hodoX2>-300 ", 50, -10., 10., "X( HodoX1 - HodoX2 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoY1_vs_hodoY2", "cluster_pos_hodoY1 - cluster_pos_hodoY2", " nFibres_hodoSmallY >0 && cluster_pos_hodoY1 >-300  && cluster_pos_hodoY2 >-300", 50, -10., 10., "Y( HodoY1 - HodoY2 ) [mm]"  );



  // corrected:

  drawSingleVariable( outputdir, tree, "wc_x_corr", "wc_x_corr", "", 50, -20., 20., "X( WC ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_y_corr", "wc_y_corr", "", 50, -20., 20., "Y( WC ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoX1_corr", "cluster_pos_corr_hodoX1", " cluster_pos_hodoX1 >-300 && nFibres_hodoSmallX >0", 50, -20., 20., "X( hodoX1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoY1_corr", "cluster_pos_corr_hodoY1", " cluster_pos_hodoY1 >-300 && nFibres_hodoSmallY >0", 50, -20., 20., "Y( hodoY1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoX2_corr", "cluster_pos_corr_hodoX2", " cluster_pos_hodoX2 >-300 &&nFibres_hodoSmallX >0", 50, -20., 20., "X( hodoX2 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoY2_corr", "cluster_pos_corr_hodoY2", "cluster_pos_hodoY2 >-300 && nFibres_hodoSmallY >0", 50, -20., 20., "Y( hodoY2 ) [mm]"  );



  drawSingleVariable( outputdir, tree, "wc_x_vs_hodoX1_corr", "wc_x_corr - cluster_pos_corr_hodoX1", " nFibres_hodoSmallX >0 && cluster_pos_hodoX1>-300 &&wc_x_corr>-400 && nTDCHits[0]>0 && nTDCHits[1]>0 && nTDCHits[0]<3 && nTDCHits[1]<3 && (nTDCHits[0]+nTDCHits[1])<4", 50, -10., 10., "X( WC - HodoX1 ) [mm]"  );

  drawSingleVariable( outputdir, tree, "wc_y_vs_hodoY1_corr", "wc_y_corr -  cluster_pos_corr_hodoY1", "nFibres_hodoSmallY >0  && cluster_pos_hodoY1>-300 &&wc_y_corr>-400 && nTDCHits[2]>0 && nTDCHits[3]>0 && nTDCHits[2]<3 && nTDCHits[3]<3 && (nTDCHits[2]+nTDCHits[3])<4 ", 50, -10., 10., "Y( WC - HodoY1 ) [mm]"  );

  drawSingleVariable( outputdir, tree, "wc_x_vs_hodoX2_corr", "wc_x_corr - cluster_pos_corr_hodoX2", "nFibres_hodoSmallX >0 && cluster_pos_corr_hodoX2>-300 &&wc_x_corr>-400 && nTDCHits[0]>0 && nTDCHits[1]>0 && nTDCHits[0]<3 && nTDCHits[1]<3 && (nTDCHits[0]+nTDCHits[1])<4", 50, -10., 10., "X( WC - HodoX2 ) [mm]"  );
 
 drawSingleVariable( outputdir, tree, "wc_y_vs_hodoY2_corr", "wc_y_corr -cluster_pos_corr_hodoY2", "nFibres_hodoSmallY >0 && cluster_pos_hodoY2>-300 &&wc_x_corr>-400 && nTDCHits[0]>0 && nTDCHits[1]>0 && nTDCHits[0]<3 && nTDCHits[1]<3 && (nTDCHits[0]+nTDCHits[1])<4", 50, -10., 10., "Y( WC - HodoY2 ) [mm]"  );

  drawSingleVariable( outputdir, tree, "hodoX1_vs_hodoX2_corr", "cluster_pos_corr_hodoX1 - cluster_pos_corr_hodoX2", "nFibres_hodoSmallX >0 && cluster_pos_corr_hodoX1>-300 && cluster_pos_corr_hodoX2>-300", 50, -10., 10., "X( HodoX1 - HodoX2 ) [mm]"  );
 
 drawSingleVariable( outputdir, tree, "hodoY1_vs_hodoY2_corr", "cluster_pos_corr_hodoY1 - cluster_pos_corr_hodoY2", "nFibres_hodoSmallY >0  && cluster_pos_hodoY2>-300 && cluster_pos_hodoY1>-300", 50, -10., 10., "Y( HodoY1 - HodoY2 ) [mm]"  );



  /*

  // uncorrected:

  drawSingleVariable( outputdir, tree, "wc_x", "wc_x", "", 50, -20., 20., "X( WC ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_y", "wc_y", "", 50, -20., 20., "Y( WC ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoX1", "pos_hodoX1[0]", "nClusters_hodoX1==1", 50, -20., 20., "X( hodoX1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoY1", "pos_hodoY1[0]", "nClusters_hodoY1==1", 50, -20., 20., "Y( hodoY1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoX2", "pos_hodoX2[0]", "nClusters_hodoX2==1", 50, -20., 20., "X( hodoX2 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoY2", "pos_hodoY2[0]", "nClusters_hodoY2==1", 50, -20., 20., "Y( hodoY2 ) [mm]"  );

  drawSingleVariable( outputdir, tree, "wc_x_vs_hodoX1", "wc_x - pos_hodoX1[0]", "nClusters_hodoX1==1", 50, -10., 10., "X( WC - HodoX1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_y_vs_hodoY1", "wc_y - pos_hodoY1[0]", "nClusters_hodoY1==1", 50, -10., 10., "Y( WC - HodoY1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_x_vs_hodoX2", "wc_x - pos_hodoX2[0]", "nClusters_hodoX2==1", 50, -10., 10., "X( WC - HodoX2 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_y_vs_hodoY2", "wc_y - pos_hodoY2[0]", "nClusters_hodoY2==1", 50, -10., 10., "Y( WC - HodoY2 ) [mm]"  );

  drawSingleVariable( outputdir, tree, "hodoX1_vs_hodoX2", "pos_hodoX1[0] - pos_hodoX2[0]", "nClusters_hodoX1==1 && nClusters_hodoX2==1", 50, -10., 10., "X( HodoX1 - HodoX2 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoY1_vs_hodoY2", "pos_hodoY1[0] - pos_hodoY2[0]", "nClusters_hodoY1==1 && nClusters_hodoY2==1", 50, -10., 10., "Y( HodoY1 - HodoY2 ) [mm]"  );


  // corrected:

  drawSingleVariable( outputdir, tree, "wc_x_corr", "wc_x_corr", "", 50, -20., 20., "X( WC ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_y_corr", "wc_y_corr", "", 50, -20., 20., "Y( WC ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoX1_corr", "pos_corr_hodoX1[0]", "nClusters_hodoX1==1", 50, -20., 20., "X( hodoX1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoY1_corr", "pos_corr_hodoY1[0]", "nClusters_hodoY1==1", 50, -20., 20., "Y( hodoY1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoX2_corr", "pos_corr_hodoX2[0]", "nClusters_hodoX2==1", 50, -20., 20., "X( hodoX2 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoY2_corr", "pos_corr_hodoY2[0]", "nClusters_hodoY2==1", 50, -20., 20., "Y( hodoY2 ) [mm]"  );

  drawSingleVariable( outputdir, tree, "wc_x_vs_hodoX1_corr", "wc_x_corr - pos_corr_hodoX1[0]", "nClusters_hodoX1==1", 50, -10., 10., "X( WC - HodoX1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_y_vs_hodoY1_corr", "wc_y_corr - pos_corr_hodoY1[0]", "nClusters_hodoY1==1", 50, -10., 10., "Y( WC - HodoY1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_x_vs_hodoX2_corr", "wc_x_corr - pos_corr_hodoX2[0]", "nClusters_hodoX2==1", 50, -10., 10., "X( WC - HodoX2 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_y_vs_hodoY2_corr", "wc_y_corr - pos_corr_hodoY2[0]", "nClusters_hodoY2==1", 50, -10., 10., "Y( WC - HodoY2 ) [mm]"  );

  drawSingleVariable( outputdir, tree, "hodoX1_vs_hodoX2_corr", "pos_corr_hodoX1[0] - pos_corr_hodoX2[0]", "nClusters_hodoX1==1 && nClusters_hodoX2==1", 50, -10., 10., "X( HodoX1 - HodoX2 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoY1_vs_hodoY2_corr", "pos_corr_hodoY1[0] - pos_corr_hodoY2[0]", "nClusters_hodoY1==1 && nClusters_hodoY2==1", 50, -10., 10., "Y( HodoY1 - HodoY2 ) [mm]"  );


  */



  return 0;

}



void drawSingleVariable( const std::string& outputdir, TTree* tree, const std::string& savename, const std::string& treeVar, const std::string& selection, int nbins, float xMin, float xMax, std::string axisName ) {

  if( axisName=="" ) axisName = treeVar;

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH1F* h1 = new TH1F( savename.c_str(), "", nbins, xMin, xMax );
  tree->Project( savename.c_str(), treeVar.c_str(), selection.c_str() );

  float yMax = 1.1*h1->GetMaximum();

  TH2D* h2_axes = new TH2D( "axes", "", nbins, xMin, xMax, 10, 0., yMax );
  h2_axes->SetXTitle( axisName.c_str() );
  h2_axes->SetYTitle( "Entries" );
  h2_axes->Draw("");




  float returnConst=0.;

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
    }

    TLine* lineThresh = new TLine( xMin, thresh, xMax, thresh );
    lineThresh->SetLineStyle(2);
    //    lineThresh->Draw("same");

    returnConst = total/denom;
    h1_usedBins->SetFillColor(kOrange);
    h1_usedBins->Draw("same");


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








  h1->SetFillColor(kOrange);

  h1->Draw("histo same");

  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");

  TLine* lineOne = new TLine( 0., 0., 0., yMax );
  lineOne->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/%s.eps", outputdir.c_str(), savename.c_str()) );
  c1->SaveAs( Form("%s/%s.png", outputdir.c_str(), savename.c_str()) );

  delete c1;
  delete h2_axes;
  delete h1;

}
