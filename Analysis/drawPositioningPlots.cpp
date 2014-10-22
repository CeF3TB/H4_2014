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

  std::string fileName = "analysisTrees_" + tag + "/Reco_" + runName + ".root";
  TFile* file = TFile::Open(fileName.c_str());
  if( file==0 ) {
    std::cout << "ERROR! Din't find file " << fileName << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  TTree* tree = (TTree*)file->Get("recoTree");


  std::string outputdir = "PositionPlots_" + tag + "/" + runName + "/";
  system( Form("mkdir -p %s", outputdir.c_str()) );


  drawSingleVariable( outputdir, tree, "wc_x_vs_hodoX1", "wc_x - pos_hodoX1[0]", "nClusters_hodoX1==1", 50, -20., 20., "X( WC - HodoX1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_y_vs_hodoY1", "wc_y - pos_hodoY1[0]", "nClusters_hodoY1==1", 50, -20., 20., "Y( WC - HodoY1 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_x_vs_hodoX2", "wc_x - pos_hodoX2[0]", "nClusters_hodoX2==1", 50, -20., 20., "X( WC - HodoX2 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_y_vs_hodoY2", "wc_y - pos_hodoY2[0]", "nClusters_hodoY2==1", 50, -20., 20., "Y( WC - HodoY2 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_x_vs_hodoSmallX", "wc_x - pos_hodoSmallX[0]", "nClusters_hodoSmallX==1", 50, -20., 20., "X( WC - HodoSmallX ) [mm]"  );
  drawSingleVariable( outputdir, tree, "wc_y_vs_hodoSmallY", "wc_y - pos_hodoSmallY[0]", "nClusters_hodoSmallY==1", 50, -20., 20., "Y( WC - HodoSmallY ) [mm]"  );

  drawSingleVariable( outputdir, tree, "hodoX1_vs_hodoX2", "pos_hodoX1[0] - pos_hodoX2[0]", "nClusters_hodoX1==1 && nClusters_hodoX2==1", 50, -20., 20., "X( HodoX1 - HodoX2 ) [mm]"  );
  drawSingleVariable( outputdir, tree, "hodoY1_vs_hodoY2", "pos_hodoY1[0] - pos_hodoY2[0]", "nClusters_hodoY1==1 && nClusters_hodoY2==1", 50, -20., 20., "Y( HodoY1 - HodoY2 ) [mm]"  );




  return 0;

}



void drawSingleVariable( const std::string& outputdir, TTree* tree, const std::string& savename, const std::string& treeVar, const std::string& selection, int nbins, float xMin, float xMax, std::string axisName ) {

  if( axisName=="" ) axisName = treeVar;

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH1F* h1 = new TH1F( savename.c_str(), "", nbins, xMin, xMax );
  tree->Project( savename.c_str(), treeVar.c_str(), selection.c_str() );

  float yMax = 1.2*h1->GetMaximum();

  TH2D* h2_axes = new TH2D( "axes", "", nbins, xMin, xMax, 10, 0., yMax );
  h2_axes->SetXTitle( axisName.c_str() );
  h2_axes->SetYTitle( "Entries" );
  h2_axes->Draw();

  h1->SetFillColor(kOrange);

  h1->Draw("histo same");

  TPaveText* labelTop = DrawTools::getLabelTop();
  labelTop->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/%s.eps", outputdir.c_str(), savename.c_str()) );
  c1->SaveAs( Form("%s/%s.png", outputdir.c_str(), savename.c_str()) );

  delete c1;
  delete h2_axes;
  delete h1;

}
