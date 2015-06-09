#include <iostream>
#include <string>
#include <cmath>


#include "DrawTools.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"





int main( int argc, char* argv[] ) {


  DrawTools::setStyle();


  gStyle->SetPadLeftMargin(0.1);

  std::string fileName = "analysisTrees_dev/Reco_corr3.root";
  TFile* file = TFile::Open(fileName.c_str());
  if( file==0 ) {
    std::cout << "ERROR! Din't find file " << fileName << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  TTree* tree = (TTree*)file->Get("recoTree");


  std::string outputdir = "AlignmentCorrvsRun";
  system( Form("mkdir -p %s", outputdir.c_str()) );


  TCanvas *canny = new TCanvas("canny", "", 1000, 600);
  canny->cd();





  std::string colors[7] ={"kBlack","kCyan","kBlue","kOrange","kRed","kGreen","kMagenta" };
  std::string energies[7] = {"beamEnergy==10","beamEnergy==15","beamEnergy==20","beamEnergy==50","beamEnergy==100","beamEnergy==150","beamEnergy==200"};
  int nentries=  tree->GetEntries();


  std::string selection = "nClusters_hodoX2==1 && abs(wc_x-cluster_pos_hodoX2)<100 && run>310  && cluster_pos_hodoX2[0]>-300";
  //  std::string selection = " abs(wc_x-pos_hodoX2[0])<100 && run>270 && cluster_pos_hodoX2>-300";

  TH2D* h1 = new TH2D( "h1", "h1", 400, 310, 600, 10, -4,4. );
  // tree->Project( "h1" ,"(wc_x - cluster_pos_hodoX2[0]) : run",  selection.c_str(), "profile");

  tree->Draw( "(wc_x - cluster_pos_hodoX2[0]) : run  >> h1",  selection.c_str(), "profile");


   std::string selection2 = "nClusters_hodoX2==1 && nFibres_hodoX2[0]>0 && abs(wc_x-cluster_pos_hodoX2[0])<100 && run>310";

  TH2D* h2 = new TH2D( "h2", "h2", 400, 310, 600, 10, -4,4. );
  //tree->Project( "h2" ,"(wc_x_corr - cluster_pos_corr_hodoX2[0]) : run",  selection2.c_str(),"profile");
  tree->Draw( "(wc_x_corr - cluster_pos_corr_hodoX2[0]) : run  >> h2",  selection2.c_str(), "profile");


  TH2D* h2_axes = new TH2D( "axes", "", 400, 310, 600, 10, -4, 2. );
  h2_axes->SetXTitle( "Run Number" );
  h2_axes->SetYTitle( "x(wire) - x(Hodo2) [mm]" );
  h2_axes->Draw("");
 h2_axes->GetYaxis()->SetTitleOffset(1.0);


  h1 = (TH2D*)gDirectory->Get("h1"); 
  h2 = (TH2D*)gDirectory->Get("h2");


  h1->SetLineWidth(2);
  h1->SetMarkerStyle(21);
  h1->SetMarkerColor(kBlue);
  h1->Draw("same");

 
  h2->SetLineColor(kRed); 
  h2->SetMarkerStyle(20);
  h2->SetMarkerColor(kRed);
  //   h2->SetLineWidth(2);
  h2->Draw("same");


  TPaveText* labelTop = DrawTools::getLabelTop("Electron Beam"); 
labelTop->Draw("same");
  TLine* lineZero = new TLine( 310 , 0., 600 , 0 );
  lineZero->SetLineColor(kBlack);  lineZero->SetLineWidth(2);  lineZero->Draw("same");
  TLine* linePlus = new TLine( 310 , 0.5, 600 , 0.5 );
  linePlus->SetLineColor(kBlack);  linePlus->SetLineWidth(2);
  linePlus->SetLineStyle(2);  linePlus->Draw("same");
  TLine* lineMinus = new TLine( 310 , -0.5, 600 , -0.5 );
  lineMinus->SetLineColor(kBlack); lineMinus->SetLineWidth(2);
  lineMinus->SetLineStyle(2);  lineMinus->Draw("same");

  //  TLegend* legend = new TLegend( 0.2, 0.75, 0.45, 0.9 );
  TLegend* legend = new TLegend( 0.7, 0.92-2.*0.06, 0.9, 0.92 );
  legend->SetTextSize(0.035);  legend->SetFillColor(0);
  legend->AddEntry( h1 , "Before Alignment", "P" );  
  legend->AddEntry( h2 , "After Alignment", "P" );
  legend->Draw("same");

  gPad->RedrawAxis();

  canny->SaveAs( Form("%s/XWC-HodoX2_small.eps", outputdir.c_str() ) );
  canny->SaveAs( Form("%s/XWC-HodoX2_small.png", outputdir.c_str() ) );
  canny->SaveAs( Form("%s/XWC-HodoX2_small.pdf", outputdir.c_str() ) );

 

  delete h1;
  delete h2;
  delete h2_axes;

  return 0;

}


