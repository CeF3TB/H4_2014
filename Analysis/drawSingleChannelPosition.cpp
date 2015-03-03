#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TProfile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TChain.h"
#include "TProfile2D.h"


#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooPlot.h"

#include "RooGlobalFunc.h"

#include "DrawTools.h"
using namespace RooFit;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
/// Fill Histos with E_i/Sum(E_i) and see how that behaves vs r  ///////////////
////////////////////////////////////////////////////////////////////////////////




float getRatioError( float num, float denom, float numErr, float denomErr );

float findRadius(TF1* f,int order, double fraction);

int getMaxFrac(double *frac);
int getSecondFrac(double *frac);


int main( int argc, char* argv[] ) {


  std::string tag="V01";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }
  std::cout << "-> Using tag: " << tag << std::endl;

  std::string outputdir = "RespVsR/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  DrawTools::setStyle();
//  TGaxis::SetMaxDigits(3);

  TGaxis::SetMaxDigits(3);
  //  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.12);

  

  //DATA RUNS// Let's start with the runs at 50 GeV and 950 V
    // TGraphErrors* gr_resp_vs_pos = new TGraphErrors(0);

    std::string cut = "";
    
    
    TCanvas* c1 = new TCanvas( "c1", "",0,0, 1200, 1200 );
    c1->cd();
    TH2D* h2_axes = new TH2D( "axes", "", 24 , -12, 12,10,1000, 160000  );
    //  h2_axes->SetXTitle("Position X [mm]");
    //    h2_axes->SetYTitle("Position Y [mm]");
    // h2_axes->Draw("");

    

  TLine* lineLower = new TLine( -12+2.1, -12., 12-2.1, -12 );
  lineLower->SetLineColor(kOrange);
  lineLower->SetLineWidth(3); 

  TLine* lineRight = new TLine( 12, 12.-2.1, 12, -12+2.1 );
  lineRight->SetLineColor(kOrange);
  lineRight->SetLineWidth(3); 

  TLine* lineLeft = new TLine( -12, 12.-2.1, -12, -12 +2.1);
  lineLeft->SetLineColor(kOrange);
  lineLeft->SetLineWidth(3);  
 
  TLine* lineUpper = new TLine( -12+2.1, 12., 12 -2.1, 12 );
  lineUpper->SetLineColor(kOrange);
  lineUpper->SetLineWidth(3);  

  TLine* chamferDR = new TLine( 12-2.1, -12., 12 , -12+2.1 );
  chamferDR->SetLineColor(kOrange);
  chamferDR->SetLineWidth(3); 

  TLine* chamferUR = new TLine( 12-2.1, 12., 12 , 12-2.1 );
  chamferUR->SetLineColor(kOrange);
  chamferUR->SetLineWidth(3); 

  TLine* chamferDL = new TLine( -12+2.1, -12., -12 , -12+2.1 );
  chamferDL->SetLineColor(kOrange);
  chamferDL->SetLineWidth(3); 

  TLine* chamferUL = new TLine(- 12+2.1, 12., -12 , 12-2.1 );
  chamferUL->SetLineColor(kOrange);
  chamferUL->SetLineWidth(3);
   

  TPaveText* label_top2 = new TPaveText();
  label_top2 = DrawTools::getLabelTop("50 GeV Electron Beam");


  std::map<string,TObject*> outHistos; 

  TProfile* hprofR[4];
  
  TProfile2D* fracProf[4];
  

  // int nBins = 8;
  //  int nBins = 16;
  int nBins = 32;
  //int nBins = 64;
  //  int nBins = 128;

  for(int j =0; j<4; j++){ //le lüp over ze fibres ~(^.^)~
    
    hprofR[j] = new TProfile(Form("hprofR_%d",j),"profile of values vs pos y",100 , 0, 50 , 0, 5000000);
    hprofR[j]->SetXTitle("Position R [mm]");
    hprofR[j]->SetYTitle("Response [ADC]");
    outHistos[hprofR[j]->GetName()]=(TObject*) hprofR[j];




    fracProf[j] = new TProfile2D(Form("fracProf_%d",j) ,"profile of values vs pos x&y", nBins, -16,16, nBins  ,-16,16, 0.0, 20.);
    fracProf[j]->SetXTitle("Position X [mm]");
    fracProf[j]->SetYTitle("Position Y [mm]");

}




  
  TChain ch("recoTree"); // creates a chain to process a Tree called "T"
  for(int i = 274; i<323; i++){
     ch.Add(Form("analysisTrees_chaInt_04_12_%s/Reco_%d.root", tag.c_str(),  i ));
  }

  float xTable;
  float yTable;
  float wc_x;
  float wc_y;
  float wc_x_corr;
  float wc_y_corr;
  float HVCeF3;
  float beamEnergy;

  float pos_2FibClust_hodoX1;
  float pos_2FibClust_hodoX2;
  float pos_2FibClust_hodoY1;
  float pos_2FibClust_hodoY2;
  
  float cluster_pos_corr_hodoX1;
  float cluster_pos_corr_hodoX2;
  float cluster_pos_corr_hodoY1;
  float cluster_pos_corr_hodoY2;
  
  float position_X1;
  float position_X2;
  float position_Y1;
  float position_Y2;
  
  std::vector<int> *nTDCHits;
  std::vector<float>   *cef3_chaInt;
  std::vector<float>   *cef3_maxAmpl;
  std::vector<float>   *cef3_chaInt_corr;
  std::vector<float>   *cef3_maxAmpl_corr;
  
  //Branches
  TBranch *b_wc_x;
  TBranch *b_wc_y;
  TBranch *b_wc_x_corr;
  TBranch *b_wc_y_corr;
  TBranch *b_xTable;
  TBranch *b_yTable;
  TBranch *b_HVCeF3;
  TBranch *b_beamEnergy;
  TBranch *b_cef3_chaInt;
  TBranch *b_cef3_maxAmpl;
  TBranch *b_cef3_chaInt_corr;
  TBranch *b_cef3_maxAmpl_corr;
  TBranch *b_nTDCHits;
  
  TBranch *b_pos_2FibClust_hodoX1;
  TBranch *b_pos_2FibClust_hodoX2;
  TBranch *b_pos_2FibClust_hodoY1;
  TBranch *b_pos_2FibClust_hodoY2;
  
  TBranch *b_cluster_pos_corr_hodoX1;
  TBranch *b_cluster_pos_corr_hodoX2;
  TBranch *b_cluster_pos_corr_hodoY1;
  TBranch *b_cluster_pos_corr_hodoY2;
  
  TBranch *b_position_X1;
  TBranch *b_position_X2;
  TBranch *b_position_Y1;
  TBranch *b_position_Y2;
  
  //Set Object Pointer
  cef3_chaInt = 0;
  cef3_maxAmpl = 0;
  cef3_chaInt_corr = 0;
  cef3_maxAmpl_corr = 0;
  nTDCHits = 0;
  
  //Set branch addresses
  ch.SetBranchAddress("wc_x", &wc_x, &b_wc_x);
  ch.SetBranchAddress("wc_y", &wc_y, &b_wc_y);
  ch.SetBranchAddress("wc_x_corr", &wc_x_corr, &b_wc_x_corr);
  ch.SetBranchAddress("wc_y_corr", &wc_y_corr, &b_wc_y_corr);
  ch.SetBranchAddress("cef3_chaInt", &cef3_chaInt, &b_cef3_chaInt);
  ch.SetBranchAddress("cef3_maxAmpl", &cef3_maxAmpl, &b_cef3_maxAmpl);
  ch.SetBranchAddress("cef3_chaInt_corr", &cef3_chaInt_corr, &b_cef3_chaInt_corr);
  ch.SetBranchAddress("cef3_maxAmpl_corr", &cef3_maxAmpl_corr, &b_cef3_maxAmpl_corr);
  ch.SetBranchAddress("xTable", &xTable, &b_xTable);
  ch.SetBranchAddress("yTable", &yTable, &b_yTable);
  
  ch.SetBranchAddress("beamEnergy", &beamEnergy, &b_beamEnergy);
  ch.SetBranchAddress("HVCeF3", &HVCeF3, &b_HVCeF3);
  
  ch.SetBranchAddress("pos_2FibClust_hodoX1", &pos_2FibClust_hodoX1, &b_pos_2FibClust_hodoX1);
  ch.SetBranchAddress("pos_2FibClust_hodoX2", &pos_2FibClust_hodoX2, &b_pos_2FibClust_hodoX2);
  ch.SetBranchAddress("pos_2FibClust_hodoY1", &pos_2FibClust_hodoY1, &b_pos_2FibClust_hodoY1 );
  ch.SetBranchAddress("pos_2FibClust_hodoY2", &pos_2FibClust_hodoY2, &b_pos_2FibClust_hodoY2);
  
  ch.SetBranchAddress("cluster_pos_corr_hodoX1", &cluster_pos_corr_hodoX1, &b_cluster_pos_corr_hodoX1);
  ch.SetBranchAddress("cluster_pos_corr_hodoX2", &cluster_pos_corr_hodoX2, &b_cluster_pos_corr_hodoX2);
  ch.SetBranchAddress("cluster_pos_corr_hodoY1", &cluster_pos_corr_hodoY1, &b_cluster_pos_corr_hodoY1 );
  ch.SetBranchAddress("cluster_pos_corr_hodoY2", &cluster_pos_corr_hodoY2, &b_cluster_pos_corr_hodoY2);
  
  ch.SetBranchAddress("position_X1", &position_X1, &b_position_X1);
  ch.SetBranchAddress("position_X2", &position_X2, &b_position_X2);
  ch.SetBranchAddress("position_Y1", &position_Y1, &b_position_Y1 );
  ch.SetBranchAddress("position_Y2", &position_Y2, &b_position_Y2);
  
  ch.SetBranchAddress("nTDCHits", &nTDCHits, &b_nTDCHits);
  
  

  int nentries = ch.GetEntries();
  for(int  iEntry=0; iEntry<nentries; ++iEntry ) {
    ch.GetEntry( iEntry );     
    
    if(cef3_chaInt_corr->at(0) <50 || cef3_chaInt_corr->at(0)>30000000)
      continue;
    if(cef3_chaInt_corr->at(1) <50 || cef3_chaInt_corr->at(1)>30000000)
      continue;
    if(cef3_chaInt_corr->at(2) <50 || cef3_chaInt_corr->at(2)>30000000)
      continue;
    if(cef3_chaInt_corr->at(3) <50 || cef3_chaInt_corr->at(3)>30000000)
      continue;
    
    if( abs(cluster_pos_corr_hodoX1)<20 &&  abs(cluster_pos_corr_hodoX2)<20 && abs(cluster_pos_corr_hodoY1)<20 && abs(cluster_pos_corr_hodoY2)<20  && ( nTDCHits->at(0)>0 && nTDCHits->at(1)>0 && nTDCHits->at(2)>0 && nTDCHits->at(3)>0 && (nTDCHits->at(0)+ nTDCHits->at(1)+ nTDCHits->at(2)+ nTDCHits->at(3)   )<7 &&   nTDCHits->at(0)<3 && nTDCHits->at(1)<3 && nTDCHits->at(2)<3 && nTDCHits->at(3)<3    ) ){
      
      double xPos =  0.5 * ( position_X2 + position_X1) - (xTable-194.);
      double yPos =  0.5 * ( position_Y2 + position_Y1) + (yTable -254.);


      if( ((abs(xPos)<12.) && (abs(yPos)<12.) ) ){

      for(int k=0; k<4; k++){
	fracProf[k]->Fill( -xPos, -yPos,  cef3_chaInt_corr->at(k)/(cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1) +cef3_chaInt_corr->at(2)+ cef3_chaInt_corr->at(3) ) );
      }


	hprofR[0]->Fill( sqrt((xPos-12.)*(xPos-12.) + (yPos+12.)*(yPos+12.)) , cef3_chaInt_corr->at(0)/(cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1) +cef3_chaInt_corr->at(2)+ cef3_chaInt_corr->at(3)   ) );   
	
	hprofR[1]->Fill( sqrt((xPos+12.)*(xPos+12.) + (yPos+12.)*(yPos+12.)) , cef3_chaInt_corr->at(1) /(cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1) +cef3_chaInt_corr->at(2)+ cef3_chaInt_corr->at(3)   ) );
	
	hprofR[2]->Fill( sqrt((xPos+12.)*(xPos+12.) + (yPos-12.)*(yPos-12.)) , cef3_chaInt_corr->at(2) /(cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1) +cef3_chaInt_corr->at(2)+ cef3_chaInt_corr->at(3)   ) );
	
	hprofR[3]->Fill( sqrt((xPos-12.)*(xPos-12.) + (yPos-12.)*(yPos-12.)) , cef3_chaInt_corr->at(3) /(cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1) +cef3_chaInt_corr->at(2)+ cef3_chaInt_corr->at(3)   ) );
      }
      
    }
    
  }
  
  TFile *fOut=TFile::Open("respVsR.root","RECREATE");   //THE ENTRIES LOOOP      
  for (std::map<string,TObject*>::iterator it=outHistos.begin();it!=outHistos.end();++it)
    it->second->Write();
  fOut->Close();
  
  


 for(int j=0; j<4; j++){
   //gPad->SetLogz();
    gStyle->SetPalette(51,0);
    fracProf[j]->Draw("colz");

    //   gr_resp_vs_pos->SetMarkerSize(100);
    
    label_top2->Draw("same");
    lineUpper->Draw("same");lineLower->Draw("same"); 
    lineRight->Draw("same");lineLeft->Draw("same");
    
    chamferDR->Draw("same"); chamferUR->Draw("same"); 
    chamferDL->Draw("same"); chamferUL->Draw("same"); 
    

    c1->SaveAs( Form( "%s/FracProf_%d.pdf", outputdir.c_str(), j ) );
    c1->SaveAs( Form( "%s/FracProf_%d.png", outputdir.c_str(), j ) );
    
    
    c1->Clear();
 }
























 // int nBins = 32;



  //  std::map<string,TObject*> outHistos; 

  TProfile2D* hprof2d[4];
  
  for(int j =0; j<4; j++){ //le lüp over ze fibres ~(^.^)~
    
    hprof2d[j] = new TProfile2D(Form("hprof2d_%d",j),"profile of values vs pos y",nBins , -16, 16,nBins,-16,16 , 0, 50000000);
    hprof2d[j]->SetXTitle("Position X [mm]");
    hprof2d[j]->SetYTitle("Position Y [mm]");
    //  outHistos[hprofR[j]->GetName()]=(TObject*) hprofR[j];
  }
  
  
  
  TProfile2D* hprof2d_tot = new TProfile2D("hprof2d_tot","profile of values vs pos x&y",nBins, -16,16, nBins  ,-16,16, 50, 25000000);
  hprof2d_tot->SetXTitle("Position X [mm]");
  hprof2d_tot->SetYTitle("Position Y [mm]");
  
  //	gStyle->SetOptStat(1);
  TH1F *diff_X = new TH1F("diff_X","hodo-calc distance in _X", 80, -40 ,40 );
  diff_X->SetXTitle("Difference X [mm]");
  diff_X->SetYTitle("Entries");
  
  TH1F *diff_Y = new TH1F("diff_Y","hodo-calc distance in _Y", 80 , -40 ,40 );
  diff_Y->SetXTitle("Difference Y [mm]");
  diff_Y->SetYTitle("Entries");
  
  
  std::cout << "Total Number of Entries = " << nentries << std::endl;
  //LOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOP///////////////////
  // for(int  iEntry=0; iEntry<100000; ++iEntry ) {
  for(int  iEntry=0; iEntry<nentries; ++iEntry ) {
    ch.GetEntry( iEntry );     
    
    if(cef3_chaInt_corr->at(0) <100 || cef3_chaInt_corr->at(0)>3000000)
      continue;
    if(cef3_chaInt_corr->at(1) <100 || cef3_chaInt_corr->at(1)>3000000)
     continue;
    if(cef3_chaInt_corr->at(2) <100 || cef3_chaInt_corr->at(2)>3000000)
      continue;
    if(cef3_chaInt_corr->at(3) <100 || cef3_chaInt_corr->at(3)>3000000)
      continue;
    
   double xPos =  0.5 * ( position_X2 + position_X1) - (xTable-194);
   double yPos =  0.5 * ( position_Y2 + position_Y1) + (yTable -254);
   if(abs(yPos)>12 || abs(xPos)>12 ) continue;
    
 
   if( abs(cluster_pos_corr_hodoX1)<20 &&  abs(cluster_pos_corr_hodoX2)<20 && abs(cluster_pos_corr_hodoY1)<20 && abs(cluster_pos_corr_hodoY2)<20  && ( nTDCHits->at(0)>0 && nTDCHits->at(1)>0 && nTDCHits->at(2)>0 && nTDCHits->at(3)>0 && (nTDCHits->at(0)+ nTDCHits->at(1)+ nTDCHits->at(2)+ nTDCHits->at(3)   )<7 &&   nTDCHits->at(0)<3 && nTDCHits->at(1)<3 && nTDCHits->at(2)<3 && nTDCHits->at(3)<3    ) ){
     
   double xPosition =  0.5 * ( position_X2 + position_X1) - (xTable-194);
   double yPosition =  0.5 * ( position_Y2 + position_Y1) + (yTable -254);
 
    double frac[4];
    frac[0] = cef3_chaInt_corr->at(0)/(cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1) +cef3_chaInt_corr->at(2) + cef3_chaInt_corr->at(3)   );
    frac[1] = cef3_chaInt_corr->at(1)/(cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1) +cef3_chaInt_corr->at(2) + cef3_chaInt_corr->at(3)   );
    frac[2] = cef3_chaInt_corr->at(2)/(cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1) +cef3_chaInt_corr->at(2) + cef3_chaInt_corr->at(3)   );
    frac[3] = cef3_chaInt_corr->at(3)/(cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1) +cef3_chaInt_corr->at(2) + cef3_chaInt_corr->at(3)   );


    //find the biggest fraction:
    //  int maxFracInd = getMaxFrac( frac );
    //   int secondFracInd = getSecondFrac( frac );
    

  TH2D* histo= new TH2D("histo","is this onnn",nBins, -16,16, nBins  ,-16,16);
  histo->SetXTitle("Position X [mm]");
  histo->SetYTitle("Position Y [mm]");

  TH2D* histo2= new TH2D("histo2","is this onnn",nBins, -16,16, nBins  ,-16,16);
  histo2->SetXTitle("Position X [mm]");
  histo2->SetYTitle("Position Y [mm]");

  TH2D* histo_overlap = new TH2D("histo_overlap","is this onnn",nBins, -16,16, nBins  ,-16,16);
  histo_overlap->SetXTitle("Position X [mm]");
  histo_overlap->SetYTitle("Position Y [mm]");

  double prob = 1.;
  

  for(int binX = 1; binX<nBins+1; binX++){
    for(int binY = 1; binY<nBins+1; binY++){
      float xx =  fracProf[0]->GetXaxis()->GetBinCenter(binX);
      float yy =  fracProf[0]->GetYaxis()->GetBinCenter(binY);

      if(abs(yy)>12 || abs(xx)>12 ) continue;
    
 
   prob = 0.;
      for(int index=0; index<4; index++){
   	
	//    value =  1.- (  abs(frac[maxFracInd] -fracProf[maxFracInd]->GetBinContent(binX,binY))) / (  abs(frac[maxFracInd]  +fracProf[maxFracInd]->GetBinContent(binX,binY)))  ;
	//   histo->Fill(xx,yy, value );
	//  value2 = 1.- (  abs(frac[secondFracInd] -fracProf[secondFracInd]->GetBinContent(binX,binY))) / (  abs(frac[secondFracInd]  +fracProf[secondFracInd]->GetBinContent(binX,binY)))  ;
	//  histo2->Fill(xx,yy,value2);


	
      //	prob = prob *( 1.- (  abs(frac[index] -fracProf[index]->GetBinContent(binX,binY) )  ) )/fracProf[index]->GetBinContent(binX,binY)  ;

		//      	prob = prob *( 1.- sqrt(sqrt( abs(frac[index] -fracProf[index]->GetBinContent(binX,binY)) / abs(frac[index]  +fracProf[index]->GetBinContent(binX,binY) ) )) )  ;

	prob = prob + ((frac[index] -fracProf[index]->GetBinContent(binX,binY))* (frac[index] -fracProf[index]->GetBinContent(binX,binY)) *log(frac[index])   )   ;
	
	histo2->Fill(xx,yy,prob);
      }//end of loop over indices of channels
      
      histo_overlap->Fill(xx,yy,1.- prob);
      
      
    }//end of loop over y  
  }//end of loop over x

  // histo_overlap->SetMaximum(0.01);
    histo_overlap->SetMinimum(0.9);
  int binnie =  histo_overlap->GetMaximumBin();
  int xBinnie;
  int yBinnie;
  int zBinnie;
  histo_overlap->GetBinXYZ(binnie, xBinnie, yBinnie, zBinnie);

  float x =  histo_overlap->GetXaxis()->GetBinCenter(xBinnie);
  float y =  histo_overlap->GetYaxis()->GetBinCenter(yBinnie);
  
    
      if(iEntry >5550 &&iEntry < 5650){

      histo_overlap->Draw("colz");
      //   gr_resp_vs_pos->SetMarkerSize(100);   
      label_top2->Draw("same");
      lineUpper->Draw("same");lineLower->Draw("same"); 
      lineRight->Draw("same");lineLeft->Draw("same");
      
      chamferDR->Draw("same"); chamferUR->Draw("same"); 
      chamferDL->Draw("same"); chamferUL->Draw("same"); 
      
      c1->SaveAs( Form( "%s/TurnedOn_overlap_%d.png", outputdir.c_str(),iEntry ) );
     }
    
  delete histo;
  delete histo2;
  delete histo_overlap;
  
  //fracProf[j]
  //Calculating the position from the fractions
  //   float x;
  //   float y;
  
  //Let's fill some damn histos and see if I produced BS or not 
  for(int j=0; j<4; j++){
    hprof2d[j]->Fill( x , y , cef3_chaInt_corr->at(j) , 1);
  }
  hprof2d_tot->Fill(x,y, cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1)+ cef3_chaInt_corr->at(2)+ cef3_chaInt_corr->at(3) , 1 );
  

    
  diff_X->Fill(xPosition-x);
  diff_Y->Fill(yPosition-y);


  
   }//end of selection
  
  }//End of loop over entries
  

 for(int j=0; j<4; j++){
    //gPad->SetLogz();
    gStyle->SetPalette(51,0);
    hprof2d[j]->Draw("colz");
    //   gr_resp_vs_pos->SetMarkerSize(100); 
    label_top2->Draw("same");
    lineUpper->Draw("same");lineLower->Draw("same"); 
    lineRight->Draw("same");lineLeft->Draw("same");
    
    chamferDR->Draw("same"); chamferUR->Draw("same"); 
    chamferDL->Draw("same"); chamferUL->Draw("same"); 
    
    c1->SaveAs( Form( "%s/SingleChannel_resp_vs_pos_%d.pdf", outputdir.c_str(), j ) );
    c1->SaveAs( Form( "%s/SingleChannel_resp_vs_pos_%d.png", outputdir.c_str(), j ) );   
    c1->Clear();
 }

    gStyle->SetPalette(51,0);
    hprof2d_tot->Draw("colz");
    //   gr_resp_vs_pos->SetMarkerSize(100);   
    label_top2->Draw("same");
    lineUpper->Draw("same");lineLower->Draw("same"); 
    lineRight->Draw("same");lineLeft->Draw("same");
    
    chamferDR->Draw("same"); chamferUR->Draw("same"); 
    chamferDL->Draw("same"); chamferUL->Draw("same"); 
    
    c1->SaveAs( Form( "%s/SingleChannel_resp_vs_pos_tot.pdf", outputdir.c_str() ) );
    c1->SaveAs( Form( "%s/SingleChannel_resp_vs_pos_tot.png", outputdir.c_str() ) );
    
    
    diff_X->Draw();

     c1->SaveAs( Form( "%s/DiffX.png", outputdir.c_str() ) );
 

    diff_Y->Draw();

     c1->SaveAs( Form( "%s/DiffY.png", outputdir.c_str() ) );
 


    return 0;
}















  int getMaxFrac( double *frac){
  int i = -1;
  if(frac[0]>frac[1] && frac[0]>frac[2] && frac[0] >frac[3]) i = 0;
  if(frac[1]>frac[0] && frac[1]>frac[2] && frac[1] >frac[3]) i = 1;
  if(frac[2]>frac[1] && frac[2]>frac[0] && frac[2] >frac[3]) i = 2;
  if(frac[3]>frac[1] && frac[3]>frac[2] && frac[3] >frac[0]) i = 3;
  return i;
  }



int getSecondFrac(double *frac){
  int i = getMaxFrac( frac);
  int j = -1;

  for(int k=0; k<4; k++){
    if(k==i) continue;
    int counter = 0; //counts if you get 2 smaller ones
    for(int m=0; m<4; m++){
      if(m==i || k==i) continue;
      if(frac[k] > frac[m])
	counter = counter + 1;
      if(counter == 2) break;
    }
    if(counter==2) j=k;
  }




 
  return j;
  }











float findRadius(TF1* f,int order, double fraction)
{

  TF1* polDer=new TF1(f->GetName()+TString("_der"),Form("pol%d",order),1.5,34);
  double* par=f->GetParameters();
  for(int i=0;i<=order;++i)
      polDer->SetParameter(i,par[i]);
  polDer->SetParameter(0, f->GetParameter(0) - fraction);
  //find root assuming monothic derivative
  float step=0.05;
  float x=0;
  //int oldSign= 1;
  int oldSign= -1 + 2*(polDer->Eval(x)>0);
  while(x<34)
    {
      x+=step;
      int newSign=-1 + 2*(polDer->Eval(x)>0);
      if (newSign*oldSign==-1)
	break;
    }
  //  std::cout << "FOUND ZeroPoint @ R =  " << x << " Fit VALUE " << polDer->Eval(x) << std::endl;
  if(fraction>0.55){
    return 1.;
  } else {
  return x;
  }
  delete polDer;
}






float getRatioError( float num, float denom, float numErr, float denomErr ) {

  return sqrt( numErr*numErr/(denom*denom) + denomErr*denomErr*num*num/(denom*denom*denom*denom) );

}



















 /*




  //  std::map<string,TObject*> outHistos; 

  TProfile2D* hprof2d[4];
  
  for(int j =0; j<4; j++){ //le lüp over ze fibres ~(^.^)~
    
    hprof2d[j] = new TProfile2D(Form("hprof2d_%d",j),"profile of values vs pos y",32 , -16, 16,32,-16,16 , 0, 50000000);
  hprof2d[j]->SetXTitle("Position X [mm]");
  hprof2d[j]->SetYTitle("Position Y [mm]");
  //  outHistos[hprofR[j]->GetName()]=(TObject*) hprofR[j];
}
  
  
  
  TProfile2D* hprof2d_tot = new TProfile2D("hprof2d_tot","profile of values vs pos x&y",32, -16,16, 32  ,-16,16, 50, 25000000);
  hprof2d_tot->SetXTitle("Position X [mm]");
  hprof2d_tot->SetYTitle("Position Y [mm]");
  
  
  TFitResultPtr fit_0 = hprofR[0]->Fit("pol5","S R+","",1.5,34);
  TFitResultPtr fit_1 = hprofR[1]->Fit("pol5","S R+","",1.5,34);
  TFitResultPtr fit_2 = hprofR[2]->Fit("pol5","S R+","",1.5,34);
  TFitResultPtr fit_3 = hprofR[3]->Fit("pol5","S R+","",1.5,34);
  
  
 //for coordinate transfomation (= rotation by pi/4:

float sine = TMath::Sin(-3.141592 / 4.);
float cosine = TMath::Cos(-3.141592 / 4.);

 //LOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOP///////////////////
// for(int  iEntry=0; iEntry<200000; ++iEntry ) {
  for(int  iEntry=0; iEntry<nentries; ++iEntry ) {
   ch.GetEntry( iEntry );     
   
   if(cef3_chaInt_corr->at(0) <100 || cef3_chaInt_corr->at(0)>3000000)
     continue;
   if(cef3_chaInt_corr->at(1) <100 || cef3_chaInt_corr->at(1)>3000000)
     continue;
   if(cef3_chaInt_corr->at(2) <100 || cef3_chaInt_corr->at(2)>3000000)
     continue;
   if(cef3_chaInt_corr->at(3) <100 || cef3_chaInt_corr->at(3)>3000000)
     continue;
   
   double xPos =  0.5 * ( position_X2 + position_X1) - (xTable-194);
   double yPos =  0.5 * ( position_Y2 + position_Y1) + (yTable -254);
   if(abs(yPos)>12 || abs(xPos>12) ) continue;
    
  if( abs(cluster_pos_corr_hodoX1)<20 &&  abs(cluster_pos_corr_hodoX2)<20 && abs(cluster_pos_corr_hodoY1)<20 && abs(cluster_pos_corr_hodoY2)<20  && ( nTDCHits->at(0)>0 && nTDCHits->at(1)>0 && nTDCHits->at(2)>0 && nTDCHits->at(3)>0 && (nTDCHits->at(0)+ nTDCHits->at(1)+ nTDCHits->at(2)+ nTDCHits->at(3)   )<7 &&   nTDCHits->at(0)<3 && nTDCHits->at(1)<3 && nTDCHits->at(2)<3 && nTDCHits->at(3)<3    ) ){
     
    double fraction_0 = cef3_chaInt_corr->at(0)/(cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1) +cef3_chaInt_corr->at(2)+ cef3_chaInt_corr->at(3)   );
    double fraction_1 = cef3_chaInt_corr->at(1)/(cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1) +cef3_chaInt_corr->at(2)+ cef3_chaInt_corr->at(3)   );
    double fraction_2 = cef3_chaInt_corr->at(2)/(cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1) +cef3_chaInt_corr->at(2)+ cef3_chaInt_corr->at(3)   );
    double fraction_3 = cef3_chaInt_corr->at(3)/(cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1) +cef3_chaInt_corr->at(2)+ cef3_chaInt_corr->at(3)   );



 float r_0 =  findRadius( hprofR[0]->GetFunction("pol5"), 5, fraction_0);
 float r_1 =  findRadius( hprofR[1]->GetFunction("pol5"), 5, fraction_1);
 float r_2 =  findRadius( hprofR[2]->GetFunction("pol5"), 5, fraction_2);
 float r_3 =  findRadius( hprofR[3]->GetFunction("pol5"), 5, fraction_3);

 //in the pi/4. turned coordinate system this corresponds to:
 float d_0 = 0.5*(r_0 + (sqrt(2.)*24- r_2)) - sqrt(2.)*12. ;
 float d_1 = 0.5*(r_3 + (sqrt(2.)*24- r_1)) - sqrt(2.)*12. ;


 //In the (X,Y) coordinate system
 float x = cosine * d_0 - sine * d_1;
 float y = sine * d_0 + cosine * d_1;

 //Let's fill some damn histos and see if I produced BS or not 
 for(int j=0; j<4; j++){
 hprof2d[j]->Fill( x , y , cef3_chaInt_corr->at(j) , 1);
 }
 hprof2d_tot->Fill(x,y, cef3_chaInt_corr->at(0)+ cef3_chaInt_corr->at(1)+ cef3_chaInt_corr->at(2)+ cef3_chaInt_corr->at(3) , 1 );


  }//end of selection

 }//End of loop over entries


 for(int j=0; j<4; j++){
    //gPad->SetLogz();
    gStyle->SetPalette(51,0);
    hprof2d[j]->Draw("colz");

    //   gr_resp_vs_pos->SetMarkerSize(100);
    
    label_top2->Draw("same");
    lineUpper->Draw("same");lineLower->Draw("same"); 
    lineRight->Draw("same");lineLeft->Draw("same");
    
    chamferDR->Draw("same"); chamferUR->Draw("same"); 
    chamferDL->Draw("same"); chamferUL->Draw("same"); 
    

    c1->SaveAs( Form( "%s/SingleChannel_resp_vs_pos_%d.pdf", outputdir.c_str(), j ) );
    c1->SaveAs( Form( "%s/SingleChannel_resp_vs_pos_%d.png", outputdir.c_str(), j ) );
    
    
    c1->Clear();
 }


    gStyle->SetPalette(51,0);
    hprof2d_tot->Draw("colz");

    //   gr_resp_vs_pos->SetMarkerSize(100);
    
    label_top2->Draw("same");
    lineUpper->Draw("same");lineLower->Draw("same"); 
    lineRight->Draw("same");lineLeft->Draw("same");
    
    chamferDR->Draw("same"); chamferUR->Draw("same"); 
    chamferDL->Draw("same"); chamferUL->Draw("same"); 
    

    c1->SaveAs( Form( "%s/SingleChannel_resp_vs_pos_tot.pdf", outputdir.c_str() ) );
    c1->SaveAs( Form( "%s/SingleChannel_resp_vs_pos_tot.png", outputdir.c_str() ) );
    
 */
