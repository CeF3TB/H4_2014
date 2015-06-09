#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TProfile.h"


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
/// Frontface Scan: Fit the runs at different positions and fill that mean /////
/// in a TH2D to produce a heatmap of that /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



struct FitStruct {

  double mean;
  double mean_err;

  double sigma;
  double sigma_err;

  double reso;
  double reso_err;

};


FitStruct getFitResults( const std::string& outputdir, TTree* tree, const std::string& name, float xMin, float xMax, std::string& cut, std::string& whatToProject );


float getRatioError( float num, float denom, float numErr, float denomErr );


TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree, float HV );

FitStruct addPhotoStatistics( FitStruct rs );



int main( int argc, char* argv[] ) {


  std::string tag="V01";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }
  std::cout << "-> Using tag: " << tag << std::endl;

  std::string outputdir = "RespVsXnY/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  DrawTools::setStyle();
//  TGaxis::SetMaxDigits(3);


  //  std::string whatToProjectSim = Form("cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]");
  //  std::string whatToProject = Form("cef3_chaInt[0]+cef3_chaInt[1]+cef3_chaInt[2]+cef3_chaInt[3]");

    /*
  std::vector<std::string> simulation; 
  std::vector<float> beamEnergySimulation; 
 ///////SIMULATION////////////
  simulation.push_back("Simulation50");
  beamEnergySimulation.push_back(50000.);
  simulation.push_back("Simulation100");
  beamEnergySimulation.push_back(100000.);

  simulation.push_back("Simulation10");
  beamEnergySimulation.push_back(10000.);
  simulation.push_back("Simulation15");
  beamEnergySimulation.push_back(15000.);
  simulation.push_back("Simulation20");
  beamEnergySimulation.push_back(20000.);

  simulation.push_back("Simulation150");
  beamEnergySimulation.push_back(150000.);
  simulation.push_back("Simulation200");
  beamEnergySimulation.push_back(200000.);
  
  TFile* energyfileS = TFile::Open("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/H4UseThis/Reco_Simulation50.root" );
  TTree* energytreeS = (TTree*)energyfileS->Get("recoTree");
  
  std::string SimCut = "";
  FitStruct fsSim= getFitResults( outputdir, energytreeS , simulation[0] ,  5000. ,  20000., SimCut , whatToProjectSim );
  float energyS = fsSim.mean;
    */

  //DATA RUNS// Let's start with the runs at 50 GeV and 950 V
    // TGraphErrors* gr_resp_vs_pos = new TGraphErrors(0);

    std::string cut = "";
    
    
    TCanvas* c1 = new TCanvas( "c1", "",0,0, 1200, 1200 );
    c1->cd();
    TH2D* h2_axes = new TH2D( "axes", "", 24 , -12, 12,10,1000, 160000  );
    //  h2_axes->SetXTitle("Position X [mm]");
    //    h2_axes->SetYTitle("Position Y [mm]");
    // h2_axes->Draw("");

    /*

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
  */

  TPaveText* label_top2 = new TPaveText();
  label_top2 = DrawTools::getLabelTop("50 GeV Electron Beam");


  std::map<string,TObject*> outHistos; 

  TProfile* hprofX[4];
  TProfile* hprofY[4];
  
  for(int j =0; j<4; j++){ //le lüp over ze fibres ~(^.^)~
    
    hprofX[j] = new TProfile(Form("hprofX_%d",j),"profile of values vs pos y",60 , -15, 15, 50, 500000);
    hprofX[j]->SetXTitle("Position X [mm]");
    hprofX[j]->SetYTitle("Response [ADC]");
    outHistos[hprofX[j]->GetName()]=(TObject*) hprofX[j];
    hprofY[j] = new TProfile(Form("hprofY_%d",j),"profile of values vs pos y",60 , -15, 15, 50,  500000) ;
    hprofY[j]->SetXTitle("Position Y [mm]");
    hprofY[j]->SetYTitle("Response [ADC]");
    outHistos[hprofY[j]->GetName()]=(TObject*) hprofY[j];
  }
  
  TProfile* hprofX_tot = new TProfile("hprofX_tot","profile of values vs pos y",60 , -15, 15, 50, 2500000);
  hprofX_tot->SetXTitle("Position X [mm]");
  hprofX_tot->SetYTitle("Response [ADC]");
  outHistos[hprofX_tot->GetName()]=(TObject*) hprofX_tot;
  
  TProfile* hprofY_tot = new TProfile("hprofY_tot","profile of values vs pos y",60 , -15, 15, 50, 2500000);
  hprofY_tot->SetXTitle("Position Y [mm]");
  hprofY_tot->SetYTitle("Response [ADC]");  
  outHistos[hprofY_tot->GetName()]=(TObject*) hprofY_tot;
  


  //   int i = 323;
  //  for(int i = 274; i< 276 ; i++){
  //  for(int i = 274; i< 323; i++){

  for(int i = 329; i< 354; i++){ //100 GeV


      std::string name = Form("%d", i);

      TFile* file = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%d.root", tag.c_str(),  i ));
      TTree* tree = (TTree*)file->Get("recoTree");

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

      int nClusters_hodoX1;
      int nClusters_hodoX2;
      int nClusters_hodoY1;
      int nClusters_hodoY2;
      /*
      float   pos_corr_hodoX1;
      float   pos_corr_hodoY1;
      float   pos_corr_hodoX2;
      float   pos_corr_hodoY2;
      */
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
 
      TBranch *b_nClusters_hodoX1;
      TBranch *b_nClusters_hodoX2;
      TBranch *b_nClusters_hodoY1;
      TBranch *b_nClusters_hodoY2;
      /*
      TBranch *b_pos_corr_hodoX1;
      TBranch *b_pos_corr_hodoY1;
      TBranch *b_pos_corr_hodoX2;
      TBranch *b_pos_corr_hodoY2;
      */
      //Set Object Pointer
      cef3_chaInt = 0;
      cef3_maxAmpl = 0;
      cef3_chaInt_corr = 0;
      cef3_maxAmpl_corr = 0;
      nTDCHits = 0;

      //Set branch addresses
      tree->SetBranchAddress("wc_x", &wc_x, &b_wc_x);
      tree->SetBranchAddress("wc_y", &wc_y, &b_wc_y);
      tree->SetBranchAddress("wc_x_corr", &wc_x_corr, &b_wc_x_corr);
      tree->SetBranchAddress("wc_y_corr", &wc_y_corr, &b_wc_y_corr);
      tree->SetBranchAddress("cef3_chaInt", &cef3_chaInt, &b_cef3_chaInt);
      tree->SetBranchAddress("cef3_maxAmpl", &cef3_maxAmpl, &b_cef3_maxAmpl);
      tree->SetBranchAddress("cef3_chaInt_corr", &cef3_chaInt_corr, &b_cef3_chaInt_corr);
      tree->SetBranchAddress("cef3_maxAmpl_corr", &cef3_maxAmpl_corr, &b_cef3_maxAmpl_corr);
      tree->SetBranchAddress("xTable", &xTable, &b_xTable);
      tree->SetBranchAddress("yTable", &yTable, &b_yTable);
 
      tree->SetBranchAddress("beamEnergy", &beamEnergy, &b_beamEnergy);
      tree->SetBranchAddress("HVCeF3", &HVCeF3, &b_HVCeF3);

      tree->SetBranchAddress("pos_2FibClust_hodoX1", &pos_2FibClust_hodoX1, &b_pos_2FibClust_hodoX1);
      tree->SetBranchAddress("pos_2FibClust_hodoX2", &pos_2FibClust_hodoX2, &b_pos_2FibClust_hodoX2);
      tree->SetBranchAddress("pos_2FibClust_hodoY1", &pos_2FibClust_hodoY1, &b_pos_2FibClust_hodoY1 );
      tree->SetBranchAddress("pos_2FibClust_hodoY2", &pos_2FibClust_hodoY2, &b_pos_2FibClust_hodoY2);

      tree->SetBranchAddress("cluster_pos_corr_hodoX1", &cluster_pos_corr_hodoX1, &b_cluster_pos_corr_hodoX1);
      tree->SetBranchAddress("cluster_pos_corr_hodoX2", &cluster_pos_corr_hodoX2, &b_cluster_pos_corr_hodoX2);
      tree->SetBranchAddress("cluster_pos_corr_hodoY1", &cluster_pos_corr_hodoY1, &b_cluster_pos_corr_hodoY1 );
      tree->SetBranchAddress("cluster_pos_corr_hodoY2", &cluster_pos_corr_hodoY2, &b_cluster_pos_corr_hodoY2);

      tree->SetBranchAddress("position_X1", &position_X1, &b_position_X1);
      tree->SetBranchAddress("position_X2", &position_X2, &b_position_X2);
      tree->SetBranchAddress("position_Y1", &position_Y1, &b_position_Y1 );
      tree->SetBranchAddress("position_Y2", &position_Y2, &b_position_Y2);

      tree->SetBranchAddress("nClusters_hodoX1", &nClusters_hodoX1, &b_nClusters_hodoX1);
      tree->SetBranchAddress("nClusters_hodoX2", &nClusters_hodoX2, &b_nClusters_hodoX2);
      tree->SetBranchAddress("nClusters_hodoY1", &nClusters_hodoY1, &b_nClusters_hodoY1 );
      tree->SetBranchAddress("nClusters_hodoY2", &nClusters_hodoY2, &b_nClusters_hodoY2);
      
      tree->SetBranchAddress("nTDCHits", &nTDCHits, &b_nTDCHits);

      int nentries = tree->GetEntries();
      for(int  iEntry=0; iEntry<nentries; ++iEntry ) {
	tree->GetEntry( iEntry );     
     //      std::cout << beamEnergy << std::endl;
	
	//	double deposite = 
	//double deposite_tot = cef3_chaInt->at(0) + cef3_chaInt->at(1) + cef3_chaInt->at(2) + cef3_chaInt->at(3) ;
	//double deposite = cef3_maxAmpl_corr->at(j);
	//	double deposite_tot = cef3_maxAmpl_corr->at(0) + cef3_maxAmpl_corr->at(1) + cef3_maxAmpl_corr->at(2) + cef3_maxAmpl_corr->at(3) ;
	double deposite_tot = cef3_chaInt_corr->at(0) + cef3_chaInt_corr->at(1) + cef3_chaInt_corr->at(2) + cef3_chaInt_corr->at(3) ;
		
	//	double xPos =  0.5 * ( cluster_pos_corr_hodoX2 + cluster_pos_corr_hodoX1) ;
	//	double yPos = 0.5 * ( cluster_pos_corr_hodoY2 + cluster_pos_corr_hodoY1) ;
	
	// if((pos_corr_hodoX1-pos_corr_hodoX2) < 1. && (pos_corr_hodoY1 - pos_corr_hodoY2) < 1.)
	  if(cef3_chaInt_corr->at(0) <100 || cef3_chaInt_corr->at(0)>30000000)
	     continue;
	  if(cef3_chaInt_corr->at(1) <100 || cef3_chaInt_corr->at(1)>30000000)
	     continue;
	  if(cef3_chaInt_corr->at(2) <100 || cef3_chaInt_corr->at(2)>30000000)
	     continue;
	  if(cef3_chaInt_corr->at(3) <100 || cef3_chaInt_corr->at(3)>30000000)
	     continue;
       
	if( abs(cluster_pos_corr_hodoX1)<20 &&  abs(cluster_pos_corr_hodoX2)<20 && abs(cluster_pos_corr_hodoY1)<20 && abs(cluster_pos_corr_hodoY2)<20  && ( nTDCHits->at(0)>0 && nTDCHits->at(1)>0 && nTDCHits->at(2)>0 && nTDCHits->at(3)>0 && (nTDCHits->at(0)+ nTDCHits->at(1)+ nTDCHits->at(2)+ nTDCHits->at(3)   )<7 &&   nTDCHits->at(0)<3 && nTDCHits->at(1)<3 && nTDCHits->at(2)<3 && nTDCHits->at(3)<3    ) ){
	  
       double xPos =  0.5 * ( position_X2 + position_X1) - (xTable-194);
       double yPos =  0.5 * ( position_Y2 + position_Y1) + (yTable -254);

       if(abs(yPos)<3){
 	 for (int j=0;j<4;++j)
	   hprofX[j]->Fill( -xPos , cef3_chaInt_corr->at(j) );

	 hprofX_tot->Fill( -xPos , deposite_tot);
	 
       }
       if(abs(xPos)<3){
	 for (int j=0;j<4;++j)
	   hprofY[j]->Fill( -yPos , cef3_chaInt_corr->at(j) );

	 hprofY_tot->Fill( -yPos , deposite_tot);

       }
	 
	}
	
      }

  }//THE ENTRIES LOOOP
  
  // TFile *fOut=TFile::Open("respVsXnY.root","RECREATE");         

 TFile *fOut=TFile::Open("respVsXnY_100GeV.root","RECREATE");         
 for (std::map<string,TObject*>::iterator it=outHistos.begin();it!=outHistos.end();++it)
   it->second->Write();
 fOut->Close();
 

/*
 hprofX->SetLineWidth(2);
   //gPad->SetLogy();
   //   gStyle->SetPalette(51,0);
   hprofX->Draw("");
   //   gr_resp_vs_pos->SetMarkerSize(100);
   
   label_top2->Draw("same");
 
   TLegend* leg4 = new TLegend(0.3 , 0.16, 0.7, 0.3);
 leg4->SetTextSize(0.038);
 leg4->AddEntry(hprofX ,Form("Ch %d, charge Integrated",j) , "p");
 leg4->SetFillColor(0);
 leg4->Draw("same");
  
   
   c1->SaveAs( Form( "%s/chaInt_Xresp_vs_pos_%d.pdf", outputdir.c_str(),j ) );
   c1->SaveAs( Form( "%s/chaInt_Xresp_vs_pos_%d.png", outputdir.c_str(),j ) );
   
   
   c1->Clear();

 hprofY->SetLineWidth(2);   
   //  gStyle->SetPalette(51,0);
   hprofY->Draw("profile");
   //   gr_resp_vs_pos->SetMarkerSize(100);

   TLegend* legY = new TLegend(0.3 , 0.16, 0.7, 0.3);
 legY->SetTextSize(0.038);
 legY->AddEntry(hprofY ,Form("Ch %d, charge Integrated",j) , "p");
 legY->SetFillColor(0);
 legY->Draw("same");
   
  
   label_top2->Draw("same");
   c1->SaveAs( Form( "%s/chaInt_Yresp_vs_pos_%d.pdf", outputdir.c_str(),j ) );
   c1->SaveAs( Form( "%s/chaInt_Yresp_vs_pos_%d.png", outputdir.c_str(),j ) );
   
   c1->Clear();
  }

 hprofY_tot->SetLineWidth(2);
   //  gStyle->SetPalette(51,0);
   hprofY_tot->Draw("profile");
   //   gr_resp_vs_pos->SetMarkerSize(100);
   
   label_top2->Draw("same");

   TLegend* legY_tot = new TLegend(0.3 , 0.16, 0.7, 0.3);
 legY_tot->SetTextSize(0.038);
 legY_tot->AddEntry(hprofY_tot ,"Sum, charge Integrated", "p");
 legY_tot->SetFillColor(0);
 legY_tot->Draw("same");
  

   c1->SaveAs( Form( "%s/chaInt_Yresp_vs_pos_tot.pdf", outputdir.c_str() ) );
   c1->SaveAs( Form( "%s/chaInt_Yresp_vs_pos_tot.png", outputdir.c_str() ) );
  

   c1->Clear();

   hprofX_tot->SetLineWidth(2);
   //  gStyle->SetPalette(51,0);
   hprofX_tot->Draw("profile");
   //   gr_resp_vs_pos->SetMarkerSize(100);
   
   label_top2->Draw("same"); 

   TLegend* legX_tot = new TLegend(0.3 , 0.16, 0.7, 0.3);
 legX_tot->SetTextSize(0.038);
 legX_tot->AddEntry(hprofX_tot ,"Sum, charge Integrated", "p");
 legX_tot->SetFillColor(0);
 legX_tot->Draw("same");
  
   c1->SaveAs( Form( "%s/chaInt_Xresp_vs_pos_tot.pdf", outputdir.c_str() ) );
   c1->SaveAs( Form( "%s/chaInt_Xresp_vs_pos_tot.png", outputdir.c_str() ) );

   c1->Clear();

*/











      /*

  TProfile* hprofX_sim_tot = new TProfile("hprofX_sim_tot","profile of values vs pos y",60 , -15, 15, 50, 2500000);
  hprofX_sim_tot->SetXTitle("Position X [mm]");
  hprofX_sim_tot->SetYTitle("Response [ADC]");
  
  TProfile* hprofY_sim_tot = new TProfile("hprofY_sim_tot","profile of values vs pos y",60 , -15, 15, 50, 2500000) ;
  hprofY_sim_tot->SetXTitle("Position Y [mm]");
  hprofY_sim_tot->SetYTitle("Response [ADC]");  
   

     std::string nameSim = Form("%d", 50);

      TFile* fileSim = TFile::Open(Form("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/H4UseThis1_5mmGap30mmBeam/Reco_Simulation50.root" ));
      TTree* treeSim = (TTree*)fileSim->Get("recoTree");
 
     
      TChain chain("recoTree");
      chain.Add(Form("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/H4UseThis1_5mmGap30mmBeam/Reco_Simulation50.root" ));
      


      float xPos;
      float yPos;
      float   *cef3_corr[4];
            
      //Branches
      TBranch *b_xPos;
      TBranch *b_yPos;
      TBranch *b_cef3_corr;
      
      //Set Object Pointer
      //cef3_corr = 0;
      //Set branch addresses
      treeSim->SetBranchAddress("cef3_corr", &cef3_corr, &b_cef3_corr);
      treeSim->SetBranchAddress("xPos", &xPos, &b_xPos);
      treeSim->SetBranchAddress("yPos", &yPos, &b_yPos);

 
      //THE ENTRIES LOOOP   
      int nentriesSim = treeSim->GetEntries();
      for(int  iEntry=0; iEntry<nentriesSim; ++iEntry ) {
	treeSim->GetEntry( iEntry );     
     	
	//	double deposite_tot = 0 ;
	//		double deposite_tot = cef3_corr[0] + cef3_corr[1] + cef3_corr[2] + cef3_corr[3] ;


	  
       double xPosition =  xPos;
       double yPosition = yPos;

       if(abs(yPosition)<5){
	  hprofX_sim_tot->Fill( -xPosition , deposite_tot);
       }
       if(abs(xPosition)<5){
	  hprofY_sim_tot->Fill( -yPosition , deposite_tot);
       }

	}
	
      
   
 hprofY_sim_tot->SetLineWidth(2);
   //  gStyle->SetPalette(51,0);
   hprofY_sim_tot->Draw("profile");
   //   gr_resp_vs_pos->SetMarkerSize(100);
   
   label_top2->Draw("same");

   TLegend* legY_sim_tot = new TLegend(0.3 , 0.16, 0.7, 0.3);
 legY_sim_tot->SetTextSize(0.038);
 legY_sim_tot->AddEntry(hprofY_sim_tot ,"Sum, charge Integrated", "p");
 legY_sim_tot->SetFillColor(0);
 legY_sim_tot->Draw("same");
  
   c1->SaveAs( Form( "%s/chaInt_Yresp_vs_pos_sim_tot.pdf", outputdir.c_str() ) );
   c1->SaveAs( Form( "%s/chaInt_Yresp_vs_pos_sim_tot.png", outputdir.c_str() ) );
  

   c1->Clear();
   hprofX_sim_tot->SetLineWidth(2);
 hprofX_sim_tot->Draw("profile");
  
   label_top2->Draw("same"); 

   TLegend* legX_sim_tot = new TLegend(0.3 , 0.16, 0.7, 0.3);
 legX_sim_tot->SetTextSize(0.038);
 legX_sim_tot->AddEntry(hprofX_sim_tot ,"Sum, charge Integrated", "p");
 legX_sim_tot->SetFillColor(0);
 legX_sim_tot->Draw("same");
  
   c1->SaveAs( Form( "%s/chaInt_Xresp_vs_pos_sim_tot.pdf", outputdir.c_str() ) );
   c1->SaveAs( Form( "%s/chaInt_Xresp_vs_pos_sim_tot.png", outputdir.c_str() ) );

   c1->Clear();

*/

    return 0;
}






















float getRatioError( float num, float denom, float numErr, float denomErr ) {

  return sqrt( numErr*numErr/(denom*denom) + denomErr*denomErr*num*num/(denom*denom*denom*denom) );

}


FitStruct getFitResults( const std::string& outputdir, TTree* tree, const std::string& name, float xMin, float xMax, std::string& cut, std::string& whatToProject ) {

  DrawTools::setStyle();

  std::string histoName(Form("h1_%s", name.c_str() ));
  std::string histoName2(Form("h2_%s", name.c_str() ));
  TH1D* h1;
  TF1* f1;
  
  gStyle->SetOptFit(1);
  
  h1 = new TH1D(histoName.c_str(), "", 5000, xMin, xMax);
  tree->Project( histoName.c_str(), whatToProject.c_str(), cut.c_str() );
  f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", xMin, xMax);
  
  h1->Fit( f1, "RQN" );
  
  f1->SetLineColor(kRed);

  
  double peakpos = h1->GetMean();
  double sigma = h1->GetRMS();
  //  double peakpos = f1->GetParameter(1);
  //  double sigma = f1->GetParameter(2);
  double fitmin;
  double fitmax;
  
  if( (peakpos-8*sigma) < 4000 ){ 
    fitmin = 4000;
  } else{ 
    fitmin = peakpos-8*sigma;
  }

  if(((peakpos>0) && (sigma > 2*peakpos)) ) { 
    fitmax = 2*peakpos;
  } else{ 
    fitmax = peakpos+5*sigma;
  }

  if (sigma < 500) sigma = 500;



  TH1D* h2 = new TH1D("h2", "", 200, fitmin,fitmax);
  tree->Project( "h2", whatToProject.c_str(), cut.c_str() );

  RooRealVar x("x","ADC Channels", fitmin, fitmax);
  RooDataHist data("data","dataset with x",x,Import(*h2) );
  
  RooPlot* frame;
  RooPlot* xframe = x.frame();   
  
  frame = x.frame("Title");
  data.plotOn(frame);  //this will show histogram data points on canvas
  //  data.statOn(frame);  //this will display hist stat on canvas
  
  RooRealVar meanr("meanr","Mean",peakpos,peakpos-sigma, peakpos+ sigma);
  RooRealVar width("width","#sigma",sigma *0.3, 150.0, 5.*sigma);
  RooRealVar A("A","Dist",2., 0.0, 7.0);
  RooRealVar N("N","Deg",5, 0.0, 10);

  meanr.setRange( 30000. , 1000000.);
  width.setRange(500, 22000);
  
  RooCBShape fit_fct("fit_fct","fit_fct",x,meanr,width,A,N); int ndf = 4;

  fit_fct.fitTo(data);
  
  fit_fct.plotOn(frame,LineColor(4));//this will show fit overlay on canvas
  // fit_fct.paramOn(frame); //this will display the fit parameters on canvas

  double mean = meanr.getVal();
  double meanErr = meanr.getError();

  double rms = width.getVal();
  double rmsErr = width.getError();

  double reso = 100.* rms/mean; //in percent
  double resoErr = 100.* getRatioError( rms, mean, meanErr, rmsErr );

  
  TCanvas* cans = new TCanvas("cans", "un canvas", 600,600);
  cans->cd();
  frame->Draw(); 

  TLegend* lego = new TLegend(0.6, 0.7, 0.9, 0.92);  
  lego->SetTextSize(0.038);
  lego->AddEntry(  (TObject*)0 ,Form("#mu = %.0f #pm %.0f", meanr.getVal(), meanr.getError() ), "");
  lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.0f #pm %.0f ", width.getVal(), width.getError() ), "");
  lego->AddEntry(  (TObject*)0 ,Form("#chi^{2} = %.2f / %d ", frame->chiSquare(ndf) , ndf ), "");
  lego->AddEntry(  (TObject*)0 ,Form("#sigma/#mu = %.1f #pm %.1f %s ", reso , resoErr ,"%"), "");
  lego->SetFillColor(0);
  lego->Draw("same");

  cans->SaveAs( Form( "%s/CBFit_%s.png", outputdir.c_str(), name.c_str() ) );
  

  
  FitStruct fs;
  fs.mean = mean;
  fs.mean_err = meanErr;
  fs.sigma = rms;
  fs.sigma_err = rmsErr;
  fs.reso = reso;
  fs.reso_err = resoErr;


  delete h1;
  delete f1;
  delete h2;
  delete cans;

  return fs;
}




FitStruct addPhotoStatistics( FitStruct rs ) {

  // MC response is already in MeV, 0.49 is the number of p.e. per MeV
  // RM = 90% of energy, so assume 90% of energy (0.9*491) is deposited in central channel
  float  nADC = rs.mean/171.6*3236;
  //  float  nADC = rs.resp/169.445*3224.36;
  float nPhotoElectrons = nADC/27.3;
  float poissonError = 100./sqrt( nPhotoElectrons ); // in percent
 
  float resoUnsmeared = rs.reso;
  float resoSmeared =  sqrt( rs.reso*rs.reso + poissonError*poissonError );

  rs.reso = resoSmeared;
  rs.reso_err = rs.reso_err * resoSmeared / resoUnsmeared; 

  return rs;
} 
