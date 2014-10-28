#include <iostream>
#include <string>
#include <stdlib.h>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TEllipse.h"
#include "TString.h"
#include "TGaxis.h" 
#include "TAxis.h" 

#include "DrawTools.h"

#include "TApplication.h"

//This serves the purpose to see if the intercalibration was sort of successful,


int main ( int argc, char* argv[] ) {

  TApplication* a = new TApplication("a", 0, 0);

  std::string inputDir = "./analysisTrees";
  std::string runName = "259";
  std::string tag = "default";


  if( argc == 3 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
    std::string tag_str(argv[2]);
    tag = tag_str;
  }else if (argc==4){
    std::string inputDir_str(argv[2]);
    inputDir =inputDir_str;
  } else{
    std::cout<<"Usage:"<<std::endl;
    std::cout<<"./calibrateCef3 XX tag "<<std::endl;
    exit(12345);
  }

  if(tag!="default") inputDir = inputDir + "_"+tag;


  std::string outputdir = "CeF3IntercalPlots/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  DrawTools::setStyle();
  TGaxis::SetMaxDigits(4);


  TFile* file = TFile::Open(Form("%s/Reco_%s.root", inputDir.c_str(),runName.c_str()));
  std::cout<<"opening file:"<<file->GetName();

  TTree* tree = (TTree*)file->Get("recoTree");

  TCanvas* canny = new TCanvas("canny", "CeF3 Intercalibration", 300,400);

  canny->cd();

  //  gPad->SetLogy();

  int nBins = 100;
  // float xMin = 180000.;
  //  float xMax = 501000.;
  //float xMax = 375000.;

  float xMin = 18.;
  float xMax = 100.;

  TH1D* hu0 = new TH1D("hu0","",nBins,xMin,xMax);
  hu0->SetLineColor(kBlue);
  hu0->SetLineWidth(2.);
  TH1D* hu1 = new TH1D("hu1","",nBins,xMin,xMax);
  hu1->SetLineColor(kAzure+8);
  hu1->SetLineWidth(2.);
  TH1D* hu2 = new TH1D("hu2","",nBins,xMin,xMax);
  hu2->SetLineWidth(2.);
  hu2->SetLineColor(kGreen+1);
  TH1D* hu3 = new TH1D("hu3","",nBins,xMin,xMax);
  hu3->SetLineColor(kMagenta-4);
  hu3->SetLineWidth(2.);
  
  TH1D* h0 = new TH1D("h0","",nBins,xMin,xMax);
  h0->SetLineColor(kBlue);
  h0->SetLineWidth(2.);
  TH1D* h1 = new TH1D("h1","",nBins,xMin,xMax);
  h1->SetLineColor(kAzure+8);
  h1->SetLineWidth(2.);  
  TH1D* h2 = new TH1D("h2","",nBins,xMin,xMax);  
  h2->SetLineColor(kGreen+1);
  h2->SetLineWidth(2.);  
  TH1D* h3 = new TH1D("h3","",nBins,xMin,xMax);
  h3->SetLineColor(kMagenta-4);
  h3->SetLineWidth(2.);

  TH2D* h2_axes = new TH2D( "axes", "", 10, xMin , xMax, 10, 0.,0.09 );
  h2_axes->SetXTitle("CeF_{3} Response [ADC]");
  h2_axes->SetYTitle("");
  h2_axes->Draw("");

  tree->Draw("cef3[0]>>hu0","","same");
  tree->Draw("cef3[1]>>hu1","","same");
  tree->Draw("cef3[2]>>hu2","","same");
  tree->Draw("cef3[3]>>hu3","","same");

  hu0->Scale(1./hu0->Integral() );
  hu1->Scale(1./hu1->Integral() );
  hu2->Scale(1./hu2->Integral() );
  hu3->Scale(1./hu3->Integral() );

  TLegend* legend2 = new TLegend( 0.7, 0.9, 0.9, 0.9-0.06*4 );
  legend2->SetFillStyle(0);
  legend2->SetLineColor(0);
  legend2->SetLineWidth(0);
  legend2->SetTextSize(0.035);

  legend2->AddEntry( hu0, "Channel 0" );
  legend2->AddEntry( hu1, "Channel 1");
  legend2->AddEntry( hu2, "Channel 2" );
  legend2->AddEntry( hu3, "Channel 3" );
  //legend2->AddEntry( hu0, Form("mean0= %1.f\n #pm %1.f\n",hu0->GetMean(),hu0->GetMeanError()) );
  // legend2->AddEntry( hu1,Form("mean1= %1.f\n #pm %1.f\n",hu1->GetMean(),hu1->GetMeanError()) );
  //legend2->AddEntry( hu2, Form("mean2= %1.f\n #pm %1.f\n",hu2->GetMean(),hu2->GetMeanError()) );
  //legend2->AddEntry( hu3, Form("mean3= %1.f\n #pm %1.f\n",hu3->GetMean() ,hu3->GetMeanError() ));
  legend2->Draw("same");


  canny->SaveAs(Form("%s/NonCalib_%s_%s.pdf", outputdir.c_str(),runName.c_str(), tag.c_str() ));
  canny->SaveAs(Form("%s/NonCalib_%s_%s.png", outputdir.c_str(),runName.c_str(), tag.c_str() ));

  canny->Clear();


  //  gPad->SetLogy();
  h2_axes->Draw("");

  tree->Draw("cef3_corr[0]>>h0","","same");
  tree->Draw("cef3_corr[1]>>h1","","same");
  tree->Draw("cef3_corr[2]>>h2","","same");
  tree->Draw("cef3_corr[3]>>h3","","same");

  h0->Scale(1./h0->Integral() );
  h1->Scale(1./h1->Integral() );
  h2->Scale(1./h2->Integral() );
  h3->Scale(1./h3->Integral() );

  TLegend* legend2f = new TLegend( 0.7, 0.9, 0.9, 0.9-0.06*4);
  legend2f->SetFillStyle(0);
  legend2f->SetLineColor(0);
  legend2f->SetLineWidth(0);
  legend2f->SetTextSize(0.035);

  //  legend2f->AddEntry( h0, Form("mean0= %1.f\n #pm %1.f\n",h0->GetMean(),h0->GetMeanError()) );
  //  legend2f->AddEntry( h1,Form("mean1= %1.f\n #pm %1.f\n",h1->GetMean(),h1->GetMeanError()) );
  //  legend2f->AddEntry( h2, Form("mean2= %1.f\n #pm %1.f\n",h2->GetMean(),h2->GetMeanError()) );
  //  legend2f->AddEntry( h3, Form("mean3= %1.f\n #pm %1.f\n",h3->GetMean() ,h3->GetMeanError() ));
  
  legend2f->AddEntry( h0, "Channel 0" );
  legend2f->AddEntry( h1, "Channel 1" );
  legend2f->AddEntry( h2, "Channel 2" );
  legend2f->AddEntry( h3, "Channel 3" );

  legend2f->Draw("same");

  canny->SaveAs(Form("%s/Calib_%s_%s.pdf", outputdir.c_str(),runName.c_str(), tag.c_str() ));
  canny->SaveAs(Form("%s/Calib_%s_%s.png", outputdir.c_str(),runName.c_str(), tag.c_str() ));



  /*
//////////////////////////////////////////////////////////////////////////////
//only central events 

  TCanvas* canny2 = new TCanvas("canny2","",1200,800);
  canny2->cd();


  TH1D* FibX=new TH1D("FibX","",6,0,6);
  FibX->SetLineColor(kAzure+8);
  TH1D* FibY=new TH1D("FibY","",6,0,6);
  FibY->SetLineColor(kBlue);

  tree->Draw("nFibres_hodoClustX>>FibX","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1)");
  tree->Draw("nFibres_hodoClustY>>FibY","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1)","same");

  FibX->GetXaxis()->SetTitle("nFibres");

  TLegend* legend = new TLegend( 0.6, 0.9, 0.9, 0.6 );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->SetLineWidth(0);
  legend->SetTextSize(0.035);
  legend->AddEntry("FibX", "nFibres_hodoClustX" );
  legend->AddEntry("FibY" , "nFibres_hodoClustY" );
  legend->Draw("same");


  /////////////////////////////////////////////////
  //POSITION of the cluster and cut to it to get only central events

    std::vector<std::string> runs; 
  std::vector<float> beamEnergy;
 
  runs.push_back("BTF_314_20140503-024715_beam");
  beamEnergy.push_back(98.3);
  runs.push_back("BTF_308_20140503-002534_beam");
  beamEnergy.push_back(147.4);
  runs.push_back("BTF_293_20140502-180258_beam");
  beamEnergy.push_back(196.5);
  runs.push_back("BTF_286_20140502-153528_beam");
  beamEnergy.push_back(294.8);
  runs.push_back("BTF_259_20140502-012847_beam");
  beamEnergy.push_back(491.4);
  
  TGraphErrors* gr_mean_x = new TGraphErrors(0);
  TGraphErrors* gr_mean_y = new TGraphErrors(0);

  TGraphErrors* gr_rms_x = new TGraphErrors(0);
  TGraphErrors* gr_rms_y = new TGraphErrors(0);

  TCanvas* canny3 = new TCanvas("canny3","",600,600);
  canny3->cd();
  DrawTools::setStyle();
  TGaxis::SetMaxDigits(3);


  for( unsigned i=0; i<runs.size(); ++i ) {

    TFile* file2 = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(), runs[i].c_str()));
    TTree* tree2 = (TTree*)file2->Get("recoTree");

    float energy = beamEnergy[i];
    float energyErr = 5.;


  TH2D* h2_axes2 = new TH2D( "axes", "", 18, -4.5,4.5, 10, 0., 17000 );
  h2_axes2->SetXTitle("Position [mm]");
  h2_axes2->SetYTitle("Entries");
  h2_axes2->Draw("");

  TH1D* PosX=new TH1D("PosX","",20,-5,5);
  PosX->SetLineColor(46);
  TH1D* PosY=new TH1D("PosY","",20,-5,5);
  PosY->SetLineColor(38); 
  //Clusters with only 2 Fibres:
  //  TH1D* PosX_2=new TH1D("PosX_2","",20,-5,5); //  PosX_2->SetLineColor(kMagenta+2);
  //  TH1D* PosY_2=new TH1D("PosY_2","",18,-5,5); //  PosY_2->SetLineColor(kTeal+3);

  tree2->Draw("pos_hodoClustX>>PosX","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1)","same");
  // tree2->Draw("pos_hodoClustY>>PosY","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1 && pos_hodoClustY<2.5 && -2.5<pos_hodoClustY)","same");
    tree2->Draw("pos_hodoClustY>>PosY","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1 && pos_hodoClustY<3. && -3.<pos_hodoClustY)","same");  
  // tree->Draw("pos_hodoClustX>>PosX_2","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1 && nFibres_hodoClustX< 3. && nFibres_hodoClustY< 3. )","same");   // tree->Draw("pos_hodoClustY>>PosY_2","(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1  && nFibres_hodoClustX< 3. && nFibres_hodoClustY< 3.)","same");
  PosX->Rebin();   PosY->Rebin();
  PosX->SetLineWidth(2.);  PosY->SetLineWidth(2.);

  TPaveText* label_top = new TPaveText();
  label_top = DrawTools::getLabelTop(Form("%.0f MeV Electron Beam", energy));
  label_top->Draw("same");

  TLegend* legend3 = new TLegend( 0.6, 0.9, 0.9, 0.75 );
  legend3->SetFillColor(0);
  legend3->SetFillStyle(0);
  legend3->SetLineColor(0);
  legend3->SetLineWidth(0);
  legend3->SetTextSize(0.035);  //  legend3->AddEntry("PosX","X" );
  legend3->AddEntry("PosX", Form("#mu_{X} = %.2f #pm %.2f",PosX->GetMean(),PosX->GetMeanError()), "L");  // legend3->AddEntry("PosY", "Y");
  legend3->AddEntry("PosY", Form("#mu_{Y} = %.2f #pm %.2f",PosY->GetMean(),PosY->GetMeanError()), "L");
  legend3->Draw("same");

  canny3->SaveAs(Form("%s/Cluster_Position_%.0f_%s.pdf", outputdir.c_str(), energy, tag.c_str() ));

  gr_mean_x->SetPoint( i, energy, PosX->GetMean() );
  gr_mean_x->SetPointError( i, energyErr, PosX->GetMeanError() );
  gr_mean_y->SetPoint( i, energy, PosY->GetMean() );
  gr_mean_y->SetPointError( i, energyErr, PosY->GetMeanError() );

  gr_rms_x->SetPoint( i, energy, PosX->GetRMS() );
  gr_rms_x->SetPointError( i, energyErr, PosX->GetRMSError() );
  gr_rms_y->SetPoint( i, energy, PosY->GetRMS() );
  gr_rms_y->SetPointError( i, energyErr, PosY->GetRMSError() );


  std::cout << "mean x = " << PosX->GetMean() << std::endl;

	canny3->Clear();
}



  */






    return 0;
}



