#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "Math/Polynomial.h"
#include "../CommonTools/interface/DrawTools.h"
#include "../CommonTools/src/DrawTools.cc"

#include <iostream>

float findMaximum(TF1* f,int order)
{
  TF1* polDer=new TF1(f->GetName()+TString("_der"),Form("pol%d",order-1),-9,9);
  double* par=f->GetParameters();
  for(int i=1;i<=order;++i)
      polDer->SetParameter(i-1,i*par[i]);
  //find root assuming monothic derivative
  float step=0.05;
  float x=-9;
  int oldSign=-1 + 2*(polDer->Eval(x)>0);
  while(x<9)
    {
      x+=step;
      int newSign=-1 + 2*(polDer->Eval(x)>0);
      if (newSign*oldSign==-1)
	break;
    }
  std::cout << "FOUND MAXIMUM @ " << x << " VALUE " << f->Eval(x) << std::endl;
  return f->Eval(x);
}

TProfile* rescaleProfile(TProfile *p)
{
  TFitResultPtr fit=p->Fit("pol4","R+","",-9,9);
  float max=findMaximum(p->GetFunction("pol4"),4);
  TProfile* p1=(TProfile*)p->Clone(p->GetName()+TString("_rescaled"));
  p1->Scale(1./max);
  p1->GetYaxis()->SetTitle("Normalized response [a.u.]");
  return p1;
}

void dataMCProfiles()
{

  gStyle->SetOptTitle(0);
  DrawTools::setStyle();
  TFile *_data = TFile::Open("respVsXnY_100GeV.root");
  //  TFile *_data = TFile::Open("respVsXnY.root");
  TFile *_mc = TFile::Open("~/BTFAnalysis/PositionAnalysis/OriginalSimulationData/simProfiles100.root");
  //  TFile *_mc = TFile::Open("~/BTFAnalysis/PositionAnalysis/OriginalSimulationData/simProfiles50.root");

   TProfile* mc_X=(TProfile*)_mc->Get("hx");
   TProfile* mc_Y=(TProfile*)_mc->Get("hy");

  TProfile* data_X=(TProfile*)_data->Get("hprofX_tot");
  TProfile* data_Y=(TProfile*)_data->Get("hprofY_tot");

  TProfile* data_X_rescaled=rescaleProfile(data_X);
  TProfile* data_Y_rescaled=rescaleProfile(data_Y);

  TProfile* mc_X_rescaled=rescaleProfile(mc_X);
  TProfile* mc_Y_rescaled=rescaleProfile(mc_Y);
  
  mc_X_rescaled->GetXaxis()->SetRangeUser(-12,12);
  data_X_rescaled->GetXaxis()->SetRangeUser(-12,12);
  mc_Y_rescaled->GetXaxis()->SetRangeUser(-12,12);
  data_Y_rescaled->GetXaxis()->SetRangeUser(-12,12);


  TCanvas* c1=new TCanvas("c1","c1",800,600);


  data_X_rescaled->GetYaxis()->SetTitleOffset(1.0);
  data_X_rescaled->GetXaxis()->SetTitleOffset(1.0);
  
  data_X_rescaled->Draw();
  data_X_rescaled->SetMarkerSize(1.2);
  data_X_rescaled->SetMarkerStyle(20);
  data_X_rescaled->SetMarkerColor(kBlue);
  data_X_rescaled->SetLineColor(kBlue);
  // data_X_rescaled->Fit("pol4","R0+","",-12,12);
  // data_X_rescaled->GetFunction("pol4")->SetLineWidth(2);
  // data_X_rescaled->GetFunction("pol4")->SetLineColor(kBlue);
  // data_X_rescaled->GetFunction("pol4")->Draw("SAME");

  mc_X_rescaled->Draw("SAME");
  mc_X_rescaled->SetMarkerSize(1.2);
  mc_X_rescaled->SetMarkerStyle(24);
  mc_X_rescaled->SetMarkerColor(kBlue);
  mc_X_rescaled->SetLineColor(kBlack);
  // mc_X_rescaled->Fit("pol4","R0+","",-12,12);
  // mc_X_rescaled->GetFunction("pol4")->SetLineWidth(2);
  // mc_X_rescaled->GetFunction("pol4")->SetLineColor(kBlack);
  // mc_X_rescaled->GetFunction("pol4")->Draw("SAME");




  TLegend* legX = new TLegend( 0.42, 0.45,0.75,  0.45 +0.06*2);
  legX->SetTextSize(0.04);
  legX->AddEntry(data_X_rescaled, "Data", "P");
  legX->AddEntry(mc_X_rescaled, "Simulation" , "P");
  legX->SetBorderSize(0);
  legX->SetFillColor(0);
  legX->Draw("same");

  TPaveText* lable_top;
  lable_top = DrawTools::getLabelTop("100 GeV Electron Beam");
  //  lable_top = DrawTools::getLabelTop("50 GeV Electron Beam");
  lable_top->Draw("same");

  TPaveText* label_low = new TPaveText(0.18 ,0.18,0.5,0.21, "brNDC");
  //  TPaveText* label_low = new TPaveText(0.165,0.175,0.5,0.21, "brNDC");
  label_low->SetFillColor(kWhite);
  label_low->SetTextSize(0.038);
  label_low->SetTextAlign(11); // align right
  label_low->SetTextFont(62);
  label_low->AddText( "W-CeF_{3} Single Tower");
  label_low->Draw("same");
 

  c1->SaveAs("dataMC_profile_X_corr_100GeV.png");
  c1->SaveAs("dataMC_profile_X_corr_100GeV.pdf");

  data_Y_rescaled->GetYaxis()->SetTitleOffset(1.0);
  data_Y_rescaled->GetXaxis()->SetTitleOffset(1.0);
 
  data_Y_rescaled->Draw();
  data_Y_rescaled->SetMarkerSize(1.2);
  data_Y_rescaled->SetMarkerStyle(20);
  data_Y_rescaled->SetMarkerColor(kBlue);
  data_Y_rescaled->SetLineColor(kBlue);
  mc_Y_rescaled->Draw("SAME");
  mc_Y_rescaled->SetMarkerSize(1.2);
  mc_Y_rescaled->SetMarkerStyle(24);
  mc_Y_rescaled->SetMarkerColor(kBlue);
  mc_Y_rescaled->SetLineColor(kBlack);

  TLegend* legY = new TLegend( 0.42, 0.45,0.75, 0.45 +0.06*2);
  legY->SetTextSize(0.04);
  legY->AddEntry(data_Y_rescaled, "Data", "P");
  legY->AddEntry(mc_Y_rescaled, "Simulation" , "P");
  legY->SetBorderSize(0);
  legY->SetFillColor(0);
  legY->Draw("same");
  lable_top->Draw("same");

 label_low->Draw("same");
 


  c1->SaveAs("dataMC_profile_Y_corr_100GeV.png");
  c1->SaveAs("dataMC_profile_Y_corr_100GeV.pdf");

}
