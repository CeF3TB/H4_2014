#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TProfile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "Math/Polynomial.h"
#include "../CommonTools/interface/DrawTools.h"
#include "../CommonTools/src/DrawTools.cc"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"

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

void fitFracResp()
{
  gStyle->SetOptTitle(0);
  DrawTools::setStyle();

  TFile *_file0 = TFile::Open("respVsR.root");

  TFile* out=TFile::Open("fitRespVsR.root","RECREATE");

  //px fits
  float max_x_values[4];
  float center_x_values[4];
  for(int i=0;i<4;++i)
    {
      TProfile* p=(TProfile*)_file0->Get(Form("hprofR_%d",i));
      TFitResultPtr fit=p->Fit("pol7","R+","",2,34);
      p->Write();
      
      /*
      TProfile*p1=rescaleProfile(p);
      p1->Write();
      max_x_values[i]=findMaximum(p->GetFunction("pol4"),4);
      center_x_values[i]=p->GetFunction("pol4")->Eval(0);

      //Let's draw some shit
      TCanvas* canny = new TCanvas("canny","canny",800,600);
      p->GetYaxis()->SetTitleOffset(1.0);
      p->GetXaxis()->SetTitleOffset(1.0);
      p->Draw();
      p->SetMarkerSize(1.2);
      p->SetMarkerStyle(20);
      p->SetMarkerColor(kBlue);
      p->SetLineColor(kBlue);

      fit->Draw("same");

      TLegend* leg = new TLegend( 0.4, 0.3, 0.75, 0.5);
      leg->SetTextSize(0.04);
      leg->AddEntry(p,Form("Ch %d, X",i), "P");
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->Draw("same");

      TPaveText* lable_top;
      lable_top = DrawTools::getLabelTop("50 GeV Electron Beam");
      lable_top->Draw("same");

      canny->SaveAs(Form("intercal_X_ch_%d.png",i));
      */
    }      //      

 
 out->Close();

  

}
