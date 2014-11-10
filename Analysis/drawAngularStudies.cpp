#include <stdio.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <stdlib.h> 

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"


#include "DrawTools.h"




struct ResoStruct {

  float resp;
  float resp_error;

  float reso;
  float reso_error;

};


struct LateralScanStruct {

 LateralScanStruct( float d, TTree* t_data, TTree* t_mc ) {
   offset = d;
   tree_mc = t_mc;
   tree_data = t_data;
 }

 float offset;
 TTree* tree_data;
 TTree* tree_mc;

};


ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name );
ResoStruct getResponseResolutionMCdata( const std::string& outputdir, TTree* tree, const std::string& name );
std::string getVarName( float LYSF[] );
ResoStruct getRespAndReso( TF1* f1, float energyErrorPercent );
float getRatioError( float num, float denom, float numErr, float denomErr );
ResoStruct getRespResoFromHisto( TH1D* h1 );
ResoStruct addPhotoStatistics( ResoStruct rs );
TGraph* shade(TCanvas *c1, TF1 *f1, TF1 *f2);

void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name, int niter, float nSigma );

void doSingleFitMC( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name, int niter, float nSigma );


int main( int argc, char* argv[]) {

  std::string tag = "V00";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }

  DrawTools::setStyle();
  gStyle->SetOptFit(1);

  std::string outputdir = "AngularReso";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );









  std::vector< float > angles;
  std::vector<std::string> anglesText; 
  angles.push_back(0.0);
  angles.push_back( 1.0 ); 
  angles.push_back( 2.0 );
  angles.push_back( 3. );
  angles.push_back( 4.0 );
  angles.push_back( 5.0 );
  angles.push_back( 6. );
  angles.push_back( 7. );
  angles.push_back( 8. );
  angles.push_back( 9. );


  anglesText.push_back("0.0");
  anglesText.push_back("1.0"); 
  anglesText.push_back("2.0"); 
  anglesText.push_back("3.0");
  anglesText.push_back("4.0");
  anglesText.push_back("5.0");
  anglesText.push_back( "6.0" );
  anglesText.push_back( "7.0" );
  anglesText.push_back( "8.0" );
  anglesText.push_back( "9.0" );


  std::vector< float > anglesData;
  std::vector<std::string> anglesTextData; 
  anglesData.push_back(0.0);
  anglesData.push_back( 2.5 );
  anglesData.push_back( 5.0 );
  anglesData.push_back( 7.5 );


  anglesTextData.push_back("497");
  anglesTextData.push_back("533"); 
  anglesTextData.push_back("523");
  anglesTextData.push_back( "528" );



  TGraphErrors* gr_reso = new TGraphErrors(0);
  TGraphErrors* gr_reso_ps = new TGraphErrors(0);

  TGraphErrors* gr_reso_data = new TGraphErrors(0);




  for( unsigned i=0; i< angles.size(); ++i ) {

    double intPart;
    double fracPart = modf( angles[i], &intPart );

    char angleText[100];
    sprintf( angleText, "%.0f.%.0f", floor(angles[i]), fracPart*10. );


      TFile* file = TFile::Open( Form("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/H4angleTest/Reco_Simulation%s.root", angleText ) );

    if( file==0 ) {
      std::cout << "WARNING: didn't find file for " << angleText << std::endl;
      continue;
    }
    TTree* tree = (TTree*)file->Get("recoTree");
    if( tree==0 ) {
      std::cout << "WARNING: didn't find tree in file: " << file->GetName() << std::endl;
      continue;
    }
     
    ResoStruct rs = getResponseResolutionMC( outputdir, tree, angleText );
    ResoStruct rs_ps = addPhotoStatistics( rs );

    gr_reso->SetPoint(i, angles[i], rs.reso);
    gr_reso->SetPointError(i, 0., rs.reso_error);

    gr_reso_ps->SetPoint(i, angles[i], rs_ps.reso);
    gr_reso_ps->SetPointError(i, 0., rs_ps.reso_error);

  }
 


  for( unsigned i=0; i< anglesData.size(); ++i ) {

    TFile* file = TFile::Open( Form("analysisTrees_%s/Reco_%s.root",tag.c_str(), anglesTextData[i].c_str() ) );

    if( file==0 ) {
      std::cout << "WARNING: didn't find file for " << anglesTextData[i] << std::endl;
      continue;
    }
    TTree* tree = (TTree*)file->Get("recoTree");
    if( tree==0 ) {
      std::cout << "WARNING: didn't find tree in file: " << file->GetName() << std::endl;
      continue;
    }
     
 
    std::string name = Form("%0f", anglesData[i]  );
    ResoStruct rs = getResponseResolutionMCdata( outputdir, tree, name.c_str() );

    gr_reso_data->SetPoint(i, anglesData[i], rs.reso);
    gr_reso_data->SetPointError(i, 0., rs.reso_error);

    std::cout << rs.reso << std::endl;
  }







  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, -0.5, 8.5, 10, 0., 85.);
  h2_axes->SetXTitle("Beam Impact Angle [deg]");
  h2_axes->SetYTitle("CeF_{3} Energy Resolution [%]");
  h2_axes->Draw("same");

  gr_reso_ps->SetMarkerStyle(20);
  gr_reso_ps->SetMarkerColor(38);
  gr_reso_ps->SetMarkerSize(1.4);
  gr_reso_ps->Draw("P same");

  TF1* fun= new TF1("fun", "[0]+[1]*x*x+[2]*x*x*x*x",0,8.5);
  gr_reso_ps->Fit("fun","RN");
  fun->SetLineWidth(1.5);
  fun->SetFillColor(38);
  fun->Draw("same");
 
  gr_reso_data->SetMarkerStyle(20);
  gr_reso_data->SetMarkerColor(46);
  gr_reso_data->SetMarkerSize(1.4);
  gr_reso_data->Draw("P same");


  //  TLegend* legend = new TLegend( 0.2, 0.75, 0.45, 0.9 );
  TLegend* legend = new TLegend( 0.2, 0.9-4.*0.06, 0.5, 0.9 );
  legend->SetTextSize(0.035);
  legend->SetFillColor(0);
  legend->AddEntry( gr_reso_data, "Data", "P" );
  legend->AddEntry( gr_reso_ps, "MC w/o stuff", "P" );
  legend->Draw("same");


  TPaveText* labelTop = DrawTools::getLabelTop("50 GeV Electron Beam");
  labelTop->Draw("same");

  c1->SaveAs(Form("%s/reso_vs_angle.pdf",outputdir.c_str()) );








  return 0;

}

/////////////////////////////////THE END.... or NOT//////////////////////




ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree,  const std::string& name ) {

  TH1D* h1 = new TH1D( name.c_str(), "", 200, 5., 20000. );

  tree->Project( name.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]");

  TF1* f1 = new TF1( Form("f1_%s", name.c_str()), "gaus", 5., 20000.);

  doSingleFitMC( h1, f1, outputdir, name,15.,1.0);


  ResoStruct rs = getRespAndReso( f1, 0. );

  return rs;

}


void doSingleFitMC( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name, int niter, float nSigma ) {

  h1->Fit( f1, "RQN" );

  for( unsigned iter=0; iter<niter; iter++ ) {

    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    float fitMin = mean - nSigma*sigma;
    float fitMax = mean + nSigma*sigma;
    f1->SetRange( fitMin, fitMax );
    if( iter==(niter-1) )
      h1->Fit( f1, "RQ+" );
    else
      h1->Fit( f1, "RQN" );
  }
  TCanvas* c1 = new TCanvas("cX", "", 600, 600);
  c1->cd();
  h1->Draw();
  c1->SaveAs( Form("%s/fit_%s.png", outputdir.c_str(), name.c_str()) );
  delete c1;
}





ResoStruct getRespAndReso( TF1* f1, float energyErrorPercent ) {

  float energy = 50000.;
  float energyErr = energyErrorPercent*energy;

  float mean = f1->GetParameter(1);
  float meanErr = f1->GetParError(1);

  float rms = f1->GetParameter(2);
  float rmsErr = f1->GetParError(2);

  float reso = 100.* rms/mean; //in percent
  float resoErr = 100.* getRatioError( rms, mean, meanErr, rmsErr );

  ResoStruct rs;
  rs.resp = mean;
  rs.resp_error = meanErr;
  rs.reso = reso;
  rs.reso_error = resoErr;


  return rs;

}


float getRatioError( float num, float denom, float numErr, float denomErr ) {

  return sqrt( numErr*numErr/(denom*denom) + denomErr*denomErr*num*num/(denom*denom*denom*denom) );

}



ResoStruct addPhotoStatistics( ResoStruct rs ) {

  // MC response is already in MeV, 0.49 is the number of p.e. per MeV
  // RM = 90% of energy, so assume 90% of energy (0.9*491) is deposited in central channel
  //  float nADC = rs.resp/185.*3200.;
  //float  nADC = rs.resp/177.*3224.36;
  //float  nADC = rs.resp/169.445*3224.36;
   float  nADC = rs.resp/168.1*3224.36;
  float nPhotoElectrons = nADC/27.3;
  //  float nPhotoElectrons = nADC/35.;

  float poissonError = 100./sqrt( nPhotoElectrons ); // in percent

  float resoUnsmeared = rs.reso;
  
  float resoSmeared =  sqrt( rs.reso*rs.reso + poissonError*poissonError );
  
  rs.reso = resoSmeared;

  rs.reso_error = rs.reso_error * resoSmeared / resoUnsmeared; 

  return rs;

} 






ResoStruct getResponseResolutionMCdata( const std::string& outputdir, TTree* tree,  const std::string& name ) {

  TH1D* h1 = new TH1D( name.c_str(), "",400, 0., 400. );

  tree->Project( name.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]");
 
  TF1* f1 = new TF1( Form("f1_%s", name.c_str()), "gaus", 10., 800.);

  doSingleFit( h1, f1, outputdir, name, 10., 1.4 );

  ResoStruct rs = getRespAndReso( f1, 0. );

  return rs;
}



void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name, int niter, float nSigma ) {

  h1->Fit( f1, "RQN" );

  for( unsigned iter=0; iter<niter; iter++ ) {

    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    float fitMin = mean - nSigma*sigma*0.6;
    float fitMax = mean + nSigma*sigma;
    f1->SetRange( fitMin, fitMax );
    if( iter==(niter-1) )
  
 h1->Fit( f1, "RQ+" );
    else
         h1->Fit( f1, "RQN" );
  }


  TCanvas* c1 = new TCanvas("cX", "", 600, 600);
  c1->cd();

  h1->Draw();

  c1->SaveAs( Form("%s/fit_%s.png", outputdir.c_str(), name.c_str()) );

  delete c1;

}
