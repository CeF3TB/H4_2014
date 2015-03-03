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


#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooDataHist.h"
#include "RooCBShape.h"

//#include "TAxis.h"
#include "RooPlot.h"

#include "RooGlobalFunc.h"


#include "DrawTools.h"

using namespace RooFit;




struct FitStruct {

  double mean;
  double mean_err;

  double sigma;
  double sigma_err;

  double reso;
  double reso_err;

};



struct ResoStruct {

  float resp;
  float resp_error;

  float reso;
  float reso_error;

};

FitStruct getFitResults( const std::string& outputdir, TTree* tree, const std::string& name, float xMin, float xMax, std::string& cut, std::string& whatToProject );
FitStruct addPhotoStatistics( FitStruct rs );
float getRatioError( float num, float denom, float numErr, float denomErr );


ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name, float energySim );
ResoStruct getRespAndReso( TF1* f1, float energySim );
ResoStruct addPhotoStatistics( ResoStruct rs );

void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name );
TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree );



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




    std::string whatToProjectSim = Form("cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]");

    //  std::string whatToProject = Form("(cef3_chaInt_corr[0]+cef3_chaInt_corr[1]+cef3_chaInt_corr[2]+cef3_chaInt_corr[3])");

    std::string whatToProject = Form("(cef3_maxAmpl_corr[0]+cef3_maxAmpl_corr[1]+cef3_maxAmpl_corr[2]+cef3_maxAmpl_corr[3])");





  std::vector< float > angles;
  std::vector<std::string> anglesText; 
  angles.push_back(0.0);
  angles.push_back( 1.0 ); 
  angles.push_back( 2.0 );
  angles.push_back( 2.5 );
  angles.push_back( 3. );
  angles.push_back( 4.0 );
  angles.push_back( 5.0 );
  angles.push_back( 6. );
  angles.push_back( 7. );
  angles.push_back( 7.5);
  angles.push_back( 8. );
  angles.push_back( 9. );
  angles.push_back( 10. );
  angles.push_back( 15. );


  anglesText.push_back("0.0");
  anglesText.push_back("1.0"); 
  anglesText.push_back("2.0"); 
  anglesText.push_back("2.5"); 
  anglesText.push_back("3.0");
  anglesText.push_back("4.0");
  anglesText.push_back("5.0");
  anglesText.push_back( "6.0" );
  anglesText.push_back( "7.0" );
  anglesText.push_back( "7.5" );
  anglesText.push_back( "8.0" );
  anglesText.push_back( "9.0" );
  anglesText.push_back( "10.0" );
  anglesText.push_back( "15.0" );


  std::vector< float > anglesData;
  std::vector<std::string> anglesTextData; 
  anglesData.push_back(-0.1);
  anglesData.push_back(0.0);
  anglesData.push_back( 2.5 );
  anglesData.push_back( 5.0 );
  anglesData.push_back( 7.5 );


  anglesTextData.push_back("323");
  anglesTextData.push_back("497");
  anglesTextData.push_back("533"); 
  anglesTextData.push_back("523");
  anglesTextData.push_back( "528" );


  //  std::string SimCut = "";
  std::string SimCut = "abs(xPos)< 3 && abs(yPos)< 3";

  //Correction: Subract difference between 600V & 950V in quadrature

  TFile* file950 = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_323.root", tag.c_str() ) );
  TTree* tree950 = (TTree*)file950->Get("recoTree");
  std::string cut = "abs(0.5* (position_X1+position_X2))< 3 && abs( 0.5* (position_Y1+position_Y2))< 3 && (wc_x_corr-cluster_pos_corr_hodoX2)<4 && (wc_y_corr-cluster_pos_corr_hodoY2)< 4  && abs( (position_X1-position_X2))<1.5 &&abs( (position_Y1-position_Y2))<1.5  ";

  FitStruct fs950 = getFitResults( outputdir, tree950,  "950HV" ,80 , 5000, cut, whatToProject );
  float reso_950 = fs950.reso;
  float reso_err_950 = fs950.reso_err;


  TFile* file600 = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_497.root", tag.c_str() ));
  TTree* tree600 = (TTree*)file600->Get("recoTree");
  FitStruct fs600 = getFitResults( outputdir, tree600,  "600HV" , 80 , 5000, cut, whatToProject );
  float reso_600 = fs600.reso;
  float reso_err_600 = fs600.reso_err;


  float corr = sqrt( reso_600*reso_600 - reso_950*reso_950);

  std::cout << " At 600V = " << reso_600 << std::endl;
  std::cout << " At 950V = " << reso_950 << std::endl;
  std::cout << " Correction = " << corr << std::endl;

  TGraphErrors* gr_reso = new TGraphErrors(0);
  TGraphErrors* gr_reso_ps = new TGraphErrors(0);

  TGraphErrors* gr_reso_data = new TGraphErrors(0);
  TGraphErrors* gr_reso_data_corr = new TGraphErrors(0);

  TGraphErrors* gr_reso_dev = new TGraphErrors(0);
  TGraphErrors* gr_reso_dev_data = new TGraphErrors(0);
  TGraphErrors* gr_reso_dev_data_corr = new TGraphErrors(0);



  TGraphErrors* gr_reso_g = new TGraphErrors(0);
  TGraphErrors* gr_reso_ps_g = new TGraphErrors(0);

  TGraphErrors* gr_reso_data_g = new TGraphErrors(0);
  TGraphErrors* gr_reso_data_corr_g = new TGraphErrors(0);

  TGraphErrors* gr_reso_dev_g = new TGraphErrors(0);
  TGraphErrors* gr_reso_dev_data_g = new TGraphErrors(0);
  TGraphErrors* gr_reso_dev_data_corr_g = new TGraphErrors(0);



  double zeroValue_sim;
  double zeroValue_data;
  double zeroValue_data_corr;

  double zeroValue_sim_g;
  double zeroValue_data_g;
  double zeroValue_data_corr_g;

  for( unsigned i=0; i< angles.size(); ++i ) {

    double intPart;
    double fracPart = modf( angles[i], &intPart );
    char angleText[100];
    sprintf( angleText, "%.0f.%.0f", floor(angles[i]), fracPart*10. );

    TFile* file = TFile::Open( Form("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/H4AngularUseThis/Reco_Simulation%s.root", angleText ) );

    if( file==0 ) {
      std::cout << "WARNING: didn't find file for " << angleText << std::endl;
      continue;
    }
    TTree* tree = (TTree*)file->Get("recoTree");
    if( tree==0 ) {
      std::cout << "WARNING: didn't find tree in file: " << file->GetName() << std::endl;
      continue;
    }

    FitStruct fsSim= getFitResults( outputdir, tree , angleText , 1000 , 20000, SimCut , whatToProjectSim );
    float energySim = fsSim.mean;
    FitStruct fsSim_ps = addPhotoStatistics( fsSim );
     
    //    ResoStruct rs = getResponseResolutionMC( outputdir, tree, angleText );
    //    ResoStruct rs_ps = addPhotoStatistics( rs );

    gr_reso->SetPoint(i, angles[i], fsSim.reso);
    gr_reso->SetPointError(i, 0., fsSim.reso_err);

    gr_reso_ps->SetPoint(i, angles[i], fsSim_ps.reso);
    gr_reso_ps->SetPointError(i, 0., fsSim_ps.reso_err);

    if(i==0){
      zeroValue_sim = fsSim.reso;
    }
    gr_reso_dev->SetPoint(i,angles[i], fsSim.reso / zeroValue_sim);



    //and now the same with an iterative gaussian...


    ResoStruct rs = getResponseResolutionMC( outputdir, tree, angleText, angles[i] );
    ResoStruct rs_ps = addPhotoStatistics( rs );

    gr_reso_g->SetPoint( i, angles[i], rs.reso );
    gr_reso_g->SetPointError( i,0, rs.reso_error);
 
    gr_reso_ps_g->SetPoint( i, angles[i], rs_ps.reso );
    gr_reso_ps_g->SetPointError( i,0,  rs_ps.reso_error );
    
    if(i==0){
      zeroValue_sim_g = rs.reso;
    }
    gr_reso_dev_g->SetPoint(i,angles[i], rs.reso / zeroValue_sim_g);



  }//End simulation loop
 


  for( unsigned i=0; i< anglesData.size(); ++i ) {

    TFile* file = TFile::Open( Form("analysisTrees_chaInt_04_12_%s/Reco_%s.root",tag.c_str(), anglesTextData[i].c_str() ) );

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

    //   ResoStruct rs = getResponseResolutionMCdata( outputdir, tree, name.c_str() );

    FitStruct fs  = getFitResults( outputdir, tree ,  name.c_str() , 80 , 5000  , cut, whatToProject );


    gr_reso_data->SetPoint(i, anglesData[i], fs.reso);
    gr_reso_data->SetPointError(i, 0., fs.reso_err);
    
    //Error of corrected resolution (funfunfun)
    float corr_error = sqrt( fs.reso_err *fs.reso_err * 2* fs.reso* 2* fs.reso /(fs.reso*fs.reso - corr*corr)  + reso_err_950 *reso_err_950  * 2* reso_950* 2* reso_950 /(fs.reso*fs.reso - corr*corr)  + reso_err_600 *reso_err_600  * 2* reso_600* 2* reso_600 /(fs.reso*fs.reso - corr*corr) );

    gr_reso_data_corr->SetPoint(i, anglesData[i], sqrt(fs.reso*fs.reso - corr*corr));
    gr_reso_data_corr->SetPointError(i, 0., corr_error);
    
    if(i==1){
      zeroValue_data = fs.reso;
      zeroValue_data_corr = sqrt(fs.reso*fs.reso - corr*corr);
    }
    gr_reso_dev_data->SetPoint(i,anglesData[i], fs.reso / zeroValue_data);
    gr_reso_dev_data_corr->SetPoint(i,anglesData[i], sqrt(fs.reso*fs.reso - corr*corr) / zeroValue_data_corr);



    //Samesame but different: Iterative Gaussian Fit
    TF1* thisFunc = fitSingleElectronPeak( outputdir,  name , tree );
    
    float mean_g = thisFunc->GetParameter(1);
    float meanErr_g = thisFunc->GetParError(1);

    float rms_g = thisFunc->GetParameter(2);
    float rmsErr_g = thisFunc->GetParError(2);

    float reso_g = 100.* rms_g/mean_g ;
    float resoErr_g = 100* getRatioError( rms_g, mean_g, rmsErr_g, meanErr_g);

    if(i==1){
      zeroValue_data_g = reso_g;
      zeroValue_data_corr_g = sqrt(reso_g*reso_g - corr*corr);
    }

    gr_reso_data_g->SetPoint( i, anglesData[i] , reso_g );
    gr_reso_data_g->SetPointError( i, 0, resoErr_g );

    gr_reso_dev_data_g->SetPoint(i, anglesData[i], reso_g / zeroValue_data_g );

    std::cout << std::endl;
    std::cout << fs.reso << std::endl;
    std::cout << std::endl;

    std::cout << reso_g << std::endl;
  }







  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, -0.5, 8.5, 10, 0., 16.);
  h2_axes->SetXTitle("Beam Impact Angle [deg]");
  h2_axes->SetYTitle("CeF_{3} Energy Resolution [%]");
  h2_axes->Draw("same");

  gr_reso_ps->SetMarkerStyle(21);
  gr_reso_ps->SetMarkerColor(kBlack);
  gr_reso_ps->SetMarkerSize(1.4);
  gr_reso_ps->Draw("P same");

  gr_reso_ps_g->SetMarkerStyle(22);
  gr_reso_ps_g->SetMarkerColor(kBlack);
  gr_reso_ps_g->SetMarkerSize(1.4);
  gr_reso_ps_g->Draw("P same");

  // TF1* fun= new TF1("fun", "[0]+[1]*x*x+[2]*x*x*x*x",0,8);
  TF1* fun= new TF1("fun", "[0]+[1]*x*x+[2]*x*x*x*x+[3]*x*x*x*x*x*x",0,8);
  gr_reso_ps->Fit("fun","RN");
  fun->SetLineWidth(1.5);
  //  fun->Draw("same");
 
  gr_reso_data->SetMarkerStyle(20);
  gr_reso_data->SetMarkerColor(46);
  gr_reso_data->SetMarkerSize(1.6);
  gr_reso_data->Draw("P same");
  
  gr_reso_data_g->SetMarkerStyle(22);
  gr_reso_data_g->SetMarkerColor(42);
  gr_reso_data_g->SetMarkerSize(1.6);
  gr_reso_data_g->Draw("P same");
  
  //Corrected data
  gr_reso_data_corr->SetMarkerStyle(24);
  gr_reso_data_corr->SetMarkerColor(38);
  gr_reso_data_corr->SetMarkerSize(1.6);
  gr_reso_data_corr->Draw("P same");

  //  TLegend* legend = new TLegend( 0.2, 0.75, 0.45, 0.9 );
  TLegend* legend = new TLegend( 0.2, 0.9-4.*0.06, 0.5, 0.9 );
  legend->SetTextSize(0.035);
  legend->SetFillColor(0);
  legend->AddEntry( gr_reso_data, "Data CB", "P" );
  legend->AddEntry( gr_reso_data_g, "Data it.Gaussian", "P" );
  legend->AddEntry( gr_reso_data_corr, "Data CB Corrected", "P" );
  legend->AddEntry( gr_reso_ps, "MC", "P" );
  legend->AddEntry( gr_reso_ps_g, "MC it.Gaussian", "P" );
  legend->Draw("same");

  TPaveText* labelTop = DrawTools::getLabelTop("50 GeV Electron Beam");
  labelTop->Draw("same");

  c1->SaveAs(Form("%s/reso_vs_angle.eps",outputdir.c_str()) );
  c1->SaveAs(Form("%s/reso_vs_angle.png",outputdir.c_str()) );







  c1->Clear();
  c1->cd();
  //Deviation from 0° measurement:

  TH2D* h2_axes2 = new TH2D("axes2", "", 10, -0.5, 8.5, 10, 0., 10.);
  h2_axes2->SetXTitle("Beam Impact Angle [deg]");
  h2_axes2->SetYTitle("Normalized CeF_{3} Energy Resolution");
  h2_axes2->Draw("");

  gr_reso_dev->SetMarkerStyle(21);
  gr_reso_dev->SetMarkerColor(kBlack);
  gr_reso_dev->SetMarkerSize(1.4);
  gr_reso_dev->Draw("P same");

  gr_reso_dev_g->SetMarkerStyle(22);
  gr_reso_dev_g->SetMarkerColor(kBlack);
  gr_reso_dev_g->SetMarkerSize(1.4);
  gr_reso_dev_g->Draw("P same");

  TF1* fun_dev= new TF1("fun_dev", "[0]+[1]*x*x+[2]*x*x*x*x",0,8);
  //TF1* fun_dev= new TF1("fun_dev", "[0]+[1]*x*x+[2]*x*x*x*x+[3]*x*x*x*x*x*x",0,8);
  gr_reso_dev->Fit("fun_dev","RN");
  fun_dev->SetLineWidth(1.5);
  fun_dev->SetFillColor(38);
  // fun_dev->Draw("same");
 
 

  
  gr_reso_dev_data->SetMarkerStyle(20);
  gr_reso_dev_data->SetMarkerColor(46);
  gr_reso_dev_data->SetMarkerSize(1.6);
  gr_reso_dev_data->Draw("P same");
  
  //Corrected data
  gr_reso_dev_data_corr->SetMarkerStyle(24);
  gr_reso_dev_data_corr->SetMarkerColor(38);
  gr_reso_dev_data_corr->SetMarkerSize(1.6);
  gr_reso_dev_data_corr->Draw("P same");



  gr_reso_dev_data_g->SetMarkerStyle(22);
  gr_reso_dev_data_g->SetMarkerColor(42);
  gr_reso_dev_data_g->SetMarkerSize(1.6);
  gr_reso_dev_data_g->Draw("P same");

  //  TLegend2* legend2 = new TLegend2( 0.2, 0.75, 0.45, 0.9 );
  TLegend* legend2 = new TLegend( 0.2, 0.9-4.*0.06, 0.5, 0.9 );
  legend2->SetTextSize(0.035);
  legend2->SetFillColor(0);
  legend2->AddEntry( gr_reso_dev_data, "Data", "P" );
  legend2->AddEntry( gr_reso_dev_data_g, "Data it. Gaussian", "P" );
  legend2->AddEntry( gr_reso_dev_data_corr, "Data Corrected", "P" );
  legend2->AddEntry( gr_reso_ps, "MC", "P" );
  legend2->AddEntry( gr_reso_ps_g, "MC it.Gaussian", "P" );
  legend2->Draw("same");

  
  labelTop->Draw("same");

  c1->SaveAs(Form("%s/reso_dev_vs_angle.eps",outputdir.c_str()) );
  c1->SaveAs(Form("%s/reso_dev_vs_angle.png",outputdir.c_str()) );




  std::cout << " At 600V = " << reso_600 << std::endl;
  std::cout << " At 950V = " << reso_950 << std::endl;
  std::cout << " Correction = " << corr << std::endl;




  return 0;

}

/////////////////////////////////THE END.... or NOT//////////////////////




































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

  
  double peakpos = f1->GetParameter(1);
  double sigma = f1->GetParameter(2);
  double fitmin;
  double fitmax;
  
 
  if( (peakpos-5*sigma) < 35 ){ 
    fitmin = 35;
  } else{ 
    fitmin = peakpos-5*sigma;
  }

  if(((peakpos>0) && (sigma > 2*peakpos)) ) { 
    fitmax = peakpos*2;
  } else{ 
    fitmax = peakpos+4*sigma;
  }
   
  // fitmax = peakpos+4*sigma;

  if(fitmin <50) fitmin =50;
  // if (sigma < 10) sigma = 10;


  TH1D* h2 = new TH1D("h2", "", 200, fitmin,fitmax);
  tree->Project( "h2", whatToProject.c_str(), cut.c_str() );

  RooRealVar x("x","ADC Channels", fitmin, fitmax);
  RooDataHist data("data","dataset with x",x,Import(*h2) );
  
  RooPlot* frame;
  RooPlot* xframe = x.frame();   
  
  frame = x.frame("Title");
  data.plotOn(frame);  //this will show histogram data points on canvas
  //  data.statOn(frame);  //this will display hist stat on canvas
  
  RooRealVar meanr("meanr","Mean",peakpos,peakpos-sigma, peakpos+1.5*sigma);
  RooRealVar width("width","#sigma",sigma, 20.0, 5.*sigma);
  RooRealVar A("A","Dist",1.1, 0.0, 15.0);
  RooRealVar N("N","Deg",3, 0.0, 150);

  meanr.setRange(100. , 25000.);
  width.setRange( 5, 11000);
  
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

  TLegend* lego = new TLegend(0.65, 0.7, 0.9, 0.92);  
  lego->SetTextSize(0.038);
  lego->AddEntry(  (TObject*)0 ,Form("#mu = %.0f #pm %.0f", meanr.getVal(), meanr.getError() ), "");
  lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.0f #pm %.0f ", width.getVal(), width.getError() ), "");
  lego->AddEntry(  (TObject*)0 ,Form("#chi^{2} = %.2f / %d ", frame->chiSquare(ndf) , ndf ), "");
  lego->AddEntry(  (TObject*)0 ,Form("#sigma/#mu = %.2f #pm  %.2f ", reso , resoErr ), "");


  lego->AddEntry(  (TObject*)0 ,Form("A = %.2f #pm  %.2f ", A.getVal(), A.getError() ), "");
  lego->AddEntry(  (TObject*)0 ,Form("N = %.2f #pm  %.2f ", N.getVal(), N.getError() ), "");

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
  //   float nADC = rs.resp/169.445*3224.36;
  //  float nADC = rs.resp/173.097*3235.76;
  // float  nADC = rs.mean/171.6*3236;
  float  nADC = rs.mean/169.445*3224.36; //without hodo cut
   float nPhotoElectrons = nADC/27.3;
  //float nADC = rs.resp/185.*3200.;
  //float nPhotoElectrons = nADC/35.;

  float poissonError = 100./sqrt( nPhotoElectrons ); // in percent
 
  float resoUnsmeared = rs.reso;
  float resoSmeared =  sqrt( rs.reso*rs.reso + poissonError*poissonError );

  rs.reso = resoSmeared;
  rs.reso_err = rs.reso_err * resoSmeared / resoUnsmeared; 

  return rs;

} 























ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name, float energySim ) {

  gStyle->SetOptFit(1);

  TH1D* h1 = new TH1D( name.c_str(), "", 500, 1000, 23000 );
  tree->Project( name.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]","abs(xPos)<3 && abs(yPos)<3" );


 TF1* f1 = new TF1( Form("f1_%s", name.c_str()), "gaus",  500, 22000 );

 f1->SetLineColor(kRed);
 doSingleFit( h1, f1, outputdir, name );
  ResoStruct rs = getRespAndReso( f1, energySim );



  return rs;

}



ResoStruct getRespAndReso( TF1* f1, float energySim) {

  // float energy = energySim;
  float energyErrS = 0;

  float meanS = f1->GetParameter(1);
  float meanErrS = f1->GetParError(1);

  float rmsS = f1->GetParameter(2);
  float rmsErrS = f1->GetParError(2);

  // float resoS = 100.* rmsS/energySim; //in percent
  float resoS = 100.* rmsS/meanS; //in percent
  float resoErrS = 100.* getRatioError( rmsS, meanS, meanErrS, rmsErrS );

  ResoStruct rs;
  rs.resp = meanS;
  rs.resp_error = meanErrS;
  rs.reso = resoS;
  rs.reso_error = resoErrS;

  return rs;

}


ResoStruct addPhotoStatistics( ResoStruct rs ) {

  // MC response is already in MeV, 0.49 is the number of p.e. per MeV
  // RM = 90% of energy, so assume 90% of energy (0.9*491) is deposited in central channel
  //   float nADC = rs.resp/169.445*3224.36;
  //  float nADC = rs.resp/173.097*3235.76;
  float  nADC = rs.resp/171.6*3236;
  //  float  nADC = rs.resp/169.445*3224.36;
   float nPhotoElectrons = nADC/27.3;
  //float nADC = rs.resp/185.*3200.;
  //float nPhotoElectrons = nADC/35.;

  float poissonError = 100./sqrt( nPhotoElectrons ); // in percent
 
  float resoUnsmeared = rs.reso;
  float resoSmeared =  sqrt( rs.reso*rs.reso + poissonError*poissonError );

  rs.reso = resoSmeared;
  rs.reso_error = rs.reso_error * resoSmeared / resoUnsmeared; 

  return rs;

} 






TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree ) {

  std::string histoName(Form("h1_%s", name.c_str() ));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 500, 0, 7500.);

    std::string whatToProject = Form("(cef3_maxAmpl_corr[0]+cef3_maxAmpl_corr[1]+cef3_maxAmpl_corr[2]+cef3_maxAmpl_corr[3])");
   std::string cut = "abs(0.5* (position_X1+position_X2))< 3 && abs( 0.5* (position_Y1+position_Y2))< 3 && (wc_x_corr-cluster_pos_corr_hodoX2)<4 && (wc_y_corr-cluster_pos_corr_hodoY2)< 4  && abs( (position_X1-position_X2))<1.5 &&abs( (position_Y1-position_Y2))<1.5  ";

   tree->Project( histoName.c_str(), whatToProject.c_str() , cut.c_str() );

   float peakpos = h1->GetMean();
   float width = h1->GetRMS();

  std::string histoName2(Form("h2_%s", name.c_str() ));
  TH1D* h2 = new TH1D(histoName2.c_str(), "", 500, peakpos-3*width, peakpos+3*width);

   tree->Project( histoName2.c_str(), whatToProject.c_str() , cut.c_str() );


  TF1* f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 20., 20000.);
  //  TF1* f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 500., 1800000.);
  //  TF1* f1 = new TF1( Form("gaus_%d", i), "gaus", 400., 1200.);
      f1->SetParameter(0, 3000.);
   f1->SetParameter(1, 200.);
   f1->SetParameter(2, 10.);

   f1->SetLineColor(kRed);

   doSingleFit( h2, f1, outputdir, name );

  return f1;

}







void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name ) {

  h1->Fit( f1, "RQN" );

  int niter = 4.;
  float nSigma =1 ;

  for( unsigned iter=0; iter<niter; iter++ ) {

    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    float fitMin = mean - nSigma*sigma*0.7;
    float fitMax = mean + nSigma*sigma * 3;
    f1->SetRange( fitMin, fitMax );
    if( iter==(niter-1) )
      h1->Fit( f1, "RQN" );
    else
      h1->Fit( f1, "RQ+" );
  }


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h1->Draw();

    c1->SaveAs( Form("%s/fit_%s.png", outputdir.c_str(), name.c_str()) );
  delete c1;

}
