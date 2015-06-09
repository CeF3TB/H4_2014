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
#include "TGraphErrors.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"


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


FitStruct getFitResults( const std::string& outputdir, TTree* tree, const std::string& name, float xMin, float xMax, std::string& cut, std::string& whatToProject );


float getRatioError( float num, float denom, float numErr, float denomErr );


FitStruct addPhotoStatistics( FitStruct rs );


int main( int argc, char* argv[] ) {


  std::string tag="V01";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }
  std::cout << "-> Using tag: " << tag << std::endl;

  std::string outputdir = "HVEnergyScan/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  DrawTools::setStyle();
  TGaxis::SetMaxDigits(3);


    std::string whatToProjectSim = Form("cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]");
    std::string whatToProject = Form("cef3_chaInt_bare_corr[0]+cef3_chaInt_bare_corr[1]+cef3_chaInt_bare_corr[2]+cef3_chaInt_bare_corr[3]");
   std::string whatToProjectCorr1 = Form("cef3_chaInt_corr1[0]+cef3_chaInt_corr1[1]+cef3_chaInt_corr1[2]+cef3_chaInt_corr1[3]"); 
   std::string whatToProjectCorr2 = Form("cef3_chaInt_corr2[0]+cef3_chaInt_corr2[1]+cef3_chaInt_corr2[2]+cef3_chaInt_corr2[3]");
   std::string whatToProjectMaxAmpCorr2 = Form("cef3_maxAmp_corr2[0]+cef3_maxAmp_corr2[1]+cef3_maxAmp_corr2[2]+cef3_maxAmp_corr2[3]");


 
  std::vector<std::string> simulation; 
  std::vector<float> beamEnergySimulation; 

  //////////////////DATA/////////////

  /////////DATA taken at 950V on CeF3//////////
  std::vector<std::string> runs950; 
  std::vector<float> beamEnergy950;
  runs950.push_back("323"); //274 a previous one
  beamEnergy950.push_back(50000);

  runs950.push_back("329");
  beamEnergy950.push_back(100000); 

  runs950.push_back("363");
  beamEnergy950.push_back(20000);

  runs950.push_back("356");  //  runs.push_back("392");
  beamEnergy950.push_back(150000); 
  runs950.push_back("377");   //runs950.push_back("387");
  beamEnergy950.push_back(200000);

  TFile* file950 = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%s.root", tag.c_str(), runs950[0].c_str()));
  TTree* tree950 = (TTree*)file950->Get("recoTree");
  // std::string cut = "";
  std::string cut = "abs(0.5* (position_X1+position_X2))< 3 && abs( 0.5* (position_Y1+position_Y2))< 3 && (wc_x_corr-cluster_pos_corr_hodoX2)<4 && (wc_y_corr-cluster_pos_corr_hodoY2)< 4  && abs( (position_X1-position_X2))<1.5 &&abs( (position_Y1-position_Y2))<1.5  ";

  FitStruct fs950 = getFitResults( outputdir, tree950,  "950HV" , 50000,  900000, cut, whatToProject );
  float adcEnergyC950 = fs950.mean;
  float energyC950= beamEnergy950[0];
 

  /////////DATA at 700 Volts ////////////
  std::vector<std::string> runs700; 
  std::vector<float> beamEnergy700;
  runs700.push_back("398");
  beamEnergy700.push_back(50000);

  runs700.push_back("418");;
  beamEnergy700.push_back(100000);
  runs700.push_back("431");
  beamEnergy700.push_back(10000);
  runs700.push_back("457");
  beamEnergy700.push_back(15000);
  runs700.push_back("455");
  beamEnergy700.push_back(20000);
  runs700.push_back("407");
  beamEnergy700.push_back(150000); 
  runs700.push_back("428");
  beamEnergy700.push_back(200000);

  TFile* file700 = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%s.root", tag.c_str(),runs700[0].c_str()));
  TTree* tree700 = (TTree*)file700->Get("recoTree");

  FitStruct fs700 = getFitResults( outputdir, tree700,  "700HV" , 500,  50000000, cut, whatToProject );
  float adcEnergyC700 = fs700.mean;
  float energyC700= beamEnergy700[0];
 
  ///////////////////// 600 V////////////////////////
  std::vector<std::string> runs600; 
  std::vector<float> beamEnergy600;
  runs600.push_back("497");
  beamEnergy600.push_back(50000);

  runs600.push_back("503");;
  beamEnergy600.push_back(100000);

  runs600.push_back("513");
  beamEnergy600.push_back(15000);
  runs600.push_back("508");
  beamEnergy600.push_back(20000);

  runs600.push_back("481");
  beamEnergy600.push_back(150000); 
  runs600.push_back("487");
  beamEnergy600.push_back(200000);

  TFile* file600 = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%s.root", tag.c_str(),runs600[0].c_str()));
  TTree* tree600 = (TTree*)file600->Get("recoTree");

  FitStruct fs600 = getFitResults( outputdir, tree600,  "600HV" , 500,  50000, cut, whatToProject );
  float adcEnergyC600 = fs600.mean;
  float energyC600= beamEnergy600[0];
 



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
  
  TFile* energyfileS = TFile::Open("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/H4UseThis1_5mmGap10mmBeam/Reco_Simulation50.root" );
  TTree* energytreeS = (TTree*)energyfileS->Get("recoTree");
  
    std::string SimCut = "abs(xPos)< 3 && abs(yPos)< 3";
    // std::string SimCut = "";
  FitStruct fsSim= getFitResults( outputdir, energytreeS , simulation[0] ,  5000. ,  20000., SimCut , whatToProjectSim );
  float energyS = fsSim.mean;


  //  ResoStruct energyrsS = getResponseResolutionMC( outputdir, energytreeS,simulation[0],beamEnergySimulation[0]);
  
  //float energyS = energyrsS.resp;




  float xMax = beamEnergy950[beamEnergy950.size()-1] +5000;

  TGraphErrors* gr_resp_vs_energy = new TGraphErrors(0);
  TGraphErrors* gr_reso_vs_energy = new TGraphErrors(0);
  TGraphErrors* gr_dev = new TGraphErrors(0);
 
  TGraphErrors* gr_resp_vs_energy_simul = new TGraphErrors(0);
  TGraphErrors* gr_reso_vs_energy_simul = new TGraphErrors(0);
  TGraphErrors* gr_dev_simul = new TGraphErrors(0);

  TGraphErrors* gr_resp_vs_energy700 = new TGraphErrors(0);
  TGraphErrors* gr_reso_vs_energy700 = new TGraphErrors(0);
  TGraphErrors* gr_dev700 = new TGraphErrors(0);

  TGraphErrors* gr_resp_vs_energy600 = new TGraphErrors(0);
  TGraphErrors* gr_reso_vs_energy600 = new TGraphErrors(0);
  TGraphErrors* gr_dev600 = new TGraphErrors(0);


 
  float yMax = 0;

  for( unsigned i=0; i<runs950.size(); ++i ) {
    
    TFile* file = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%s.root", tag.c_str(),  runs950[i].c_str() ));
    TTree* tree = (TTree*)file->Get("recoTree");

    //    TF1* thisFunc = fitSingleElectronPeak( outputdir, runs950[i], tree, 950. );
    
    FitStruct fs  = getFitResults( outputdir, tree ,  runs950[i].c_str() , 20000, 19000000   , cut, whatToProject );

    float energy = beamEnergy950[i];
    float energyErr = 0.005*energy;
 
    float mean = fs.mean;
    float meanErr = fs.mean_err;

    float reso = fs.reso;
    float resoErr = fs.reso_err;

    std::cout << "mean at "<< energy << " = " << mean << std::endl;
    gr_resp_vs_energy->SetPoint( i, energy/1000., mean /adcEnergyC950 * adcEnergyC600 );
    gr_resp_vs_energy->SetPointError( i, energyErr/1000., meanErr);

    gr_reso_vs_energy->SetPoint( i, energy/1000., reso );
    gr_reso_vs_energy->SetPointError( i, energyErr/1000., resoErr );

    gr_dev->SetPoint(i, energy/1000. , mean / energy / ( adcEnergyC950 / energyC950) );
    gr_dev->SetPointError(i, 0, meanErr / energy / ( adcEnergyC950 / energyC950) );

 }


  for( unsigned i=0; i<runs700.size(); ++i ) {
    
    TFile* file = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%s.root", tag.c_str(),  runs700[i].c_str()));
    TTree* tree = (TTree*)file->Get("recoTree");
    
    //  TF1* thisFunc = fitSingleElectronPeak( outputdir, runs700[i], tree , 700.);
  
    FitStruct fs  = getFitResults( outputdir, tree ,  runs700[i].c_str() , 100, 19000000   , cut, whatToProject );

    float energy = beamEnergy700[i];
    float energyErr = 0.005*energy;
 
    float mean = fs.mean;
    float meanErr = fs.mean_err;

    float reso = fs.reso;
    float resoErr = fs.reso_err;


    gr_resp_vs_energy700->SetPoint( i, energy/1000., mean /adcEnergyC700 * adcEnergyC600 );
    gr_resp_vs_energy700->SetPointError( i, energyErr/1000., meanErr);

    gr_reso_vs_energy700->SetPoint( i, energy/1000., reso );
    gr_reso_vs_energy700->SetPointError( i, energyErr/1000., resoErr );

    gr_dev700->SetPoint(i, energy/1000. , mean / energy / ( adcEnergyC700 / energyC700) );
    gr_dev700->SetPointError(i, 0, meanErr / energy / ( adcEnergyC700 / energyC700) );

 }




  for( unsigned i=0; i<runs600.size(); ++i ) {
    
    TFile* file = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%s.root", tag.c_str(),  runs600[i].c_str()));
    TTree* tree = (TTree*)file->Get("recoTree");
    
    FitStruct fs  = getFitResults( outputdir, tree ,  runs600[i].c_str() , 100, 19000000   , cut, whatToProject );

    float energy = beamEnergy600[i];
    float energyErr = 0.005*energy;
 
    float mean = fs.mean;
    float meanErr = fs.mean_err;

    float reso = fs.reso;
    float resoErr = fs.reso_err;

    gr_resp_vs_energy600->SetPoint( i, energy/1000., mean );
    gr_resp_vs_energy600->SetPointError( i, energyErr/1000., meanErr);

    gr_reso_vs_energy600->SetPoint( i, energy/1000., reso );
    gr_reso_vs_energy600->SetPointError( i, energyErr/1000., resoErr );

    gr_dev600->SetPoint(i, energy/1000. , mean / energy / ( adcEnergyC600 / energyC600) );
    gr_dev600->SetPointError(i, 0, meanErr / energy / ( adcEnergyC600 / energyC600) );

 }




 for( unsigned i=0; i<simulation.size(); ++i ) {

    ///////////////// (1x1) Shashlik ("real setup") //////////////////////////////
    TFile* fileS = TFile::Open(Form("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/H4UseThis1_5mmGap10mmBeam/Reco_%s.root", simulation[i].c_str()));

    TTree* treeS = (TTree*)fileS->Get("recoTree");

    FitStruct fsSim= getFitResults( outputdir, treeS , simulation[i] , beamEnergySimulation[i]/10. ,  beamEnergySimulation[i], SimCut , whatToProjectSim );
    float energySim = fsSim.mean;


   FitStruct fsSim_ps = addPhotoStatistics( fsSim );

    // ResoStruct rs = getResponseResolutionMC( outputdir, treeS, simulation[i], beamEnergySimulation[i] );
    //ResoStruct rs_ps = addPhotoStatistics( rs );

    gr_resp_vs_energy_simul->SetPoint( i, beamEnergySimulation[i]/1000., energySim* adcEnergyC600/energyS );
    gr_resp_vs_energy_simul->SetPointError( i,0, fsSim.mean_err * adcEnergyC600/energyS);
 
    gr_reso_vs_energy_simul->SetPoint( i, beamEnergySimulation[i]/1000., fsSim_ps.reso );
    gr_reso_vs_energy_simul->SetPointError( i,0,  fsSim_ps.reso_err );
    
    
    std::cout << "reso = " << fsSim_ps.reso << std::endl;

    yMax = fsSim.mean* adcEnergyC600/energyS + 20.;


    gr_dev_simul->SetPoint(i, beamEnergySimulation[i] /1000. , fsSim.mean / beamEnergySimulation[i] / ( energyS / energyC600 ) );
    gr_dev_simul->SetPointError(i, 0, fsSim.mean_err / beamEnergySimulation[i] / ( energyS / energyC600) );

  }




  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();


  TH2D* h2_axes = new TH2D( "axes", "", 10, 0., xMax/1000.+1, 10, 0., yMax );
  h2_axes->SetXTitle("Beam Energy [GeV]");
  h2_axes->SetYTitle("CeF_{3} Response [ADC]");
  h2_axes->Draw("");

  gr_resp_vs_energy->SetMarkerStyle(20);
  gr_resp_vs_energy->SetMarkerSize(1.6);
  gr_resp_vs_energy->SetMarkerColor(46);


  TF1* f1_line = new TF1("line", "[0] + [1]*x", 0., xMax/1000.+5 );
  gr_resp_vs_energy->Fit(f1_line,"RN");
  f1_line->SetLineWidth(1.);
  f1_line->SetLineColor(46);
  // f1_line->Draw("L same");


  gr_resp_vs_energy700->SetMarkerStyle(21);
  gr_resp_vs_energy700->SetMarkerSize(1.6);
  gr_resp_vs_energy700->SetMarkerColor(38);
  gr_resp_vs_energy700->Draw("p same");

  /* TF1* f1_line = new TF1("line", "[0] + [1]*x", 0., xMax/1000.+5 );
  gr_resp_vs_energy->Fit(f1_line,"RN");
  f1_line->SetLineWidth(1.);
  f1_line->SetLineColor(46);
  f1_line->Draw("L same");
  */

  /////////// 600 V ///////////////
  gr_resp_vs_energy600->SetMarkerStyle(21);
  gr_resp_vs_energy600->SetMarkerSize(1.6);
  gr_resp_vs_energy600->SetMarkerColor(41);
  gr_resp_vs_energy600->Draw("p same");
 /*
  TF1* f1_line = new TF1("line", "[0] + [1]*x", 0., xMax/1000.+5 );
  gr_resp_vs_energy->Fit(f1_line,"RN");
  f1_line->SetLineWidth(1.);
  f1_line->SetLineColor(46);
  f1_line->Draw("L same");
 */

  gr_resp_vs_energy->Draw("p same");

  gr_resp_vs_energy_simul->SetMarkerStyle(20);
  gr_resp_vs_energy_simul->SetMarkerSize(1.6);
  gr_resp_vs_energy_simul->SetMarkerColor(38);
  //  gr_resp_vs_energy_simul->Draw("p same");


  TF1* f1_lines = new TF1("lines", "[0] + [1]*x", 0., xMax/1000.+5 );
  gr_resp_vs_energy_simul->Fit(f1_lines,"RN");
  f1_lines->SetLineWidth(1.);
  f1_lines->SetLineColor(kBlack);
  f1_lines->Draw("L same");


  TLegend* leg0 = new TLegend(0.65, 0.4, 0.9, 0.2);
  // TLegend* leg0 = new TLegend(0.45, 0.2, 0.9, 0.4);  
leg0->SetTextSize(0.038);


  leg0->AddEntry(gr_resp_vs_energy,"Data 950 V","p");
  leg0->AddEntry(gr_resp_vs_energy700,"Data 700 V","p");
  leg0->AddEntry(gr_resp_vs_energy600,"Data 600 V","p");
  // leg0->AddEntry(f1_line,Form("Offset = %.1f\n #pm %.2f\n",f1_line->GetParameter(0), f1_line->GetParError(0) ),"L");
  //leg0->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.0f\n / %d",f1_line->GetChisquare(), f1_line->GetNDF() ), "");

  leg0->AddEntry(f1_lines,"MC","L");
  /*
  leg0->AddEntry(gr_resp_vs_energy_simul,"MC (1x1)","p");
  leg0->AddEntry(f1_lines,Form("Offset = %.1f\n #pm %.1f\n",f1_lines->GetParameter(0), f1_lines->GetParError(0) ),"L");
  leg0->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.0f\n / %d",f1_lines->GetChisquare(), f1_lines->GetNDF() ), "");
  */

  leg0->SetFillColor(0);
  leg0->Draw("same");

  TPaveText* label_top2 = new TPaveText();
  label_top2 = DrawTools::getLabelTop("Electron Beam");
  label_top2->Draw("same");

  c1->SaveAs( Form( "%s/resp_vs_energy.pdf", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/resp_vs_energy.png", outputdir.c_str() ) );
    
  c1->Clear();







 ///////////////////////RESOLUTION///////////////////////////////
  TH2D* h2_axes2 = new TH2D( "axes", "", 100, -0.0, 205. , 10, 0.0, 5 );
 h2_axes2->SetXTitle("Beam Energy [GeV]");
 h2_axes2->SetYTitle("Energy Resolution [%]");
 h2_axes2->Draw("");
 

 // Data (1x1)
 gr_reso_vs_energy->SetMarkerStyle(20);
 gr_reso_vs_energy->SetMarkerSize(1.6);
 gr_reso_vs_energy->SetMarkerColor(46);
 gr_reso_vs_energy->Draw("p same");
 
 TF1 *fun= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",0.9, xMax/1000.+5.);
 fun->SetParameter(1,1.);
 gr_reso_vs_energy->Fit(fun,"RN");
 fun->SetLineWidth(1.);
 fun->SetLineColor(46);
 // fun->Draw("L same");

 // Data 700 V
 gr_reso_vs_energy700->SetMarkerStyle(20);
 gr_reso_vs_energy700->SetMarkerSize(1.6);
 gr_reso_vs_energy700->SetMarkerColor(38);
 // gr_reso_vs_energy700->Draw("p same");

 TF1 *fun7= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",0.9, xMax/1000.+5.);
 gr_reso_vs_energy700->Fit(fun7,"RN");
 fun7->SetLineWidth(1.);
 fun7->SetLineColor(38);
 // fun7->Draw("L same");


 // Data 600 V
 gr_reso_vs_energy600->SetMarkerStyle(20);
 gr_reso_vs_energy600->SetMarkerSize(1.6);
 gr_reso_vs_energy600->SetMarkerColor(41);
 // gr_reso_vs_energy600->Draw("p same");

 TF1 *fun6= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",1, xMax/1000.+5.);
 fun6->SetParameter(1, 1.);

 fun6->SetParameter(0, 12.);
 gr_reso_vs_energy600->Fit(fun6,"RN");
 fun6->SetLineWidth(1.);
 fun6->SetLineColor(41);
 // fun6->Draw("L same");


 // MC (1x1)
 gr_reso_vs_energy_simul->SetMarkerStyle(21);
 gr_reso_vs_energy_simul->SetMarkerSize(1.4);
 gr_reso_vs_energy_simul->SetMarkerColor(kBlack);
 gr_reso_vs_energy_simul->Draw("p same");
 
 TF1 *fun1= new TF1("fun1",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",0.005, xMax/1000.+5.);
 // TF1 *fun1= new TF1("fun1",  "sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0.02, xMax/1000.);
 //  fun1->SetParameter(2,0.01);
 gr_reso_vs_energy_simul->Fit(fun1,"RN");

 fun1->SetMarkerStyle(21);
 fun1->SetMarkerSize(1.4);
 fun1->SetLineWidth(1.);
 fun1->SetLineColor(kBlack);
 fun1->Draw("L same");




 TLegend* leg4 = new TLegend(0.5, 0.9-0.06*4 +0.1, 0.9, 0.9);
 // TLegend* leg4 = new TLegend(0.3, 0.9-0.06*10 +0.1, 0.78, 0.9);
 leg4->SetTextSize(0.038);

 leg4->AddEntry(gr_reso_vs_energy,"Data 950 V","p");
// leg4->AddEntry(fun,Form("S = %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun->GetParameter(0), (fun->GetParError(0)),"%" ),"L");
// leg4->AddEntry( (TObject*)0,Form("C =   %.2f\n #pm %.2f\n %s",(fun->GetParameter(1)), (fun->GetParError(1)),"%" ),"");
// leg4->AddEntry( (TObject*)0,Form("N = %.2f\n #pm %.2f\n %s / E [GeV]",(fun1->GetParameter(2)), (fun1->GetParError(2)),"%" ),"");
// leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun->GetChisquare(), fun->GetNDF() ), "");
//

// leg4->AddEntry(gr_reso_vs_energy700,"Data 700 HV","p");
 // leg4->AddEntry(fun7,Form("S = %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun7->GetParameter(0), (fun7->GetParError(0)),"%" ),"L");
// leg4->AddEntry( (TObject*)0,Form("C =   %.2f\n #pm %.2f\n %s",(fun7->GetParameter(1)), (fun7->GetParError(1)),"%" ),"");
// leg4->AddEntry( (TObject*)0,Form("N =  %.2f\n #pm %.2f\n %s / E [GeV]",(fun7->GetParameter(2)), (fun7->GetParError(2)),"%" ),"");
// leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun7->GetChisquare(), fun7->GetNDF() ), "");
//

// leg4->AddEntry(gr_reso_vs_energy600,"Data 600 HV","p");
// leg4->AddEntry(fun6,Form("S = %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun6->GetParameter(0), (fun6->GetParError(0)),"%" ),"L");
// leg4->AddEntry( (TObject*)0,Form("C =   %.2f\n #pm %.2f\n %s",(fun6->GetParameter(1)), (fun6->GetParError(1)),"%" ),"");
// leg4->AddEntry( (TObject*)0,Form("N = %.2f\n #pm %.2f\n %s / E [GeV]",(fun6->GetParameter(2)), (fun6->GetParError(2)),"%" ),"");
// leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun6->GetChisquare(), fun6->GetNDF() ), "");


 leg4->AddEntry(fun1,"MC (1x1)","LP");
 // leg4->AddEntry(fun1,Form("S = %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun1->GetParameter(0), (fun1->GetParError(0)),"%" ),"L");
 // leg4->AddEntry( (TObject*)0,Form("C =   %.2f\n #pm %.2f\n %s",(fun1->GetParameter(1)), (fun1->GetParError(1)),"%" ),"");
 // leg4->AddEntry( (TObject*)0,Form("N =   %.2f\n #pm %.2f\n %s / E [GeV]",(fun1->GetParameter(2)), (fun1->GetParError(2)),"%" ),"");
 // leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun1->GetChisquare(), fun1->GetNDF() ), "");

 leg4->SetFillColor(0);
 leg4->Draw("same");
 
 TPaveText* label_low = new TPaveText(0.165,0.175,0.5,0.21, "brNDC");
 label_low->SetFillColor(kWhite);
 label_low->SetTextSize(0.038);
 label_low->SetTextAlign(11); // align right
 label_low->SetTextFont(62);
 label_low->AddText( "W-CeF_{3} Single Tower");
 // label_low->Draw("same");
 
  label_top2->Draw("same");

 
 c1->SaveAs( Form( "%s/resolution.pdf", outputdir.c_str() ) );
 c1->SaveAs( Form( "%s/resolution.png", outputdir.c_str() ) );
 c1->SaveAs( Form( "%s/resolution.eps", outputdir.c_str() ) );












  c1->Clear();

  //////////////////DEVIATION////////////////////////////

  TH2D* h2_axes3 = new TH2D( "axes", "", 10, 0., xMax/1000.+1, 10, 0., 1.6 );
  h2_axes3->SetXTitle("Beam Energy [GeV]");
  h2_axes3->SetYTitle("CeF_{3} Response / Beam Energy ");
  h2_axes3->Draw("");

  gr_dev->SetMarkerStyle(21);
  gr_dev->SetMarkerSize(1.6);
  gr_dev->SetMarkerColor(46);
  gr_dev->Draw("p same");

  gr_dev700->SetMarkerStyle(21);
  gr_dev700->SetMarkerSize(1.6);
  gr_dev700->SetMarkerColor(38);
  gr_dev700->Draw("p same");

  gr_dev600->SetMarkerStyle(20);
  gr_dev600->SetMarkerSize(1.6);
  gr_dev600->SetMarkerColor(41);
  gr_dev600->Draw("p same");

  gr_dev_simul->SetMarkerStyle(20);
  gr_dev_simul->SetMarkerSize(1.6);
  gr_dev_simul->SetMarkerColor(kBlack);
  gr_dev_simul->Draw("p same");



  TLegend* legd0 = new TLegend(0.55, 0.2 + 0.06 *4 , 0.9, 0.2);
  legd0->SetTextSize(0.038);

  legd0->AddEntry(gr_dev,"Data 950 V","p");
  legd0->AddEntry(gr_dev700,"Data 700 V","p");
  legd0->AddEntry(gr_dev600,"Data 600 V","p");
  legd0->AddEntry(gr_dev_simul,"MC (1x1)","p");
  
  legd0->SetFillColor(0);
  legd0->Draw("same");

  label_top2->Draw("same");

  c1->SaveAs( Form( "%s/deviation.pdf", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/deviation.png", outputdir.c_str() ) );






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

  
  double peakpos = f1->GetParameter(1);
  double sigma = f1->GetParameter(2);
  double fitmin;
  double fitmax;
  
  if( (peakpos-5*sigma) < 1000 ){ 
    fitmin = 1000;
  } else{ 
    fitmin = peakpos-5*sigma;
  }

  if(((peakpos>0) && (sigma > 2*peakpos)) ) { 
    fitmax = peakpos*2;
  } else{ 
    fitmax = peakpos+5*sigma;
  }

  if (sigma < 100) sigma = 100;



  TH1D* h2 = new TH1D("h2", "", 200, fitmin,fitmax);
  tree->Project( "h2", whatToProject.c_str(), cut.c_str() );

  RooRealVar x("x","ADC Channel", fitmin, fitmax);
  RooDataHist data("data","dataset with x",x,Import(*h2) );
  
  RooPlot* frame;
  RooPlot* xframe = x.frame();   
  
  frame = x.frame("Title");
  data.plotOn(frame);  //this will show histogram data points on canvas
  //  data.statOn(frame);  //this will display hist stat on canvas
  
  RooRealVar meanr("meanr","Mean",peakpos,peakpos-2*sigma, peakpos+2*sigma);
  RooRealVar width("width","#sigma",sigma, 100.0, 5.*sigma);
  RooRealVar A("A","Dist",1., 0.0, 5.0);
  RooRealVar N("N","Deg",3, 0.0, 15);

  meanr.setRange(1000. , 3000000.);
  width.setRange(100, 1000000);
  
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
  float  nADC = rs.mean/171.6*3236;
  //  float  nADC = rs.resp/169.445*3224.36;
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
