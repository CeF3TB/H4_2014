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

#include "DrawTools.h"


struct ResoStruct {

  float resp;
  float resp_error;

  float reso;
  float reso_error;

};

ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name, float energySim );
ResoStruct getRespAndReso( TF1* f1, float energySim );
float getRatioError( float num, float denom, float numErr, float denomErr );
ResoStruct addPhotoStatistics( ResoStruct rs );

void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name );
TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree, float HV );


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

  TFile* file950 = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(), runs950[0].c_str()));
  TTree* tree950 = (TTree*)file950->Get("recoTree");
 
  // for CeF3:
  TF1* energyfunc950 = fitSingleElectronPeak( outputdir, runs950[0], tree950, 950 );
  float energyC950= beamEnergy950[0];
  float adcEnergyC950 = energyfunc950->GetParameter(1);

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

  TFile* file700 = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(),runs700[0].c_str()));
  TTree* tree700 = (TTree*)file700->Get("recoTree");

  TF1* energyfunc700 = fitSingleElectronPeak( outputdir, runs700[0], tree700, 700. );

  float energyC700= beamEnergy700[0];
  float adcEnergyC700 = energyfunc700->GetParameter(1);

 
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

  TFile* file600 = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(),runs600[0].c_str()));
  TTree* tree600 = (TTree*)file600->Get("recoTree");

  TF1* energyfunc600 = fitSingleElectronPeak( outputdir, runs600[0], tree600, 600 );

  float energyC600= beamEnergy600[0];
  float adcEnergyC600 = energyfunc600->GetParameter(1);



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
  
  TFile* energyfileS = TFile::Open("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/H4new/Reco_Simulation50.root" );
  TTree* energytreeS = (TTree*)energyfileS->Get("recoTree");
  
  ResoStruct energyrsS = getResponseResolutionMC( outputdir, energytreeS,simulation[0],beamEnergySimulation[0]);
  
  float energyS = energyrsS.resp;



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
    
    TFile* file = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(),  runs950[i].c_str()));
    TTree* tree = (TTree*)file->Get("recoTree");
    
    TF1* thisFunc = fitSingleElectronPeak( outputdir, runs950[i], tree, 950. );
    
    float energy = beamEnergy950[i];
    float energyErr = 0.005*energy;
 
    float mean = thisFunc->GetParameter(1);
    float meanErr = thisFunc->GetParError(1);

    float rms = thisFunc->GetParameter(2);
    float rmsErr = thisFunc->GetParError(2);

    float reso = 100.* rms/mean ;
    float resoErr = 100* getRatioError( rms, mean, rmsErr, meanErr);

    std::cout << "mean at "<< energy << " = " << mean << std::endl;
    gr_resp_vs_energy->SetPoint( i, energy/1000., mean /adcEnergyC950 * adcEnergyC600 );
    gr_resp_vs_energy->SetPointError( i, energyErr/1000., meanErr);

    gr_reso_vs_energy->SetPoint( i, energy/1000., reso );
    gr_reso_vs_energy->SetPointError( i, energyErr/1000., resoErr );

    gr_dev->SetPoint(i, energy/1000. , mean / energy / ( adcEnergyC950 / energyC950) );
    gr_dev->SetPointError(i, 0, meanErr / energy / ( adcEnergyC950 / energyC950) );

 }


  for( unsigned i=0; i<runs700.size(); ++i ) {
    
    TFile* file = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(),  runs700[i].c_str()));
    TTree* tree = (TTree*)file->Get("recoTree");
    
    TF1* thisFunc = fitSingleElectronPeak( outputdir, runs700[i], tree , 700.);
    
    float energy = beamEnergy700[i];
    float energyErr = 0.005*energy;
 
    float mean = thisFunc->GetParameter(1);
    float meanErr = thisFunc->GetParError(1);

    float rms = thisFunc->GetParameter(2);
    float rmsErr = thisFunc->GetParError(2);

    float reso = 100.* rms/mean ;
    float resoErr = 100* getRatioError( rms, mean, rmsErr, meanErr);

    gr_resp_vs_energy700->SetPoint( i, energy/1000., mean /adcEnergyC700 * adcEnergyC600 );
    gr_resp_vs_energy700->SetPointError( i, energyErr/1000., meanErr);

    gr_reso_vs_energy700->SetPoint( i, energy/1000., reso );
    gr_reso_vs_energy700->SetPointError( i, energyErr/1000., resoErr );

    gr_dev700->SetPoint(i, energy/1000. , mean / energy / ( adcEnergyC700 / energyC700) );
    gr_dev700->SetPointError(i, 0, meanErr / energy / ( adcEnergyC700 / energyC700) );

 }




  for( unsigned i=0; i<runs600.size(); ++i ) {
    
    TFile* file = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(),  runs600[i].c_str()));
    TTree* tree = (TTree*)file->Get("recoTree");
    
    TF1* thisFunc = fitSingleElectronPeak( outputdir, runs600[i], tree , 600.);
    
    float energy = beamEnergy600[i];
    float energyErr = 0.005*energy;
 
    float mean = thisFunc->GetParameter(1);
    float meanErr = thisFunc->GetParError(1);

    float rms = thisFunc->GetParameter(2);
    float rmsErr = thisFunc->GetParError(2);

    float reso = 100.* rms/mean ;
    float resoErr = 100* getRatioError( rms, mean, rmsErr, meanErr);

    gr_resp_vs_energy600->SetPoint( i, energy/1000., mean );
    gr_resp_vs_energy600->SetPointError( i, energyErr/1000., meanErr);

    gr_reso_vs_energy600->SetPoint( i, energy/1000., reso );
    gr_reso_vs_energy600->SetPointError( i, energyErr/1000., resoErr );

    gr_dev600->SetPoint(i, energy/1000. , mean / energy / ( adcEnergyC600 / energyC600) );
    gr_dev600->SetPointError(i, 0, meanErr / energy / ( adcEnergyC600 / energyC600) );

 }




 for( unsigned i=0; i<simulation.size(); ++i ) {

    ///////////////// (1x1) Shashlik ("real setup") //////////////////////////////
    TFile* fileS = TFile::Open(Form("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/H4new/Reco_%s.root", simulation[i].c_str()));

    TTree* treeS = (TTree*)fileS->Get("recoTree");

    ResoStruct rs = getResponseResolutionMC( outputdir, treeS, simulation[i], beamEnergySimulation[i] );
    ResoStruct rs_ps = addPhotoStatistics( rs );

    gr_resp_vs_energy_simul->SetPoint( i, beamEnergySimulation[i]/1000., rs.resp* adcEnergyC600/energyS );
    gr_resp_vs_energy_simul->SetPointError( i,0, rs.resp_error * adcEnergyC600/energyS);
 
    gr_reso_vs_energy_simul->SetPoint( i, beamEnergySimulation[i]/1000., rs_ps.reso );
    gr_reso_vs_energy_simul->SetPointError( i,0,  rs_ps.reso_error );
    
    
    std::cout << "reso = " << rs_ps.reso << std::endl;

    yMax = rs.resp* adcEnergyC600/energyS + 20.;


    gr_dev_simul->SetPoint(i, beamEnergySimulation[i] /1000. , rs.resp / beamEnergySimulation[i] / ( energyS / energyC600 ) );
    gr_dev_simul->SetPointError(i, 0, rs.resp_error / beamEnergySimulation[i] / ( energyS / energyC600) );

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
  
  c1->Clear();







 ///////////////////////RESOLUTION///////////////////////////////
  TH2D* h2_axes2 = new TH2D( "axes", "", 100, -0.0, 205. , 10, 0.0, 30 );
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
 fun->Draw("L same");

 // Data 700 V
 gr_reso_vs_energy700->SetMarkerStyle(20);
 gr_reso_vs_energy700->SetMarkerSize(1.6);
 gr_reso_vs_energy700->SetMarkerColor(38);
 gr_reso_vs_energy700->Draw("p same");

 TF1 *fun7= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",0.9, xMax/1000.+5.);
 gr_reso_vs_energy700->Fit(fun7,"RN");
 fun7->SetLineWidth(1.);
 fun7->SetLineColor(38);
 fun7->Draw("L same");


 // Data 600 V
 gr_reso_vs_energy600->SetMarkerStyle(20);
 gr_reso_vs_energy600->SetMarkerSize(1.6);
 gr_reso_vs_energy600->SetMarkerColor(41);
 gr_reso_vs_energy600->Draw("p same");

 TF1 *fun6= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",1, xMax/1000.+5.);
 fun6->SetParameter(1, 1.);

 fun6->SetParameter(0, 12.);
 gr_reso_vs_energy600->Fit(fun6,"RN");
 fun6->SetLineWidth(1.);
 fun6->SetLineColor(41);
 fun6->Draw("L same");


 // MC (1x1)
 gr_reso_vs_energy_simul->SetMarkerStyle(20);
 gr_reso_vs_energy_simul->SetMarkerSize(1.6);
 gr_reso_vs_energy_simul->SetMarkerColor(kBlack);
 gr_reso_vs_energy_simul->Draw("p same");
 
 TF1 *fun1= new TF1("fun1",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",0.005, xMax/1000.+5.);
 // TF1 *fun1= new TF1("fun1",  "sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0.02, xMax/1000.);
 //  fun1->SetParameter(2,0.01);
 gr_reso_vs_energy_simul->Fit(fun1,"RN");
 fun1->SetLineWidth(1.);
 fun1->SetLineColor(kBlack);
 fun1->Draw("L same");




 TLegend* leg4 = new TLegend(0.3, 0.9-0.06*10 +0.1, 0.78, 0.9);
 leg4->SetTextSize(0.038);

 leg4->AddEntry(gr_reso_vs_energy,"Data 950 V","p");
 leg4->AddEntry(fun,Form("S = %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun->GetParameter(0), (fun->GetParError(0)),"%" ),"L");
 leg4->AddEntry( (TObject*)0,Form("C =   %.2f\n #pm %.2f\n %s",(fun->GetParameter(1)), (fun->GetParError(1)),"%" ),"");
 leg4->AddEntry( (TObject*)0,Form("N = %.2f\n #pm %.2f\n %s / E [GeV]",(fun1->GetParameter(2)), (fun1->GetParError(2)),"%" ),"");
 leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun->GetChisquare(), fun->GetNDF() ), "");


 leg4->AddEntry(gr_reso_vs_energy700,"Data 700 HV","p");
 leg4->AddEntry(fun7,Form("S = %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun7->GetParameter(0), (fun7->GetParError(0)),"%" ),"L");
 leg4->AddEntry( (TObject*)0,Form("C =   %.2f\n #pm %.2f\n %s",(fun7->GetParameter(1)), (fun7->GetParError(1)),"%" ),"");
 leg4->AddEntry( (TObject*)0,Form("N =  %.2f\n #pm %.2f\n %s / E [GeV]",(fun7->GetParameter(2)), (fun7->GetParError(2)),"%" ),"");
 leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun7->GetChisquare(), fun7->GetNDF() ), "");




 leg4->AddEntry(gr_reso_vs_energy600,"Data 600 HV","p");
 leg4->AddEntry(fun6,Form("S = %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun6->GetParameter(0), (fun6->GetParError(0)),"%" ),"L");
 leg4->AddEntry( (TObject*)0,Form("C =   %.2f\n #pm %.2f\n %s",(fun6->GetParameter(1)), (fun6->GetParError(1)),"%" ),"");
 leg4->AddEntry( (TObject*)0,Form("N = %.2f\n #pm %.2f\n %s / E [GeV]",(fun6->GetParameter(2)), (fun6->GetParError(2)),"%" ),"");
 leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun6->GetChisquare(), fun6->GetNDF() ), "");



 leg4->AddEntry(fun1,"MC (1x1)","L");
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












  c1->Clear();





  //////////////////DEVIATION////////////////////////////

  TH2D* h2_axes3 = new TH2D( "axes", "", 10, 0., xMax/1000.+1, 10, 0., 1.4 );
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






  return 0;

}















ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name, float energySim ) {

  gStyle->SetOptFit(1);

  TH1D* h1 = new TH1D( name.c_str(), "", 500, energySim/4.-10, energySim/2.5 +100 );
  tree->Project( name.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]" );


 TF1* f1 = new TF1( Form("f1_%s", name.c_str()), "gaus",  energySim/4.-10, energySim/2.5+10 );

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


float getRatioError( float num, float denom, float numErr, float denomErr ) {

  return sqrt( numErr*numErr/(denom*denom) + denomErr*denomErr*num*num/(denom*denom*denom*denom) );

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
























TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree, float HV ) {

  std::string histoName(Form("h1_%s", name.c_str() ));
  TH1D* h1;
  TF1* f1;

  if(HV == 950. ){

    h1 = new TH1D(histoName.c_str(), "", 500, 1000, 7000 );
  tree->Project( histoName.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]");
  f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 1000., 7000);
 
  } else  if(HV == 700. ){

    h1 = new TH1D(histoName.c_str(), "", 3500,50, 1700 );
  tree->Project( histoName.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]");
  f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 50., 1700);
     f1->SetParameter(0, 10.);
  }else{

    h1 = new TH1D(histoName.c_str(), "", 500, 0, 600 );
  tree->Project( histoName.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]");
  f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 0., 600);
 
  }

  // std::string histoName(Form("h1_%s", name.c_str() ));
  // TH1D* h1 = new TH1D(histoName.c_str(), "", 5100, 0., 10000 );
  //  TH1D* h1 = new TH1D(histoName.c_str(), "", 500, 0.,180000.);

 //  TF1* f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 500., 1800000.);
  //  TF1* f1 = new TF1( Form("gaus_%d", i), "gaus", 400., 1200.);
  //     f1->SetParameter(0, 60.);
   f1->SetParameter(1, 1000.);
   f1->SetParameter(2, 100.);

   f1->SetLineColor(kRed);

    doSingleFit( h1, f1, outputdir, name );

  return f1;

}







void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name ) {

  h1->Fit( f1, "RQN" );

  int niter = 4.;
  float nSigma =1.5 ;

  for( unsigned iter=0; iter<niter; iter++ ) {

    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    float fitMin = mean - nSigma*sigma*0.8;
    float fitMax = mean + nSigma*sigma;
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
