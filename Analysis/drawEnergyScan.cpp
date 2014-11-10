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
TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree );


int main( int argc, char* argv[] ) {

  std::string tag="V01";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }
  std::cout << "-> Using tag: " << tag << std::endl;

  std::string outputdir = "EnergyScan/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  DrawTools::setStyle();
  TGaxis::SetMaxDigits(3);

  std::vector<std::string> runs; 
  std::vector<float> beamEnergy;
 
  std::vector<std::string> simulation; 
  std::vector<float> beamEnergySimulation; 

  //////////////////DATA/////////////

  /*
  runs.push_back("274");
  beamEnergy.push_back(50000);
  runs.push_back("329");
  //  runs.push_back("392");
  beamEnergy.push_back(100000); 
  runs.push_back("363");
  beamEnergy.push_back(20000);
  runs.push_back("355");
  //  runs.push_back("392");
  beamEnergy.push_back(150000); 
  //runs.push_back("387");
  runs.push_back("377");
  beamEnergy.push_back(200000);
  */

  /*
  runs.push_back("418");;
  beamEnergy.push_back(100000);

  runs.push_back("431");
  beamEnergy.push_back(10000);
  runs.push_back("457");
  beamEnergy.push_back(15000);
  runs.push_back("455");
  beamEnergy.push_back(20000);
  runs.push_back("398");
  beamEnergy.push_back(50000);
  runs.push_back("407");
  beamEnergy.push_back(150000); 
  runs.push_back("428");
  beamEnergy.push_back(200000);
  */

  runs.push_back("503");;
  beamEnergy.push_back(100000);

  runs.push_back("513");
  beamEnergy.push_back(15000);
  runs.push_back("508");
  beamEnergy.push_back(20000);
  runs.push_back("497");
  beamEnergy.push_back(50000);
  runs.push_back("481");
  beamEnergy.push_back(150000); 
  runs.push_back("487");
  beamEnergy.push_back(200000);



  ///////SIMULATION/////////////
  simulation.push_back("Simulation10");
  beamEnergySimulation.push_back(10000.);
  simulation.push_back("Simulation20");
  beamEnergySimulation.push_back(20000.);

  simulation.push_back("Simulation50");
  beamEnergySimulation.push_back(50000.);

  simulation.push_back("Simulation100");
  beamEnergySimulation.push_back(100000.);

  simulation.push_back("Simulation150");
  beamEnergySimulation.push_back(150000.);
  
  simulation.push_back("Simulation200");
  beamEnergySimulation.push_back(200000.);
  

  float xMax = beamEnergy[beamEnergy.size()-1] +5000;

  TGraphErrors* gr_resp_vs_energy = new TGraphErrors(0);
  TGraphErrors* gr_reso_vs_energy = new TGraphErrors(0);

  TGraphErrors* gr_resp_vs_energy_simul = new TGraphErrors(0);
  TGraphErrors* gr_reso_vs_energy_simul = new TGraphErrors(0);

  TGraphErrors* gr_dev = new TGraphErrors(0);
  TGraphErrors* gr_dev_simul = new TGraphErrors(0);



  TFile* fileMean = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(),runs[0].c_str()));
  TTree* treeMean = (TTree*)fileMean->Get("recoTree");
 
  //and for CeF3:
  TF1* energyfuncC = fitSingleElectronPeak( outputdir, runs[0], treeMean );

  float energyC= beamEnergy[0];

  float adcEnergyC = energyfuncC->GetParameter(1);

  //for the simulated Cef3:
  TFile* energyfileS = TFile::Open("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/newnewH4/Reco_Simulation100.root");
  TTree* energytreeS = (TTree*)energyfileS->Get("recoTree");
  
  ResoStruct energyrsS = getResponseResolutionMC( outputdir, energytreeS,"100", 100000);
  
  float energyS = energyrsS.resp;

 
  float yMax = 0;

  for( unsigned i=0; i<runs.size(); ++i ) {
    
    TFile* file = TFile::Open(Form("analysisTrees_%s/Reco_%s.root", tag.c_str(),  runs[i].c_str()));
    TTree* tree = (TTree*)file->Get("recoTree");
    
    TF1* thisFunc = fitSingleElectronPeak( outputdir, runs[i], tree );
    
    float energy = beamEnergy[i];
    float energyErr = 0.005*energy;
 
    float mean = thisFunc->GetParameter(1);
    float meanErr = thisFunc->GetParError(1);

    float rms = thisFunc->GetParameter(2);
    float rmsErr = thisFunc->GetParError(2);

    float reso = 100.* rms/mean ;
    float resoErr = 100* getRatioError( rms, mean, rmsErr, meanErr);

    std::cout << "mean at "<< energy << " = " << mean << std::endl;
    gr_resp_vs_energy->SetPoint( i, energy/1000., mean );
    gr_resp_vs_energy->SetPointError( i, energyErr/1000., meanErr);

    gr_reso_vs_energy->SetPoint( i, energy/1000., reso );
    gr_reso_vs_energy->SetPointError( i, energyErr/1000., resoErr );

    gr_dev->SetPoint(i, energy/1000. , mean / energy / ( adcEnergyC / energyC) );
    gr_dev->SetPointError(i, 0, meanErr / energy / ( adcEnergyC / energyC) );

 }


 for( unsigned i=0; i<simulation.size(); ++i ) {

    ///////////////// (1x1) Shashlik ("real setup") //////////////////////////////
    TFile* fileS = TFile::Open(Form("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/newnewH4/Reco_%s.root", simulation[i].c_str()));

    TTree* treeS = (TTree*)fileS->Get("recoTree");

    ResoStruct rs = getResponseResolutionMC( outputdir, treeS, simulation[i], beamEnergySimulation[i] );
    ResoStruct rs_ps = addPhotoStatistics( rs );

    gr_resp_vs_energy_simul->SetPoint( i, beamEnergySimulation[i]/1000., rs.resp* adcEnergyC/energyS );
    gr_resp_vs_energy_simul->SetPointError( i,0, rs.resp_error * adcEnergyC/energyS);
 
    gr_reso_vs_energy_simul->SetPoint( i, beamEnergySimulation[i]/1000., rs_ps.reso );
    gr_reso_vs_energy_simul->SetPointError( i,0,  rs_ps.reso_error );
    
    
    std::cout << "reso = " << rs_ps.reso << std::endl;

    yMax = rs.resp* adcEnergyC/energyS + 20.;


    gr_dev_simul->SetPoint(i, beamEnergySimulation[i] /1000. , rs.resp / beamEnergySimulation[i] / ( energyS / energyC ) );
    gr_dev_simul->SetPointError(i, 0, rs.resp_error / beamEnergySimulation[i] / ( energyS / energyC) );

  }




  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();


  TH2D* h2_axes = new TH2D( "axes", "", 10, 0., xMax/1000.+1, 10, 0., yMax );
  h2_axes->SetXTitle("Beam Energy [GeV]");
  h2_axes->SetYTitle("CeF_{3} Response [ADC]");
  h2_axes->Draw("");

  gr_resp_vs_energy->SetMarkerStyle(21);
  gr_resp_vs_energy->SetMarkerSize(1.6);
  gr_resp_vs_energy->SetMarkerColor(46);
  gr_resp_vs_energy->Draw("p same");

  TF1* f1_line = new TF1("line", "[0] + [1]*x", 0., xMax/1000.+5 );
  gr_resp_vs_energy->Fit(f1_line,"RN");
  f1_line->SetLineWidth(1.);
  f1_line->SetLineColor(46);
  f1_line->Draw("L same");

  gr_resp_vs_energy_simul->SetMarkerStyle(20);
  gr_resp_vs_energy_simul->SetMarkerSize(1.6);
  gr_resp_vs_energy_simul->SetMarkerColor(38);
  gr_resp_vs_energy_simul->Draw("p same");


  TF1* f1_lines = new TF1("lines", "[0] + [1]*x", 0., xMax/1000.+5 );
  gr_resp_vs_energy_simul->Fit(f1_lines,"RN");
  f1_lines->SetLineWidth(1.);
  f1_lines->SetLineColor(38);
  f1_lines->Draw("L same");


  TLegend* leg0 = new TLegend(0.55, 0.5, 0.9, 0.2);
  // TLegend* leg0 = new TLegend(0.45, 0.2, 0.9, 0.4);  
leg0->SetTextSize(0.038);


  leg0->AddEntry(gr_resp_vs_energy,"Data","p");
  leg0->AddEntry(f1_line,Form("Offset = %.1f\n #pm %.2f\n",f1_line->GetParameter(0), f1_line->GetParError(0) ),"L");
  leg0->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.0f\n / %d",f1_line->GetChisquare(), f1_line->GetNDF() ), "");

  leg0->AddEntry(gr_resp_vs_energy_simul,"MC (1x1)","p");
  leg0->AddEntry(f1_lines,Form("Offset = %.1f\n #pm %.1f\n",f1_lines->GetParameter(0), f1_lines->GetParError(0) ),"L");
  leg0->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.0f\n / %d",f1_lines->GetChisquare(), f1_lines->GetNDF() ), "");


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
 
 //TF1 *fun= new TF1("fun", "sqrt([0]*[0]/(x))",50, xMax);
 TF1 *fun= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]/(x*x))",0.1, xMax/1000.+5.);
 gr_reso_vs_energy->Fit(fun,"RN");
 fun->SetLineWidth(1.);
 fun->SetLineColor(46);
 fun->Draw("L same");



 // MC (1x1)
 gr_reso_vs_energy_simul->SetMarkerStyle(20);
 gr_reso_vs_energy_simul->SetMarkerSize(1.6);
 gr_reso_vs_energy_simul->SetMarkerColor(38);
 gr_reso_vs_energy_simul->Draw("p same");
 
 TF1 *fun1= new TF1("fun1",  "sqrt([0]*[0]/x+[1]*[1])",0.0005, xMax/1000.+5.);
 // TF1 *fun1= new TF1("fun1",  "sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",0.02, xMax/1000.);
 //  fun1->SetParameter(2,0.01);
 gr_reso_vs_energy_simul->Fit(fun1,"RN");
 fun1->SetLineWidth(1.);
 fun1->SetLineColor(38);
 fun1->Draw("L same");




 TLegend* leg4 = new TLegend(0.3, 0.9-0.06*8, 0.78, 0.9);
 leg4->SetTextSize(0.038);

 leg4->AddEntry(gr_reso_vs_energy,"Data","p");
 leg4->AddEntry(fun,Form("S = %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun->GetParameter(0), (fun->GetParError(0)),"%" ),"L");
 leg4->AddEntry( (TObject*)0,Form("C =   %.2f\n #pm %.2f\n %s",(fun->GetParameter(1)), (fun->GetParError(1)),"%" ),"");
 // leg4->AddEntry( (TObject*)0,Form("N =   %.2f\n #pm %.2f\n %s / E [GeV]",(fun1->GetParameter(2)), (fun1->GetParError(2)),"%" ),"");
 leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun->GetChisquare(), fun->GetNDF() ), "");





 leg4->AddEntry(gr_reso_vs_energy_simul,"MC (1x1) H4","p");
 leg4->AddEntry(fun1,Form("S = %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun1->GetParameter(0), (fun1->GetParError(0)),"%" ),"L");
 leg4->AddEntry( (TObject*)0,Form("C =   %.2f\n #pm %.2f\n %s",(fun1->GetParameter(1)), (fun1->GetParError(1)),"%" ),"");
 // leg4->AddEntry( (TObject*)0,Form("N =   %.2f\n #pm %.2f\n %s / E [GeV]",(fun1->GetParameter(2)), (fun1->GetParError(2)),"%" ),"");
 leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun1->GetChisquare(), fun1->GetNDF() ), "");



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

  gr_dev_simul->SetMarkerStyle(20);
  gr_dev_simul->SetMarkerSize(1.6);
  gr_dev_simul->SetMarkerColor(38);
  gr_dev_simul->Draw("p same");



  TLegend* legd0 = new TLegend(0.55, 0.5, 0.9, 0.2);
  legd0->SetTextSize(0.038);

  legd0->AddEntry(gr_dev,"Data","p");
  legd0->AddEntry(gr_dev_simul,"MC (1x1)","p");
  
  legd0->SetFillColor(0);
  legd0->Draw("same");

  label_top2->Draw("same");

  c1->SaveAs( Form( "%s/deviation.pdf", outputdir.c_str() ) );






  return 0;

}















ResoStruct getResponseResolutionMC( const std::string& outputdir, TTree* tree, const std::string& name, float energySim ) {

  gStyle->SetOptFit(1);

  TH1D* h1 = new TH1D( name.c_str(), "", 500, energySim/3.-10, energySim/2.3+100 );
  tree->Project( name.c_str(),"cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]" );


 TF1* f1 = new TF1( Form("f1_%s", name.c_str()), "gaus",  energySim/6.-100, energySim/2.+10 );

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
























TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree ) {

  std::string histoName(Form("h1_%s", name.c_str() ));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 500, 0.,600.);
  //  TH1D* h1 = new TH1D(histoName.c_str(), "", 500, 0.,180000.);

  tree->Project( histoName.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]");

  TF1* f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 20., 600.);
  //  TF1* f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", 500., 1800000.);
  //  TF1* f1 = new TF1( Form("gaus_%d", i), "gaus", 400., 1200.);
      f1->SetParameter(0, 300.);
   f1->SetParameter(1, 200.);
   f1->SetParameter(2, 10.);

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
