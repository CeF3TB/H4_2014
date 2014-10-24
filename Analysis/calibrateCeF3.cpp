#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h> 

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h" 
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLine.h"

#include "DrawTools.h"

#include "TApplication.h"


void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name );
TF1* fitSingleElectronPeak( const std::string& outputdir, int i, TTree* tree );
TF1* fitSingleElectronPeakCentral( const std::string& outputdir, int i, TTree* tree );
TF1* checkTotalResolution( const std::string& outputdir, TTree* tree );
void checkIntercalibration(std::vector<float> constant, std::vector<float> const_uncert ,const std::string& outputdir, const std::string& runName, const std::string& tag);
float sumVector( std::vector<float> v );
bool savePlots=true;
bool checkIntercal = true;


int main( int argc, char* argv[] ) {

  TApplication* a = new TApplication("a", 0, 0);

  DrawTools::setStyle();


  std::string inputDir = "./analysisTrees";
  std::string runName = "259";
  std::string tag = "V01";

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
    std::cout<<"./calibrateCef3 runNr tag [inputDir]"<<std::endl;
    exit(12345);
  }

  if(tag!="default") inputDir = inputDir + "_"+tag;



  TFile* file = TFile::Open(Form("%s/Reco_%s.root", inputDir.c_str(),runName.c_str()));
  std::cout<<"opening file:"<<file->GetName();


  TTree* tree = (TTree*)file->Get("recoTree");

  std::string outputdir = "CeF3Calibration/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  std::string ofsName = outputdir +Form( "/constants_%s_%s.txt", runName.c_str(), tag.c_str() );
  ofstream ofs(ofsName.c_str());
  std::string ofsNameU = outputdir + Form( "/constants_uncert_%s_%s.txt", runName.c_str(), tag.c_str() );
  ofstream ofsU(ofsNameU.c_str());


  int nentries = tree->GetEntries();

  std::vector<float> cef3_calibration;
  std::vector<float> cef3_calib_uncert;


  for( unsigned i=0; i<4; ++i ) {
    TF1* f1 = fitSingleElectronPeak( outputdir, i, tree );
    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    float mean_err = f1->GetParError(1);
    std::cout << std::endl;
    std::cout << "Channel " << i << std::endl;
    std::cout << "  Mean       : " << mean << std::endl;
    std::cout << "  Sigma      : " << sigma << std::endl;
    std::cout << "  Resolution : " << sigma/mean << std::endl;

    cef3_calibration.push_back(mean);
    cef3_calib_uncert.push_back(mean_err);
  }


  float cef3CalibrationAverage = sumVector(cef3_calibration)/cef3_calibration.size();

  for(unsigned i=0; i<cef3_calibration.size(); ++i ){
    cef3_calib_uncert[i] = abs(cef3CalibrationAverage-cef3_calibration[i]/4.)/(cef3_calibration[i]*cef3_calibration[i])*cef3_calib_uncert[i];
    ofsU << cef3_calib_uncert[i] << std::endl;

    cef3_calibration[i] = cef3CalibrationAverage/cef3_calibration[i];
    ofs << cef3_calibration[i] << std::endl;
  }



  ofs.close();
  ofsU.close();


  if(checkIntercal == true){
    checkIntercalibration(cef3_calibration, cef3_calib_uncert, outputdir, runName, tag);
  }

  std::cout << "-> Saved constants in: " << ofsName << std::endl;


  checkTotalResolution( outputdir, tree );




  return 0;

}



TF1* fitSingleElectronPeak( const std::string& outputdir, int i, TTree* tree ) {

  std::string histoName(Form("h1_%d", i));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 1000, 0., 390000.);
  //tree->Project( histoName.c_str(), Form("cef3_corr[%d]", i), "");
  tree->Project( histoName.c_str(), Form("cef3_corr[%d]", i));

  //tree->Project( histoName.c_str(), Form("cef3_corr[%d]", i), "(isSingleEle_scintFront==1 && nHodoClustersX==1 && nHodoClustersY==1 )");

  TF1* f1 = new TF1( Form("gaus_%d", i), "gaus", 5000., 390000.);
  //  TF1* f1 = new TF1( Form("gaus_%d", i), "gaus", 400., 1200.);
    f1->SetParameter(0, 300000.);
    f1->SetParameter(1, 250000.);
    f1->SetParameter(2, 15000.);

  doSingleFit( h1, f1, outputdir, Form("%d", i) );

  return f1;

}




TF1* checkTotalResolution( const std::string& outputdir, TTree* tree ) {

  std::string histoName("h1_tot");
  TH1D* h1 = new TH1D(histoName.c_str(), "", 2000, 200000.,1500000.);
  tree->Project( histoName.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "");

  // tree->Project( histoName.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(scintFront>500. && scintFront<2000. && nHodoClustersX==1 && nHodoClustersY==1)");

  //  TF1* f1 = new TF1("gaus_tot", "gaus", 1600., 4800.);
  TF1* f1 = new TF1("gaus_tot", "gaus", 200000.,1500000. );
  // f1->SetParameter(0, 3000000.);
  //f1->SetParameter(1, 3000000.);
  //f1->SetParameter(2, 60000.);

  f1->SetLineColor(kRed);

  doSingleFit( h1, f1, outputdir, "tot" );

  std::cout << std::endl;
  std::cout << std::endl;
  float mean  = f1->GetParameter(1);
  float sigma = f1->GetParameter(2);
  float reso = sigma/mean;
  std::cout << "Total " << std::endl;
  std::cout << "  Mean       : " << mean << std::endl;
  std::cout << "  Sigma      : " << sigma << std::endl;
  std::cout << "  Resolution : " << sigma/mean << std::endl;
  std::cout << "Corresponds to a stochastic term of: " << 100.*reso*sqrt(0.5) << " %" << std::endl;

  return f1;

}


void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name ) {

  h1->Fit( f1, "RQN" );

  int niter = 4.;
  float nSigma =1.5 ;

  for( unsigned iter=0; iter<niter; iter++ ) {

    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    float fitMin = mean - nSigma*sigma;
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

  if(savePlots){
    c1->SaveAs( Form("%s/fit_%s.png", outputdir.c_str(), name.c_str()) );
  }

  delete c1;

}

float sumVector( std::vector<float> v ) {

  float sum=0.;
  for( unsigned i=0; i<v.size(); ++i ) sum += v[i];

  return sum;

}





void checkIntercalibration(std::vector<float> constant, std::vector<float> const_uncert, const std::string& outputdir, const std::string& runName, const std::string& tag){

  TCanvas* canny = new TCanvas("canny", "",200,200);

  int n=4;
  float vv[4] = { 0.05,1.05,2.05,3.05 };
  std::vector<float> x(&vv[0], &vv[0]+4);
  std::vector<float> xerr (4,0);

  TGraphErrors* graf = new TGraphErrors(n,&(x[0]),&(constant[0]),&(xerr[0]),&(const_uncert[0]));
  graf->SetMarkerStyle(20);
  graf->SetMarkerSize(0.5);
  graf->SetMarkerColor(38);



  float vv2[4] = { -0.05,0.95,1.95,2.95 };
  std::vector<float> x2(&vv2[0], &vv2[0]+4);
  float constantBTF[4] = {1.10023,0.9082,1.02656,0.984356};
  float constantBTF_uncert[4] = {0.0012271,0.00130331,0.00118187,0.00121551};

  TGraphErrors* graf2 = new TGraphErrors(n,&(x2[0]),&(constantBTF[0]),&(xerr[0]),&(constantBTF_uncert[0]));
  graf2->SetMarkerStyle(21);
  graf2->SetMarkerSize(0.5);
  graf2->SetMarkerColor(46);





  TMultiGraph *multi= new TMultiGraph();
  multi->Add(graf);
  multi->Add(graf2);
  multi->SetTitle(" ;Channel Nr.; Correction Factor");
  multi->Draw("AP");



  multi->GetYaxis()->SetRangeUser(0.85,1.15);
  //multi->GetYaxis()->SetRangeUser(0.9,1.1);
  multi->Draw("AP");
  canny->Update();


  TLine* lin = new TLine(-0.2,1.,3.2,1.);
  lin->SetLineColor(kRed);
  lin->Draw();

  TLegend* leg = new TLegend(0.6, 0.9-0.06*2, 0.9, 0.9); 
  leg->AddEntry(graf,"H4","P");
  leg->AddEntry(graf2,"Frascati","P");
  leg->SetFillColor(0);
  leg->Draw("same");

  if(savePlots==1){  canny->SaveAs( Form( "%s/corrPlot_%s_%s.pdf", outputdir.c_str(), runName.c_str(), tag.c_str()  ));
}

  delete canny;

}
