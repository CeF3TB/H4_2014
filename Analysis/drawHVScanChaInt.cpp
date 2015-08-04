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
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/FitMethodFunction.h"
#include "HFitInterface.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
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

  double chi2;
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

void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name );
TF1* fitSingleElectronPeak( const std::string& outputdir, const std::string& name, TTree* tree );


// definition of shared parameter
// 20 GeV
int ipar20GeV[2] = { 0,   //noise term
		     1    //constant 20GeV
};
// 50GeV
int ipar50GeV[2] = { 0, // common noise term
		     2 // constant term 50GeV
};
// 100GeV
int ipar100GeV[2] = { 0, // common noise term
		      3 // constant term 100GeV
};
// 150GeV
int ipar150GeV[2] = { 0, // common noise term
		      4 // constant term 150GeV
};
/*
// 200GeV
int ipar200GeV[2] = { 0, // common noise term
		      4 // constant term 200GeV
};
*/

class GlobalChi2 : public ROOT::Math::FitMethodFunction { 

public:
  GlobalChi2( int dim, int npoints, ROOT::Math::FitMethodFunction & f1,   ROOT::Math::FitMethodFunction & f2,  ROOT::Math::FitMethodFunction & f3,  ROOT::Math::FitMethodFunction & f4) :
      ROOT::Math::FitMethodFunction(dim,npoints),
      fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3), fChi2_4(&f4) {}

   ROOT::Math::IMultiGenFunction * Clone() const { 
      // copy using default copy-ctor
      // i.e. function pointer will be copied (and not the functions)
      return new GlobalChi2(*this);   
   }

   double  DataElement( const double *par, unsigned int ipoint, double *g = 0) const { 
      // implement evaluation of single chi2 element
      double p1[2];
      for (int i = 0; i < 2; ++i) p1[i] = par[ipar20GeV[i] ];
      double p2[2]; 
      for (int i = 0; i < 2; ++i) p2[i] = par[ipar50GeV[i] ];
      double p3[2]; 
      for (int i = 0; i < 2; ++i) p3[i] = par[ipar100GeV[i] ];
      double p4[2]; 
      for (int i = 0; i < 2; ++i) p4[i] = par[ipar150GeV[i] ];
      //     double p5[2]; 
      //     for (int i = 0; i < 2; ++i) p5[i] = par[ipar150GeV[i] ];

      double g1[2];  double g2[2];  double g3[2];  double g4[2]; // double g5[2];
      double value = 0; 

     
      if (g != 0) {  for (int i = 0; i < 4; ++i) g[i] = 0;       }
      //      if (g != 0) {  for (int i = 0; i < 4; ++i) g[i] = 0;       }
      
      if (ipoint < fChi2_1->NPoints() ) {
	
	if (g != 0) { 
	  value = fChi2_1->DataElement(p1, ipoint, g1);
	  // update gradient values
	  for (int i= 0; i < 2; ++i) g[ipar20GeV[i]] = g1[i];         
	} else 
	  // no need to compute gradient in this case
	  value = fChi2_1->DataElement(p1, ipoint);
      }

      else if(ipoint > (fChi2_1->NPoints())-1 && ipoint < (fChi2_1->NPoints()+fChi2_2->NPoints()) )   
	{ //second/////////////////////////////////////
	unsigned int jpoint = ipoint - fChi2_1->NPoints();
	assert ( jpoint < fChi2_2->NPoints() );
	if ( g != 0) { 
	  value =  fChi2_2->DataElement(p2, jpoint, g2);
	  // update gradient values
	  for (int i= 0; i < 2; ++i) g[ipar50GeV[i]] = g2[i];
	} else 
	  // no need to compute gradient in this case
	  value =  fChi2_2->DataElement(p2, jpoint);
      }
      else if(ipoint > (fChi2_1->NPoints() +fChi2_2->NPoints() -1) && ipoint < (fChi2_1->NPoints()+fChi2_2->NPoints()+fChi2_3->NPoints() ) )   
	{  //third//////////////////////////////////////
	unsigned int kpoint = ipoint - fChi2_2->NPoints() - fChi2_1->NPoints() ;

	//	std::cout << "kpoint = " << kpoint << std::endl;

       	assert (kpoint < fChi2_3->NPoints()  );
	if ( g != 0) { 
	  value =  fChi2_3->DataElement(p3, kpoint, g3);
	  // update gradient values
	  for (int i= 0; i < 2; ++i) g[ipar100GeV[i]] = g3[i];
	}
	else 
	  // no need to compute gradient in this case
	  value =  fChi2_3->DataElement(p3, kpoint);
      }

      else
	//      else if(ipoint > (fChi2_1->NPoints() +fChi2_2->NPoints() +fChi2_3->NPoints() -1) && ipoint < (fChi2_1->NPoints()+fChi2_2->NPoints()+fChi2_3->NPoints() +fChi2_4->NPoints() ) )   
	{  //fourth//////////////////////////////////////
	unsigned int lpoint = ipoint - fChi2_3->NPoints() - fChi2_2->NPoints() - fChi2_1->NPoints() ;
       	assert (lpoint < fChi2_4->NPoints()  );
	if ( g != 0) { 
	  value =  fChi2_4->DataElement(p4, lpoint, g4);
	  // update gradient values
	  for (int i= 0; i < 2; ++i) g[ipar100GeV[i]] = g4[i];
	}
	else 
	  // no need to compute gradient in this case
	  value =  fChi2_4->DataElement(p4, lpoint);
      }

      /*
    else{  //fifth//////////////////////////////////////
	unsigned int mpoint = ipoint - fChi2_4->NPoints() - fChi2_3->NPoints()- fChi2_2->NPoints() - fChi2_1->NPoints() ;
     	assert (mpoint < fChi2_5->NPoints()  );
	if ( g != 0) { 
	  value =  fChi2_5->DataElement(p5, mpoint, g5);
	  // update gradient values
	  for (int i= 0; i < 2; ++i) g[ipar100GeV[i]] = g5[i];
	}
	else 
	  // no need to compute gradient in this case
	  value =  fChi2_5->DataElement(p5, mpoint);
      }
      */

      return value; 
   }
  
    
   // needed if want to use Fumili or Fumili2
   virtual Type_t Type() const { return ROOT::Math::FitMethodFunction::kLeastSquare; }

private:
   // parameter vector is first background (in common 1 and 2) and then is signal (only in 2)
   virtual double DoEval (const double *par) const {

      double p1[2];
      for (int i = 0; i < 2; ++i) p1[i] = par[ipar20GeV[i] ];
      double p2[2]; 
      for (int i = 0; i < 2; ++i) p2[i] = par[ipar50GeV[i] ];
      double p3[2]; 
      for (int i = 0; i < 2; ++i) p3[i] = par[ipar100GeV[i] ];
      double p4[2]; 
      for (int i = 0; i < 2; ++i) p4[i] = par[ipar150GeV[i] ];
      //   double p5[2]; 
      //     for (int i = 0; i < 2; ++i) p5[i] = par[ipar150GeV[i] ];

      return (*fChi2_1)(p1) + (*fChi2_2)(p2)  + (*fChi2_3)(p3) + (*fChi2_4)(p4);
   } 

   const  ROOT::Math::FitMethodFunction * fChi2_1;
   const  ROOT::Math::FitMethodFunction * fChi2_2;
   const  ROOT::Math::FitMethodFunction * fChi2_3;
   const  ROOT::Math::FitMethodFunction * fChi2_4;
  //  const  ROOT::Math::FitMethodFunction * fChi2_5;
};






int main( int argc, char* argv[] ) {

  std::string tag="V01";
  if( argc>1 ) {
    std::string tag_str(argv[1]);
    tag = tag_str;
  }
  std::cout << "-> Using tag: " << tag << std::endl;

  std::string outputdir = "HVScanChaInt/";
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  DrawTools::setStyle();   TGaxis::SetMaxDigits(3);

  std::string whatToProjectSim=Form("cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]");
  std::string whatToProject = Form("cef3_chaInt_corr[0]+cef3_chaInt_corr[1]+cef3_chaInt_corr[2]+cef3_chaInt_corr[3]");
  //  std::string whatToProject_bare = Form("cef3_chaInt_bare[0]+cef3_chaInt_bare_corr[1]+cef3_chaInt_bare_corr[2]+cef3_chaInt_bare_corr[3]");
  
  
  std::vector<std::string> simulation; 
  std::vector<float> beamEnergySimulation; 
  
  //////////////////DATA/////////////
  /////////DATA taken at 950V on CeF3//////////
  std::vector<std::string> runs950; 
  std::vector<float> beamEnergy950;
  runs950.push_back("323"); //274 a previous one
  beamEnergy950.push_back(50000);
  runs950.push_back("363");
  beamEnergy950.push_back(20000);
  runs950.push_back("329");
  beamEnergy950.push_back(100000); 
  runs950.push_back("356");  
  beamEnergy950.push_back(150000); 
  runs950.push_back("377");   
  beamEnergy950.push_back(200000);

  TFile* file950 = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%s.root", tag.c_str(), runs950[0].c_str()));
  TTree* tree950 = (TTree*)file950->Get("recoTree");

  std::string cut = "abs(0.5*(position_X1+position_X2))< 3 && abs( 0.5* (position_Y1+position_Y2))< 3 && (wc_x_corr-cluster_pos_corr_hodoX2)<4 && (wc_y_corr-cluster_pos_corr_hodoY2)< 4  && abs( (position_X1-position_X2))<1.5 &&abs((position_Y1-position_Y2))<1.5";

  FitStruct fs950 = getFitResults( outputdir, tree950,  "950HV" ,  500, 140000000, cut, whatToProject );
  float adcEnergyC950 = fs950.mean;  float energyC950= beamEnergy950[0];
 

  /////////200 GeV "HV SCAN" /////////////////////////////////////////
  std::vector<std::string> runs200GeV; 
  std::vector<float> beamEnergy;
  runs200GeV.push_back("487"); //600 V
  beamEnergy.push_back(200000);
  runs200GeV.push_back("428"); // 700 V
  beamEnergy.push_back(200000);
  runs200GeV.push_back("377"); //950 V 
  beamEnergy950.push_back(200000);

  std::vector<std::string> simulation_noBGO; 
 std::vector<float> beamEnergySimulation_noBGO; 

  ///////SIMULATION w/o BGO////////////
  simulation_noBGO.push_back("Simulation50");
  simulation_noBGO.push_back("Simulation20");
  simulation_noBGO.push_back("Simulation100");
  simulation_noBGO.push_back("Simulation150");
  simulation_noBGO.push_back("Simulation200");
   


  ///////SIMULATION////////////
  simulation.push_back("Simulation50");
  beamEnergySimulation.push_back(50000.);
  //  simulation.push_back("Simulation10");
  //  beamEnergySimulation.push_back(10000.);
  //  simulation.push_back("Simulation15");
  //  beamEnergySimulation.push_back(15000.);
  simulation.push_back("Simulation20");
  beamEnergySimulation.push_back(20000.);
  simulation.push_back("Simulation100");
  beamEnergySimulation.push_back(100000.);
  simulation.push_back("Simulation150");
  beamEnergySimulation.push_back(150000.);
  simulation.push_back("Simulation200");
  beamEnergySimulation.push_back(200000.);
  
  TFile* energyfileS = TFile::Open("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/H4UseThis1_5mmGap10mmBeam/Reco_Simulation50.root" );
  //  TFile* energyfileS = TFile::Open("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/H4UseThis1_5mmGap10mmBeam/Reco_Simulation50.root" );
  TTree* energytreeS = (TTree*)energyfileS->Get("recoTree");
 
  std::string SimCut = "abs(xPos)< 3 && abs(yPos)< 3";
  FitStruct fsSim= getFitResults( outputdir, energytreeS , simulation[0] ,  20000. , 350000, SimCut , whatToProjectSim );
  float energyS = fsSim.mean;


  //Taking shit out of the while cond... fFs//
  TCanvas* c1 = new TCanvas( "c1", "", 600, 600 );
  c1->cd();

  TH2D* h2_axes_all2 = new TH2D( "axes_all2", "", 10, 10000 ,900000, 10, 0, 9 );
  h2_axes_all2->SetXTitle("CeF_{3} Response [ADC]");
  h2_axes_all2->SetYTitle("Energy Resolution [%]");


  TF1 *fun_rere20= new TF1("fun_rere20",  "sqrt([0]*[0]/(x*x) +[1]*[1])",  20000, 700000);
  fun_rere20->SetParameter(0, 500000); fun_rere20->SetParameter(1,2);
  //  gr_reso_vs_resp_20GeV->Fit(fun_rere20,"RN");
  fun_rere20->SetLineColor(kBlue+3);

  TF1 *fun_rere50= new TF1("fun_rere50", "sqrt([0]*[0]/(x*x)+[1]*[1]) ",  20000, 700000);
  fun_rere50->SetParameter(0, 500000);
  fun_rere50->SetParameter(1,1.0); 
  //  gr_reso_vs_resp_50GeV->Fit(fun_rere50,"RN");
  fun_rere50->SetLineColor(kAzure);

  TF1 *fun_rere100= new TF1("fun_rere100","sqrt([0]*[0]/(x*x)+[1]*[1])" ,  20000, 700000);
  fun_rere100->SetParameter(0,  500000 );
  fun_rere100->SetParameter(1, 1);
  // gr_reso_vs_resp_100GeV->Fit(fun_rere100,"RN");
  fun_rere100->SetLineColor(kAzure+7);

  TF1 *fun_rere150= new TF1("fun_rere150", "sqrt([0]*[0]/(x*x) +[1]*[1])",  20000, 700000);
  fun_rere150->SetParameter(0, 500000);
  fun_rere150->SetParameter(1, 1.0);
  //  gr_reso_vs_resp_150GeV->Fit(fun_rere150,"RN");
  fun_rere150->SetLineColor(kCyan-6);
  //  fun_rere150->SetLineColor(kCyan-9);

  TF1 *fun_rere200= new TF1("fun_rere200", "sqrt([0]*[0]/(x*x) +[1]*[1])",  20000, 700000);
  fun_rere200->SetParameter(0, 470000);
  fun_rere200->SetParameter(1, 1.0);
  // gr_reso_vs_resp_200GeV->Fit(fun_rere200,"RN");
  fun_rere200->SetLineColor(kBlack);




  
  double globalChiSquare = 10;
  bool doTheBreak=0;
  //  double fermi = 0.0;
  //  double fermi_50 = 0.0;

  // while (globalChiSquare > 1.0){

  //NO, this is not beautiful, but it works and I don't want to spend an hour making it nicer
  for(int i=20; i<50; i++){
    for(int j=17; j < i; j++){
      double fermi = i/100.;
      double fermi_50 = j/100.;

      TGraphErrors* gr_reso_vs_HV = new TGraphErrors(0);
      TGraphErrors* gr_resp_vs_HV = new TGraphErrors(0);
  
      TGraphErrors* gr_reso_vs_HV_50GeV = new TGraphErrors(0);
      TGraphErrors* gr_resp_vs_HV_50GeV = new TGraphErrors(0);
  
      TGraphErrors* gr_reso_vs_HV_100GeV = new TGraphErrors(0);
      TGraphErrors* gr_resp_vs_HV_100GeV = new TGraphErrors(0);
  
      TGraphErrors* gr_reso_vs_HV_150GeV = new TGraphErrors(0);
      TGraphErrors* gr_resp_vs_HV_150GeV = new TGraphErrors(0);
  
      TGraphErrors* gr_resp_vs_energy = new TGraphErrors(0);
      TGraphErrors* gr_reso_vs_energy = new TGraphErrors(0);
      TGraphErrors* gr_dev = new TGraphErrors(0);
  
      TGraphErrors* gr_resp_vs_energy_simul_noBGO = new TGraphErrors(0);
      TGraphErrors* gr_reso_vs_energy_simul_noBGO = new TGraphErrors(0);
 
      TGraphErrors* gr_resp_vs_energy_simul = new TGraphErrors(0);
      TGraphErrors* gr_reso_vs_energy_simul = new TGraphErrors(0);
      TGraphErrors* gr_dev_simul = new TGraphErrors(0);

      TGraphErrors* gr_resp_vs_energy_simul_ideal = new TGraphErrors(0);
      TGraphErrors* gr_reso_vs_energy_simul_ideal = new TGraphErrors(0);
      TGraphErrors* gr_dev_simul_ideal = new TGraphErrors(0);
  
      TGraphErrors* gr_reso_vs_resp_20GeV = new TGraphErrors(0);
      TGraphErrors* gr_reso_vs_resp_50GeV = new TGraphErrors(0);
      TGraphErrors* gr_reso_vs_resp_100GeV = new TGraphErrors(0);
      TGraphErrors* gr_reso_vs_resp_150GeV = new TGraphErrors(0);
      TGraphErrors* gr_reso_vs_resp_200GeV = new TGraphErrors(0);

      TGraphErrors* gr_reso_vs_energy_stat = new TGraphErrors(0);
  
      TGraphErrors* gr_reso_corr = new TGraphErrors(0);

      float xMax = 205000;
      float yMax = 2350000;

  


      for( unsigned i=0; i<runs950.size(); ++i ) {
    
	TFile* file = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%s.root", tag.c_str(),  runs950[i].c_str() ));
	TTree* tree = (TTree*)file->Get("recoTree");

  

	FitStruct fs  = getFitResults( outputdir, tree ,  runs950[i].c_str() , 500, 140000000, cut, whatToProject );
	if(i==1) FitStruct fs  = getFitResults( outputdir, tree ,  Form("%d",i) , 50, 9900000, cut, whatToProject );
 

	float energy = beamEnergy950[i];
	float energyErr = 0.005*energy;
	float mean = fs.mean;
	float meanErr = fs.mean_err;
	float reso = fs.reso;
	float resoErr = fs.reso_err;
    
	std::cout << "mean at "<< energy << " = " << mean << std::endl;
	gr_resp_vs_energy->SetPoint( i, energy/1000., mean /adcEnergyC950 * adcEnergyC950);
	gr_resp_vs_energy->SetPointError( i, energyErr/1000., meanErr);
    
	gr_reso_vs_energy->SetPoint( i, energy/1000., reso );
	gr_reso_vs_energy->SetPointError( i, energyErr/1000., resoErr );
    
	gr_dev->SetPoint(i, energy/1000. , mean / energy / ( adcEnergyC950 / energyC950) );
	gr_dev->SetPointError(i, 0, meanErr / energy / ( adcEnergyC950 / energyC950) );

     
	if(i == 0) {   
	  gr_reso_vs_resp_50GeV->SetPoint(465-459+1, mean, reso);   
	  gr_reso_vs_resp_50GeV->SetPointError(465-459+1, meanErr, sqrt(resoErr* resoErr +fermi_50*fermi_50) );
	  //      gr_reso_vs_resp_50GeV->SetPointError(465-459+1, meanErr, sqrt(resoErr* resoErr +fermi*fermi) );
	  //    gr_reso_vs_resp_50GeV->SetPointError(465-459+1, meanErr, resoErr );
	}
	/*
	  if(i == 1) {
	  gr_reso_vs_resp_20GeV->SetPoint( 453 - 432, mean, reso);
	  gr_reso_vs_resp_20GeV->SetPointError( 453 - 432, meanErr, sqrt(resoErr* resoErr +fermi *fermi) );	  }
	*/

	/*
	//Anyway out of range//
	if(i == 2) {
	gr_reso_vs_resp_100GeV->SetPoint(472 - 466, mean, reso);
	gr_reso_vs_resp_100GeV->SetPointError(472 - 466, meanErr, resoErr);	}

	if(i == 3){
	gr_reso_vs_resp_150GeV->SetPoint(480 - 474, mean, reso);
	gr_reso_vs_resp_150GeV->SetPointError(480 - 474, meanErr, resoErr);	}

	//  gr_reso_corr->SetPoint(i, energy/1000. , reso);
	// gr_reso_corr->SetPointError(i, energyErr/1000. , resoErr); // commenting it out
	*/

	delete file;
      }




      for( unsigned i=0; i<runs200GeV.size(); ++i ) {
	TFile* file = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%s.root", tag.c_str(),  runs200GeV[i].c_str() ));
	TTree* tree = (TTree*)file->Get("recoTree");
	FitStruct fs  = getFitResults( outputdir, tree ,  runs200GeV[i].c_str() , 500, 14000000, cut, whatToProject );
    
	float energy = beamEnergy[i];
	//    float energyErr = 0.005*energy;
	float mean = fs.mean;
	float meanErr = fs.mean_err;
	float reso = fs.reso;
	float resoErr = fs.reso_err;

	gr_reso_vs_resp_200GeV->SetPoint(i , mean , reso);
	gr_reso_vs_resp_200GeV->SetPointError(i, meanErr,  resoErr);

	delete file;
      }





      for( unsigned i=0; i<simulation.size(); ++i ) {
	///////////////// (1x1) Calorimeter //////////////////////////////
	TFile* fileS = TFile::Open(Form("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/H4UseThis1_5mmGap10mmBeam/Reco_%s.root", simulation[i].c_str()));
   
	TTree* treeS = (TTree*)fileS->Get("recoTree");
   
	FitStruct fsSim= getFitResults( outputdir, treeS , simulation[i] , beamEnergySimulation[i]/10. ,  beamEnergySimulation[i]*1.5, SimCut , whatToProjectSim );
	float energySim = fsSim.mean;
   
	FitStruct fsSim_ps = addPhotoStatistics( fsSim );

	gr_resp_vs_energy_simul->SetPoint(i,beamEnergySimulation[i]/1000., energySim*adcEnergyC950/energyS );
	gr_resp_vs_energy_simul->SetPointError( i,0, fsSim.mean_err * adcEnergyC950/energyS);
   
	gr_reso_vs_energy_simul->SetPoint( i, beamEnergySimulation[i]/1000., fsSim_ps.reso );
	gr_reso_vs_energy_simul->SetPointError( i,0,  fsSim_ps.reso_err );
   
	gr_dev_simul->SetPoint(i, beamEnergySimulation[i]/1000. , fsSim.mean/beamEnergySimulation[i]/(energyS/energyC950) );
	gr_dev_simul->SetPointError(i, 0, fsSim.mean_err/beamEnergySimulation[i]/(energyS/energyC950) );

	delete fileS;
      }




      ///////////////// (5x5) Calo  ///////////////////  
      for( unsigned i=0; i<simulation.size(); ++i ) {
	TFile* fileS = TFile::Open(Form("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/H4Ideal/Reco_%s.root", simulation[i].c_str()));
   
	TTree* treeS = (TTree*)fileS->Get("recoTree");
   
	FitStruct fsSim= getFitResults( outputdir, treeS , simulation[i] , beamEnergySimulation[i]/10. ,  beamEnergySimulation[i]*1.5, SimCut , whatToProjectSim );
	float energySim = fsSim.mean;
   
	FitStruct fsSim_ps = addPhotoStatistics( fsSim );

	gr_resp_vs_energy_simul_ideal->SetPoint(i,beamEnergySimulation[i]/1000., energySim*adcEnergyC950/energyS );
	gr_resp_vs_energy_simul_ideal->SetPointError( i,0, fsSim.mean_err * adcEnergyC950/energyS);
   
	gr_reso_vs_energy_simul_ideal->SetPoint( i, beamEnergySimulation[i]/1000., fsSim_ps.reso );
	gr_reso_vs_energy_simul_ideal->SetPointError( i,0,  fsSim_ps.reso_err );
   
	gr_dev_simul_ideal->SetPoint(i, beamEnergySimulation[i]/1000. , fsSim.mean/beamEnergySimulation[i]/(energyS/energyC950) );
	gr_dev_simul_ideal->SetPointError(i, 0, fsSim.mean_err/beamEnergySimulation[i]/(energyS/energyC950) );

	delete fileS;
      }
 
 


      //Resolution vs HV at 20GeV
      for( unsigned i=432; i < 454; ++i ) { 
	TFile* file = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%d.root", tag.c_str(),  i ));
	TTree* tree = (TTree*)file->Get("recoTree");

	float HV;
	tree->SetBranchAddress("HVCeF3",&HV);
	tree->GetEntry(5);
	std::cout << HV << std::endl;
	FitStruct fs  = getFitResults( outputdir, tree ,  Form("%d",i) , 50, 9900000, cut, whatToProject );
	if(i==449) fs  = getFitResults( outputdir, tree ,  Form("%d",i) , 5000, 9900000, cut, whatToProject );
  
	if(fs.reso > 10. || fs.chi2 > 1.1) continue;
	float mean = fs.mean;
	float meanErr = fs.mean_err; 
	float reso = fs.reso;
	float resoErr = fs.reso_err;
   
	gr_reso_vs_HV->SetPoint(i, HV, reso);
	gr_reso_vs_HV->SetPointError(i, 0, resoErr);
	gr_resp_vs_HV->SetPoint(i, HV, mean);
	gr_resp_vs_HV->SetPointError(i, 0, meanErr);
   
	//  gr_reso_vs_resp_20GeV1->SetPoint(i, mean, reso);
	gr_reso_vs_resp_20GeV->SetPoint(i, mean, reso);
	//   gr_reso_vs_resp_20GeV->SetPointError(i, meanErr, resoErr);
	gr_reso_vs_resp_20GeV->SetPointError(i, meanErr, sqrt(resoErr*resoErr + fermi*fermi));

	delete file;
      }


      //Resolution vs HV at 50GeV
      for( unsigned i=459; i < 466; ++i ) {
	TFile* file = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%d.root", tag.c_str(),  i ));
	TTree* tree = (TTree*)file->Get("recoTree");
	float HV;
	tree->SetBranchAddress("HVCeF3",&HV);
	tree->GetEntry(5);
	std::cout << HV << std::endl;

	FitStruct fs  = getFitResults( outputdir, tree ,  Form("%d",i) ,250, 70000000, cut, whatToProject ); 
	if(i==465) fs  = getFitResults( outputdir, tree ,  Form("%d",i) , 10000,  90000000, cut, whatToProject );
	if(fs.reso > 10. || fs.chi2 > 1.1) continue; 
	float mean = fs.mean;
	float meanErr = fs.mean_err; 
	float reso = fs.reso;
	float resoErr = fs.reso_err;
   
	gr_reso_vs_HV_50GeV->SetPoint(i, HV, reso);
	gr_reso_vs_HV_50GeV->SetPointError(i, 0, resoErr);
	gr_resp_vs_HV_50GeV->SetPoint(i, HV, mean);
	gr_resp_vs_HV_50GeV->SetPointError(i, 0, meanErr);

	//   gr_reso_vs_resp_50GeV1->SetPoint(i, mean, reso);

	gr_reso_vs_resp_50GeV->SetPoint(i, mean, reso);
	/*
	  if( resoErr > fermi){
	  gr_reso_vs_resp_50GeV->SetPointError(i, meanErr, resoErr);
	  }else{
	  gr_reso_vs_resp_50GeV->SetPointError(i, meanErr, sqrt(resoErr*resoErr + fermi*fermi) );
	  }
	*/
	gr_reso_vs_resp_50GeV->SetPointError(i, meanErr, sqrt(resoErr*resoErr + fermi_50*fermi_50) );
	//   gr_reso_vs_resp_50GeV->SetPointError(i, meanErr, sqrt(resoErr*resoErr + fermi*fermi) );
	//gr_reso_vs_resp_50GeV->SetPointError(i, meanErr, resoErr);

	delete file;
      }


      //Resolution vs HV at 100GeV
      for( unsigned i=466; i < 473; ++i ) {
	TFile* file = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%d.root", tag.c_str(),  i ));
	TTree* tree = (TTree*)file->Get("recoTree");
	float HV;
	tree->SetBranchAddress("HVCeF3",&HV);
	tree->GetEntry(5);
	std::cout << HV << std::endl;
	FitStruct fs  = getFitResults( outputdir, tree ,  Form("%d",i) , 25, 19000000, cut, whatToProject );  
	if(fs.reso < 0 || fs.reso > 10. || fs.chi2 > 1.1) continue;
	float mean = fs.mean;
	float meanErr = fs.mean_err; 
	float reso = fs.reso;
	float resoErr = fs.reso_err;
   
	gr_reso_vs_HV_100GeV->SetPoint(i, HV, reso);
	gr_reso_vs_HV_100GeV->SetPointError(i, 0, resoErr);
	gr_resp_vs_HV_100GeV->SetPoint(i, HV, mean);
	gr_resp_vs_HV_100GeV->SetPointError(i, 0, meanErr);

	//   gr_reso_vs_resp_100GeV1->SetPoint(i, mean, reso);
 
	gr_reso_vs_resp_100GeV->SetPoint(i, mean, reso);
	gr_reso_vs_resp_100GeV->SetPointError(i, meanErr, resoErr );
	/*
	  if( resoErr > fermi){
	  gr_reso_vs_resp_100GeV->SetPointError(i, meanErr, resoErr);
	  }else{
	  //gr_reso_vs_resp_100GeV->SetPointError(i, meanErr, fermi/5. );
	  gr_reso_vs_resp_100GeV->SetPointError(i, meanErr, sqrt( resoErr* resoErr + fermi*fermi/25.) );
	  }
	*/

	delete file;
      }



      //Resolution vs HV at 150GeV
      for( unsigned i=474; i < 481; ++i ) {
	TFile* file = TFile::Open(Form("analysisTrees_chaInt_04_12_%s/Reco_%d.root", tag.c_str(),  i ));
	TTree* tree = (TTree*)file->Get("recoTree");
	float HV;
	tree->SetBranchAddress("HVCeF3",&HV);
	tree->GetEntry(5);
	std::cout << HV << std::endl;
	FitStruct fs  = getFitResults( outputdir, tree ,  Form("%d",i) ,50 , 140000000, cut, whatToProject );  
	if(fs.reso > 10. || fs.chi2 > 1.1) continue;
	float mean = fs.mean;
	float meanErr = fs.mean_err; 
	float reso = fs.reso;
	float resoErr = fs.reso_err;
   
	gr_reso_vs_HV_150GeV->SetPoint(i, HV, reso);
	gr_reso_vs_HV_150GeV->SetPointError(i, 0, resoErr);
	gr_resp_vs_HV_150GeV->SetPoint(i, HV, mean);
	gr_resp_vs_HV_150GeV->SetPointError(i, 0, meanErr);

	//   gr_reso_vs_resp_150GeV1->SetPoint(i, mean, reso);
 
	gr_reso_vs_resp_150GeV->SetPoint(i, mean, reso);
	// gr_reso_vs_resp_150GeV->SetPointError(i, meanErr, sqrt(resoErr*resoErr + fermi*fermi) );
	gr_reso_vs_resp_150GeV->SetPointError(i, meanErr, resoErr);

	delete file;
      }



      //SIMULATION WIHTOUT BGO
      for( unsigned i=0; i<simulation_noBGO.size(); ++i ) {
	///////////////// (1x1) Calorimeter //////////////////////////////
	TFile* fileS = TFile::Open(Form("/home/myriam/BTFAnalysis/PositionAnalysis/OriginalSimulationData/noBGO/Reco_%s.root", simulation[i].c_str()));
	TTree* treeS = (TTree*)fileS->Get("recoTree");
   
	FitStruct fsSim= getFitResults( outputdir, treeS , simulation[i] , beamEnergySimulation[i]/10. ,  beamEnergySimulation[i]*1.5, SimCut , whatToProjectSim );
	float energySim = fsSim.mean;
   
	FitStruct fsSim_ps = addPhotoStatistics( fsSim );

	gr_resp_vs_energy_simul_noBGO->SetPoint(i,beamEnergySimulation[i]/1000., energySim*adcEnergyC950/energyS );
	gr_resp_vs_energy_simul_noBGO->SetPointError( i,0, fsSim.mean_err * adcEnergyC950/energyS);
   
	gr_reso_vs_energy_simul_noBGO->SetPoint( i, beamEnergySimulation[i]/1000., fsSim_ps.reso );
	gr_reso_vs_energy_simul_noBGO->SetPointError( i,0,  fsSim_ps.reso_err );
   
	delete fileS;
      }



 
      //PLOTTING THINGS///////////////////

      TH2D* h2_axes = new TH2D( "axes", "", 10, 0., xMax/1000.+1, 10, 0., yMax );
      h2_axes->SetXTitle("Beam Energy [GeV]");
      h2_axes->SetYTitle("CeF_{3} Response [ADC]");
      h2_axes->Draw("");

      gr_resp_vs_energy->SetMarkerStyle(20);
      gr_resp_vs_energy->SetMarkerSize(1.6);
      gr_resp_vs_energy->SetMarkerColor(46);
      gr_resp_vs_energy->Draw("p same");

      TF1* f1_line = new TF1("line", "[0] + [1]*x", 0., xMax/1000.+5 );
      gr_resp_vs_energy->Fit(f1_line,"RN");
      f1_line->SetLineWidth(1.);
      f1_line->SetLineColor(46);
      // f1_line->Draw("L same");

      gr_resp_vs_energy_simul->SetMarkerStyle(21);
      gr_resp_vs_energy_simul->SetMarkerSize(1.3);
      gr_resp_vs_energy_simul->SetMarkerColor(kBlack);
      gr_resp_vs_energy_simul->Draw("p same");

      TF1* f1_lines = new TF1("lines", "[0] + [1]*x", 0., xMax/1000.+5 );
      gr_resp_vs_energy_simul->Fit(f1_lines,"RN");
      f1_lines->SetLineWidth(1.);
      f1_lines->SetLineColor(kBlack);
      f1_lines->Draw("L same");

      TLegend* leg0 = new TLegend(0.6, 0.4-3*0.06, 0.9, 0.4);
      leg0->SetTextSize(0.038);
      leg0->AddEntry(gr_resp_vs_energy,"Data ","p");
      //leg0->AddEntry(f1_line,Form("Offset = %.1f\n #pm %.2f\n",f1_line->GetParameter(0), f1_line->GetParError(0) ),"L");
      //leg0->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.0f\n / %d",f1_line->GetChisquare(), f1_line->GetNDF() ), "");
      leg0->AddEntry(f1_lines,"MC","L");

      leg0->SetFillColor(0);
      leg0->Draw("same");

      TPaveText* label_top2 = new TPaveText();
      label_top2 = DrawTools::getLabelTop("Electron Beam");
      label_top2->Draw("same");

      c1->SaveAs( Form( "%s/resp_vs_energy.eps", outputdir.c_str() ) );
      c1->SaveAs( Form( "%s/resp_vs_energy.png", outputdir.c_str() ) );
    
      c1->Clear();




      gStyle->SetOptFit(0);
  
      gPad->SetLeftMargin(0.11);
      gPad->SetRightMargin(0.09);
      gPad->SetBottomMargin(0.15);
      gStyle->SetTitleYOffset(1.01); // => 1.15 if exponents
  
      c1->Clear();

      //////////////////Reso vs Response for correction///////////////////////////

      h2_axes_all2->Draw("");

      gr_reso_vs_resp_150GeV->SetMarkerStyle(34);
      gr_reso_vs_resp_150GeV->SetMarkerSize(1.6);
      gr_reso_vs_resp_150GeV->SetMarkerColor(kCyan-6);
      //  gr_reso_vs_resp_150GeV->SetMarkerColor(kCyan-9);

      gr_reso_vs_resp_100GeV->SetMarkerStyle(29);
      gr_reso_vs_resp_100GeV->SetMarkerSize(1.6);
      gr_reso_vs_resp_100GeV->SetMarkerColor(kAzure+7);

      gr_reso_vs_resp_50GeV->SetMarkerStyle(20);
      gr_reso_vs_resp_50GeV->SetMarkerSize(1.4);
      gr_reso_vs_resp_50GeV->SetMarkerColor(kAzure);
 
      gr_reso_vs_resp_20GeV->SetMarkerStyle(21);
      gr_reso_vs_resp_20GeV->SetMarkerSize(1.3);
      gr_reso_vs_resp_20GeV->SetMarkerColor(kBlue+3);

      gr_reso_vs_resp_200GeV->SetMarkerStyle(20);
      gr_reso_vs_resp_200GeV->SetMarkerSize(1.6);
      gr_reso_vs_resp_200GeV->SetMarkerColor(kBlack);





      //SIMULTANEOUUUS//
      ROOT::Math::WrappedMultiTF1 wf20GeV(*fun_rere20,1);
      ROOT::Math::WrappedMultiTF1 wf50GeV(*fun_rere50,1);
      ROOT::Math::WrappedMultiTF1 wf100GeV(*fun_rere100,1);
      ROOT::Math::WrappedMultiTF1 wf150GeV(*fun_rere150,1);
      ROOT::Math::WrappedMultiTF1 wf200GeV(*fun_rere200,1);
      ROOT::Fit::DataOptions opt; 

      // set the data range
      ROOT::Fit::DataRange range20GeV; 
      // range20GeV.SetRange(  200, 2000000 ); //Full range//
      // range20GeV.SetRange(  60000, 700000 );//use this
      range20GeV.SetRange(  90, 700000 );//use this (nice ones)
      // range20GeV.SetRange(  90000, 600000 );
      ROOT::Fit::BinData data20GeV(opt,range20GeV); 
      data20GeV.Opt().fCoordErrors = false; 
      ROOT::Fit::FillData(data20GeV, gr_reso_vs_resp_20GeV);
  
      ROOT::Fit::DataRange range50GeV; 
      // range50GeV.SetRange(  200, 2000000 ); //full range
      // range50GeV.SetRange(  30000, 2000000 );
      range50GeV.SetRange(  10, 700000 ); //use this (nice fits)
      //  range50GeV.SetRange(  30000, 700000 ); //use this
      // range50GeV.SetRange(  100000, 800000 );
      ROOT::Fit::BinData data50GeV(opt,range50GeV);
      data50GeV.Opt().fCoordErrors = false; 
      ROOT::Fit::FillData(data50GeV, gr_reso_vs_resp_50GeV );
  
      ROOT::Fit::DataRange range100GeV; 
      // range100GeV.SetRange( 200, 2000000 ); //full range
      //  range100GeV.SetRange(  60000, 2000000 );
      range100GeV.SetRange(  20, 700000 );//use this nice
      //range100GeV.SetRange(  60000, 700000 );//use this
      //  range100GeV.SetRange(  60000, 700000 ); //one more point
      //  range100GeV.SetRange(  200000, 1000000 );
      ROOT::Fit::BinData data100GeV(opt,range100GeV); 
      data100GeV.Opt().fCoordErrors = false; 
      ROOT::Fit::FillData(data100GeV, gr_reso_vs_resp_100GeV );
  
      ROOT::Fit::DataRange range150GeV; 
      // range150GeV.SetRange(  200, 2000000 );//full range
      //  range150GeV.SetRange(  30000, 1200000 ); //two more points
      range150GeV.SetRange(  20, 700000 );//use this
      // range150GeV.SetRange(  100000, 700000 );//use this
      //  range150GeV.SetRange(  100000, 700000 );//one more point
      // range150GeV.SetRange(  200000, 1200000 ); //original
      ROOT::Fit::BinData data150GeV(opt,range150GeV); 
      data150GeV.Opt().fCoordErrors = false; 
      ROOT::Fit::FillData(data150GeV, gr_reso_vs_resp_150GeV );
  
      ROOT::Fit::DataRange range200GeV; 
      //  range200GeV.SetRange(  1000000, 2000000 );
      range200GeV.SetRange(  400000, 1000000 );
      ROOT::Fit::BinData data200GeV(opt,range200GeV);
      data200GeV.Opt().fCoordErrors = false; 
      ROOT::Fit::FillData(data200GeV, gr_reso_vs_resp_200GeV );


      ROOT::Fit::Chi2Function chi2_20GeV(data20GeV, wf20GeV);
      ROOT::Fit::Chi2Function chi2_50GeV(data50GeV, wf50GeV);
      ROOT::Fit::Chi2Function chi2_100GeV(data100GeV, wf100GeV);
      ROOT::Fit::Chi2Function chi2_150GeV(data150GeV, wf150GeV);
      ROOT::Fit::Chi2Function chi2_200GeV(data200GeV, wf200GeV);
  

      ROOT::Fit::Fitter fitter;
      const int Npar = 5; 
      double par0[Npar] = { 460000, 2.9 , 1.9, 1.23 , 1.1 };

      //  const int Npar = 6; 
      //  double par0[Npar] = { 460000, 2.9 , 1.9, 1.23 , 1.2,  0.87};

      // create before the parameter settings in order to fix or set range on them
      fitter.Config().SetParamsSettings(Npar,par0);
      fitter.Config().ParSettings(0).SetLimits(100000, 650000);
      // fitter.Config().ParSettings(0).SetLimits(450000, 650000);
      fitter.Config().ParSettings(1).SetLimits( 2.5,3.5); 
      fitter.Config().ParSettings(2).SetLimits(1,3);
      fitter.Config().ParSettings(3).SetLimits(0.5,2);
      fitter.Config().ParSettings(4).SetLimits(0.5,2);
  
      fitter.Config().MinimizerOptions().SetPrintLevel(3);

      // fitter.Config().MinimizerOptions().SetMinimizerType("GSLMultiMin");
      // fitter.Config().MinimizerOptions().SetMinimizerType("GSLMultiFit");
      // fitter.Config().MinimizerOptions().SetMinimizerType("Genetic");
      // fitter.Config().MinimizerOptions().SetMinimizerType("Fumili2");
      fitter.Config().MinimizerOptions().SetMinimizerType("Minuit2");
  
      //fit FCN function directly (specify optionally data size and flag to indicate that is a chi2 fit)
      GlobalChi2 globalChi2(Npar, (data20GeV.Size()+data50GeV.Size()+data100GeV.Size()+data150GeV.Size() ) ,chi2_20GeV, chi2_50GeV, chi2_100GeV, chi2_150GeV );

      fitter.FitFCN(globalChi2, 0, data20GeV.Size() + data50GeV.Size() + data100GeV.Size() + data150GeV.Size() , true);
      //  fitter.FitFCN(globalChi2, 0, data20GeV.Size() + data50GeV.Size() + data100GeV.Size() + data150GeV.Size() , true);
      ROOT::Fit::FitResult result = fitter.Result();

      std::cout << std::endl;
      std::cout << std::endl;
      result.Print(std::cout);
      std::cout << std::endl;
      std::cout << std::endl;  
      std::cout << std::endl;
      std::cout << result.Chi2() / result.Ndf() << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "----------- END of first fit results ---------------- " << std::endl;



      //  ipar200GeV[0] = 0.9 ;
      fun_rere20->SetFitResult( result, ipar20GeV );
      fun_rere20->Draw("L same");
      fun_rere50->SetFitResult( result, ipar50GeV );
      fun_rere50->Draw("L same");
      fun_rere100->SetFitResult( result, ipar100GeV );
      fun_rere100->Draw("L same");
      fun_rere150->SetFitResult( result, ipar150GeV );
      fun_rere150->Draw("L same");
      //  fun_rere200->SetFitResult( result, ipar200GeV );
      // fun_rere200->Draw("L same");

      gr_reso_vs_resp_20GeV ->Draw("p same");
      gr_reso_vs_resp_50GeV ->Draw("p same");
      gr_reso_vs_resp_100GeV->Draw("p same");
      gr_reso_vs_resp_150GeV->Draw("p same");
  


      TLegend* leggie_all_resp = new TLegend(0.45, 0.92 - 0.06 *8 , 0.85, 0.92);
      leggie_all_resp->SetTextSize(0.038);
      leggie_all_resp->AddEntry(gr_reso_vs_resp_20GeV,"Data 20GeV","p");
      leggie_all_resp->AddEntry( fun_rere20,Form("N =  %.00f\n #pm %.0f\n kADC",(fun_rere20->GetParameter(0))/1000., (fun_rere20->GetParError(0))/1000. ),"");
      leggie_all_resp->AddEntry( (TObject*)0,Form("C = %.2f\n #pm %.2f\n %s",(fun_rere20->GetParameter(1)), (fun_rere20->GetParError(1)),"%" ),"");
  
      leggie_all_resp->AddEntry(gr_reso_vs_resp_50GeV,"Data 50GeV","p");
      leggie_all_resp->AddEntry( fun_rere50,Form("N =  %.00f\n #pm %.0f\n kADC",(fun_rere50->GetParameter(0))/1000., (fun_rere50->GetParError(0))/1000. ),"");
      leggie_all_resp->AddEntry( (TObject*)0,Form("C = %.2f\n #pm %.2f\n %s",(fun_rere50->GetParameter(1)), (fun_rere50->GetParError(1)) ,"%"),"");
  
      leggie_all_resp->AddEntry(gr_reso_vs_resp_100GeV,"Data 100GeV","p");
      leggie_all_resp->AddEntry( fun_rere100,Form("N =  %.00f\n #pm %.0f\n kADC",(fun_rere100->GetParameter(0))/1000., (fun_rere100->GetParError(0))/1000. ),"");
      leggie_all_resp->AddEntry( (TObject*)0,Form("C = %.2f\n #pm %.2f\n %s",(fun_rere100->GetParameter(1)), (fun_rere100->GetParError(1)),"%" ),"");
  
      leggie_all_resp->AddEntry(gr_reso_vs_resp_150GeV,"Data 150GeV","p");
      leggie_all_resp->AddEntry( fun_rere150,Form("N =  %.00f\n #pm %.0f\n kADC",(fun_rere150->GetParameter(0))/1000., (fun_rere150->GetParError(0))/1000.),"");
      leggie_all_resp->AddEntry( (TObject*)0,Form("C  = %.2f\n #pm %.2f\n %s",(fun_rere150->GetParameter(1)), (fun_rere150->GetParError(1)) ,"%"),"");
      /*
	leggie_all_resp->AddEntry(gr_reso_vs_resp_200GeV,"Data 200GeV","p");
	leggie_all_resp->AddEntry( fun_rere200,Form("N =  %.00f\n #pm %.0f\n ADC",(fun_rere200->GetParameter(0)), (fun_rere200->GetParError(0)) ),"");
	leggie_all_resp->AddEntry( (TObject*)0,Form("C  = %.2f\n #pm %.2f\n %s",(fun_rere200->GetParameter(1)), (fun_rere200->GetParError(1)),"%" ),"");
      */

      leggie_all_resp->AddEntry( (TObject*)0, Form("#chi^{2}/NDF = %.0f/%d = %.3f", result.Chi2() , result.Ndf(), result.Chi2()/ result.Ndf()  ), "" );
      leggie_all_resp->SetFillColor(0);
      leggie_all_resp->Draw("same");
  
  
  
  
      TPaveText* label_top =  new TPaveText(0.3,0.953,0.945,0.975, "brNDC");
      label_top->SetFillColor(kWhite);
      label_top->SetTextSize(0.038);
      label_top->SetTextAlign(31); // align right
      label_top->SetTextFont(62);
      label_top->AddText("Electron Beam");
      label_top->Draw("same");
  
  
      // label_top2->Draw("same");
  
      c1->SaveAs( Form( "%s/reso_vs_resp_all_%d_%d.eps", outputdir.c_str(), int(fermi*100), int(fermi_50*100) ) );
      c1->SaveAs( Form( "%s/reso_vs_resp_all_%d_%d.png", outputdir.c_str(), int(fermi*100), int(fermi_50*100) ) );
      c1->SaveAs( Form( "%s/reso_vs_resp_all_%d_%d.pdf", outputdir.c_str(), int(fermi*100), int(fermi_50*100) ) );




      c1->Clear();

      //NEAT VERSION FOR FRANNY ///////////////////////////////////////////
      TH2D* h2_axes_neat2 = new TH2D( "axes_neat2", "", 10, 10000 , 900000, 10, 0, 7 );
      h2_axes_neat2->SetXTitle("CeF_{3} Response [ADC]");
      h2_axes_neat2->SetYTitle("#sigma_{R} [%]");
      h2_axes_neat2->Draw("");

      fun_rere20->Draw("L same");
      fun_rere50->Draw("L same");
      fun_rere100->Draw("L same");
      fun_rere150->Draw("L same");

      gr_reso_vs_resp_20GeV ->Draw("p same");
      gr_reso_vs_resp_50GeV ->Draw("p same");
      gr_reso_vs_resp_100GeV->Draw("p same");
      gr_reso_vs_resp_150GeV->Draw("p same");
  
      TLegend* leggie_neat_resp = new TLegend(0.65, 0.9 - 0.06 *4 , 0.85, 0.9);
      leggie_neat_resp->SetTextSize(0.038);
      leggie_neat_resp->AddEntry(gr_reso_vs_resp_20GeV,"20GeV","p");
      leggie_neat_resp->AddEntry(gr_reso_vs_resp_50GeV,"50GeV","p");
      leggie_neat_resp->AddEntry(gr_reso_vs_resp_100GeV,"100GeV","p");
      leggie_neat_resp->AddEntry(gr_reso_vs_resp_150GeV,"150GeV","p");
      leggie_neat_resp->SetFillColor(0);
      leggie_neat_resp->Draw("same");
  


      TLine* line = new TLine( 700000, 0  , 700000, 7 );//x1y1x2y2
      line ->SetLineColor(kGray+2);
      line ->SetLineWidth(3);
      line ->SetLineStyle(2);
      //  line->Draw("same"); 

      TPaveText* label_low2 = new TPaveText(0.123,0.175,0.5,0.23, "brNDC");
      label_low2->SetFillColor(kWhite);
      label_low2->SetTextSize(0.038);
      label_low2->SetTextAlign(11); // align right
      label_low2->SetTextFont(62);
      label_low2->AddText( "W-CeF_{3} Single Tower");
      label_low2->Draw("same");
 
      label_top->Draw("same");

      gPad->RedrawAxis();
 
      c1->SaveAs( Form( "%s/reso_vs_resp_all_neat_%d.eps", outputdir.c_str(), int(fermi*100) ) );
      c1->SaveAs( Form( "%s/reso_vs_resp_all_neat_%d.png", outputdir.c_str(), int(fermi*100) ) );
      c1->SaveAs( Form( "%s/reso_vs_resp_all_neat_%d.pdf", outputdir.c_str(), int(fermi*100) ) );
  
    c1->SaveAs( Form( "%s/reso_vs_resp_all_neat_%d_%d.pdf", outputdir.c_str(), int(fermi*100), int(fermi_50*100) ) );
  
    c1->SaveAs( Form( "%s/reso_vs_resp_all_neat_%d_%d.png", outputdir.c_str(), int(fermi*100), int(fermi_50*100) ) );



      c1->Clear();


  
      double reso200 = fun_rere200->GetParameter(1);
      double reso200_err = fun_rere200->GetParError(1);
      double reso150 = fun_rere150->GetParameter(1);
      double reso150_err = fun_rere150->GetParError(1);
      double reso100 = fun_rere100->GetParameter(1);
      double reso100_err = fun_rere100->GetParError(1);
      double reso50 = fun_rere50->GetParameter(1);
      double reso50_err = fun_rere50->GetParError(1);
      double reso20 = fun_rere20->GetParameter(1);
      double reso20_err = fun_rere20->GetParError(1);
  
  













  //Systematic errors: Set Value 700k, vary first up to 900k
  // set the data range
  ROOT::Fit::DataRange range20GeV2; 
  range20GeV2.SetRange(  20, 900000 );
  ROOT::Fit::BinData data20GeV2(opt,range20GeV2); 
  data20GeV2.Opt().fCoordErrors = false; 
  ROOT::Fit::FillData(data20GeV2, gr_reso_vs_resp_20GeV);
  
  ROOT::Fit::DataRange range50GeV2; 
  range50GeV2.SetRange(  10, 900000 ); 
  ROOT::Fit::BinData data50GeV2(opt,range50GeV2);
  data50GeV2.Opt().fCoordErrors = false; 
  ROOT::Fit::FillData(data50GeV2, gr_reso_vs_resp_50GeV );
  
  ROOT::Fit::DataRange range100GeV2; 
  range100GeV2.SetRange(  20, 900000 );
  ROOT::Fit::BinData data100GeV2(opt,range100GeV2); 
  data100GeV2.Opt().fCoordErrors = false; 
  ROOT::Fit::FillData(data100GeV2, gr_reso_vs_resp_100GeV );
  
  ROOT::Fit::DataRange range150GeV2; 
  range150GeV2.SetRange(  20, 900000 );//one more point diff
  ROOT::Fit::BinData data150GeV2(opt,range150GeV2); 
  data150GeV2.Opt().fCoordErrors = false; 
  ROOT::Fit::FillData(data150GeV2, gr_reso_vs_resp_150GeV );
  
  ROOT::Fit::DataRange range200GeV2; 
  range200GeV2.SetRange(  400000, 1000000 );
  ROOT::Fit::BinData data200GeV2(opt,range200GeV2);
  data200GeV2.Opt().fCoordErrors = false; 
  ROOT::Fit::FillData(data200GeV2, gr_reso_vs_resp_200GeV );


  ROOT::Fit::Chi2Function chi2_20GeV2(data20GeV2, wf20GeV);
  ROOT::Fit::Chi2Function chi2_50GeV2(data50GeV2, wf50GeV);
  ROOT::Fit::Chi2Function chi2_100GeV2(data100GeV2, wf100GeV);
  ROOT::Fit::Chi2Function chi2_150GeV2(data150GeV2, wf150GeV);
  ROOT::Fit::Chi2Function chi2_200GeV2(data200GeV2, wf200GeV);
  

  ROOT::Fit::Fitter fitter2;
  fitter2.Config().SetParamsSettings(Npar,par0);
  fitter2.Config().MinimizerOptions().SetPrintLevel(3);
  fitter2.Config().ParSettings(0).SetLimits(450000, 650000);
  fitter2.Config().ParSettings(1).SetLimits( 2.5,3.5); 
  fitter2.Config().ParSettings(2).SetLimits(1,3);
  fitter2.Config().ParSettings(3).SetLimits(0.5,2);
  fitter2.Config().ParSettings(4).SetLimits(0.5,2);
  
  fitter2.Config().MinimizerOptions().SetPrintLevel(3);

  fitter2.Config().MinimizerOptions().SetMinimizerType("Minuit2");
  //  fitter2.Config().MinimizerOptions().SetMinimizerType("GSLMultiFit");

  GlobalChi2 globalChi22(Npar, (data20GeV2.Size()+data50GeV2.Size()+data100GeV2.Size()+data150GeV2.Size() ) ,chi2_20GeV2, chi2_50GeV2, chi2_100GeV2, chi2_150GeV2);

  fitter2.FitFCN(globalChi22, 0, data20GeV2.Size() + data50GeV2.Size() + data100GeV2.Size() + data150GeV2.Size() , true);
  ROOT::Fit::FitResult result2 = fitter2.Result();


  std::cout << std::endl;  std::cout << std::endl;
  std::cout << "------------- Second Result -------------- " << std::endl;
  result2.Print(std::cout);
  std::cout << std::endl;  std::cout << std::endl;


  fun_rere20->SetFitResult( result2, ipar20GeV );
  fun_rere50->SetFitResult( result2, ipar50GeV );
  fun_rere100->SetFitResult( result2, ipar100GeV );
  fun_rere150->SetFitResult( result2, ipar150GeV );
  //  fun_rere200->SetFitResult( result2, ipar200GeV );

  float reso200_2 = fun_rere200->GetParameter(1);
  float reso150_2 = fun_rere150->GetParameter(1);
  float reso100_2 = fun_rere100->GetParameter(1);
  float reso50_2 = fun_rere50->GetParameter(1);
  float reso20_2 = fun_rere20->GetParameter(1);
 

  h2_axes_all2->Draw("");

  fun_rere20->SetFitResult( result2, ipar20GeV );
  fun_rere20->Draw("L same");
  fun_rere50->SetFitResult( result2, ipar50GeV );
  fun_rere50->Draw("L same");
  fun_rere100->SetFitResult( result2, ipar100GeV );
  fun_rere100->Draw("L same");
  fun_rere150->SetFitResult( result2, ipar150GeV );
  fun_rere150->Draw("L same");
  // fun_rere200->SetFitResult( result2, ipar200GeV );
  // fun_rere200->Draw("L same");

  //  gr_reso_vs_resp_200GeV->Draw("p same");
  gr_reso_vs_resp_20GeV ->Draw("p same");
  gr_reso_vs_resp_50GeV ->Draw("p same");
  gr_reso_vs_resp_100GeV->Draw("p same");
  gr_reso_vs_resp_150GeV->Draw("p same");
  
  TLegend* leg_up = new TLegend(0.45, 0.92 - 0.06 *8 , 0.85, 0.92);
  leg_up->SetTextSize(0.038);
  leg_up->AddEntry(gr_reso_vs_resp_20GeV,"Data 20GeV","p");
  leg_up->AddEntry( fun_rere20,Form("N =  %.00f\n #pm %.0f\n kADC",(fun_rere20->GetParameter(0))/1000., (fun_rere20->GetParError(0))/1000. ),"");
  leg_up->AddEntry( (TObject*)0,Form("C = %.2f\n #pm %.2f\n %s",(fun_rere20->GetParameter(1)), (fun_rere20->GetParError(1)),"%" ),"");

  leg_up->AddEntry(gr_reso_vs_resp_50GeV,"Data 50GeV","p");
  leg_up->AddEntry( fun_rere50,Form("N =  %.00f\n #pm %.0f\n kADC",(fun_rere50->GetParameter(0))/1000., (fun_rere50->GetParError(0))/1000. ),"");
  leg_up->AddEntry( (TObject*)0,Form("C = %.2f\n #pm %.2f\n %s",(fun_rere50->GetParameter(1)), (fun_rere50->GetParError(1)) ,"%"),"");

  leg_up->AddEntry(gr_reso_vs_resp_100GeV,"Data 100GeV","p");
  leg_up->AddEntry( fun_rere100,Form("N =  %.00f\n #pm %.0f\n kADC",(fun_rere100->GetParameter(0))/1000., (fun_rere100->GetParError(0))/1000. ),"");
  leg_up->AddEntry( (TObject*)0,Form("C = %.2f\n #pm %.2f\n %s",(fun_rere100->GetParameter(1)), (fun_rere100->GetParError(1)),"%" ),"");

  leg_up->AddEntry(gr_reso_vs_resp_150GeV,"Data 150GeV","p");
  leg_up->AddEntry( fun_rere150,Form("N =  %.00f\n #pm %.0f\n kADC",(fun_rere150->GetParameter(0))/1000., (fun_rere150->GetParError(0))/1000.),"");
  leg_up->AddEntry( (TObject*)0,Form("C  = %.2f\n #pm %.2f\n %s",(fun_rere150->GetParameter(1)), (fun_rere150->GetParError(1)) ,"%"),"");

  leg_up->SetFillColor(0);
  leg_up->Draw("same");




  label_top->Draw("same");

  c1->SaveAs( Form( "%s/reso_vs_resp_all_up.eps", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/reso_vs_resp_all_up.png", outputdir.c_str() ) );
 c1->SaveAs( Form( "%s/reso_vs_resp_all_up.pdf", outputdir.c_str() ) );













  //Systematic errors: Set Value 700k, then down to 400k
  // set the data range
  ROOT::Fit::DataRange range20GeV3; 
  range20GeV3.SetRange(  20, 500000 );
  ROOT::Fit::BinData data20GeV3(opt,range20GeV3); 
  data20GeV3.Opt().fCoordErrors = false; 
  ROOT::Fit::FillData(data20GeV3, gr_reso_vs_resp_20GeV);
  
  ROOT::Fit::DataRange range50GeV3; 
  range50GeV3.SetRange(  10, 500000 ); 
  ROOT::Fit::BinData data50GeV3(opt,range50GeV3);
  data50GeV3.Opt().fCoordErrors = false; 
  ROOT::Fit::FillData(data50GeV3, gr_reso_vs_resp_50GeV );
  
  ROOT::Fit::DataRange range100GeV3; 
  range100GeV3.SetRange(  20, 500000 );
  ROOT::Fit::BinData data100GeV3(opt,range100GeV3); 
  data100GeV3.Opt().fCoordErrors = false; 
  ROOT::Fit::FillData(data100GeV3, gr_reso_vs_resp_100GeV );
  
  ROOT::Fit::DataRange range150GeV3; 
  range150GeV3.SetRange(  20, 500000 );//one more point diff
  ROOT::Fit::BinData data150GeV3(opt,range150GeV3); 
  data150GeV3.Opt().fCoordErrors = false; 
  ROOT::Fit::FillData(data150GeV3, gr_reso_vs_resp_150GeV );
  
  ROOT::Fit::DataRange range200GeV3; 
  range200GeV3.SetRange(  400000, 1000000 );
  ROOT::Fit::BinData data200GeV3(opt,range200GeV3);
  data200GeV3.Opt().fCoordErrors = false; 
  ROOT::Fit::FillData(data200GeV3, gr_reso_vs_resp_200GeV );


  ROOT::Fit::Chi2Function chi2_20GeV3(data20GeV3, wf20GeV);
  ROOT::Fit::Chi2Function chi2_50GeV3(data50GeV3, wf50GeV);
  ROOT::Fit::Chi2Function chi2_100GeV3(data100GeV3, wf100GeV);
  ROOT::Fit::Chi2Function chi2_150GeV3(data150GeV3, wf150GeV);
  ROOT::Fit::Chi2Function chi2_200GeV3(data200GeV3, wf200GeV);
  

  ROOT::Fit::Fitter fitter3;
  fitter3.Config().SetParamsSettings(Npar,par0);
  fitter3.Config().MinimizerOptions().SetPrintLevel(3);
  fitter3.Config().ParSettings(0).SetLimits(450000, 650000);
  fitter3.Config().ParSettings(1).SetLimits( 2.5,3.5); 
  fitter3.Config().ParSettings(2).SetLimits(1,3);
  fitter3.Config().ParSettings(3).SetLimits(0.5,2);
  fitter3.Config().ParSettings(4).SetLimits(0.5,2);
  
  fitter3.Config().MinimizerOptions().SetPrintLevel(3);

  fitter3.Config().MinimizerOptions().SetMinimizerType("Minuit2");
  //  fitter3.Config().MinimizerOptions().SetMinimizerType("GSLMultiFit");

  GlobalChi2 globalChi23(Npar, (data20GeV3.Size()+data50GeV3.Size()+data100GeV3.Size()+data150GeV3.Size() ) ,chi2_20GeV3, chi2_50GeV3, chi2_100GeV3, chi2_150GeV3);

  fitter3.FitFCN(globalChi23, 0, data20GeV3.Size() + data50GeV3.Size() + data100GeV3.Size() + data150GeV3.Size() , true);
  ROOT::Fit::FitResult result3 = fitter3.Result();


  std::cout << std::endl;  std::cout << std::endl;
  std::cout << "------------- THIRD Result -------------- " << std::endl;
  result3.Print(std::cout);
  std::cout << std::endl;  std::cout << std::endl;


  fun_rere20->SetFitResult( result3, ipar20GeV );
  fun_rere50->SetFitResult( result3, ipar50GeV );
  fun_rere100->SetFitResult( result3, ipar100GeV );
  fun_rere150->SetFitResult( result3, ipar150GeV );
  // fun_rere200->SetFitResult( result3, ipar200GeV );

  float reso200_3 = fun_rere200->GetParameter(1);
  float reso150_3 = fun_rere150->GetParameter(1);
  float reso100_3 = fun_rere100->GetParameter(1);
  float reso50_3 = fun_rere50->GetParameter(1);
  float reso20_3 = fun_rere20->GetParameter(1);


  h2_axes_all2->Draw("");

  fun_rere20->SetFitResult( result3, ipar20GeV );
  fun_rere20->Draw("L same");
  fun_rere50->SetFitResult( result3, ipar50GeV );
  fun_rere50->Draw("L same");
  fun_rere100->SetFitResult( result3, ipar100GeV );
  fun_rere100->Draw("L same");
  fun_rere150->SetFitResult( result3, ipar150GeV );
  fun_rere150->Draw("L same");
  //  fun_rere200->SetFitResult( result3, ipar200GeV );
  // fun_rere200->Draw("L same");

  //  gr_reso_vs_resp_200GeV->Draw("p same");
  gr_reso_vs_resp_20GeV ->Draw("p same");
  gr_reso_vs_resp_50GeV ->Draw("p same");
  gr_reso_vs_resp_100GeV->Draw("p same");
  gr_reso_vs_resp_150GeV->Draw("p same");
  
  
  TLegend* leg_low = new TLegend(0.45, 0.92 - 0.06 *8 , 0.85, 0.92);
  leg_low->SetTextSize(0.038);
  leg_low->AddEntry(gr_reso_vs_resp_20GeV,"Data 20GeV","p");
 leg_low->AddEntry( fun_rere20,Form("N =  %.00f\n #pm %.0f\n kADC",(fun_rere20->GetParameter(0))/1000., (fun_rere20->GetParError(0))/1000. ),"");
 leg_low->AddEntry( (TObject*)0,Form("C = %.2f\n #pm %.2f\n %s",(fun_rere20->GetParameter(1)), (fun_rere20->GetParError(1)),"%" ),"");

  leg_low->AddEntry(gr_reso_vs_resp_50GeV,"Data 50GeV","p");
  leg_low->AddEntry( fun_rere50,Form("N =  %.00f\n #pm %.0f\n kADC",(fun_rere50->GetParameter(0))/1000., (fun_rere50->GetParError(0))/1000. ),"");
 leg_low->AddEntry( (TObject*)0,Form("C = %.2f\n #pm %.2f\n %s",(fun_rere50->GetParameter(1)), (fun_rere50->GetParError(1)) ,"%"),"");

  leg_low->AddEntry(gr_reso_vs_resp_100GeV,"Data 100GeV","p");
 leg_low->AddEntry( fun_rere100,Form("N =  %.00f\n #pm %.0f\n kADC",(fun_rere100->GetParameter(0))/1000., (fun_rere100->GetParError(0))/1000. ),"");
 leg_low->AddEntry( (TObject*)0,Form("C = %.2f\n #pm %.2f\n %s",(fun_rere100->GetParameter(1)), (fun_rere100->GetParError(1)),"%" ),"");

  leg_low->AddEntry(gr_reso_vs_resp_150GeV,"Data 150GeV","p");
  leg_low->AddEntry( fun_rere150,Form("N =  %.00f\n #pm %.0f\n kADC",(fun_rere150->GetParameter(0))/1000., (fun_rere150->GetParError(0))/1000.),"");
 leg_low->AddEntry( (TObject*)0,Form("C  = %.2f\n #pm %.2f\n %s",(fun_rere150->GetParameter(1)), (fun_rere150->GetParError(1)) ,"%"),"");

  leg_low->SetFillColor(0);
  leg_low->Draw("same");

  label_top->Draw("same");

  c1->SaveAs( Form( "%s/reso_vs_resp_all_low.eps", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/reso_vs_resp_all_low.png", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/reso_vs_resp_all_low.pdf", outputdir.c_str() ) );











  //stat only 
  gr_reso_vs_energy_stat->SetPoint(0, 20, reso20);
  gr_reso_vs_energy_stat->SetPoint(1, 50, reso50);
  gr_reso_vs_energy_stat->SetPoint(2, 100, reso100);
  gr_reso_vs_energy_stat->SetPoint(3, 150, reso150);
  //  gr_reso_vs_energy_stat->SetPoint(4, 200, reso200);
  
  
  gr_reso_vs_energy_stat->SetPointError(0, 0,sqrt(reso20_err*reso20_err -fermi*fermi));
  gr_reso_vs_energy_stat->SetPointError(1, 0,sqrt(reso50_err*reso50_err - fermi_50*fermi_50));
  //  gr_reso_vs_energy_stat->SetPointError(2, 0,sqrt(reso100_err*reso100_err -fermi*fermi));
  // gr_reso_vs_energy_stat->SetPointError(3, 0,sqrt(reso150_err*reso150_err - fermi*fermi));
  gr_reso_vs_energy_stat->SetPointError(2, 0, reso100_err);
  gr_reso_vs_energy_stat->SetPointError(3, 0, reso150_err);
 
  
  /*
  gr_reso_vs_energy_stat->SetPointError(0, 0, reso20_err);
  gr_reso_vs_energy_stat->SetPointError(1, 0, reso50_err);
  gr_reso_vs_energy_stat->SetPointError(2, 0, reso100_err);
  gr_reso_vs_energy_stat->SetPointError(3, 0, reso150_err);
  //  gr_reso_corr->SetPointError(4, 0, reso200_err);
  */

  //maximum difference to mean
  float max_20 = TMath::Max( TMath::Abs(reso20 - reso20_2),TMath::Abs(reso20 - reso20_3) );
  float max_50 = TMath::Max( TMath::Abs(reso50 - reso50_2),TMath::Abs(reso50 - reso50_3) );
  float max_100 =TMath::Max(TMath::Abs(reso100-reso100_2),TMath::Abs(reso100 - reso100_3));
  float max_150 =TMath::Max(TMath::Abs(reso150-reso150_2),TMath::Abs(reso150 - reso150_3));


  reso20_err = sqrt(reso20_err*reso20_err + max_20*max_20);
  reso50_err = sqrt(reso50_err*reso50_err +  max_50*max_50);
  reso100_err = sqrt(reso100_err*reso100_err +  max_100*max_100 );
  reso150_err = sqrt(reso150_err*reso150_err +  max_150*max_150 );
  // reso200_err = sqrt(reso200_err*reso200_err + max_20*max_20 );

  gr_reso_corr->SetPoint(0, 20, reso20);
  gr_reso_corr->SetPoint(1, 50, reso50);
  gr_reso_corr->SetPoint(2, 100, reso100);
  gr_reso_corr->SetPoint(3, 150, reso150);
  //  gr_reso_corr->SetPoint(4, 200, reso200);
   
  gr_reso_corr->SetPointError(0, 0, reso20_err);
  gr_reso_corr->SetPointError(1, 0, reso50_err);
  gr_reso_corr->SetPointError(2, 0, reso100_err);
  gr_reso_corr->SetPointError(3, 0, reso150_err);
  //  gr_reso_corr->SetPointError(4, 0, reso200_err);
  
 



  std::cout << max_20 << std::endl;
  std::cout << max_20/reso20*100. << std::endl;
  std::cout << max_50 << std::endl;
  std::cout << max_50/reso50*100. << std::endl;
  std::cout << max_100 << std::endl;
  std::cout << max_100/reso100*100. << std::endl;
  std::cout << max_150 << std::endl;
  std::cout << max_150/reso150*100. << std::endl;


















  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  gStyle->SetTitleYOffset(1.4); // => 1.15 if exponents
  
 ///////////////////////RESOLUTION///////////////////////////////
  TH2D* h2_axes2 = new TH2D( "axes", "", 100, -0.0, xMax/1000 , 10, 0.0, 5 );
  h2_axes2->SetXTitle("Beam Energy [GeV]");
  h2_axes2->SetYTitle("Energy Resolution [%]");
  h2_axes2->Draw("");
  
  // Data (1x1)
  gr_reso_vs_energy->SetMarkerStyle(20);
  gr_reso_vs_energy->SetMarkerSize(1.6);
  gr_reso_vs_energy->SetMarkerColor(kRed);
  // gr_reso_vs_energy->Draw("p same");
  
  TF1 *fun= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1])",0.9, xMax/1000.+5.);
  // TF1 *fun= new TF1("fun",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",0.9, xMax/1000.);
  fun->SetParameter(2,0.1);
  fun->SetParameter(1,1.1);
  gr_reso_vs_energy->Fit(fun,"RN");
  fun->SetLineWidth(1.);
  fun->SetLineColor(46);
  //fun->Draw("L same");
  
 
  // MC (1x1)
  gr_reso_vs_energy_simul->SetMarkerStyle(24);
  gr_reso_vs_energy_simul->SetMarkerSize(1.6);
  gr_reso_vs_energy_simul->SetMarkerColor(kBlue);

  // MC (5x5)
  gr_reso_vs_energy_simul_ideal->SetMarkerStyle(21);
  gr_reso_vs_energy_simul_ideal->SetMarkerSize(1.2);
  gr_reso_vs_energy_simul_ideal->SetMarkerColor(kBlack);
  // gr_reso_vs_energy_simul_ideal->Draw("p same");
 
  TF1 *fun1= new TF1("fun1","sqrt([0]*[0]/x+[1]*[1])",5.5, xMax/1000.+5.);
  // TF1 *fun1= new TF1("fun1",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))", 5.5, xMax/1000.+5.);
  fun1->SetMarkerSize(1.6); fun1->SetMarkerStyle(24);
  fun1->SetMarkerColor(kBlue);
  fun1->SetParameter(1,0.5); fun1->SetParameter(0, 10);
  gr_reso_vs_energy_simul->Fit(fun1,"RN");
  fun1->SetLineWidth(1.);
  fun1->SetLineColor(kBlue);
  // fun1->Draw("L same");
  
  //  TF1 *fun_ideal= new TF1("fun_ideal","sqrt([0]*[0]/x)",0.05, xMax/1000.+5.);
  TF1 *fun_ideal= new TF1("fun_ideal","sqrt([0]*[0]/x+[1]*[1])",0.05, xMax/1000.+5.);
  // TF1 *fun_ideal= new TF1("fun_ideal",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",0.005, xMax/1000.+5.);
  fun_ideal->SetMarkerSize(1.2); fun_ideal->SetMarkerStyle(21);
  fun_ideal->SetMarkerColor(kBlack);
  fun_ideal->SetParameter(1,0.5); fun_ideal->SetParameter(0, 10);
  gr_reso_vs_energy_simul_ideal->Fit(fun_ideal,"RN");
  fun_ideal->SetLineWidth(1.);
  fun_ideal->SetLineColor(kBlack);
  fun_ideal->SetLineStyle(2);
  fun_ideal->Draw("L same");
  
  
  gr_reso_corr->SetMarkerStyle(20);
  gr_reso_corr->SetMarkerSize(1.6);
  gr_reso_corr->SetMarkerColor(kBlue);
  // gr_reso_corr->Draw("p same");
  
  
  
  TF1 *fun_corr= new TF1("fun_corr",  "sqrt([0]*[0]/x+[1]*[1])",0.9, xMax/1000.+5);
  // TF1 *fun_corr= new TF1("fun_corr",  "sqrt([0]*[0]/x+[1]*[1]+ [2]*[2]/(x*x))",0.9, xMax/1000.+5);
  fun_corr->SetParameter(0,12.3);
  fun_corr->SetParameter(1,0.5);
  // fun_corr->SetParameter(2,0.02);
  gr_reso_corr->Fit(fun_corr,"RN");
  fun_corr->SetLineWidth(1.);
  fun_corr->SetLineColor(43);
  // fun_corr->Draw("L same");
  fun1->Draw("L same");
  
  gr_reso_vs_energy_stat->SetMarkerStyle(20);
  gr_reso_vs_energy_stat->SetMarkerSize(1.6);
  gr_reso_vs_energy_stat->SetMarkerColor(kBlue);
  //  gr_reso_vs_energy_stat->Draw("p same");
  

  // MC (1x1) NOBGO
  gr_reso_vs_energy_simul_noBGO->SetMarkerStyle(24);
  gr_reso_vs_energy_simul_noBGO->SetMarkerSize(1.6);
  gr_reso_vs_energy_simul_noBGO->SetMarkerColor(kBlack);
  gr_reso_vs_energy_simul_noBGO->Draw("p same");

  TF1 *fun_noBGO= new TF1("fun_noBGO","sqrt([0]*[0]/x+[1]*[1])",5.5, xMax/1000.+5.);
  fun_noBGO->SetMarkerSize(1.6); fun_noBGO->SetMarkerStyle(24);
  fun_noBGO->SetMarkerColor(kBlack);
  fun_noBGO->SetParameter(1,0.5); fun_noBGO->SetParameter(0, 10);
  gr_reso_vs_energy_simul_noBGO->Fit(fun_noBGO,"RN");
  fun_noBGO->SetLineWidth(1.);
  fun_noBGO->SetLineColor(kBlack);
  fun_noBGO->Draw("L same");


  gr_reso_vs_energy_simul->Draw("p same");
    

  TLegend* leg4 = new TLegend(0.32, 0.92-0.06*7 , 0.9, 0.92);
  leg4->SetTextSize(0.038);
  leg4->AddEntry(gr_reso_vs_energy,"Data (Uncorrected)","p");
  leg4->AddEntry(gr_reso_corr,"Data (Corrected)","p");
  
  leg4->AddEntry(fun_corr,Form("S = %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun_corr->GetParameter(0), (fun_corr->GetParError(0)),"%" ),"");
  leg4->AddEntry( (TObject*)0,Form("C =   %.2f\n #pm %.2f\n %s",(fun_corr->GetParameter(1)), (fun_corr->GetParError(1)),"%" ),"");
  // leg4->AddEntry( (TObject*)0,Form("N =  %.2f\n #pm %.2f\n GeV",(fun_corr->GetParameter(2)), (fun_corr->GetParError(2)) ),"");
  //  leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun_corr->GetChisquare(), fun_corr->GetNDF() ), "");

  leg4->AddEntry(gr_reso_vs_energy_simul ,"MC (1x1)","P");
  leg4->AddEntry((TObject*)0 ,Form("S = %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun1->GetParameter(0), (fun1->GetParError(0)),"%" ),"");
  leg4->AddEntry( (TObject*)0 ,Form("C =   %.2f\n #pm %.2f\n %s",(fun1->GetParameter(1)), (fun1->GetParError(1)),"%" ),"");
  // leg4->AddEntry( (TObject*)0 ,Form("N =   %.2f\n #pm %.2f\n GeV",(fun1->GetParameter(2)), (fun1->GetParError(2)) ),"");
  // leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun1->GetChisquare(), fun1->GetNDF() ), "");

  leg4->AddEntry(gr_reso_vs_energy_simul_noBGO ,"MC (1x1) no BGOs","P");
  leg4->AddEntry((TObject*)0 ,Form("S = %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun_noBGO->GetParameter(0), (fun_noBGO->GetParError(0)),"%" ),"");
  leg4->AddEntry( (TObject*)0 ,Form("C =   %.2f\n #pm %.2f\n %s",(fun_noBGO->GetParameter(1)), (fun_noBGO->GetParError(1)),"%" ),"");


 
  leg4->AddEntry( fun_ideal ,"MC (5x5)","L");
 /*
  leg4->AddEntry((TObject*)0 ,Form("S =   %.2f\n #pm %.2f\n %s / #sqrt{E [GeV]}",fun_ideal->GetParameter(0), (fun_ideal->GetParError(0)),"%" ),"");
  leg4->AddEntry( (TObject*)0 ,Form("C =   %.2f\n #pm %.2f\n %s",(fun_ideal->GetParameter(1)), (fun_ideal->GetParError(1)),"%" ),"");
  //leg4->AddEntry( (TObject*)0 ,Form("N =   %.2f\n #pm %.2f\n GeV",(fun_ideal->GetParameter(2)), (fun_ideal->GetParError(2)) ),"");
   leg4->AddEntry( (TObject*)0, Form("#chi^{2} / NDF = %.2f\n / %d",fun_ideal->GetChisquare(), fun_ideal->GetNDF() ), "");
  */
  leg4->SetFillColor(0);
  leg4->Draw("same");
 
  TPaveText* label_low = new TPaveText(0.165,0.175,0.5,0.21, "brNDC");
  label_low->SetFillColor(kWhite);
  label_low->SetTextSize(0.038);
  label_low->SetTextAlign(11); // align right
  label_low->SetTextFont(62);
  label_low->AddText( "W-CeF_{3} Single Tower");
  label_low->Draw("same");
 
  label_top2->Draw("same");

  

 
  c1->SaveAs( Form( "%s/resolution.eps", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/resolution.png", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/resolution.pdf", outputdir.c_str() ) );
















  ///////////////////////RESOLUTION cleaned up ///////////////////////////////
  TH2D* h2_axes_neat = new TH2D( "axes_neat", "", 100, -0.0, 165 , 4, 0.0, 4 );
  h2_axes_neat->SetXTitle("Beam Energy [GeV]");
  h2_axes_neat->SetYTitle("Energy Resolution [%]");
  h2_axes_neat->Draw("");
 
  // MC (1x1)
  gr_reso_vs_energy_simul->Draw("p same");
  // Data (1x1) CORRECTED
  fun_ideal->Draw("L same");

  gr_reso_corr->Draw("p same");
   
  //gr_reso_vs_energy_stat->Draw("p same");

  gr_reso_vs_energy_simul_ideal->Draw("p same");
 
  fun1->Draw("L same");

  TLegend* leg_neat = new TLegend(0.42, 0.92-0.06*5 , 0.9, 0.92);
  leg_neat->SetTextSize(0.038);
  leg_neat->AddEntry(gr_reso_corr,"Data Noise Subtracted","p");
  leg_neat->AddEntry( fun1,"MC (1x1)","LP");
  // leg_neat->AddEntry(gr_reso_vs_energy_simul ,"MC (1x1)","P");
  leg_neat->AddEntry( fun_ideal ,"MC (5x5)","LP");
  leg_neat->AddEntry((TObject*)0 ,Form("S =  %.2f\n%s / #sqrt{E [GeV]}",fun_ideal->GetParameter(0),"%" ),"");
  leg_neat->AddEntry( (TObject*)0 ,Form("C =  %.2f\n%s",(fun_ideal->GetParameter(1)) ,"%" ),"");

  leg_neat->SetFillColor(0);
  leg_neat->Draw("same");
 
  label_low->Draw("same");
 
  label_top2->Draw("same");
 
  c1->SaveAs( Form( "%s/resolution_neat.eps", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/resolution_neat.png", outputdir.c_str() ) );
  c1->SaveAs( Form( "%s/resolution_neat.pdf", outputdir.c_str() ) );

  c1->SaveAs( Form( "%s/resolution_neat_%d.eps", outputdir.c_str(), int(fermi*100) ) );
  c1->SaveAs( Form( "%s/resolution_neat_%d.png", outputdir.c_str(), int(fermi*100) ) );
  c1->SaveAs( Form( "%s/resolution_neat_%d.pdf", outputdir.c_str(), int(fermi*100) ) );



  double x, with, without;

  gr_resp_vs_energy_simul->GetPoint(1, x,with);


  gr_resp_vs_energy_simul_noBGO->GetPoint(1, x,without);



  std::cout << "With BGOs = " << with /  adcEnergyC950 *energyS  << std::endl;

  std::cout << "Without BGOs = " << without  /  adcEnergyC950*energyS   << std::endl;

  std::cout << "With / Without BGOs = " << with/without << std::endl;



  globalChiSquare = result.Chi2() / result.Ndf();


  if(globalChiSquare<1.0005)
    {doTheBreak=1;
      break;}
  // fermi += 0.05;


  delete gr_reso_vs_resp_20GeV ;
  delete gr_reso_vs_resp_50GeV ;
  delete gr_reso_vs_resp_100GeV ;
  delete gr_reso_vs_resp_150GeV ;
  delete gr_reso_vs_resp_200GeV ;
  delete gr_reso_vs_energy_stat ;
       
  delete gr_reso_vs_HV ;
  delete gr_resp_vs_HV ;
  
  delete gr_reso_vs_HV_50GeV ;
  delete gr_resp_vs_HV_50GeV ;
   
  delete gr_reso_vs_HV_100GeV ;
  delete gr_resp_vs_HV_100GeV ;
  
  delete gr_reso_vs_HV_150GeV ;
  delete gr_resp_vs_HV_150GeV ;
  
  delete gr_resp_vs_energy ;
  delete gr_reso_vs_energy ;
  delete gr_dev ;
  
  delete gr_resp_vs_energy_simul ;
  delete gr_reso_vs_energy_simul ;
  delete gr_dev_simul ;

  delete gr_resp_vs_energy_simul_ideal ;
  delete gr_reso_vs_energy_simul_ideal ;
  delete gr_dev_simul_ideal ;

  }  

 

    if(doTheBreak==1) break;
  }



 

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
  
  h1 = new TH1D(histoName.c_str(), "", 10000, xMin, xMax);
  tree->Project( histoName.c_str(), whatToProject.c_str(), cut.c_str() );
  f1 = new TF1( Form("gaus_%s", name.c_str()), "gaus", xMin, xMax);
  
  h1->Fit( f1, "RQN" );
  
  f1->SetLineColor(kRed);

  
  double peakpos = f1->GetParameter(1);
  double sigma = f1->GetParameter(2);
  double fitmin;
  double fitmax;

  if(sigma < 0.005 * peakpos) sigma = 0.005*peakpos;
  //  fitmax = peakpos + 3 * sigma;  //  fitmax = peakpos+2.7*sigma;
  //  fitmin = peakpos - 4 * sigma;  //  fitmin = peakpos - 1.3 *sigma;

  fitmax = int(peakpos + 4.8 * sigma ) ;  //  fitmax = peakpos+2.7*sigma;
  fitmin = int(peakpos - 6 * sigma ) ;  //  fitmin = peakpos - 1.3 *sigma;

  if(peakpos > 1000000) {
   fitmin = fitmin -15000; 
  fitmin = int(fitmin/10000) *10000;

  fitmax = fitmin+ int((fitmax-fitmin)/20000)*20000; 
  }

  /*
  if(peakpos < 500000 && peakpos >400000 ) {
     fitmax = fitmax - 3000; 
    fitmin = fitmin -6000; 
  fitmin = int(fitmin/1000) *1000;

  fitmax = fitmin+ int((fitmax-fitmin)/2000)*2000; 
  }
  */

  //  fitmax = fitmax -9000; fitmax = int(fitmax/10000) *10000;

  //  fitmax = peakpos + 4 * sigma;  //  fitmax = peakpos+2.7*sigma;
  //  fitmin = peakpos - 5 * sigma;  //  fitmin = peakpos - 1.3 *sigma;
  if(fitmin <5) fitmin =5;
  if(fitmax < 100) fitmax = 100;

  int nBins = 200;
  // int nBins = 60.;
   TH1D* h2 = new TH1D("h2", "", nBins, fitmin,fitmax);
   //   TH1D* h2 = new TH1D("h2", "", 200, fitmin,fitmax);
  //  TH1D* h2 = new TH1D("h2", "", int((fitmax-fitmin)/3000), fitmin,fitmax);
  tree->Project( "h2", whatToProject.c_str(), cut.c_str() );

  RooRealVar x("x","ADC Channel", fitmin, fitmax);
  RooDataHist data("data","dataset with x",x,Import(*h2) );
  
  RooPlot* frame;
  RooPlot* xframe = x.frame();   
  
  frame = x.frame("Title");
   data.plotOn(frame);  //this will show histogram data points on canvas
  // data.plotOn(frame, Binning( int((fitmax-fitmin)/500.) , fitmin,fitmax));
  //  data.statOn(frame);  //this will display hist stat on canvas
  
  RooRealVar meanr("meanr","Mean",peakpos,peakpos-sigma, peakpos+sigma);
  RooRealVar width("width","#sigma",sigma, 5.0, 15.*sigma);
  RooRealVar A("A","Dist",1., 0.0, 5.0);
  RooRealVar N("N","Deg",3, 0.0, 15);

  meanr.setRange( 55., 2500000.);
  width.setRange( 30, 100000);
  
  RooCBShape fit_fct("fit_fct","fit_fct",x,meanr,width,A,N); //int ndf = 4;

  fit_fct.fitTo(data);
  
  fit_fct.plotOn(frame,LineColor(4));//this will show fit overlay on canvas
  // fit_fct.paramOn(frame); //this will display the fit parameters on canvas
    
  double mean = meanr.getVal();
  double meanErr = meanr.getError();
  double rms = width.getVal();
  double rmsErr = width.getError();
  double reso = 100.* rms/mean; //in percent
  double resoErr = 100.* getRatioError( rms, mean, meanErr, rmsErr );
  double chiSquare = frame->chiSquare();
  
  TCanvas* cans = new TCanvas("cans", "un canvas", 600,600);
  cans->cd();
  /*
  TH2D* h2_axes = new TH2D( "axes", "", 5 , fitmin, fitmax, 5, 0, h2->GetBinContent(h2->GetMaximumBin())*1.1 );

  h2_axes->GetXaxis()->SetNdivisions(5, 0, 0, kFALSE);

  h2_axes->SetXTitle("ADC Channels");
  h2_axes->SetYTitle(Form("Events / (%.0f ADC Channels)", (fitmax-fitmin)/200.  )   );
  h2_axes->Draw();
  */ 
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.09);
  gPad->SetBottomMargin(0.15);
  gStyle->SetTitleYOffset(0.7); // => 1.15 if exponents
  
  frame->GetYaxis()->SetTitle(Form("Events / (%.0f ADC Channel)", (fitmax-fitmin)/nBins  ) );
  frame->GetXaxis()->SetNdivisions(5, 5, 0, kTRUE);
  //  frame->GetXaxis()->SetNdivisions(5, 5, 0, kFALSE);
  frame->Draw(""); 

  /*
  TLegend* lego = new TLegend(0.65, 0.7, 0.9, 0.92);  
  lego->SetTextSize(0.038);
  lego->AddEntry(  (TObject*)0 ,Form("#mu = %.0f #pm %.0f", meanr.getVal(), meanr.getError() ), "");
  lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.0f #pm %.0f ", width.getVal(), width.getError() ), "");
  lego->AddEntry(  (TObject*)0 ,Form("#chi^{2} = %.2f ", frame->chiSquare() ), "");
  lego->AddEntry(  (TObject*)0 ,Form("#sigma/#mu = %.2f #pm  %.2f ", reso , resoErr ), "");
  lego->AddEntry(  (TObject*)0 ,Form("A = %.1f #pm  %.1f ", A.getVal() , A.getError() ), "");
  lego->AddEntry(  (TObject*)0 ,Form("N = %.1f #pm  %.1f ", N.getVal() , N.getError() ), "");
  lego->SetFillColor(0);  lego->Draw("same");
  */


  TPaveText* label_fit = new TPaveText(0.19,0.89-0.06*3,0.4,0.89, "brNDC");
  label_fit->SetFillColor(kWhite);
  label_fit->SetTextSize(0.038);
  label_fit->SetTextAlign(10); // align right
  label_fit->SetTextFont(62);
  label_fit->AddText("W-CeF_{3} Single Tower");

  std::string mu_str = Form("#mu = (%.1f#pm%.1f)#upoint10^{3}", meanr.getVal()/1000., meanr.getError()/1000. );
  label_fit->AddText(mu_str.c_str());
  std::string sigma_str = Form("#sigma = (%.1f#pm%.1f)#upoint10^{3}", width.getVal()/1000., width.getError()/1000. );
  label_fit->AddText(sigma_str.c_str());
  /*
  std::string reso_str =Form("#sigma/#mu = %.2f #pm  %.2f %s ", reso , resoErr, "%" );
  label_fit->AddText(reso_str.c_str());
  std::string A_str =  Form("A = %.1f #pm  %.1f ", A.getVal() , A.getError() );
  label_fit->AddText(A_str.c_str());
  std::string N_str =  Form("N = %.1f #pm  %.1f ", N.getVal() , N.getError());
  label_fit->AddText(N_str.c_str());
  std::string chi_str = Form("#chi^{2} = %.2f ", chiSquare );
  label_fit->AddText(chi_str.c_str());
  */
  label_fit->Draw("same");

  TPaveText* label_top =  new TPaveText(0.3,0.953,0.945,0.975, "brNDC");
  label_top->SetFillColor(kWhite);
  label_top->SetTextSize(0.038);
  label_top->SetTextAlign(31); // align right
  label_top->SetTextFont(62);
  if(name == "329" || name == "471")
  label_top->AddText("100 GeV Electron Beam");
  label_top->Draw("same");


  cans->SaveAs( Form( "%s/CBFit_%s.png", outputdir.c_str(), name.c_str() ) );
  cans->SaveAs( Form( "%s/CBFit_%s.pdf", outputdir.c_str(), name.c_str() ) );
  

  FitStruct fs;
  fs.mean = mean;
  fs.mean_err = meanErr;
  fs.sigma = rms;
  fs.sigma_err = rmsErr;
  fs.reso = reso;
  fs.reso_err = resoErr;
  fs.chi2 = chiSquare;


  delete h1;
  delete f1;
  delete h2;
  delete cans;

  return fs;
}




FitStruct addPhotoStatistics( FitStruct rs ) {

  // MC response is already in MeV, 0.49 is the number of p.e. per MeV
  // RM = 90% of energy, so assume 90% of energy (0.9*491) is deposited in central channel
  // float  nADC = rs.mean/171.6*3236;
  float  nADC = rs.mean/169.445*3224.36; //without hodo cut
  float nPhotoElectrons = nADC/27.3;
  
  float poissonError = 100./sqrt( nPhotoElectrons ); // in percent
  
  float resoUnsmeared = rs.reso;
  float resoSmeared =  sqrt( rs.reso*rs.reso + poissonError*poissonError );
  
  rs.reso = resoSmeared;
  rs.reso_err = rs.reso_err * resoSmeared / resoUnsmeared; 

  return rs;
} 

























