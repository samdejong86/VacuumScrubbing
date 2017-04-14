#include <iostream>
#include <vector>
#include <sstream>

//some root classes
#include "TROOT.h" 
#include "TGraph.h"   //A Graph
#include "TDatime.h"  //Time stamp
#include "TMatrixD.h" //Matrix
#include "TString.h"  //string
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"  
#include "TDirectory.h"
#include "TPad.h"
#include "TMath.h"   //some math operations
#include "TSystem.h"
#include "TAxis.h"
#include "TColor.h"  //colours

//custom colours in RGB space                 R         G         B
TColor *CustomRed   = new TColor(1002, 210./256,  38./256,  48./256);
TColor *CustomBlue  = new TColor(1003,  60./256, 125./256, 196./256);

using namespace std;

#include "fitter.h"


TString dataDir = "../v2_All/";
TString IntCurDir = "../../superKEKBData/";

void ScrubOneDay(TDatime date, double &HER, double &LER, double &HERInt, double &LERInt){

  stringstream ss;

  ss<<date.GetMonth();
  TString month = ss.str();
  ss.str("");
  
  ss<<date.GetDay();
  TString day = ss.str();
  ss.str("");
  
  if(day.Length()==1) day = "0"+day;
  
  ss<<date.GetDate();
  TString d = ss.str();
  ss.str("");
  
  TString IntCurfile = IntCurDir+"IntegratedBeamCurrent_"+d+".root";
  TString filename = dataDir+"BEAST_2016-0"+month+"-"+day+".root";
  
  TFile *f_IntCur = new TFile(IntCurfile);
  TGraph *gr_IntCurHER = (TGraph*)f_IntCur->Get("grHERcurIntegral");
  TGraph *gr_IntCurLER = (TGraph*)f_IntCur->Get("grLERcurIntegral");


  TFile *f = new TFile(filename);
  if(f->IsZombie()) return;
  if(! f->GetListOfKeys()->Contains("tout")){
    cout<<"no tout\n";
    f->Close();
    return;
  }
  
  TTree *tout = (TTree*)f->Get("tout");

  UInt_t ts;
  vector<double> *he3 = new vector<double>;
  vector<double> *LERCurrent = new vector<double>;
  vector<double> *HERCurrent = new vector<double>;
  vector<double> *LERPressure = new vector<double>;
  vector<double> *HERPressure = new vector<double>;

  tout->SetBranchAddress("ts", &ts);
  tout->SetBranchAddress("HE3_rate", &he3);
  tout->SetBranchAddress("SKB_LER_current", &LERCurrent);
  tout->SetBranchAddress("SKB_HER_current", &HERCurrent);
  tout->SetBranchAddress("SKB_LER_pressure_average", &LERPressure);
  tout->SetBranchAddress("SKB_HER_pressure_average", &HERPressure);

  int nentries = tout->GetEntries();
    
  TMatrixD ydata(1,1);
  TMatrixD X(1,2);
  int ctr=0;
  
  //vector<UInt_t> t;
  vector<double> LC;
  vector<double> HC;
  
  for(int i=0; i<nentries; i++){
    if(i>nentries) break;
    tout->GetEntry(i);
    
    //data quality cuts
    if(LERCurrent->size()==0) continue;
    if(HERCurrent->size()==0) continue;
    if(HERPressure->size()==0) continue;
    if(LERPressure->size()==0) continue;
    if(he3->size()<4) continue;
    
    if(LERCurrent->at(0)<=10) continue;
    if(HERCurrent->at(0)<=10) continue;
    if(LERPressure->at(0)<=0) continue;
    if(HERPressure->at(0)<=0) continue;
    
    if(he3->at(0)==0||he3->at(1)==0||he3->at(2)==0||he3->at(3)==0) continue;
    
    double HERCURPRE = HERCurrent->at(0)*HERPressure->at(0);
    double LERCURPRE = LERCurrent->at(0)*LERPressure->at(0);
      
    //check for nan's
    if(HERCURPRE!=HERCURPRE) continue;
    if(LERCURPRE!=LERCURPRE) continue;
    
    
    double he3Ave = (he3->at(0)+he3->at(1)+he3->at(2)+he3->at(3))*0.25;
    
    ydata[ctr][0] = he3Ave;
    ydata.ResizeTo(ydata.GetNrows()+1, 1);
      
    X[ctr][0] = HERCURPRE;
    X[ctr][1] = LERCURPRE;
    
    LC.push_back(LERCurrent->at(0));
    HC.push_back(HERCurrent->at(0));
    
    X.ResizeTo(X.GetNrows()+1, X.GetNcols());
    ctr++;
    
  }
  
  
  ydata.ResizeTo(ydata.GetNrows()-1, 1);
  X.ResizeTo(X.GetNrows()-1, X.GetNcols());
  
  bool goodfit;
  
  TMatrixD soln = fitter(ydata, X, goodfit);
  if(!goodfit) return;  
  TMatrixD appliedSoln = MatrixMultiply(X, soln);
  TMatrixD Var = uncertainty(ydata, appliedSoln, X);
    
  double AverageHER = 0;
  double AverageLER = 0;
  double n=0;

  for(int i=0; i<(int)HC.size(); i=i+100){
    AverageHER+=soln[0][0]*X[i][0]/(HC[i]*HC[i]);
    AverageLER+=soln[1][0]*X[i][1]/(LC[i]*LC[i]);
    n++; 
  }

  HER = AverageHER/n;
  LER = AverageLER/n;
  HERInt = gr_IntCurHER->Eval(ts);
  LERInt = gr_IntCurLER->Eval(ts);

  f->Close();
  f_IntCur->Close();
}


int studyDays[16] = {20160405,
		     20160407,
		     20160412,
		     20160419,
		     20160422,
		     20160516,
		     20160517,
		     20160518,
		     20160523,
		     20160524,
		     20160525,		   
		     20160529,
		     20160608,
		     20160621,
		     20160622,
		     20160628};



void VacuumScrubbingFitter(int startdate, int enddate){

  TDatime start(startdate, 0);
  TDatime end(enddate,0);


  stringstream ss;

  TGraph *gr_HER = new TGraph();
  TGraph *gr_LER = new TGraph();

  int LCtr=0;
  int HCtr=0;

  double conversionFctr = 1/3.6e6;

  while(start.GetDate()<end.GetDate()){
   
    cout<<start.GetDate()<<endl;
    
    bool goodDay=true;
    for(int i=0; i<16; i++){
      if(start.GetDate()==studyDays[i]) goodDay=false;
    }

    if(!goodDay){
      cout<<"Machine study day, skipping\n";
      start = start.Convert()+86400;
      continue;
    }

    double HER=-1;
    double LER=-1;
    double IntHER, IntLER;

    ScrubOneDay(start, HER, LER, IntHER, IntLER);

    if(HER!=-1){
      gr_HER->SetPoint(HCtr, IntHER*conversionFctr, HER);
      HCtr++;
    }
    if(LER!=-1){
      gr_LER->SetPoint(LCtr, IntLER*conversionFctr, LER);
      LCtr++;
    }

    start = start.Convert()+86400;
  }



  TCanvas *c1 = new TCanvas();
  gr_HER->GetXaxis()->SetTitle("HER Integrated Current (A#timeshr)");
  gr_HER->GetYaxis()->SetTitle("Rate_{HER}/I_{HER}^{2}  (Hz/mA^{2})");
  
  gr_HER->SetMarkerColor(1002);
  gr_HER->SetLineColor(1002);
  gr_HER->SetMarkerStyle(20);
  gr_HER->Draw("AP");
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1,1);

  c1->Print("HER.pdf");
  c1->Print("HER.C");

  TCanvas *c2 = new TCanvas();
  gr_LER->GetXaxis()->SetTitle("LER Integrated Current (A#timeshr)");
  gr_LER->GetYaxis()->SetTitle("Rate_{LER}/I_{LER}^{2} (Hz/mA^{2})");

  gr_LER->SetMarkerColor(1003);
  gr_LER->SetLineColor(1003);
  gr_LER->SetMarkerStyle(20);
  gr_LER->Draw("AP");
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1,1);

  c2->Print("LER.pdf");
  c2->Print("LER.C");

}

int main(int argc, char* argv[]){

  VacuumScrubbingFitter(20160226, 20160629);

}
