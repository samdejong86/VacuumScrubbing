

void PowerFitter(TString filename){

  TFile *f = new TFile(filename);
  TGraph *gr = (TGraph*)f->Get("VaccumScrub");

  TGraph *logGr = new TGraph();

  double max=-999;

  int n=gr->GetN();
  for(int i=0; i<n; i++){
    double intCur;
    double rate;
    gr->GetPoint(i, intCur, rate);

    if(rate>max) max=rate;

    logGr->SetPoint(i, log(intCur), log(rate));


  }
  //logGr->SetMarkerStyle(20);
  TFitResultPtr result = logGr->Fit("pol1", "S");

  double A=exp(result->Value(0));
  double Aerror = exp(result->Value(0))*result->ParError(0);
  double k=-result->Value(1);
  double errork=result->ParError(1);

  cout<<A<<"\t"<<Aerror<<endl;
  cout<<k<<"\t"<<errork<<endl;


  TF1* Fn = new TF1("Fn", "[0]*x^(-[1])", 0, 700);
  Fn->SetParameter(0, A);
  Fn->SetParameter(1, k);


  Fn->SetParName(0, "a (Hz mA^{-1})");
  Fn->SetParName(1, "k");


  gr->SetMarkerStyle(20);

  if(filename.Contains("HER")){
    gr->SetMarkerColor(1002);
    Fn->SetLineColor(1002);
  } else {
    gr->SetMarkerColor(1003);
    Fn->SetLineColor(1003);
  }

  //TFitResultPtr finalResult = gr->Fit(Fn, "SQ");



  TString Ring = "LER";
  if(filename.Contains("HER")) Ring="HER";


  TFile *dPDIFile = new TFile(Ring+"DPDI.root");
  TGraph *grdpdi = (TGraph*)dPDIFile->Get(Ring+"DPDI");

  TGraph *grDPDILog = new TGraph();

  double min=999;

  int n=grdpdi->GetN();
  for(int i=0; i<n; i++){
    double intCur;
    double dpdi;
    grdpdi->GetPoint(i, intCur, dpdi);

    if(dpdi<min) min=dpdi;

    grDPDILog->SetPoint(i, log(intCur), log(dpdi));

  }

  //grDPDILog->Draw("AP");
  
  TFitResultPtr result = grDPDILog->Fit("pol1", "S");
  
  double A=exp(result->Value(0));
  double Aerror = exp(result->Value(0))*result->ParError(0);
  double k=-result->Value(1);
  double errork=result->ParError(1);
  
  TF1* Fn2 = new TF1("Fn2", "[0]*x^(-[1])", 0, 700);
  Fn2->SetParameter(0, A);
  Fn2->SetParameter(1, k);

  cout<<A<<"\t"<<Aerror<<endl;
  cout<<k<<"\t"<<errork<<endl;

  gr->SetMaximum(max*1.5);
  gr->SetMinimum(min/1.5);

  gr->Draw("AP");
  Fn->Draw("same");

  grdpdi->SetMarkerStyle(20);
  //grdpdi->Fit(Fn2,"S");
  grdpdi->Draw("sameP");
  Fn2->Draw("same");
  
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1,1);

}
