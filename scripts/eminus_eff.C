#include <iostream>
#include <TMinuit.h>
#import "TH1.h"
#import "TFile.h"
#import "TF1.h"
#import "TFitter.h"

//-----------------------------------------------------------------------------
double trk_eff(double P) {
  /*
    function to parametrize the efficiency, when adding statistic you may want to change
   the parameters twiking the Fit.
  */

   // 1  p0           1.26213e-01   3.65006e-03  -1.05833e-07   1.66691e-02
   // 2  p1          -3.87269e-03   8.10695e-05  -8.35686e-11   6.80404e-01
   // 3  p2           2.96238e-05   1.24173e-06   5.35688e-11   2.78023e+01

  double c[3] = {1.26213e-01, -3.87269e-03, 2.96238e-05 }; 
  double f = c[0]+c[1]*P+c[2]*P*P;
  if (f < 0) f = 0;
  if (P < 65) f = 0;

  return f;
}

void efficiencies()
{
//-----------------------------------------------------------------------------
  TH1F * r_0   = gh1("eminus_gun_stnmaker.val2_stn.hist","Validation2","gen_0/p");
  // TH1F * r_1   = gh1("eminus_gun_stnmaker.val2_stn.hist","Validation2","gen_1/p");
  TH1F * r_2   = gh1("eminus_gun_stnmaker.val2_stn.hist","Validation2","gen_2/p");
  TH1F * mc    = gh1("eminus_gun_stnmaker.val2_stn.hist","Validation2","evt_0/mce");
  
  TH1F * eff_0 = (TH1F *) r_0->Clone(); eff_0->Reset();
  // TH1F * eff_1 = (TH1F *) r_1->Clone(); eff_1->Reset();
  TH1F * eff_2 = (TH1F *) r_2->Clone(); eff_2->Reset();
  THStack * s   = new THStack("s","Stack of efficiencies");

  eff_0->Divide(r_0,mc);
  eff_0->SetLineColor(kBlue);

  // eff_1->Divide(r_1,mc);
  // eff_1->SetLineColor(kRed);

  eff_2->Divide(r_2,mc);
  eff_2->SetLineColor(kGreen);
 
//-----------------------------------------------------------------------------
   //fit to get the parametrization for the efficiency
  TF1 * f_eff = new TF1("f_eff","[0]+[1]*x+[2]*x**2",60,110);
  f_eff->SetParameters(8e-2,-3e-3,2e-5);
  eff_2->Fit("f_eff","L","",60,110);


  s->Add(eff_0);
  //s->Add(eff_1);
  s->Add(eff_2);

  TCanvas * c_eff = new TCanvas("c_eff","",400,400);
  c_eff->cd();

  s->Draw("nostack");

//-----------------------------------------------------------------------------
//create cz_histo 
  Int_t N=0;
  Int_t nbin=1061;
  std::ifstream cz_file;
 
  Double_t  * x;
  x = new Double_t [nbin];
  Double_t  * y;
  y = new Double_t [nbin];


  Int_t i=0;
  cz_file.open("ConditionsService/data/czarnecki_Al.tbl");
  
    while(true){
      cz_file>>x[i]>>y[i];
      if (cz_file.eof()) break;
      i=i+1;
    }
  
  cz_file.close();
  N=i;
  cout<<"loaded"<<endl;

  TH1F *cz=new TH1F("cz","cz_histogram",nbin,0,106);
  cz->SetLineColor(kBlue);
  for(Int_t j=0; j<nbin; j++){
    int ib = cz->GetXaxis()->FindBin(x[j]);
    cz->SetBinContent(ib,y[j]);
  }

  TCanvas * c_cz = new TCanvas("c_cz","",400,400);
  c_cz->cd();
  cz->Draw("");
  
  Double_t cz_integ = cz->Integral(0,1061)*cz->GetBinWidth(5);
  cout<<"whole Integral of cz = "<<cz_integ<<endl;



//-----------------------------------------------------------------------------
 //use the czarnecki file to evaluate the DIO reconstruction eff  
  TH1F * cz_eff = (TH1F *) cz->Clone(); //to have cz_eff to be like cz
  THStack * DIO_s   = new THStack("DIO_s","Stack of efficiencies for DIO");

  for (int i=0; i<nbin; i++) {
    double p = cz_eff->GetBinCenter(i);
    double y = cz->GetBinContent(i)*trk_eff(p); //multiply cz by the eff(p)
    cz_eff->SetBinContent(i,y);
  }

//-----------------------------------------------------------------------------
  //whole Integral and from a to b
  Double_t cz_eff_integ = cz_eff->Integral(0,1061)*cz_eff->GetBinWidth(5);
  cout<<"whole Integral of cz_eff = "<<cz_eff_integ<<endl;
  
  Int_t a=80,b=110;
  TAxis * x_ax = cz_eff->GetXaxis();
  cz_eff_integ = cz_eff->Integral(x_ax->FindBin(a),x_ax->FindBin(b))*cz_eff->GetBinWidth(5);

  cout<<"Integral of cz_eff from "<<a<<" = "<<cz_eff_integ<<endl;
  cz_eff->SetLineColor(kRed);
  DIO_s->Add(cz);
  DIO_s->Add(cz_eff);

  TH1F * cz_eff_copy = (TH1F *) cz_eff -> Clone();
  cz_eff_copy->GetXaxis()->SetRange(x_ax->FindBin(a),x_ax->FindBin(110));
 
  TCanvas * c_DIO_eff = new TCanvas("c_DIO_eff","",400,400);
  c_DIO_eff->Divide(2,1);
  c_DIO_eff->SetLogy();
  c_DIO_eff->cd(1);
  DIO_s->Draw("nostack");
  c_DIO_eff->cd(2);
  cz_eff_copy->Draw();
  // f_eff->Draw("SAME"); // per vedere l'efficienza come faccio? moltiplico tutto?

  
  
  

  
}
