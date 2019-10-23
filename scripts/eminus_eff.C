#include <iostream>
#include <TMinuit.h>
#import "TH1.h"
#import "TFile.h"
#import "TF1.h"
#import "TFitter.h"
//-----------------------------------------------------------------------------
Double_t eff_pdf(Double_t *x,Double_t *par){
  Double_t val=0;
  val=par[0]/(1+exp(1+par[1]*(x[0]-par[2])));
  return val;
 
}

//-----------------------------------------------------------------------------
double trk_eff(double P) {
  /*
    function to parametrize the efficiency, when adding statistic you may want to change
   the parameters twiking the Fit.
   1  p0           2.16312e-01   7.48075e-02   3.02994e-05  -5.98481e-03
   2  p1          -2.42336e-01   1.40564e-01   7.33587e-05   7.01581e-03
   3  p2           8.56444e+01   4.10006e+00   1.83116e-03  -1.79631e-05
   */

  double par[3] = {2.2e-01, -2e-2, 8.56e1 }; 
  double val=par[0]/(1+exp(1+par[1]*(P-par[2])));
  if (val < 0) val = 0;
  if (P < 65) val = 0;

  return val;
}

void efficiencies()
{
//-----------------------------------------------------------------------------
  // TH1F * r_0   = gh1("eminus_gun_stnmaker.xxx.000003.val2_stn.hist","Validation2","gen_0/p");
  // TH1F * r_1   = gh1("eminus_gun_stnmaker.xxx.000003.val2_stn.hist","Validation2","gen_1/p");
  TH1F * r_2   = gh1("/mu2e/app/users/bvitali/summer/eminus_gun_stnmaker.xxx.000003.val2_stn.hist","Validation2","gen_2/p");
  //TH1F * mc    = gh1("eminus_gun_stnmaker.xxx.000003.val2_stn.hist","Validation2","evt_0/mce");
  TH1F * mc = new TH1F("mc","particelle generate",1000,0,200);
  for(int i=300; i<550; i++){
    mc->SetBinContent(i,400);
    mc->SetBinError(i,0);
  }
  //mc->Scale(5); //was created for 20k event. now are 100k (x5)
 
  // TH1F * eff_0 = (TH1F *) r_0->Clone(); eff_0->Reset();
  // TH1F * eff_1 = (TH1F *) r_1->Clone(); eff_1->Reset();
  TH1F * eff_2 = (TH1F *) r_2->Clone(); eff_2->Reset();
  
  
  /*
  eff_0->Divide(r_0,mc);
  eff_0->SetLineColor(kRed);

  eff_1->Divide(r_1,mc);
  eff_1->SetLineColor(kGreen);
  */

  eff_2->Divide(r_2,mc);
  eff_2->SetLineColor(kBlue);
  eff_2->SetTitle("Reconstruction efficiency for electrons");
  eff_2->GetXaxis()->SetTitle("E_{e} [MeV]");
  
  int rbin = 5;
  eff_2->Rebin(rbin);
  eff_2->Scale(1./rbin);

//-----------------------------------------------------------------------------
  //fit to get the parametrization for the efficiency
  TF1 * f_eff = new TF1("f_eff",eff_pdf,60,110,3);
  f_eff->SetParameters( -0.3,0.3,100);
  eff_2->Fit("f_eff","L","",60,110);
  eff_2->GetYaxis()->SetTitle("Efficiency #varepsilon_{e}");
  /*
  THStack * s   = new THStack("s","Stack of efficiencies");
  s->Add(eff_0); 
  s->Add(eff_1);
  s->Add(eff_2);
  
  TCanvas * c_eff = new TCanvas("c_eff","",400,400);
  c_eff->cd();

  s->Draw("nostack");
  */
 
  

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
  cz->SetTitle("Spectrum of the DIOs");
  cz->GetXaxis()->SetTitle("E_{e} [MeV]");
  cz->GetYaxis()->SetTitle("#frac{1}{#Gamma_{0}} #frac{d#Gamma}{dE_{e}} [MeV^{-1}]");

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
  cz_eff->SetTitle("Spectrum of reconstructed  DIOs");
  THStack * DIO_s   = new THStack("DIO_s","Stack of efficiencies for DIO");

  for (int i=0; i<nbin; i++) {
    double p = cz_eff->GetBinCenter(i);
    double y = cz->GetBinContent(i)*f_eff->Eval(p);//trk_eff(p); //multiply cz by the eff(p)
    cz_eff->SetBinContent(i,y);
  }

//-----------------------------------------------------------------------------
  //whole Integral and from a to b
  Double_t cz_eff_integ = cz_eff->Integral(0,1061)*cz_eff->GetBinWidth(5);
  cout<<"whole Integral of cz_eff = "<<cz_eff_integ<<endl;
  
  Int_t a=80,b=110;
  TAxis * x_ax = cz_eff->GetXaxis();
  cz_eff_integ = cz_eff->Integral(x_ax->FindBin(a),x_ax->FindBin(b))*cz_eff->GetBinWidth(5);
  
//-----------------------------------------------------------------------------
  cout<<"Integral of cz_eff from "<<a<<" = "<<cz_eff_integ<<endl;
  cz_eff->SetLineColor(kRed);
  DIO_s->Add(cz);
  DIO_s->Add(cz_eff);

  TH1F * cz_eff_copy = (TH1F *) cz_eff -> Clone();
  cz_eff_copy->GetXaxis()->SetRange(x_ax->FindBin(a-23),x_ax->FindBin(b+23));
 
  TCanvas * c_DIO_eff = new TCanvas("c_DIO_eff","",400,400);
  c_DIO_eff->SetLogy();
  c_DIO_eff->cd();
  cz_eff_copy->Draw();
}
