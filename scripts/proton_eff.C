#include <iostream>
#include <TMinuit.h>
#import "TH1.h"
#import "TFile.h"
#import "TF1.h"
#import "TFitter.h"

void percento( double x )
{
  /*
    at which energy should be cut if we want to rejected the fraction x of low energy StrawHits ?
  */

  TFile * f = new TFile("trackerMCCheck_read.hist");
  
  TH1F * ehit     = (TH1F *)f->Get("TrackerMCCheck/ehit");
  TH1F * ppr      = (TH1F *)f->Get("TrackerMCCheck/ppr");

  TCanvas * c = new TCanvas("c","Canvas THStack",800,500);
  c->cd();
  ehit->Draw();

  double whole_integ = 0;
  double partial_integ = 0;
  int bin_dx = 1;

  whole_integ = ehit->Integral();   
  std::cout<<"whole integral of ehit ="<<whole_integ<<std::endl;
  
  while(partial_integ < whole_integ*x/100.){
    bin_dx += 1;
    partial_integ = ehit->Integral(0,bin_dx);
  }
    
  TAxis *  xaxis = ehit->GetXaxis();
  double x_bin_dx = xaxis->GetBinLowEdge(bin_dx);
  std::cout<<"wholepartial integral of ehit from 0 to "<<x_bin_dx<<" (bin number "<<bin_dx<<") ="<<partial_integ<<std::endl;
  
  //x will not be exactly the requested fraction due to binning
  x = (double) partial_integ/whole_integ*100;
  std::cout<<"the treshold to cut the "<<x<<"%"<<" is "<< x_bin_dx<<std::endl;
}


double EjectedProtonSpectrum(double p){
  //taken from GMC
  //
  //   Ed Hungerford  Houston University May 17 1999
  //   Rashid Djilkibaev New York University (modified) May 18 1999
  //
  //   e - proton kinetic energy (MeV)
  //   p - proton Momentum (MeV/c)
  //
  //   Generates a proton spectrum similar to that observed in
  //   u capture in Si.  JEPT 33(1971)11 and PRL 20(1967)569

  double e = p*p/(2*938.3); //momentum->energy

  //these numbers are in MeV!!!!
  static const double emn = 1.4; // replacing par1 from GMC
  static const double par2 = 1.3279;
  static const double par3=17844.0;
  static const double par4=.32218;
  static const double par5=100.;
  static const double par6=10.014;
  static const double par7=1050.;
  static const double par8=5.103;
  
  double spectrumWeight;
  
  if (e >= 20)
    {
      spectrumWeight=par5*exp(-(e-20.)/par6);
    }
  
  else if(e >= 8.0 && e <= 20.0)
    {
      spectrumWeight=par7*exp(-(e-8.)/par8);
    }
  else if (e > emn)
    {
      double xw=(1.-emn/e);
      double xu=std::pow(xw,par2);
      double xv=par3*exp(-par4*e);
      spectrumWeight=xv*xu;
    }
  else
    {
      spectrumWeight = 0.;
    }
  return spectrumWeight;
}


void proton_eff(){
  /*
    Let's find the efficiency for proton reconstruction
  */
  TH1F * Ppr   = gh1("protonMDC2018.hist","Validation2","gen_0/ppr");
  Ppr->GetXaxis()->SetTitle("p [MeV/c]");
  Ppr->GetYaxis()->SetTitle("Number of Trk");
  TH1F * Precopr    = gh1("protonMDC2018.hist","Validation2","trk_0/precopr");
  Precopr->GetXaxis()->SetTitle("p [MeV/c]");
  TH1F * eff = new TH1F("eff","Efficiency of proton reconstruction",2000,0,400);
  eff->GetXaxis()->SetTitle("p [MeV/c]");
  eff->GetYaxis()->SetTitle("Efficiency #varepsilon_{p}");

  //create a draw the 'genp' spectrum trough the ejected proton spectrum (10M events)
  TH1F * genp = new TH1F("genp","Generated protons",2000,0,400);
  genp->GetXaxis()->SetTitle("p [MeV/c]");
  TF1 * spectrum = new TF1("spectrum","EjectedProtonSpectrum(x) ",0,400);
  TCanvas * c_spectrum = new TCanvas("c_spectrum","EjectedProtonSpectrum",400,400);
  c_spectrum->cd();
  spectrum->Draw();
  genp->FillRandom("spectrum",10000000);
  TCanvas * c_genp = new TCanvas("c_genp","Histo of generated particles",400,400);
  c_genp->cd();
  genp->Draw();

  eff->Divide(Precopr,genp);// eff->Divide(Precopr,Ppr);

  //draw true momenta of genp and trk
  int rbin = 25;  // <------------------------- change here the n of rebin
  Ppr->Rebin(rbin);
  Ppr->Scale(1./rbin);
  Precopr->Rebin(rbin);
  Precopr->Scale(1./rbin);
  TCanvas * c = new TCanvas("c","",400,400);
  c->Divide(2,1);
  c->cd(1);
  Ppr->Draw("hist");
  c->cd(2);
  Precopr->Draw("hist");

  //Draw the eff_histogram
  TCanvas * c_eff = new TCanvas("c_eff","",400,400);
  c_eff->cd();
  eff->DrawCopy("");

  // rebinning the eff_histogram
  rbin = 50;  // <------------------------- change here the n of rebin
  eff->Rebin(rbin);
  eff->Scale(1./rbin);
  TCanvas * c_eff2 = new TCanvas("c_eff2","",400,400);
  c_eff2->cd();
  eff->Draw("hist");
 
}
