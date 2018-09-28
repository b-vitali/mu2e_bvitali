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

