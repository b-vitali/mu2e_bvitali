//// bvitali, March 2020
// In order to update the spectra for ejected protons and deuterons,
// I made this script to generate the weights' tables.
// The parametrization is the one developed by Pasha for the general meeting
// in Feb 2020, after the release of the new AlCap-TWIST analysis
////
#include "TF1.h"
#include "TH1.h"
#include <stdlib.h>
#include <string>
#include <iostream>
using namespace std;

//new parametrization by Pasha (General Meeting Feb 2020:  Doc-31745-v1)
double EjectedProtonSpectrum_2020(double *p_poit, double *par){
  double m = par[0];
  double choice = par[1];
  double e;
  if(choice == 0.){ e = p_poit[0];}
  else if(choice == 1.){
    double p = p_poit[0]; // p - Momentum (MeV/c)
    e = sqrt(p*p + m*m)-m; //momentum->kinetic energy rel (MeV/c)
  }
  else cout<<"nither energy or momentum?"<<endl;
  double par0,par1,par2,par3,par4,par5;

  //one function for proton and deuterons but different par values
  if(m==938.3){//protons
    par0 = 0.0099;
    par1 = 0.505;
    par2 = 1.8;
    par3 = 7.755;
    par4 = 3.4;
    par5 = 5.9;
  }

  else if(m==1875.6){//deuterons
    par0 = 0.0011;
    par1 = 0.5;
    par2 = 6.5;
    par3 = 7.755;
    par4 = 2.5;
    par5 = 7.6;
  }

  double f;
  if (e <= 0) {
    f = 0;
  }
  else if (e < par3) {
    f = pow(e/(1+par1*e),par2)*exp(-e/par4);
  }
  else {
    double c4 = pow(par3/(1+par1*par3),par2)*exp(-par3/par4);
    f = c4*exp(-(e-par3)/par5);
  }
  return f;
}

void make_spectra(){  
  double step = 0.1 ;
  TF1 * spectrum   = new TF1("spectrum","EjectedProtonSpectrum_2020",-0.1,1000.1,2);
  char choose;

  double par[2];
  
  //choose the particle
  cout<<"proton [p] or deuteron [d] spectrum?"<<endl;
  cin>>choose;

  //choose the variable
  cout<<"energy [0] or momentum [1] as variable?"<<endl;
  cin>>par[1];
  
  //Name the file accordingly to your choice
  TString out_name;

  if(choose == 'p') {
    par[0]=938.3;
    spectrum->SetParameters(par);
    out_name = "ejected_protons_";
  }
  else if(choose == 'd'){
    par[0]=1875.6;
    spectrum->SetParameters(par);
    out_name = "ejected_deuterons_";
  }  
  else {cout<<"Something is off"<<endl; return;}

  if(par[1] == 0)out_name = out_name + "energy_weights.tbl";
  else if(par[1] == 1)out_name = out_name + "momentum_weights.tbl";
  else cout<<"nither energy or momentum?"<<endl;
  
  //in order to keep the spectrum normalized to 1
  double integral = spectrum -> Integral(0,1000);
  cout<<integral<<endl;

  ofstream output;
  output.open(out_name);

  double x=0., y=0.;

  for(Int_t n = 0; n<1000/step; n+=1){
    y = (spectrum -> Eval(x))/integral;
    cout<<x<<" "<<y<<endl;
    output << x<<" "<<y<<endl;
    x += step;

  }
  output.close();  
  
}

//voit to make a histogram out of a weights' table
void read_txt(){
  //defined like this because we want the values to be the center of the bins
  TH1F * h = new TH1F("h","h",10000,-.05,999.95); 
  Int_t bin;

  string read;
  double x, y;
  Int_t line=0;

  cout<<"what table do you want to make a histogram?"<<endl;
  string read_file;
  cin>>read_file;
  ifstream input;
  input.open(read_file);
  if(input.is_open()){
    while(!input.eof()){
      std::getline(input,read);
      cout<<read<<endl;
      if(isdigit(read[0])||read[0]=='-'){
	stringstream(read) >> x >> y;

	bin=h->GetXaxis()->FindBin(x);
	h->SetBinContent(bin,y);
      }
      else if (read == "skip all"){cout<<"['skip all' at line: "<<line<<"] "<<endl; break;}
      else {std::cout<<"[Comment at line: "<<line<<"] "<<read<<endl;}
      line+=1;
    }
  }
  h->Draw("");
}


/////////////////////////////////////////////
//---- from here some additional stuff ----//
/////////////////////////////////////////////

//support void just to plot the functions
void see_proton(){
  TF1 * spectrum   = new TF1("spectrum","EjectedProtonSpectrum_2020",0,500,2);
  double par[2]={938.3,0};
  spectrum->SetParameters(par);
  spectrum->Draw();
}

void see_deuteron(){
  TF1 * spectrum   = new TF1("spectrum","EjectedProtonSpectrum_2020",0,1000,2);
  double par[2]={1875.6,0};
  spectrum->SetParameters(par);
  spectrum->Draw();
}


//old parametrization. Not used but here for comparison
double EjectedProtonSpectrum(double *p_poit, double *m_poit){
  //taken from GMC
  //
  //   Ed Hungerford  Houston University May 17 1999
  //   Rashid Djilkibaev New York University (modified) May 18 1999
  //
  //   e - Kinetic energy (MeV)
  //   p - Momentum (MeV/c)
  //
  //   Generates a proton spectrum similar to that observed in
  //   u capture in Si.  JEPT 33(1971)11 and PRL 20(1967)569
  double p = p_poit[0];
  double m = m_poit[0];

  double e = sqrt(p*p + m*m)-m; //momentum->energy rel

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
