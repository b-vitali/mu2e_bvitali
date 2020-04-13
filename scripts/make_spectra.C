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
//the funci
double EjectedProtonSpectrum_2020(double *x_poit, double *par){
  double m = par[0];
  double choice = par[1];
  double e;
  double jacob;

  if(choice == 0.){ e = x_poit[0]; jacob = 1;}
  else if(choice == 1.){
    double p = x_poit[0]; // p - Momentum (MeV/c)
    e = sqrt(p*p + m*m)-m; //momentum->kinetic energy rel (MeV/c)
    jacob = p / sqrt(p*p + m*m); //change of variable->jacobiano
  }
  else cout<<"choice = " << choice << "nither energy or momentum?"<<endl;
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
  f = f*jacob*par0;
  return f;
}

void make_spectra(){  
  double step = 0.1, limit ;
  TF1 * spectrum;
  char choose;
  double par[2];
  double m_p=938.3, m_d=1875.6;
  double p_max = 1000.; //kEnergy: 249.9 MeV deuterons,  432.9 MeV protons
  double e_max;

  //choose the particle
  cout<<"proton [p] or deuteron [d] spectrum?"<<endl;
  cin>>choose;

  //choose the variable
  cout<<"energy [0] or momentum [1] as variable?"<<endl;
  cin>>par[1];

  //different function to have the same limits
  if(par[1]==1) limit = p_max;
  else if(par[1]==0) {
    if(choose=='p') limit = sqrt(p_max*p_max+m_p*m_p)-m_p;
    else if(choose=='d') limit = sqrt(p_max*p_max+m_d*m_d)-m_d;
  }
  cout<<limit<<endl;
  spectrum = new TF1("spectrum","EjectedProtonSpectrum_2020",0.,limit,2);  

  //Name the file accordingly to your choice
  TString out_name;
  if(choose == 'p') {
    par[0]=m_p;
    spectrum->SetParameters(par);
    out_name = "ejected_protons_";
  }
  else if(choose == 'd'){
    par[0]=m_d;
    spectrum->SetParameters(par);
    out_name = "ejected_deuterons_";
  }  
  else {cout<<"Something is off"<<endl; return;}
  
  if(par[1] == 0)out_name = out_name + "energy_weights.tbl";
  else if(par[1] == 1)out_name = out_name + "momentum_weights.tbl";
  else cout<<"nither energy or momentum?"<<endl;
  
  cout<<"we are here"<<endl;

  //in order to keep the spectrum normalized to 1
  double integral = spectrum -> Integral(0,limit);
  cout<<integral<<endl;

  ofstream output;
  output.open(out_name);

  double x=0., y=0.;

  x=0;
  for(Int_t n = 0; n<limit/step; n+=1){
    y = (spectrum -> Eval(x))/integral;
    if(n%500==0){ cout<<x<<" "<<y<<endl; }
    output << x<<" "<<y<<endl;
    x += step;

  }
  output.close();
  cout<<"limit "<<limit<<endl;

}

//void to make a histogram out of a weights' table
void read_txt(){
  //defined like this because we want the values to be the center of the bins
  TH1F * h = new TH1F("h","h",10000,-.05,999.95); 
  Int_t bin;

  //TH1F * h_k = new TH1F("k","k",10000,-.05,999.95); 

  double m;
  char choose;
  cout<<"proton [p] or deuteron [d]?"<<endl;
  cin>>choose;
  if(choose == 'p') m = 938.3; 
  else if (choose == 'd') m = 1875.6;

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
	
	//x = sqrt(x*x + m*m)-m;
	//bin=h_k->GetXaxis()->FindBin(x);
	//h_k->AddBinContent(bin,y);
      }
      else if (read == "skip all"){cout<<"['skip all' at line: "<<line<<"] "<<endl; break;}
      else {std::cout<<"[Comment at line: "<<line<<"] "<<read<<endl;}
      line+=1;
    }
  }
  TCanvas * c = new TCanvas();
  c->cd();
  h->Draw("");
  
  //TCanvas * ck = new TCanvas();
  //ck->cd();
  //h_k->Draw();
}
