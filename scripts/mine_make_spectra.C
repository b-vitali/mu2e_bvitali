//// bvitali, March 2020
// In order to update the spectra for ejected protons and deuterons,
// I made this script to generate the weights' tables.
// The parametrization is the one developed by Pasha for the general meeting
// in Feb 2020, after the release of the new AlCap-TWIST analysis
// 
// this is my version and is full of usless stuff. the right script should be in
// CoditionServices/data with the table this script can produce.
//
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
  if(m==938.3){    //protons
    par0 = 0.01;   // 0.009927 ± 0.0009604
    par1 = 0.505;  // 0.5048 ± 0.008764
    par2 = 1.8;    // 1.763 ± 0.6477
    par3 = 7.755;  // 7.755 ± 0.004043
    par4 = 3.4;    // 3.445 ± 0.3426
    par5 = 5.87;   // 5.869 ± 0.1362
  }

  else if(m==1875.6){ //deuterons
    par0 = 0.0011;    // 0.0011 ± 0.0003352
    par1 = 0.5;       // 0.5 ± 0
    par2 = 6.5;       // 6.543 ± 0.783
    par3 = 7.755;     // 7.755 ± 0
    par4 = 2.5;       // 2.492 ± 0.2597
    par5 = 7.6;       // 7.597 ± 0.2495
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
    y = (spectrum -> Eval(x))/integral; // /integral
    if(n%500==0){ cout<<x<<" "<<y<<endl; }
    output << x<<" "<<y<<endl;
    x += step;

  }
  output.close();
  cout<<"limit "<<limit<<endl;

}

//voit to make a histogram out of a weights' table
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


/////////////////////////////////////////////
//---- from here some additional stuff ----//
/////////////////////////////////////////////

//support void just to plot the functions and compare them
double kEmax=130;

TF1 * spectrum_p;
void see_proton(){
  spectrum_p   = new TF1("spectrum_p","EjectedProtonSpectrum_2020",0,kEmax,2);
  double par[2]={938.3,0};
  spectrum_p->SetParameters(par);
  //spectrum_p->SetNormalized(true);
  spectrum_p->SetNpx(100000);
  spectrum_p->Draw();
}

TF1 * spectrum_d;
void see_deuteron(){
  spectrum_d   = new TF1("spectrum_d","EjectedProtonSpectrum_2020",0,kEmax,2);
  double par[2]={1875.6,0};
  spectrum_d->SetParameters(par);
  //spectrum_d->SetNormalized(true);
  spectrum_d->SetLineColor(kBlue);
  spectrum_d->SetNpx(100000);
  spectrum_d->Draw();
}

TF1 * spectrum_h;
void see_h(){
  spectrum_h   = new TF1("spectrum_h","EjectedProtonSpectrum(x)",0,kEmax);
  //spectrum_h->SetNormalized(true);
  spectrum_h->SetLineColor(kGreen);  
  spectrum_h->SetNpx(100000);
  spectrum_h->Draw();
}

void comparison(){
  see_proton();
  see_deuteron();
  see_h();

  double integral_p = spectrum_p -> Integral(0,kEmax);
  double integral_d = spectrum_d -> Integral(0,kEmax);
  double integral_hp = spectrum_h -> Integral(0,kEmax);
  double integral_hd = spectrum_h -> Integral(0,kEmax);

  cout<<"TF1 Integral "<<integral_p<<" "<<integral_d<<" "<<integral_hp<<" "<<integral_hd<<endl;

  int N = 10000;
  double x[N];
  double y_p[N],y_d[N],y_hp[N],y_hd[N];
  for(int i = 0; i<N; i++){
    x[i] =i*kEmax/N;
    y_p[i] =  (spectrum_p -> Eval(x[i]));// /integral_p * 0.045;
    y_d[i] =  (spectrum_d -> Eval(x[i]));// /integral_d * 0.018;
    y_hp[i] =  (spectrum_h -> Eval(x[i]));// /integral_hp * 0.05;
    y_hd[i] =  (spectrum_h -> Eval(x[i]));// / integral_hd *0.05*0.5; //p:d=2:1 ?
  }

  double scale_p, scale_d;
  double sum_p, sum_d, sum_hp, sum_hd, sum_d_reco, sum_hd_reco;
  //rescale p>3.4 -> 0.032 ; d>4.5 -> 0.011
  double x_p = 3.4;
  double x_d = 4.5;
  for(int i=0;i<N;i++){
    if(x[i]>x_p) sum_p = sum_p + (x[i+1]-x[i])*(y_p[i+1]+y_p[i])/2;
    if(x[i]>x_p) sum_hp = sum_hp + (x[i+1]-x[i])*(y_hp[i+1]+y_hp[i])/2;

    if(x[i]>x_d) sum_d = sum_d + (x[i+1]-x[i])*(y_d[i+1]+y_d[i])/2;
    if(x[i]>x_d) sum_hd = sum_hd + (x[i+1]-x[i])*(y_hd[i+1]+y_hd[i])/2;

  } 
  cout<<"integral above 3.4p and 4.5d "<<sum_p<<" "<<sum_d<<" ( "<<sum_hp<<" "<<sum_hd<<")"<<endl;

  for(int j = 0; j<N; j++){
    y_p[j] = y_p[j] * 0.032 / sum_p;
    y_d[j] = y_d[j] * 0.012 / sum_d;
    y_hp[j] = y_hp[j] * 0.032 / sum_hp;
    y_hd[j] = y_hd[j] * 0.012 / sum_hd;
 }

  sum_p=0; sum_d=0;
  for(int i=0;i<N;i++){
    if(x[i]>x_p) sum_p = sum_p + (x[i+1]-x[i])*(y_p[i+1]+y_p[i])/2;
    if(x[i]>x_d) sum_d = sum_d + (x[i+1]-x[i])*(y_d[i+1]+y_d[i])/2;
  } 

  cout<<"rescaled "<<sum_p<<" "<<sum_d<<endl;
  sum_p=0; sum_d=0; sum_hp=0; sum_hd=0;
  for(int i=0;i<N;i++){
    sum_p = sum_p + (x[i+1]-x[i])*(y_p[i+1]+y_p[i])/2;
    sum_d = sum_d + (x[i+1]-x[i])*(y_d[i+1]+y_d[i])/2;
    sum_hp = sum_hp + (x[i+1]-x[i])*(y_hp[i+1]+y_hp[i])/2;
    sum_hd = sum_hd + (x[i+1]-x[i])*(y_hd[i+1]+y_hd[i])/2;
    
    if(x[i]>10 && x[i+1]<65){
      sum_d_reco = sum_d_reco + (x[i+1]-x[i])*(y_d[i+1]+y_d[i])/2;
      sum_hd_reco = sum_hd_reco + (x[i+1]-x[i])*(y_hd[i+1]+y_hd[i])/2;
    }
  }

  cout<<"final integrals "<<sum_p<<" "<<sum_d<<" "<<sum_hp<<" "<<sum_hd<<endl;
  cout<<"reco interval for d "<<sum_d_reco<<" "<<sum_hd_reco<<endl;

  TGraph * gr_p = new TGraph(N,x,y_p);
  gr_p->SetMarkerColor(kRed);
  gr_p->SetMarkerStyle(7);
  TGraph * gr_d = new TGraph(N,x,y_d);
  gr_d->SetMarkerColor(kBlue);
  gr_d->SetMarkerStyle(7);
  TGraph * gr_hp = new TGraph(N,x,y_hp);
  gr_hp->SetMarkerColor(kRed-7);
  gr_hp->SetMarkerStyle(7);
  TGraph * gr_hd = new TGraph(N,x,y_hd);
  gr_hd->SetMarkerColor(kBlue-7);
  gr_hd->SetMarkerStyle(7);

  TMultiGraph * mg = new TMultiGraph();
  mg->SetTitle("Comparison with Hungerford spectra");
  mg->Add(gr_p);
  mg->Add(gr_d);
  mg->Add(gr_hp);
  mg->Add(gr_hd);

  TCanvas * c_mg = new TCanvas();
  c_mg->cd();
  mg->Draw("AP");
  //spectrum_p->Draw("same");
  //spectrum_d->Draw("same");

  TLegend * l = new TLegend(0.6,0.7,0.9,0.9);
  l->AddEntry(gr_p,"gr_p","p");
  l->AddEntry(gr_d,"gr_d","p");
  l->AddEntry(gr_hp,"gr_hp","p");
  l->AddEntry(gr_hd,"gr_hd","p");

  double p[4] = {50,100,150,500};
  for(int i = 0; i<4; i++){
    double m = 938.3;
    double x = sqrt(p[i]*p[i]+m*m)-m;
    TLine * line = new TLine(x,0,x,0.008);
    line->SetLineStyle(i+1);
    line->SetLineColor(kRed);
    line->DrawLine(x,0,x,0.008);
    m = 1875.6;
    x = sqrt(p[i]*p[i]+m*m)-m;
    line->SetLineColor(kBlue);
    line->DrawLine(x,0,x,0.008);
    l->AddEntry(line,TString::Format("lines p=%5.1f MeV/c",p[i]),"l");
  }
  l->Draw();

}

void comparison_shape(){
  see_proton();
  see_deuteron();
  see_h();

  TLegend * l = new TLegend(0.6,0.7,0.9,0.9);
  l->AddEntry(spectrum_p,"new proton","l");
  l->AddEntry(spectrum_d,"new deuteron","l");
  l->AddEntry(spectrum_h,"h proton","l");


  spectrum_h->SetNormalized(true);
  spectrum_p->SetNormalized(true);
  spectrum_d->SetNormalized(true);

  spectrum_h->Draw();
  spectrum_p->Draw("same");
  spectrum_d->Draw("same");
  l->Draw();
  
  double p[4] = {50,100,150,500};
  for(int i = 0; i<4; i++){
    double m = 938.3;
    double x = sqrt(p[i]*p[i]+m*m)-m;
    TLine * line = new TLine(x,0,x,0.16);
    line->SetLineStyle(i+1);
    line->SetLineColor(kRed);
    line->DrawLine(x,0,x,0.16);
    m = 1875.6;
    x = sqrt(p[i]*p[i]+m*m)-m;
    line->SetLineColor(kBlue);
    line->DrawLine(x,0,x,0.16);
    l->AddEntry((TObject*)0,TString::Format("lines p=%5.1f MeV/c",p[i]),"");
  }
  l->Draw();

  double integ=0;
  integ = spectrum_p->Integral(0,1.4);
  l->AddEntry((TObject*)0,TString::Format("new_p_>Integ(0,1.4) = %.4f",integ),"");
  l->Draw();
}


//h parametrization. Not used but here for comparison
double EjectedProtonSpectrum(double e){
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
