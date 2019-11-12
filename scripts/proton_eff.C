#include <iostream>
#include <TMinuit.h>
#import "TH1.h"
#import "TFile.h"
#import "TF1.h"
#import "TFitter.h"

void generate_spectrum(TH1F* h, int N, double m){
  TF1 * spectrum   = new TF1("spectrum","EjectedProtonSpectrum",0,500,1);
  spectrum->SetParameter(0,m);
  //TCanvas * c_spectrum = new TCanvas("c_spectrum","EjectedProtonSpectrum",400,400);
  //c_spectrum->cd();
  //spectrum->Draw();
  h->FillRandom("spectrum",N);
  TCanvas * c_genp = new TCanvas("c_genp","Histo of generated particles",400,400);
  c_genp->cd();
  h->DrawCopy("histo");
};

void proton_eff(){//using the spectrum
  /*
    Let's find the efficiency for proton reconstruction

    The first lines are used to decide the type of file and the name
    
    Ppr is the momentum of the generated protons
    Precopr is the generated momentum of the reconstructed protons
    genp is the spectrum of protons generated via EjectedProtonSpectrum()
    
    proton_gun/pgun2/proton_gun.2.SHminE0_L


  */

  TString mdc_file  = "/mu2e/data/users/bvitali/MDC2018e/MDC2018e.SHminE002_L.hist";
  TString pgun_file = "/mu2e/data/users/bvitali/proton_gun/pgun2/proton_gun.2.SHminE0.hist";
  TString file;

  bool single = false ;                                         //<------------------------ decide if pgun or mdc
  int choice;
  printf("By default is it proton_gun? %d \n",single);
  printf("Do you want to change? (y=1 / n=0) \n");
  cin>>choice;
  if(choice==1) single = !single;
  if(single)  file = pgun_file;
  else        file = mdc_file;
  cout<<"Using the file:"<<file<<endl;

  //Protons
  TH1F * Ppr       = gh1(file,"Validation2","gen_3/ppr");      // <----------------- the input .hist 
  Ppr->GetXaxis()->SetTitle("p [MeV/c]");
  
  TH1F * Precopr   = gh1(file,"Validation2","trk_71/precopr");  // <---------------- the input .hist
  Precopr->GetXaxis()->SetTitle("p [MeV/c]");

  //create and draw the 'genp' spectrum of ejected proton spectrum (10M events)
  TH1F * genp      = new TH1F("genp","Generated protons",2000,0,500);
  genp->GetXaxis()->SetTitle("p [MeV/c]");
  generate_spectrum(genp,10000000,938.3);

  //create and draw the 'genp' spectrum of ejected Deuton spectrum (10M events)
  TH1F * genp_deu     = new TH1F("genp_deu","Generated deutons",2000,0,500);
  genp_deu->GetXaxis()->SetTitle("p [MeV/c]");
  generate_spectrum(genp_deu,10000000,1875.6);
  

  //create and draw the KineticEnergy spectrum of ejected proton spectrum (10M events) 
  if(false){      //<-----------------------Set to TRUE if you want the Kinetic Energy
    TH1F * genk = new TH1F("genk","Generated protons (ke)",2000,0,130);
    genk->GetXaxis()->SetTitle("E [MeV/c]");
    TF1 * spectrum2 = new TF1("spectrum2","EjectedProtonSpectrum2(x) ",0,130);
    TCanvas * c_spectrum2 = new TCanvas("c_spectrum2","EjectedProtonSpectrum2",400,400);
    c_spectrum2->cd();
    spectrum2->Draw();
    genk->FillRandom("spectrum2",10000000);
    TCanvas * c_genk = new TCanvas("c_genk","Histo of generated particles (ke)",400,400);
    c_genk->cd();
    genk->Scale(0.03/10000000);
    genk->Draw();
  }

  //rebin and draw momenta
  int rbin = 50;  // <------------------------- change here the n of rebin
  Ppr->Rebin(rbin);
  Precopr->Rebin(rbin);
  
  TCanvas * c2 = new TCanvas("c2","",400,400);
  c2->Divide(2,1);
  c2->cd(1);
  Ppr->Draw("hist");
  c2->cd(2);
  Precopr->Draw("hist");

  /*
    The Efficiency has 2 definitions
    - recostructed / saved       (easy but not usefull)
    - recostructed / generated   (the number you can multiply by something)
  */

  //to evaluate uncertenties in Devide(A,B)
  Precopr->Sumw2(); 
  Ppr->Sumw2();
  genp->Sumw2();
  genp_deu->Sumw2();


  //Efficiency definded as reconstructed/saved
  TH1F * eff = new TH1F("eff","Efficiency of proton reconstruction",2000/rbin,0,500);
  eff->GetXaxis()->SetTitle("p [MeV/c]");
  eff->GetYaxis()->SetTitle("Efficiency #varepsilon_{p}");
  TCanvas * c_eff = new TCanvas("c_eff","",400,400);
  c_eff->cd();
  eff->Divide(Precopr,Ppr,1,1,"b");
  eff->Draw("histo E");

  //Efficiency definded as reconstructed/generated
  TH1F * eff2 = (TH1F *) Precopr->Clone("eff2");
  eff2->SetTitle("Efficiency of proton reconstruction genp");
  eff2->GetXaxis()->SetTitle("p [MeV/c]");
  eff2->GetYaxis()->SetTitle("Efficiency #varepsilon_{p}");
 
  if(single){
    //efficiency as reconstructed/generated for signle proton event.
    int n_generated_protons = 100000;
    double min_p = 100;
    double max_p = 500;

    double_t rescale = (double_t) (n_generated_protons*( 500 / (max_p-min_p) ) / (2000/rbin)) ; 
    eff2->Scale(1/rescale);
    TCanvas * c_eff2 = new TCanvas("c_eff2","",400,400);
    c_eff2->cd();
    eff2->DrawCopy("histo E");

    //if done with an histo and the binomial error is the same
    /*
    TH1F * flat             = new TH1F("flat","Generated protons",2000,0,500);
    TF1 * flat_spectrum = new TF1("flat_spectrum","flat_spectrum",0,500,2);
    flat_spectrum->SetParameters(100,500);
    flat->GetXaxis()->SetTitle("p [MeV/c]");
    flat->FillRandom("flat_spectrum",100000*100); //100k in pgun2, non 10k!
    flat->Rebin(rbin);
    flat->Scale(1./100);
    flat->Sumw2();

    eff2->Divide(Precopr,flat,1,1,"b");
    TCanvas * c_eff2_flathist = new TCanvas("c_eff2_flathist","",400,400);
    c_eff2_flathist->cd();
    eff2->Draw("histo E");
    */
  }

  else {
    //using the luminosity as rescaling factor
    TH2F * istlum_vs_ntrk   = gh2(file,"Validation2","evt_0/istlum_vs_ntrk");  // <-------------change here "MDC2018e.SHminE002.hist    

    //change 0.0019 to 0.0018 and 3.8e7 to 3.9e7 (talking with Andy)
    double number_event     = (double) istlum_vs_ntrk->GetEntries();
    double mean_lum = 3.8e7; //istlum_vs_ntrk->GetMeanY();
    printf("mean %f\n",mean_lum);
    double prob_p_p = 0.0019*0.61*0.05; //stopped_mu_per_proton_on_target * nuclear_cap * proto_ejection (PASHA)
    double number_protons_genp =  mean_lum * number_event * prob_p_p; 
    double_t rescale = (double_t)(number_protons_genp/genp->GetEntries());
    printf("number of generated protons = %f",number_protons_genp);

    genp->Rebin(rbin);
    printf("rescale %f\n",rescale);
    genp->Scale(rescale);

    TCanvas * c_genp_scaled= new TCanvas("c_genp_scaled","can",400,400);
    c_genp_scaled->cd();
    genp->SetTitle("Generated protons (rescaled)");
    genp->Draw("histo E");


    eff2->Divide(Precopr,genp,1,1,"b");
    TCanvas * c_eff2 = new TCanvas("c_eff2","",400,400);
    c_eff2->cd();
    eff2->Draw("histo E");

    //Deuterons
    TH1F * Pdeu      = gh1(file,"Validation2","gen_0/pdeu");  
    Pdeu->GetXaxis()->SetTitle("p [MeV/c]");
    TH1F * DeuPreco  = gh1(file,"Validation2","trk_72/precopr");
    DeuPreco->GetXaxis()->SetTitle("p [MeV/c]");
  
    Pdeu->Rebin(rbin);
    DeuPreco->Rebin(rbin);

    TCanvas * c3 = new TCanvas("c3","",400,400);
    c3->Divide(2,1);
    c3->cd(1);
    Pdeu->Draw("hist");
    c3->cd(2);
    DeuPreco->Draw("hist");

    Pdeu->Sumw2();
    DeuPreco->Sumw2();

    //Efficiency definded as reconstructed/saved DEUTON
    TH1F * eff_deu = new TH1F("eff_deu","Efficiency of deuton reconstruction",2000/rbin,0,500);
    eff_deu->GetXaxis()->SetTitle("p [MeV/c]");
    eff_deu->GetYaxis()->SetTitle("Efficiency #varepsilon_{p}");
    TCanvas * c_eff_deu = new TCanvas("c_eff_deu","",400,400);
    c_eff_deu->cd();
    eff_deu->Divide(DeuPreco,Pdeu,1,1,"b");
    eff_deu->Draw("histo E");

    TH1F * eff2_deu = (TH1F *) DeuPreco->Clone("eff2_deu");
    eff2_deu->SetTitle("Efficiency of deuterons reconstruction genp");
    eff2_deu->GetXaxis()->SetTitle("p [MeV/c]");
    eff2_deu->GetYaxis()->SetTitle("Efficiency #varepsilon_{p}");

    genp_deu->SetTitle("Generated deuternos (rescaled)");
    genp_deu->Rebin(rbin);
    printf("rescale %f\n",rescale);
    genp_deu->Scale(rescale*0.5);//<--------deut:prot eject?-------BERTACCHI suggested 0.1. check (PASHA 0.5)
    TCanvas * c_genp_deu_scaled= new TCanvas("c_genp_deu_scaled","can2",400,400);
    c_genp_deu_scaled->cd();
    genp_deu->Draw("histo E");

    eff2_deu->Divide(DeuPreco,genp_deu,1,1,"b");
    TCanvas * c_eff2_deu = new TCanvas("c_eff2_deu","",400,400);
    c_eff2_deu->cd();
    eff2_deu->Draw("histo E");
    
    /*
    //fraction of the spectra
    TH1F * deu_p = (TH1F *) genp_deu->Clone("deu_p");
    deu_p->SetTitle("Deuteron spectrum / proton spectrum");
    deu_p->GetXaxis()->SetTitle("p [MeV/c]");
    deu_p->GetYaxis()->SetTitle("Efficiency #varepsilon_{p}");
    deu_p->Divide(genp);
    deu_p->Draw("histo E");
    */
  }
}

double flat_spectrum(double* p, double* par){
  double w = 0;
  if(p[0]>=par[0] && p[0]<=par[1]) w=1;
  return w;
}

void GenpCheck(){
  // SINGLE PROTON EVENT
  TH1F * genp_ppr         = gh1("/mu2e/data/users/bvitali/proton_gun/proton_gun.2.SHminE0_L.hist","Validation2","gen_0/ppr"); //gen3 geometry, gen0 all
  TH1F * flat             = new TH1F("flat","Generated protons",2000,0,500);
  TF1 * flat_spectrum = new TF1("flat_spectrum","flat_spectrum",0,500,2);
  flat_spectrum->SetParameters(100,500);
  flat_spectrum->Draw();
  flat->GetXaxis()->SetTitle("p [MeV/c]");
  flat->FillRandom("flat_spectrum",100000*100); //100k in pgun2, non 10k!
  flat->Scale(1./100);
 
  genp_ppr->Sumw2();
  flat->Sumw2();

  int rbin = 50;
  genp_ppr->Rebin(rbin);
  flat->Rebin(rbin);

  TCanvas * c_1= new TCanvas("c_1","c_1",400,400);
  c_1->cd();
  flat->DrawCopy("hist");

  TCanvas * c_2= new TCanvas("c_2","c_2",400,400);
  c_2->cd();
  genp_ppr->DrawCopy("hist");


  TH1F * eff = (TH1F *) flat->Clone("eff");
  eff->SetTitle("Efficiency of saved protons");
  eff->GetXaxis()->SetTitle("p [MeV/c]");
  eff->GetYaxis()->SetTitle("Efficiency #varepsilon_{p}");
  eff->Divide(genp_ppr,flat);
  
  TCanvas * c_3= new TCanvas("c_3","c_3",400,400);
  c_3->cd();
  eff->Draw("hist E");
  
  //MDC2018e
  TH1F * genp_ppr_mdc      = gh1("/mu2e/data/users/bvitali/MDC2018e/MDC2018e.ff5f.SHminE002_L.hist","Validation2","gen_0/ppr");
  
  TH1F * eject = new TH1F("eject","Generated ejected protons",2000,0,500);
  eject->GetXaxis()->SetTitle("p [MeV/c]");
  TF1 * spectrum   = new TF1("spectrum","EjectedProtonSpectrum",0,500,1);
  spectrum->SetParameter(0,938.3);
  eject->FillRandom("spectrum",10000000);

  double mean_lum = 3.8e7;
  printf("mean %f\n",mean_lum);
  double prob_p_p = 0.0019*0.61*0.05; //stopped_mu_per_proton_on_target * nuclear_cap * proto_ejection
  double number_protons_genp =  mean_lum * 500 * prob_p_p; 
  double_t rescale = (double_t)(number_protons_genp/eject->GetEntries());
  printf("number of generated protons = %f",number_protons_genp);
  printf("rescale %f\n",rescale);

  eject->Scale(rescale);
  
  genp_ppr_mdc->Sumw2();
  eject->Sumw2();
  genp_ppr_mdc->Rebin(rbin);
  eject->Rebin(rbin);

  TCanvas * c_4= new TCanvas("c_4","c_4",400,400);
  c_4->cd();
  eject->DrawCopy("histo"); 

  TCanvas * c_5= new TCanvas("c_5","c_5",400,400);
  c_5->cd();
  genp_ppr_mdc->DrawCopy("hist e");

  TH1F * eff_mdc = (TH1F *) eject->Clone("eff_mdc");
  eff_mdc->SetTitle("Efficiency of saved protons MDC");
  eff_mdc->GetXaxis()->SetTitle("p [MeV/c]");
  eff_mdc->GetYaxis()->SetTitle("Efficiency #varepsilon_{p}");
  eff_mdc->Divide(genp_ppr_mdc,eject,1,1,"b");
  
  TCanvas * c_6= new TCanvas("c_6","c_6",400,400);
  c_6->cd();
  eff_mdc->Draw("hist E");

  /*
  //DEUTERON
  TH1F * genp_pdeu_mdc      = gh1("/mu2e/data/users/bvitali/MDC2018e/MDC2018e.ff5f.SHminE002_L.hist","Validation2","gen_0/pdeu");//$$$$$$$
  
  TH1F * eject_deu = new TH1F("eject_deu","Generated ejected deuterons",2000,0,500);
  eject_deu->GetXaxis()->SetTitle("p [MeV/c]");
  TF1 * spectrum_deu   = new TF1("spectrum_deu","EjectedProtonSpectrum",0,500,1);
  spectrum_deu->SetParameter(0,938.3*2);
  eject_deu->FillRandom("spectrum_deu",10000000);

  double prob_p_d = 0.0019*0.61*0.025; //stopped_mu_per_proton_on_target * nuclear_cap * deuteron_ejection
  double number_deu_genp =  mean_lum * 500 * prob_p_d; 
  double_t rescale_deu = (double_t)(number_deu_genp/eject->GetEntries());
  printf("number of generated deuterons = %f",number_deu_genp);
  printf("rescale %f\n",rescale_deu);

  eject_deu->Scale(rescale_deu);
  
  genp_pdeu_mdc->Sumw2();
  eject_deu->Sumw2();
  genp_pdeu_mdc->Rebin(rbin);
  eject_deu->Rebin(rbin);

  TCanvas * c_4_deu= new TCanvas("c_4_deu","c_4_deu",400,400);
  c_4_deu->cd();
  eject_deu->DrawCopy("histo"); 

  TCanvas * c_5_deu= new TCanvas("c_5_deu","c_5_deu",400,400);
  c_5_deu->cd();
  genp_pdeu_mdc->DrawCopy("hist e");

  TH1F * eff_mdc_deu = (TH1F *) eject_deu->Clone("eff_mdc_deu");
  eff_mdc_deu->SetTitle("Efficiency of saved deuteron MDC");
  eff_mdc_deu->GetXaxis()->SetTitle("p [MeV/c]");
  eff_mdc_deu->GetYaxis()->SetTitle("Efficiency #varepsilon_{p}");
  eff_mdc_deu->Divide(genp_pdeu_mdc,eject_deu,1,1,"b");
  
  TCanvas * c_6_deu= new TCanvas("c_6_deu","c_6_deu",400,400);
  c_6_deu->cd();
  eff_mdc_deu->Draw("hist E");
  */
}


/*
  Momentum spectrum for both protons and deutons (different mass as parameter)
*/
double EjectedProtonSpectrum(double *p_poit, double *m_poit){
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
  double p = p_poit[0];
  double m = m_poit[0];

  // double e = p*p/(2*m); //momentum->energy
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

/*
  Energy spectrum for protons (copied from EjectedProtonSpectrum)
*/
double EjectedProtonSpectrum2(double e){
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


/*
  At which energy should be cut if we want to rejected the fraction x of low energy StrawHits ?
*/
void percento( double x )
{
  

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
