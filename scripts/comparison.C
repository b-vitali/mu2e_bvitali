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


double EjectedProtonSpectrum_new(double *p_poit, double *m_poit){
  //   e - proton kinetic energy (MeV)
  //   p - proton Momentum (MeV/c)
  //   Pasha at the general meeting Feb 2020:  Doc-31745-v1
  double p = p_poit[0];
  double m = m_poit[0];
  
  double e = sqrt(p*p + m*m)-m; //momentum->kinetic energy rel
  
  double par0,par1,par2,par3,par4,par5;

  if(m>800 && m<1100){
    par0 = 0.0099;
    par1 = 0.505;
    par2 = 1.8;
    par3 = 7.755;
    par4 = 3.4;
    par5 = 5.9;
  }

  else if(m>1800 && m<2200){
    par0 = 0.0011;
    par1 = 0.5;
    par2 = 6.5;
    par3 = 7.755;
    par4 = 2.5;
    par5 = 7.6;
  }

  double spectrumWeight;
  
  if (e >= par3)
    {
      spectrumWeight=exp(-e*par5);
    }
  
  else if (e < par3 && e>0)
    {
      double xarg=(e / (1 + par1*e));
      double xpot=std::pow(xarg,par2);
      double xexp=exp(-par4*e);
      spectrumWeight=par0*xpot*xexp;
    }
  else
    {
      spectrumWeight = 0.;
    }

  return spectrumWeight;
}

void make_txt(){
  TF1 * spectrum   = new TF1("spectrum","EjectedProtonSpectrum",0,500,1);
  spectrum->SetParameter(0,938.3);
  //spectrum->Draw();

  TH1F * h = new TH1F("h","h",5000,0,500);
  h->FillRandom("spectrum",300000);
  h->Scale(1./h->Integral());
  h->Draw("");
  
  h->Fill(1000);
  h->Fill(-1000);
  h->Fill(-100);

  ofstream output;
  output.open("make_txt.txt");
  for(Int_t i=0; i< h->GetNbinsX()+2; i++){
    output << h->GetBinCenter(i)<<" "<<h->GetBinContent(i)<<endl;
  }
  output.close();
}

void read_txt(){
  
  TH1F * h = new TH1F("h","h",5000,0,500); //
  Int_t bin;//

  string test;
  Float_t x, y;
  Int_t line=0;

  ifstream input;
  input.open("make_txt.txt");
  if(input.is_open()){
    while(!input.eof()){
      std::getline(input,test);

      if(isdigit(test[0])||test[0]=='-'){
	stringstream(test) >> x >> y;

	bin=h->GetXaxis()->FindBin(x);//
	h->SetBinContent(bin,y);//
      }
      else if (test == "skip all"){cout<<"['skip all' at line: "<<line<<"] "<<endl; break;}
      else {std::cout<<"[Comment at line: "<<line<<"] "<<test<<endl;}
      line+=1;
    }
  }
  h->Draw("");//
}


void difference(){
  TF1 * spectrum   = new TF1("spectrum","EjectedProtonSpectrum",0,500,1);
  TF1 * spectrum_new   = new TF1("spectrum_new","EjectedProtonSpectrum_new",0,500,1);
  TF1 * spectrum_deu   = new TF1("spectrum_deu","EjectedProtonSpectrum",0,500,1);
  TF1 * spectrum_deu_new   = new TF1("spectrum_deu_new","EjectedProtonSpectrum_new",0,500,1);

  spectrum->SetParameter(0,938.3);
  spectrum_new->SetParameter(0,938.3);
  spectrum_deu->SetParameter(0,1876.6);
  spectrum_deu_new->SetParameter(0,1876.6);

  TH1F * h_spectrum = new TH1F("h_spectrum","",100,0,500);
  h_spectrum->FillRandom("spectrum",100000);
  h_spectrum->Scale(1/h_spectrum->GetEntries());
  h_spectrum->SetLineColor(kRed);
  h_spectrum->SetLineStyle(2);
  TH1F * h_spectrum_deu = new TH1F("h_spectrum_deu","",100,0,500);
  h_spectrum_deu->FillRandom("spectrum_deu",100000);
  h_spectrum_deu->Scale(1/h_spectrum_deu->GetEntries());
  h_spectrum_deu->SetLineStyle(2);
  TH1F * h_spectrum_new = new TH1F("h_spectrum_new","",100,0,500);
  h_spectrum_new->FillRandom("spectrum_new",100000);
  h_spectrum_new->Scale(1/h_spectrum_new->GetEntries());
  h_spectrum_new->SetLineColor(kRed);
  TH1F * h_spectrum_deu_new = new TH1F("h_spectrum_deu_new","",100,0,500);
  h_spectrum_deu_new->FillRandom("spectrum_deu_new",100000);
  h_spectrum_deu_new->Scale(1/h_spectrum_deu_new->GetEntries());

  h_spectrum->Draw("histo");
  h_spectrum_deu->Draw("histo same");
  h_spectrum_new->Draw("histo same");
  h_spectrum_deu_new->Draw("histo same");



}
