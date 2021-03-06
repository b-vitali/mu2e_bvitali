//// bvitali, March 2020
// In order to update the spectra for ejected protons and deuterons,
// I made this script to generate the weights' tables.
// The parametrization is the one developed by Pasha for the general meeting
// in Feb 2020, after the release of the new AlCap-TWIST analysis
////

//new parametrization by Pasha (General Meeting Feb 2020:  Doc-31745-v1)
double EjectedProtonSpectrum_2020(double *p_poit, double *m_poit){

  double p = p_poit[0]; // p - Momentum (MeV/c)
  double m = m_poit[0];  
  double e = sqrt(p*p + m*m)-m; //momentum->kinetic energy rel (MeV/c)
  
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

//FillRandom an histogram with the proper parametrization,
//normalize it to 1 and save the weights in a .tbl
void make_proton_deuteron_weights(){
  TF1 * spectrum   = new TF1("spectrum","EjectedProtonSpectrum_2020",0,1000,1);
  char choose;
  int n;

  //Do you whant the proton or deuteron spectrum? how many?
  cout<<"proton [p] or deuteron [d] spectrum?"<<endl;
  cin>>choose;
  cout<<"how many should I generate?"<<endl;
  cin>>n;

  //Name the file accordingly to your choice
  TString out_name;
  if(choose == 'p') {
    spectrum->SetParameter(0,938.3);
    out_name = "ejected_protons_weights.tbl";
  }
  else if(choose == 'd'){
    spectrum->SetParameter(0,1875.6);
    out_name = "ejected_deuterons_weights.tbl";
  }  
  else {cout<<"Something is off"<<endl; return;}

  //create the histo, fill it and scale it
  TH1F * h = new TH1F("h","h",10000,-.05,999.05);
  h->FillRandom("spectrum",n);
  h->Scale(1./h->Integral());
  h->Draw("");
  spectrum->Draw("same");

  //save the weights on file
  //I wanted to save underflows but negative energy created problems (i from 1) 
  ofstream output;
  output.open(out_name);
  for(Int_t i=1; i< h->GetNbinsX()+2; i++){
    output << h->GetBinCenter(i)<<" "<<h->GetBinContent(i)<<endl;
  }
  output.close();
}

void test(){  
  Float_t step = 0.10000 ;
  TF1 * spectrum   = new TF1("spectrum","EjectedProtonSpectrum_2020",-0.1,1000.1,1);
  char choose;

  //Do you whant the proton or deuteron spectrum? how many?
  cout<<"proton [p] or deuteron [d] spectrum?"<<endl;
  cin>>choose;

  //Name the file accordingly to your choice
  TString out_name;
  if(choose == 'p') {
    spectrum->SetParameter(0,938.3);
    out_name = "ejected_protons_weights.tbl";
  }
  else if(choose == 'd'){
    spectrum->SetParameter(0,1875.6);
    out_name = "ejected_deuterons_weights.tbl";
  }  
  
  else {cout<<"Something is off"<<endl; return;}
  
  Float_t integral=spectrum->Integral(0,1000);
  cout<<integral<<endl;

  ofstream output;
  output.open(out_name);

  Float_t x=0., y=0.;

  for(Int_t n = 0; n<1000/step; n+=1){
    x += step;
    y = (spectrum -> Eval(x))/integral;
    cout<<x<<" "<<y<<endl;
    output << x<<" "<<y<<endl;
  }
  output.close();  
}

//voit to make a histogram out of a weights' table
//just 2 columns of numbers but can manage some commenting
void read_txt(){
  //defined like this because we want the values to be the center of the bins
  TH1F * h = new TH1F("h","h",10000,-.05,999.05); 
  Int_t bin;

  string read;
  Float_t x, y;
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
  TF1 * spectrum   = new TF1("spectrum","EjectedProtonSpectrum_2020",0,500,1);
  spectrum->SetParameter(0,938.3);
  spectrum->Draw();
}

void see_deuteron(){
  TF1 * spectrum   = new TF1("spectrum","EjectedProtonSpectrum_2020",0,1000,1);
  spectrum->SetParameter(0,1875.6);
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

//simple void to see how the spectra changed
void difference(){
  int n = 1e7;
  TF1 * spectrum   = new TF1("spectrum","EjectedProtonSpectrum",0,500,1);
  TF1 * spectrum_new   = new TF1("spectrum_new","EjectedProtonSpectrum_2020",0,500,1);
  TF1 * spectrum_deu   = new TF1("spectrum_deu","EjectedProtonSpectrum",0,500,1);
  TF1 * spectrum_deu_new   = new TF1("spectrum_deu_new","EjectedProtonSpectrum_2020",0,500,1);

  spectrum->SetParameter(0,938.3);
  spectrum_new->SetParameter(0,938.3);
  spectrum_deu->SetParameter(0,1875.6);
  spectrum_deu_new->SetParameter(0,1875.6);

  TH1F * h_spectrum = new TH1F("h_spectrum","",100,0,500);
  h_spectrum->FillRandom("spectrum",n);
  h_spectrum->Scale(1/h_spectrum->GetEntries());
  h_spectrum->SetLineColor(kRed);
  h_spectrum->SetLineStyle(2);
  TH1F * h_spectrum_deu = new TH1F("h_spectrum_deu","",100,0,500);
  h_spectrum_deu->FillRandom("spectrum_deu",n);
  h_spectrum_deu->Scale(1/h_spectrum_deu->GetEntries());
  h_spectrum_deu->SetLineStyle(2);
  TH1F * h_spectrum_new = new TH1F("h_spectrum_new","",100,0,500);
  h_spectrum_new->FillRandom("spectrum_new",n);
  h_spectrum_new->Scale(1/h_spectrum_new->GetEntries());
  h_spectrum_new->SetLineColor(kRed);
  TH1F * h_spectrum_deu_new = new TH1F("h_spectrum_deu_new","",100,0,500);
  h_spectrum_deu_new->FillRandom("spectrum_deu_new",n);
  h_spectrum_deu_new->Scale(1/h_spectrum_deu_new->GetEntries());

  h_spectrum->Draw("histo");
  h_spectrum_deu->Draw("histo same");
  h_spectrum_new->Draw("histo same");
  h_spectrum_deu_new->Draw("histo same");
}
