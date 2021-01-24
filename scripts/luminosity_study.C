#include <iostream>
#include <TMinuit.h>
#import "TH1.h"
#import "TFile.h"
#import "TF1.h"
#import "TFitter.h"
#include <TStyle.h>
#include "TMatrixD.h"




double lum_ave = 3.9e7;

TH2F * sum_histos(){
  TString file0  = "/mu2e/data/personal/bvitali/proton1M/mdc2018/chunks/reco_2qsh_slim.000.hist";
  TH2F * lum0       = gh2(file0,"Validation2","evt_0/istlum_vs_ntrk"); 
  TH2F * lum = (TH2F *)lum0->Clone("lum");

  TString file_tmp;
  TH2F * lum_tmp;
  for(int i = 1; i< 7+1; i++){
    file_tmp  = "/mu2e/data/personal/bvitali/proton1M/mdc2018/chunks/reco_2qsh_slim.";
    if(i<10) {file_tmp += "00", file_tmp += std::to_string(i);}

    file_tmp += ".hist";
    std::cout<<"open: "<<file_tmp<<std::endl;
    lum_tmp  = gh2(file_tmp,"Validation2","evt_0/istlum_vs_ntrk");  
  
    lum->Add(lum_tmp);
  }
  return lum;
}

double Matrices4(TFitResultPtr r, double x){
  
    TMatrixD cor = r->GetCorrelationMatrix();
    TMatrixD cov = r->GetCovarianceMatrix();
    cov.Print();
    cor.Print();

    TMatrixD M(4,4);
    TMatrixD X_tmp(4,1);

    for(int i = 0; i<3; i++){
      X_tmp(i,0)=r->Parameter(i);
      for(int j = 0; j<3; j++){
	if(i==3 || j==3) M(3,3)=0;
	else M(i,j)=cov(i,j);
      }
    }
    X_tmp(3,0)=x;
    M(3,3)=sqrt(x)/500.;
    M.Print();
    X_tmp.Print();

    TMatrixD X(4,1);

    X(0,0)=1;
    X(1,0)=X_tmp(3,0);
    X(2,0)=X_tmp(3,0)*X_tmp(3,0);
    X(3,0)=X_tmp(1,0)+X_tmp(2,0);
   
    X.Print();

    TMatrixD MX(4,1);
    MX.Mult(M,X);
    MX.Print();

    TMatrixD XT(1,4);

    XT.Transpose(X);
    XT.Print();
    TMatrixD XTMX(1,1);
    XTMX.Mult(XT,MX);
    XTMX.Print();

    double uncert = sqrt(XTMX(0,0));
    return uncert;
}


double Matrices(TFitResultPtr r, double x){
  
    TMatrixD cor = r->GetCorrelationMatrix();
    TMatrixD cov = r->GetCovarianceMatrix();
    cov.Print();
    cor.Print();

    TMatrixD M(3,3);
    TMatrixD X_tmp(4,1);

    for(int i = 0; i<3; i++){
      X_tmp(i,0)=r->Parameter(i);
      for(int j = 0; j<3; j++){
        M(i,j)=cov(i,j);
      }
    }
    X_tmp(3,0)=x;
    M.Print();
    X_tmp.Print();

    TMatrixD X(3,1);

    X(0,0)=1;
    X(1,0)=X_tmp(3,0);
    X(2,0)=X_tmp(3,0)*X_tmp(3,0);
   
    X.Print();

    TMatrixD MX(3,1);
    MX.Mult(M,X);
    MX.Print();

    TMatrixD XT(1,3);

    XT.Transpose(X);
    XT.Print();
    TMatrixD XTMX(1,1);
    XTMX.Mult(XT,MX);
    XTMX.Print();

    double uncert = sqrt(XTMX(0,0));
    return uncert;
}

void contrario(TGraphErrors * gr, bool b){
  double fitLow = 0;
  double fitHigh = 20;

  double conf = 0.9;
  double conf_draw = 0.9;

  TCanvas * c_gr_contrario = new TCanvas("c_gr_contrario", "", 700, 500);
  TGraphErrors * gr_contrario = new TGraphErrors;;
  c_gr_contrario->SetGrid();

  TPad * pad = new TPad("pad","",0.02,0.02,0.97,0.7);
  TPad * pad_res = new TPad("pad_res","",0.02,0.7,0.97,0.97);

  
  pad->SetGrid();
  pad->Draw();
  pad->cd();
  double tmp_x, tmp_y;
  for (int i = 0; i<gr->GetN(); i++){
    gr->GetPoint(i,tmp_x,tmp_y);
    gr_contrario->SetPoint(i, tmp_y, tmp_x);
    gr_contrario->SetPointError(i, gr->GetErrorY(i), gr->GetErrorX(i));
  }
  gr_contrario->SetTitle();
  gr_contrario->InitPolynom();
  TF1 * f = new TF1("f", "pol2",0,20);
  TFitResultPtr r = gr_contrario->Fit("f", "EMS","",fitLow,fitHigh);

  if(b){
    double x = 4.5;
    double uncert;
    uncert = Matrices(r,x);
    std::cout<< f->Eval(x)<<" "<<uncert<<std::endl;
  }

  TGraphErrors * ge = new TGraphErrors(gr_contrario->GetN());
  ge->SetTitle();
  ge->SetName("ge");

  for (int i=0; i<gr_contrario->GetN(); i++) ge->SetPoint(i, gr_contrario->GetX()[i], 0);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(ge,conf_draw);
  ge->SetFillColor(4);
  ge->SetFillStyle(3002);
  ge->Draw("a3");
  ge->GetXaxis()->SetLimits(1,9);
  gr_contrario->Draw("psame");
  gr_contrario->GetXaxis()->SetLimits(1,9);

  c_gr_contrario->cd();
  pad_res->SetGrid();
  pad_res->Draw();
  pad_res->cd();

  float pool_i;
  TGraph * gr_pool = new TGraph;
  int j = 0;
  for (int i = 0; i<gr_contrario->GetN(); i++){
    gr_contrario->GetPoint(i,tmp_x,tmp_y);
    if(tmp_x<fitLow || tmp_x>fitHigh) continue;
    pool_i = ( tmp_y - f->Eval(tmp_x) )/gr_contrario->GetErrorY(i);
    gr_pool->SetPoint(j, tmp_x, pool_i);
    j++;
  }
  gr_pool->SetTitle();
  gr_pool->GetYaxis()->SetLabelSize(0.1);
  gr_pool->GetXaxis()->SetLabelSize(0);
  gr_pool->GetXaxis()->SetLimits(1,9);
  gr_pool->SetMarkerStyle(7);
  gr_pool->Draw("AP");


  int N = 30;
  double *x = new double [N];
  double *x_up = new double [N];
  double *x_down = new double [N];

  double err[N], err_up[N], err_down[N];  // error on the function at point x0
  double y[N];

  for(int k =0; k<N; k++){
    x[k]= 0.5+0.5*(double)k;
  }

  for(int k =0; k<N; k++){
    x_up[k]= x[k]+sqrt(x[k])/sqrt(1/0.0017);
    x_down[k]= x[k]-sqrt(x[k])/sqrt(1/0.0017);
  }
  r->GetConfidenceIntervals(N,1,1, x, err, conf, false);
  r->GetConfidenceIntervals(N,1,1, x_up, err_up, conf, false);
  r->GetConfidenceIntervals(N,1,1, x_down, err_down, conf, false);
  
  for(int k =0; k<N; k++){
    y[k]=f->Eval(x_down[k]);
    cout << " function value at " << x_down[k] << " = " << y[k] << " +/- " << err_down[k] << " err % " << err_down[k]/y[k]<<endl;

    y[k]=f->Eval(x[k]);
    cout << " function value at " << x[k] << " = " << y[k] << " +/- " << err[k] << " err % " << err[k]/y[k]<<endl;

    y[k]=f->Eval(x_up[k]);
    cout << " function value at " << x_up[k] << " = " << y[k] << " +/- " << err_up[k] << " err % " << err_up[k]/y[k]<<endl;


    cout << f->Eval(x[k])<< " +- " <<  (f->Eval(x_up[k])+err_up[k]) - (f->Eval(x_down[k])+err_down[k])<<endl;
  }
  
  TGraph * gr_perc = new TGraph;
  gr_perc->SetMarkerStyle(7);
  double tmp_perc;
  for(int k =0; k<N; k++){
    y[k]=f->Eval(x[k]);
    tmp_perc =  ( f->Eval(x_up[k]) - f->Eval(x_down[k]) )/2;
    cout << " Check: " << x[k] << " = " << y[k] << " +/- " << err[k] << " err % " << err[k]/y[k]<<"  from x "<<tmp_perc<<" % "<< tmp_perc/y[k]<<endl;

    //tmp_perc = (f->Eval(x_up[k]) - f->Eval(x_down[k])+ err_up[k]+err_down[k]) /(2*y[k]);
    //tmp_perc =  ( tmp_perc )/y[k] + err[k]/y[k];
    //tmp_perc =  sqrt( pow(( f->Eval(x_up[k]) - f->Eval(x_down[k]) ),2) + pow(err[k],2) )/y[k];
    tmp_perc = Matrices(r,x[k]); cout<<"matrix" <<y[k]<<" "<<tmp_perc<<endl; tmp_perc = tmp_perc / y[k];

    gr_perc->SetPoint(k,x[k],tmp_perc);
  }
  TCanvas * c_gr_perc = new TCanvas("c_gr_perc", "", 700, 500);
  c_gr_perc->SetGrid();
  gr_perc->Draw("AP");
  
}

void study(){
  TH2F* lum = sum_histos();
  new TCanvas;
  lum->Draw();

  new TCanvas;
  lum->ProfileY()->Draw();

  float percentage = 0.1;
  float PocheBriciole = 200;
  float y1 = 10e6;
  float y2 = y1*(1+percentage);
  int bin1 = lum->GetYaxis()->FindFixBin(y1);
  int bin_tmp = lum->GetYaxis()->FindFixBin(y2);
  TString name_tmp;
  std::cout<< y1 << " " << y2 <<std::endl;

  std::cout<<  lum->GetYaxis()->GetLast() <<std::endl;

  bool ok=false;
  
  TCanvas * c_compact = new TCanvas("c_compact", "", 700, 500);
 
  c_compact->Divide(4,3);
  int i = 1;

  TGraphErrors * gr = new TGraphErrors;
  gr->SetTitle("slices mean");
  TGraphErrors * gr_eff = new TGraphErrors;
  gr_eff->SetTitle("Efficiency");

  float prob = (1.6e-3)*(0.61)*(0.05+0.025);
  float mean, dmean;
  for(int bin2 = bin_tmp; bin2 <= lum->GetYaxis()->GetLast(); bin2 ++){
    y2 = lum->GetYaxis()->GetBinCenter(bin2);
    bin1 = lum->GetYaxis()->FindBin(y1);
    if( bin2< bin_tmp+4 ) std::cout<< y1 << " " << y2 <<std::endl;
    if( (y2-y1)/(y2+y1) > percentage && lum->Integral(1,20,bin1,bin2)> PocheBriciole ) ok = true;
    if(bin2 == lum->GetYaxis()->GetLast()) ok = true;
    if(ok){
      name_tmp = "projection_";
      name_tmp +=  std::to_string(y1/10e6);   
      name_tmp +=  "-";   
      name_tmp +=  std::to_string(y2/10e6); 
      //new TCanvas;
      TCanvas * c_tmp = new TCanvas("c_tmp", "", 700, 500);
      TH1F* h_tmp = (TH1F*)lum->ProjectionX(name_tmp,bin1,bin2,"");
      h_tmp->SetTitle(";Number of tracks;POT");

      /*
      TFitResultPtr r_tmp = h_tmp->Fit("gaus","S");
      mean = r_tmp->Parameter(1);
      dmean = r_tmp->ParError(1);
      */

      mean = h_tmp->GetMean(); 
      dmean = h_tmp->GetRMS() / h_tmp->GetEntries();

      c_compact->cd(i);
      h_tmp->DrawCopy();
      i++;

      std::cout<< mean << " " << dmean <<std::endl;

      gr->SetPoint(gr->GetN(), (y2+y1)/2, mean);
      gr->SetPointError(gr->GetN()-1, (y2-y1)/sqrt(12), dmean);

      //gr->SetPoint(gr->GetN(), (y2+y1)/2/lum_ave, mean);
      //gr->SetPointError(gr->GetN()-1, (y2-y1)/sqrt(12)/lum_ave, dmean);

      gr_eff->SetPoint(gr->GetN(), (y2+y1)/2, mean/((y2+y1)/2 * prob));
      gr_eff->SetPointError(gr->GetN(), (y2-y1)/sqrt(12),dmean/((y2+y1)/2 * prob)  );

      y1 = y2;
      bin1 = lum->GetYaxis()->FindBin(y1);
      ok = false;
    }
  }

  TCanvas * c_gr = new TCanvas("c_gr", "", 700, 500);
  c_gr->SetGrid();
  gr->Fit("pol2", "Q");
  gr->Draw("ap");


  TCanvas * c_gr_eff = new TCanvas("c_gr_eff", "", 700, 500);
  c_gr_eff->SetGrid();
  gr_eff->Draw("AP");
  
  contrario(gr, false);
  
}
