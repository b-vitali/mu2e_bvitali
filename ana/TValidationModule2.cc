//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
// 
//  3  : events with set C tracks and 70mm < |dx|  < 90 mm
//  4  : events with DpF > 1 MeV : obviously, misreconstructed ones
//  5  : events with N(tracks) > 1
//  6  : events trk_41 with 0.8< E/P < 1.1 - tracks missed by CalPatRec
//  7  : events (muo) with LogLHRCal >   20
//  8  : events (ele) with LogLHRCal < - 20
//  9  : events (muo) with 0.42 < E/P < 0.46
// 10  : events (muo) with Set C track with ECL > 80 MeV
// 28  : Set C DEM tracks with E/P > 1.1
// 29  : TRK_19 (Set C DEM tracks with a cluster) and LLHR(cal) < 0
// 31  : EVT_6 events with ce_costh > 0.8 
// 32  : TRK_1 events with chi2tcm > 100. 
// 33  : DU < -80mm - study edge effects
// 34  : EVT_7: events with E_CL > 60 and no tracks (makes sense only for single CE events)
// 35  : TRK_1: events with P > 106 MeV/c - misreconstruction
// 36  : TRK_23 events with P < 80: odd misidentified muons - turned out to be DIO electrons
// 37  : TRK_26 LLHR_CAL > 5
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
// #include "CalorimeterGeom/inc/HexMap.hh"

#include "ana/TValidationModule2.hh"

ClassImp(TValidationModule2)
//-----------------------------------------------------------------------------
TValidationModule2::TValidationModule2(const char* name, const char* title):
  TStnModule(name,title)
{
  fPtMin  = 1.;
  fTrackNumber.Set(100);

  fFillDioHist     = 1;

  fMinT0 = 0; // do not cut on time by default

  fTrackID = new TStnTrackID();
  fLogLH   = new TEmuLogLH();
//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
//-----------------------------------------------------------------------------
  fPdgCode       = 2212;
  fGeneratorCode = 28;			// conversionGun, 28:StoppedParticleReactionGun
  fTrackBlockName = "TrackBlock";
}

//-----------------------------------------------------------------------------
TValidationModule2::~TValidationModule2() {
}


//-----------------------------------------------------------------------------
void TValidationModule2::BookTimeClusterHistograms   (TimeClusterHist_t*   Hist, const char* Folder){
  
    HBook1F(Hist->fNHits         ,"nhits"      ,Form("%s: number of StrawHits; nSH"              ,Folder), 150,   0,   150,Folder);
    HBook1F(Hist->fNComboHits    ,"ncombohits" ,Form("%s: number of ComboHits; nCH"              ,Folder), 150,   0,   150,Folder);
    HBook1F(Hist->fT0            ,"t0"         ,Form("%s: t_{0}[ns]"       ,Folder), 800, 400,  1700,Folder);
    HBook1F(Hist->fClusterEnergy ,"clusterE"   ,Form("%s: cluster energy; E [MeV]      ",Folder), 400,   0,  200,Folder);  

    HBook2F(Hist->fNComboHitsVsP ,"ncombohits_vs_p" ,Form("%s: number of ComboHits vs generated p; p [MeV/c]; nCH",Folder), 2400, 0,600, 200, 0, 200, Folder);              //
    HBook2F(Hist->fFracSHVsP ,"fracSH_vs_p" ,Form("%s: Fraction of hits in the TimeCLuster vs generated p; p [MeV/c]; sh(TC)/sh(TOT)",Folder), 2400, 0,600, 110,0,1.1, Folder);              //

}
//-----------------------------------------------------------------------------
void TValidationModule2::BookHelixHistograms   (HelixHist_t*   Hist, const char* Folder){
  
    HBook1F(Hist->fNHits         ,"nhits"      ,Form("%s: number of StrawHits; nSH"            ,Folder), 150,   0,   150,Folder);
    HBook1F(Hist->fClusterTime   ,"clusterTime",Form("%s: cluster time; t_{cluster}[ns]"   ,Folder), 800, 400,  1700,Folder);
    HBook1F(Hist->fClusterEnergy ,"clusterE"   ,Form("%s: cluster energy; E [MeV]      "   ,Folder), 400,   0,  200,Folder);
    HBook1F(Hist->fRadius        ,"radius"     ,Form("%s: curvature radius; r [mm]"        ,Folder), 600,   0,   600,Folder);
    HBook1F(Hist->fMom           ,"p"          ,Form("%s: momentum; p [MeV/c]"             ,Folder), 2400,  0,   600,Folder);
    HBook1F(Hist->fPt            ,"pT"         ,Form("%s: transverse momentum; pT [MeV/c]" ,Folder), 600,   0,   150,Folder);
    HBook1F(Hist->fLambda        ,"lambda"     ,Form("%s: lambda; -#lambda [mm/rad]"       ,Folder), 2500,  0,   2500,Folder);
    HBook1F(Hist->fT0            ,"t0"         ,Form("%s: t0; t0[ns]"                      ,Folder), 100,   0,    10,Folder);
    HBook1F(Hist->fT0Err         ,"t0err"      ,Form("%s: t0err; t0err [ns]"               ,Folder), 100,   0,    10,Folder);
    HBook1F(Hist->fD0            ,"d0"         ,Form("%s: D0; d0 [mm]"                     ,Folder), 1600,   -400,    400,Folder);
    HBook2F(Hist->fChi2XYNDof    ,"chi2ndofxy" ,Form("%s: Chi2/NDof XY vs generate p; p [MeV/c]; #chi2/dof"    ,Folder), 2400, 0,600, 150, 0,   50,Folder);
    HBook2F(Hist->fChi2PhiZNDof    ,"chi2ndofphiz" ,Form("%s: Chi2/NDof #Phi Z vs generate p; p [MeV/c]; #chi2/dof" ,Folder), 2400, 0,600, 150, 0,   50,Folder);

    HBook2F(Hist->fLambdaVsP ,"lambda_vs_p" ,Form("%s:Lambda vs generated p; p [MeV/c]; -#lambda [mm/rad]",Folder), 2400, 0,600, 2500, 0, 2500, Folder);            
    HBook2F(Hist->fNRotVsP ,"nrot_vs_p" ,Form("%s: rotations vs generated p; p [MeV/c]; Rot.",Folder), 2400, 0,600, 1000, 0, 5, Folder);    
    HBook2F(Hist->fRadiusVsP,"radius_vs_p" ,Form("%s: radius vs generated p; p [MeV/c]; Rot.",Folder), 2400, 0,600, 600, 100, 700, Folder);    

    HBook2F(Hist->fFracSHVsP ,"fracsh_vs_p" ,Form("%s: fraction of SH vs generated p; p [MeV/c]; nSH active/nSH",Folder), 2400, 0,600, 200, 0, 1.1, Folder);         //
    HBook2F(Hist->fPrecoVsP ,"preco_vs_p" ,Form("%s: reconstructed vs generated p; p [MeV/c]; p [MeV/c]",Folder), 2400, 0,600, 2400, 0,600, Folder);                           //



}

//-----------------------------------------------------------------------------
void TValidationModule2::BookTrackSeedHistograms   (TrackSeedHist_t*   Hist, const char* Folder){
  
    HBook1F(Hist->fNHits         ,"nhits"    ,Form("%s: # of straw hits"              ,Folder), 150,   0,   150,Folder);
    HBook1F(Hist->fClusterTime ,"clusterTime",Form("%s: cluster time; t_{cluster}[ns]",Folder), 800, 400,  1700,Folder);
    HBook1F(Hist->fClusterEnergy ,"clusterE" ,Form("%s: cluster energy; E [MeV]      ",Folder), 400,   0,  200,Folder);
    HBook1F(Hist->fRadius        ,"radius"   ,Form("%s: curvature radius; r [mm]"     ,Folder), 600,   0,   600,Folder);
    HBook1F(Hist->fMom           ,"p"        ,Form("%s: momentum; p [MeV/c]"          ,Folder), 1000,   0,   1000,Folder);
    HBook1F(Hist->fPt            ,"pT"       ,Form("%s: pT; pT [MeV/c]"               ,Folder), 600,   0,   600,Folder);
    HBook1F(Hist->fTanDip        ,"tanDip"   ,Form("%s: tanDip; tanDip"               ,Folder), 300,   0,     3,Folder);
    HBook1F(Hist->fChi2          ,"chi2"   ,Form("%s: #chi^{2}-XY; #chi^{2}/ndof"     ,Folder), 100,   0,    10,Folder);
    HBook1F(Hist->fFitCons      ,"FitCons" ,Form("%s: Fit consistency; Fit-cons"      ,Folder), 100,   0,    1, Folder);
    HBook1F(Hist->fD0            ,"d0"       ,Form("%s: D0; d0 [mm]"                  ,Folder), 1600,   -400,    400,Folder);

  
}

//-----------------------------------------------------------------------------
void TValidationModule2::BookClusterHistograms(ClusterHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fVaneID ,"vane_id",Form("%s: Vane ID"       ,Folder), 10, 0,  10,Folder);
  HBook1F(Hist->fEnergy ,"energy" ,Form("%s: Cluster Energy",Folder),150, 0, 300,Folder);
  HBook1F(Hist->fT0     ,"t0"     ,Form("%s: cluster T0"    ,Folder),200, 0,2000,Folder);
  HBook1F(Hist->fRow    ,"row"    ,Form("%s: cluster Row"   ,Folder),200, 0, 200,Folder);
  HBook1F(Hist->fCol    ,"col"    ,Form("%s: cluster column",Folder),200, 0, 200,Folder);
  HBook1F(Hist->fX      ,"x"      ,Form("%s: cluster X"     ,Folder),200, -5000,5000,Folder);
  HBook1F(Hist->fY      ,"y"      ,Form("%s: cluster Y"     ,Folder),200,-1000,1000,Folder);
  HBook1F(Hist->fZ      ,"z"      ,Form("%s: cluster Z"     ,Folder),200, 11500,13500,Folder);
  HBook1F(Hist->fR      ,"r"      ,Form("%s: cluster Radius",Folder),100, 0,  1000,Folder);
  HBook1F(Hist->fYMean  ,"ymean"  ,Form("%s: cluster YMean" ,Folder),400,-400,400,Folder);
  HBook1F(Hist->fZMean  ,"zmean"  ,Form("%s: cluster ZMean" ,Folder),400,-400,400,Folder);
  HBook1F(Hist->fSigY   ,"sigy"   ,Form("%s: cluster SigY"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fSigZ   ,"sigz"   ,Form("%s: cluster SigZ"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fSigR   ,"sigr"   ,Form("%s: cluster SigR"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fNCr0   ,"ncr0"   ,Form("%s: cluster NCR[0]",Folder),100, 0,100,Folder);
  HBook1F(Hist->fNCr1   ,"ncr1"   ,Form("%s: cluster NCR[1]",Folder),100, 0,100,Folder);
  HBook1F(Hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fSigE1  ,"sige1"   ,Form("%s: SigmaE/Etot"  ,Folder),200, 0, 10,Folder);
  HBook1F(Hist->fSigE2  ,"sige2"   ,Form("%s: SigmaE/Emean" ,Folder),200, 0, 10,Folder);
}

//-----------------------------------------------------------------------------
void TValidationModule2::BookGenpHistograms(GenpHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP      ,"p"       ,Form("%s: momentum; p [MeV/c] "     ,Folder),1000,     0, 1000,Folder);
  HBook1F(Hist->fPdgCode[0],"pdg_code_0",Form("%s: PDG Code[0]"     ,Folder),4000, -100, 3000,Folder);
  HBook1F(Hist->fPdgCode[1],"pdg_code_1",Form("%s: PDG Code[1]"     ,Folder),500, -2500, 2500,Folder);
  HBook1F(Hist->fGenID  ,"gen_id"  ,Form("%s: Generator ID" ,Folder), 100,     0, 100,Folder);
  HBook1F(Hist->fZ0     ,"z0"      ,Form("%s: Z0"           ,Folder), 500,  5400, 6400,Folder);
  HBook1F(Hist->fT0     ,"t0"      ,Form("%s: T0"           ,Folder), 200,     0, 2000,Folder);
  HBook1F(Hist->fR0     ,"r"       ,Form("%s: R0"           ,Folder), 100,     0,  100,Folder);
  HBook1F(Hist->fCosTh  ,"cos_th"  ,Form("%s: Cos(Theta)"   ,Folder), 200,   -1.,   1.,Folder);

  HBook1F(Hist->fPpr      ,"ppr"        ,Form("%s:generated p for protons; p [MeV/c]",Folder), 2400,  0, 600,Folder);                      //
  HBook1F(Hist->fPprCheck ,"pprcheck"   ,Form("%s:True momentum CHECK of protons",Folder), 2000,  0, 50000,Folder);              //
  HBook1F(Hist->fPdeu     ,"pdeu"       ,Form("%s:generated p for deuteria; p [MeV/c]",Folder), 2400,  0, 600,Folder);                     //
  HBook1F(Hist->fPdeuCheck,"pdeucheck"  ,Form("%s:True momentum CHECK of deuteria",Folder), 2000,  0, 50000,Folder);             //
  HBook2F(Hist->fPrecoVsNSH ,"preco_vs_nsh" ,Form("%s:number of StrawHits vs generated p; p [MeV/c]; nSH",Folder), 2400, 0,600, 150, 0, 150, Folder); //%%
}

//-----------------------------------------------------------------------------
void TValidationModule2::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP[0]       ,"p"        ,Form("%s: Track P(Z1)"       ,Folder), 600,  35  ,600. ,Folder);
  HBook1F(Hist->fP[1]       ,"p_1"      ,Form("%s: Track P(total)[1]" ,Folder), 100, 104.5,105.5,Folder);
  HBook1F(Hist->fP[2]       ,"p_2"      ,Form("%s: Track P(total)[1]" ,Folder),1000,   0  ,600. ,Folder);
  HBook1F(Hist->fP0         ,"p0"       ,Form("%s: Track P(Z0)"       ,Folder),1000,   0  ,600. ,Folder);
  HBook1F(Hist->fP2         ,"p2"       ,Form("%s: Track P(z=-1540)"  ,Folder),1000,   0  ,400. ,Folder);
  HBook1D(Hist->fPDio       ,"pdio"     ,Form("%s: Track P(DIO WT)"   ,Folder), 400,  70  ,150. ,Folder);
  Hist->fPDio->Sumw2(kTRUE);

  HBook1F(Hist->fFitMomErr  ,"momerr"   ,Form("%s: Track FitMomError" ,Folder), 200,   0  ,  2.5 ,Folder);
  HBook1F(Hist->fPFront     ,"pf"       ,Form("%s: Track P(front)   " ,Folder), 400,   0  ,  110. ,Folder);
  HBook1F(Hist->fDpFront    ,"dpf"      ,Form("%s: Track P-P(front) " ,Folder), 200,  -5. ,  50. ,Folder);
  HBook1F(Hist->fDpFront0   ,"dp0f"     ,Form("%s: Track P0-P(front)" ,Folder), 200,  -5. ,  50. ,Folder);
  HBook1F(Hist->fDpFront2   ,"dp2f"     ,Form("%s: Track P2-P(front)" ,Folder), 200,  -5. ,  50. ,Folder);
  HBook1F(Hist->fPStOut     ,"pstout"   ,Form("%s: Track P(ST_Out)  " ,Folder), 400,   0  ,  110. ,Folder);
  HBook1F(Hist->fDpFSt      ,"dpfst"    ,Form("%s: Track Pf-Psto"     ,Folder), 200,  -5  ,  5. ,Folder);
  HBook2F(Hist->fDpFVsZ1    ,"dpf_vs_z1",Form("%s: Track DPF Vs Z1"   ,Folder), 200, -2000.,0,200,-5.,5,Folder);

  HBook1F(Hist->fPt         ,"pt"       ,Form("%s: Track Pt"          ,Folder), 600, 50,150,Folder);
  HBook1F(Hist->fCosTh      ,"costh"    ,Form("%s: Track cos(theta)"  ,Folder), 100,-1,1,Folder);
  HBook1F(Hist->fChi2       ,"chi2"     ,Form("%s: Track chi2 total"  ,Folder), 200, 0,400,Folder);
  HBook1F(Hist->fNDof       ,"ndof"     ,Form("%s: Number of DOF"     ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fChi2Dof    ,"chi2d"    ,Form("%s: track chi2/N(dof)" ,Folder), 600, 0, 20,Folder);
  HBook1F(Hist->fNActive    ,"nactv"    ,Form("%s: N(active)"         ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fT0         ,"t0"       ,Form("%s: track T0"          ,Folder), 200, 0,2000,Folder);
  HBook1F(Hist->fT0Err      ,"t0err"    ,Form("%s: track T0Err"       ,Folder), 100, 0,  20,Folder);
  HBook1F(Hist->fQ          ,"q"        ,Form("%s: track Q"           ,Folder),   4,-2,   2,Folder);
  HBook1F(Hist->fFitCons[0] ,"fcon"     ,Form("%s: track fit cons [0]",Folder), 200, 0,   2,Folder);
  HBook1F(Hist->fFitCons[1] ,"fcon1"    ,Form("%s: track fit cons [1]",Folder), 1000, 0,   0.2,Folder);
  HBook1F(Hist->fD0         ,"d0"       ,Form("%s: track D0      "    ,Folder), 200,-270, 270,Folder);
  HBook1F(Hist->fZ0         ,"z0"       ,Form("%s: track Z0      "    ,Folder), 200,-2070,2070,Folder);
  HBook1F(Hist->fTanDip     ,"tdip"     ,Form("%s: track tan(dip)"    ,Folder), 200, 0.0 ,3.0,Folder);
  HBook1F(Hist->fResid      ,"resid"    ,Form("%s: hit residuals"     ,Folder), 500,-0.5 ,0.5,Folder);
  HBook1F(Hist->fAlgMask    ,"alg"      ,Form("%s: algorithm mask"    ,Folder),  10,  0, 10,Folder);

  HBook1F(Hist->fDt         ,"dt"       ,Form("%s: track delta(T)"    ,Folder), 200,-40  ,20 ,Folder);
  HBook1F(Hist->fChi2Match  ,"chi2tcm"  ,Form("%s: chi2(t-c match)"   ,Folder), 250,  -200  ,350 ,Folder);

  HBook1F(Hist->fDx         ,"dx"       ,Form("%s: track delta(X)"    ,Folder), 200,-700 ,600,Folder);
  HBook1F(Hist->fDy         ,"dy"       ,Form("%s: track delta(Y)"    ,Folder), 200,-700 ,600,Folder);
  HBook1F(Hist->fDz         ,"dz"       ,Form("%s: track delta(Z)"    ,Folder), 200,-500 ,250,Folder);
  HBook1F(Hist->fDu         ,"du"       ,Form("%s: track-cluster DU)" ,Folder), 250,-500 ,300,Folder);
  HBook1F(Hist->fDv         ,"dv"       ,Form("%s: track-cluster DV)" ,Folder), 200,-200 ,200,Folder);
  HBook2F(Hist->fDvVsDu     ,"dv_vs_du" ,Form("%s: Track Dv Vs Du"    ,Folder), 100, -250,250,100,-100.,100,Folder);
  HBook1F(Hist->fPath       ,"path"     ,Form("%s: track sdisk"       ,Folder),  50,   0 ,500,Folder);
  HBook2F(Hist->fDuVsPath   ,"du_vs_path",Form("%s: Track Du Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDucVsPath  ,"duc_vs_path",Form("%s: T-C Duc Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDvVsPath   ,"dv_vs_path",Form("%s: T-C  Dv Vs Path"  ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDvcVsPath  ,"dvc_vs_path",Form("%s: T-C Dvc Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDtVsPath   ,"dt_vs_path",Form("%s: T-C DT Vs Path"   ,Folder),  50,   0 ,500,100,  -5.,  5.,Folder);
  HBook2F(Hist->fDuVsTDip   ,"du_vs_tdip",Form("%s: Track Du Vs TDip" ,Folder), 100, 0.5 ,1.5,200,-200.,200.,Folder);
  HBook2F(Hist->fDvVsTDip   ,"dv_vs_tdip",Form("%s: Track Dv Vs TDip" ,Folder), 100, 0.5 ,1.5,200,-200.,200.,Folder);

  HBook1F(Hist->fZ1         ,"z1"       ,Form("%s: track Z1      "    ,Folder), 200,-2000,2000,Folder);
  HBook1F(Hist->fNClusters  ,"ncl"      ,Form("%s: track N(clusters)" ,Folder),  10, 0   , 10,Folder);
  HBook1F(Hist->fVaneID     ,"vid"      ,Form("%s: track vane ID"     ,Folder),  10,-5   ,  5,Folder);
  HBook1F(Hist->fXCal       ,"xcal"     ,Form("%s: track XCal"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fYCal       ,"ycal"     ,Form("%s: track YCal"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fZCal       ,"zcal"     ,Form("%s: track ZCal"        ,Folder), 200, 1500,3500,Folder);
  HBook1F(Hist->fXTrk       ,"xtrk"     ,Form("%s: track XTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fYTrk       ,"ytrk"     ,Form("%s: track YTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fRTrk       ,"rtrk"     ,Form("%s: track RTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fZTrk       ,"ztrk"     ,Form("%s: track ZTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fECl        ,"ecl"      ,Form("%s: cluster E"         ,Folder), 300, 0   ,150,Folder);
  HBook1F(Hist->fEClEKin    ,"ecl_ekin" ,Form("%s: cluster E/Ekin(mu)",Folder), 200, 0   ,2,Folder);
  HBook1F(Hist->fEp         ,"ep"       ,Form("%s: track E/P"         ,Folder), 300, 0   ,1.5,Folder);
  HBook2F(Hist->fEpVsPath   ,"ep_vs_path",Form("%s: E/P Vs Path"      ,Folder),  50,   0 ,500,150,  0.,  1.5,Folder);
  HBook2F(Hist->fNHVsStation,"nh_vs_st" ,Form("%s: N(hits) Vs Station",Folder),  40, 0,40,10,-0.5,9.5,Folder);
  HBook2F(Hist->fNHVsNSt    ,"nh_vs_nst",Form("%s: N(hits) Vs NSt"    ,Folder),  10,-0.5,9.5,40,-0.5,39.5,Folder);

  HBook1F(Hist->fRSlope     ,"rslope"   ,Form("%s: Res Slope"         ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fXSlope     ,"xslope"   ,Form("%s: Res/Sig Slope"     ,Folder), 200,-20 , 20,Folder);

  HBook2F(Hist->fEpVsDt     ,"ep_vs_dt" ,Form("%s: E/P vs Dt"         ,Folder), 200, -10, 10,150,0.,1.5,Folder);
  HBook1F(Hist->fEleLogLHCal,"ele_llh_c",Form("%s: ELE Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
  HBook1F(Hist->fMuoLogLHCal,"muo_llh_c",Form("%s: MUO Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
  HBook1F(Hist->fLogLHRCal  ,"llhr_cal" ,Form("%s: LogLH(e/m) Cal"    ,Folder), 200,-100,100,Folder);
  HBook1F(Hist->fLogLHRDeDx ,"llhr_dedx",Form("%s: LogLH(e/m) De/Dx"  ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fLogLHRXs   ,"llhr_xs"  ,Form("%s: LogLH(e/m) XSlope" ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fLogLHRTrk  ,"llhr_trk" ,Form("%s: LogLH(e/m) Trk"    ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fLogLHR     ,"llhr"     ,Form("%s: LogLH(e/m)"        ,Folder), 200,-100 ,100,Folder);

  HBook1F(Hist->fPdgCode    ,"pdg"      ,Form("%s: track PDG code"    ,Folder),2350  ,-100,2250,Folder);
  HBook1F(Hist->fFrGH       ,"fgh"      ,Form("%s: Fraction Goog Hits",Folder), 100, -5,5,Folder);

  HBook2F(Hist->fNEPlVsNHPl ,"nep_vs_nhp",Form("%s: Track NEXP vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
  HBook2F(Hist->fNDPlVsNHPl ,"ndp_vs_nhp",Form("%s: Track NDIF vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
  HBook2F(Hist->fChi2dVsNDPl,"chi2d_vs_ndp",Form("%s: Track Chi2/Dof vs NDP",Folder), 30, 0,30,100,0.,10,Folder);
  HBook2F(Hist->fDpFVsNDPl  ,"dpf_vs_ndp"  ,Form("%s: Track DpF vs NDP",Folder)     , 30, 0,30,100,-5,5,Folder);

  HBook1F(Hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);
  /*
  HBook1F(Hist->fDeuPreco  ,"deupreco"    ,Form("%s:True momentum of reco_deutons",Folder), 2000,  0, 500,Folder);                        //Deuton
  HBook2F(Hist->fDeuPrecoVsP ,"deupreco_vs_p" ,Form("%s:Reco_mom vs true for deutons",Folder), 2000, 0,500, 2000, 0,500, Folder);         //Deuton
  */

  /*
  HBook1F(Hist->fPrecoprAll  ,"precoprAll"    ,Form("%s: generated p of reco_particles All ; p [MeV/c]",Folder), 2400,  0, 600,Folder);                          //
  HBook2F(Hist->fPrecoVsPAll ,"preco_vs_pAll" ,Form("%s: reconstructed vs generated p All; p [MeV/c]; p [MeV/c]",Folder), 2400, 0,600, 2400, 0,600, Folder);                           //
  HBook2F(Hist->fChi2dVsPAll ,"chi2d_vs_pAll" ,Form("%s:Chi2/dof vs generated momentum All; p [MeV/c]; #chi2/dof",Folder), 2400, 0,600, 200, 0,50, Folder);                    //
  */

  HBook1F(Hist->fPrecopr  ,"precopr"    ,Form("%s: generated p of reco_protons; p [MeV/c]",Folder), 2400,  0, 600,Folder);                          //
  HBook2F(Hist->fPrecoVsP ,"preco_vs_p" ,Form("%s: reconstructed vs generated p; p [MeV/c]; p [MeV/c]",Folder), 2400, 0,600, 2400, 0,600, Folder);                           //
  HBook2F(Hist->fChi2dVsP ,"chi2d_vs_p" ,Form("%s:Chi2/dof vs generated momentum; p [MeV/c]; #chi2/dof",Folder), 2400, 0,600, 200, 0,50, Folder);                    //

  HBook2F(Hist->fPrecoVsTanTh ,"preco_vs_tanth" ,Form("%s: TanTh vs generated p; p [MeV/c] ",Folder), 2400, 0,600, 2000, 0 ,10, Folder);                    //
  HBook2F(Hist->fPrecoVsSHE ,"preco_vs_shE" ,Form("%s: SHE vs generated p; p [MeV/c]; SHE [keV] ",Folder), 2400, 0,600, 200, 0, 0.015, Folder);                        //
  HBook2F(Hist->fTanThVsSHE ,"tanth_vs_shE" ,Form("%s: TanTh vs SHE",Folder), 2000, 0, 10, 200, 0, 0.015, Folder);                         //
  HBook2F(Hist->fPrecoVsNSH ,"preco_vs_nsh" ,Form("%s: generated p vs number of SH; p [MeV/c]; nSH",Folder), 2400, 0,600, 150, 0, 150, Folder);                 //
  HBook2F(Hist->fFracSHVsPreco ,"fracsh_vs_preco" ,Form("%s: fraction of SH vs generated p; p [MeV/c]; nSH active/nSH",Folder), 2400, 0,600, 200, 0, 1.1, Folder);         //

  HBook2F(Hist->fLastTrkZVsPreco ,"lasttrkz_vs_preco" ,Form("%s: z last SH vs generated p; p[MeV/c] ;station",Folder),2400, 0,600, 34, -2, 32, Folder);     //
  HBook2F(Hist->fFirstTrkZVsPreco ,"firsttrkz_vs_preco" ,Form("%s: z first SH vs generated p; p[MeV/c] ;station",Folder),2400, 0,600, 34, -2, 32, Folder);  //

  HBook2F(Hist->fUnusedSHVsPreco ,"unusedsh_vs_preco" ,Form("%s:Unused SH vs true P",Folder), 2400, 0,600, 2000, 0,100000, Folder);       //
  HBook2F(Hist->fUsedSHVsPreco ,"usedsh_vs_preco" ,Form("%s:Used SH vs true P",Folder),  2400, 0,600, 2000, 0,100000, Folder);             //

}

//-----------------------------------------------------------------------------
void TValidationModule2::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
  HBook1F(Hist->fEleMom    ,"ce_mom"   ,Form("%s: Conversion Electron Momentum"    ,Folder),1000,  0,200,Folder);
  HBook1D(Hist->fDioMom    ,"dio_mom"  ,Form("%s: DIO momentum"                    ,Folder),1000, 50,150,Folder);
  HBook1F(Hist->fRv         ,"rv"      ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 1000,Folder);
  HBook1F(Hist->fZv         ,"zv"      ,Form("%s: Z(Vertex)"                       ,Folder), 300, 0,15000,Folder);
  HBook1F(Hist->fNClusters ,"ncl"      ,Form("%s: Number of Reconstructed Clusters",Folder),200,0,200,Folder);
  HBook1F(Hist->fNTracks   ,"ntrk"     ,Form("%s: Number of Reconstructed Tracks, ntrk"  ,Folder),50,0,50,Folder);

  HBook1F(Hist->fNTracksP   ,"ntrk_p"     ,Form("%s: Number of Deuteron Reconstructed Tracks; ntrk"  ,Folder),50,0,50,Folder);
  HBook1F(Hist->fNTracksD   ,"ntrk_d"     ,Form("%s: Number of Proton Reconstructed Tracks; ntrk"  ,Folder),50,0,50,Folder);
  HBook1F(Hist->fNTracksO   ,"ntrk_o"     ,Form("%s: Number of Other Reconstructed Tracks; ntrk"  ,Folder),50,0,50,Folder);
  HBook1F(Hist->fNTracksCut   ,"ntrk_cut"     ,Form("%s: Number of Reconstructed Tracks (Cut), ntrk"  ,Folder),50,0,50,Folder);
  HBook1F(Hist->fNTracksPCut   ,"ntrk_p_cut"     ,Form("%s: Number of Deuteron Reconstructed Tracks (Cut); ntrk"  ,Folder),50,0,50,Folder);
  HBook1F(Hist->fNTracksDCut   ,"ntrk_d_cut"     ,Form("%s: Number of Proton Reconstructed Tracks (Cut); ntrk"  ,Folder),50,0,50,Folder);
  HBook1F(Hist->fNTracksOCut   ,"ntrk_o_cut"     ,Form("%s: Number of Other Reconstructed Tracks (Cut); ntrk"  ,Folder),50,0,50,Folder);


  HBook1F(Hist->fNGoodTracks  ,"ngtrk" ,Form("%s: Number of Good          Tracks"  ,Folder),100,0,100,Folder);
  HBook1F(Hist->fNStrawHits[0],"nsh_0" ,Form("%s: Number of Straw Hits [0]"        ,Folder),250,0,250,Folder);
  HBook1F(Hist->fNStrawHits[1],"nsh_1" ,Form("%s: Number of Straw Hits [1]"        ,Folder),250,0,5000,Folder);
  HBook1F(Hist->fNGoodSH   ,"nsh50"    ,Form("%s: N(SH) +/-50"                     ,Folder),300,0,1500,Folder);
  HBook1F(Hist->fDtClT     ,"dt_clt"   ,Form("%s: DT(cluster-track)"               ,Folder),100,-100,100,Folder);
  HBook1F(Hist->fDtClS     ,"dt_cls"   ,Form("%s: DT(cluster-straw hit)"           ,Folder),200,-200,200,Folder);
  HBook1F(Hist->fSHTime    ,"shtime"   ,Form("%s: Straw Hit Time"                  ,Folder),400,0,2000,Folder);
  HBook1F(Hist->fEMax      ,"emax"     ,Form("%s: Max cluster energy"              ,Folder),150,0,150,Folder);
  HBook1F(Hist->fNHyp      ,"nhyp"     ,Form("%s: N(fit hypotheses)"               ,Folder),5,0,5,Folder);
  HBook1F(Hist->fBestHyp[0],"bfh0"     ,Form("%s: Best Fit Hyp[0](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(Hist->fBestHyp[1],"bfh1"     ,Form("%s: Best Fit Hyp[1](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  
  HBook1F(Hist->fNTimeClusters ,"ntimecl"      ,Form("%s: Number of Time Clusters",Folder),200,0,200,Folder);
  HBook2F(Hist->fNTimeClustersVsMom ,"ntimecl_vs_mom"      ,Form("%s: Number of Time Clusters VS Mom",Folder),2400,0,600,20,0,20,Folder);

  HBook1F(Hist->fNGenp     ,"ngenp"    ,Form("%s: N(Gen Particles)"                ,Folder),500,0,500,Folder);
  HBook1F(Hist->fQSH_p     ,"qsh_p"      ,Form("%s: Charge of SH from proton"              ,Folder),500,0,0.014,Folder);  
  HBook1F(Hist->fQSH_d     ,"qsh_d"      ,Form("%s: Charge of SH from deuteron"              ,Folder),500,0,0.014,Folder);                       
  HBook1F(Hist->fQSH_e     ,"qsh_e"      ,Form("%s: Charge of SH from electron"              ,Folder),500,0,0.014,Folder);
  HBook1F(Hist->fQSH       ,"qsh"      ,Form("%s: Charge for each SH"              ,Folder),500,0,0.014,Folder);                          // 
  HBook2F(Hist->fNSHVsPreco ,"nsh_vs_preco" ,Form("%s:Number of SH  vs true(?) momentum",Folder), 2400, 0,600, 200,0,200, Folder);        //
  HBook2F(Hist->fLumVsNTrk ,"istlum_vs_ntrk" ,Form("%s:Istant Luminosity vs N of Tracks",Folder), 20, 0,20, 2000, 1e6,100*1e6, Folder);   //
  HBook2F(Hist->fLumVsNTrkCut ,"istlum_vs_ntrk_cut" ,Form("%s:Istant Luminosity vs N of Tracks (cut)",Folder), 20, 0,20, 2000, 1e6,100*1e6, Folder);   //

  HBook2F(Hist->fLastZVsPreco ,"lastz_vs_preco" ,Form("%s:LastZ vs generated momentum",Folder), 2400, 0,600, 34, -2, 32, Folder);           //
  HBook2F(Hist->fFirstZVsPreco ,"firstz_vs_preco" ,Form("%s:FirstZ vs generated momentum",Folder),2400, 0,600, 34, -2, 32, Folder);        //
  HBook2F(Hist->fZVsPreco ,"z_vs_preco" ,Form("%s:FirstZ vs generated momentum",Folder),2400, 0,600, 34, -2, 32, Folder);        //

  HBook2F(Hist->fStdTimeVsMom     ,"stdtime_vs_mom"    ,Form("%s: cluster time std vs mom",Folder),2400, 0,600,500,0,50,Folder);
  HBook2F(Hist->fWidthTimeVsMom     ,"widthtime_vs_mom"    ,Form("%s: cluster time width vs mom",Folder),2400, 0,600,500,0,50,Folder);


}

//-----------------------------------------------------------------------------
void TValidationModule2::BookSimpHistograms(SimpHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fPdgCode   ,"pdg"         ,Form("%s: PDG code"                     ,Folder),200,-10000,10000,Folder);
  HBook1F(Hist->fNStrawHits,"nsth"        ,Form("%s: n straw hits"                 ,Folder),200,   0,200,Folder);
  HBook1F(Hist->fMomTargetEnd    ,"ptarg" ,Form("%s: Mom after Stopping Target" ,Folder),2400, 0,600,Folder);
  HBook1F(Hist->fMomTrackerFront ,"pfront",Form("%s: Mom at the Tracker Front"  ,Folder),2400, 0,600,Folder);

  HBook1F(Hist->fPvd      ,"pvd"       ,Form("%s: Momentum_VD"     ,Folder), 2400,  0, 600,Folder);
  HBook1F(Hist->fCosTh  ,"cos_th"  ,Form("%s: Cos(Theta)_VD"   ,Folder), 200,   -1.,   1.,Folder);
  HBook2F(Hist->fPvdVsNSH ,"pvd_vs_nsh" ,Form("%s:Momentum_VD vs number of SH",Folder), 2400, 0,600, 150, 0, 150, Folder); //%%
  HBook2F(Hist->fPvdVsCosTh ,"pvd_vs_cos_th" ,Form("%s:Momentum_VD vs Cos(Theta)",Folder), 2400, 0,600, 400, -1, 1, Folder); //%%
  HBook2F(Hist->fCosThVsNSH ,"cos_th_vs_nsh" ,Form("%s:Cos(Theta)_VD vs number of SH",Folder), 400, -1, 1, 150, 0, 150, Folder); //%%
  HBook2F(Hist->fPvdVsPgen ,"pvd_vs_pgen" ,Form("%s: VD vs generated momentum",Folder), 2400, 0,600,2400, 0,600, Folder); //%%
  HBook1F(Hist->fPvdEndT      ,"pvd_end"       ,Form("%s: Momentum_End Tracker"     ,Folder), 2400,  0, 600,Folder);

}
//_____________________________________________________________________________
void TValidationModule2::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");

//-----------------------------------------------------------------------------
// book crystal histograms
//-----------------------------------------------------------------------------
  HBook1F(fHist.fCrystalR[0],"rc_0"     ,Form("disk [0] crystal radius"),100,0,1000,"Hist");
  HBook1F(fHist.fCrystalR[1],"rc_1"     ,Form("disk [1] crystal radius"),100,0,1000,"Hist");

//--------------------------------------------------------------------------------
// book timecluster histograms
//--------------------------------------------------------------------------------
  int book_timecluster_histset[kNTimeClusterHistSets];
  for (int i=0; i<kNTimeClusterHistSets; ++i)  book_timecluster_histset[i] = 0;

  book_timecluster_histset[0] = 1;   // all events
  book_timecluster_histset[1] = 1;   // timeclusters with NHits>5
  book_timecluster_histset[2] = 1;   // timeclusters with NHits>10
  book_timecluster_histset[3] = 1;   // timeclusters with NHits>15
  book_timecluster_histset[4] = 1;   // No helices


   for (int i=0; i<kNTrackSeedHistSets; i++) {
    if (book_timecluster_histset[i] != 0) {
      sprintf(folder_name,"timecluster_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTimeCluster[i] = new TimeClusterHist_t;
      BookTimeClusterHistograms(fHist.fTimeCluster[i],Form("Hist/%s",folder_name));
    }
  }


//--------------------------------------------------------------------------------
// book trackSeed histograms
//--------------------------------------------------------------------------------
  int book_trackSeed_histset[kNTrackSeedHistSets];
  for (int i=0; i<kNTrackSeedHistSets; ++i)  book_trackSeed_histset[i] = 0;

  book_trackSeed_histset[0] = 1;   // events with at least one trackSeed
  book_trackSeed_histset[1] = 1;   // events with at least one trackSeed with p > 150 MeV/c
  book_trackSeed_histset[2] = 1;   // events with at least one trackSeed with nhits >= 10
  book_trackSeed_histset[3] = 1;   // events with at least one trackSeed with nhits >= 10 and chi2<5

   for (int i=0; i<kNTrackSeedHistSets; i++) {
    if (book_trackSeed_histset[i] != 0) {
      sprintf(folder_name,"trkseed_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrackSeed[i] = new TrackSeedHist_t;
      BookTrackSeedHistograms(fHist.fTrackSeed[i],Form("Hist/%s",folder_name));
    }
  }
  
//--------------------------------------------------------------------------------
// book helix histograms
//--------------------------------------------------------------------------------
  int book_helix_histset[kNHelixHistSets];
  for (int i=0; i<kNHelixHistSets; ++i)  book_helix_histset[i] = 0;

  book_helix_histset[0] = 1;   // events with at least one helix
  book_helix_histset[1] = 1;   // events with at least one helix with p > 250 MeV/c
  book_helix_histset[2] = 1;   // events with at least one helix with p < 250 MeV/c
  book_helix_histset[3] = 1;   // events with at least one helix with nhits >= 10
  book_helix_histset[4] = 1;   // events with at least one helix with nhits >= 10 and chi2<5 (XY and PhiZ)
  book_helix_histset[5] = 1;   // events with no tracks
  book_helix_histset[6] = 1;   // events protons
  book_helix_histset[7] = 1;   // events deuteron


   for (int i=0; i<kNHelixHistSets; i++) {
    if (book_helix_histset[i] != 0) {
      sprintf(folder_name,"helix_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fHelix[i] = new HelixHist_t;
      BookHelixHistograms(fHist.fHelix[i],Form("Hist/%s",folder_name));
    }
  }
  

//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;	        // events with a reconstructed track
  book_event_histset[ 2] = 1;	        // events without reconstructed tracks
  book_event_histset[ 6] = 1;	        // events with tracks passing "Set C" cuts

  book_event_histset[10] = 0;	        // 
  book_event_histset[11] = 1;	        // selection cuts
  book_event_histset[12] = 1;	        // 
  book_event_histset[13] = 1;	        // 
  book_event_histset[14] = 1;	        // 
  book_event_histset[15] = 1;	        // 
  book_event_histset[16] = 1;	        // 
  book_event_histset[17] = 1;	        // 
  book_event_histset[18] = 1;	        // 
  book_event_histset[19] = 0;	        // 
  book_event_histset[20] = 0;	        // 
  book_event_histset[21] = 0;	        // 
  book_event_histset[22] = 0;	        // 
  book_event_histset[23] = 0;	        // 
					// TrkPatRec tracks
  book_event_histset[24] = 1;	        // events with at least one reco track
  book_event_histset[25] = 1;	        // 
  book_event_histset[26] = 1;	        // 
  book_event_histset[27] = 1;	        // 
  book_event_histset[28] = 1;	        // 


  book_event_histset[40] = 1;	        // p or d



  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book simp histograms
//-----------------------------------------------------------------------------
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[ 0] = 1;		// all events
  
  for (int i=0; i<kNSimpHistSets; i++) {
    if (book_simp_histset[i] != 0) {
      sprintf(folder_name,"sim_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fSimp[i] = new SimpHist_t;
      BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track histograms
//-----------------------------------------------------------------------------
  int book_track_histset[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

  book_track_histset[  0] = 1;		// all tracks
  book_track_histset[  1] = 1;		// all tracks passing Set C cuts 

  book_track_histset[  6] = 1;		// all tracks E/P > 0.4

  book_track_histset[ 11] = 1;		// all tracks with P > 150

  book_track_histset[ 14] = 1;		// tracks with fcons < 1.e-2
  book_track_histset[ 15] = 1;		// tracks intersecting the 1st disk
  book_track_histset[ 16] = 1;		// tracks intersecting the 2nd disk
  book_track_histset[ 17] = 1;		// tracks with no calorimeter intersections
  
  book_track_histset[ 20] = 1;		// tracks with Nhits >= 20
  book_track_histset[ 21] = 1;		// tracks with Nhits >= 20 and chi/Ndof < 5

  book_track_histset[ 28] = 1;		// NO Set C tracks, E/P > 1.1
  book_track_histset[ 29] = 1;		// NO Set C tracks, E/P > 0, P > 100 *precursor for TRK_25*

  book_track_histset[ 30] = 1;		// tracks with Nhits >= 25
  book_track_histset[ 31] = 1;		// tracks with Nhits >= 25 and chi/Ndof < 3

  book_track_histset[ 40] = 1;		// all tracks, alg_mask = 1
  book_track_histset[ 41] = 1;		// Set "C" tracks, alg_mask = 1
  book_track_histset[ 42] = 1;		// Set "C" tracks, alg_mask = 1, T > 700

  book_track_histset[ 50] = 1;		// all tracks, alg_mask = 2
  book_track_histset[ 51] = 1;		// Set "C" tracks, alg_mask = 2
  book_track_histset[ 52] = 1;		// Set "C" tracks, alg_mask = 2, T > 700

  book_track_histset[ 60] = 1;		// all tracks, alg_mask = 3
  book_track_histset[ 61] = 1;		// Set "C" tracks, alg_mask = 3
  book_track_histset[ 62] = 1;		// Set "C" tracks, alg_mask = 3, T > 700

  //
  book_track_histset[ 63] = 1;          // all tracks with P>270                        <---------------
  book_track_histset[ 64] = 1;          // all tracks with P<270                        <---------------
  book_track_histset[ 71] = 1;          // all track Protons                            <---------------
  book_track_histset[ 72] = 1;          // all track Deuterons                          <---------------
  //
  
  for (int i=0; i<kNTrackHistSets; i++) {
    if (book_track_histset[i] != 0) {
      sprintf(folder_name,"trk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrack[i] = new TrackHist_t;
      BookTrackHistograms(fHist.fTrack[i],Form("Hist/%s",folder_name));
    }
  
  }
//-----------------------------------------------------------------------------
// book cluster histograms
//-----------------------------------------------------------------------------
  int book_cluster_histset[kNClusterHistSets];
  for (int i=0; i<kNClusterHistSets; i++) book_cluster_histset[i] = 0;

  book_cluster_histset[0] = 1;		// all clusters
  book_cluster_histset[1] = 1;		// clusters in events with the reconstructed e-
  book_cluster_histset[2] = 1;		// clusters in events with the track passing SetC cuts
  book_cluster_histset[3] = 1;		// cNO lusters in events w/track passing SetC cuts and |dt|<2.5ns 
  book_cluster_histset[4] = 1;		// clusters > 10 MeV
  book_cluster_histset[5] = 1;		// clusters > 60 MeV
  book_cluster_histset[6] = 1;		// clusters disk#0
  book_cluster_histset[7] = 1;		// clusters disk#1

  for (int i=0; i<kNClusterHistSets; i++) {
    if (book_cluster_histset[i] != 0) {
      sprintf(folder_name,"cls_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCluster[i] = new ClusterHist_t;
      BookClusterHistograms(fHist.fCluster[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book Genp histograms
//-----------------------------------------------------------------------------
  int book_genp_histset[kNGenpHistSets];
  for (int i=0; i<kNGenpHistSets; i++) book_genp_histset[i] = 0;

  book_genp_histset[0] = 1;		// all particles
  book_genp_histset[1] = 1;		// events with tracks
  book_genp_histset[2] = 1;		// events with good tracks
  book_genp_histset[3] = 1;		// At least 1 SH
  book_genp_histset[4] = 1;		// At least 5 SH
  book_genp_histset[5] = 1;		// At least 20 SH



  for (int i=0; i<kNGenpHistSets; i++) {
    if (book_genp_histset[i] != 0) {
      sprintf(folder_name,"gen_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fGenp[i] = new GenpHist_t;
      BookGenpHistograms(fHist.fGenp[i],Form("Hist/%s",folder_name));
    }
  }

}

//-----------------------------------------------------------------------------
// need MC truth branch
//-----------------------------------------------------------------------------
void TValidationModule2::FillEventHistograms(EventHist_t* Hist) {
  double            cos_th(-2.), dio_wt(-1.), xv(0.), yv(0.), rv(0.), zv(0.), p(0);
  TLorentzVector    mom;

  if (fParticle) {
    fParticle->Momentum(mom);
    p      = mom.P();

    cos_th = mom.Pz()/p;
    dio_wt = TStntuple::DioWeightAl(p);

    xv = fParticle->Vx()+3904.;
    yv = fParticle->Vy();
    rv = sqrt(xv*xv+yv*yv);
    zv = fParticle->Vz();
  }

  Hist->fEleMom->Fill(p);
  Hist->fDioMom->Fill(p,dio_wt);
  Hist->fEleCosTh->Fill(cos_th);
  Hist->fRv->Fill(rv);
  Hist->fZv->Fill(zv);

  Hist->fNClusters->Fill(fNClusters);
  Hist->fNTracks->Fill  (fNTracks[0]);
  Hist->fNGoodTracks->Fill(fNGoodTracks);
  Hist->fNStrawHits[0]->Fill(fNStrawHits);
  Hist->fNStrawHits[1]->Fill(fNStrawHits);

  Hist->fNTimeClusters->Fill(fNTimeClusters[0]);
  Hist->fNTimeClustersVsMom->Fill(p,fNTimeClusters[0]);
  Hist->fNSHVsPreco->Fill(p,fNStrawHits);


  Hist->fLumVsNTrk->Fill(fTrackBlock->NTracks(),GetHeaderBlock()->InstLum());
  int ntrkchi2=0;
  int ntrk_p = 0;
  int ntrk_d = 0;
  int ntrk_o = 0;
  int ntrk_p_cut = 0;
  int ntrk_d_cut = 0;
  int ntrk_o_cut = 0;
  bool nice;
  for(int i=0;i<fTrackBlock->NTracks();i++) {
    nice = false;
    TStnTrack* trk = fTrackBlock->Track(i);
    if(trk->Chi2Dof()<5) {ntrkchi2=ntrkchi2+1; nice = true;}//---------------------------------
    if(trk->fPdgCode == 2212) {ntrk_p = ntrk_p + 1; if(nice) ntrk_p_cut = ntrk_p_cut+1; }
    else if(trk->fPdgCode == 1000010020) {ntrk_d = ntrk_d + 1; if(nice) ntrk_d_cut = ntrk_d_cut+1;}
    else { ntrk_o = ntrk_o + 1;if(nice) ntrk_o_cut = ntrk_o_cut+1;}

  }
  Hist->fNTracksCut->Fill(ntrkchi2);
  Hist->fLumVsNTrkCut->Fill(ntrkchi2,GetHeaderBlock()->InstLum()); // ------------------------------------------------------------fNTracks[21]
  Hist->fNTracksP->Fill(ntrk_p); //protons
  Hist->fNTracksD->Fill(ntrk_d); //deuterons
  Hist->fNTracksO->Fill(ntrk_o); //others?

  Hist->fNTracksPCut->Fill(ntrk_p_cut); //protons
  Hist->fNTracksDCut->Fill(ntrk_d_cut); //deuterons
  Hist->fNTracksOCut->Fill(ntrk_o_cut); //others?


  double emax   = -1;
  double t0_cls = -1;
  double dt     = 9999.;

  TStnCluster* cluster(0);
  if (fNClusters > 0) cluster = fClusterBlock->Cluster(0);

  TStnTrack* track(0);
  if (fNTracks[0] > 0) track = fTrackBlock->Track(0);

  if (cluster) {
    emax   = cluster->Energy();
    t0_cls = cluster->Time();
  }

  double t0_trk = -1;
  if (track) {
    t0_trk = track->fT0;
  }

  if (track && cluster) {
    dt = t0_cls-t0_trk;
  }

  Hist->fDtClT->Fill(dt);
  Hist->fEMax->Fill(emax);

  //  SH Energy and max&min Station
  TStrawHitData*  sh;
  int n_good_hits = 0;
  int minStation = 30;
  int maxStation = -2;

  float fSTDtime_sum=0;
  float stdtime=0;
  float tfirst=2000;
  float tlast=0;
  for (int i=0; i<fNStrawHits; i++ ) {
    sh = fStrawDataBlock->Hit(i);
    Hist->fZVsPreco->Fill(p,sh->Station());  // bastiano proton eff

    dt = t0_cls-sh->Time() + 15;
    Hist->fDtClS->Fill(dt);
    Hist->fSHTime->Fill(sh->Time());

    Hist->fQSH->Fill(sh->Energy());


    if(sh->PdgCode()==2212 && sh->GeneratorCode()==28){                 //JJ
      Hist->fQSH_p->Fill(sh->Energy()); //proton SH chargesh->Energy()
    }
    else if(sh->PdgCode()==1000010020 && sh->GeneratorCode()==28){                 //JJ
      Hist->fQSH_d->Fill(sh->Energy()); //deuteron SH chargesh->Energy()
    }
    else if(sh->PdgCode()==11){   
      Hist->fQSH_e->Fill(sh->Energy());
      //if(sh->Energy()>0.0055) sh->Print();
    } 

    if(sh->Station()>maxStation) maxStation=sh->Station();
    if(sh->Station()<minStation) minStation=sh->Station();
    if (fabs(dt+15.)< 50) n_good_hits += 1;
    fSTDtime_sum=fSTDtime_sum+sh->Time();  

    if(sh->Time()>tlast) tlast = sh->Time();
    if(sh->Time()<tfirst) tfirst = sh->Time();
  }
  fSTDtime_sum=fSTDtime_sum/fNStrawHits;  
  for (int i=0; i<fNStrawHits; i++ ) { 
        sh = fStrawDataBlock->Hit(i);
	stdtime=stdtime+(sh->Time()-fSTDtime_sum)*(sh->Time()-fSTDtime_sum);
     }
  stdtime=sqrt(stdtime/(fNStrawHits-1)); 
  Hist->fStdTimeVsMom->Fill(p,stdtime);
  Hist->fWidthTimeVsMom->Fill(p,tlast-tfirst);

  Hist->fLastZVsPreco->Fill(p,maxStation);
  Hist->fFirstZVsPreco->Fill(p,minStation);
  Hist->fNGoodSH->Fill(n_good_hits);

  Hist->fNHyp->Fill(fNHyp);
  Hist->fBestHyp[0]->Fill(fBestHyp[0]);
  Hist->fBestHyp[1]->Fill(fBestHyp[1]);

  Hist->fNGenp->Fill(fNGenp);  
}

//--------------------------------------------------------------------------------
// function to fill TrasckSeedHit block
//--------------------------------------------------------------------------------
void TValidationModule2::FillTrackSeedHistograms(TrackSeedHist_t*   Hist, TStnTrackSeed*    TrkSeed){
  
  int         nhits    = TrkSeed->NHits      ();
  double      clusterT = TrkSeed->ClusterTime();
  double      clusterE = TrkSeed->ClusterEnergy();
  
  double      mm2MeV   = 3/10.;
  double      pT       = TrkSeed->Pt();
  double      radius   = pT/mm2MeV;

  double      tanDip   = TrkSeed->TanDip();  
  double      p        = pT/std::cos( std::atan(tanDip));
  

  Hist->fNHits      ->Fill(nhits);	 
  Hist->fClusterTime->Fill(clusterT);
  Hist->fClusterEnergy->Fill(clusterE);
  
  Hist->fRadius     ->Fill(radius);    
  Hist->fMom        ->Fill(p);	 
  Hist->fPt         ->Fill(pT);	 
  Hist->fTanDip     ->Fill(tanDip);    
  
  Hist->fChi2       ->Fill(TrkSeed->Chi2());
  Hist->fFitCons    ->Fill(TrkSeed->FitCons());
  Hist->fD0         ->Fill(TrkSeed->D0());

}

//--------------------------------------------------------------------------------
// function to fill Helix block
//--------------------------------------------------------------------------------
void TValidationModule2::FillTimeClusterHistograms(TimeClusterHist_t*   Hist, TStnTimeCluster*    TimeCluster){
  
  int         nhits      = TimeCluster->NHits      ();
  int         ncombohits = TimeCluster->NComboHits ();

  double      time       = TimeCluster->T0();
  double      clusterE   = TimeCluster->ClusterEnergy();

  Hist->fNHits         ->Fill(nhits);	 
  Hist->fNComboHits    ->Fill(ncombohits);	 
  Hist->fT0            ->Fill(time);
  Hist->fClusterEnergy ->Fill(clusterE);

  if(GetHeaderBlock()->InstLum() == -1 && fSimpBlock){
    TSimParticle* simp = fSimpBlock->Particle(0);
    double p  = simp->StartMom()->P();
    float frac = (float) nhits/( fStrawDataBlock->NHits() );
    std::cout<< nhits<<" "<<fStrawDataBlock->NHits()<<" "<<frac<<std::endl;
    if(frac < 0.4) GetHeaderBlock()->Print();
    Hist->fNComboHitsVsP ->Fill(p,ncombohits);
    Hist->fFracSHVsP->Fill(p,frac);
  }
  

}

//--------------------------------------------------------------------------------
// function to fill Helix block
//--------------------------------------------------------------------------------
void TValidationModule2::FillHelixHistograms(HelixHist_t*   Hist, TStnHelix*    Helix){
  
  int         nhits    = Helix->NHits      ();
  double      clusterT = Helix->ClusterTime();
  double      clusterE = Helix->ClusterEnergy();
  
  double      radius   = Helix->Radius();

  double      lambda   = Helix->Lambda();  
  double      tanDip   = lambda/radius;
  double      mm2MeV   = 3/10.;
  double      pT       = radius*mm2MeV;
  double      p        = pT/std::cos( std::atan(tanDip));
  
  double      nRot     = 3270/(abs(lambda)*TMath::TwoPi()); //

  Hist->fNHits         ->Fill(nhits);	 
  Hist->fClusterTime   ->Fill(clusterT);
  Hist->fClusterEnergy ->Fill(clusterE);
  
  Hist->fRadius        ->Fill(radius);    
  Hist->fMom           ->Fill(p);	 
  Hist->fPt            ->Fill(pT);	 
  Hist->fLambda        ->Fill(-lambda); // -lambda because lambda is negative for protons    
  //  Hist->fAlg           ->Fill(Helix->AlgorithmID()); // does this exist?
  Hist->fD0            ->Fill(Helix->D0());

  double p_sim;
  if(fSimpBlock){
    TSimParticle* simp = fSimpBlock->Particle(0);
    p_sim  = simp->StartMom()->P();
    

  }

  Hist->fChi2XYNDof    ->Fill(p_sim,Helix->Chi2XY());
  Hist->fChi2PhiZNDof  ->Fill(p_sim,Helix->Chi2ZPhi());

  Hist->fLambdaVsP     ->Fill(p_sim,-lambda);
  Hist->fRadiusVsP     ->Fill(p_sim,radius);    
  Hist->fNRotVsP       ->Fill(p_sim,nRot);

  float frac = (float) nhits/( fStrawDataBlock->NHits() );
  Hist->fFracSHVsP     ->Fill(p_sim, frac);
  Hist->fPrecoVsP->Fill(p_sim,p);

}

//-----------------------------------------------------------------------------
void TValidationModule2::FillClusterHistograms(ClusterHist_t* Hist, TStnCluster* Cluster) {
  int   row, col;
  float  x, y, z, r;

  row = Cluster->Ix1();
  col = Cluster->Ix2();

  x   = Cluster->fX+3904.;
  y   = Cluster->fY;
  z   = Cluster->fZ;
  r   = sqrt(x*x+y*y);

  if ((row < 0) || (row > 9999)) row = -9999;
  if ((col < 0) || (col > 9999)) col = -9999;

  Hist->fVaneID->Fill(Cluster->DiskID());
  Hist->fEnergy->Fill(Cluster->Energy());
  Hist->fT0->Fill(Cluster->Time());
  Hist->fRow->Fill(row);
  Hist->fCol->Fill(col);
  Hist->fX->Fill(x);
  Hist->fY->Fill(y);
  Hist->fZ->Fill(z);
  Hist->fR->Fill(r);

  Hist->fYMean->Fill(Cluster->fYMean);
  Hist->fZMean->Fill(Cluster->fZMean);
  Hist->fSigY->Fill(Cluster->fSigY);
  Hist->fSigZ->Fill(Cluster->fSigZ);
  Hist->fSigR->Fill(Cluster->fSigR);
  Hist->fNCr0->Fill(Cluster->fNCrystals);
  Hist->fNCr1->Fill(Cluster->fNCr1);
  Hist->fFrE1->Fill(Cluster->fFrE1);
  Hist->fFrE2->Fill(Cluster->fFrE2);
  Hist->fSigE1->Fill(Cluster->fSigE1);
  Hist->fSigE2->Fill(Cluster->fSigE2);
}

//-----------------------------------------------------------------------------
void TValidationModule2::FillGenpHistograms(GenpHist_t* Hist, TGenParticle* Genp) {
  int    gen_id;
  float  p, cos_th, z0, t0, r0, x0, y0;

  TLorentzVector mom, v;

  Genp->Momentum(mom);
  //  Genp->ProductionVertex(v);

  p      = mom.P();
  cos_th = mom.CosTheta();

  x0     = Genp->Vx()+3904.;
  y0     = Genp->Vy();

  z0     = Genp->Vz();
  t0     = Genp->T();
  r0     = sqrt(x0*x0+y0*y0);
  gen_id = Genp->GetStatusCode();

  Hist->fPdgCode[0]->Fill(Genp->GetPdgCode());
  Hist->fPdgCode[1]->Fill(Genp->GetPdgCode());
  Hist->fGenID->Fill(gen_id);
  Hist->fZ0->Fill(z0);
  Hist->fT0->Fill(t0);
  Hist->fR0->Fill(r0);
  Hist->fP->Fill(p);
  Hist->fCosTh->Fill(cos_th);
 
  if(Genp->GetPdgCode() == 2212){ //If it is proton
    if (Genp->GetStatusCode()==28) Hist->fPpr->Fill(p);
    Hist->fPprCheck->Fill(p); 
  }

  if(Genp->GetPdgCode() == 1000010020){ //If it is deuton
    if (Genp->GetStatusCode()==28) Hist->fPdeu->Fill(p);    
    Hist->fPdeuCheck->Fill(p); 
  }
  Hist->fPrecoVsNSH->Fill(p,fNStrawHits); //%%
}

//-----------------------------------------------------------------------------
void TValidationModule2::FillSimpHistograms(SimpHist_t* Hist, TSimParticle* Simp) {

  Hist->fPdgCode->Fill(Simp->fPdgCode);
  Hist->fMomTargetEnd->Fill(Simp->fMomTargetEnd);
  Hist->fMomTrackerFront->Fill(Simp->fMomTrackerFront);
  Hist->fNStrawHits->Fill(Simp->fNStrawHits);

  //Trying to use the VirtualDetectors
  if(fStepPointMCBlock) {
    //std::cout<<"there seems to be something"<<std::endl;
    //std::cout<<fStepPointMCBlock->NStepPoints()<<std::endl;
    TStepPointMC* fSpmc;
    Float_t p;
    Float_t cos_th;
    
    for(int j = 0; j<fStepPointMCBlock->NStepPoints(); j++){
      fSpmc = fStepPointMCBlock->StepPointMC(j);

      if(fSpmc->VolumeID()==15){ //15 end tracker (DataProducts/inc/VirtualDetectorId.hh)
	p = fSpmc->fMom.Mag();
	Hist->fPvdEndT->Fill(p);
      }


      if(fSpmc->VolumeID()==10){ //10,13,14 (DataProducts/inc/VirtualDetectorId.hh) fSpmc->VolumeID()==13 || fSpmc->VolumeID()==14  
	//fSpmc->Print();
	p = fSpmc->fMom.Mag();
	cos_th = fSpmc->fMom.Z()/p;
	Hist->fPvd->Fill(p);
	Hist->fCosTh->Fill(cos_th);
	Hist->fPvdVsNSH->Fill(p ,fNStrawHits);
	Hist->fPvdVsCosTh->Fill(p,cos_th);
	Hist->fCosThVsNSH->Fill(cos_th,fNStrawHits);
	
	// looking at the difference with genp
	std::cout<<fGenpBlock->NParticles()<<std::endl;

	if(fGenpBlock->NParticles()==1){ 
	  TLorentzVector momgenp;
	  fGenpBlock->Particle(0)->Momentum(momgenp);
	  Hist->fPvdVsPgen->Fill(momgenp.P(),p);// %%%
	}
      }
    }
  }


}

//-----------------------------------------------------------------------------
// for DIO : ultimately, one would need to renormalize the distribution
//-----------------------------------------------------------------------------
void TValidationModule2::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track) {

  TLorentzVector  mom;
  double          r;
  int             itrk; 
  TrackPar_t*     tp;
					// pointer to local track parameters
  itrk = Track->Number();
  tp   = fTrackPar+itrk;

  Hist->fP[0]->Fill (Track->fP);
  Hist->fP[1]->Fill (Track->fP);
  Hist->fP[2]->Fill (Track->fP);
  Hist->fP0->  Fill (Track->fP0);
  Hist->fP2->  Fill (Track->fP2);

  Hist->fPDio->Fill(Track->fP,tp->fDioWt);

  Hist->fFitMomErr->Fill(Track->fFitMomErr);

  Hist->fPt    ->Fill(Track->fPt    );
  Hist->fPFront->Fill(Track->fPFront);
  Hist->fPStOut->Fill(Track->fPStOut);
					// dp: Tracker-only resolution

  Hist->fDpFront ->Fill(tp->fDpF);
  Hist->fDpFront0->Fill(tp->fDp0);
  Hist->fDpFront2->Fill(tp->fDp2);
  Hist->fDpFSt   ->Fill(tp->fDpFSt);
  Hist->fDpFVsZ1 ->Fill(Track->fZ1,tp->fDpF);

  Hist->fCosTh->Fill(Track->Momentum()->CosTheta());
  Hist->fChi2->Fill (Track->fChi2);
  Hist->fNDof->Fill(Track->NActive()-5.);
  Hist->fChi2Dof->Fill(Track->fChi2/(Track->NActive()-5.));
  Hist->fNActive->Fill(Track->NActive());
  Hist->fT0->Fill(Track->fT0);
  Hist->fT0Err->Fill(Track->fT0Err);
  //  printf("TValidationModule2::FillTrackHistograms: track charge is not defined yet\n");
  Hist->fQ->Fill(-1);
  Hist->fFitCons[0]->Fill(Track->fFitCons);
  Hist->fFitCons[1]->Fill(Track->fFitCons);

  Hist->fD0->Fill(Track->fD0);
  Hist->fZ0->Fill(Track->fZ0);
  Hist->fTanDip->Fill(Track->fTanDip);
  Hist->fAlgMask->Fill(Track->AlgMask());

  //  int nh, nst_with_nh[10];
					// 2014-04-29: currently not saved

  //  for (int i=0; i<10; i++) nst_with_nh[i] = 0;

//   for (int i=0; i<40; i++) {
//     Hist->fNHVsStation->Fill(i,Track->fNHPerStation[i]);
//     nh = Track->fNHPerStation[i];
//     if ((nh >= 0) && (nh < 10)) {
//       nst_with_nh[nh] += 1;
//     }
//     else {
//       printf(">>> ERROR : nh = %20i, IGNORE \n",nh);
//     }
//  }

//   for (int i=0; i<10; i++) {
//     Hist->fNHVsNSt->Fill(i,nst_with_nh[i]);
//   }
  //-----------------------------------------------------------------------------
  // track-cluster matching part: 
  // - for residuals, determine intersection with the most energetic cluster
  // - for track -only parameters use intersection with lowest trajectory length
  //-----------------------------------------------------------------------------
  TStnTrack::InterData_t*    vt = Track->fVMinS;  // track-only
  //  TStnTrack::InterData_t*    vr = Track->fVMaxEp; // residuals

  if (vt) {
    Hist->fVaneID->Fill(vt->fID  );
    Hist->fXTrk->Fill  (vt->fXTrk);
    Hist->fYTrk->Fill  (vt->fYTrk);

    r = sqrt(vt->fXTrk*vt->fXTrk+vt->fYTrk*vt->fYTrk);
    Hist->fRTrk->Fill  (r);

    Hist->fZTrk->Fill  (vt->fZTrk);
  }
  else {
//-----------------------------------------------------------------------------
// fill histograms with numbers easy to recognize as dummy
//-----------------------------------------------------------------------------
    Hist->fVaneID->Fill(-1.);
    Hist->fXTrk->Fill  (999.);
    Hist->fYTrk->Fill  (999.);
    Hist->fRTrk->Fill  (999.);
    Hist->fZTrk->Fill  (-1. );
  }

//-----------------------------------------------------------------------------
// there is an inconsistency in the SIMP block filling - in Mu2e offline 
// the particle momentumis is kept in MeV/c, while the PDG mass  -in GeV/c^2..
// thus the energy is screwed up... kludge around
// assign muon mass
//-----------------------------------------------------------------------------
  double ekin(-1.);
  if (fSimp) {
    double p, m;
    //    p    = fSimp->fStartMom.P();
    p = Track->fP;
    m    = 105.658; // in MeV
    ekin = sqrt(p*p+m*m)-m;
  }

  Hist->fECl->Fill(tp->fEcl);
  Hist->fEClEKin->Fill(tp->fEcl/ekin);
  Hist->fEp->Fill(tp->fEp);
  Hist->fEpVsPath->Fill(tp->fPath,tp->fEp);

  Hist->fDx->Fill(tp->fDx);
  Hist->fDy->Fill(tp->fDy);
  Hist->fDz->Fill(tp->fDz);

  Hist->fDt->Fill(tp->fDt);
  Hist->fChi2Match->Fill(tp->fChi2Match);

  Hist->fDu->Fill    (tp->fDu);
  Hist->fDv->Fill    (tp->fDv);
  Hist->fDvVsDu->Fill(tp->fDu,tp->fDv);

  Hist->fPath->Fill(tp->fPath);
  Hist->fDuVsPath->Fill(tp->fPath,tp->fDu);

  //  double cu[4] = { -60.3049, -0.749111,    0.00522242,  -7.52018e-06};
  double cu[4] = { -59.5174, -0.541226, 0.00414309, -5.84989e-06 };
  double cv[4] = {  6.44161, 0.0722353, -0.000653084, 1.14054e-06};

  double x = tp->fPath;
  double corr_u = cu[0]+cu[1]*x+cu[2]*x*x+cu[3]*x*x*x;
  double corr_v = cv[0]+cv[1]*x+cv[2]*x*x+cv[3]*x*x*x;

  //  double duc = tp->fDu-0.34*(tp->fPath-350.);
  double duc = tp->fDu-corr_u;
  double dvc = tp->fDv-corr_v;
  
  Hist->fDucVsPath->Fill(tp->fPath,duc);
  Hist->fDvcVsPath->Fill(tp->fPath,dvc);

  Hist->fDvVsPath->Fill(tp->fPath,tp->fDv);
  Hist->fDtVsPath->Fill(tp->fPath,tp->fDt);

  Hist->fDuVsTDip->Fill(Track->fTanDip,tp->fDu);
  Hist->fDvVsTDip->Fill(Track->fTanDip,tp->fDv);

  Hist->fZ1->Fill(Track->fZ1);

  int ncl = Track->NClusters();
  Hist->fNClusters->Fill(ncl);

  Hist->fRSlope->Fill(Track->RSlope());
  Hist->fXSlope->Fill(Track->XSlope());

  double llhr_dedx, llhr_xs, llhr_cal, llhr_trk, llhr;

  Hist->fEleLogLHCal->Fill(Track->EleLogLHCal());
  Hist->fMuoLogLHCal->Fill(Track->MuoLogLHCal());

  llhr_cal = Track->LogLHRCal();
  Hist->fLogLHRCal->Fill(llhr_cal);

  llhr_dedx = Track->LogLHRDeDx();
  llhr_xs   = Track->LogLHRXs();
  llhr_trk  = Track->LogLHRTrk();
  llhr      = llhr_cal+llhr_trk;

  Hist->fEpVsDt->Fill(tp->fDt,tp->fEp);
  Hist->fLogLHRDeDx->Fill(llhr_dedx);
  Hist->fLogLHRXs->Fill(llhr_xs);
  Hist->fLogLHRTrk->Fill(llhr_trk);
  Hist->fLogLHR->Fill(llhr);

  Hist->fPdgCode->Fill(Track->fPdgCode);
  Hist->fFrGH->Fill(Track->fNGoodMcHits/(Track->NActive()+1.e-5));

  Hist->fNEPlVsNHPl->Fill(tp->fNEPl,tp->fNHPl);
  Hist->fNDPlVsNHPl->Fill(tp->fNDPl,tp->fNHPl);
  Hist->fChi2dVsNDPl->Fill(tp->fNDPl,Track->Chi2Dof());
  Hist->fDpFVsNDPl  ->Fill(tp->fNDPl,tp->fDpF);

  float        fre1(-1), fre2(-1);
  int          icl;
  TStnCluster* cl;

  if (Track->fVMinS) {
    icl = Track->fVMinS->fClusterIndex;
    if (icl >= 0) {
      cl = fClusterBlock->Cluster(icl);
      fre1 = cl->fFrE1;
      fre2 = cl->fFrE2;
    }
  }

  Hist->fFrE1->Fill(fre1);
  Hist->fFrE2->Fill(fre2);
   
//<-------------------------------------------------------------------------------------------------------------------------||
//<---------------------------$$$-Checking the proton reconstruction-$$$----------------------------------------------------||
//<--------------------------------------------------------------------------------------------------------bvitali-Oct-2019-||
/*
This are the hist that bvitali used to improve and check the efficency of proton reconstruction
First part is for every type of event, then it's split between SingleParticle and NoPrimary.
You enter both if it is 
*/  
/*
  if (Track->fPdgCode == 1000010020) {  
    //import the MC particle associated to the track
    TSimParticle* simp = fSimpBlock->FindParticle(Track->fPartID);
    float p_genp  = simp->StartMom()->P();
    
    //True momentum of the reconstructed particle
    Hist->fDeuPreco->Fill(p_genp);

    //True momentum VS reconstructed momentum
    Hist->fDeuPrecoVsP->Fill(p_genp,Track->fP);    
  }
*/  


/*

      //True momentum of the reconstructed particle
    Hist->fPrecoprAll->Fill(p_genp);

    //True momentum VS reconstructed momentum
    Hist->fPrecoVsPAll->Fill(p_genp,Track->fP);    

    //Chi2/dof vs momentum. Is the dip at 250 linked to chi2?
    Hist->fChi2dVsPAll->Fill(p_genp,Track->fChi2/(Track->NActive()-5.));
*/
    

  //Check if the particle is a proton
  if (Track->fPdgCode == 2212|| Track->fPdgCode == 1000010020) {       // WHAT ABOUT DEUTONS    || Track->fPdgCode == 1000010020

      //import the MC particle associated to the track
  TSimParticle* simp = fSimpBlock->FindParticle(Track->fPartID);
  float p_genp  = simp->StartMom()->P();

  //True momentum of the reconstructed particle
    Hist->fPrecopr->Fill(p_genp);

    //True momentum VS reconstructed momentum
    Hist->fPrecoVsP->Fill(p_genp,Track->fP);    

    //Chi2/dof vs momentum. Is the dip at 250 linked to chi2?
    Hist->fChi2dVsP->Fill(p_genp,Track->fChi2/(Track->NActive()-5.));

    TStrawHitData*  sh = NULL;
    TTrackStrawHitData*  tsh = NULL;
    double tan = simp->StartMom()->Pz()/simp->StartMom()->P(); ////////////////COSENO 
    
    //How are p and the angle linked?
    Hist->fPrecoVsTanTh->Fill(p_genp,tan); 


    //SINGLE PROTON EVENT
    if(GetHeaderBlock()->InstLum() == -1 && fNStrawHits > 0){
      
      //Number of hits and fraction of 'active hits' Vs True momentum (with SHBlock and not TrackStrawHitBlock)  
      double fractionActive = (double)Track->NActive()/(double)fNStrawHits;
      Hist->fPrecoVsNSH->Fill(p_genp,fNStrawHits);   
      Hist->fFracSHVsPreco->Fill(p_genp,fractionActive);  

      int nTrackHits = -1; 
      bool used = false;
      int nonActive = 0;
      
      //variables to check initial and final point of the track
      int minTrkStation = 30;
      int maxTrkStation = -2;
      
      //loop on all the SH in the SHBlock (we are in a 'Single proton event')
      for (int k=0; k<fNStrawHits; k++ ) {
	sh = fStrawDataBlock->Hit(k);
	Hist->fPrecoVsSHE->Fill(p_genp,sh->Energy());
	Hist->fTanThVsSHE->Fill(tan,sh->Energy()); 
	
	used = false;
	nonActive=0;

	// If there is something in TrackSHBlock
	if(fTrackSHBlock->NTracks()!=0){
	  //there should be just 1 track, the number 0
	  nTrackHits = fTrackSHBlock->NTrackHits(0);
	  //if(nTrackHits != Track->NHits() && k == 0) printf("SHBlock %i : Track %i \n",nTrackHits, Track->NHits());
	  if(nTrackHits>fNStrawHits ||Track->NHits()> fNStrawHits) printf("Something is wrong\n");

	  //loop over the TrackSH: check if the sh is used and take the max Station of the sh
	  //is done like this because there is no 'Station()' in TrackSH
	  for(int i = 0; i <  nTrackHits; i++){
	    nonActive = 0;
	    tsh = fTrackSHBlock->Hit(0,i);
	    //printf("TrackStrawHitData Index = %i \n", tsh->Index());
	    if(tsh->Active()!=1) nonActive = nonActive+1;
	    if(tsh->Index()==sh->Index()) {
	      used=true;
	      Hist->fUsedSHVsPreco->Fill(p_genp,tsh->Index());
	      
	      if(sh->Station()>maxTrkStation) maxTrkStation=sh->Station();
	      if(sh->Station()<minTrkStation) minTrkStation=sh->Station();
	    }
	  }
	  if(!used) Hist->fUnusedSHVsPreco->Fill(p_genp,sh->Index()); 
	}
	else printf("No TrackSHBlock \n");
      }
      Hist->fLastTrkZVsPreco->Fill(p_genp,maxTrkStation);
      Hist->fFirstTrkZVsPreco->Fill(p_genp,minTrkStation);
    }
   
    //NON SIGNLE PARTICLE EVENTS
    else if(GetHeaderBlock()->InstLum() !=-1 && fTrackSHBlock->NTracks()!=0 ){
      //Number of hits and fraction of 'active hits' Vs True momentum  
      double fractionActive = (double)Track->NActive()/(double)Track->NHits();
      Hist->fPrecoVsNSH->Fill(p_genp,Track->NHits());   
      Hist->fFracSHVsPreco->Fill(p_genp,fractionActive);  
      
      int which_track = Track->Number();
      for(int k=0;k<fTrackSHBlock->NTrackHits(which_track);k++){
	tsh = fTrackSHBlock->Hit(which_track,k);
	Hist->fPrecoVsSHE->Fill(p_genp,tsh->Energy());
	Hist->fTanThVsSHE->Fill(tan,tsh->Energy()); 
      }
    }
  }
  
  else {printf("non proton or deuton track, look in pdio, at 120 \n"); Hist->fPDio->Fill(120);} //Random just to see if there are non proton/deuton tracks
//<--------------------------------------------------------------------------------------------------------bvitali-Oct-2019-||
//<-------------------------------------------------------------------------------------------------------------------------||  
//<-------------------------------------------------------------------------------------------------------------------------||  
} 

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TValidationModule2::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("TimeClusterBlockDpP" ,"TStnTimeClusterBlock",&fTimeClusterBlock);
  RegisterDataBlock("TrackSeedBlockDpP","TStnTrackSeedBlock",&fTrackSeedBlock);
  RegisterDataBlock("HelixBlock"    ,"TStnHelixBlock"    ,&fHelixBlock); //^_^
  RegisterDataBlock(fTrackBlockName.Data() ,"TStnTrackBlock"    ,&fTrackBlock  );
  RegisterDataBlock("ClusterBlock"  ,"TStnClusterBlock"  ,&fClusterBlock);
  RegisterDataBlock("CalDataBlock"  ,"TCalDataBlock"     ,&fCalDataBlock);
  RegisterDataBlock("StrawDataBlock","TStrawDataBlock"   ,&fStrawDataBlock);
  RegisterDataBlock("GenpBlock"     ,"TGenpBlock"        ,&fGenpBlock);
  RegisterDataBlock("SimpBlock"     ,"TSimpBlock"        ,&fSimpBlock);
  RegisterDataBlock("TrackHitBlock", "TTrackStrawHitBlock"   ,&fTrackSHBlock); //<<-----------TrackHitBlock-------------||||||||||||||
  RegisterDataBlock("SpmcBlockVDet","TStepPointMCBlock",&fStepPointMCBlock);   //<<-----------StepPointMC---------------||||||||||||||


//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// initialize likelihood histograms
//-----------------------------------------------------------------------------
  // fTrackID->SetMinT0(fMinT0);
  fTrackID->SetMaxMomErr(0.3); // MDC2018
  fTrackID->SetMaxT0Err (2.);

//-----------------------------------------------------------------------------
// PID initialization: read the likelihood templates
//-----------------------------------------------------------------------------
  fLogLH->Init("v5_7_0");

  return 0;
}


//_____________________________________________________________________________
int TValidationModule2::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
void TValidationModule2::FillHistograms() {

  int          alg_mask, nsh, nactive;
  float        pfront, ce_pitch, reco_pitch, fcons, t0, sigt, sigp, p; 
  
  // Possible usefull variables not used for now.
  /*
  double       cos_th (-2.);
  cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();
  TStnCluster  *cl0;
  if (fNClusters > 0) {
    double  cl_e(-1.);
    int     disk_id(-1);
    cl0     = fClusterBlock->Cluster(0);
    cl_e    = cl0->Energy();
    disk_id = cl0->DiskID();
  }
  */

  //-----------------------------------------------------------------------------
  // event histograms
  //-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);

  if (fNTracks[0]> 0) FillEventHistograms(fHist.fEvent[1]);
  else                FillEventHistograms(fHist.fEvent[2]);

  if (fNGoodTracks > 0) {
    FillEventHistograms(fHist.fEvent[6]); 

    TLorentzVector    mom;
    
    fParticle->Momentum(mom);

    double p, cos_th;

    p      = mom.P();
    cos_th = mom.Pz()/p;

    if (GetDebugBit(31) && (cos_th > 0.8)) {
      GetHeaderBlock()->Print(Form(" bit:031 cos_th = %10.3f p = %10.3f ntrk = %5i",
				   cos_th, p, fNTracks[0]));
    }
  }

  TLorentzVector    mom;
  fParticle->Momentum(mom);
  p      = mom.P();

  if(fParticle->GetPdgCode() == 2212 || fParticle->GetPdgCode() == 1000010020) FillEventHistograms(fHist.fEvent[40]);       //////////////
 

//-----------------------------------------------------------------------------
// Dave's ladder for all tracks
// 1. N(straw hits) > 20
//-----------------------------------------------------------------------------
  if (fSimp) {
    nsh    = fSimp->NStrawHits();
    pfront = fSimp->fMomTrackerFront;
  }
  else {
    nsh    = -1;
    pfront = -1.e6;
  }
  
  if (nsh >= 20) {
    FillEventHistograms(fHist.fEvent[11]);
    if (pfront > 100.) {
      FillEventHistograms(fHist.fEvent[12]);
      
      ce_pitch = 0.7; // kludge
      if ((ce_pitch > 0.577) && (ce_pitch < 1.)) {
	FillEventHistograms(fHist.fEvent[13]);

	if (fNTracks[0] > 0) {
	  FillEventHistograms(fHist.fEvent[14]);

					// here we have a track reconstructed

	  TStnTrack* trk = fTrackBlock->Track(0);

	  fcons = trk->fFitCons;
	  t0    = trk->T0();
	  reco_pitch = trk->fTanDip;
	  sigp       = trk->fFitMomErr;
	  sigt       = trk->fT0Err;
	  nactive    = trk->fNActive;
	  p          = trk->fP;
					// fit quality
	  if ((nactive > 25) && (fcons > 2.e-3) && (sigp < 0.25) && (sigt < 1.0))  {
	    FillEventHistograms(fHist.fEvent[15]);
	    if (t0 > 700) {
	      FillEventHistograms(fHist.fEvent[16]);
	      if ((reco_pitch > 0.577) && (reco_pitch < 1.)) {
		FillEventHistograms(fHist.fEvent[17]);
		if (p > 103.5) {
		  FillEventHistograms(fHist.fEvent[18]);
		}
	      }
	    }
	  }

	  alg_mask = trk->AlgMask();

	  if ((alg_mask == 1) || (alg_mask == 3)) {
//-----------------------------------------------------------------------------
// track reconstructed with TrkPatRec 
//-----------------------------------------------------------------------------
	    FillEventHistograms(fHist.fEvent[24]);
	    if ((nactive > 25) && (fcons > 2.e-3) && (sigp < 0.25) && (sigt < 1.0))  {
	      FillEventHistograms(fHist.fEvent[25]);
	      if (t0 > 700) {
		FillEventHistograms(fHist.fEvent[26]);
		if ((reco_pitch > 0.577) && (reco_pitch < 1.)) {
		  FillEventHistograms(fHist.fEvent[27]);
		  if (p > 103.5) {
		    FillEventHistograms(fHist.fEvent[28]);
		  }
		}
	      }
	    }
	  }
	  else if (alg_mask == 2) {
//-----------------------------------------------------------------------------
// track reconstructed with CalPatRec, but not with TrkPatRec
//-----------------------------------------------------------------------------
//	    int x=0;
	  }
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// the same ladder for TrkPatRec tracks 
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Simp histograms
//-----------------------------------------------------------------------------
  if (fSimp) {
    FillSimpHistograms(fHist.fSimp[0],fSimp);
  }

//--------------------------------------------------------------------------------
// TimeCluster histograms
//--------------------------------------------------------------------------------
  TStnTimeCluster* tCluster;
  for (int i=0; i<fNTimeClusters[0]; ++i){
    
    tCluster = fTimeClusterBlock->TimeCluster(i);
    
    FillTimeClusterHistograms(fHist.fTimeCluster[0], tCluster);
    
    int         nhits    = tCluster->NHits();
    if (nhits >= 5 )    FillTimeClusterHistograms(fHist.fTimeCluster[1], tCluster);
    if (nhits >= 10 )   FillTimeClusterHistograms(fHist.fTimeCluster[2], tCluster);
    if (nhits >= 15 )   FillTimeClusterHistograms(fHist.fTimeCluster[3], tCluster);
    if (fNHelices[0]<1) FillTimeClusterHistograms(fHist.fTimeCluster[4], tCluster); //Only if there are NO helices
  }  

//--------------------------------------------------------------------------------
// trackseed histograms
//--------------------------------------------------------------------------------
  TStnTrackSeed* trkSeed;
  for (int i=0; i<fNTrackSeeds[0]; ++i){
    
    trkSeed = fTrackSeedBlock->TrackSeed(i);
    
    FillTrackSeedHistograms(fHist.fTrackSeed[0], trkSeed);
    
    int         nhits    = trkSeed->NHits();
    //    double      mm2MeV   = 3/10.;
    double      p        = trkSeed->P();
    //    double      radius   = trkSeed->Pt()/mm2MeV;
    //    double      tanDip   = trkSeed->TanDip();  
    
    double      chi2      = trkSeed->Chi2();
    //    double     fitCons = trkSeed->Fitcons();
    
   
    if (p > 150.) {
      FillTrackSeedHistograms(fHist.fTrackSeed[1], trkSeed);
    }

    if ( nhits>=10) {
      FillTrackSeedHistograms(fHist.fTrackSeed[2], trkSeed);
    }
    
    if ( (chi2 < 5) /*&& (chi2zphi < 4)*/ && (nhits>=10)){
      FillTrackSeedHistograms(fHist.fTrackSeed[3], trkSeed);
    }
  
  }


//--------------------------------------------------------------------------------
// helix histograms
//--------------------------------------------------------------------------------
  TStnHelix* helix;
  for (int i=0; i<fNHelices[0]; ++i){
    
    helix = fHelixBlock->Helix(i);
    
    FillHelixHistograms(fHist.fHelix[0], helix);
    
    int         nhits    = helix->NHits();
    double      radius   = helix->Radius();
    
    double      lambda   = helix->Lambda();
    double      tanDip   = lambda/radius;
    double      mm2MeV   = 3/10.;
    double      p        = radius*mm2MeV/std::cos( std::atan(tanDip));
    
    double      chi2xy   = helix->Chi2XY();
    double      chi2zphi = helix->Chi2ZPhi();
    
   
    if (p >= 250.) {
      FillHelixHistograms(fHist.fHelix[1], helix);
    }
    
    if (p < 250) {
      FillHelixHistograms(fHist.fHelix[2], helix);
    }
    
    if ( nhits>=10) {
      FillHelixHistograms(fHist.fHelix[3], helix);
    }
    
    if ( (chi2xy < 5) && (chi2zphi < 5) && (nhits>=10)){
      FillHelixHistograms(fHist.fHelix[4], helix);
    }

    if ( fNTracks[0] < 1 ) {
      FillHelixHistograms(fHist.fHelix[5], helix);  //Events with no tracks
    }

    if(helix->PDGMother1() == 2212){
      FillHelixHistograms(fHist.fHelix[6], helix);  //it is a proton
    }

    if(helix->PDGMother1() == 1000010020){
      FillHelixHistograms(fHist.fHelix[7], helix);  //it is a deuteron 1000010020
    }
  
  }
//-----------------------------------------------------------------------------
// track histograms
//-----------------------------------------------------------------------------
  TStnTrack*   trk;
  TrackPar_t*  tp;

  for (int i=0; i<fNTracks[0]; ++i ) {
    trk = fTrackBlock->Track(i);
    tp  = fTrackPar+i;

    FillTrackHistograms(fHist.fTrack[0],trk);

    if (trk->fIDWord == 0) {
					// track passes selection "C" 
      FillTrackHistograms(fHist.fTrack[1],trk);
    }
//-----------------------------------------------------------------------------
// TRK 6 : events with a track and a cluster E/P > 0.4
//-----------------------------------------------------------------------------
    if (trk->Ep() > 0.4) {
      FillTrackHistograms(fHist.fTrack[6],trk);
    }

//-----------------------------------------------------------------------------
// TRK 11: tracks with P > 150 MeV
//-----------------------------------------------------------------------------
    if (trk->P() > 150) FillTrackHistograms(fHist.fTrack[11],trk);
//-----------------------------------------------------------------------------
// TRK 14: tracks with fcon < 1e-2
//-----------------------------------------------------------------------------

    if (trk->fFitCons < 1.e-2) FillTrackHistograms(fHist.fTrack[14],trk);

    TStnTrack::InterData_t*    vt = trk->fVMinS;  // track-only
//-----------------------------------------------------------------------------
// TRK 15: tracks which have intersection with the 1st disk
// TRK 16: tracks which have intersection with the 2nd disk
// TRK 17: tracks which do not have intersections with the calorimeter
//-----------------------------------------------------------------------------
    if (vt) {
      if      (vt->fID == 0) {
	FillTrackHistograms(fHist.fTrack[15],trk);
      }
      else if (vt->fID == 1) {
	FillTrackHistograms(fHist.fTrack[16],trk);
      }
    }
    else {
      FillTrackHistograms(fHist.fTrack[17],trk);
    }


//-----------------------------------------------------------------------------
// TRK 20: tracks with >= 20 hits
// TRK 21: tracks with >= 20 hits and Chi2/Ndof < 3
//-----------------------------------------------------------------------------
    if (trk->fNActive >= 20) {
      FillTrackHistograms(fHist.fTrack[20],trk);
      if (trk->Chi2Dof() < 5) {
	FillTrackHistograms(fHist.fTrack[21],trk);
      }
    }
//-----------------------------------------------------------------------------
// TRK 28 : events with a "Set C" track and a cluster E/P > 1.1
//-----------------------------------------------------------------------------
    if (trk->fIDWord == 0) {
      if (trk->Ep() > 1.1) {
	FillTrackHistograms(fHist.fTrack[28],trk);

	if (GetDebugBit(28)) {
	  GetHeaderBlock()->Print(Form(" bit:028 LLHR(CAL) = %10.3f ep = %10.3f dt = %10.3f",
				       trk->LogLHRCal(), trk->Ep(), trk->Dt()));
	}
      }
    }
//-----------------------------------------------------------------------------
// TRK 29 : events with a "Set C" track, E/P>0 and P>100 : precursor for TRK 25
// in effect startign from 2014-06-17 10:31am
//-----------------------------------------------------------------------------
    if (trk->fIDWord == 0) {
      if ((tp->fEp > 0) && (trk->fP > 100.)) {
	FillTrackHistograms(fHist.fTrack[29],trk);
      }
    }
//-----------------------------------------------------------------------------
// TRK 30: tracks with >= 25 hits
// TRK 31: tracks with >= 25 hits and Chi2/Ndof < 3
//-----------------------------------------------------------------------------
    if (trk->fNActive >= 25) {
      FillTrackHistograms(fHist.fTrack[30],trk);
      if (trk->Chi2Dof() < 3) {
	FillTrackHistograms(fHist.fTrack[31],trk);
      }
    }
//-----------------------------------------------------------------------------
// TRK 63: tracks with P>=250                                                      <-----------------
// TRK 64: tracks with 250>P                                                       <-----------------
//-----------------------------------------------------------------------------
    if (trk->fP >= 250) {
      FillTrackHistograms(fHist.fTrack[63],trk);
      GetAna()->GetHeaderBlock()->Print(Form("63: momentum: %12.5f",trk->fP));
    }
    else if(trk->fP < 250){
      FillTrackHistograms(fHist.fTrack[64],trk);
      GetAna()->GetHeaderBlock()->Print(Form("64: momentum: %12.5f",trk->fP));
    }
//-----------------------------------------------------------------------------
// TRK 71: Just Protons                                                           <-----------------
// TRK 72: Just Deuterons                                                         <-----------------
//-----------------------------------------------------------------------------
    if (trk->fPdgCode == 2212) {
      FillTrackHistograms(fHist.fTrack[71],trk);
      GetAna()->GetHeaderBlock()->Print(Form("71: momentum: %12.5f",trk->fP));
    }
    else if(trk->fPdgCode == 1000010020){
      FillTrackHistograms(fHist.fTrack[72],trk);
      GetAna()->GetHeaderBlock()->Print(Form("72: momentum: %12.5f",trk->fP));
    }

//-----------------------------------------------------------------------------
// split tracks by the algorithm mask: 1 , 2 , or 3
//-----------------------------------------------------------------------------
    alg_mask = trk->AlgMask();
    if      (alg_mask == 1) {
//-----------------------------------------------------------------------------
// TrkPatRec-only tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[40],trk);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[41],trk);
	// print run/event numbers :
	if (GetDebugBit(6)) {
	  double ep = trk->Ep();
	  if ((ep > 0.8) && (ep < 1.1)) {
	    GetHeaderBlock()->Print(Form(" bit:006 trk_41: track E/P = %8.3f",trk->Ep()));
	  }
	}
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[42],trk);
      }
    }
    else if (alg_mask == 2) {
//-----------------------------------------------------------------------------
// CalPatRec-only tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[50],trk);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[51],trk);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[52],trk);
      }
    }
    else if (alg_mask == 3) {
//-----------------------------------------------------------------------------
// TrkPatRec+CalPatRec tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[60],trk);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[61],trk);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[62],trk);
      }
    }
}
//-----------------------------------------------------------------------------
// cluster histograms
//-----------------------------------------------------------------------------
  TStnCluster* cl;
  int id;
  for (int i=0; i<fNClusters; ++i ) {
    cl = fClusterBlock->Cluster(i);
    id = cl->DiskID();
    FillClusterHistograms(fHist.fCluster[0],cl);

    if (fNTracks[0]     >  0 ) FillClusterHistograms(fHist.fCluster[1],cl);
    if (fNGoodTracks    >  0 ) FillClusterHistograms(fHist.fCluster[2],cl);
    if (fNMatchedTracks >  0 ) FillClusterHistograms(fHist.fCluster[3],cl);
    if (cl->Energy()    > 10.) FillClusterHistograms(fHist.fCluster[4],cl);
    if (cl->Energy()    > 60.) FillClusterHistograms(fHist.fCluster[5],cl);

    if      (id == 0         ) FillClusterHistograms(fHist.fCluster[6],cl);
    else if (id == 1         ) FillClusterHistograms(fHist.fCluster[7],cl);
  }
//-----------------------------------------------------------------------------
// fill GENP histograms
// GEN_0: all particles
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  for (int i=0; i<fNGenp; i++) {
    genp = fGenpBlock->Particle(i);
    FillGenpHistograms(fHist.fGenp[0],genp);
    if (fNTracks[0] > 0) {
      FillGenpHistograms(fHist.fGenp[1],genp);
      if (fNGoodTracks > 0) {
	FillGenpHistograms(fHist.fGenp[2],genp);
      }
    }
    //if there is at least a SH  GetHeaderBlock()->InstLum()==-1
    if(fNStrawHits>0) FillGenpHistograms(fHist.fGenp[3],genp);
    if(fNStrawHits>4) FillGenpHistograms(fHist.fGenp[4],genp);
    if(fNStrawHits>19) FillGenpHistograms(fHist.fGenp[5],genp);

  }
}



//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TValidationModule2::Event(int ientry) {

  double                xs, p;
  TEmuLogLH::PidData_t  dat;
  TStnTrack*            track;
  int                   id_word;
  TLorentzVector        mom;

  fTrackBlock  ->GetEntry(ientry);
  fHelixBlock  ->GetEntry(ientry);
  fClusterBlock->GetEntry(ientry);
  fStrawDataBlock->GetEntry(ientry);
  fTrackSHBlock->GetEntry(ientry);
  fTimeClusterBlock->GetEntry(ientry);
  fCalDataBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fStepPointMCBlock->GetEntry(ientry); //%%%

//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fNGenp    = fGenpBlock->NParticles();

  TGenParticle* genp;
  int           pdg_code, generator_code;

  fParticle = NULL;
  for (int i=fNGenp-1; i>=0; i--) {
    genp           = fGenpBlock->Particle(i);
    pdg_code       = genp->GetPdgCode();
    generator_code = genp->GetStatusCode();
    if ((abs(pdg_code) == fPdgCode) && (generator_code == fGeneratorCode)) {
      fParticle = genp;
      break;
    }
  }
					// may want to revisit the definition of fSimp
  fSimp     = fSimpBlock->Particle(0);


  if (fParticle) {
    fParticle->Momentum(mom);
					// this is a kludge, to be removed at the next 
					// ntupling 
  //  fEleE     = fParticle->Energy();
    p         = mom.P();
  }
  else {
    p = 0.;
  }
  fEleE     = sqrt(p*p+0.511*0.511);


  fNTimeClusters[0] = fTimeClusterBlock->NTimeClusters();
  
  fNTracks[0] = fTrackBlock->NTracks();
  fNClusters  = fClusterBlock->NClusters();
  fNCalHits   = fCalDataBlock->NHits();
  fNStrawHits = fStrawDataBlock->NHits();
  fNHelices[0]= fHelixBlock->NHelices();
  fNTrackSeeds[0] = fTrackSeedBlock->NTrackSeeds();

  fNHyp       = -1;
  fBestHyp[0] = -1;
  fBestHyp[1] = -1;

  fNGoodTracks    = 0;
  fNMatchedTracks = 0;

  fNTracks[0] = fTrackBlock->NTracks();
  if (fNTracks[0] == 0) fTrack = 0;
  else                  fTrack = fTrackBlock->Track(0);

  int ntrk = fNTracks[0];

  TrackPar_t*   tp;

  for (int itrk=0; itrk<ntrk; itrk++) {
					// assume less 20 tracks
    tp             = fTrackPar+itrk;

    track          = fTrackBlock->Track(itrk);
    id_word        = fTrackID->IDWord(track);
    track->fIDWord = id_word;
    if (id_word == 0) {
      fNGoodTracks += 1;
      if ((track->fVMaxEp != NULL) && (fabs(track->fVMaxEp->fDt) < 2.5)) {
	fNMatchedTracks += 1;
      }
    }
//-----------------------------------------------------------------------------
// process hit masks
//-----------------------------------------------------------------------------
    int i1, i2, n1(0) ,n2(0), ndiff(0);
    int nbits = track->fHitMask.GetNBits();
    for (int i=0; i<nbits; i++) {
      i1 = track->HitMask()->GetBit(i);
      i2 = track->ExpectedHitMask()->GetBit(i);
      n1 += i1;
      n2 += i2;
      if (i1 != i2) ndiff += 1;
    }
//-----------------------------------------------------------------------------
// define additional parameters
//-----------------------------------------------------------------------------
    tp->fNHPl = n1;
    tp->fNEPl = n2;
    tp->fNDPl = ndiff;

    tp->fDpF   = track->fP     -track->fPFront;
    tp->fDp0   = track->fP0    -track->fPFront;
    tp->fDp2   = track->fP2    -track->fPFront;
    tp->fDpFSt = track->fPFront-track->fPStOut;

    if (fFillDioHist == 0) tp->fDioWt = 1.;
    else                   tp->fDioWt = TStntuple::DioWeightAl(fEleE);
//-----------------------------------------------------------------------------
// track residuals
//-----------------------------------------------------------------------------
    TStnTrack::InterData_t*  vr = track->fVMaxEp; 
    double    nx, ny;

    tp->fEcl       = -1.e6;
    tp->fEp        = -1.e6;

    tp->fDu        = -1.e6;
    tp->fDv        = -1.e6;
    tp->fDx        = -1.e6;
    tp->fDy        = -1.e6;
    tp->fDz        = -1.e6;
    tp->fDt        = -1.e6;

    tp->fChi2Match = -1.e6;
    tp->fPath      = -1.e6;

    if (vr) {
      tp->fEcl = vr->fEnergy;
      tp->fEp  = tp->fEcl/track->fP;

      tp->fDx  = vr->fDx;
      tp->fDy  = vr->fDy;
      tp->fDz  = vr->fDz;
//-----------------------------------------------------------------------------
// v4_2_4: correct by additional 0.22 ns - track propagation by 6 cm
//-----------------------------------------------------------------------------
      tp->fDt  = vr->fDt - 0.22; // - 1.;

      nx  = vr->fNxTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);
      ny  = vr->fNyTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);

      tp->fDu        = vr->fDx*nx+vr->fDy*ny;
      tp->fDv        = vr->fDx*ny-vr->fDy*nx;
      tp->fChi2Match = vr->fChi2Match;
      tp->fPath      = vr->fPath;
    }

    if ((tp->fEp > 0) && (track->fEp > 0) && (fabs(tp->fEp-track->fEp) > 1.e-6)) {
      GetHeaderBlock()->Print(Form(" TValidationModule2 ERROR: tp->fEp = %10.5f  track->fEp = %10.5f",tp->fEp,track->fEp));
    }
//-----------------------------------------------------------------------------
// PID likelihoods
//-----------------------------------------------------------------------------
    dat.fDt   = tp->fDt;
    dat.fEp   = tp->fEp;
    dat.fPath = tp->fPath;
      
    xs = track->XSlope();

    track->fEleLogLHCal = fLogLH->LogLHCal(&dat,11);
    track->fMuoLogLHCal = fLogLH->LogLHCal(&dat,13);

    double llhr_cal = track->fEleLogLHCal-track->fMuoLogLHCal;

    if (GetDebugBit(7)) {
      if ((id_word == 0) && (llhr_cal > 20)) {
	GetHeaderBlock()->Print(Form("bit:007: dt = %10.3f ep = %10.3f",track->Dt(),tp->fEp));
      }
    }

    if (GetDebugBit(8)) {
      if ((id_word == 0) && (llhr_cal < -20)) {
	GetHeaderBlock()->Print(Form("bit:008: p = %10.3f dt = %10.3f ep = %10.3f",
				     track->P(),track->Dt(),tp->fEp));
      }
    }

    track->fLogLHRXs    = fLogLH->LogLHRXs(xs);
  }

  fNClusters = fClusterBlock->NClusters();
  if (fNClusters == 0) fCluster = 0;
  else                 fCluster = fClusterBlock->Cluster(0);

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TValidationModule2::Debug() {

  TStnTrack* trk;
  TrackPar_t* tp;
  int ntrk = fTrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = fTrackBlock->Track(itrk);
    tp  = &fTrackPar[itrk];
//-----------------------------------------------------------------------------
// bit 3: Set C tracks with large DX : 70mm < |DX| < 90mm
//-----------------------------------------------------------------------------
    if (GetDebugBit(3) == 1) {
      if (trk->fIDWord == 0) {
	TStnTrack::InterData_t*    vr = trk->fVMaxEp; // residuals
	if ((vr && (fabs(vr->fDx) > 70) && (fabs(vr->fDx) < 90))) {
	  GetHeaderBlock()->Print(Form("large DX: %f",vr->fDx));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 4: tracks with DpF > 1MeV - positive tail...
//-----------------------------------------------------------------------------
    if (GetDebugBit(4) == 1) {
      if (tp->fDpF > 1.) {
	GetHeaderBlock()->Print(Form("pF pRec, fDpf = %10.3f  %10.3f  %10.3f",
				     trk->fPFront, trk->Momentum()->P(),tp->fDpF));
      }
    }
//-----------------------------------------------------------------------------
// bit 9: Set C tracks with DpF > 1MeV - positive tail...
//-----------------------------------------------------------------------------
    if (GetDebugBit(9) == 1) {
      double ep = trk->Ep();
      if (trk->fIDWord == 0) { 
	if (((ep > 0.42) && (ep < 0.46)) || ((ep > 0.35) && (ep < 0.39))) {
	  GetHeaderBlock()->Print(Form("bit:009 ep = %10.3f e = %10.3f p = %10.3f",
				       trk->fEp,trk->fEp*trk->fP,trk->fP));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 10: Set C tracks with Ecl > 80
//-----------------------------------------------------------------------------
    if (GetDebugBit(10) == 1) {
      double ecl = trk->ClusterE();
      if (trk->fIDWord == 0) { 
	if (ecl > 60) {
	  GetHeaderBlock()->Print(Form("bit:010 e = %10.3f p = %10.3f",
				       ecl,trk->fP));
	}
      }
    }
  }

//-----------------------------------------------------------------------------
// bit 5: events with N(tracks) > 1
//-----------------------------------------------------------------------------
   if (GetDebugBit(5) == 1) {
    int ntrk = fTrackBlock->NTracks();
    if (ntrk > 1) {
      GetHeaderBlock()->Print(Form("NTracks = %i5",ntrk));
    }
  }
}

//_____________________________________________________________________________
int TValidationModule2::EndJob() {
  int ntrk = fTrackBlock->NTracks();
  printf("----- NTracks = %i\n",ntrk);
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TValidationModule2::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}

