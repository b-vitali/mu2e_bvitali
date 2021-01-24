///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_ana_TValidationModule2_hh
#define Stntuple_ana_TValidationModule2_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"

#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TStnCrystal.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

#include "Stntuple/obj/TTrackStrawHitBlock.hh"
#include "Stntuple/obj/TTrackStrawHitData.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh" //


class TValidationModule2: public TStnModule {
public:

  struct TrackPar_t {
    int     fNHPl;
    int     fNEPl;
    int     fNDPl;
    float   fDpF ;    // tracker-only resolution
    float   fDp0 ;
    float   fDp2 ;
    float   fDpFSt;
    double  fDioWt;

    double  fEcl;
    double  fEp;
    double  fDx;
    double  fDy;
    double  fDz;
    double  fDt;
    double  fDu;			// rotated residuals
    double  fDv;
    double  fChi2Match;
    double  fPath;
  };
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct ClusterHist_t {
    TH1F*    fVaneID;
    TH1F*    fEnergy;
    TH1F*    fT0;
    TH1F*    fRow;
    TH1F*    fCol;
    TH1F*    fX;
    TH1F*    fY;
    TH1F*    fZ;
    TH1F*    fR;
    TH1F*    fNCr0;			// all clustered
    TH1F*    fNCr1;			// above 1MeV
    TH1F*    fYMean;
    TH1F*    fZMean;
    TH1F*    fSigY;
    TH1F*    fSigZ;
    TH1F*    fSigR;
    TH1F*    fFrE1;
    TH1F*    fFrE2;
    TH1F*    fSigE1;
    TH1F*    fSigE2;
  };

  struct EventHist_t {
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;
    TH1F*    fEleMom;
    TH1D*    fDioMom;
    TH1F*    fEleCosTh;
    TH1F*    fNClusters;
    TH1F*    fNTracks;

    TH1F*    fNTracksP;                 // bastiano Proton reco
    TH1F*    fNTracksD;                 // bastiano Deuteron reco
    TH1F*    fNTracksO;                 // bastiano others reco

    TH1F*    fNTracksCut;
    TH1F*    fNTracksPCut;                 // bastiano Proton reco
    TH1F*    fNTracksDCut;                 // bastiano Deuteron reco
    TH1F*    fNTracksOCut;                 // bastiano others reco

    TH2F*    fEffVsLum;                // number of trks and luminosity


    TH1F*    fNGoodTracks;
    TH1F*    fNStrawHits[2];
    TH1F*    fNGoodSH;
    TH1F*    fDtClT;
    TH1F*    fEMax;			// energy of the first reco cluster
    TH1F*    fDtClS;
    TH1F*    fSHTime;
    TH1F*    fNHyp;
    TH1F*    fBestHyp[2];		// [0]: by chi2, [1]: by fit consistency
    TH1F*    fNTimeClusters;            // bastiano number of timeclusters
    TH2F*    fNTimeClustersVsMom;       // bastiano number of timeclusters
    
    TH2F*    fTCTDiscanceVsMom;          // bastiano number of timeclusters
    TH2F*    fTCZDiscanceVsMom;

    TH1F*    fNGenp;                    // N(particles in GENP block)
    TH1F*    fQSH;                      // SH charge
    TH1F*    fQSH_p;                    // bastiano SH charge from protons
    TH1F*    fQSH_d;                    // bastiano SH charge from deuteron
    TH1F*    fQSH_e;                    // JJ SH charge from electro
    TH2F*    fNSHVsPreco;               // bastiano proton eff
    
    TH2F*    fLumVsNSH;
    TH2F*    fLumVsNTrk;                // number of trks and luminosity

    TH2F*    fLumVsNTrkCut;             // bastiano proton eff
    TH1F*    fLum;                      // luminosity
    TH1F*    fLumRel;                   // relative luminosity


    TH2F*    fLastZVsPreco;             // bastiano proton eff
    TH2F*    fFirstZVsPreco;            // bastiano proton eff
    TH2F*    fZVsPreco;                 // bastiano proton eff
    TH2F*    fStdTimeVsMom;             // bastiano proton eff
    TH2F*    fWidthTimeVsMom;
  };

  struct TimeClusterHist_t {
    TH1F*    fNHits;
    TH1F*    fNComboHits;
    TH1F*    fT0;
    TH1F*    fClusterEnergy;
    TH2F*    fNComboHitsVsP;
    TH2F*    fFracSHVsP;
  
  };


  struct TrackSeedHist_t {
    TH1F*    fNHits;	 
    TH1F*    fClusterTime;
    TH1F*    fClusterEnergy;
    TH1F*    fRadius;
    TH1F*    fMom;
    TH1F*    fPt;
    TH1F*    fTanDip;   
    TH1F*    fChi2;
    TH1F*    fFitCons;
    TH1F*    fD0;
  };


  struct HelixHist_t {
    TH1F*    fNHits;	 
    TH1F*    fClusterTime;
    TH1F*    fClusterEnergy;
    TH1F*    fRadius;    // fabs(1/omega)
    TH1F*    fMom;
    TH1F*    fPt;
    TH1F*    fLambda;   //dz/dphi
    TH1F*    fT0;
    TH1F*    fT0Err;
    TH1F*    fD0;
    TH1F*    fAlg;
    TH2F*    fChi2XYNDof;
    TH2F*    fChi2PhiZNDof;
    TH2F*    fLambdaVsP;          // bvitali proton_reco
    TH2F*    fRadiusVsP;          // bvitali proton_reco
    TH2F*    fNRotVsP;            // bvitali proton_reco
    TH2F*    fFracSHVsP;          // bvitali proton_reco
    TH2F*    fPrecoVsP;           // bvitali proton_reco
  };

  struct TrackHist_t {
    TH1F*    fP[3];			// total momentum, 3 hists with different binning
    TH1F*    fP0;
    TH1F*    fP2;
    TH1F*    fPt;
    TH1D*    fPDio;                     // momentum dist weighted with the DIO weight
    TH1F*    fFitMomErr;
    TH1F*    fPFront;
    TH1F*    fDpFront;
    TH1F*    fDpFront0;
    TH1F*    fDpFront2;
    TH2F*    fDpFVsZ1;
    TH1F*    fPStOut;
    TH1F*    fDpFSt;			// P(TT_Hollow) - P(ST_Out)
    TH1F*    fCosTh;
    TH1F*    fChi2;
    TH1F*    fNDof;
    TH1F*    fChi2Dof;
    TH1F*    fNActive;
    TH1F*    fT0;
    TH1F*    fT0Err;
    TH1F*    fQ;
    TH1F*    fFitCons[2];		// fit consistency (0 to 1)
    TH1F*    fD0;
    TH1F*    fZ0;
    TH1F*    fTanDip;
    TH1F*    fResid;
    TH1F*    fAlgMask;
					// matching histograms
    TH1F*    fNClusters;
    TH1F*    fVaneID;
    TH1F*    fXCal;
    TH1F*    fYCal;
    TH1F*    fZCal;
    TH1F*    fXTrk;
    TH1F*    fYTrk;
    TH1F*    fZTrk;
    TH1F*    fRTrk;
    TH1F*    fDt;			// track-cluster residuals
    TH1F*    fChi2Match;
    TH1F*    fDt_eMinus;
    TH1F*    fDt_ePlus;
    TH1F*    fDt_muMinus;
    TH1F*    fDt_muPlus;
    TH1F*    fDx;
    TH1F*    fDy;
    TH1F*    fDz;
    TH1F*    fDu;
    TH1F*    fDv;
    TH2F*    fDvVsDu;
    TH1F*    fPath;
    TH2F*    fDuVsPath;
    TH2F*    fDucVsPath;
    TH2F*    fDvVsPath;
    TH2F*    fDvcVsPath;
    TH2F*    fDtVsPath;
    TH2F*    fDuVsTDip;
    TH2F*    fDvVsTDip;
    TH1F*    fZ1;
    TH1F*    fECl;
    TH1F*    fEClEKin;
    TH1F*    fEp;
    TH2F*    fEpVsPath;
    TH1F*    fEp_eMinus;
    TH1F*    fEp_ePlus;
    TH1F*    fEp_muMinus;
    TH1F*    fEp_muPlus;
    TH2F*    fNHVsStation;
    TH2F*    fNHVsNSt;

    TH1F*    fRSlope;
    TH1F*    fXSlope;
					// likelihoods
    TH2F*    fEpVsDt;
    TH1F*    fEleLogLHCal;
    TH1F*    fMuoLogLHCal;
    TH1F*    fLogLHRCal;
    TH1F*    fLogLHRDeDx;
    TH1F*    fLogLHRXs;
    TH1F*    fLogLHRTrk;
    TH1F*    fLogLHR;
					// MC truth
    TH1F*    fPdgCode;	                // PDG code of the particle produced most hits
    TH1F*    fFrGH;			// fraction of hits produced by the particle

    TH2F*    fNEPlVsNHPl;
    TH2F*    fNDPlVsNHPl;
    TH2F*    fChi2dVsNDPl;
    TH2F*    fDpFVsNDPl;

    TH1F*    fFrE1;
    TH1F*    fFrE2;
    /*
    TH1F*    fDeuPreco;                 // bastiano deuton eff     
    TH2F*    fDeuPrecoVsP;              // bastiano deuton eff
    */

    //TH1F*    fPrecoprAll;                  // bastiano proton eff	          
    //TH2F*    fPrecoVsPAll;                 // bastiano proton eff
    //TH2F*    fChi2dVsPAll;                 // bastiano proton eff

    TH1F*    fPrecopr;                  // bastiano proton eff	          
    TH2F*    fPrecoVsP;                 // bastiano proton eff
    TH2F*    fChi2dVsP;                 // bastiano proton eff
    TH2F*    fPrecoVsTanTh;             // bastiano proton eff
    TH2F*    fPrecoVsSHE;               // bastiano proton eff
    TH2F*    fTanThVsSHE;               // bastiano proton eff
    TH2F*    fPrecoVsNSH;               // bastiano proton eff
    TH2F*    fFracSHVsPreco;            // bastiano proton eff
    TH2F*    fLastTrkZVsPreco;          // bastiano proton eff
    TH2F*    fFirstTrkZVsPreco;         // bastiano proton eff
    TH2F*    fUnusedSHVsPreco;          // bastiano proton eff
    TH2F*    fUsedSHVsPreco;            // bastiano proton eff
  };

  struct GenpHist_t {
    TH1F*    fPdgCode[2];		// same distribution in different scale
    TH1F*    fGenID;			// 
    TH1F*    fZ0;			// 
    TH1F*    fT0;			// 
    TH1F*    fR0;			// 
    TH1F*    fP;			// 
    TH1F*    fCosTh;			// 
    
    TH1F*    fPpr;                      //bastiano proton eff
    TH1F*    fPprCheck;                 //bastiano proton eff
    TH1F*    fPdeu;                     //bastiano proton eff
    TH1F*    fPdeuCheck;                //bastiano proton eff
    TH2F*    fPrecoVsNSH;               //bastiano proton eff

  };
					// histograms for the simulated CE
  struct SimpHist_t {
    TH1F*    fPdgCode;
    TH1F*    fMomTargetEnd;
    TH1F*    fMomTrackerFront;
    TH1F*    fNStrawHits;

    TH1F*    fPvd;                      //bastiano proton eff
    TH1F*    fCosTh;                    //bastiano proton eff
    TH2F*    fPvdVsNSH;                 //bastiano proton eff
    TH2F*    fPvdVsCosTh;               //bastiano proton eff
    TH2F*    fCosThVsNSH;               //bastiano proton eff
    TH2F*    fPvdVsPgen;                //bastiano proton eff
    TH1F*    fPvdEndT;                  //bastiano proton eff

  };

//-----------------------------------------------------------------------------
//  fTrackHist[ 0]: all the tracks
//  fTrackHist[ 1]: all the tracks Pt > PtMin and within the fiducial
//  fTrackHist[ 2]: [1]+ matched to MC
//  fTrackHist[ 3]: [1]+ not matched to MC
//  fTrackHist[ 4]: [2]+ inside  the jet
//  fTrackHist[ 5]: [3]+ inside  the jet
//  fTrackHist[ 6]: [2]+ outside the jet
//  fTrackHist[ 7]: [3]+ outside the jet
//  fTrackHist[ 8]:
//  fTrackHist[ 9]:
//  fTrackHist[10]:
//  fTrackHist[11]: tracks with pt.10 inside the COT
//
//  fTrackEffHist[0]
//  fTrackEffHist[1]
//  fTrackEffHist[2]
//  fTrackEffHist[3]
//-----------------------------------------------------------------------------
  enum { kNEventHistSets   = 100 };
  enum { kNHelixHistSets   = 100 };
  enum { kNTrackSeedHistSets = 100 };
  enum { kNTrackHistSets   = 400 };
  enum { kNTimeClusterHistSets = 100};
  enum { kNClusterHistSets = 100 };
  enum { kNCaloHistSets    = 100 };
  enum { kNGenpHistSets    = 100 };
  enum { kNSimpHistSets    = 100 };

  struct Hist_t {
    TH1F*            fCrystalR  [2];	          // crystal radius
    EventHist_t*     fEvent     [kNEventHistSets];
    TrackHist_t*     fTrack     [kNTrackHistSets];
    TimeClusterHist_t* fTimeCluster  [kNTimeClusterHistSets];
    TrackSeedHist_t* fTrackSeed [kNTrackSeedHistSets];
    HelixHist_t*     fHelix     [kNHelixHistSets];
    ClusterHist_t*   fCluster   [kNClusterHistSets];
    GenpHist_t*      fGenp      [kNGenpHistSets];
    SimpHist_t*      fSimp      [kNSimpHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTimeClusterBlock*  fTimeClusterBlock;
  TStnTrackSeedBlock*    fTrackSeedBlock;
  TStnHelixBlock*        fHelixBlock;
  TStnTrackBlock*        fTrackBlock;
  TStnClusterBlock*      fClusterBlock;
  TCalDataBlock*         fCalDataBlock;
  TStrawDataBlock*       fStrawDataBlock;
  TGenpBlock*            fGenpBlock;
  TSimpBlock*            fSimpBlock;
  TTrackStrawHitBlock*   fTrackSHBlock;
  TStepPointMCBlock*     fStepPointMCBlock; //??

					// additional track parameters (assume ntracks < 20)
  TrackPar_t        fTrackPar[20];
					// histograms filled
  Hist_t            fHist;
					// cut values
  double            fPtMin;

  TGenParticle*     fParticle;		// electron or muon
  int               fPdgCode;		// determines which one
  int               fGeneratorCode;      

  TSimParticle*     fSimp;
  double            fEleE;		// electron energy

  int               fCalorimeterType;

  int               fNClusters;
  int               fNTimeClusters [5];
  int               fNTrackSeeds[5];
  int               fNHelices[5];
  int               fNTracks[10];
  int               fNGoodTracks;
  int               fNMatchedTracks;
  int               fNStrawHits;
  int               fNCalHits;
  int               fNGenp;		// N(generated particles)

  int               fNHyp;
  int               fBestHyp[10];
  int               fFillDioHist;
					// fTrackNumber[i]: track number, 
					// corresponding to OBSP particle #i
					// or -1
  TStnArrayI        fTrackNumber;

  TStnTrack*        fTrack;
  TStnCluster*      fCluster;

  TStnTrackID*      fTrackID;
  TEmuLogLH*        fLogLH;

  double            fMinT0;
  
  TString           fTrackBlockName;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TValidationModule2(const char* name="Validation2", const char* title="Validation");
  ~TValidationModule2();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
  TStnTrackBlock*    GetTrackBlock  () { return fTrackBlock;   }
  TStnClusterBlock*  GetClusterBlock() { return fClusterBlock; }

  TStnTrackID*       GetTrackID     () { return fTrackID; }
  TEmuLogLH*         GetLogLH       () { return fLogLH; }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void               SetFillDioHist  (int YesNo) { fFillDioHist   = YesNo; }
  void               SetPdgCode      (int Code ) { fPdgCode       = Code ; }
  void               SetGeneratorCode(int Code ) { fGeneratorCode = Code ; }

  void               SetTrackBlockName(const char* Name ) { fTrackBlockName = Name ; }
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  int     BeginJob();
  int     BeginRun();
  int     Event   (int ientry);
  int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookClusterHistograms (ClusterHist_t* Hist, const char* Folder);
  void    BookGenpHistograms    (GenpHist_t*    Hist, const char* Folder);
  void    BookEventHistograms   (EventHist_t*   Hist, const char* Folder);
  void    BookSimpHistograms    (SimpHist_t*    Hist, const char* Folder);
  void    BookTrackHistograms   (TrackHist_t*   Hist, const char* Folder);

  void    BookTimeClusterHistograms (TimeClusterHist_t* Hist, const char* Folder);
  void    BookHelixHistograms       (HelixHist_t*       Hist, const char* Folder);
  void    BookTrackSeedHistograms   (TrackSeedHist_t*   Hist, const char* Folder);

  void    FillEventHistograms    (EventHist_t* Hist);
  void    FillClusterHistograms  (ClusterHist_t* Hist, TStnCluster*  Cluster);
  void    FillGenpHistograms     (GenpHist_t*    Hist, TGenParticle* Genp   );
  void    FillSimpHistograms     (SimpHist_t*    Hist, TSimParticle* Simp   );
  void    FillTrackHistograms    (TrackHist_t*   Hist, TStnTrack*    Trk    );

  void    FillTrackSeedHistograms  (TrackSeedHist_t*   Hist, TStnTrackSeed*    TrkSeed);
  void    FillHelixHistograms      (HelixHist_t*       Hist, TStnHelix*        Helix  );
  void    FillTimeClusterHistograms(TimeClusterHist_t* Hist, TStnTimeCluster*  TimeCluster);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TValidationModule2,0)
};

#endif
