//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat May  7 16:24:47 2016 by ROOT version 6.04/14
// from TTree CollectionTree/xAOD event tree
// found on file: data15_13TeV_h012pre2_scalar_or_graviton.root
//////////////////////////////////////////////////////////

#ifndef HGammaMxAOD_h
#define HGammaMxAOD_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;
using std::vector;

class HGammaMxAOD {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   /*
   const Int_t kMaxHGamPhotonsAux = 1;
   const Int_t kMaxHGamElectronsAux = 1;
   const Int_t kMaxHGamAntiKt4EMTopoJetsAux = 1;
   const Int_t kMaxHGamMuonsAux = 1;
   const Int_t kMaxHGamMET_Reference_AntiKt4EMTopoAux = 1;
   const Int_t kMaxHGamMuonsInJetsAux = 1;
   const Int_t kMaxHGamEventInfoAux = 1;
   const Int_t kMaxEventInfoAux = 1;
   */

   // Declaration of leaf types
   //DataVector<xAOD::Photon_v1> *HGamPhotons;
   //xAOD::AuxContainerBase *HGamPhotonsAux_;
   vector<float>   *HGamPhotonsAuxDyn_topoetcone20;
   vector<float>   *HGamPhotonsAuxDyn_ptcone40;
   vector<float>   *HGamPhotonsAuxDyn_topoetcone40;
   vector<float>   *HGamPhotonsAuxDyn_eta_s2;
   vector<char>    *HGamPhotonsAuxDyn_isTight;
   vector<unsigned int> *HGamPhotonsAuxDyn_isEMTight;
   vector<float>   *HGamPhotonsAuxDyn_maxEcell_energy;
   vector<int>     *HGamPhotonsAuxDyn_maxEcell_gain;
   vector<float>   *HGamPhotonsAuxDyn_pt;
   vector<float>   *HGamPhotonsAuxDyn_relEreso;
   vector<float>   *HGamPhotonsAuxDyn_eta;
   vector<int>     *HGamPhotonsAuxDyn_conversionType;
   vector<float>   *HGamPhotonsAuxDyn_phi;
   vector<float>   *HGamPhotonsAuxDyn_m;
   vector<float>   *HGamPhotonsAuxDyn_conversionRadius;
   vector<float>   *HGamPhotonsAuxDyn_ratioE1E2;
   vector<char>    *HGamPhotonsAuxDyn_isIsoFixedCutLoose;
   vector<char>    *HGamPhotonsAuxDyn_isIsoFixedCutTight;
   vector<char>    *HGamPhotonsAuxDyn_isIsoFixedCutTightCaloOnly;
   vector<float>   *HGamPhotonsAuxDyn_E0_raw;
   vector<char>    *HGamPhotonsAuxDyn_isIsoFixedCutLooseCaloOnly;
   vector<float>   *HGamPhotonsAuxDyn_E1_raw;
   vector<float>   *HGamPhotonsAuxDyn_E2_raw;
   vector<float>   *HGamPhotonsAuxDyn_E3_raw;
   vector<float>   *HGamPhotonsAuxDyn_ptcone20_original;
   vector<float>   *HGamPhotonsAuxDyn_ptcone20_corr;
   vector<float>   *HGamPhotonsAuxDyn_ptcone20;
   //DataVector<xAOD::Electron_v1> *HGamElectrons;
   //xAOD::AuxContainerBase *HGamElectronsAux_;
   //DataVector<xAOD::Jet_v1> *HGamAntiKt4EMTopoJets;
   //xAOD::AuxContainerBase *HGamAntiKt4EMTopoJetsAux_;
   vector<char>    *HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_77;
   vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_Jvt;
   vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_DetectorEta;
   vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_pt;
   vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_eta;
   vector<char>    *HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_60;
   vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_phi;
   vector<char>    *HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_85;
   vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_m;
   vector<char>    *HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_70;
   //DataVector<xAOD::Muon_v1> *HGamMuons;
   //xAOD::AuxContainerBase *HGamMuonsAux_;
   //xAOD::MissingETContainer_v1 *HGamMET_Reference_AntiKt4EMTopo;
   //xAOD::AuxContainerBase *HGamMET_Reference_AntiKt4EMTopoAux_;
   vector<ULong64_t> *HGamMET_Reference_AntiKt4EMTopoAuxDyn_source;
   vector<double>  *HGamMET_Reference_AntiKt4EMTopoAuxDyn_sumet;
   vector<double>  *HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpx;
   vector<double>  *HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpy;
   vector<string>  *HGamMET_Reference_AntiKt4EMTopoAuxDyn_name;
   //DataVector<xAOD::Muon_v1> *HGamMuonsInJets;
   //xAOD::AuxContainerBase *HGamMuonsInJetsAux_;
   //xAOD::EventInfo_v1 *HGamEventInfo;
   //xAOD::AuxInfoBase *HGamEventInfoAux_;
   Float_t         HGamEventInfoAuxDyn_yybb_m_yybb;
   Float_t         HGamEventInfoAuxDyn_yybb_weight;
   Int_t           HGamEventInfoAuxDyn_yybb_bTagCat;
   Int_t           HGamEventInfoAuxDyn_yybb_cutFlow;
   Float_t         HGamEventInfoAuxDyn_yAbs_yy;
   Float_t         HGamEventInfoAuxDyn_m_yy_resolution;
   Float_t         HGamEventInfoAuxDyn_pTt_yy;
   Int_t           HGamEventInfoAuxDyn_NLoosePhotons;
   Float_t         HGamEventInfoAuxDyn_m_yy;
   Char_t          HGamEventInfoAuxDyn_passMeyCut;
   Int_t           HGamEventInfoAuxDyn_N_j_btag;
   Float_t         HGamEventInfoAuxDyn_pT_yy;
   Int_t           HGamEventInfoAuxDyn_N_j_btag30;
   Float_t         HGamEventInfoAuxDyn_pT_y1;
   Float_t         HGamEventInfoAuxDyn_m_yy_hardestVertex;
   Float_t         HGamEventInfoAuxDyn_pT_y2;
   Float_t         HGamEventInfoAuxDyn_m_yy_zCommon;
   Float_t         HGamEventInfoAuxDyn_E_y1;
   Char_t          HGamEventInfoAuxDyn_isPassedPreselection;
   Float_t         HGamEventInfoAuxDyn_E_y2;
   Char_t          HGamEventInfoAuxDyn_isPassedPID;
   Float_t         HGamEventInfoAuxDyn_pT_hard;
   Char_t          HGamEventInfoAuxDyn_isPassedIsolation;
   Float_t         HGamEventInfoAuxDyn_cosTS_yy;
   Char_t          HGamEventInfoAuxDyn_isPassedRelPtCuts;
   Char_t          HGamEventInfoAuxDyn_isPassedMassCut;
   Char_t          HGamEventInfoAuxDyn_isPassedIsolationMoriond;
   Float_t         HGamEventInfoAuxDyn_Dy_y_y;
   Char_t          HGamEventInfoAuxDyn_isPassedMoriond;
   Int_t           HGamEventInfoAuxDyn_N_e;
   Char_t          HGamEventInfoAuxDyn_isPassedIsolationLowHighMyy;
   Int_t           HGamEventInfoAuxDyn_N_mu;
   Char_t          HGamEventInfoAuxDyn_isPassedIsolationLowHighMyyMoriond;
   Int_t           HGamEventInfoAuxDyn_N_j;
   Char_t          HGamEventInfoAuxDyn_isPassedRelPtCutsLowHighMyy;
   Char_t          HGamEventInfoAuxDyn_isPassedLowHighMyy;
   Int_t           HGamEventInfoAuxDyn_N_j_central;
   Char_t          HGamEventInfoAuxDyn_isPassedLowHighMyyMoriond;
   Int_t           HGamEventInfoAuxDyn_N_j_central30;
   Char_t          HGamEventInfoAuxDyn_isPassedIsolationExotic;
   Float_t         HGamEventInfoAuxDyn_pT_j1;
   Char_t          HGamEventInfoAuxDyn_isPassedlPtCutsExotic;
   Float_t         HGamEventInfoAuxDyn_pT_j2;
   Char_t          HGamEventInfoAuxDyn_isPassedExotic;
   Float_t         HGamEventInfoAuxDyn_pT_jj;
   Float_t         HGamEventInfoAuxDyn_m_jj;
   Float_t         HGamEventInfoAuxDyn_Dy_j_j;
   Float_t         HGamEventInfoAuxDyn_Dphi_j_j;
   Float_t         HGamEventInfoAuxDyn_Dphi_yy_jj;
   Float_t         HGamEventInfoAuxDyn_m_ee;
   Float_t         HGamEventInfoAuxDyn_m_mumu;
   Float_t         HGamEventInfoAuxDyn_DRmin_y_j;
   Float_t         HGamEventInfoAuxDyn_DR_y_y;
   Float_t         HGamEventInfoAuxDyn_Zepp;
   Float_t         HGamEventInfoAuxDyn_cosTS_yyjj;
   Float_t         HGamEventInfoAuxDyn_met_TST;
   Float_t         HGamEventInfoAuxDyn_sumet_TST;
   Float_t         HGamEventInfoAuxDyn_phi_TST;
   Char_t          HGamEventInfoAuxDyn_isPassedBasic;
   Char_t          HGamEventInfoAuxDyn_isPassed;
   Char_t          HGamEventInfoAuxDyn_isPassedJetEventClean;
   Int_t           HGamEventInfoAuxDyn_cutFlow;
   Float_t         HGamEventInfoAuxDyn_weightInitial;
   Float_t         HGamEventInfoAuxDyn_weight;
   Float_t         HGamEventInfoAuxDyn_vertexWeight;
   Float_t         HGamEventInfoAuxDyn_pileupWeight;
   Float_t         HGamEventInfoAuxDyn_weightCatCoup_dev;
   Int_t           HGamEventInfoAuxDyn_catCoup_dev;
   Float_t         HGamEventInfoAuxDyn_weightCatCoup_Moriond2016;
   Int_t           HGamEventInfoAuxDyn_catCoup_Moriond2016;
   Int_t           HGamEventInfoAuxDyn_numberOfPrimaryVertices;
   Float_t         HGamEventInfoAuxDyn_selectedVertexZ;
   Float_t         HGamEventInfoAuxDyn_hardestVertexZ;
   Float_t         HGamEventInfoAuxDyn_zCommon;
   Float_t         HGamEventInfoAuxDyn_eventShapeDensity;
   Float_t         HGamEventInfoAuxDyn_mu;
   //xAOD::EventInfo_v1 *EventInfo;
 //xAOD::EventAuxInfo_v1 *EventInfoAux_;
 //xAOD::EventAuxInfo_v1 *EventInfoAux_xAOD__AuxInfoBase;
   UInt_t          EventInfoAux_runNumber;
   ULong64_t       EventInfoAux_eventNumber;
   UInt_t          EventInfoAux_lumiBlock;
   UInt_t          EventInfoAux_timeStamp;
   UInt_t          EventInfoAux_timeStampNSOffset;
   UInt_t          EventInfoAux_bcid;
   UInt_t          EventInfoAux_detectorMask0;
   UInt_t          EventInfoAux_detectorMask1;
   UInt_t          EventInfoAux_detectorMask2;
   UInt_t          EventInfoAux_detectorMask3;
 //vector<pair<string,string> > EventInfoAux_detDescrTags;
   UInt_t          EventInfoAux_eventTypeBitmask;
   UInt_t          EventInfoAux_statusElement;
   UInt_t          EventInfoAux_extendedLevel1ID;
   UShort_t        EventInfoAux_level1TriggerType;
   vector<string>  EventInfoAux_streamTagNames;
   vector<string>  EventInfoAux_streamTagTypes;
   vector<char>    EventInfoAux_streamTagObeysLumiblock;
   Float_t         EventInfoAux_actualInteractionsPerCrossing;
   Float_t         EventInfoAux_averageInteractionsPerCrossing;
   UInt_t          EventInfoAux_pixelFlags;
   UInt_t          EventInfoAux_sctFlags;
   UInt_t          EventInfoAux_trtFlags;
   UInt_t          EventInfoAux_larFlags;
   UInt_t          EventInfoAux_tileFlags;
   UInt_t          EventInfoAux_muonFlags;
   UInt_t          EventInfoAux_forwardDetFlags;
   UInt_t          EventInfoAux_coreFlags;
   UInt_t          EventInfoAux_backgroundFlags;
   UInt_t          EventInfoAux_lumiFlags;
   Float_t         EventInfoAux_beamPosX;
   Float_t         EventInfoAux_beamPosY;
   Float_t         EventInfoAux_beamPosZ;
   Float_t         EventInfoAux_beamPosSigmaX;
   Float_t         EventInfoAux_beamPosSigmaY;
   Float_t         EventInfoAux_beamPosSigmaZ;
   Float_t         EventInfoAux_beamPosSigmaXY;
   Float_t         EventInfoAux_beamTiltXZ;
   Float_t         EventInfoAux_beamTiltYZ;
   UInt_t          EventInfoAux_beamStatus;
   Char_t          EventInfoAuxDyn_passTrig_HLT_g35_loose_g25_loose;
   Char_t          EventInfoAuxDyn_passTrig_HLT_g35_medium_g25_medium;
   Char_t          EventInfoAuxDyn_passTrig_HLT_2g50_loose;
   Float_t         HGamEventInfoAuxDyn_Dy_yy_jj;
   vector<float>   *HGamMuonsInJetsAuxDyn_charge;
   vector<float>   *HGamMuonsInJetsAuxDyn_ptvarcone20;
   vector<float>   *HGamMuonsInJetsAuxDyn_pt;
   vector<unsigned short> *HGamMuonsInJetsAuxDyn_muonType;
   vector<float>   *HGamMuonsInJetsAuxDyn_eta;
   vector<float>   *HGamMuonsInJetsAuxDyn_phi;
   vector<char>    *HGamMuonsInJetsAuxDyn_passIPCut;
   vector<float>   *HGamMuonsInJetsAuxDyn_topoetcone20;
   vector<float>   *HGamElectronsAuxDyn_charge;
   vector<float>   *HGamElectronsAuxDyn_ptvarcone20;
   vector<float>   *HGamElectronsAuxDyn_pt;
   vector<float>   *HGamElectronsAuxDyn_eta;
   vector<float>   *HGamElectronsAuxDyn_phi;
   vector<float>   *HGamElectronsAuxDyn_eta_s2;
   vector<float>   *HGamElectronsAuxDyn_m;
   vector<char>    *HGamElectronsAuxDyn_isTight;
   vector<float>   *HGamElectronsAuxDyn_topoetcone20;
   vector<float>   *HGamMuonsAuxDyn_charge;
   vector<float>   *HGamMuonsAuxDyn_ptvarcone20;
   vector<float>   *HGamMuonsAuxDyn_pt;
   vector<char>    *HGamMuonsAuxDyn_passIPCut;
   vector<unsigned short> *HGamMuonsAuxDyn_muonType;
   vector<float>   *HGamMuonsAuxDyn_eta;
   vector<float>   *HGamMuonsAuxDyn_phi;
   vector<float>   *HGamMuonsAuxDyn_topoetcone20;

   // List of branches
   //TBranch        *b_HGamPhotons;   //!
   //TBranch        *b_HGamPhotonsAux_;   //!
   TBranch        *b_HGamPhotonsAuxDyn_topoetcone20;   //!
   TBranch        *b_HGamPhotonsAuxDyn_ptcone40;   //!
   TBranch        *b_HGamPhotonsAuxDyn_topoetcone40;   //!
   TBranch        *b_HGamPhotonsAuxDyn_eta_s2;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isTight;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isEMTight;   //!
   TBranch        *b_HGamPhotonsAuxDyn_maxEcell_energy;   //!
   TBranch        *b_HGamPhotonsAuxDyn_maxEcell_gain;   //!
   TBranch        *b_HGamPhotonsAuxDyn_pt;   //!
   TBranch        *b_HGamPhotonsAuxDyn_relEreso;   //!
   TBranch        *b_HGamPhotonsAuxDyn_eta;   //!
   TBranch        *b_HGamPhotonsAuxDyn_conversionType;   //!
   TBranch        *b_HGamPhotonsAuxDyn_phi;   //!
   TBranch        *b_HGamPhotonsAuxDyn_m;   //!
   TBranch        *b_HGamPhotonsAuxDyn_conversionRadius;   //!
   TBranch        *b_HGamPhotonsAuxDyn_ratioE1E2;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isIsoFixedCutLoose;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isIsoFixedCutTight;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isIsoFixedCutTightCaloOnly;   //!
   TBranch        *b_HGamPhotonsAuxDyn_E0_raw;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isIsoFixedCutLooseCaloOnly;   //!
   TBranch        *b_HGamPhotonsAuxDyn_E1_raw;   //!
   TBranch        *b_HGamPhotonsAuxDyn_E2_raw;   //!
   TBranch        *b_HGamPhotonsAuxDyn_E3_raw;   //!
   TBranch        *b_HGamPhotonsAuxDyn_ptcone20_original;   //!
   TBranch        *b_HGamPhotonsAuxDyn_ptcone20_corr;   //!
   TBranch        *b_HGamPhotonsAuxDyn_ptcone20;   //!
   //TBranch        *b_HGamElectrons;   //!
   //TBranch        *b_HGamElectronsAux_;   //!
   //TBranch        *b_HGamAntiKt4EMTopoJets;   //!
   //TBranch        *b_HGamAntiKt4EMTopoJetsAux_;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_77;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_Jvt;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_DetectorEta;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_pt;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_eta;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_60;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_phi;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_85;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_m;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_70;   //!
   //TBranch        *b_HGamMuons;   //!
   //TBranch        *b_HGamMuonsAux_;   //!
   TBranch        *b_HGamMET_Reference_AntiKt4EMTopo;   //!
   TBranch        *b_HGamMET_Reference_AntiKt4EMTopoAux_;   //!
   TBranch        *b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_source;   //!
   TBranch        *b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_sumet;   //!
   TBranch        *b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpx;   //!
   TBranch        *b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpy;   //!
   TBranch        *b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_name;   //!
   //TBranch        *b_HGamMuonsInJets;   //!
   //TBranch        *b_HGamMuonsInJetsAux_;   //!
   //TBranch        *b_HGamEventInfo;   //!
   //TBranch        *b_HGamEventInfoAux_;   //!
   TBranch        *b_HGamEventInfoAuxDyn_yybb_m_yybb;   //!
   TBranch        *b_HGamEventInfoAuxDyn_yybb_weight;   //!
   TBranch        *b_HGamEventInfoAuxDyn_yybb_bTagCat;   //!
   TBranch        *b_HGamEventInfoAuxDyn_yybb_cutFlow;   //!
   TBranch        *b_HGamEventInfoAuxDyn_yAbs_yy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_m_yy_resolution;   //!
   TBranch        *b_HGamEventInfoAuxDyn_pTt_yy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_NLoosePhotons;   //!
   TBranch        *b_HGamEventInfoAuxDyn_m_yy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_passMeyCut;   //!
   TBranch        *b_HGamEventInfoAuxDyn_N_j_btag;   //!
   TBranch        *b_HGamEventInfoAuxDyn_pT_yy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_N_j_btag30;   //!
   TBranch        *b_HGamEventInfoAuxDyn_pT_y1;   //!
   TBranch        *b_HGamEventInfoAuxDyn_m_yy_hardestVertex;   //!
   TBranch        *b_HGamEventInfoAuxDyn_pT_y2;   //!
   TBranch        *b_HGamEventInfoAuxDyn_m_yy_zCommon;   //!
   TBranch        *b_HGamEventInfoAuxDyn_E_y1;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedPreselection;   //!
   TBranch        *b_HGamEventInfoAuxDyn_E_y2;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedPID;   //!
   TBranch        *b_HGamEventInfoAuxDyn_pT_hard;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedIsolation;   //!
   TBranch        *b_HGamEventInfoAuxDyn_cosTS_yy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedRelPtCuts;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedMassCut;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedIsolationMoriond;   //!
   TBranch        *b_HGamEventInfoAuxDyn_Dy_y_y;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedMoriond;   //!
   TBranch        *b_HGamEventInfoAuxDyn_N_e;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedIsolationLowHighMyy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_N_mu;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedIsolationLowHighMyyMoriond;   //!
   TBranch        *b_HGamEventInfoAuxDyn_N_j;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedRelPtCutsLowHighMyy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedLowHighMyy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_N_j_central;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedLowHighMyyMoriond;   //!
   TBranch        *b_HGamEventInfoAuxDyn_N_j_central30;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedIsolationExotic;   //!
   TBranch        *b_HGamEventInfoAuxDyn_pT_j1;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedlPtCutsExotic;   //!
   TBranch        *b_HGamEventInfoAuxDyn_pT_j2;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedExotic;   //!
   TBranch        *b_HGamEventInfoAuxDyn_pT_jj;   //!
   TBranch        *b_HGamEventInfoAuxDyn_m_jj;   //!
   TBranch        *b_HGamEventInfoAuxDyn_Dy_j_j;   //!
   TBranch        *b_HGamEventInfoAuxDyn_Dphi_j_j;   //!
   TBranch        *b_HGamEventInfoAuxDyn_Dphi_yy_jj;   //!
   TBranch        *b_HGamEventInfoAuxDyn_m_ee;   //!
   TBranch        *b_HGamEventInfoAuxDyn_m_mumu;   //!
   TBranch        *b_HGamEventInfoAuxDyn_DRmin_y_j;   //!
   TBranch        *b_HGamEventInfoAuxDyn_DR_y_y;   //!
   TBranch        *b_HGamEventInfoAuxDyn_Zepp;   //!
   TBranch        *b_HGamEventInfoAuxDyn_cosTS_yyjj;   //!
   TBranch        *b_HGamEventInfoAuxDyn_met_TST;   //!
   TBranch        *b_HGamEventInfoAuxDyn_sumet_TST;   //!
   TBranch        *b_HGamEventInfoAuxDyn_phi_TST;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedBasic;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassed;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedJetEventClean;   //!
   TBranch        *b_HGamEventInfoAuxDyn_cutFlow;   //!
   TBranch        *b_HGamEventInfoAuxDyn_weightInitial;   //!
   TBranch        *b_HGamEventInfoAuxDyn_weight;   //!
   TBranch        *b_HGamEventInfoAuxDyn_vertexWeight;   //!
   TBranch        *b_HGamEventInfoAuxDyn_pileupWeight;   //!
   TBranch        *b_HGamEventInfoAuxDyn_weightCatCoup_dev;   //!
   TBranch        *b_HGamEventInfoAuxDyn_catCoup_dev;   //!
   TBranch        *b_HGamEventInfoAuxDyn_weightCatCoup_Moriond2016;   //!
   TBranch        *b_HGamEventInfoAuxDyn_catCoup_Moriond2016;   //!
   TBranch        *b_HGamEventInfoAuxDyn_numberOfPrimaryVertices;   //!
   TBranch        *b_HGamEventInfoAuxDyn_selectedVertexZ;   //!
   TBranch        *b_HGamEventInfoAuxDyn_hardestVertexZ;   //!
   TBranch        *b_HGamEventInfoAuxDyn_zCommon;   //!
   TBranch        *b_HGamEventInfoAuxDyn_eventShapeDensity;   //!
   TBranch        *b_HGamEventInfoAuxDyn_mu;   //!
   //TBranch        *b_EventInfo;   //!
   TBranch        *b_EventInfoAux_runNumber;   //!
   TBranch        *b_EventInfoAux_eventNumber;   //!
   TBranch        *b_EventInfoAux_lumiBlock;   //!
   TBranch        *b_EventInfoAux_timeStamp;   //!
   TBranch        *b_EventInfoAux_timeStampNSOffset;   //!
   TBranch        *b_EventInfoAux_bcid;   //!
   TBranch        *b_EventInfoAux_detectorMask0;   //!
   TBranch        *b_EventInfoAux_detectorMask1;   //!
   TBranch        *b_EventInfoAux_detectorMask2;   //!
   TBranch        *b_EventInfoAux_detectorMask3;   //!
   TBranch        *b_EventInfoAux_eventTypeBitmask;   //!
   TBranch        *b_EventInfoAux_statusElement;   //!
   TBranch        *b_EventInfoAux_extendedLevel1ID;   //!
   TBranch        *b_EventInfoAux_level1TriggerType;   //!
   TBranch        *b_EventInfoAux_streamTagNames;   //!
   TBranch        *b_EventInfoAux_streamTagTypes;   //!
   TBranch        *b_EventInfoAux_streamTagObeysLumiblock;   //!
   TBranch        *b_EventInfoAux_actualInteractionsPerCrossing;   //!
   TBranch        *b_EventInfoAux_averageInteractionsPerCrossing;   //!
   TBranch        *b_EventInfoAux_pixelFlags;   //!
   TBranch        *b_EventInfoAux_sctFlags;   //!
   TBranch        *b_EventInfoAux_trtFlags;   //!
   TBranch        *b_EventInfoAux_larFlags;   //!
   TBranch        *b_EventInfoAux_tileFlags;   //!
   TBranch        *b_EventInfoAux_muonFlags;   //!
   TBranch        *b_EventInfoAux_forwardDetFlags;   //!
   TBranch        *b_EventInfoAux_coreFlags;   //!
   TBranch        *b_EventInfoAux_backgroundFlags;   //!
   TBranch        *b_EventInfoAux_lumiFlags;   //!
   TBranch        *b_EventInfoAux_beamPosX;   //!
   TBranch        *b_EventInfoAux_beamPosY;   //!
   TBranch        *b_EventInfoAux_beamPosZ;   //!
   TBranch        *b_EventInfoAux_beamPosSigmaX;   //!
   TBranch        *b_EventInfoAux_beamPosSigmaY;   //!
   TBranch        *b_EventInfoAux_beamPosSigmaZ;   //!
   TBranch        *b_EventInfoAux_beamPosSigmaXY;   //!
   TBranch        *b_EventInfoAux_beamTiltXZ;   //!
   TBranch        *b_EventInfoAux_beamTiltYZ;   //!
   TBranch        *b_EventInfoAux_beamStatus;   //!
   TBranch        *b_EventInfoAuxDyn_passTrig_HLT_g35_loose_g25_loose;   //!
   TBranch        *b_EventInfoAuxDyn_passTrig_HLT_g35_medium_g25_medium;   //!
   TBranch        *b_EventInfoAuxDyn_passTrig_HLT_2g50_loose;   //!
   TBranch        *b_HGamEventInfoAuxDyn_Dy_yy_jj;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_charge;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_ptvarcone20;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_pt;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_muonType;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_eta;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_phi;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_passIPCut;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_topoetcone20;   //!
   TBranch        *b_HGamElectronsAuxDyn_charge;   //!
   TBranch        *b_HGamElectronsAuxDyn_ptvarcone20;   //!
   TBranch        *b_HGamElectronsAuxDyn_pt;   //!
   TBranch        *b_HGamElectronsAuxDyn_eta;   //!
   TBranch        *b_HGamElectronsAuxDyn_phi;   //!
   TBranch        *b_HGamElectronsAuxDyn_eta_s2;   //!
   TBranch        *b_HGamElectronsAuxDyn_m;   //!
   TBranch        *b_HGamElectronsAuxDyn_isTight;   //!
   TBranch        *b_HGamElectronsAuxDyn_topoetcone20;   //!
   TBranch        *b_HGamMuonsAuxDyn_charge;   //!
   TBranch        *b_HGamMuonsAuxDyn_ptvarcone20;   //!
   TBranch        *b_HGamMuonsAuxDyn_pt;   //!
   TBranch        *b_HGamMuonsAuxDyn_passIPCut;   //!
   TBranch        *b_HGamMuonsAuxDyn_muonType;   //!
   TBranch        *b_HGamMuonsAuxDyn_eta;   //!
   TBranch        *b_HGamMuonsAuxDyn_phi;   //!
   TBranch        *b_HGamMuonsAuxDyn_topoetcone20;   //!

   HGammaMxAOD(TTree *tree=0, TString tag="");
   virtual ~HGammaMxAOD();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, TString tag);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HGammaMxAOD_cxx
HGammaMxAOD::HGammaMxAOD(TTree *tree, TString tag) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data15_13TeV_h012pre2_scalar_or_graviton.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data15_13TeV_h012pre2_scalar_or_graviton.root");
      }
      f->GetObject("CollectionTree",tree);

   }
   Init(tree, tag);
}

HGammaMxAOD::~HGammaMxAOD()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HGammaMxAOD::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HGammaMxAOD::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HGammaMxAOD::Init(TTree *tree, TString tag)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   //HGamPhotons = 0;
   //HGamPhotonsAux_ = 0;
   HGamPhotonsAuxDyn_topoetcone20 = 0;
   HGamPhotonsAuxDyn_ptcone40 = 0;
   HGamPhotonsAuxDyn_topoetcone40 = 0;
   HGamPhotonsAuxDyn_eta_s2 = 0;
   HGamPhotonsAuxDyn_isTight = 0;
   HGamPhotonsAuxDyn_isEMTight = 0;
   HGamPhotonsAuxDyn_maxEcell_energy = 0;
   HGamPhotonsAuxDyn_maxEcell_gain = 0;
   HGamPhotonsAuxDyn_pt = 0;
   HGamPhotonsAuxDyn_relEreso = 0;
   HGamPhotonsAuxDyn_eta = 0;
   HGamPhotonsAuxDyn_conversionType = 0;
   HGamPhotonsAuxDyn_phi = 0;
   HGamPhotonsAuxDyn_m = 0;
   HGamPhotonsAuxDyn_conversionRadius = 0;
   HGamPhotonsAuxDyn_ratioE1E2 = 0;
   HGamPhotonsAuxDyn_isIsoFixedCutLoose = 0;
   HGamPhotonsAuxDyn_isIsoFixedCutTight = 0;
   HGamPhotonsAuxDyn_isIsoFixedCutTightCaloOnly = 0;
   HGamPhotonsAuxDyn_E0_raw = 0;
   HGamPhotonsAuxDyn_isIsoFixedCutLooseCaloOnly = 0;
   HGamPhotonsAuxDyn_E1_raw = 0;
   HGamPhotonsAuxDyn_E2_raw = 0;
   HGamPhotonsAuxDyn_E3_raw = 0;
   HGamPhotonsAuxDyn_ptcone20_original = 0;
   HGamPhotonsAuxDyn_ptcone20_corr = 0;
   HGamPhotonsAuxDyn_ptcone20 = 0;
   //HGamElectrons = 0;
   //HGamElectronsAux_ = 0;
   //HGamAntiKt4EMTopoJets = 0;
   //HGamAntiKt4EMTopoJetsAux_ = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_77 = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_Jvt = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_DetectorEta = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_pt = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_eta = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_60 = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_phi = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_85 = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_m = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_70 = 0;
   //HGamMuons = 0;
   //HGamMuonsAux_ = 0;
   //HGamMET_Reference_AntiKt4EMTopo = 0;
   //HGamMET_Reference_AntiKt4EMTopoAux_ = 0;
   HGamMET_Reference_AntiKt4EMTopoAuxDyn_source = 0;
   HGamMET_Reference_AntiKt4EMTopoAuxDyn_sumet = 0;
   HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpx = 0;
   HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpy = 0;
   HGamMET_Reference_AntiKt4EMTopoAuxDyn_name = 0;
   //HGamMuonsInJets = 0;
   //HGamMuonsInJetsAux_ = 0;
   //HGamEventInfo = 0;
   //HGamEventInfoAux_ = 0;
   //EventInfo = 0;
   HGamMuonsInJetsAuxDyn_charge = 0;
   HGamMuonsInJetsAuxDyn_ptvarcone20 = 0;
   HGamMuonsInJetsAuxDyn_pt = 0;
   HGamMuonsInJetsAuxDyn_muonType = 0;
   HGamMuonsInJetsAuxDyn_eta = 0;
   HGamMuonsInJetsAuxDyn_phi = 0;
   HGamMuonsInJetsAuxDyn_passIPCut = 0;
   HGamMuonsInJetsAuxDyn_topoetcone20 = 0;
   HGamElectronsAuxDyn_charge = 0;
   HGamElectronsAuxDyn_ptvarcone20 = 0;
   HGamElectronsAuxDyn_pt = 0;
   HGamElectronsAuxDyn_eta = 0;
   HGamElectronsAuxDyn_phi = 0;
   HGamElectronsAuxDyn_eta_s2 = 0;
   HGamElectronsAuxDyn_m = 0;
   HGamElectronsAuxDyn_isTight = 0;
   HGamElectronsAuxDyn_topoetcone20 = 0;
   HGamMuonsAuxDyn_charge = 0;
   HGamMuonsAuxDyn_ptvarcone20 = 0;
   HGamMuonsAuxDyn_pt = 0;
   HGamMuonsAuxDyn_passIPCut = 0;
   HGamMuonsAuxDyn_muonType = 0;
   HGamMuonsAuxDyn_eta = 0;
   HGamMuonsAuxDyn_phi = 0;
   HGamMuonsAuxDyn_topoetcone20 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   //fChain->SetBranchAddress("HGamPhotons", &HGamPhotons, &b_HGamPhotons);
   //fChain->SetBranchAddress("HGamPhotonsAux.", &HGamPhotonsAux_, &b_HGamPhotonsAux_);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.topoetcone20", &HGamPhotonsAuxDyn_topoetcone20, &b_HGamPhotonsAuxDyn_topoetcone20);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.ptcone40", &HGamPhotonsAuxDyn_ptcone40, &b_HGamPhotonsAuxDyn_ptcone40);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.topoetcone40", &HGamPhotonsAuxDyn_topoetcone40, &b_HGamPhotonsAuxDyn_topoetcone40);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.eta_s2", &HGamPhotonsAuxDyn_eta_s2, &b_HGamPhotonsAuxDyn_eta_s2);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isTight", &HGamPhotonsAuxDyn_isTight, &b_HGamPhotonsAuxDyn_isTight);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isEMTight", &HGamPhotonsAuxDyn_isEMTight, &b_HGamPhotonsAuxDyn_isEMTight);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.maxEcell_energy", &HGamPhotonsAuxDyn_maxEcell_energy, &b_HGamPhotonsAuxDyn_maxEcell_energy);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.maxEcell_gain", &HGamPhotonsAuxDyn_maxEcell_gain, &b_HGamPhotonsAuxDyn_maxEcell_gain);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.pt", &HGamPhotonsAuxDyn_pt, &b_HGamPhotonsAuxDyn_pt);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.relEreso", &HGamPhotonsAuxDyn_relEreso, &b_HGamPhotonsAuxDyn_relEreso);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.eta", &HGamPhotonsAuxDyn_eta, &b_HGamPhotonsAuxDyn_eta);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.conversionType", &HGamPhotonsAuxDyn_conversionType, &b_HGamPhotonsAuxDyn_conversionType);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.phi", &HGamPhotonsAuxDyn_phi, &b_HGamPhotonsAuxDyn_phi);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.m", &HGamPhotonsAuxDyn_m, &b_HGamPhotonsAuxDyn_m);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.conversionRadius", &HGamPhotonsAuxDyn_conversionRadius, &b_HGamPhotonsAuxDyn_conversionRadius);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.ratioE1E2", &HGamPhotonsAuxDyn_ratioE1E2, &b_HGamPhotonsAuxDyn_ratioE1E2);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isIsoFixedCutLoose", &HGamPhotonsAuxDyn_isIsoFixedCutLoose, &b_HGamPhotonsAuxDyn_isIsoFixedCutLoose);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isIsoFixedCutTight", &HGamPhotonsAuxDyn_isIsoFixedCutTight, &b_HGamPhotonsAuxDyn_isIsoFixedCutTight);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isIsoFixedCutTightCaloOnly", &HGamPhotonsAuxDyn_isIsoFixedCutTightCaloOnly, &b_HGamPhotonsAuxDyn_isIsoFixedCutTightCaloOnly);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.E0_raw", &HGamPhotonsAuxDyn_E0_raw, &b_HGamPhotonsAuxDyn_E0_raw);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isIsoFixedCutLooseCaloOnly", &HGamPhotonsAuxDyn_isIsoFixedCutLooseCaloOnly, &b_HGamPhotonsAuxDyn_isIsoFixedCutLooseCaloOnly);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.E1_raw", &HGamPhotonsAuxDyn_E1_raw, &b_HGamPhotonsAuxDyn_E1_raw);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.E2_raw", &HGamPhotonsAuxDyn_E2_raw, &b_HGamPhotonsAuxDyn_E2_raw);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.E3_raw", &HGamPhotonsAuxDyn_E3_raw, &b_HGamPhotonsAuxDyn_E3_raw);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.ptcone20_original", &HGamPhotonsAuxDyn_ptcone20_original, &b_HGamPhotonsAuxDyn_ptcone20_original);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.ptcone20_corr", &HGamPhotonsAuxDyn_ptcone20_corr, &b_HGamPhotonsAuxDyn_ptcone20_corr);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.ptcone20", &HGamPhotonsAuxDyn_ptcone20, &b_HGamPhotonsAuxDyn_ptcone20);
   //fChain->SetBranchAddress("HGamElectrons", &HGamElectrons, &b_HGamElectrons);
   //fChain->SetBranchAddress("HGamElectronsAux.", &HGamElectronsAux_, &b_HGamElectronsAux_);
   //fChain->SetBranchAddress("HGamAntiKt4EMTopoJets", &HGamAntiKt4EMTopoJets, &b_HGamAntiKt4EMTopoJets);
   //fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAux.", &HGamAntiKt4EMTopoJetsAux_, &b_HGamAntiKt4EMTopoJetsAux_);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.MV2c20_FixedCutBEff_77", &HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_77, &b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_77);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.Jvt", &HGamAntiKt4EMTopoJetsAuxDyn_Jvt, &b_HGamAntiKt4EMTopoJetsAuxDyn_Jvt);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.DetectorEta", &HGamAntiKt4EMTopoJetsAuxDyn_DetectorEta, &b_HGamAntiKt4EMTopoJetsAuxDyn_DetectorEta);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.pt", &HGamAntiKt4EMTopoJetsAuxDyn_pt, &b_HGamAntiKt4EMTopoJetsAuxDyn_pt);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.eta", &HGamAntiKt4EMTopoJetsAuxDyn_eta, &b_HGamAntiKt4EMTopoJetsAuxDyn_eta);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.MV2c20_FixedCutBEff_60", &HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_60, &b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_60);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.phi", &HGamAntiKt4EMTopoJetsAuxDyn_phi, &b_HGamAntiKt4EMTopoJetsAuxDyn_phi);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.MV2c20_FixedCutBEff_85", &HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_85, &b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_85);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.m", &HGamAntiKt4EMTopoJetsAuxDyn_m, &b_HGamAntiKt4EMTopoJetsAuxDyn_m);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.MV2c20_FixedCutBEff_70", &HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_70, &b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_FixedCutBEff_70);
   //fChain->SetBranchAddress("HGamMuons", &HGamMuons, &b_HGamMuons);
   //fChain->SetBranchAddress("HGamMuonsAux.", &HGamMuonsAux_, &b_HGamMuonsAux_);
   //fChain->SetBranchAddress("HGamMET_Reference_AntiKt4EMTopo", &HGamMET_Reference_AntiKt4EMTopo, &b_HGamMET_Reference_AntiKt4EMTopo);
   //fChain->SetBranchAddress("HGamMET_Reference_AntiKt4EMTopoAux.", &HGamMET_Reference_AntiKt4EMTopoAux_, &b_HGamMET_Reference_AntiKt4EMTopoAux_);
   fChain->SetBranchAddress("HGamMET_Reference_AntiKt4EMTopoAuxDyn.source", &HGamMET_Reference_AntiKt4EMTopoAuxDyn_source, &b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_source);
   fChain->SetBranchAddress("HGamMET_Reference_AntiKt4EMTopoAuxDyn.sumet", &HGamMET_Reference_AntiKt4EMTopoAuxDyn_sumet, &b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_sumet);
   fChain->SetBranchAddress("HGamMET_Reference_AntiKt4EMTopoAuxDyn.mpx", &HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpx, &b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpx);
   fChain->SetBranchAddress("HGamMET_Reference_AntiKt4EMTopoAuxDyn.mpy", &HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpy, &b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpy);
   fChain->SetBranchAddress("HGamMET_Reference_AntiKt4EMTopoAuxDyn.name", &HGamMET_Reference_AntiKt4EMTopoAuxDyn_name, &b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_name);
   //fChain->SetBranchAddress("HGamMuonsInJets", &HGamMuonsInJets, &b_HGamMuonsInJets);
   //fChain->SetBranchAddress("HGamMuonsInJetsAux.", &HGamMuonsInJetsAux_, &b_HGamMuonsInJetsAux_);
   //fChain->SetBranchAddress("HGamEventInfo", &HGamEventInfo, &b_HGamEventInfo);
   //fChain->SetBranchAddress("HGamEventInfoAux.", &HGamEventInfoAux_, &b_HGamEventInfoAux_);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.yybb_m_yybb", &HGamEventInfoAuxDyn_yybb_m_yybb, &b_HGamEventInfoAuxDyn_yybb_m_yybb);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.yybb_weight", &HGamEventInfoAuxDyn_yybb_weight, &b_HGamEventInfoAuxDyn_yybb_weight);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.yybb_bTagCat", &HGamEventInfoAuxDyn_yybb_bTagCat, &b_HGamEventInfoAuxDyn_yybb_bTagCat);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.yybb_cutFlow", &HGamEventInfoAuxDyn_yybb_cutFlow, &b_HGamEventInfoAuxDyn_yybb_cutFlow);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.yAbs_yy", &HGamEventInfoAuxDyn_yAbs_yy, &b_HGamEventInfoAuxDyn_yAbs_yy);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.m_yy_resolution", &HGamEventInfoAuxDyn_m_yy_resolution, &b_HGamEventInfoAuxDyn_m_yy_resolution);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.pTt_yy", &HGamEventInfoAuxDyn_pTt_yy, &b_HGamEventInfoAuxDyn_pTt_yy);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.NLoosePhotons", &HGamEventInfoAuxDyn_NLoosePhotons, &b_HGamEventInfoAuxDyn_NLoosePhotons);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.m_yy", &HGamEventInfoAuxDyn_m_yy, &b_HGamEventInfoAuxDyn_m_yy);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.passMeyCut", &HGamEventInfoAuxDyn_passMeyCut, &b_HGamEventInfoAuxDyn_passMeyCut);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.N_j_btag", &HGamEventInfoAuxDyn_N_j_btag, &b_HGamEventInfoAuxDyn_N_j_btag);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.pT_yy", &HGamEventInfoAuxDyn_pT_yy, &b_HGamEventInfoAuxDyn_pT_yy);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.N_j_btag30", &HGamEventInfoAuxDyn_N_j_btag30, &b_HGamEventInfoAuxDyn_N_j_btag30);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.pT_y1", &HGamEventInfoAuxDyn_pT_y1, &b_HGamEventInfoAuxDyn_pT_y1);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.m_yy_hardestVertex", &HGamEventInfoAuxDyn_m_yy_hardestVertex, &b_HGamEventInfoAuxDyn_m_yy_hardestVertex);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.pT_y2", &HGamEventInfoAuxDyn_pT_y2, &b_HGamEventInfoAuxDyn_pT_y2);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.m_yy_zCommon", &HGamEventInfoAuxDyn_m_yy_zCommon, &b_HGamEventInfoAuxDyn_m_yy_zCommon);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.E_y1", &HGamEventInfoAuxDyn_E_y1, &b_HGamEventInfoAuxDyn_E_y1);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedPreselection", &HGamEventInfoAuxDyn_isPassedPreselection, &b_HGamEventInfoAuxDyn_isPassedPreselection);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.E_y2", &HGamEventInfoAuxDyn_E_y2, &b_HGamEventInfoAuxDyn_E_y2);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedPID", &HGamEventInfoAuxDyn_isPassedPID, &b_HGamEventInfoAuxDyn_isPassedPID);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.pT_hard", &HGamEventInfoAuxDyn_pT_hard, &b_HGamEventInfoAuxDyn_pT_hard);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedIsolation", &HGamEventInfoAuxDyn_isPassedIsolation, &b_HGamEventInfoAuxDyn_isPassedIsolation);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.cosTS_yy", &HGamEventInfoAuxDyn_cosTS_yy, &b_HGamEventInfoAuxDyn_cosTS_yy);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedRelPtCuts", &HGamEventInfoAuxDyn_isPassedRelPtCuts, &b_HGamEventInfoAuxDyn_isPassedRelPtCuts);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedMassCut", &HGamEventInfoAuxDyn_isPassedMassCut, &b_HGamEventInfoAuxDyn_isPassedMassCut);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedIsolationMoriond", &HGamEventInfoAuxDyn_isPassedIsolationMoriond, &b_HGamEventInfoAuxDyn_isPassedIsolationMoriond);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.Dy_y_y", &HGamEventInfoAuxDyn_Dy_y_y, &b_HGamEventInfoAuxDyn_Dy_y_y);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedMoriond", &HGamEventInfoAuxDyn_isPassedMoriond, &b_HGamEventInfoAuxDyn_isPassedMoriond);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.N_e", &HGamEventInfoAuxDyn_N_e, &b_HGamEventInfoAuxDyn_N_e);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedIsolationLowHighMyy", &HGamEventInfoAuxDyn_isPassedIsolationLowHighMyy, &b_HGamEventInfoAuxDyn_isPassedIsolationLowHighMyy);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.N_mu", &HGamEventInfoAuxDyn_N_mu, &b_HGamEventInfoAuxDyn_N_mu);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedIsolationLowHighMyyMoriond", &HGamEventInfoAuxDyn_isPassedIsolationLowHighMyyMoriond, &b_HGamEventInfoAuxDyn_isPassedIsolationLowHighMyyMoriond);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.N_j", &HGamEventInfoAuxDyn_N_j, &b_HGamEventInfoAuxDyn_N_j);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedRelPtCutsLowHighMyy", &HGamEventInfoAuxDyn_isPassedRelPtCutsLowHighMyy, &b_HGamEventInfoAuxDyn_isPassedRelPtCutsLowHighMyy);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedLowHighMyy", &HGamEventInfoAuxDyn_isPassedLowHighMyy, &b_HGamEventInfoAuxDyn_isPassedLowHighMyy);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.N_j_central", &HGamEventInfoAuxDyn_N_j_central, &b_HGamEventInfoAuxDyn_N_j_central);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedLowHighMyyMoriond", &HGamEventInfoAuxDyn_isPassedLowHighMyyMoriond, &b_HGamEventInfoAuxDyn_isPassedLowHighMyyMoriond);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.N_j_central30", &HGamEventInfoAuxDyn_N_j_central30, &b_HGamEventInfoAuxDyn_N_j_central30);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedIsolationExotic", &HGamEventInfoAuxDyn_isPassedIsolationExotic, &b_HGamEventInfoAuxDyn_isPassedIsolationExotic);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.pT_j1", &HGamEventInfoAuxDyn_pT_j1, &b_HGamEventInfoAuxDyn_pT_j1);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedlPtCutsExotic", &HGamEventInfoAuxDyn_isPassedlPtCutsExotic, &b_HGamEventInfoAuxDyn_isPassedlPtCutsExotic);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.pT_j2", &HGamEventInfoAuxDyn_pT_j2, &b_HGamEventInfoAuxDyn_pT_j2);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedExotic", &HGamEventInfoAuxDyn_isPassedExotic, &b_HGamEventInfoAuxDyn_isPassedExotic);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.pT_jj", &HGamEventInfoAuxDyn_pT_jj, &b_HGamEventInfoAuxDyn_pT_jj);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.m_jj", &HGamEventInfoAuxDyn_m_jj, &b_HGamEventInfoAuxDyn_m_jj);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.Dy_j_j", &HGamEventInfoAuxDyn_Dy_j_j, &b_HGamEventInfoAuxDyn_Dy_j_j);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.Dphi_j_j", &HGamEventInfoAuxDyn_Dphi_j_j, &b_HGamEventInfoAuxDyn_Dphi_j_j);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.Dphi_yy_jj", &HGamEventInfoAuxDyn_Dphi_yy_jj, &b_HGamEventInfoAuxDyn_Dphi_yy_jj);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.m_ee", &HGamEventInfoAuxDyn_m_ee, &b_HGamEventInfoAuxDyn_m_ee);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.m_mumu", &HGamEventInfoAuxDyn_m_mumu, &b_HGamEventInfoAuxDyn_m_mumu);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.DRmin_y_j", &HGamEventInfoAuxDyn_DRmin_y_j, &b_HGamEventInfoAuxDyn_DRmin_y_j);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.DR_y_y", &HGamEventInfoAuxDyn_DR_y_y, &b_HGamEventInfoAuxDyn_DR_y_y);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.Zepp", &HGamEventInfoAuxDyn_Zepp, &b_HGamEventInfoAuxDyn_Zepp);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.cosTS_yyjj", &HGamEventInfoAuxDyn_cosTS_yyjj, &b_HGamEventInfoAuxDyn_cosTS_yyjj);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.met_TST", &HGamEventInfoAuxDyn_met_TST, &b_HGamEventInfoAuxDyn_met_TST);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.sumet_TST", &HGamEventInfoAuxDyn_sumet_TST, &b_HGamEventInfoAuxDyn_sumet_TST);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.phi_TST", &HGamEventInfoAuxDyn_phi_TST, &b_HGamEventInfoAuxDyn_phi_TST);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedBasic", &HGamEventInfoAuxDyn_isPassedBasic, &b_HGamEventInfoAuxDyn_isPassedBasic);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassed", &HGamEventInfoAuxDyn_isPassed, &b_HGamEventInfoAuxDyn_isPassed);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedJetEventClean", &HGamEventInfoAuxDyn_isPassedJetEventClean, &b_HGamEventInfoAuxDyn_isPassedJetEventClean);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.cutFlow", &HGamEventInfoAuxDyn_cutFlow, &b_HGamEventInfoAuxDyn_cutFlow);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.weightInitial", &HGamEventInfoAuxDyn_weightInitial, &b_HGamEventInfoAuxDyn_weightInitial);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.weight", &HGamEventInfoAuxDyn_weight, &b_HGamEventInfoAuxDyn_weight);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.vertexWeight", &HGamEventInfoAuxDyn_vertexWeight, &b_HGamEventInfoAuxDyn_vertexWeight);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.pileupWeight", &HGamEventInfoAuxDyn_pileupWeight, &b_HGamEventInfoAuxDyn_pileupWeight);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.weightCatCoup_dev", &HGamEventInfoAuxDyn_weightCatCoup_dev, &b_HGamEventInfoAuxDyn_weightCatCoup_dev);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.catCoup_dev", &HGamEventInfoAuxDyn_catCoup_dev, &b_HGamEventInfoAuxDyn_catCoup_dev);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.weightCatCoup_Moriond2016", &HGamEventInfoAuxDyn_weightCatCoup_Moriond2016, &b_HGamEventInfoAuxDyn_weightCatCoup_Moriond2016);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.catCoup_Moriond2016", &HGamEventInfoAuxDyn_catCoup_Moriond2016, &b_HGamEventInfoAuxDyn_catCoup_Moriond2016);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.numberOfPrimaryVertices", &HGamEventInfoAuxDyn_numberOfPrimaryVertices, &b_HGamEventInfoAuxDyn_numberOfPrimaryVertices);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.selectedVertexZ", &HGamEventInfoAuxDyn_selectedVertexZ, &b_HGamEventInfoAuxDyn_selectedVertexZ);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.hardestVertexZ", &HGamEventInfoAuxDyn_hardestVertexZ, &b_HGamEventInfoAuxDyn_hardestVertexZ);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.zCommon", &HGamEventInfoAuxDyn_zCommon, &b_HGamEventInfoAuxDyn_zCommon);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.eventShapeDensity", &HGamEventInfoAuxDyn_eventShapeDensity, &b_HGamEventInfoAuxDyn_eventShapeDensity);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.mu", &HGamEventInfoAuxDyn_mu, &b_HGamEventInfoAuxDyn_mu);
   //fChain->SetBranchAddress("EventInfo", &EventInfo, &b_EventInfo);
   fChain->SetBranchAddress("EventInfoAux.runNumber", &EventInfoAux_runNumber, &b_EventInfoAux_runNumber);
   fChain->SetBranchAddress("EventInfoAux.eventNumber", &EventInfoAux_eventNumber, &b_EventInfoAux_eventNumber);
   fChain->SetBranchAddress("EventInfoAux.lumiBlock", &EventInfoAux_lumiBlock, &b_EventInfoAux_lumiBlock);
   fChain->SetBranchAddress("EventInfoAux.timeStamp", &EventInfoAux_timeStamp, &b_EventInfoAux_timeStamp);
   fChain->SetBranchAddress("EventInfoAux.timeStampNSOffset", &EventInfoAux_timeStampNSOffset, &b_EventInfoAux_timeStampNSOffset);
   fChain->SetBranchAddress("EventInfoAux.bcid", &EventInfoAux_bcid, &b_EventInfoAux_bcid);
   fChain->SetBranchAddress("EventInfoAux.detectorMask0", &EventInfoAux_detectorMask0, &b_EventInfoAux_detectorMask0);
   fChain->SetBranchAddress("EventInfoAux.detectorMask1", &EventInfoAux_detectorMask1, &b_EventInfoAux_detectorMask1);
   fChain->SetBranchAddress("EventInfoAux.detectorMask2", &EventInfoAux_detectorMask2, &b_EventInfoAux_detectorMask2);
   fChain->SetBranchAddress("EventInfoAux.detectorMask3", &EventInfoAux_detectorMask3, &b_EventInfoAux_detectorMask3);
   fChain->SetBranchAddress("EventInfoAux.eventTypeBitmask", &EventInfoAux_eventTypeBitmask, &b_EventInfoAux_eventTypeBitmask);
   fChain->SetBranchAddress("EventInfoAux.statusElement", &EventInfoAux_statusElement, &b_EventInfoAux_statusElement);
   fChain->SetBranchAddress("EventInfoAux.extendedLevel1ID", &EventInfoAux_extendedLevel1ID, &b_EventInfoAux_extendedLevel1ID);
   fChain->SetBranchAddress("EventInfoAux.level1TriggerType", &EventInfoAux_level1TriggerType, &b_EventInfoAux_level1TriggerType);
   fChain->SetBranchAddress("EventInfoAux.streamTagNames", &EventInfoAux_streamTagNames, &b_EventInfoAux_streamTagNames);
   fChain->SetBranchAddress("EventInfoAux.streamTagTypes", &EventInfoAux_streamTagTypes, &b_EventInfoAux_streamTagTypes);
   fChain->SetBranchAddress("EventInfoAux.streamTagObeysLumiblock", &EventInfoAux_streamTagObeysLumiblock, &b_EventInfoAux_streamTagObeysLumiblock);
   fChain->SetBranchAddress("EventInfoAux.actualInteractionsPerCrossing", &EventInfoAux_actualInteractionsPerCrossing, &b_EventInfoAux_actualInteractionsPerCrossing);
   fChain->SetBranchAddress("EventInfoAux.averageInteractionsPerCrossing", &EventInfoAux_averageInteractionsPerCrossing, &b_EventInfoAux_averageInteractionsPerCrossing);
   fChain->SetBranchAddress("EventInfoAux.pixelFlags", &EventInfoAux_pixelFlags, &b_EventInfoAux_pixelFlags);
   fChain->SetBranchAddress("EventInfoAux.sctFlags", &EventInfoAux_sctFlags, &b_EventInfoAux_sctFlags);
   fChain->SetBranchAddress("EventInfoAux.trtFlags", &EventInfoAux_trtFlags, &b_EventInfoAux_trtFlags);
   fChain->SetBranchAddress("EventInfoAux.larFlags", &EventInfoAux_larFlags, &b_EventInfoAux_larFlags);
   fChain->SetBranchAddress("EventInfoAux.tileFlags", &EventInfoAux_tileFlags, &b_EventInfoAux_tileFlags);
   fChain->SetBranchAddress("EventInfoAux.muonFlags", &EventInfoAux_muonFlags, &b_EventInfoAux_muonFlags);
   fChain->SetBranchAddress("EventInfoAux.forwardDetFlags", &EventInfoAux_forwardDetFlags, &b_EventInfoAux_forwardDetFlags);
   fChain->SetBranchAddress("EventInfoAux.coreFlags", &EventInfoAux_coreFlags, &b_EventInfoAux_coreFlags);
   fChain->SetBranchAddress("EventInfoAux.backgroundFlags", &EventInfoAux_backgroundFlags, &b_EventInfoAux_backgroundFlags);
   fChain->SetBranchAddress("EventInfoAux.lumiFlags", &EventInfoAux_lumiFlags, &b_EventInfoAux_lumiFlags);
   fChain->SetBranchAddress("EventInfoAux.beamPosX", &EventInfoAux_beamPosX, &b_EventInfoAux_beamPosX);
   fChain->SetBranchAddress("EventInfoAux.beamPosY", &EventInfoAux_beamPosY, &b_EventInfoAux_beamPosY);
   fChain->SetBranchAddress("EventInfoAux.beamPosZ", &EventInfoAux_beamPosZ, &b_EventInfoAux_beamPosZ);
   fChain->SetBranchAddress("EventInfoAux.beamPosSigmaX", &EventInfoAux_beamPosSigmaX, &b_EventInfoAux_beamPosSigmaX);
   fChain->SetBranchAddress("EventInfoAux.beamPosSigmaY", &EventInfoAux_beamPosSigmaY, &b_EventInfoAux_beamPosSigmaY);
   fChain->SetBranchAddress("EventInfoAux.beamPosSigmaZ", &EventInfoAux_beamPosSigmaZ, &b_EventInfoAux_beamPosSigmaZ);
   fChain->SetBranchAddress("EventInfoAux.beamPosSigmaXY", &EventInfoAux_beamPosSigmaXY, &b_EventInfoAux_beamPosSigmaXY);
   fChain->SetBranchAddress("EventInfoAux.beamTiltXZ", &EventInfoAux_beamTiltXZ, &b_EventInfoAux_beamTiltXZ);
   fChain->SetBranchAddress("EventInfoAux.beamTiltYZ", &EventInfoAux_beamTiltYZ, &b_EventInfoAux_beamTiltYZ);
   fChain->SetBranchAddress("EventInfoAux.beamStatus", &EventInfoAux_beamStatus, &b_EventInfoAux_beamStatus);
   fChain->SetBranchAddress("EventInfoAuxDyn.passTrig_HLT_g35_loose_g25_loose", &EventInfoAuxDyn_passTrig_HLT_g35_loose_g25_loose, &b_EventInfoAuxDyn_passTrig_HLT_g35_loose_g25_loose);
   fChain->SetBranchAddress("EventInfoAuxDyn.passTrig_HLT_g35_medium_g25_medium", &EventInfoAuxDyn_passTrig_HLT_g35_medium_g25_medium, &b_EventInfoAuxDyn_passTrig_HLT_g35_medium_g25_medium);
   fChain->SetBranchAddress("EventInfoAuxDyn.passTrig_HLT_2g50_loose", &EventInfoAuxDyn_passTrig_HLT_2g50_loose, &b_EventInfoAuxDyn_passTrig_HLT_2g50_loose);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.Dy_yy_jj", &HGamEventInfoAuxDyn_Dy_yy_jj, &b_HGamEventInfoAuxDyn_Dy_yy_jj);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.charge", &HGamMuonsInJetsAuxDyn_charge, &b_HGamMuonsInJetsAuxDyn_charge);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.ptvarcone20", &HGamMuonsInJetsAuxDyn_ptvarcone20, &b_HGamMuonsInJetsAuxDyn_ptvarcone20);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.pt", &HGamMuonsInJetsAuxDyn_pt, &b_HGamMuonsInJetsAuxDyn_pt);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.muonType", &HGamMuonsInJetsAuxDyn_muonType, &b_HGamMuonsInJetsAuxDyn_muonType);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.eta", &HGamMuonsInJetsAuxDyn_eta, &b_HGamMuonsInJetsAuxDyn_eta);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.phi", &HGamMuonsInJetsAuxDyn_phi, &b_HGamMuonsInJetsAuxDyn_phi);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.passIPCut", &HGamMuonsInJetsAuxDyn_passIPCut, &b_HGamMuonsInJetsAuxDyn_passIPCut);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.topoetcone20", &HGamMuonsInJetsAuxDyn_topoetcone20, &b_HGamMuonsInJetsAuxDyn_topoetcone20);
   fChain->SetBranchAddress("HGamElectronsAuxDyn.charge", &HGamElectronsAuxDyn_charge, &b_HGamElectronsAuxDyn_charge);
   fChain->SetBranchAddress("HGamElectronsAuxDyn.ptvarcone20", &HGamElectronsAuxDyn_ptvarcone20, &b_HGamElectronsAuxDyn_ptvarcone20);
   fChain->SetBranchAddress("HGamElectronsAuxDyn.pt", &HGamElectronsAuxDyn_pt, &b_HGamElectronsAuxDyn_pt);
   fChain->SetBranchAddress("HGamElectronsAuxDyn.eta", &HGamElectronsAuxDyn_eta, &b_HGamElectronsAuxDyn_eta);
   fChain->SetBranchAddress("HGamElectronsAuxDyn.phi", &HGamElectronsAuxDyn_phi, &b_HGamElectronsAuxDyn_phi);
   fChain->SetBranchAddress("HGamElectronsAuxDyn.eta_s2", &HGamElectronsAuxDyn_eta_s2, &b_HGamElectronsAuxDyn_eta_s2);
   fChain->SetBranchAddress("HGamElectronsAuxDyn.m", &HGamElectronsAuxDyn_m, &b_HGamElectronsAuxDyn_m);
   fChain->SetBranchAddress("HGamElectronsAuxDyn.isTight", &HGamElectronsAuxDyn_isTight, &b_HGamElectronsAuxDyn_isTight);
   fChain->SetBranchAddress("HGamElectronsAuxDyn.topoetcone20", &HGamElectronsAuxDyn_topoetcone20, &b_HGamElectronsAuxDyn_topoetcone20);
   fChain->SetBranchAddress("HGamMuonsAuxDyn.charge", &HGamMuonsAuxDyn_charge, &b_HGamMuonsAuxDyn_charge);
   fChain->SetBranchAddress("HGamMuonsAuxDyn.ptvarcone20", &HGamMuonsAuxDyn_ptvarcone20, &b_HGamMuonsAuxDyn_ptvarcone20);
   fChain->SetBranchAddress("HGamMuonsAuxDyn.pt", &HGamMuonsAuxDyn_pt, &b_HGamMuonsAuxDyn_pt);
   fChain->SetBranchAddress("HGamMuonsAuxDyn.passIPCut", &HGamMuonsAuxDyn_passIPCut, &b_HGamMuonsAuxDyn_passIPCut);
   fChain->SetBranchAddress("HGamMuonsAuxDyn.muonType", &HGamMuonsAuxDyn_muonType, &b_HGamMuonsAuxDyn_muonType);
   fChain->SetBranchAddress("HGamMuonsAuxDyn.eta", &HGamMuonsAuxDyn_eta, &b_HGamMuonsAuxDyn_eta);
   fChain->SetBranchAddress("HGamMuonsAuxDyn.phi", &HGamMuonsAuxDyn_phi, &b_HGamMuonsAuxDyn_phi);
   fChain->SetBranchAddress("HGamMuonsAuxDyn.topoetcone20", &HGamMuonsAuxDyn_topoetcone20, &b_HGamMuonsAuxDyn_topoetcone20);
   Notify();
}

Bool_t HGammaMxAOD::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HGammaMxAOD::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HGammaMxAOD::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HGammaMxAOD_cxx
