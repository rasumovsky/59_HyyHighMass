//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar  8 11:22:22 2016 by ROOT version 6.04/14
// from TTree toy/toy
// found on file: toy_mu0.root
//////////////////////////////////////////////////////////

#ifndef ToyTree_h
#define ToyTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;
using std::vector;

class ToyTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           seed;
   Int_t           bestFitUpdate;
   Double_t        numEvents;
   vector<double>  *numEventsPerCate;
   vector<double>  *nllPerRetry;
   Double_t        profiledPOIVal;
   Bool_t          convergedMu0;
   Bool_t          convergedMu1;
   Bool_t          convergedMuFree;
   Double_t        nllMu0;
   Double_t        nllMu1;
   Double_t        nllMuFree;
   Double_t        llrL1L0;
   Double_t        llrL0Lfree;
   Double_t        llrL1Lfree;
   vector<string>  *namesNP;
   vector<double>  *valuesNPMu0;
   vector<double>  *valuesNPMu1;
   vector<double>  *valuesNPMuFree;
   vector<string>  *namesGlobs;
   vector<double>  *valuesGlobsMu1;
   vector<double>  *valuesGlobsMu0;
   vector<double>  *valuesGlobsMuFree;
   vector<string>  *namesPoIs;
   vector<double>  *valuesPoIsMu1;
   vector<double>  *valuesPoIsMu0;
   vector<double>  *valuesPoIsMuFree;

   // List of branches
   TBranch        *b_seed;   //!
   TBranch        *b_bestFitUpdate;   //!
   TBranch        *b_numEvents;   //!
   TBranch        *b_numEventsPerCate;   //!
   TBranch        *b_nllPerRetry;   //!
   TBranch        *b_profiledPOIVal;   //!
   TBranch        *b_convergedMu0;   //!
   TBranch        *b_convergedMu1;   //!
   TBranch        *b_convergedMuFree;   //!
   TBranch        *b_nllMu0;   //!
   TBranch        *b_nllMu1;   //!
   TBranch        *b_nllMuFree;   //!
   TBranch        *b_llrL1L0;   //!
   TBranch        *b_llrL0Lfree;   //!
   TBranch        *b_llrL1Lfree;   //!
   TBranch        *b_namesNP;   //!
   TBranch        *b_valuesNPMu0;   //!
   TBranch        *b_valuesNPMu1;   //!
   TBranch        *b_valuesNPMuFree;   //!
   TBranch        *b_namesGlobs;   //!
   TBranch        *b_valuesGlobsMu1;   //!
   TBranch        *b_valuesGlobsMu0;   //!
   TBranch        *b_valuesGlobsMuFree;   //!
   TBranch        *b_namesPoIs;   //!
   TBranch        *b_valuesPoIsMu1;   //!
   TBranch        *b_valuesPoIsMu0;   //!
   TBranch        *b_valuesPoIsMuFree;   //!

   ToyTree(TTree *tree=0);
   virtual ~ToyTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ToyTree_cxx
ToyTree::ToyTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("toy_mu0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("toy_mu0.root");
      }
      f->GetObject("toy",tree);

   }
   Init(tree);
}

ToyTree::~ToyTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ToyTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ToyTree::LoadTree(Long64_t entry)
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

void ToyTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   numEventsPerCate = 0;
   nllPerRetry = 0;
   namesNP = 0;
   valuesNPMu0 = 0;
   valuesNPMu1 = 0;
   valuesNPMuFree = 0;
   namesGlobs = 0;
   valuesGlobsMu1 = 0;
   valuesGlobsMu0 = 0;
   valuesGlobsMuFree = 0;
   namesPoIs = 0;
   valuesPoIsMu1 = 0;
   valuesPoIsMu0 = 0;
   valuesPoIsMuFree = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("seed", &seed, &b_seed);
   fChain->SetBranchAddress("bestFitUpdate", &bestFitUpdate, &b_bestFitUpdate);
   fChain->SetBranchAddress("numEvents", &numEvents, &b_numEvents);
   fChain->SetBranchAddress("numEventsPerCate", &numEventsPerCate, &b_numEventsPerCate);
   fChain->SetBranchAddress("nllPerRetry", &nllPerRetry, &b_nllPerRetry);
   fChain->SetBranchAddress("profiledPOIVal", &profiledPOIVal, &b_profiledPOIVal);
   fChain->SetBranchAddress("convergedMu0", &convergedMu0, &b_convergedMu0);
   fChain->SetBranchAddress("convergedMu1", &convergedMu1, &b_convergedMu1);
   fChain->SetBranchAddress("convergedMuFree", &convergedMuFree, &b_convergedMuFree);
   fChain->SetBranchAddress("nllMu0", &nllMu0, &b_nllMu0);
   fChain->SetBranchAddress("nllMu1", &nllMu1, &b_nllMu1);
   fChain->SetBranchAddress("nllMuFree", &nllMuFree, &b_nllMuFree);
   fChain->SetBranchAddress("llrL1L0", &llrL1L0, &b_llrL1L0);
   fChain->SetBranchAddress("llrL0Lfree", &llrL0Lfree, &b_llrL0Lfree);
   fChain->SetBranchAddress("llrL1Lfree", &llrL1Lfree, &b_llrL1Lfree);
   fChain->SetBranchAddress("namesNP", &namesNP, &b_namesNP);
   fChain->SetBranchAddress("valuesNPMu0", &valuesNPMu0, &b_valuesNPMu0);
   fChain->SetBranchAddress("valuesNPMu1", &valuesNPMu1, &b_valuesNPMu1);
   fChain->SetBranchAddress("valuesNPMuFree", &valuesNPMuFree, &b_valuesNPMuFree);
   fChain->SetBranchAddress("namesGlobs", &namesGlobs, &b_namesGlobs);
   fChain->SetBranchAddress("valuesGlobsMu1", &valuesGlobsMu1, &b_valuesGlobsMu1);
   fChain->SetBranchAddress("valuesGlobsMu0", &valuesGlobsMu0, &b_valuesGlobsMu0);
   fChain->SetBranchAddress("valuesGlobsMuFree", &valuesGlobsMuFree, &b_valuesGlobsMuFree);
   fChain->SetBranchAddress("namesPoIs", &namesPoIs, &b_namesPoIs);
   fChain->SetBranchAddress("valuesPoIsMu1", &valuesPoIsMu1, &b_valuesPoIsMu1);
   fChain->SetBranchAddress("valuesPoIsMu0", &valuesPoIsMu0, &b_valuesPoIsMu0);
   fChain->SetBranchAddress("valuesPoIsMuFree", &valuesPoIsMuFree, &b_valuesPoIsMuFree);
   Notify();
}

Bool_t ToyTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ToyTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ToyTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ToyTree_cxx
