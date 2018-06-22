#include "TString.h"
#include "TFile.h"
#include "TH1.h"

void removeOldSV(TString treeName, TString inputFile, TString outputFile) {
  TFile* oldFile = TFile::Open(inputFile);
  TFile* newFile = TFile::Open(outputFile, "recreate");

  TH1F*  nevents = (TH1F*)oldFile->Get("nevents");
  TTree* tree  = (TTree*)oldFile->Get(treeName);
  tree->SetBranchStatus("*",1);
  tree->SetBranchStatus("m_sv", 0);
  tree->SetBranchStatus("pt_sv", 0);
  tree->SetBranchStatus("*_sv_*", 0);
  
  TTree* newTree  = tree->CloneTree();
  nevents->Write();
  newTree->Write();
  delete oldFile;
  delete newFile;
}
  
