#include "PhysicsTools/FWLite/interface/CommandLineParser.h" 
#include "TFile.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TKey.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include <math.h> 
#include "TMath.h" 
#include <limits>
#include "TSystem.h"
#include <vector>
#include <string>

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

//If recoilType 0 then don't do recoil
//              FIXME amc@nlo is not ready yet!!! 1 then aMC@NLO DY and W+Jets MC samples
//                1 is not longer an option
//              2 MG5 DY and W+Jets MC samples or Higgs MC samples
//
//If doES       0 does not apply any ES shifts
//              1 applies ES shifts to TT channel, no effect on other channels
//
//If isWJets    0 no shift in number of jets used for recoil corrections
//              1 shifts njets + 1 for recoil corrections
//
//If metType    1 use mvamet
//        -1 use pf met

ClassicSVfit svfitAlgorithm;
bool tylerCode = false;

void copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) ;
void readdir(TDirectory *dir, optutl::CommandLineParser parser,  char TreeToUse[], int recoilType, int doES, int isWJets, int metType, double tesSize) ;
void CopyFile(const char *fname, optutl::CommandLineParser parser);
void CopyDir(TDirectory *source,optutl::CommandLineParser parser);

void runSVFit(std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons,
              double measuredMETx, double measuredMETy,
              TMatrixD &covMET, float num, float &svFitMass, float& svFitPt, float &svFitEta,
              float &svFitPhi, float &svFitMET, float &svFitTransverseMass
);

int main (int argc, char* argv[]) 
{
  optutl::CommandLineParser parser ("Sets Event Weights in the ntuple");
  parser.addOption("branch",optutl::CommandLineParser::kString,"Branch","__svFit__");
  parser.addOption("newFile",optutl::CommandLineParser::kString,"newFile","newFile.root");
  parser.addOption("inputFile",optutl::CommandLineParser::kString,"input File");
  parser.addOption("newOutputFile",optutl::CommandLineParser::kDouble,"New Output File",0.0);
  parser.addOption("recoilType",optutl::CommandLineParser::kDouble,"recoilType",0.0);
  parser.addOption("doES",optutl::CommandLineParser::kDouble,"doES",0.0);
  parser.addOption("isWJets",optutl::CommandLineParser::kDouble,"isWJets",0.0);
  parser.addOption("metType",optutl::CommandLineParser::kDouble,"metType",-1.0); // 1 = mvamet, -1 = pf met
  parser.addOption("tesSize",optutl::CommandLineParser::kDouble,"tesSize",0.012); // Default TES = 1.2%
  parser.addOption("numEvents",optutl::CommandLineParser::kInteger,"numEvents",-1);
  parser.parseArguments (argc, argv);
  
  std::cout << "EXTRA COMMANDS:"
        << "\n --- numEvents: " << parser.integerValue("numEvents")
        << "\n --- recoilType: " << parser.doubleValue("recoilType")
        << "\n --- doES: " << parser.doubleValue("doES")
        << "\n --- isWJets: " << parser.doubleValue("isWJets")
        << "\n --- metType: " << parser.doubleValue("metType")
        << "\n --- tesSize: " << parser.doubleValue("tesSize") << std::endl;
  
  // Make sure a proper Met Type is chosen
  assert (parser.doubleValue("metType") == 1.0 || parser.doubleValue("metType") == -1.0);
  // Temp: require pfMet
  assert (parser.doubleValue("metType") == -1.0);
  
  // No DiTauMass constraint
  svfitAlgorithm.setDiTauMassConstraint(-1.0);
  
  char TreeToUse[80]="first";
  
  TFile *fProduce;//= new TFile(parser.stringValue("newFile").c_str(),"UPDATE");
  
  if(parser.doubleValue("newOutputFile")>0){
    TFile *f = new TFile(parser.stringValue("inputFile").c_str(),"READ");
    std::cout<<"Creating new outputfile"<<std::endl;
    std::string newFileName = parser.stringValue("newFile");
    
    fProduce = new TFile(newFileName.c_str(),"RECREATE");
    copyFiles(parser, f, fProduce);//new TFile(parser.stringValue("inputFile").c_str()+"SVFit","UPDATE");
    fProduce = new TFile(newFileName.c_str(),"UPDATE");
    std::cout<<"listing the directories================="<<std::endl;
    fProduce->ls();
    readdir(fProduce,parser,TreeToUse,parser.doubleValue("recoilType"),parser.doubleValue("doES"),
            parser.doubleValue("isWJets"),parser.doubleValue("metType"),parser.doubleValue("tesSize"));
    
    fProduce->Close();
    f->Close();
  }
  else{
    TFile *f = new TFile(parser.stringValue("inputFile").c_str(),"UPDATE");
    readdir(f,parser,TreeToUse,parser.doubleValue("recoilType"),parser.doubleValue("doES"),
            parser.doubleValue("isWJets"),parser.doubleValue("metType"),parser.doubleValue("tesSize"));
    f->Close();
  }
  
  
} 


void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int recoilType, int doES, int isWJets, int metType, double tesSize) 
{
  
  
  std::string recoilFileName = "HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2016BtoH.root";
  if(recoilType == 1) { //amc@nlo
    std::cout << "Alexei no long specified MG vs. AMC@NLO, so use recoilType = 2" << std::endl;
    return; }
  if(recoilType == 2 && metType == 1) { // mva met (Alexei no long specified MG vs. AMC@NLO)
    std::cout << "Alexei does not provide full 2016 data recoil corrections for Mva Met\n\n" << std::endl;
    std::cout << "Using ICHEP Mva Met corrections\n\n" << std::endl;
    recoilFileName = "HTT-utilities/RecoilCorrections/data/MvaMET_2016BCD.root";}
  if(recoilType == 2 && metType == -1) { // pf met (Alexei no long specified MG vs. AMC@NLO)
    recoilFileName = "HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2016BtoH.root";}
  
  classic_svFit::MeasuredTauLepton::kDecayType decayType1 = classic_svFit::MeasuredTauLepton::kUndefinedDecayType;
  classic_svFit::MeasuredTauLepton::kDecayType decayType2 = classic_svFit::MeasuredTauLepton::kUndefinedDecayType; 
  
  // Both masses should depend on decay mode and particle?
  std::string channel = "x";
  float mass1 = 0.;
  float mass2 = 0.;
  int finalState = 0;

  TDirectory *dirsav = gDirectory;
  TKey *key;
  dir->cd();      

  std::vector<TString> processedNames;
  
  TIter next(dir->GetListOfKeys());
  while ((key = (TKey*)next())) {
    printf("Found key=%s \n",key->GetName());
    
    TObject *obj = key->ReadObj();
    if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
      std::cout << "This is a directory, diving in!" << std::endl;
      // zero the processedNames vector, to allow processing trees with duplicate names in separate directories
      processedNames.clear();

      dir->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      sprintf(TreeToUse,"%s",key->GetName());
      readdir(subdir,parser,TreeToUse,parser.doubleValue("recoilType"),parser.doubleValue("doES"),
          parser.doubleValue("isWJets"),parser.doubleValue("metType"),parser.doubleValue("tesSize"));
      
      dirsav->cd();
    }
    else if(obj->IsA()->InheritsFrom(TTree::Class())) {
      // check  if this tree was already processed
      std::vector<TString>::const_iterator it = find(processedNames.begin(), processedNames.end(), key->GetName());
      if ( it != processedNames.end() ) {
    std::cout << "This tree was already processed, skipping..." <<  std::endl;
    continue;
      }
      std::cout << "This is the tree! Start processing" << std::endl;
      processedNames.push_back(key->GetName());
     
    if((std::string(TreeToUse).find("eemt")!= std::string::npos) ||
                  (std::string(TreeToUse).find("mmmt")!= std::string::npos) ||
                  (parser.stringValue("inputFile").find("_eemt.root") != std::string::npos) ||
                  (parser.stringValue("inputFile").find("_mmmt.root") != std::string::npos)
              ){
             decayType1 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
             decayType2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
             mass1 = 0.105658;
             mass2 = 0;
             channel = "mt";
             std::cout << "Identified channel mt and using kappa = 4" << std::endl;
             svfitAlgorithm.addLogM_fixed(true, 4);
             if (parser.stringValue("inputFile").find("_eemt.root") != std::string::npos) finalState = 10;
             if (parser.stringValue("inputFile").find("_mmmt.root") != std::string::npos) finalState = 20;
    }

    else if((std::string(TreeToUse).find("eeet")!= std::string::npos) ||
                  (std::string(TreeToUse).find("emmt")!= std::string::npos) ||
                  (parser.stringValue("inputFile").find("_eeet.root") != std::string::npos) ||
                  (parser.stringValue("inputFile").find("_emmt.root") != std::string::npos)
              ){
             decayType1 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
             decayType2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
             mass1 = 0.00051100;
             mass2 = 0;
             channel = "et";
             std::cout << "Identified channel et and using kappa = 4" << std::endl;
             svfitAlgorithm.addLogM_fixed(true, 4);
             if (parser.stringValue("inputFile").find("_eeet.root") != std::string::npos) finalState = 11;
             if (parser.stringValue("inputFile").find("_emmt.root") != std::string::npos) finalState = 21;
    }

    else if((std::string(TreeToUse).find("eeem")!= std::string::npos) ||
                  (std::string(TreeToUse).find("emmm")!= std::string::npos) ||
                  (parser.stringValue("inputFile").find("_eeem.root") != std::string::npos) ||
                  (parser.stringValue("inputFile").find("_emmm.root") != std::string::npos)
              ){
             decayType1 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
             decayType2 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
             mass1 = 0.00051100;
             mass2 = 0.105658;
             channel = "em";
             std::cout << "Identified channel em and using kappa = 3" << std::endl;
             svfitAlgorithm.addLogM_fixed(true, 3);
             if (parser.stringValue("inputFile").find("_eeem.root") != std::string::npos) finalState = 12;
             if (parser.stringValue("inputFile").find("_emmm.root") != std::string::npos) finalState = 22;
    }
    else if((std::string(TreeToUse).find("eett")!= std::string::npos) ||
                  (std::string(TreeToUse).find("mmtt")!= std::string::npos) ||
                  (parser.stringValue("inputFile").find("_eett.root") != std::string::npos) ||
                  (parser.stringValue("inputFile").find("_mmtt.root") != std::string::npos)
              ){
             decayType1 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
             decayType2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
             mass1 = 0.13957;
             mass2 = 0.13957;
             channel = "tt";
             std::cout << "Identified channel tt and using kappa = 5" << std::endl;
             svfitAlgorithm.addLogM_fixed(true, 5);
             if (parser.stringValue("inputFile").find("_eett.root") != std::string::npos) finalState = 13;
             if (parser.stringValue("inputFile").find("_mmtt.root") != std::string::npos) finalState = 23;
    }
    else{
             std::cout<<"TreeToUse "<< std::string(TreeToUse)<<" does not match muTauEvent or eleTauEvent... Skipping!!"<<std::endl;
    }


    // Make sure we found a finalState code
    assert (
         finalState == 10.0
         || finalState == 11.0
         || finalState == 12.0
         || finalState == 13.0
         || finalState == 20.0
         || finalState == 21.0
         || finalState == 22.0
         || finalState == 23.0
    );





 
      TTree *t = (TTree*)obj;
      float svFitMass = -10;
      float svFitPt = -10;
      float svFitEta = -10;
      float svFitPhi = -10;
      float svFitMET = -10;
      float svFitTransverseMass = -10;

      float metcorr_ex = -10; // corrected met px (float)
      float metcorr_ey = -10;  // corrected met py (float)
      float metcorrUncUp_ex = -10; // corrUncUpected met px (float)
      float metcorrUncUp_ey = -10;  // corrUncUpected met py (float)
      float metcorrUncDown_ex = -10; // corrUncDownected met px (float)
      float metcorrUncDown_ey = -10;  // corrUncDownected met py (float)
      float metcorrClusteredUp_ex = -10; // corrClusteredUpected met px (float)
      float metcorrClusteredUp_ey = -10;  // corrClusteredUpected met py (float)
      float metcorrClusteredDown_ex = -10; // corrClusteredDownected met px (float)
      float metcorrClusteredDown_ey = -10;  // corrClusteredDownected met py (float)

      // For saving
      float metcor = -10; // corrected metcor
      float metcorphi = -10; // corrected metcorphi
      float metcorClusteredDown = -10;   
      float metcorphiClusteredDown = -10;
      float metcorClusteredUp = -10;     
      float metcorphiClusteredUp = -10;  
      float metcorUncDown = -10;         
      float metcorphiUncDown = -10;      
      float metcorUncUp = -10;           
      float metcorphiUncUp = -10;        

      TBranch *newBranch1 = t->Branch("m_sv", &svFitMass, "m_sv/F");
      TBranch *newBranch2 = t->Branch("pt_sv", &svFitPt, "pt_sv/F");
      TBranch *newBranch3 = t->Branch("eta_sv", &svFitEta, "eta_sv/F");
      TBranch *newBranch4 = t->Branch("phi_sv", &svFitPhi, "phi_sv/F");
      TBranch *newBranch5 = t->Branch("met_sv", &svFitMET, "met_sv/F");
      TBranch *newBranch6 = t->Branch("mt_sv", &svFitTransverseMass, "mt_sv/F");

      TBranch *newBranch7 = t->Branch("metcorr_ex", &metcorr_ex, "metcorr_ex/F");
      TBranch *newBranch8 = t->Branch("metcorr_ey", &metcorr_ey, "metcorr_ey/F");
      TBranch *newBranch9 = t->Branch("metcor", &metcor, "metcor/F");
      TBranch *newBranch10 = t->Branch("metcorphi", &metcorphi, "metcorphi/F");

      // If doing ES shifts, we need extra ouput branches
      float svFitMass_UP = -10;
      float svFitPt_UP = -10;
      float svFitEta_UP = -10;
      float svFitPhi_UP = -10;
      float svFitMET_UP = -10;
      float svFitTransverseMass_UP = -10;
      float svFitMass_DOWN = -10;
      float svFitPt_DOWN = -10;
      float svFitEta_DOWN = -10;
      float svFitPhi_DOWN = -10;
      float svFitMET_DOWN = -10;
      float svFitTransverseMass_DOWN = -10;

      float svFitMassEES_UP = -10;
      float svFitPtEES_UP = -10;
      float svFitEtaEES_UP = -10;
      float svFitPhiEES_UP = -10;
      float svFitMETEES_UP = -10;
      float svFitTransverseMassEES_UP = -10;
      float svFitMassEES_DOWN = -10;
      float svFitPtEES_DOWN = -10;
      float svFitEtaEES_DOWN = -10;
      float svFitPhiEES_DOWN = -10;
      float svFitMETEES_DOWN = -10;
      float svFitTransverseMassEES_DOWN = -10;

      float svFitMass_DM0_UP = -10;
      float svFitPt_DM0_UP = -10;
      float svFitEta_DM0_UP = -10;
      float svFitPhi_DM0_UP = -10;
      float svFitMET_DM0_UP = -10;
      float svFitTransverseMass_DM0_UP = -10;
      float svFitMass_DM0_DOWN = -10;
      float svFitPt_DM0_DOWN = -10;
      float svFitEta_DM0_DOWN = -10;
      float svFitPhi_DM0_DOWN = -10;
      float svFitMET_DM0_DOWN = -10;
      float svFitTransverseMass_DM0_DOWN = -10;

      float svFitMass_DM1_UP = -10;
      float svFitPt_DM1_UP = -10;
      float svFitEta_DM1_UP = -10;
      float svFitPhi_DM1_UP = -10;
      float svFitMET_DM1_UP = -10;
      float svFitTransverseMass_DM1_UP = -10;
      float svFitMass_DM1_DOWN = -10;
      float svFitPt_DM1_DOWN = -10;
      float svFitEta_DM1_DOWN = -10;
      float svFitPhi_DM1_DOWN = -10;
      float svFitMET_DM1_DOWN = -10;
      float svFitTransverseMass_DM1_DOWN = -10;

      float svFitMass_DM10_UP = -10;
      float svFitPt_DM10_UP = -10;
      float svFitEta_DM10_UP = -10;
      float svFitPhi_DM10_UP = -10;
      float svFitMET_DM10_UP = -10;
      float svFitTransverseMass_DM10_UP = -10;
      float svFitMass_DM10_DOWN = -10;
      float svFitPt_DM10_DOWN = -10;
      float svFitEta_DM10_DOWN = -10;
      float svFitPhi_DM10_DOWN = -10;
      float svFitMET_DM10_DOWN = -10;
      float svFitTransverseMass_DM10_DOWN = -10;

      float svFitMass_UncMet_UP = -10;
      float svFitPt_UncMet_UP = -10;
      float svFitEta_UncMet_UP = -10;
      float svFitPhi_UncMet_UP = -10;
      float svFitMET_UncMet_UP = -10;
      float svFitTransverseMass_UncMet_UP = -10;
      float svFitMass_UncMet_DOWN = -10;
      float svFitPt_UncMet_DOWN = -10;
      float svFitEta_UncMet_DOWN = -10;
      float svFitPhi_UncMet_DOWN = -10;
      float svFitMET_UncMet_DOWN = -10;
      float svFitTransverseMass_UncMet_DOWN = -10;

      float svFitMass_ClusteredMet_UP = -10;
      float svFitPt_ClusteredMet_UP = -10;
      float svFitEta_ClusteredMet_UP = -10;
      float svFitPhi_ClusteredMet_UP = -10;
      float svFitMET_ClusteredMet_UP = -10;
      float svFitTransverseMass_ClusteredMet_UP = -10;
      float svFitMass_ClusteredMet_DOWN = -10;
      float svFitPt_ClusteredMet_DOWN = -10;
      float svFitEta_ClusteredMet_DOWN = -10;
      float svFitPhi_ClusteredMet_DOWN = -10;
      float svFitMET_ClusteredMet_DOWN = -10;
      float svFitTransverseMass_ClusteredMet_DOWN = -10;
      
                                                                                  
      TBranch *newBranch11 = t->Branch("m_sv_UP", &svFitMass_UP, "m_sv_UP/F");
      TBranch *newBranch12 = t->Branch("pt_sv_UP", &svFitPt_UP, "pt_sv_UP/F");
      TBranch *newBranch13 = t->Branch("eta_sv_UP", &svFitEta_UP, "eta_sv_UP/F");
      TBranch *newBranch14 = t->Branch("phi_sv_UP", &svFitPhi_UP, "phi_sv_UP/F");
      TBranch *newBranch15 = t->Branch("met_sv_UP", &svFitMET_UP, "met_sv_UP/F");
      TBranch *newBranch16 = t->Branch("mt_sv_UP", &svFitTransverseMass_UP, "mt_sv_UP/F");

      TBranch *newBranch17 = t->Branch("m_sv_DOWN", &svFitMass_DOWN, "m_sv_DOWN/F");
      TBranch *newBranch18 = t->Branch("pt_sv_DOWN", &svFitPt_DOWN, "pt_sv_DOWN/F");
      TBranch *newBranch19 = t->Branch("eta_sv_DOWN", &svFitEta_DOWN, "eta_sv_DOWN/F");
      TBranch *newBranch20 = t->Branch("phi_sv_DOWN", &svFitPhi_DOWN, "phi_sv_DOWN/F");
      TBranch *newBranch21 = t->Branch("met_sv_DOWN", &svFitMET_DOWN, "met_sv_DOWN/F");
      TBranch *newBranch22 = t->Branch("mt_sv_DOWN", &svFitTransverseMass_DOWN, "mt_sv_DOWN/F");

      TBranch *newBranch23 = t->Branch("m_sv_DM0_UP", &svFitMass_DM0_UP, "m_sv_DM0_UP/F");
      TBranch *newBranch24 = t->Branch("pt_sv_DM0_UP", &svFitPt_DM0_UP, "pt_sv_DM0_UP/F");
      TBranch *newBranch25 = t->Branch("eta_sv_DM0_UP", &svFitEta_DM0_UP, "eta_sv_DM0_UP/F");
      TBranch *newBranch26 = t->Branch("phi_sv_DM0_UP", &svFitPhi_DM0_UP, "phi_sv_DM0_UP/F");
      TBranch *newBranch27 = t->Branch("met_sv_DM0_UP", &svFitMET_DM0_UP, "met_sv_DM0_UP/F");
      TBranch *newBranch28 = t->Branch("mt_sv_DM0_UP", &svFitTransverseMass_DM0_UP, "mt_sv_DM0_UP/F");

      TBranch *newBranch29 = t->Branch("m_sv_DM0_DOWN", &svFitMass_DM0_DOWN, "m_sv_DM0_DOWN/F");
      TBranch *newBranch30 = t->Branch("pt_sv_DM0_DOWN", &svFitPt_DM0_DOWN, "pt_sv_DM0_DOWN/F");
      TBranch *newBranch31 = t->Branch("eta_sv_DM0_DOWN", &svFitEta_DM0_DOWN, "eta_sv_DM0_DOWN/F");
      TBranch *newBranch32 = t->Branch("phi_sv_DM0_DOWN", &svFitPhi_DM0_DOWN, "phi_sv_DM0_DOWN/F");
      TBranch *newBranch33 = t->Branch("met_sv_DM0_DOWN", &svFitMET_DM0_DOWN, "met_sv_DM0_DOWN/F");
      TBranch *newBranch34 = t->Branch("mt_sv_DM0_DOWN", &svFitTransverseMass_DM0_DOWN, "mt_sv_DM0_DOWN/F");

      TBranch *newBranch35 = t->Branch("m_sv_DM1_UP", &svFitMass_DM1_UP, "m_sv_DM1_UP/F");
      TBranch *newBranch36 = t->Branch("pt_sv_DM1_UP", &svFitPt_DM1_UP, "pt_sv_DM1_UP/F");
      TBranch *newBranch37 = t->Branch("eta_sv_DM1_UP", &svFitEta_DM1_UP, "eta_sv_DM1_UP/F");
      TBranch *newBranch38 = t->Branch("phi_sv_DM1_UP", &svFitPhi_DM1_UP, "phi_sv_DM1_UP/F");
      TBranch *newBranch39 = t->Branch("met_sv_DM1_UP", &svFitMET_DM1_UP, "met_sv_DM1_UP/F");
      TBranch *newBranch40 = t->Branch("mt_sv_DM1_UP", &svFitTransverseMass_DM1_UP, "mt_sv_DM1_UP/F");

      TBranch *newBranch41 = t->Branch("m_sv_DM1_DOWN", &svFitMass_DM1_DOWN, "m_sv_DM1_DOWN/F");
      TBranch *newBranch42 = t->Branch("pt_sv_DM1_DOWN", &svFitPt_DM1_DOWN, "pt_sv_DM1_DOWN/F");
      TBranch *newBranch43 = t->Branch("eta_sv_DM1_DOWN", &svFitEta_DM1_DOWN, "eta_sv_DM1_DOWN/F");
      TBranch *newBranch44 = t->Branch("phi_sv_DM1_DOWN", &svFitPhi_DM1_DOWN, "phi_sv_DM1_DOWN/F");
      TBranch *newBranch45 = t->Branch("met_sv_DM1_DOWN", &svFitMET_DM1_DOWN, "met_sv_DM1_DOWN/F");
      TBranch *newBranch46 = t->Branch("mt_sv_DM1_DOWN", &svFitTransverseMass_DM1_DOWN, "mt_sv_DM1_DOWN/F");

      TBranch *newBranch47 = t->Branch("m_sv_DM10_UP", &svFitMass_DM10_UP, "m_sv_DM10_UP/F");
      TBranch *newBranch48 = t->Branch("pt_sv_DM10_UP", &svFitPt_DM10_UP, "pt_sv_DM10_UP/F");
      TBranch *newBranch49 = t->Branch("eta_sv_DM10_UP", &svFitEta_DM10_UP, "eta_sv_DM10_UP/F");
      TBranch *newBranch50 = t->Branch("phi_sv_DM10_UP", &svFitPhi_DM10_UP, "phi_sv_DM10_UP/F");
      TBranch *newBranch51 = t->Branch("met_sv_DM10_UP", &svFitMET_DM10_UP, "met_sv_DM10_UP/F");
      TBranch *newBranch52 = t->Branch("mt_sv_DM10_UP", &svFitTransverseMass_DM10_UP, "mt_sv_DM10_UP/F");

      TBranch *newBranch53 = t->Branch("m_sv_DM10_DOWN", &svFitMass_DM10_DOWN, "m_sv_DM10_DOWN/F");
      TBranch *newBranch54 = t->Branch("pt_sv_DM10_DOWN", &svFitPt_DM10_DOWN, "pt_sv_DM10_DOWN/F");
      TBranch *newBranch55 = t->Branch("eta_sv_DM10_DOWN", &svFitEta_DM10_DOWN, "eta_sv_DM10_DOWN/F");
      TBranch *newBranch56 = t->Branch("phi_sv_DM10_DOWN", &svFitPhi_DM10_DOWN, "phi_sv_DM10_DOWN/F");
      TBranch *newBranch57 = t->Branch("met_sv_DM10_DOWN", &svFitMET_DM10_DOWN, "met_sv_DM10_DOWN/F");
      TBranch *newBranch58 = t->Branch("mt_sv_DM10_DOWN", &svFitTransverseMass_DM10_DOWN, "mt_sv_DM10_DOWN/F");

      TBranch *newBranch59 = t->Branch("m_sv_UncMet_UP", &svFitMass_UncMet_UP, "m_sv_UncMet_UP/F");
      TBranch *newBranch60 = t->Branch("pt_sv_UncMet_UP", &svFitPt_UncMet_UP, "pt_sv_UncMet_UP/F");
      TBranch *newBranch61 = t->Branch("eta_sv_UncMet_UP", &svFitEta_UncMet_UP, "eta_sv_UncMet_UP/F");
      TBranch *newBranch62 = t->Branch("phi_sv_UncMet_UP", &svFitPhi_UncMet_UP, "phi_sv_UncMet_UP/F");
      TBranch *newBranch63 = t->Branch("met_sv_UncMet_UP", &svFitMET_UncMet_UP, "met_sv_UncMet_UP/F");
      TBranch *newBranch64 = t->Branch("mt_sv_UncMet_UP", &svFitTransverseMass_UncMet_UP, "mt_sv_UncMet_UP/F");

      TBranch *newBranch65 = t->Branch("m_sv_UncMet_DOWN", &svFitMass_UncMet_DOWN, "m_sv_UncMet_DOWN/F");
      TBranch *newBranch66 = t->Branch("pt_sv_UncMet_DOWN", &svFitPt_UncMet_DOWN, "pt_sv_UncMet_DOWN/F");
      TBranch *newBranch67 = t->Branch("eta_sv_UncMet_DOWN", &svFitEta_UncMet_DOWN, "eta_sv_UncMet_DOWN/F");
      TBranch *newBranch68 = t->Branch("phi_sv_UncMet_DOWN", &svFitPhi_UncMet_DOWN, "phi_sv_UncMet_DOWN/F");
      TBranch *newBranch69 = t->Branch("met_sv_UncMet_DOWN", &svFitMET_UncMet_DOWN, "met_sv_UncMet_DOWN/F");
      TBranch *newBranch70 = t->Branch("mt_sv_UncMet_DOWN", &svFitTransverseMass_UncMet_DOWN, "mt_sv_UncMet_DOWN/F");

      TBranch *newBranch71 = t->Branch("m_sv_ClusteredMet_UP", &svFitMass_ClusteredMet_UP, "m_sv_ClusteredMet_UP/F");
      TBranch *newBranch72 = t->Branch("pt_sv_ClusteredMet_UP", &svFitPt_ClusteredMet_UP, "pt_sv_ClusteredMet_UP/F");
      TBranch *newBranch73 = t->Branch("eta_sv_ClusteredMet_UP", &svFitEta_ClusteredMet_UP, "eta_sv_ClusteredMet_UP/F");
      TBranch *newBranch74 = t->Branch("phi_sv_ClusteredMet_UP", &svFitPhi_ClusteredMet_UP, "phi_sv_ClusteredMet_UP/F");
      TBranch *newBranch75 = t->Branch("met_sv_ClusteredMet_UP", &svFitMET_ClusteredMet_UP, "met_sv_ClusteredMet_UP/F");
      TBranch *newBranch76 = t->Branch("mt_sv_ClusteredMet_UP", &svFitTransverseMass_ClusteredMet_UP, "mt_sv_ClusteredMet_UP/F");

      TBranch *newBranch77 = t->Branch("m_sv_ClusteredMet_DOWN", &svFitMass_ClusteredMet_DOWN, "m_sv_ClusteredMet_DOWN/F");
      TBranch *newBranch78 = t->Branch("pt_sv_ClusteredMet_DOWN", &svFitPt_ClusteredMet_DOWN, "pt_sv_ClusteredMet_DOWN/F");
      TBranch *newBranch79 = t->Branch("eta_sv_ClusteredMet_DOWN", &svFitEta_ClusteredMet_DOWN, "eta_sv_ClusteredMet_DOWN/F");
      TBranch *newBranch80 = t->Branch("phi_sv_ClusteredMet_DOWN", &svFitPhi_ClusteredMet_DOWN, "phi_sv_ClusteredMet_DOWN/F");
      TBranch *newBranch81 = t->Branch("met_sv_ClusteredMet_DOWN", &svFitMET_ClusteredMet_DOWN, "met_sv_ClusteredMet_DOWN/F");
      TBranch *newBranch82 = t->Branch("mt_sv_ClusteredMet_DOWN", &svFitTransverseMass_ClusteredMet_DOWN, "mt_sv_ClusteredMet_DOWN/F");
    
      TBranch *newBranch83 = t->Branch("metcorClusteredDown",    &metcorClusteredDown,   "metcorClusteredDown/F");
      TBranch *newBranch84 = t->Branch("metcorphiClusteredDown", &metcorphiClusteredDown,"metcorphiClusteredDown/F");
      TBranch *newBranch85 = t->Branch("metcorClusteredUp",      &metcorClusteredUp,     "metcorClusteredUp/F");
      TBranch *newBranch86 = t->Branch("metcorphiClusteredUp",   &metcorphiClusteredUp,  "metcorphiClusteredUp/F");
      TBranch *newBranch87 = t->Branch("metcorUncDown",          &metcorUncDown,         "metcorUncDown/F");
      TBranch *newBranch89 = t->Branch("metcorphiUncDown",       &metcorphiUncDown,      "metcorphiUncDown/F");
      TBranch *newBranch90 = t->Branch("metcorUncUp",            &metcorUncUp,           "metcorUncUp/F");
      TBranch *newBranch91 = t->Branch("metcorphiUncUp",         &metcorphiUncUp,        "metcorphiUncUp/F");

      TBranch *newBranch92 = t->Branch("m_sv_EES_UP", &svFitMassEES_UP, "m_sv_EES_UP/F");
      TBranch *newBranch93 = t->Branch("pt_sv_EES_UP", &svFitPtEES_UP, "pt_sv_EES_UP/F");
      TBranch *newBranch94 = t->Branch("eta_sv_EES_UP", &svFitEtaEES_UP, "eta_sv_EES_UP/F");
      TBranch *newBranch95 = t->Branch("phi_sv_EES_UP", &svFitPhiEES_UP, "phi_sv_EES_UP/F");
      TBranch *newBranch96 = t->Branch("met_sv_EES_UP", &svFitMETEES_UP, "met_sv_EES_UP/F");
      TBranch *newBranch97 = t->Branch("mt_sv_EES_UP", &svFitTransverseMassEES_UP, "mt_sv_EES_UP/F");

      TBranch *newBranch98 = t->Branch("m_sv_EES_DOWN", &svFitMassEES_DOWN, "m_sv_EES_DOWN/F");
      TBranch *newBranch99 = t->Branch("pt_sv_EES_DOWN", &svFitPtEES_DOWN, "pt_sv_EES_DOWN/F");
      TBranch *newBranch100 = t->Branch("eta_sv_EES_DOWN", &svFitEtaEES_DOWN, "eta_sv_EES_DOWN/F");
      TBranch *newBranch101 = t->Branch("phi_sv_EES_DOWN", &svFitPhiEES_DOWN, "phi_sv_EES_DOWN/F");
      TBranch *newBranch102 = t->Branch("met_sv_EES_DOWN", &svFitMETEES_DOWN, "met_sv_EES_DOWN/F");
      TBranch *newBranch103 = t->Branch("mt_sv_EES_DOWN", &svFitTransverseMassEES_DOWN, "mt_sv_EES_DOWN/F");

    
    
      unsigned long long evt;
      int run, lumi;
      float pt1;
      float eta1;
      float phi1;
      float m1;
      float gen_match_1;
      float pt2;
      float eta2;
      float phi2;
      float m2;
      float gen_match_2;
      float decayMode=-999.;
      float decayMode2;
      //float mvaCovMatrix00;
      //float mvaCovMatrix10;
      //float mvaCovMatrix01;
      //float mvaCovMatrix11;
      float pfCovMatrix00;
      float pfCovMatrix10;
      float pfCovMatrix01;
      float pfCovMatrix11;
      //float mvamet_ex, // uncorrected mva met px (float)
      //  mvamet_ey, // uncorrected mva met py (float)
      float  genPx=-999.    , // generator Z/W/Higgs px (float)
        genPy =-999.   , // generator Z/W/Higgs py (float)
        visPx =-999.   , // generator visible Z/W/Higgs px (float)
        visPy =-999.   ; // generator visible Z/W/Higgs py (float)

      float njets =-999.   ;  // number of jets (hadronic jet multiplicity) (int)

      // define MET
      double measuredMETx = 0.;
      double measuredMETy = 0.;
      //float mvamet;
      //float mvametphi;
      float pfmet;
      float pfmetphi;
      TLorentzVector TMet(0,0,0,0);
      // define MET covariance
      TMatrixD covMET(2, 2);

      // MET Uncertainties
      float uncMetPtUp;
      float uncMetPtDown;
      float uncMetPhiUp;
      float uncMetPhiDown;
      float clusteredMetPtUp;
      float clusteredMetPtDown;
      float clusteredMetPhiUp;
      float clusteredMetPhiDown;
      double uncMetUpMETx = 0.;
      double uncMetUpMETy = 0.;
      double uncMetDownMETx = 0.;
      double uncMetDownMETy = 0.;
      double clusteredMetUpMETx = 0.;
      double clusteredMetUpMETy = 0.;
      double clusteredMetDownMETx = 0.;
      double clusteredMetDownMETy = 0.;
      
      std::basic_string<char> l1 = "x";
      std::basic_string<char> l2 = "x";
      // ZEE
      if (finalState==10) {l1 = "m"; l2 = "t";};
      if (finalState==11) {l1 = "e3"; l2 = "t";};
      if (finalState==12) {l1 = "e3"; l2 = "m";};
      if (finalState==13) {l1 = "t1"; l2 = "t2";};
      // ZMM
      if (finalState==20) {l1 = "m3"; l2 = "t";};
      if (finalState==21) {l1 = "e"; l2 = "t";};
      if (finalState==22) {l1 = "e"; l2 = "m3";};
      if (finalState==23) {l1 = "t1"; l2 = "t2";};

      t->SetBranchAddress("evt",&evt);
      t->SetBranchAddress("run",&run);
      t->SetBranchAddress("lumi",&lumi);
      //t->SetBranchAddress( (l1+"Pt").c_str(), &pt1);
      t->SetBranchAddress( "shiftedPt_3", &pt1); // Has Tau Energy Corrections alreayd applied
      t->SetBranchAddress( (l1+"Eta").c_str(), &eta1);
      t->SetBranchAddress( (l1+"Phi").c_str(), &phi1);
      t->SetBranchAddress( (l1+"Mass").c_str(), &m1);
      t->SetBranchAddress( (l1+"ZTTGenMatching").c_str(), &gen_match_1);
      t->SetBranchAddress( (l1+"DecayMode").c_str(), &decayMode);
      //t->SetBranchAddress( (l2+"Pt").c_str(), &pt2);
      t->SetBranchAddress( "shiftedPt_4", &pt2); // Has Tau Energy Corrections alreayd applied
      t->SetBranchAddress( (l2+"Eta").c_str(), &eta2);
      t->SetBranchAddress( (l2+"Phi").c_str(), &phi2);
      t->SetBranchAddress( (l2+"Mass").c_str(), &m2);
      t->SetBranchAddress( (l2+"ZTTGenMatching").c_str(), &gen_match_2);
      t->SetBranchAddress( (l2+"DecayMode").c_str(), &decayMode2);
      t->SetBranchAddress( "shiftedMET", &pfmet);
      t->SetBranchAddress( "shiftedMETPhi", &pfmetphi);

      // Use shifted METs from addPreSVFitCorrectionsZH.py
      t->SetBranchAddress( "shiftedUncMETUp", &uncMetPtUp);
      t->SetBranchAddress( "shiftedUncMETPhiUp", &uncMetPhiUp);
      t->SetBranchAddress( "shiftedUncMETDown", &uncMetPtDown);
      t->SetBranchAddress( "shiftedUncMETPhiDown", &uncMetPhiDown);
      t->SetBranchAddress( "shiftedClustMETUp", &clusteredMetPtUp);
      t->SetBranchAddress( "shiftedClustMETPhiUp", &clusteredMetPhiUp);
      t->SetBranchAddress( "shiftedClustMETDown", &clusteredMetPtDown);
      t->SetBranchAddress( "shiftedClustMETPhiDown", &clusteredMetPhiDown);

      //// Check if we have new met naming or old
      //std::basic_string<char> testString1 = "type1_pfMetPt_JetEnUp";
      //TObject* brNew = t->GetListOfBranches()->FindObject( testString1.c_str() );
      //std::cout << "testString: " << testString1 << " : " << brNew << std::endl;
      //std::basic_string<char> testString2 = "type1_pfMet_shiftedPt_JetEnUp";
      //TObject* brOld = t->GetListOfBranches()->FindObject( testString2.c_str() );
      //std::cout << "testString: " << testString2 << " : " << brOld << std::endl;
      //// New Files
      //if (brNew != 0) {
      //    t->SetBranchAddress( "type1_pfMetPt_UnclusteredEnUp", &uncMetPtUp);
      //    t->SetBranchAddress( "type1_pfMetPhi_UnclusteredEnUp", &uncMetPhiUp);
      //    t->SetBranchAddress( "type1_pfMetPt_UnclusteredEnDown", &uncMetPtDown);
      //    t->SetBranchAddress( "type1_pfMetPhi_UnclusteredEnDown", &uncMetPhiDown);
      //    t->SetBranchAddress( "type1_pfMetPt_JetEnUp", &clusteredMetPtUp);
      //    t->SetBranchAddress( "type1_pfMetPhi_JetEnUp", &clusteredMetPhiUp);
      //    t->SetBranchAddress( "type1_pfMetPt_JetEnDown", &clusteredMetPtDown);
      //    t->SetBranchAddress( "type1_pfMetPhi_JetEnDown", &clusteredMetPhiDown);
      //}
      //// Old Files - only use old if new isn't present
      //if (brOld != 0 && brNew == 0) {
      //    t->SetBranchAddress("type1_pfMet_shiftedPt_UnclusteredEnUp",&uncMetPtUp);
      //    t->SetBranchAddress("type1_pfMet_shiftedPhi_UnclusteredEnUp",&uncMetPhiUp);
      //    t->SetBranchAddress("type1_pfMet_shiftedPt_UnclusteredEnDown",&uncMetPtDown);
      //    t->SetBranchAddress("type1_pfMet_shiftedPhi_UnclusteredEnDown",&uncMetPhiDown);
      //    t->SetBranchAddress("type1_pfMet_shiftedPt_JetEnUp",&clusteredMetPtUp);
      //    t->SetBranchAddress("type1_pfMet_shiftedPhi_JetEnUp",&clusteredMetPhiUp);
      //    t->SetBranchAddress("type1_pfMet_shiftedPt_JetEnDown",&clusteredMetPtDown);
      //    t->SetBranchAddress("type1_pfMet_shiftedPhi_JetEnDown",&clusteredMetPhiDown);
      //}
      //t->SetBranchAddress("t1_t2_MvaMetCovMatrix00",&mvaCovMatrix00);
      //t->SetBranchAddress("t1_t2_MvaMetCovMatrix01",&mvaCovMatrix01);
      //t->SetBranchAddress("t1_t2_MvaMetCovMatrix10",&mvaCovMatrix10);
      //t->SetBranchAddress("t1_t2_MvaMetCovMatrix11",&mvaCovMatrix11);
      //t->SetBranchAddress("t1_t2_MvaMet",&mvamet);
      //t->SetBranchAddress("t1_t2_MvaMetPhi",&mvametphi);
      t->SetBranchAddress("jetVeto30", &njets);
      // Recoil variables below
      t->SetBranchAddress( "genpX", &genPx);
      t->SetBranchAddress( "genpY", &genPy);
      t->SetBranchAddress( "vispX", &visPx);
      t->SetBranchAddress( "vispY", &visPy);
      // FOR PF MET ANALYSIS
      t->SetBranchAddress("metcov00",&pfCovMatrix00);
      t->SetBranchAddress("metcov01",&pfCovMatrix01);
      t->SetBranchAddress("metcov10",&pfCovMatrix10);
      t->SetBranchAddress("metcov11",&pfCovMatrix11);

      // Recoil variables below
      t->SetBranchAddress( "genpX", &genPx);
      t->SetBranchAddress( "genpY", &genPy);
      t->SetBranchAddress( "vispX", &visPx);
      t->SetBranchAddress( "vispY", &visPy);
      // FOR PF MET ANALYSIS
      t->SetBranchAddress("metcov00",&pfCovMatrix00);
      t->SetBranchAddress("metcov01",&pfCovMatrix01);
      t->SetBranchAddress("metcov10",&pfCovMatrix10);
      t->SetBranchAddress("metcov11",&pfCovMatrix11);

      // use this RooT file when running on aMC@NLO DY and W+Jets MC samples
      RecoilCorrector* recoilCorrector = new RecoilCorrector(recoilFileName);
      if (metType == 1) std::cout<<"MetType: MvaMet"<<std::endl;
      if (metType == -1) std::cout<<"MetType: PF Met"<<std::endl;
      std::cout<<"recoiltype "<<recoilType<<" recoilFileName "<<recoilFileName<<std::endl;
      
      printf("Found tree -> weighting\n");
    
      double tesUP = 1.0 + tesSize;
      double tesDOWN = 1.0 - tesSize;

      int nevents = t->GetEntries();
      if ( parser.integerValue("numEvents") != -1 ) nevents = parser.integerValue("numEvents");
      for(Int_t i=0;i<nevents;++i){
         t->GetEntry(i);

         //Recoil Correction time
         // Correct WJets recoil for faked lepton / extra jet
         int recoilNJets;
         if(isWJets) {
            recoilNJets = njets + 1;
            std::cout << " - njets: " << njets << " recoilNJets: " << recoilNJets << std::endl;
         }
         else recoilNJets = njets;

         //// Using PF Met or Mva Met?
         //if (metType == 1) { // 1 = Mva Met
         //    std::cout << "Don't use MVA MET Right now!" << std::endl;
         //    TMet.SetPtEtaPhiM(mvamet,0,mvametphi,0);
         //    measuredMETx = mvamet*TMath::Cos(mvametphi);
         //    measuredMETy = mvamet*TMath::Sin(mvametphi);

         //    covMET[0][0] =  mvaCovMatrix00;
         //    covMET[1][0] =  mvaCovMatrix10;
         //    covMET[0][1] =  mvaCovMatrix01;
         //    covMET[1][1] =  mvaCovMatrix11;
         //} // mva met
         if (metType == -1) { // -1 = PF Met
             TMet.SetPtEtaPhiM(pfmet,0,pfmetphi,0);
             measuredMETx = pfmet*TMath::Cos(pfmetphi);
             measuredMETy = pfmet*TMath::Sin(pfmetphi);
             // Shifted METs
             uncMetUpMETx         = uncMetPtUp*TMath::Cos(uncMetPhiUp);
             uncMetUpMETy         = uncMetPtUp*TMath::Sin(uncMetPhiUp);
             uncMetDownMETx       = uncMetPtDown*TMath::Cos(uncMetPhiDown);
             uncMetDownMETy       = uncMetPtDown*TMath::Sin(uncMetPhiDown);
             clusteredMetUpMETx   = clusteredMetPtUp*TMath::Cos(clusteredMetPhiUp);
             clusteredMetUpMETy   = clusteredMetPtUp*TMath::Sin(clusteredMetPhiUp);
             clusteredMetDownMETx = clusteredMetPtDown*TMath::Cos(clusteredMetPhiDown);
             clusteredMetDownMETy = clusteredMetPtDown*TMath::Sin(clusteredMetPhiDown);

             covMET[0][0] =  pfCovMatrix00;
             covMET[1][0] =  pfCovMatrix10;
             covMET[0][1] =  pfCovMatrix01;
             covMET[1][1] =  pfCovMatrix11;
         } // pf met

         // Do recoil corrections if requested
         if(recoilType != 0){
            // Alexie shows that corrections via quantile mapping provide best results
            // for mva met and pf met, so
            // use that as the defaul.  People can switch if they want below
            //
            // RecoilCorrector::Correct == Quantile Mapping
            // RecoilCorrector::CorrectByMeanResolution == By Mean Res

            //recoilCorrector->Correct(
            recoilCorrector->CorrectByMeanResolution( // This method is faster (Alexei)
                              measuredMETx, // uncorrected mva met px (float)
                              measuredMETy, // uncorrected mva met py (float)
                              genPx, // generator Z/W/Higgs px (float)
                              genPy, // generator Z/W/Higgs py (float)
                              visPx, // generator visible Z/W/Higgs px (float)
                              visPy, // generator visible Z/W/Higgs py (float)
                              recoilNJets,  // number of jets (hadronic jet multiplicity) (int)
                              metcorr_ex, // corrected met px (float)
                              metcorr_ey  // corrected met py (float)
                              );
            // Shifted MET unc Up
            recoilCorrector->CorrectByMeanResolution( // This method is faster (Alexei)
                              uncMetUpMETx, // uncorrected mva met px (float)
                              uncMetUpMETy, // uncorrected mva met py (float)
                              genPx, // generator Z/W/Higgs px (float)
                              genPy, // generator Z/W/Higgs py (float)
                              visPx, // generator visible Z/W/Higgs px (float)
                              visPy, // generator visible Z/W/Higgs py (float)
                              recoilNJets,  // number of jets (hadronic jet multiplicity) (int)
                              metcorrUncUp_ex, // corrected met px (float)
                              metcorrUncUp_ey  // corrected met py (float)
                              );
            // Shifted MET unc Down
            recoilCorrector->CorrectByMeanResolution( // This method is faster (Alexei)
                              uncMetDownMETx, // uncorrected mva met px (float)
                              uncMetDownMETy, // uncorrected mva met py (float)
                              genPx, // generator Z/W/Higgs px (float)
                              genPy, // generator Z/W/Higgs py (float)
                              visPx, // generator visible Z/W/Higgs px (float)
                              visPy, // generator visible Z/W/Higgs py (float)
                              recoilNJets,  // number of jets (hadronic jet multiplicity) (int)
                              metcorrUncDown_ex, // corrected met px (float)
                              metcorrUncDown_ey  // corrected met py (float)
                              );
            // Shifted MET clustered Up
            recoilCorrector->CorrectByMeanResolution( // This method is faster (Alexei)
                              clusteredMetUpMETx, // uncorrected mva met px (float)
                              clusteredMetUpMETy, // uncorrected mva met py (float)
                              genPx, // generator Z/W/Higgs px (float)
                              genPy, // generator Z/W/Higgs py (float)
                              visPx, // generator visible Z/W/Higgs px (float)
                              visPy, // generator visible Z/W/Higgs py (float)
                              recoilNJets,  // number of jets (hadronic jet multiplicity) (int)
                              metcorrClusteredUp_ex, // corrected met px (float)
                              metcorrClusteredUp_ey  // corrected met py (float)
                              );
            // Shifted MET clustered Down
            recoilCorrector->CorrectByMeanResolution( // This method is faster (Alexei)
                              clusteredMetDownMETx, // uncorrected mva met px (float)
                              clusteredMetDownMETy, // uncorrected mva met py (float)
                              genPx, // generator Z/W/Higgs px (float)
                              genPy, // generator Z/W/Higgs py (float)
                              visPx, // generator visible Z/W/Higgs px (float)
                              visPy, // generator visible Z/W/Higgs py (float)
                              recoilNJets,  // number of jets (hadronic jet multiplicity) (int)
                              metcorrClusteredDown_ex, // corrected met px (float)
                              metcorrClusteredDown_ey  // corrected met py (float)
                              );
           std::cout << " - MEASURED:  met_ex: " << measuredMETx << "  met_ey: " << measuredMETy << std::endl;
           std::cout << " - CORRECTED: met_ex: " << metcorr_ex << "  met_ey: " << metcorr_ey << std::endl;
        }
        else{
           metcorr_ex = measuredMETx;
           metcorr_ey = measuredMETy;
           metcorrUncUp_ex = uncMetUpMETx;
           metcorrUncUp_ey = uncMetUpMETy;
           metcorrUncDown_ex = uncMetDownMETx;
           metcorrUncDown_ey = uncMetDownMETy;
           metcorrClusteredUp_ex = clusteredMetUpMETx;
           metcorrClusteredUp_ey = clusteredMetUpMETy;
           metcorrClusteredDown_ex = clusteredMetDownMETx;
           metcorrClusteredDown_ey = clusteredMetDownMETy;
        }

        metcor = TMath::Sqrt( metcorr_ex*metcorr_ex + metcorr_ey*metcorr_ey);
        metcorphi = TMath::ATan2( metcorr_ey, metcorr_ex );
        std::cout << " - metcor "<<metcor<<" metcorphi "<<metcorphi<<std::endl;

        std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;

      // LLMT only receives TES (not doing MES)
      if(channel=="mt"){

         mass2 = m2;
         measuredTauLeptons.push_back(
           classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1));
        
        measuredTauLeptons.push_back(
           classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2, decayMode2)); 
        std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1 << " mass1 " << mass1 << " pt2: "<< pt2<< " mass2: "<< mass2 <<std::endl;        
        //modified
        runSVFit(measuredTauLeptons, metcorr_ex, metcorr_ey, covMET, 0, svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass);
        std::cout<<"finished runningSVFit LLMT: "<< svFitMass <<std::endl;





        if(doES) {

          // Fill all with nominal values
          svFitMass_UP=svFitMass;
          svFitPt_UP=svFitPt;
          svFitEta_UP=svFitEta;
          svFitPhi_UP=svFitPhi;
          svFitMET_UP=svFitMET;
          svFitTransverseMass_UP=svFitTransverseMass;
          svFitMass_DOWN=svFitMass;
          svFitPt_DOWN=svFitPt;
          svFitEta_DOWN=svFitEta;
          svFitPhi_DOWN=svFitPhi;
          svFitMET_DOWN=svFitMET;
          svFitTransverseMass_DOWN=svFitTransverseMass;

          svFitMass_DM0_UP=svFitMass;
          svFitPt_DM0_UP=svFitPt;
          svFitEta_DM0_UP=svFitEta;
          svFitPhi_DM0_UP=svFitPhi;
          svFitMET_DM0_UP=svFitMET;
          svFitTransverseMass_DM0_UP=svFitTransverseMass;
          svFitMass_DM0_DOWN=svFitMass;
          svFitPt_DM0_DOWN=svFitPt;
          svFitEta_DM0_DOWN=svFitEta;
          svFitPhi_DM0_DOWN=svFitPhi;
          svFitMET_DM0_DOWN=svFitMET;
          svFitTransverseMass_DM0_DOWN=svFitTransverseMass;

          svFitMass_DM1_UP=svFitMass;
          svFitPt_DM1_UP=svFitPt;
          svFitEta_DM1_UP=svFitEta;
          svFitPhi_DM1_UP=svFitPhi;
          svFitMET_DM1_UP=svFitMET;
          svFitTransverseMass_DM1_UP=svFitTransverseMass;
          svFitMass_DM1_DOWN=svFitMass;
          svFitPt_DM1_DOWN=svFitPt;
          svFitEta_DM1_DOWN=svFitEta;
          svFitPhi_DM1_DOWN=svFitPhi;
          svFitMET_DM1_DOWN=svFitMET;
          svFitTransverseMass_DM1_DOWN=svFitTransverseMass;

          svFitMass_DM10_UP=svFitMass;
          svFitPt_DM10_UP=svFitPt;
          svFitEta_DM10_UP=svFitEta;
          svFitPhi_DM10_UP=svFitPhi;
          svFitMET_DM10_UP=svFitMET;
          svFitTransverseMass_DM10_UP=svFitTransverseMass;
          svFitMass_DM10_DOWN=svFitMass;
          svFitPt_DM10_DOWN=svFitPt;
          svFitEta_DM10_DOWN=svFitEta;
          svFitPhi_DM10_DOWN=svFitPhi;
          svFitMET_DM10_DOWN=svFitMET;
          svFitTransverseMass_DM10_DOWN=svFitTransverseMass;

          // Sort the DM based shifts below
          // Shift Tau Up
          if (gen_match_2==5){
             float  ES_UP_scale=tesUP; // for real taus
             double pt2_UP = pt2 * ES_UP_scale;
             double metcorr_ex_UP, metcorr_ey_UP;
             double dx2_UP, dy2_UP;
             dx2_UP = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_UP_scale ) - 1.);
             dy2_UP = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_UP_scale ) - 1.);
             metcorr_ex_UP = metcorr_ex + dx2_UP;
             metcorr_ey_UP = metcorr_ey + dy2_UP;
             
             std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUP;
             measuredTauLeptonsUP.push_back(
              classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1));
             measuredTauLeptonsUP.push_back(
              classic_svFit::MeasuredTauLepton(decayType2,  pt2_UP, eta2, phi2,  mass2, decayMode2));

             runSVFit(measuredTauLeptonsUP, metcorr_ex_UP, metcorr_ey_UP, covMET, 0, svFitMass_UP, svFitPt_UP, svFitEta_UP, svFitPhi_UP, svFitMET_UP, svFitTransverseMass_UP);

             // Shift Tau Down
             float ES_DOWN_scale=tesDOWN; // tau
             double pt2_DOWN;
             pt2_DOWN = pt2 * ES_DOWN_scale;
             double metcorr_ex_DOWN, metcorr_ey_DOWN;
             double dx2_DOWN, dy2_DOWN;
             dx2_DOWN = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_DOWN_scale ) - 1.);
             dy2_DOWN = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_DOWN_scale ) - 1.);
             metcorr_ex_DOWN = metcorr_ex + dx2_DOWN;
             metcorr_ey_DOWN = metcorr_ey + dy2_DOWN;

             std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDOWN;          
             measuredTauLeptonsDOWN.push_back(
              classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1));
             measuredTauLeptonsDOWN.push_back(
              classic_svFit::MeasuredTauLepton(decayType2,  pt2_DOWN, eta2, phi2,  mass2, decayMode2));

             runSVFit(measuredTauLeptonsDOWN, metcorr_ex_DOWN, metcorr_ey_DOWN, covMET, 0, svFitMass_DOWN, svFitPt_DOWN, svFitEta_DOWN, svFitPhi_DOWN, svFitMET_DOWN, svFitTransverseMass_DOWN);

             if (decayMode2 == 0) {
                svFitMass_DM0_UP=svFitMass_UP;
                svFitPt_DM0_UP=svFitPt_UP;
                svFitEta_DM0_UP=svFitEta_UP;
                svFitPhi_DM0_UP=svFitPhi_UP;
                svFitMET_DM0_UP=svFitMET_UP;
                svFitTransverseMass_DM0_UP=svFitTransverseMass_UP;
                svFitMass_DM0_DOWN=svFitMass_DOWN;
                svFitPt_DM0_DOWN=svFitPt_DOWN;
                svFitEta_DM0_DOWN=svFitEta_DOWN;
                svFitPhi_DM0_DOWN=svFitPhi_DOWN;
                svFitMET_DM0_DOWN=svFitMET_DOWN;
                svFitTransverseMass_DM0_DOWN=svFitTransverseMass_DOWN;
             }
             if (decayMode2 == 1) {
                svFitMass_DM1_UP=svFitMass_UP;
                svFitPt_DM1_UP=svFitPt_UP;
                svFitEta_DM1_UP=svFitEta_UP;
                svFitPhi_DM1_UP=svFitPhi_UP;
                svFitMET_DM1_UP=svFitMET_UP;
                svFitTransverseMass_DM1_UP=svFitTransverseMass_UP;
                svFitMass_DM1_DOWN=svFitMass_DOWN;
                svFitPt_DM1_DOWN=svFitPt_DOWN;
                svFitEta_DM1_DOWN=svFitEta_DOWN;
                svFitPhi_DM1_DOWN=svFitPhi_DOWN;
                svFitMET_DM1_DOWN=svFitMET_DOWN;
                svFitTransverseMass_DM1_DOWN=svFitTransverseMass_DOWN;
             }
             if (decayMode2 == 10) {
                svFitMass_DM10_UP=svFitMass_UP;
                svFitPt_DM10_UP=svFitPt_UP;
                svFitEta_DM10_UP=svFitEta_UP;
                svFitPhi_DM10_UP=svFitPhi_UP;
                svFitMET_DM10_UP=svFitMET_UP;
                svFitTransverseMass_DM10_UP=svFitTransverseMass_UP;
                svFitMass_DM10_DOWN=svFitMass_DOWN;
                svFitPt_DM10_DOWN=svFitPt_DOWN;
                svFitEta_DM10_DOWN=svFitEta_DOWN;
                svFitPhi_DM10_DOWN=svFitPhi_DOWN;
                svFitMET_DM10_DOWN=svFitMET_DOWN;
                svFitTransverseMass_DM10_DOWN=svFitTransverseMass_DOWN;
             }
          }
        } // end doES
      } // LLMT

      // For LLET do TES & EES
      else if(channel=="et"){

         mass2 = m2;
         measuredTauLeptons.push_back(
           classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1));
        
        measuredTauLeptons.push_back(
           classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2, decayMode2)); 
        std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1 << " mass1 " << mass1 << " pt2: "<< pt2<< " mass2: "<< mass2 <<std::endl;        
        //modified
        runSVFit(measuredTauLeptons, metcorr_ex, metcorr_ey, covMET, 0, svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass);
        std::cout<<"finished runningSVFit LLET: "<< svFitMass <<std::endl;



        if(doES) {

          // Fill all with nominal values
          svFitMass_UP=svFitMass;
          svFitPt_UP=svFitPt;
          svFitEta_UP=svFitEta;
          svFitPhi_UP=svFitPhi;
          svFitMET_UP=svFitMET;
          svFitTransverseMass_UP=svFitTransverseMass;
          svFitMass_DOWN=svFitMass;
          svFitPt_DOWN=svFitPt;
          svFitEta_DOWN=svFitEta;
          svFitPhi_DOWN=svFitPhi;
          svFitMET_DOWN=svFitMET;
          svFitTransverseMass_DOWN=svFitTransverseMass;

          svFitMass_DM0_UP=svFitMass;
          svFitPt_DM0_UP=svFitPt;
          svFitEta_DM0_UP=svFitEta;
          svFitPhi_DM0_UP=svFitPhi;
          svFitMET_DM0_UP=svFitMET;
          svFitTransverseMass_DM0_UP=svFitTransverseMass;
          svFitMass_DM0_DOWN=svFitMass;
          svFitPt_DM0_DOWN=svFitPt;
          svFitEta_DM0_DOWN=svFitEta;
          svFitPhi_DM0_DOWN=svFitPhi;
          svFitMET_DM0_DOWN=svFitMET;
          svFitTransverseMass_DM0_DOWN=svFitTransverseMass;

          svFitMass_DM1_UP=svFitMass;
          svFitPt_DM1_UP=svFitPt;
          svFitEta_DM1_UP=svFitEta;
          svFitPhi_DM1_UP=svFitPhi;
          svFitMET_DM1_UP=svFitMET;
          svFitTransverseMass_DM1_UP=svFitTransverseMass;
          svFitMass_DM1_DOWN=svFitMass;
          svFitPt_DM1_DOWN=svFitPt;
          svFitEta_DM1_DOWN=svFitEta;
          svFitPhi_DM1_DOWN=svFitPhi;
          svFitMET_DM1_DOWN=svFitMET;
          svFitTransverseMass_DM1_DOWN=svFitTransverseMass;

          svFitMass_DM10_UP=svFitMass;
          svFitPt_DM10_UP=svFitPt;
          svFitEta_DM10_UP=svFitEta;
          svFitPhi_DM10_UP=svFitPhi;
          svFitMET_DM10_UP=svFitMET;
          svFitTransverseMass_DM10_UP=svFitTransverseMass;
          svFitMass_DM10_DOWN=svFitMass;
          svFitPt_DM10_DOWN=svFitPt;
          svFitEta_DM10_DOWN=svFitEta;
          svFitPhi_DM10_DOWN=svFitPhi;
          svFitMET_DM10_DOWN=svFitMET;
          svFitTransverseMass_DM10_DOWN=svFitTransverseMass;

          // Shift Tau Up
          if (gen_match_2==5){
             float ES_UP_scale=1.0; // this value is for jet -> tau fakes
             if (gen_match_2<5) ES_UP_scale=1.015; // for gen matched ele/muon
             if (gen_match_2==5) ES_UP_scale=tesUP; // for real taus
             double pt2_UP;
             pt2_UP = pt2 * ES_UP_scale;
             double metcorr_ex_UP, metcorr_ey_UP;
             double dx2_UP, dy2_UP;
             dx2_UP = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_UP_scale ) - 1.);
             dy2_UP = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_UP_scale ) - 1.);
             metcorr_ex_UP = metcorr_ex + dx2_UP;
             metcorr_ey_UP = metcorr_ey + dy2_UP;
             
             std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUP;
             measuredTauLeptonsUP.push_back(
              classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1));
             measuredTauLeptonsUP.push_back(
              classic_svFit::MeasuredTauLepton(decayType2,  pt2_UP, eta2, phi2,  mass2, decayMode2));

             runSVFit(measuredTauLeptonsUP, metcorr_ex_UP, metcorr_ey_UP, covMET, 0, svFitMass_UP, svFitPt_UP, svFitEta_UP, svFitPhi_UP, svFitMET_UP, svFitTransverseMass_UP);

             // Shift Tau Down
             float ES_DOWN_scale=1.0; // jet
             if (gen_match_2==5) ES_DOWN_scale=tesDOWN; // tau
             if (gen_match_2<5) ES_DOWN_scale=0.985;  // elec/mu
             double pt2_DOWN;
             pt2_DOWN = pt2 * ES_DOWN_scale;
             double metcorr_ex_DOWN, metcorr_ey_DOWN;
             double dx2_DOWN, dy2_DOWN;
             dx2_DOWN = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_DOWN_scale ) - 1.);
             dy2_DOWN = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_DOWN_scale ) - 1.);
             metcorr_ex_DOWN = metcorr_ex + dx2_DOWN;
             metcorr_ey_DOWN = metcorr_ey + dy2_DOWN;

             std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDOWN;          
             measuredTauLeptonsDOWN.push_back(
              classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1));
             measuredTauLeptonsDOWN.push_back(
              classic_svFit::MeasuredTauLepton(decayType2,  pt2_DOWN, eta2, phi2,  mass2, decayMode2));

             runSVFit(measuredTauLeptonsDOWN, metcorr_ex_DOWN, metcorr_ey_DOWN, covMET, 0, svFitMass_DOWN, svFitPt_DOWN, svFitEta_DOWN, svFitPhi_DOWN, svFitMET_DOWN, svFitTransverseMass_DOWN);

             if (decayMode2 == 0) {
                svFitMass_DM0_UP=svFitMass_UP;
                svFitPt_DM0_UP=svFitPt_UP;
                svFitEta_DM0_UP=svFitEta_UP;
                svFitPhi_DM0_UP=svFitPhi_UP;
                svFitMET_DM0_UP=svFitMET_UP;
                svFitTransverseMass_DM0_UP=svFitTransverseMass_UP;
                svFitMass_DM0_DOWN=svFitMass_DOWN;
                svFitPt_DM0_DOWN=svFitPt_DOWN;
                svFitEta_DM0_DOWN=svFitEta_DOWN;
                svFitPhi_DM0_DOWN=svFitPhi_DOWN;
                svFitMET_DM0_DOWN=svFitMET_DOWN;
                svFitTransverseMass_DM0_DOWN=svFitTransverseMass_DOWN;
             }
             if (decayMode2 == 1) {
                svFitMass_DM1_UP=svFitMass_UP;
                svFitPt_DM1_UP=svFitPt_UP;
                svFitEta_DM1_UP=svFitEta_UP;
                svFitPhi_DM1_UP=svFitPhi_UP;
                svFitMET_DM1_UP=svFitMET_UP;
                svFitTransverseMass_DM1_UP=svFitTransverseMass_UP;
                svFitMass_DM1_DOWN=svFitMass_DOWN;
                svFitPt_DM1_DOWN=svFitPt_DOWN;
                svFitEta_DM1_DOWN=svFitEta_DOWN;
                svFitPhi_DM1_DOWN=svFitPhi_DOWN;
                svFitMET_DM1_DOWN=svFitMET_DOWN;
                svFitTransverseMass_DM1_DOWN=svFitTransverseMass_DOWN;
             }
             if (decayMode2 == 10) {
                svFitMass_DM10_UP=svFitMass_UP;
                svFitPt_DM10_UP=svFitPt_UP;
                svFitEta_DM10_UP=svFitEta_UP;
                svFitPhi_DM10_UP=svFitPhi_UP;
                svFitMET_DM10_UP=svFitMET_UP;
                svFitTransverseMass_DM10_UP=svFitTransverseMass_UP;
                svFitMass_DM10_DOWN=svFitMass_DOWN;
                svFitPt_DM10_DOWN=svFitPt_DOWN;
                svFitEta_DM10_DOWN=svFitEta_DOWN;
                svFitPhi_DM10_DOWN=svFitPhi_DOWN;
                svFitMET_DM10_DOWN=svFitMET_DOWN;
                svFitTransverseMass_DM10_DOWN=svFitTransverseMass_DOWN;
             }
          }

          // Now EES
          // Only shift Electrons
          // 1% shift in barrel and 2.5% shift in endcap
          // applied to electrons in EMu and ETau channel
          float etaBarrelEndcap  = 1.479;
          float ES_UP_scale;
          if (fabs(eta1) < etaBarrelEndcap) ES_UP_scale = 1.01;
          else ES_UP_scale = 1.025;
          double pt1_UP = pt1 * ES_UP_scale;
          std::cout << "E eta: " << eta1 << " ees SF: " << ES_UP_scale << std::endl;
          double metcorr_ex_UP, metcorr_ey_UP;
          double dx1_UP, dy1_UP;
          dx1_UP = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_UP_scale ) - 1.);
          dy1_UP = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_UP_scale ) - 1.);
          metcorr_ex_UP = metcorr_ex + dx1_UP;
          metcorr_ey_UP = metcorr_ey + dy1_UP;
          std::cout << "px1 " << pt1 * TMath::Cos( phi1 ) << "  met px1 cor " << dx1_UP <<std::endl;
          std::cout << "py1 " << pt1 * TMath::Sin( phi1 ) << "  met py1 cor " << dy1_UP <<std::endl;
          std::cout << "pt1 " << pt1 << "  pt1_up " << pt1_UP <<std::endl;
          std::cout << "metcor_ex " << metcorr_ex << " ees: " << metcorr_ex_UP << std::endl;
          std::cout << "metcor_ey " << metcorr_ey << " ees: " << metcorr_ey_UP << std::endl;
          
          std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUP;
          measuredTauLeptonsUP.push_back(
           classic_svFit::MeasuredTauLepton(decayType1, pt1_UP, eta1,  phi1, mass1));
          measuredTauLeptonsUP.push_back(
           classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2, mass2));

          std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1_UP << " mass1 " << mass1 << " pt2: "<< pt2 << " mass2: "<< mass2 <<std::endl;        
          runSVFit(measuredTauLeptonsUP, metcorr_ex_UP, metcorr_ey_UP, covMET, 0, svFitMassEES_UP, svFitPtEES_UP, svFitEtaEES_UP, svFitPhiEES_UP, svFitMETEES_UP, svFitTransverseMassEES_UP);
          std::cout<<"finished runningSVFit in LLET EES Up: "<<svFitMassEES_UP<<std::endl;


          // EES Down
          float ES_DOWN_scale;
          if (fabs(eta1) < etaBarrelEndcap) ES_DOWN_scale = 0.99;
          else ES_DOWN_scale = 0.975;
          std::cout << "E eta: " << eta1 << " ees SF: " << ES_DOWN_scale << std::endl;
          double pt1_DOWN;
          pt1_DOWN = pt1 * ES_DOWN_scale;
          double metcorr_ex_DOWN, metcorr_ey_DOWN;
          double dx1_DOWN, dy1_DOWN;
          dx1_DOWN = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_DOWN_scale ) - 1.);
          dy1_DOWN = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_DOWN_scale ) - 1.);
          metcorr_ex_DOWN = metcorr_ex + dx1_DOWN;
          metcorr_ey_DOWN = metcorr_ey + dy1_DOWN;
          std::cout << "px1 " << pt1 * TMath::Cos( phi1 ) << "  met px1 cor " << dx1_DOWN <<std::endl;
          std::cout << "py1 " << pt1 * TMath::Sin( phi1 ) << "  met py1 cor " << dy1_DOWN <<std::endl;
          std::cout << "metcor_ex " << metcorr_ex << " ees: " << metcorr_ex_DOWN << std::endl;
          std::cout << "metcor_ey " << metcorr_ey << " ees: " << metcorr_ey_DOWN << std::endl;

          std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDOWN;
          measuredTauLeptonsDOWN.push_back(
           classic_svFit::MeasuredTauLepton(decayType1, pt1_DOWN, eta1,  phi1, mass1));
          measuredTauLeptonsDOWN.push_back(
           classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2));

          runSVFit(measuredTauLeptonsDOWN, metcorr_ex_DOWN, metcorr_ey_DOWN, covMET, 0, svFitMassEES_DOWN, svFitPtEES_DOWN, svFitEtaEES_DOWN, svFitPhiEES_DOWN, svFitMETEES_DOWN, svFitTransverseMassEES_DOWN);
          std::cout<<"finished runningSVFit in LLET EES Down: "<<svFitMassEES_DOWN<<std::endl;
        } // end doES
      } // LLET

      // For LLEM only do EES (no MES)
      else if(channel=="em"){
        // define lepton four vectors
        measuredTauLeptons.push_back(
         classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1));
        measuredTauLeptons.push_back(
         classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2, mass2));

        std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1 << " mass1 " << mass1 << " pt2: "<< pt2<< " mass2: "<< mass2 <<std::endl;        
        runSVFit(measuredTauLeptons, metcorr_ex, metcorr_ey, covMET, 0, svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass);
        std::cout<<"finished runningSVFit LLEM: "<< svFitMass <<std::endl;



        if(doES) {

          // Only shift Electrons
          // 1% shift in barrel and 2.5% shift in endcap
          // applied to electrons in emu channel
          float etaBarrelEndcap  = 1.479;
          float ES_UP_scale;
          if (fabs(eta1) < etaBarrelEndcap) ES_UP_scale = 1.01;
          else ES_UP_scale = 1.025;
          double pt1_UP = pt1 * ES_UP_scale;
          std::cout << "E eta: " << eta1 << " ees SF: " << ES_UP_scale << std::endl;
          double metcorr_ex_UP, metcorr_ey_UP;
          double dx1_UP, dy1_UP;
          dx1_UP = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_UP_scale ) - 1.);
          dy1_UP = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_UP_scale ) - 1.);
          metcorr_ex_UP = metcorr_ex + dx1_UP;
          metcorr_ey_UP = metcorr_ey + dy1_UP;
          std::cout << "px1 " << pt1 * TMath::Cos( phi1 ) << "  met px1 cor " << dx1_UP <<std::endl;
          std::cout << "py1 " << pt1 * TMath::Sin( phi1 ) << "  met py1 cor " << dy1_UP <<std::endl;
          std::cout << "pt1 " << pt1 << "  pt1_up " << pt1_UP <<std::endl;
          std::cout << "metcor_ex " << metcorr_ex << " ees: " << metcorr_ex_UP << std::endl;
          std::cout << "metcor_ey " << metcorr_ey << " ees: " << metcorr_ey_UP << std::endl;
          
          std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUP;
          measuredTauLeptonsUP.push_back(
           classic_svFit::MeasuredTauLepton(decayType1, pt1_UP, eta1,  phi1, mass1));
          measuredTauLeptonsUP.push_back(
           classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2, mass2));

          std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1_UP << " mass1 " << mass1 << " pt2: "<< pt2 << " mass2: "<< mass2 <<std::endl;        
          runSVFit(measuredTauLeptonsUP, metcorr_ex_UP, metcorr_ey_UP, covMET, 0, svFitMassEES_UP, svFitPtEES_UP, svFitEtaEES_UP, svFitPhiEES_UP, svFitMETEES_UP, svFitTransverseMassEES_UP);
          std::cout<<"finished runningSVFit in EMu EES Up"<<std::endl;


          // EES Down
          float ES_DOWN_scale;
          if (fabs(eta1) < etaBarrelEndcap) ES_DOWN_scale = 0.99;
          else ES_DOWN_scale = 0.975;
          std::cout << "E eta: " << eta1 << " ees SF: " << ES_DOWN_scale << std::endl;
          double pt1_DOWN;
          pt1_DOWN = pt1 * ES_DOWN_scale;
          double metcorr_ex_DOWN, metcorr_ey_DOWN;
          double dx1_DOWN, dy1_DOWN;
          dx1_DOWN = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_DOWN_scale ) - 1.);
          dy1_DOWN = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_DOWN_scale ) - 1.);
          metcorr_ex_DOWN = metcorr_ex + dx1_DOWN;
          metcorr_ey_DOWN = metcorr_ey + dy1_DOWN;
          std::cout << "px1 " << pt1 * TMath::Cos( phi1 ) << "  met px1 cor " << dx1_DOWN <<std::endl;
          std::cout << "py1 " << pt1 * TMath::Sin( phi1 ) << "  met py1 cor " << dy1_DOWN <<std::endl;
          std::cout << "metcor_ex " << metcorr_ex << " ees: " << metcorr_ex_DOWN << std::endl;
          std::cout << "metcor_ey " << metcorr_ey << " ees: " << metcorr_ey_DOWN << std::endl;

          std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDOWN;
          measuredTauLeptonsDOWN.push_back(
           classic_svFit::MeasuredTauLepton(decayType1, pt1_DOWN, eta1,  phi1, mass1));
          measuredTauLeptonsDOWN.push_back(
           classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2));

          runSVFit(measuredTauLeptonsDOWN, metcorr_ex_DOWN, metcorr_ey_DOWN, covMET, 0, svFitMassEES_DOWN, svFitPtEES_DOWN, svFitEtaEES_DOWN, svFitPhiEES_DOWN, svFitMETEES_DOWN, svFitTransverseMassEES_DOWN);


        } // end doES
      } // eMu



      else if(channel=="tt"){

         mass1 = m1;
         mass2 = m2;

        // Add Tau of higest Pt first
        if (pt1 > pt2) {
          measuredTauLeptons.push_back(
           classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1, decayMode));
        
          measuredTauLeptons.push_back(
           classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2, decayMode2));
        }
        else {
          measuredTauLeptons.push_back(
           classic_svFit::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2, decayMode2));

          measuredTauLeptons.push_back(
           classic_svFit::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1, decayMode));
        }

        std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1 << " mass1 " << mass1 << " pt2: "<< pt2<< " mass2: "<< mass2 <<std::endl;        
        runSVFit(measuredTauLeptons, metcorr_ex, metcorr_ey, covMET, 0, svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass);
        std::cout<<"finished runningSVFit LLTT: "<< svFitMass <<std::endl;


        if(doES) {

          // Fill all with nominal values
          svFitMass_UP=svFitMass;
          svFitPt_UP=svFitPt;
          svFitEta_UP=svFitEta;
          svFitPhi_UP=svFitPhi;
          svFitMET_UP=svFitMET;
          svFitTransverseMass_UP=svFitTransverseMass;
          svFitMass_DOWN=svFitMass;
          svFitPt_DOWN=svFitPt;
          svFitEta_DOWN=svFitEta;
          svFitPhi_DOWN=svFitPhi;
          svFitMET_DOWN=svFitMET;
          svFitTransverseMass_DOWN=svFitTransverseMass;

          svFitMass_DM0_UP=svFitMass;
          svFitPt_DM0_UP=svFitPt;
          svFitEta_DM0_UP=svFitEta;
          svFitPhi_DM0_UP=svFitPhi;
          svFitMET_DM0_UP=svFitMET;
          svFitTransverseMass_DM0_UP=svFitTransverseMass;
          svFitMass_DM0_DOWN=svFitMass;
          svFitPt_DM0_DOWN=svFitPt;
          svFitEta_DM0_DOWN=svFitEta;
          svFitPhi_DM0_DOWN=svFitPhi;
          svFitMET_DM0_DOWN=svFitMET;
          svFitTransverseMass_DM0_DOWN=svFitTransverseMass;

          svFitMass_DM1_UP=svFitMass;
          svFitPt_DM1_UP=svFitPt;
          svFitEta_DM1_UP=svFitEta;
          svFitPhi_DM1_UP=svFitPhi;
          svFitMET_DM1_UP=svFitMET;
          svFitTransverseMass_DM1_UP=svFitTransverseMass;
          svFitMass_DM1_DOWN=svFitMass;
          svFitPt_DM1_DOWN=svFitPt;
          svFitEta_DM1_DOWN=svFitEta;
          svFitPhi_DM1_DOWN=svFitPhi;
          svFitMET_DM1_DOWN=svFitMET;
          svFitTransverseMass_DM1_DOWN=svFitTransverseMass;

          svFitMass_DM10_UP=svFitMass;
          svFitPt_DM10_UP=svFitPt;
          svFitEta_DM10_UP=svFitEta;
          svFitPhi_DM10_UP=svFitPhi;
          svFitMET_DM10_UP=svFitMET;
          svFitTransverseMass_DM10_UP=svFitTransverseMass;
          svFitMass_DM10_DOWN=svFitMass;
          svFitPt_DM10_DOWN=svFitPt;
          svFitEta_DM10_DOWN=svFitEta;
          svFitPhi_DM10_DOWN=svFitPhi;
          svFitMET_DM10_DOWN=svFitMET;
          svFitTransverseMass_DM10_DOWN=svFitTransverseMass;


          //***************************************************************************
          //********************** Tau DM0 shifted up *********************************
          //***************************************************************************

          if ((gen_match_2==5 && decayMode2==0) or (gen_match_1==5 && decayMode==0)){
             std::cout << "DM0 UP    ---  ";
             float ES_UP_scale1 = 1.0;
             float ES_UP_scale2 = 1.0;
             if(gen_match_1==5 && decayMode==0) ES_UP_scale1 = tesUP;
             if(gen_match_2==5 && decayMode2==0) ES_UP_scale2 = tesUP;
             double pt1_UP, pt2_UP;
             pt1_UP = pt1 * ES_UP_scale1;
             pt2_UP = pt2 * ES_UP_scale2;
             double metcorr_ex_UP, metcorr_ey_UP;
             double dx1_UP, dy1_UP, dx2_UP, dy2_UP;
             dx1_UP = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_UP_scale1 ) - 1.);
             dy1_UP = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_UP_scale1 ) - 1.);
             dx2_UP = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_UP_scale2 ) - 1.);
             dy2_UP = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_UP_scale2 ) - 1.);
             metcorr_ex_UP = metcorr_ex + dx1_UP + dx2_UP;
             metcorr_ey_UP = metcorr_ey + dy1_UP + dy2_UP;
          
             std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUP;
             // Add Tau of higest Pt first
             if (pt1 > pt2) {
               measuredTauLeptonsUP.push_back(
                classic_svFit::MeasuredTauLepton(decayType1, pt1_UP, eta1,  phi1, mass1, decayMode));
               measuredTauLeptonsUP.push_back(
                classic_svFit::MeasuredTauLepton(decayType2,  pt2_UP, eta2, phi2,  mass2, decayMode2));
             }
             else {
               measuredTauLeptonsUP.push_back(
                classic_svFit::MeasuredTauLepton(decayType2,  pt2_UP, eta2, phi2,  mass2, decayMode2));
               measuredTauLeptonsUP.push_back(
                classic_svFit::MeasuredTauLepton(decayType1, pt1_UP, eta1,  phi1, mass1, decayMode));
             }

             runSVFit(measuredTauLeptonsUP, metcorr_ex_UP, metcorr_ey_UP, covMET, 0, svFitMass_DM0_UP, svFitPt_DM0_UP, svFitEta_DM0_UP, svFitPhi_DM0_UP, svFitMET_DM0_UP, svFitTransverseMass_DM0_UP);
          }

          //***************************************************************************
          //********************** Tau DM1 shifted up *********************************
          //***************************************************************************

          if ((decayMode==1 && gen_match_1==5) or (decayMode2==1 && gen_match_2==5)){
             std::cout << "DM1 UP    ---  ";
             float ES_UP_scale1 = 1.0;
             float ES_UP_scale2 = 1.0;
             if (decayMode==1 && gen_match_1==5) ES_UP_scale1 = tesUP;
             if (decayMode2==1 && gen_match_2==5) ES_UP_scale2 = tesUP;
             double pt1_UP, pt2_UP;
             pt1_UP = pt1 * ES_UP_scale1;
             pt2_UP = pt2 * ES_UP_scale2;
             double metcorr_ex_UP, metcorr_ey_UP;
             double dx1_UP, dy1_UP, dx2_UP, dy2_UP;
             dx1_UP = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_UP_scale1 ) - 1.);
             dy1_UP = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_UP_scale1 ) - 1.);
             dx2_UP = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_UP_scale2 ) - 1.);
             dy2_UP = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_UP_scale2 ) - 1.);
             metcorr_ex_UP = metcorr_ex + dx1_UP + dx2_UP;
             metcorr_ey_UP = metcorr_ey + dy1_UP + dy2_UP;
             std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUP;
             // Add Tau of higest Pt first
             if (pt1 > pt2) {
               measuredTauLeptonsUP.push_back(
                classic_svFit::MeasuredTauLepton(decayType1, pt1_UP, eta1,  phi1, mass1, decayMode));
               measuredTauLeptonsUP.push_back(
                classic_svFit::MeasuredTauLepton(decayType2,  pt2_UP, eta2, phi2,  mass2, decayMode2));
             }
             else {
               measuredTauLeptonsUP.push_back(
                classic_svFit::MeasuredTauLepton(decayType2,  pt2_UP, eta2, phi2,  mass2, decayMode2));
               measuredTauLeptonsUP.push_back(
                classic_svFit::MeasuredTauLepton(decayType1, pt1_UP, eta1,  phi1, mass1, decayMode));
             }

             runSVFit(measuredTauLeptonsUP, metcorr_ex_UP, metcorr_ey_UP, covMET, 0, svFitMass_DM1_UP, svFitPt_DM1_UP, svFitEta_DM1_UP, svFitPhi_DM1_UP, svFitMET_DM1_UP, svFitTransverseMass_DM1_UP);
          }

          //***************************************************************************
          //********************* Tau DM10 shifted up *********************************
          //***************************************************************************

          if ((decayMode2==10 && gen_match_2==5) or (decayMode==10 && gen_match_1==5)){
             std::cout << "DM10 UP    ---  ";
             float ES_UP_scale1 = 1.0;
             float ES_UP_scale2 = 1.0;
             if(decayMode==10 && gen_match_1==5) ES_UP_scale1 = tesUP;
             if(decayMode2==10 && gen_match_2==5) ES_UP_scale2 = tesUP;
             double pt1_UP, pt2_UP;
             pt1_UP = pt1 * ES_UP_scale1;
             pt2_UP = pt2 * ES_UP_scale2;
             double metcorr_ex_UP, metcorr_ey_UP;
             double dx1_UP, dy1_UP, dx2_UP, dy2_UP;
             dx1_UP = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_UP_scale1 ) - 1.);
             dy1_UP = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_UP_scale1 ) - 1.);
             dx2_UP = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_UP_scale2 ) - 1.);
             dy2_UP = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_UP_scale2 ) - 1.);
             metcorr_ex_UP = metcorr_ex + dx1_UP + dx2_UP;
             metcorr_ey_UP = metcorr_ey + dy1_UP + dy2_UP;
          
             std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUP;
          
             // Add Tau of higest Pt first
             if (pt1 > pt2) {
               measuredTauLeptonsUP.push_back(
                classic_svFit::MeasuredTauLepton(decayType1, pt1_UP, eta1,  phi1, mass1, decayMode));
               measuredTauLeptonsUP.push_back(
                classic_svFit::MeasuredTauLepton(decayType2,  pt2_UP, eta2, phi2,  mass2, decayMode2));
             }
             else {
               measuredTauLeptonsUP.push_back(
                classic_svFit::MeasuredTauLepton(decayType2,  pt2_UP, eta2, phi2,  mass2, decayMode2));
               measuredTauLeptonsUP.push_back(
                classic_svFit::MeasuredTauLepton(decayType1, pt1_UP, eta1,  phi1, mass1, decayMode));
             }

             runSVFit(measuredTauLeptonsUP, metcorr_ex_UP, metcorr_ey_UP, covMET, 0, svFitMass_DM10_UP, svFitPt_DM10_UP, svFitEta_DM10_UP, svFitPhi_DM10_UP, svFitMET_DM10_UP, svFitTransverseMass_DM10_UP);
          }


          //*****************************************************
          //************* Tau DM0 shifted down  *****************
          //*****************************************************

          if ((decayMode==0 && gen_match_1==5) or (decayMode2==0 && gen_match_2==5)){
             std::cout << "DM0 DOWN  ---  ";
             float ES_DOWN_scale1 = 1.0;
             float ES_DOWN_scale2 = 1.0;
             if (decayMode==0 && gen_match_1==5) ES_DOWN_scale1 = tesDOWN;
             if (decayMode2==0 && gen_match_2==5) ES_DOWN_scale2 = tesDOWN;
             double pt1_DOWN, pt2_DOWN;
             pt1_DOWN = pt1 * ES_DOWN_scale1;
             pt2_DOWN = pt2 * ES_DOWN_scale2;
             double metcorr_ex_DOWN, metcorr_ey_DOWN;
             double dx1_DOWN, dy1_DOWN, dx2_DOWN, dy2_DOWN;
             dx1_DOWN = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_DOWN_scale1 ) - 1.);
             dy1_DOWN = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_DOWN_scale1 ) - 1.);
             dx2_DOWN = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_DOWN_scale2 ) - 1.);
             dy2_DOWN = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_DOWN_scale2 ) - 1.);
             metcorr_ex_DOWN = metcorr_ex + dx1_DOWN + dx2_DOWN;
             metcorr_ey_DOWN = metcorr_ey + dy1_DOWN + dy2_DOWN;

             std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDOWN;
             if (pt1 > pt2) {
               measuredTauLeptonsDOWN.push_back(
                classic_svFit::MeasuredTauLepton(decayType1, pt1_DOWN, eta1,  phi1, mass1, decayMode));
               measuredTauLeptonsDOWN.push_back(
                classic_svFit::MeasuredTauLepton(decayType2,  pt2_DOWN, eta2, phi2,  mass2, decayMode2));
             }
             else {
               measuredTauLeptonsDOWN.push_back(
                classic_svFit::MeasuredTauLepton(decayType2,  pt2_DOWN, eta2, phi2,  mass2, decayMode2));
               measuredTauLeptonsDOWN.push_back(
                classic_svFit::MeasuredTauLepton(decayType1, pt1_DOWN, eta1,  phi1, mass1, decayMode));
             }

             runSVFit(measuredTauLeptonsDOWN, metcorr_ex_DOWN, metcorr_ey_DOWN, covMET, 0, svFitMass_DM0_DOWN, svFitPt_DM0_DOWN, svFitEta_DM0_DOWN, svFitPhi_DM0_DOWN, svFitMET_DM0_DOWN, svFitTransverseMass_DM0_DOWN);
          }

          //*****************************************************
          //************** Tau DM1 shifted down *****************
          //*****************************************************

          if ((decayMode==1 && gen_match_1==5) or (decayMode2==1 && gen_match_2==5)){
             std::cout << "DM1 DOWN  ---  ";
             float ES_DOWN_scale1 = 1.0;
             float ES_DOWN_scale2 = 1.0;
             if (decayMode==1 && gen_match_1==5) ES_DOWN_scale1 = tesDOWN;
             if (decayMode2==1 && gen_match_2==5) ES_DOWN_scale2 = tesDOWN;
             double pt1_DOWN, pt2_DOWN;
             pt1_DOWN = pt1 * ES_DOWN_scale1;
             pt2_DOWN = pt2 * ES_DOWN_scale2;
             double metcorr_ex_DOWN, metcorr_ey_DOWN;
             double dx1_DOWN, dy1_DOWN, dx2_DOWN, dy2_DOWN;
             dx1_DOWN = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_DOWN_scale1 ) - 1.);
             dy1_DOWN = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_DOWN_scale1 ) - 1.);
             dx2_DOWN = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_DOWN_scale2 ) - 1.);
             dy2_DOWN = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_DOWN_scale2 ) - 1.);
             metcorr_ex_DOWN = metcorr_ex + dx1_DOWN + dx2_DOWN;
             metcorr_ey_DOWN = metcorr_ey + dy1_DOWN + dy2_DOWN;

             std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDOWN;

             if (pt1 > pt2) {
               measuredTauLeptonsDOWN.push_back(
                classic_svFit::MeasuredTauLepton(decayType1, pt1_DOWN, eta1,  phi1, mass1, decayMode));
               measuredTauLeptonsDOWN.push_back(
                classic_svFit::MeasuredTauLepton(decayType2,  pt2_DOWN, eta2, phi2,  mass2, decayMode2));
             }
             else {
               measuredTauLeptonsDOWN.push_back(
                classic_svFit::MeasuredTauLepton(decayType2,  pt2_DOWN, eta2, phi2,  mass2, decayMode2));
               measuredTauLeptonsDOWN.push_back(
                classic_svFit::MeasuredTauLepton(decayType1, pt1_DOWN, eta1,  phi1, mass1, decayMode));
             }

             runSVFit(measuredTauLeptonsDOWN, metcorr_ex_DOWN, metcorr_ey_DOWN, covMET, 0, svFitMass_DM1_DOWN, svFitPt_DM1_DOWN, svFitEta_DM1_DOWN, svFitPhi_DM1_DOWN, svFitMET_DM1_DOWN, svFitTransverseMass_DM1_DOWN);
          }

          //*****************************************************
          //************* Tau DM10 shifted down *****************
          //*****************************************************

          if ((decayMode==10 && gen_match_1==5) or (decayMode2==10 && gen_match_2==5)){
             std::cout << "DM10 DOWN  ---  ";
             float ES_DOWN_scale1 = 1.0;
             float ES_DOWN_scale2 = 1.0;
             if (decayMode==10 && gen_match_1==5) ES_DOWN_scale1 = tesDOWN;
             if (decayMode2==10 && gen_match_2==5) ES_DOWN_scale2 = tesDOWN;
             double pt1_DOWN, pt2_DOWN;
             pt1_DOWN = pt1 * ES_DOWN_scale1;
             pt2_DOWN = pt2 * ES_DOWN_scale2;
             double metcorr_ex_DOWN, metcorr_ey_DOWN;
             double dx1_DOWN, dy1_DOWN, dx2_DOWN, dy2_DOWN;
             dx1_DOWN = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_DOWN_scale1 ) - 1.);
             dy1_DOWN = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_DOWN_scale1 ) - 1.);
             dx2_DOWN = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_DOWN_scale2 ) - 1.);
             dy2_DOWN = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_DOWN_scale2 ) - 1.);
             metcorr_ex_DOWN = metcorr_ex + dx1_DOWN + dx2_DOWN;
             metcorr_ey_DOWN = metcorr_ey + dy1_DOWN + dy2_DOWN;

             std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDOWN;
             if (pt1 > pt2) {
               measuredTauLeptonsDOWN.push_back(
                classic_svFit::MeasuredTauLepton(decayType1, pt1_DOWN, eta1,  phi1, mass1, decayMode));
               measuredTauLeptonsDOWN.push_back(
                classic_svFit::MeasuredTauLepton(decayType2,  pt2_DOWN, eta2, phi2,  mass2, decayMode2));
             }
             else {
               measuredTauLeptonsDOWN.push_back(
                classic_svFit::MeasuredTauLepton(decayType2,  pt2_DOWN, eta2, phi2,  mass2, decayMode2));
               measuredTauLeptonsDOWN.push_back(
                classic_svFit::MeasuredTauLepton(decayType1, pt1_DOWN, eta1,  phi1, mass1, decayMode));
             }

             runSVFit(measuredTauLeptonsDOWN, metcorr_ex_DOWN, metcorr_ey_DOWN, covMET, 0, svFitMass_DM10_DOWN, svFitPt_DM10_DOWN, svFitEta_DM10_DOWN, svFitPhi_DM10_DOWN, svFitMET_DM10_DOWN, svFitTransverseMass_DM10_DOWN);
          }
        }// Do ES (TT)

      } // Double Hadronic (TT)


      // Do the MET Uncertainties
      // The values for measuredTauLeptons have already been set
      if (doES) { 

        // Corrected MET values for saving
        // Will be re-corrected with TEC if running tautau channel
        metcorClusteredDown = TMath::Sqrt( metcorrClusteredDown_ex*metcorrClusteredDown_ex + metcorrClusteredDown_ey*metcorrClusteredDown_ey);
        metcorphiClusteredDown = TMath::ATan2( metcorrClusteredDown_ey, metcorrClusteredDown_ex );

        metcorClusteredUp = TMath::Sqrt( metcorrClusteredUp_ex*metcorrClusteredUp_ex + metcorrClusteredUp_ey*metcorrClusteredUp_ey);
        metcorphiClusteredUp = TMath::ATan2( metcorrClusteredUp_ey, metcorrClusteredUp_ex );

        metcorUncDown = TMath::Sqrt( metcorrUncDown_ex*metcorrUncDown_ex + metcorrUncDown_ey*metcorrUncDown_ey);
        metcorphiUncDown = TMath::ATan2( metcorrUncDown_ey, metcorrUncDown_ex );

        metcorUncUp = TMath::Sqrt( metcorrUncUp_ex*metcorrUncUp_ex + metcorrUncUp_ey*metcorrUncUp_ey);
        metcorphiUncUp = TMath::ATan2( metcorrUncUp_ey, metcorrUncUp_ex );

        runSVFit(measuredTauLeptons, metcorrUncUp_ex, metcorrUncUp_ey, covMET, 0, svFitMass_UncMet_UP, svFitPt_UncMet_UP, svFitEta_UncMet_UP, svFitPhi_UncMet_UP, svFitMET_UncMet_UP, svFitTransverseMass_UncMet_UP);
        runSVFit(measuredTauLeptons, metcorrUncDown_ex, metcorrUncDown_ey, covMET, 0, svFitMass_UncMet_DOWN, svFitPt_UncMet_DOWN, svFitEta_UncMet_DOWN, svFitPhi_UncMet_DOWN, svFitMET_UncMet_DOWN, svFitTransverseMass_UncMet_DOWN);
        runSVFit(measuredTauLeptons, metcorrClusteredUp_ex, metcorrClusteredUp_ey, covMET, 0, svFitMass_ClusteredMet_UP, svFitPt_ClusteredMet_UP, svFitEta_ClusteredMet_UP, svFitPhi_ClusteredMet_UP, svFitMET_ClusteredMet_UP, svFitTransverseMass_ClusteredMet_UP);
        runSVFit(measuredTauLeptons, metcorrClusteredDown_ex, metcorrClusteredDown_ey, covMET, 0, svFitMass_ClusteredMet_DOWN, svFitPt_ClusteredMet_DOWN, svFitEta_ClusteredMet_DOWN, svFitPhi_ClusteredMet_DOWN, svFitMET_ClusteredMet_DOWN, svFitTransverseMass_ClusteredMet_DOWN);
        std::cout << "MET Unclustered Energy Up   ---  "<< svFitMass_UncMet_UP << std::endl;
        std::cout << "MET Unclustered Energy Down ---  "<< svFitMass_UncMet_DOWN << std::endl;
        std::cout << "MET Clustered Energy Up     ---  "<< svFitMass_ClusteredMet_UP << std::endl;
        std::cout << "MET Clustered Energy Down   ---  "<< svFitMass_ClusteredMet_DOWN << std::endl;
        std::cout<< "Shifted MET Summary:\nmetcorr_ex " << metcorr_ex << "\n --- metcorrUncUp_ex " << metcorrUncUp_ex << " metcorrUncDown_ex " << metcorrUncDown_ex
        << " metcorrClusteredUp_ex " << metcorrClusteredUp_ex << " metcorrClusteredDown_ex " << metcorrClusteredDown_ex << std::endl;
        std::cout<< "metcorr_ey " << metcorr_ey << "\n --- metcorrUncUp_ey " << metcorrUncUp_ey << " metcorrUncDown_ey " << metcorrUncDown_ey
        << " metcorrClusteredUp_ey " << metcorrClusteredUp_ey << " metcorrClusteredDown_ey " << metcorrClusteredDown_ey << std::endl;

     }



      
      std::cout << "\n\n" << std::endl;
      //std::cout << "\n\nex: " << metcorr_ex << "   ey: " << metcorr_ey <<  " phi: " << metcorphi<<"\n"<<std::endl; 
      newBranch1->Fill();
      newBranch2->Fill();
      newBranch3->Fill();
      newBranch4->Fill();
      newBranch5->Fill();
      newBranch6->Fill();
      newBranch7->Fill();
      newBranch8->Fill();
      newBranch9->Fill();
      newBranch10->Fill();

      newBranch11->Fill();
      newBranch12->Fill();
      newBranch13->Fill();
      newBranch14->Fill();
      newBranch15->Fill();
      newBranch16->Fill();
      newBranch17->Fill();
      newBranch18->Fill();
      newBranch19->Fill();
      newBranch20->Fill();
      newBranch21->Fill();
      newBranch22->Fill();

      newBranch23->Fill();
      newBranch24->Fill();
      newBranch25->Fill();
      newBranch26->Fill();
      newBranch27->Fill();
      newBranch28->Fill();
      newBranch29->Fill();
      newBranch30->Fill();
      newBranch31->Fill();
      newBranch32->Fill();
      newBranch33->Fill();
      newBranch34->Fill();

      newBranch35->Fill();
      newBranch36->Fill();
      newBranch37->Fill();
      newBranch38->Fill();
      newBranch39->Fill();
      newBranch40->Fill();
      newBranch41->Fill();
      newBranch42->Fill();
      newBranch43->Fill();
      newBranch44->Fill();
      newBranch45->Fill();
      newBranch46->Fill();

      newBranch47->Fill();
      newBranch48->Fill();
      newBranch49->Fill();
      newBranch50->Fill();
      newBranch51->Fill();
      newBranch52->Fill();
      newBranch53->Fill();
      newBranch54->Fill();
      newBranch55->Fill();
      newBranch56->Fill();
      newBranch57->Fill();
      newBranch58->Fill();

      newBranch59->Fill();
      newBranch60->Fill();
      newBranch61->Fill();
      newBranch62->Fill();
      newBranch63->Fill();
      newBranch64->Fill();
      newBranch65->Fill();
      newBranch66->Fill();
      newBranch67->Fill();
      newBranch68->Fill();
      newBranch69->Fill();
      newBranch70->Fill();

      newBranch71->Fill();
      newBranch72->Fill();
      newBranch73->Fill();
      newBranch74->Fill();
      newBranch75->Fill();
      newBranch76->Fill();
      newBranch77->Fill();
      newBranch78->Fill();
      newBranch79->Fill();
      newBranch80->Fill();
      newBranch81->Fill();
      newBranch82->Fill();

      newBranch83->Fill();
      newBranch84->Fill();
      newBranch85->Fill();
      newBranch86->Fill();
      newBranch87->Fill();
      newBranch89->Fill();
      newBranch90->Fill();
      newBranch91->Fill();

      newBranch92->Fill();
      newBranch93->Fill();
      newBranch94->Fill();
      newBranch95->Fill();
      newBranch96->Fill();
      newBranch97->Fill();
      newBranch98->Fill();
      newBranch99->Fill();
      newBranch100->Fill();
      newBranch101->Fill();
      newBranch102->Fill();
      newBranch103->Fill();

      }
      dir->cd();
      t->Write("",TObject::kOverwrite);
      delete  t;
    } // if the iterator of the key is a TTree
  }
}

void runSVFit(std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons, 
          double measuredMETx, double measuredMETy, TMatrixD &covMET, float num, 
          float &svFitMass, float& svFitPt, float &svFitEta, float &svFitPhi, float &svFitMET, float &svFitTransverseMass){

  svfitAlgorithm.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
  if ( svfitAlgorithm.isValidSolution()) {

    svFitMass = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->getMass();
    svFitPt = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->getPt();
    svFitEta = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->getEta();
    svFitPhi = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->getPhi();
    svFitTransverseMass = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->getTransverseMass();

    const classic_svFit::LorentzVector ditau = static_cast<classic_svFit::HistogramAdapterDiTau*>(svfitAlgorithm.getHistogramAdapter())->getP4();
    svFitMET = (ditau - (measuredTauLeptons[0].p4() + measuredTauLeptons[1].p4())).Pt();

    //TLorentzVector testTau1, testTau2;
    //testTau1.SetPtEtaPhiM(measuredTauLeptons[0].p4().Pt(), measuredTauLeptons[0].p4().Eta(), measuredTauLeptons[0].p4().Phi(), measuredTauLeptons[0].p4().M());
    //testTau2.SetPtEtaPhiM(measuredTauLeptons[1].p4().Pt(), measuredTauLeptons[1].p4().Eta(), measuredTauLeptons[1].p4().Phi(), measuredTauLeptons[1].p4().M());

    //TLorentzVector ditau;
    //ditau.SetPtEtaPhiM(svFitPt, svFitEta, svFitPhi, svFitMass);
  }

}

//Thank you Renee Brun :)
void CopyDir(TDirectory *source, optutl::CommandLineParser parser) {
  //copy all objects and subdirs of directory source as a subdir of the current directory
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir; 
  if(source->GetName()!=parser.stringValue("inputFile")){
    adir = savdir->mkdir(source->GetName());
    std::cout<<"Source name is not outputfile name"<<std::endl;
    adir->cd();    
  }
  else{
    //adir = savdir->mkdir("input");
    adir->cd();    
  }

  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {
      source->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      adir->cd();
      CopyDir(subdir,parser);
      adir->cd();
    } else if (cl->InheritsFrom(TTree::Class())) {
      TTree *T = (TTree*)source->Get(key->GetName());
      adir->cd();
      TTree *newT = T->CloneTree(-1,"fast");
      newT->Write();
    } else {
      source->cd();
      TObject *obj = key->ReadObj();
      adir->cd();
      obj->Write();
      delete obj;
    }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}
void CopyFile(const char *fname, optutl::CommandLineParser parser) {
  //Copy all objects and subdirs of file fname as a subdir of the current directory
  TDirectory *target = gDirectory;
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    printf("Cannot copy file: %s\n",fname);
    target->cd();
    return;
  }
  target->cd();
  CopyDir(f,parser);
  delete f;
  target->cd();
}
void copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) 
{
  //prepare files to be copied
  if(gSystem->AccessPathName(parser.stringValue("inputFile").c_str())) {
    gSystem->CopyFile("hsimple.root", parser.stringValue("inputFile").c_str());
  }

  fNew->cd();
  CopyFile(parser.stringValue("inputFile").c_str(),parser);
  fNew->ls();
  fNew->Close();

}

