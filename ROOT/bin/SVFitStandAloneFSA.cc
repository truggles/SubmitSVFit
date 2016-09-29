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

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

//If recoilType 0 then don't do recoil
//              1 then aMC@NLO DY and W+Jets MC samples
//              2 MG5 DY and W+Jets MC samples or Higgs MC samples
//
//If doES      0 does not apply any ES shifts
//              1 applies ES shifts to TT channel, no effect on other channels
//
//If isWJets    0 no shift in number of jets used for recoil corrections
//              1 shifts njets + 1 for recoil corrections

void copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) ;
void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int recoilType, int doES, int isWJets) ;
void CopyFile(const char *fname, optutl::CommandLineParser parser);
void CopyDir(TDirectory *source,optutl::CommandLineParser parser);
void runSVFit(std::vector<svFitStandalone::MeasuredTauLepton> & measuredTauLeptons, TFile * inputFile_visPtResolution, double measuredMETx, double measuredMETy, TMatrixD &covMET, float num, float &svFitMass, float& svFitPt, float &svFitEta, float &svFitPhi, float &svFitMET, float &svFitTransverseMass);

int main (int argc, char* argv[]) 
{
   optutl::CommandLineParser parser ("Sets Event Weights in the ntuple");
   parser.addOption("branch",optutl::CommandLineParser::kString,"Branch","__svFit__");
   parser.addOption("newFile",optutl::CommandLineParser::kString,"newFile","newFile");
   parser.addOption("newOutputFile",optutl::CommandLineParser::kDouble,"New Output File",0.0);
   parser.addOption("recoilType",optutl::CommandLineParser::kDouble,"recoilType",0.0);
   parser.addOption("doES",optutl::CommandLineParser::kDouble,"doES",0.0);
   parser.addOption("isWJets",optutl::CommandLineParser::kDouble,"isWJets",0.0);

   parser.parseArguments (argc, argv);

   std::cout << "EXTRA COMMANDS:"
    << "\n --- recoilType: " << parser.doubleValue("recoilType")
    << "\n --- doES: " << parser.doubleValue("doES")
    << "\n --- isWJets: " << parser.doubleValue("isWJets") << std::endl;
   
   char TreeToUse[80]="first" ;

   TFile *fProduce;//= new TFile(parser.stringValue("newFile").c_str(),"UPDATE");

   if(parser.doubleValue("newOutputFile")>0){
   TFile *f = new TFile(parser.stringValue("outputFile").c_str(),"READ");
     std::cout<<"Creating new outputfile"<<std::endl;
     std::string newFileName = parser.stringValue("newFile");

     fProduce = new TFile(newFileName.c_str(),"RECREATE");
     copyFiles(parser, f, fProduce);//new TFile(parser.stringValue("outputFile").c_str()+"SVFit","UPDATE");
     fProduce = new TFile(newFileName.c_str(),"UPDATE");
     std::cout<<"listing the directories================="<<std::endl;
     fProduce->ls();
     readdir(fProduce,parser,TreeToUse,parser.doubleValue("recoilType"),parser.doubleValue("doES"),parser.doubleValue("isWJets"));
     fProduce->Close();
     f->Close();
   }
   else{
     TFile *f = new TFile(parser.stringValue("outputFile").c_str(),"UPDATE");
     readdir(f,parser,TreeToUse,parser.doubleValue("recoilType"),parser.doubleValue("doES"),parser.doubleValue("isWJets"));
     f->Close();
   }


} 


void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int recoilType, int doES, int isWJets) 
{
  std::string recoilFileName = "HTT-utilities/RecoilCorrections/data/recoilMvaMEt_76X_newTraining.root";
  if(recoilType == 1) //amc@nlo
    recoilFileName = "HTT-utilities/RecoilCorrections/data/recoilMvaMEt_76X_newTraining.root";
  if(recoilType == 2) //MG5
    recoilFileName = "HTT-utilities/RecoilCorrections/data/recoilMvaMEt_76X_newTraining_MG5.root";

  TDirectory *dirsav = gDirectory;
  TIter next(dir->GetListOfKeys());
  TKey *key;
  char stringA[80]="first";
  dir->cd();      
  while ((key = (TKey*)next())) {
    printf("Found key=%s \n",key->GetName());

    TObject *obj = key->ReadObj();
    if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
      dir->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      sprintf(TreeToUse,"%s",key->GetName());
      readdir(subdir,parser,TreeToUse,parser.doubleValue("recoilType"),parser.doubleValue("doES"),parser.doubleValue("isWJets"));
      dirsav->cd();
    }
    else if(obj->IsA()->InheritsFrom(TTree::Class())) {

      TTree *t = (TTree*)obj;
      float svFitMass = -10;
      float svFitPt = -10;
      float svFitEta = -10;
      float svFitPhi = -10;
      float svFitMET = -10;
      float svFitTransverseMass = -10;

      float mvametcorr_ex = -10; // corrected met px (float)
      float mvametcorr_ey = -10;  // corrected met py (float)
      float mvamet = -10; // corrected mvamet
      float mvametphi = -10; // corrected mvametphi

      //TBranch *newBranch = t->Branch(parser.stringValue("branch").c_str(),&svFitMass,(parser.stringValue("branch")+"/F").c_str());
      TBranch *newBranch1 = t->Branch("m_sv", &svFitMass, "m_sv/F");
      TBranch *newBranch2 = t->Branch("p_sv", &svFitPt, "pt_sv/F");
      TBranch *newBranch3 = t->Branch("eta_sv", &svFitEta, "eta_sv/F");
      TBranch *newBranch4 = t->Branch("phi_sv", &svFitPhi, "phi_sv/F");
      TBranch *newBranch5 = t->Branch("met_sv", &svFitMET, "met_sv/F");
      TBranch *newBranch6 = t->Branch("mt_sv", &svFitTransverseMass, "mt_sv/F");

      TBranch *newBranch7 = t->Branch("mvametcorr_ex", &mvametcorr_ex, "mvametcorr_ex/F");
      TBranch *newBranch8 = t->Branch("mvametcorr_ey", &mvametcorr_ey, "mvametcorr_ey/F");
      TBranch *newBranch9 = t->Branch("mvamet", &mvamet, "mvamet/F");
      TBranch *newBranch10 = t->Branch("mvametphi", &mvametphi, "mvametphi/F");

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
      //if(doES) {
      TBranch *newBranch11 = t->Branch("m_sv_UP", &svFitMass_UP, "m_sv_UP/F");
      TBranch *newBranch12 = t->Branch("p_sv_UP", &svFitPt_UP, "pt_sv_UP/F");
      TBranch *newBranch13 = t->Branch("eta_sv_UP", &svFitEta_UP, "eta_sv_UP/F");
      TBranch *newBranch14 = t->Branch("phi_sv_UP", &svFitPhi_UP, "phi_sv_UP/F");
      TBranch *newBranch15 = t->Branch("met_sv_UP", &svFitMET_UP, "met_sv_UP/F");
      TBranch *newBranch16 = t->Branch("mt_sv_UP", &svFitTransverseMass_UP, "mt_sv_UP/F");

      TBranch *newBranch17 = t->Branch("m_sv_DOWN", &svFitMass_DOWN, "m_sv_DOWN/F");
      TBranch *newBranch18 = t->Branch("p_sv_DOWN", &svFitPt_DOWN, "pt_sv_DOWN/F");
      TBranch *newBranch19 = t->Branch("eta_sv_DOWN", &svFitEta_DOWN, "eta_sv_DOWN/F");
      TBranch *newBranch20 = t->Branch("phi_sv_DOWN", &svFitPhi_DOWN, "phi_sv_DOWN/F");
      TBranch *newBranch21 = t->Branch("met_sv_DOWN", &svFitMET_DOWN, "met_sv_DOWN/F");
      TBranch *newBranch22 = t->Branch("mt_sv_DOWN", &svFitTransverseMass_DOWN, "mt_sv_DOWN/F");
      //}

      unsigned int evt, run, lumi;
      int evt2, run2, lumi2;
      evt2=0; run2=0; lumi2=0;
      float pt1;
      float eta1;
      float phi1;
      float pt2;
      float eta2;
      float phi2;
      float m2;
      float decayMode;
      float decayMode2;
      float covMatrix00;
      float covMatrix10;
      float covMatrix01;
      float covMatrix11;
      //float mvamet_ex, // uncorrected mva met px (float)
      //  mvamet_ey, // uncorrected mva met py (float)
      float  genPx    , // generator Z/W/Higgs px (float)
        genPy    , // generator Z/W/Higgs py (float)
        visPx    , // generator visible Z/W/Higgs px (float)
        visPy    , // generator visible Z/W/Higgs py (float)
        njets    ;  // number of jets (hadronic jet multiplicity) (int)

      // define MET
      double measuredMETx;
      double measuredMETy;
      float met;
      float metphi;
      TLorentzVector TMet(0,0,0,0);
      // define MET covariance
      TMatrixD covMET(2, 2);
      //ele/mu variables
      svFitStandalone::kDecayType decayType1 = svFitStandalone::kUndefinedDecayType; //svFitStandalone::kTauToElecDecay
      svFitStandalone::kDecayType decayType2 = svFitStandalone::kUndefinedDecayType; //svFitStandalone::kTauToElecDecay
      // Both masses should depend on decay mode and particle?
      float mass1;
      float mass2;
      std::string channel = "x";
      if(std::string(TreeToUse).find("muTauEvent")!= std::string::npos){
    decayType1 = svFitStandalone::kTauToMuDecay;
    decayType2 = svFitStandalone::kTauToHadDecay; 
    mass1 = 0.105658;
    mass2 = 0;
    channel = "mt";
      }
      else if(std::string(TreeToUse).find("eleTauEvent")!= std::string::npos){
    std::cout<<"eleTauTree"<<std::endl;
    decayType1 = svFitStandalone::kTauToElecDecay;
    decayType2 = svFitStandalone::kTauToHadDecay;
    mass1 = 0.00051100;
    mass2 = 0;
    channel = "et";
      }
      //else if(parser.stringValue("outputFile").find("_em.root") != std::string::npos){
      else if(std::string(TreeToUse).find("em")!= std::string::npos){
        //std::cout<< parser.stringValue("outputFile").c_str() << std::endl;
        std::cout<< "EMu sample" <<std::endl;
    decayType1 = svFitStandalone::kTauToElecDecay;
    decayType2 = svFitStandalone::kTauToMuDecay;
    mass1 = 0.00051100;
    mass2 = 0;
    channel = "em";
      }
      //else if(parser.stringValue("outputFile").find("_tt.root") != std::string::npos){
      else if((std::string(TreeToUse).find("tt")!= std::string::npos) ||
		(parser.stringValue("outputFile").find("_tt.root") != std::string::npos)){
        //std::cout<< parser.stringValue("outputFile").c_str() << std::endl;
        std::cout<< "Double Hadronic sample" <<std::endl;
    decayType1 = svFitStandalone::kTauToHadDecay;
    decayType2 = svFitStandalone::kTauToHadDecay;
    mass1 = 0.13957;
    mass2 = 0.13957;
    channel = "tt";
      }
      else{
    std::cout<<"TreeToUse "<< std::string(TreeToUse)<<" does not match muTauEvent or eleTauEvent... Skipping!!"<<std::endl;
    continue;
      }
      
      //MET, MET Covariance, lepton objects,
      TBranch *pt1branch;
      //TBranch *newBranch = t->Branch(parser.stringValue("branch").c_str(),&svFitMass,(parser.stringValue("branch")+"/F").c_str());
      //newBranch->SetFile("svfitout.root:/input/muTauEventTree");

      if(channel=="tt" || channel=="em"){
        t->SetBranchAddress("evt",&evt2);
        t->SetBranchAddress("run",&run2);
        t->SetBranchAddress("lumi",&lumi2);
        t->SetBranchAddress("jetVeto30RecoilZTT", &njets);
      }
      if(channel=="tt") {
        t->SetBranchAddress("t1Pt",&pt1,&pt1branch);
        t->SetBranchAddress("t1Eta",&eta1);
        t->SetBranchAddress("t1Phi",&phi1);
        t->SetBranchAddress("t1Mass",&mass1);
        t->SetBranchAddress("t2Pt",&pt2);
        t->SetBranchAddress("t2Eta",&eta2);
        t->SetBranchAddress("t2Phi",&phi2);
        t->SetBranchAddress("t2Mass",&mass2);
        t->SetBranchAddress("t1DecayMode",&decayMode);
        t->SetBranchAddress("t2DecayMode",&decayMode2);
        t->SetBranchAddress("t1_t2_MvaMetCovMatrix00",&covMatrix00);
        t->SetBranchAddress("t1_t2_MvaMetCovMatrix01",&covMatrix01);
        t->SetBranchAddress("t1_t2_MvaMetCovMatrix10",&covMatrix10);
        t->SetBranchAddress("t1_t2_MvaMetCovMatrix11",&covMatrix11);
        t->SetBranchAddress("t1_t2_MvaMet",&met);
        t->SetBranchAddress("t1_t2_MvaMetPhi",&metphi);
      }
      else if(channel=="em") {
        t->SetBranchAddress("ePt",&pt1,&pt1branch);
        t->SetBranchAddress("eEta",&eta1);
        t->SetBranchAddress("ePhi",&phi1);
        t->SetBranchAddress("mPt",&pt2);
        t->SetBranchAddress("mEta",&eta2);
        t->SetBranchAddress("mPhi",&phi2);
        t->SetBranchAddress("e_m_MvaMetCovMatrix00",&covMatrix00);
        t->SetBranchAddress("e_m_MvaMetCovMatrix01",&covMatrix01);
        t->SetBranchAddress("e_m_MvaMetCovMatrix10",&covMatrix10);
        t->SetBranchAddress("e_m_MvaMetCovMatrix11",&covMatrix11);
        t->SetBranchAddress("e_m_MvaMet",&met);
        t->SetBranchAddress("e_m_MvaMetPhi",&metphi);
      }
      else {
        t->SetBranchAddress("EVENT",&evt);
        t->SetBranchAddress("run",&run);
        t->SetBranchAddress("lumi",&lumi);
        t->SetBranchAddress("pt_1",&pt1,&pt1branch);
        t->SetBranchAddress("eta_1",&eta1);
        t->SetBranchAddress("phi_1",&phi1);
        t->SetBranchAddress("pt_2",&pt2);
        t->SetBranchAddress("eta_2",&eta2);
        t->SetBranchAddress("phi_2",&phi2);
        t->SetBranchAddress("m_2",&m2);
        t->SetBranchAddress("tauDecayMode",&decayMode);
        t->SetBranchAddress("t2DecayMode",&decayMode2);
        t->SetBranchAddress("mvacov00",&covMatrix00);
        t->SetBranchAddress("mvacov01",&covMatrix01);
        t->SetBranchAddress("mvacov10",&covMatrix10);
        t->SetBranchAddress("mvacov11",&covMatrix11);
        t->SetBranchAddress("mvamet",&met);
        t->SetBranchAddress("mvametphi",&metphi);
        t->SetBranchAddress( "njets", &njets);
      }
      t->SetBranchAddress( "genpX", &genPx);
      t->SetBranchAddress( "genpY", &genPy);
      t->SetBranchAddress( "vispX", &visPx);
      t->SetBranchAddress( "vispY", &visPy);
      // use this RooT file when running on aMC@NLO DY and W+Jets MC samples
      RecoilCorrector* recoilMvaMetCorrector = new RecoilCorrector(recoilFileName);
      std::cout<<"recoiltype "<<recoilType<<" recoilFileName "<<recoilFileName<<std::endl;

      printf("Found tree -> weighting\n");

      // Only open this once!
      edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
      TH1::AddDirectory(false);  
      TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());  

      for(Int_t i=0;i<t->GetEntries();++i)
    {
      t->GetEntry(i);

      if(channel=="tt" || channel=="em"){
        evt = (unsigned int)evt2;
        run = (unsigned int)run2;
        lumi = (unsigned int)lumi2;
      }

      TMet.SetPtEtaPhiM(met,0,metphi,0);
      measuredMETx = met*TMath::Cos(metphi);
      measuredMETy = met*TMath::Sin(metphi);

      //Recoil Correction time
      // Correct WJets recoil for faked lepton / extra jet
      int recoilNJets;
      if(isWJets) {
        recoilNJets = njets + 1;
        std::cout << "njets: " << njets << " recoilNJets: " << recoilNJets << std::endl;
      }
      else recoilNJets = njets;


      if(recoilType != 0){

        //RecoilCorrector recoilMvaMetCorrector(recoilFileName);
        //recoilMvaMetCorrector->CorrectByMeanResolution(measuredMETx, // uncorrected mva met px (float)
        //                      measuredMETy, // uncorrected mva met py (float)
        //                      genPx, // generator Z/W/Higgs px (float)
        //                      genPy, // generator Z/W/Higgs py (float)
        //                      visPx, // generator visible Z/W/Higgs px (float)
        //                      visPy, // generator visible Z/W/Higgs py (float)
        //                      //recoilNJets,  // number of jets (hadronic jet multiplicity) (int)
        //                      njets,  // number of jets (hadronic jet multiplicity) (int)
        //                      mvametcorr_ex, // corrected met px (float)
        //                      mvametcorr_ey  // corrected met py (float)
        //                      );
        //                      //);}
        //std::cout << "njets       mvamet_ex: " << mvametcorr_ex << "  mvamet_ey: " << mvametcorr_ey << std::endl;

        recoilMvaMetCorrector->CorrectByMeanResolution(measuredMETx, // uncorrected mva met px (float)
                              measuredMETy, // uncorrected mva met py (float)
                              genPx, // generator Z/W/Higgs px (float)
                              genPy, // generator Z/W/Higgs py (float)
                              visPx, // generator visible Z/W/Higgs px (float)
                              visPy, // generator visible Z/W/Higgs py (float)
                              recoilNJets,  // number of jets (hadronic jet multiplicity) (int)
                              mvametcorr_ex, // corrected met px (float)
                              mvametcorr_ey  // corrected met py (float)
                              );
        std::cout << "recoilNJets mvamet_ex: " << mvametcorr_ex << "  mvamet_ey: " << mvametcorr_ey << std::endl;
      }
      else{
        mvametcorr_ex = measuredMETx;
        mvametcorr_ey = measuredMETy;
      }

      mvamet = TMath::Sqrt( mvametcorr_ex*mvametcorr_ex + mvametcorr_ey*mvametcorr_ey);
      mvametphi = TMath::ATan2( mvametcorr_ey, mvametcorr_ex );
      std::cout << " - mvamet "<<mvamet<<" mvametphi "<<mvametphi<<std::endl;

      if(channel=="et" || channel=="mt"){
        mass2 = m2;
        if(decayMode==0)
        mass2 = 0.13957;
      }
      //mods needed for tau tau here

      covMET[0][0] =  covMatrix00;
      covMET[1][0] =  covMatrix10;
      covMET[0][1] =  covMatrix01;
      covMET[1][1] =  covMatrix11;
      //std::cout<<"getting decay mode "<< decayMode<<std::endl;
      if((channel!="tt"&&channel!="em")&&(decayMode==0||decayMode==1||decayMode==10)){
        // define lepton four vectors
        std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
        // check if electron or muon
        measuredTauLeptons.push_back(
         svFitStandalone::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1)
                     ); // tau -> electron decay (Pt, eta, phi, mass)
        
        measuredTauLeptons.push_back(
          svFitStandalone::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2, decayMode)
                     ); // tau -> 1prong0pi0 hadronic decay (Pt, eta, phi, mass, pat::Tau.decayMode())
        std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1 << " mass1 " << mass1 << " pt2: "<< pt2<< " mass2: "<< mass2 <<std::endl;        
        //modified
        runSVFit(measuredTauLeptons, inputFile_visPtResolution, mvametcorr_ex, mvametcorr_ey, covMET, 0, svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass);
        std::cout<<"finished runningSVFit"<<std::endl;
      } // eTau / muTau



      else if(channel=="em"){
        // define lepton four vectors
        std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
        measuredTauLeptons.push_back(
         svFitStandalone::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1));
        measuredTauLeptons.push_back(
         svFitStandalone::MeasuredTauLepton(decayType2,  pt2, eta2, phi2, mass2));

        std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1 << " mass1 " << mass1 << " pt2: "<< pt2<< " mass2: "<< mass2 <<std::endl;        
        //modified
        runSVFit(measuredTauLeptons, inputFile_visPtResolution, mvametcorr_ex, mvametcorr_ey, covMET, 0, svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass);
        std::cout<<"finished running non-EES SVFit in EMu"<<std::endl;
        if(doES) {

          // Only shift Electrons
          // 1% shift in barrel and 2.5% shift in endcap
          // applied to electrons in emu channel
          float etaBarrelEndcap  = 1.479;
          float ES_UP_scale;
          if (abs(eta1) < etaBarrelEndcap) ES_UP_scale = 1.01;
          else ES_UP_scale = 1.025;
          double pt1_UP = pt1 * ES_UP_scale;
          std::cout << "E eta: " << eta1 << " ees SF: " << ES_UP_scale << std::endl;
          double mvametcorr_ex_UP, mvametcorr_ey_UP;
          double dx1_UP, dy1_UP;
          dx1_UP = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_UP_scale ) - 1.);
          dy1_UP = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_UP_scale ) - 1.);
          mvametcorr_ex_UP = mvametcorr_ex + dx1_UP;
          mvametcorr_ey_UP = mvametcorr_ey + dy1_UP;
          std::cout << "px1 " << pt1 * TMath::Cos( phi1 ) << "  met px1 cor " << dx1_UP <<std::endl;
          std::cout << "py1 " << pt1 * TMath::Sin( phi1 ) << "  met py1 cor " << dy1_UP <<std::endl;
          std::cout << "pt1 " << pt1 << "  pt1_up " << pt1_UP <<std::endl;
          std::cout << "mvamet_ex_i " << mvametcorr_ex << " ees: " << mvametcorr_ex_UP << std::endl;
          std::cout << "mvamet_ey_i " << mvametcorr_ey << " ees: " << mvametcorr_ey_UP << std::endl;
          
          std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptonsUP;
          measuredTauLeptonsUP.push_back(
           svFitStandalone::MeasuredTauLepton(decayType1, pt1_UP, eta1,  phi1, mass1));
          measuredTauLeptonsUP.push_back(
           svFitStandalone::MeasuredTauLepton(decayType2,  pt2, eta2, phi2, mass2));

          std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1_UP << " mass1 " << mass1 << " pt2: "<< pt2 << " mass2: "<< mass2 <<std::endl;        
          runSVFit(measuredTauLeptonsUP, inputFile_visPtResolution, mvametcorr_ex_UP, mvametcorr_ey_UP, covMET, 0, svFitMass_UP, svFitPt_UP, svFitEta_UP, svFitPhi_UP, svFitMET_UP, svFitTransverseMass_UP);
          std::cout<<"finished runningSVFit in EMu EES Up"<<std::endl;


          // EES Down
          float ES_DOWN_scale;
          if (abs(eta1) < etaBarrelEndcap) ES_DOWN_scale = 0.99;
          else ES_DOWN_scale = 0.975;
          std::cout << "E eta: " << eta1 << " ees SF: " << ES_DOWN_scale << std::endl;
          double pt1_DOWN;
          pt1_DOWN = pt1 * ES_DOWN_scale;
          double mvametcorr_ex_DOWN, mvametcorr_ey_DOWN;
          double dx1_DOWN, dy1_DOWN;
          dx1_DOWN = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_DOWN_scale ) - 1.);
          dy1_DOWN = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_DOWN_scale ) - 1.);
          mvametcorr_ex_DOWN = mvametcorr_ex + dx1_DOWN;
          mvametcorr_ey_DOWN = mvametcorr_ey + dy1_DOWN;
          std::cout << "px1 " << pt1 * TMath::Cos( phi1 ) << "  met px1 cor " << dx1_DOWN <<std::endl;
          std::cout << "py1 " << pt1 * TMath::Sin( phi1 ) << "  met py1 cor " << dy1_DOWN <<std::endl;
          std::cout << "mvamet_ex_i " << mvametcorr_ex << " ees: " << mvametcorr_ex_DOWN << std::endl;
          std::cout << "mvamet_ey_i " << mvametcorr_ey << " ees: " << mvametcorr_ey_DOWN << std::endl;

          std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptonsDOWN;
          measuredTauLeptonsDOWN.push_back(
           svFitStandalone::MeasuredTauLepton(decayType1, pt1_DOWN, eta1,  phi1, mass1));
          measuredTauLeptonsDOWN.push_back(
           svFitStandalone::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2));

          std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1_DOWN << " mass1 " << mass1 << " pt2: "<< pt2 << " mass2: "<< mass2 <<std::endl;
          runSVFit(measuredTauLeptonsDOWN, inputFile_visPtResolution, mvametcorr_ex_DOWN, mvametcorr_ey_DOWN, covMET, 0, svFitMass_DOWN, svFitPt_DOWN, svFitEta_DOWN, svFitPhi_DOWN, svFitMET_DOWN, svFitTransverseMass_DOWN);
          std::cout<<"finished runningSVFit in EMu EES Down"<<std::endl;


        } // end doES
      } // eMu



      else if(channel=="tt"){
        // define lepton four vectors
        std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
        // Add Tau of higest Pt first
        if (pt1 > pt2) {
          measuredTauLeptons.push_back(
           svFitStandalone::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1, decayMode));
        
          measuredTauLeptons.push_back(
           svFitStandalone::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2, decayMode2));
        }
        else {
          measuredTauLeptons.push_back(
           svFitStandalone::MeasuredTauLepton(decayType2,  pt2, eta2, phi2,  mass2, decayMode2));

          measuredTauLeptons.push_back(
           svFitStandalone::MeasuredTauLepton(decayType1, pt1, eta1,  phi1, mass1, decayMode));
        }

        std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1 << " mass1 " << mass1 << " pt2: "<< pt2<< " mass2: "<< mass2 <<std::endl;        
        runSVFit(measuredTauLeptons, inputFile_visPtResolution, mvametcorr_ex, mvametcorr_ey, covMET, 0, svFitMass, svFitPt, svFitEta, svFitPhi, svFitMET, svFitTransverseMass);
        std::cout<<"finished running non-TES SVFit in TT"<<std::endl;

        if(doES) {

          // First ES UP x 1.03
          float ES_UP_scale = 1.03;
          double pt1_UP, pt2_UP;
          pt1_UP = pt1 * ES_UP_scale;
          pt2_UP = pt2 * ES_UP_scale;
          double mvametcorr_ex_UP, mvametcorr_ey_UP;
          double dx1_UP, dy1_UP, dx2_UP, dy2_UP;
          dx1_UP = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_UP_scale ) - 1.);
          dy1_UP = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_UP_scale ) - 1.);
          dx2_UP = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_UP_scale ) - 1.);
          dy2_UP = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_UP_scale ) - 1.);
          mvametcorr_ex_UP = mvametcorr_ex + dx1_UP + dx2_UP;
          mvametcorr_ey_UP = mvametcorr_ey + dy1_UP + dy2_UP;
          std::cout << "px1 " << pt1 * TMath::Cos( phi1 ) << "  met px1 cor " << dx1_UP <<std::endl;
          std::cout << "px2 " << pt2 * TMath::Cos( phi2 ) << "  met px2 cor " << dx2_UP <<std::endl;
          std::cout << "py1 " << pt1 * TMath::Sin( phi1 ) << "  met py1 cor " << dy1_UP <<std::endl;
          std::cout << "py2 " << pt2 * TMath::Sin( phi2 ) << "  met py2 cor " << dy2_UP <<std::endl;
          std::cout << "mvamet_ex_i " << mvametcorr_ex << " tes: " << mvametcorr_ex_UP << std::endl;
          std::cout << "mvamet_ey_i " << mvametcorr_ey << " tes: " << mvametcorr_ey_UP << std::endl;
          
          std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptonsUP;
          
          // Add Tau of higest Pt first
          if (pt1 > pt2) {
            measuredTauLeptonsUP.push_back(
             svFitStandalone::MeasuredTauLepton(decayType1, pt1_UP, eta1,  phi1, mass1, decayMode));
          
            measuredTauLeptonsUP.push_back(
             svFitStandalone::MeasuredTauLepton(decayType2,  pt2_UP, eta2, phi2,  mass2, decayMode2));
          }
          else {
            measuredTauLeptonsUP.push_back(
             svFitStandalone::MeasuredTauLepton(decayType2,  pt2_UP, eta2, phi2,  mass2, decayMode2));

            measuredTauLeptonsUP.push_back(
             svFitStandalone::MeasuredTauLepton(decayType1, pt1_UP, eta1,  phi1, mass1, decayMode));
          }

          std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1_UP << " mass1 " << mass1 << " pt2: "<< pt2_UP<< " mass2: "<< mass2 <<std::endl;        
          runSVFit(measuredTauLeptonsUP, inputFile_visPtResolution, mvametcorr_ex_UP, mvametcorr_ey_UP, covMET, 0, svFitMass_UP, svFitPt_UP, svFitEta_UP, svFitPhi_UP, svFitMET_UP, svFitTransverseMass_UP);
          std::cout<<"finished runningSVFit in TT TES Up"<<std::endl;


          // Second ES Down, x 0.97
          float ES_DOWN_scale = 0.97;
          double pt1_DOWN, pt2_DOWN;
          pt1_DOWN = pt1 * ES_DOWN_scale;
          pt2_DOWN = pt2 * ES_DOWN_scale;
          double mvametcorr_ex_DOWN, mvametcorr_ey_DOWN;
          double dx1_DOWN, dy1_DOWN, dx2_DOWN, dy2_DOWN;
          dx1_DOWN = pt1 * TMath::Cos( phi1 ) * (( 1. / ES_DOWN_scale ) - 1.);
          dy1_DOWN = pt1 * TMath::Sin( phi1 ) * (( 1. / ES_DOWN_scale ) - 1.);
          dx2_DOWN = pt2 * TMath::Cos( phi2 ) * (( 1. / ES_DOWN_scale ) - 1.);
          dy2_DOWN = pt2 * TMath::Sin( phi2 ) * (( 1. / ES_DOWN_scale ) - 1.);
          mvametcorr_ex_DOWN = mvametcorr_ex + dx1_DOWN + dx2_DOWN;
          mvametcorr_ey_DOWN = mvametcorr_ey + dy1_DOWN + dy2_DOWN;
          std::cout << "px1 " << pt1 * TMath::Cos( phi1 ) << "  met px1 cor " << dx1_DOWN <<std::endl;
          std::cout << "px2 " << pt2 * TMath::Cos( phi2 ) << "  met px2 cor " << dx2_DOWN <<std::endl;
          std::cout << "py1 " << pt1 * TMath::Sin( phi1 ) << "  met py1 cor " << dy1_DOWN <<std::endl;
          std::cout << "py2 " << pt2 * TMath::Sin( phi2 ) << "  met py2 cor " << dy2_DOWN <<std::endl;
          std::cout << "mvamet_ex_i " << mvametcorr_ex << " tes: " << mvametcorr_ex_DOWN << std::endl;
          std::cout << "mvamet_ey_i " << mvametcorr_ey << " tes: " << mvametcorr_ey_DOWN << std::endl;


          std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptonsDOWN;
          // Add Tau of higest Pt first
          
          if (pt1 > pt2) {
            measuredTauLeptonsDOWN.push_back(
             svFitStandalone::MeasuredTauLepton(decayType1, pt1_DOWN, eta1,  phi1, mass1, decayMode));
          
            measuredTauLeptonsDOWN.push_back(
             svFitStandalone::MeasuredTauLepton(decayType2,  pt2_DOWN, eta2, phi2,  mass2, decayMode2));
          }
          else {
            measuredTauLeptonsDOWN.push_back(
             svFitStandalone::MeasuredTauLepton(decayType2,  pt2_DOWN, eta2, phi2,  mass2, decayMode2));

            measuredTauLeptonsDOWN.push_back(
             svFitStandalone::MeasuredTauLepton(decayType1, pt1_DOWN, eta1,  phi1, mass1, decayMode));
          }

          std::cout<< "evt: "<<evt<<" run: "<<run<<" lumi: "<<lumi<< " pt1 " << pt1_DOWN << " mass1 " << mass1 << " pt2: "<< pt2_DOWN<< " mass2: "<< mass2 <<std::endl;
          runSVFit(measuredTauLeptonsDOWN, inputFile_visPtResolution, mvametcorr_ex_DOWN, mvametcorr_ey_DOWN, covMET, 0, svFitMass_DOWN, svFitPt_DOWN, svFitEta_DOWN, svFitPhi_DOWN, svFitMET_DOWN, svFitTransverseMass_DOWN);
          std::cout<<"finished runningSVFit in TT ES Down"<<std::endl;


        } // end doES
      } // Double Hadronic (TT)

      else {
        svFitMass = -100;
        svFitPt = -100;
        svFitEta = -100;
        svFitPhi = -100;
        svFitMET = -100;
        svFitTransverseMass = -100;
      }
      //std::cout << "\n\nex: " << mvametcorr_ex << "   ey: " << mvametcorr_ey <<  " phi: " << mvametphi<<"\n"<<std::endl; 
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
      //if(doES) {
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
      //} // end do ES
    }
      delete inputFile_visPtResolution;
      dir->cd();
      t->Write("",TObject::kOverwrite);
      strcpy(TreeToUse,stringA) ;

    }
  }
}

void runSVFit(std::vector<svFitStandalone::MeasuredTauLepton> & measuredTauLeptons, TFile * inputFile_visPtResolution, double measuredMETx, double measuredMETy, TMatrixD &covMET, float num, float &svFitMass, float& svFitPt, float &svFitEta, float &svFitPhi, float &svFitMET, float &svFitTransverseMass){

  //edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
  //TH1::AddDirectory(false);  
  //TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());  
  //float svFitMass = 0;
  SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, 0);
  algo.addLogM(false);  
  algo.shiftVisPt(true, inputFile_visPtResolution);
  algo.integrateMarkovChain();
  svFitMass = algo.getMass(); // return value is in units of GeV
  svFitPt = algo.pt();
  svFitEta = algo.eta();
  svFitPhi = algo.phi();
  svFitMET = algo.fittedMET().Rho();
  svFitTransverseMass = algo.transverseMass();
  if ( algo.isValidSolution() ) {
    std::cout << "found mass = " << svFitMass << std::endl;
  } else {
    std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
  }
  //delete inputFile_visPtResolution;

}

//Thank you Renee Brun :)
void CopyDir(TDirectory *source, optutl::CommandLineParser parser) {
  //copy all objects and subdirs of directory source as a subdir of the current directory
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir; 
  if(source->GetName()!=parser.stringValue("outputFile")){
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
  if(gSystem->AccessPathName(parser.stringValue("outputFile").c_str())) {
    gSystem->CopyFile("hsimple.root", parser.stringValue("outputFile").c_str());
  }

  fNew->cd();
  CopyFile(parser.stringValue("outputFile").c_str(),parser);
  fNew->ls();
  fNew->Close();

}

