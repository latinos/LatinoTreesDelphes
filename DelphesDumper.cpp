
/*c++ -O2 -lm `root-config --cflags --glibs` -L /afs/cern.ch/work/d/depoyraz/VBS/TP/CMSSW_6_2_0_SLHC19/src/Delphes_Two/Delphes -I /afs/cern.ch/work/d/depoyraz/VBS/TP/CMSSW_6_2_0_SLHC19/src/Delphes_Two/Delphes -l Delphes -o DelphesDumper DelphesDumper.cpp
  ./DelphesDumper rootfile  output.root
*/


#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"


using namespace std;

//****************************************************************************************

//****************************************************************************************
// runs over all entries in the delphes tree

 
//****************************************************************************************

struct Lepton {
  unsigned int type; //0:electron , 1:muon
  unsigned int index;
  float lpt, leta, lphi, liso, lisoDBeta, lisoRhoCorr, lsumChargedHadron, lsumNeutral, lsumChargedPU, lsumAllParticles;
  double lch;
};

struct lheParticleDescendingPt
{
  bool operator() (LHEParticle* a, LHEParticle* b)
  {
    return a->PT > b->PT;
  }
};	
struct JetDescendingPt
{
  bool operator() (Jet* a, Jet* b)
  {
    return a->PT > b->PT;
  }
};


struct leptonDescendingPt
{
  bool operator() (const Lepton& a, const Lepton& b)
  {
    return a.lpt > b.lpt;
  }
};

/* struct genParticleDescendingPt 
{
  bool operator() (GenParticle* a, GenParticle* b) 
  {     
    return a->PT > b->PT;
  }
}; 
*/  

float DeltaR(float eta1, float eta2, float phi1, float phi2)
{
  float deta = eta2 - eta1;
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  float DELTAR = sqrt(pow(dphi,2)+pow(deta,2))*1.0;
  return DELTAR;
}

float DeltaPhi(float phi1, float phi2)
{
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  return dphi;
}


//****************************************************************************************
// main
int main (int argc, char *argv[])
{
  //----------------------------------------------------------------------------------------
  //importing delphes libraries
  gSystem->Load("libDelphes");
  //----------------------------------------------------------------------------------------
  //complex object definitions
  //change it for root input
  TChain* delphesNtuples = new TChain("Delphes");
  ExRootTreeReader *delphesTree = new ExRootTreeReader(delphesNtuples);
  //new
  //TTree *LT ;//=  (TTree*)(delphesNtuples);
  TFile* outputFile = TFile::Open(argv[2],"recreate");
  TTree* easyTree = new TTree("easyDelphes","easyDelphes");
  //----------------------------------------------------------------------------------------
  // reading input files
  if(argc < 3)
    {
      cout << "ERROR: not enough info provided" << endl;
      return 0;
    }
    
	
  delphesNtuples -> Add(argv[1]);
	
  delphesNtuples -> BranchRef();
    
  Long64_t numberOfEntries = delphesTree->GetEntries();
  cout<<"##### Number of Entires in the Delphes Tree is: "<<numberOfEntries<<endl;
    
    
	
    
  //----------------------------------------------------------------------------------------
  //variable management
  //
  // -> all objects in the output tree are stored in decreasing pt order.
    
  //--------- getting objects from the delphes tree
  TClonesArray* branchLHEParticle = delphesTree->UseBranch("LHEParticles");
  TClonesArray* branchEl = delphesTree->UseBranch("Electron");
  TClonesArray* branchMu = delphesTree->UseBranch("Muon");
  TClonesArray* branchGenJet = delphesTree->UseBranch("GenJet");
  TClonesArray* branchTrackJet = delphesTree->UseBranch("TrackJet");
  TClonesArray* branchJet = delphesTree->UseBranch("JetPUID");
  TClonesArray* branchPuppiJet = delphesTree->UseBranch("PuppiJetPUID");
  //TClonesArray* branchGenParticle = delphesTree->UseBranch("GenParticles");
    
  TClonesArray* branchMET = delphesTree->UseBranch("MissingET");
  TClonesArray* branchPuppiMET = delphesTree->UseBranch("PuppiMissingET");
  TClonesArray* branchGenMET = delphesTree->UseBranch("GenMissingET");
    
  TClonesArray* branchNPU = delphesTree->UseBranch("NPU");
  TClonesArray* branchGlobalRhokt4 = delphesTree->UseBranch("GlobalRhoKt4");
  TClonesArray* branchGlobalRhoGFJ= delphesTree->UseBranch("GlobalRhoGridFastJet");
  TClonesArray* branchRhokt4 = delphesTree->UseBranch("RhoKt4");
  TClonesArray* branchRhoGFJ= delphesTree->UseBranch("RhoGridFastJet");
  TClonesArray* branchPuppiRhokt4 = delphesTree->UseBranch("PuppiRhoKt4");
  TClonesArray* branchPuppiRhoGFJ= delphesTree->UseBranch("PuppiRhoGridFastJet");
  

    

    
    
  //--------- Creating branches for the new (light) tree
    
  //--------- LHE Information
  int nlhe=4;
  int nlhes=2;
    
  float leptonLHEpt_tmp[nlhe], leptonLHEeta_tmp[nlhe], leptonLHEphi_tmp[nlhe], leptonLHEpid_tmp[nlhe], leptonLHEch_tmp[nlhe], leptonLHEm_tmp[nlhe] ;
  float neutrinoLHEpt_tmp[nlhe], neutrinoLHEeta_tmp[nlhe], neutrinoLHEphi_tmp[nlhe], neutrinoLHEpid_tmp[nlhe];
  float jetLHEPartonpt_tmp[nlhe], jetLHEPartoneta_tmp[nlhe], jetLHEPartonphi_tmp[nlhe], jetLHEPartonpid_tmp[nlhe];
  float jetLHEGluonpt_tmp[nlhe], jetLHEGluoneta_tmp[nlhe], jetLHEGluonphi_tmp[nlhe], jetLHEGluonpid_tmp[nlhe] ;
  float vbosonLHEpt_tmp[nlhe], vbosonLHEeta_tmp[nlhe], vbosonLHEphi_tmp[nlhe], vbosonLHEpid_tmp[nlhe], vbosonLHEch_tmp[nlhe],vbosonLHEm_tmp[nlhe] ;
	
	
  for(int ilhe =0; ilhe<nlhe; ilhe++){
    TString leptonLHEptStr = "leptonLHEpt"; leptonLHEptStr += (ilhe+1);
    TString leptonLHEetaStr = "leptonLHEeta"; leptonLHEetaStr += (ilhe+1);
    TString leptonLHEphiStr = "leptonLHEphi"; leptonLHEphiStr += (ilhe+1);
    TString leptonLHEpidStr = "leptonLHEpid"; leptonLHEpidStr += (ilhe+1);
    TString leptonLHEchStr = "leptonLHEch"; leptonLHEchStr += (ilhe+1);
    TString leptonLHEmStr = "leptonLHEm"; leptonLHEmStr += (ilhe+1);
    TString neutrinoLHEptStr = "neutrinoLHEpt"; neutrinoLHEptStr += (ilhe+1);
    TString neutrinoLHEetaStr = "neutrinoLHEeta"; neutrinoLHEetaStr += (ilhe+1);
    TString neutrinoLHEphiStr = "neutrinoLHEphi"; neutrinoLHEphiStr += (ilhe+1);
    TString neutrinoLHEpidStr = "neutrinoLHEpid"; neutrinoLHEpidStr += (ilhe+1);
    TString jetLHEPartonptStr = "jetLHEPartonpt"; jetLHEPartonptStr += (ilhe+1);
    TString jetLHEPartonetaStr = "jetLHEPartoneta"; jetLHEPartonetaStr += (ilhe+1);
    TString jetLHEPartonphiStr = "jetLHEPartonphi"; jetLHEPartonphiStr += (ilhe+1);
    TString jetLHEPartonpidStr = "jetLHEPartonpid"; jetLHEPartonpidStr += (ilhe+1);
        
    easyTree -> Branch(leptonLHEptStr,&leptonLHEpt_tmp[ilhe],leptonLHEptStr+"/F");
    easyTree -> Branch(leptonLHEetaStr,&leptonLHEeta_tmp[ilhe],leptonLHEetaStr+"/F");
    easyTree -> Branch(leptonLHEphiStr,&leptonLHEphi_tmp[ilhe],leptonLHEphiStr+"/F");
    easyTree -> Branch(leptonLHEpidStr,&leptonLHEpid_tmp[ilhe],leptonLHEpidStr+"/F");
    easyTree -> Branch(leptonLHEchStr,&leptonLHEch_tmp[ilhe],leptonLHEchStr+"/F");
    easyTree -> Branch(leptonLHEmStr,&leptonLHEm_tmp[ilhe],leptonLHEmStr+"/F");
    easyTree -> Branch(neutrinoLHEptStr,&neutrinoLHEpt_tmp[ilhe],neutrinoLHEptStr+"/F");
    easyTree -> Branch(neutrinoLHEetaStr,&neutrinoLHEeta_tmp[ilhe],neutrinoLHEetaStr+"/F");
    easyTree -> Branch(neutrinoLHEphiStr,&neutrinoLHEphi_tmp[ilhe],neutrinoLHEphiStr+"/F");
    easyTree -> Branch(neutrinoLHEpidStr,&neutrinoLHEpid_tmp[ilhe],neutrinoLHEpidStr+"/F");
    easyTree -> Branch(jetLHEPartonptStr,&jetLHEPartonpt_tmp[ilhe],jetLHEPartonptStr+"/F");
    easyTree -> Branch(jetLHEPartonetaStr,&jetLHEPartoneta_tmp[ilhe],jetLHEPartonetaStr+"/F");
    easyTree -> Branch(jetLHEPartonphiStr,&jetLHEPartonphi_tmp[ilhe],jetLHEPartonphiStr+"/F");
    easyTree -> Branch(jetLHEPartonpidStr,&jetLHEPartonpid_tmp[ilhe],jetLHEPartonpidStr+"/F");
        
  }
    
  for(int ilhe =0; ilhe<nlhes; ilhe++){
    TString vbosonLHEptStr = "vbosonLHEpt"; vbosonLHEptStr += (ilhe+1);
    TString vbosonLHEetaStr = "vbosonLHEeta"; vbosonLHEetaStr += (ilhe+1);
    TString vbosonLHEphiStr = "vbosonLHEphi"; vbosonLHEphiStr += (ilhe+1);
    TString vbosonLHEpidStr = "vbosonLHEpid"; vbosonLHEpidStr += (ilhe+1);
    TString vbosonLHEchStr = "vbosonLHEch"; vbosonLHEchStr += (ilhe+1);
    TString vbosonLHEmStr = "vbosonLHEm"; vbosonLHEmStr += (ilhe+1);
    TString jetLHEGluonptStr = "jetLHEGluonpt"; jetLHEGluonptStr += (ilhe+1);
    TString jetLHEGluonetaStr = "jetLHEGluoneta"; jetLHEGluonetaStr += (ilhe+1);
    TString jetLHEGluonphiStr = "jetLHEGluonphi"; jetLHEGluonphiStr += (ilhe+1);
    TString jetLHEGluonpidStr = "jetLHEGluonpid"; jetLHEGluonpidStr += (ilhe+1);
        
    easyTree -> Branch(vbosonLHEptStr,&vbosonLHEpt_tmp[ilhe],vbosonLHEptStr+"/F");
    easyTree -> Branch(vbosonLHEetaStr,&vbosonLHEeta_tmp[ilhe],vbosonLHEetaStr+"/F");
    easyTree -> Branch(vbosonLHEphiStr,&vbosonLHEphi_tmp[ilhe],vbosonLHEphiStr+"/F");
    easyTree -> Branch(vbosonLHEpidStr,&vbosonLHEpid_tmp[ilhe],vbosonLHEpidStr+"/F");
    easyTree -> Branch(vbosonLHEchStr,&vbosonLHEch_tmp[ilhe],vbosonLHEchStr+"/F");
    easyTree -> Branch(vbosonLHEmStr,&vbosonLHEm_tmp[ilhe],vbosonLHEmStr+"/F");
    easyTree -> Branch(jetLHEGluonptStr,&jetLHEGluonpt_tmp[ilhe],jetLHEGluonptStr+"/F");
    easyTree -> Branch(jetLHEGluonetaStr,&jetLHEGluoneta_tmp[ilhe],jetLHEGluonetaStr+"/F");
    easyTree -> Branch(jetLHEGluonphiStr,&jetLHEGluonphi_tmp[ilhe],jetLHEGluonphiStr+"/F");
    easyTree -> Branch(jetLHEGluonpidStr,&jetLHEGluonpid_tmp[ilhe],jetLHEGluonpidStr+"/F");
	            
  }
  
  //------Gen Particles  
  
/*  int ngen=4;
  float leptonGenpt_tmp[ngen];
  float leptonGenpid_tmp[ngen];
  float leptonGenphi_tmp[ngen];
  float leptonGeneta_tmp[ngen];
  float neutrinoGenpt_tmp[ngen];
  float neutrinoGenpid_tmp[ngen];
  float neutrinoGenphi_tmp[ngen];
  float neutrinoGeneta_tmp[ngen];
  
  	
  for(int igen = 0; igen<ngen; igen++){
    TString leptonGenptStr = "leptonGenpt"; leptonGenptStr += (igen+1);
    TString leptonGenetaStr = "leptonGeneta"; leptonGenetaStr += (igen+1);
    TString leptonGenphiStr = "leptonGenphi"; leptonGenphiStr += (igen+1);
    TString leptonGenpidStr = "leptonGenpid"; leptonGenpidStr += (igen+1);
    TString neutrinoGenptStr = "neutrinoGenpt"; neutrinoGenptStr += (igen+1);
    TString neutrinoGenetaStr = "neutrinoGeneta"; neutrinoGenetaStr += (igen+1);
    TString neutrinoGenphiStr = "neutrinoGenphi"; neutrinoGenphiStr += (igen+1);
    TString neutrinoGenpidStr = "neutrinoGenpid"; neutrinoGenpidStr += (igen+1);
    
    easyTree -> Branch(leptonGenptStr,&leptonGenpt_tmp[igen],leptonGenptStr+"/F");
    easyTree -> Branch(leptonGenetaStr,&leptonGeneta_tmp[igen],leptonGenetaStr+"/F");
    easyTree -> Branch(leptonGenphiStr,&leptonGenphi_tmp[igen],leptonGenphiStr+"/F");
    easyTree -> Branch(leptonGenpidStr,&leptonGenpid_tmp[igen],leptonGenpidStr+"/F");
    
    easyTree -> Branch(neutrinoGenptStr,&neutrinoGenpt_tmp[igen],neutrinoGenptStr+"/F");
    easyTree -> Branch(neutrinoGenetaStr,&neutrinoGeneta_tmp[igen],neutrinoGenetaStr+"/F");
    easyTree -> Branch(neutrinoGenphiStr,&neutrinoGenphi_tmp[igen],neutrinoGenphiStr+"/F");
    easyTree -> Branch(neutrinoGenpidStr,&neutrinoGenpid_tmp[igen],neutrinoGenpidStr+"/F");
  }
  
*/  
	
  //--------- GEN JETS Information
  int ngjet=4;
  float jetGenpt_tmp[ngjet], jetGeneta_tmp[ngjet], jetGenphi_tmp[ngjet],  jetGenm_tmp[ngjet] ;
  float jetGenAreaX_tmp[ngjet], jetGenAreaY_tmp[ngjet], jetGenAreaZ_tmp[ngjet], jetGenAreaT_tmp[ngjet];
	
  for(int ijet = 0; ijet<ngjet; ijet++){
    TString jetGenptStr = "jetGenpt"; jetGenptStr += (ijet+1);
    TString jetGenetaStr = "jetGeneta"; jetGenetaStr += (ijet+1);
    TString jetGenphiStr = "jetGenphi"; jetGenphiStr += (ijet+1);
    TString jetGenmStr = "jetGenm"; jetGenmStr += (ijet+1);
    TString jetGenAreaxStr = "jetGenAreaX"; jetGenAreaxStr += (ijet+1);
    TString jetGenAreayStr = "jetGenAreaY"; jetGenAreayStr += (ijet+1);
    TString jetGenAreazStr = "jetGenAreaZ"; jetGenAreazStr += (ijet+1);
    TString jetGenAreatStr = "jetGenAreaT"; jetGenAreatStr += (ijet+1);
    
    easyTree -> Branch(jetGenptStr,&jetGenpt_tmp[ijet],jetGenptStr+"/F");
    easyTree -> Branch(jetGenetaStr,&jetGeneta_tmp[ijet],jetGenetaStr+"/F");
    easyTree -> Branch(jetGenphiStr,&jetGenphi_tmp[ijet],jetGenphiStr+"/F");
    easyTree -> Branch(jetGenmStr,&jetGenm_tmp[ijet],jetGenmStr+"/F");
    easyTree -> Branch(jetGenAreaxStr,&jetGenAreaX_tmp[ijet],jetGenAreaxStr+"/F");
    easyTree -> Branch(jetGenAreayStr,&jetGenAreaY_tmp[ijet],jetGenAreayStr+"/F");
    easyTree -> Branch(jetGenAreazStr,&jetGenAreaZ_tmp[ijet],jetGenAreazStr+"/F");
    easyTree -> Branch(jetGenAreatStr,&jetGenAreaT_tmp[ijet],jetGenAreatStr+"/F");
        
  }
	
  //--------- TRACK JETS Information
  int ntjet=8;
	
  float jetTrackpt_tmp[ntjet], jetTracketa_tmp[ntjet], jetTrackphi_tmp[ntjet],  jetTrackm_tmp[ntjet] ;
  float jetTrackAreaX_tmp[ntjet], jetTrackAreaY_tmp[ntjet], jetTrackAreaZ_tmp[ntjet], jetTrackAreaT_tmp[ntjet];
	
	
  for(int ijet = 0; ijet<ntjet; ijet++){
    TString jetTrackptStr = "jetTrackpt"; jetTrackptStr += (ijet+1);
    TString jetTracketaStr = "jetTracketa"; jetTracketaStr += (ijet+1);
    TString jetTrackphiStr = "jetTrackphi"; jetTrackphiStr += (ijet+1);
    TString jetTrackmStr = "jetTrackm"; jetTrackmStr += (ijet+1);
    TString jetTrackAreaxStr = "jetTrackAreaX"; jetTrackAreaxStr += (ijet+1);
    TString jetTrackAreayStr = "jetTrackAreaY"; jetTrackAreayStr += (ijet+1);
    TString jetTrackAreazStr = "jetTrackAreaZ"; jetTrackAreazStr += (ijet+1);
    TString jetTrackAreatStr = "jetTrackAreaT"; jetTrackAreatStr += (ijet+1);
        
    easyTree -> Branch(jetTrackptStr,&jetTrackpt_tmp[ijet],jetTrackptStr+"/F");
    easyTree -> Branch(jetTracketaStr,&jetTracketa_tmp[ijet],jetTracketaStr+"/F");
    easyTree -> Branch(jetTrackphiStr,&jetTrackphi_tmp[ijet],jetTrackphiStr+"/F");
    easyTree -> Branch(jetTrackmStr,&jetTrackm_tmp[ijet],jetTrackmStr+"/F");
    easyTree -> Branch(jetTrackAreaxStr,&jetTrackAreaX_tmp[ijet],jetTrackAreaxStr+"/F");
    easyTree -> Branch(jetTrackAreayStr,&jetTrackAreaY_tmp[ijet],jetTrackAreayStr+"/F");
    easyTree -> Branch(jetTrackAreazStr,&jetTrackAreaZ_tmp[ijet],jetTrackAreazStr+"/F");
    easyTree -> Branch(jetTrackAreatStr,&jetTrackAreaT_tmp[ijet],jetTrackAreatStr+"/F");
        
  }
	
    
  //---------  JETS (JetPUID) Information
  int njet=8;
	
  float jeteta_tmp[njet], jetphi_tmp[njet], jetpt_tmp[njet],  jetmass_tmp[njet] ;
  float jetAreaX_tmp[njet], jetAreaY_tmp[njet], jetAreaZ_tmp[njet], jetAreaT_tmp[njet];
  float jetBTagAlgo_tmp[njet], jetBTagDefault_tmp[njet],jetBTagPhysics_tmp[njet], jetBTagNearest2_tmp[njet], jetBTagNearest3_tmp[njet], jetBTagHeaviest_tmp[njet] ;
  float jetFlavourAlgo_tmp[njet], jetFlavourDefault_tmp[njet],jetFlavourPhysics_tmp[njet], jetFlavourNearest2_tmp[njet], jetFlavourNearest3_tmp[njet], jetFlavourHeaviest_tmp[njet];
  float jetptD_tmp[njet], jetptDNe_tmp[njet], jetptDCh_tmp[njet];
  float jetnNeutral_tmp[njet], jetnCharged_tmp[njet], jetneuEMfrac_tmp[njet], jetneuHadfrac_tmp[njet];
  float jetbetaClassic_tmp[njet],jetbetaClassicStar_tmp[njet], jetbeta_tmp[njet], jetbetaStar_tmp[njet], jetconstituents_tmp[njet], jetaxis2_tmp[njet];
  float mjj_tmp, detajj_tmp; 
  int njet_tmp, njetid_tmp, nbjet_tmp, hardbjpb_tmp, softbjpb_tmp;

  easyTree -> Branch("mjj",&mjj_tmp,"mjj/F");
  easyTree -> Branch("detajj",&detajj_tmp,"detajj/F");
  easyTree -> Branch("njet",&njet_tmp,"njet/I");
  easyTree -> Branch("nbjet",&nbjet_tmp,"nbjet/I");
  easyTree -> Branch("hardbjpb",&hardbjpb_tmp,"hardbjpb/I");    
  easyTree -> Branch("softbjpb",&softbjpb_tmp,"softbjpb/I");    
  easyTree -> Branch("njetid",&njetid_tmp,"njetid/I");
	
  for(int ijet = 0; ijet<njet; ijet++){
    TString jetptStr = "jetpt"; jetptStr += (ijet+1);
    TString jetetaStr = "jeteta"; jetetaStr += (ijet+1);
    TString jetphiStr = "jetphi"; jetphiStr += (ijet+1);
    TString jetmStr = "jetmass"; jetmStr += (ijet+1);
    TString jetAreaxStr = "jetAreaX"; jetAreaxStr += (ijet+1);
    TString jetAreayStr = "jetAreaY"; jetAreayStr += (ijet+1);
    TString jetAreazStr = "jetAreaZ"; jetAreazStr += (ijet+1);
    TString jetAreatStr = "jetAreaT"; jetAreatStr += (ijet+1);
    TString jetBTagAlgoStr = "jetBTagAlgo"; jetBTagAlgoStr += (ijet+1);
    TString jetBTagDefaultStr = "jetBTagDefault"; jetBTagDefaultStr += (ijet+1);
    TString jetBTagPhysicsStr = "jetBTagPhysics"; jetBTagPhysicsStr += (ijet+1);
    TString jetBTagNearest2Str = "jetBTagNearest2_"; jetBTagNearest2Str += (ijet+1);
    TString jetBTagNearest3Str = "jetBTagNearest3_"; jetBTagNearest3Str += (ijet+1);
    TString jetBTagHeaviestStr = "jetBTagHeaviest_"; jetBTagHeaviestStr += (ijet+1);
    TString jetFlavourAlgoStr = "jetFlavourAlgo"; jetFlavourAlgoStr += (ijet+1);
    TString jetFlavourDefaultStr = "jetFlavourDefault"; jetFlavourDefaultStr += (ijet+1);
    TString jetFlavourPhysicsStr = "jetFlavourPhysics"; jetFlavourPhysicsStr += (ijet+1);
    TString jetFlavourNearest2Str = "jetFlavourNearest2_"; jetFlavourNearest2Str += (ijet+1);
    TString jetFlavourNearest3Str = "jetFlavourNearest3_"; jetFlavourNearest3Str += (ijet+1);
    TString jetFlavourHeaviestStr = "jetFlavourHeaviest_"; jetFlavourHeaviestStr += (ijet+1);
    TString jetptDStr = "jetptD"; jetptDStr += (ijet+1);
    TString jetptDNeStr = "jetptDNe"; jetptDNeStr += (ijet+1);
    TString jetptDChStr = "jetptDCh"; jetptDChStr += (ijet+1);
    TString jetnNeutralStr = "jetnNeutral"; jetnNeutralStr += (ijet+1);
    TString jetnChargedStr = "jetnCharged"; jetnChargedStr += (ijet+1);
    TString jetneuEMfracStr = "jetneuEMfrac"; jetneuEMfracStr += (ijet+1);
    TString jetneuHadfracStr = "jetneuHadfrac"; jetneuHadfracStr += (ijet+1);
    TString jetbetaClassicStr = "jetbetaClassic"; jetbetaClassicStr += (ijet+1);
    TString jetbetaClassicStarStr = "jetbetaClassicStar"; jetbetaClassicStarStr += (ijet+1);
    TString jetbetaStr = "jetbeta"; jetbetaStr += (ijet+1);
    TString jetbetaStarStr = "jetbetaStar"; jetbetaStarStr += (ijet+1);
    TString jetconstituentsStr = "jetconstituents"; jetconstituentsStr += (ijet+1);
    TString jetaxis2Str = "jetaxis2_"; jetaxis2Str += (ijet+1);
        
    easyTree -> Branch(jetptStr,&jetpt_tmp[ijet],jetptStr+"/F");
    easyTree -> Branch(jetetaStr,&jeteta_tmp[ijet],jetetaStr+"/F");
    easyTree -> Branch(jetphiStr,&jetphi_tmp[ijet],jetphiStr+"/F");
    easyTree -> Branch(jetmStr,&jetmass_tmp[ijet],jetmStr+"/F");
    easyTree -> Branch(jetAreaxStr,&jetAreaX_tmp[ijet],jetAreaxStr+"/F");
    easyTree -> Branch(jetAreayStr,&jetAreaY_tmp[ijet],jetAreayStr+"/F");
    easyTree -> Branch(jetAreazStr,&jetAreaZ_tmp[ijet],jetAreazStr+"/F");
    easyTree -> Branch(jetAreatStr,&jetAreaT_tmp[ijet],jetAreatStr+"/F");
    easyTree -> Branch(jetBTagAlgoStr,&jetBTagAlgo_tmp[ijet],jetBTagAlgoStr+"/F");
    easyTree -> Branch(jetBTagDefaultStr,&jetBTagDefault_tmp[ijet],jetBTagDefaultStr+"/F");
    easyTree -> Branch(jetBTagPhysicsStr,&jetBTagPhysics_tmp[ijet],jetBTagPhysicsStr+"/F");
    easyTree -> Branch(jetBTagNearest2Str,&jetBTagNearest2_tmp[ijet],jetBTagNearest2Str+"/F");
    easyTree -> Branch(jetBTagNearest3Str,&jetBTagNearest3_tmp[ijet],jetBTagNearest3Str+"/F");
    easyTree -> Branch(jetBTagHeaviestStr,&jetBTagHeaviest_tmp[ijet],jetBTagHeaviestStr+"/F");
    easyTree -> Branch(jetFlavourAlgoStr,&jetFlavourAlgo_tmp[ijet],jetFlavourAlgoStr+"/F");
    easyTree -> Branch(jetFlavourDefaultStr,&jetFlavourDefault_tmp[ijet],jetFlavourDefaultStr+"/F");
    easyTree -> Branch(jetFlavourPhysicsStr,&jetFlavourPhysics_tmp[ijet],jetFlavourPhysicsStr+"/F");
    easyTree -> Branch(jetFlavourNearest2Str,&jetFlavourNearest2_tmp[ijet],jetFlavourNearest2Str+"/F");
    easyTree -> Branch(jetFlavourNearest3Str,&jetFlavourNearest3_tmp[ijet],jetFlavourNearest3Str+"/F");
    easyTree -> Branch(jetFlavourHeaviestStr,&jetFlavourHeaviest_tmp[ijet],jetFlavourHeaviestStr+"/F");
    easyTree -> Branch(jetptDStr,&jetptD_tmp[ijet],jetptDStr+"/F");
    easyTree -> Branch(jetptDNeStr,&jetptDNe_tmp[ijet],jetptDNeStr+"/F");
    easyTree -> Branch(jetptDChStr,&jetptDCh_tmp[ijet],jetptDChStr+"/F");
    easyTree -> Branch(jetnNeutralStr,&jetnNeutral_tmp[ijet],jetnNeutralStr+"/F");
    easyTree -> Branch(jetnChargedStr,&jetnCharged_tmp[ijet],jetnChargedStr+"/F");
    easyTree -> Branch(jetneuEMfracStr,&jetneuEMfrac_tmp[ijet],jetneuEMfracStr+"/F");
    easyTree -> Branch(jetneuHadfracStr,&jetneuHadfrac_tmp[ijet],jetneuHadfracStr+"/F");
    easyTree -> Branch(jetbetaClassicStr,&jetbetaClassic_tmp[ijet],jetbetaClassicStr+"/F");
    easyTree -> Branch(jetbetaClassicStarStr,&jetbetaClassicStar_tmp[ijet],jetbetaClassicStarStr+"/F");
    easyTree -> Branch(jetbetaStr,&jetbeta_tmp[ijet],jetbetaStr+"/F");
    easyTree -> Branch(jetbetaStarStr,&jetbetaStar_tmp[ijet],jetbetaStarStr+"/F");
    easyTree -> Branch(jetconstituentsStr,&jetconstituents_tmp[ijet],jetconstituentsStr+"/F");
    easyTree -> Branch(jetaxis2Str,&jetaxis2_tmp[ijet],jetaxis2Str+"/F");
  }
        
  //---------  PUPPI JETS (JetPUID) Information
  int npujet=8;

  float jeteta_puppi_tmp[npujet], jetphi_puppi_tmp[npujet], jetpt_puppi_tmp[npujet],  jetmass_puppi_tmp[npujet] ;
  float jetAreaX_puppi_tmp[npujet], jetAreaY_puppi_tmp[npujet], jetAreaZ_puppi_tmp[npujet], jetAreaT_puppi_tmp[npujet];
  float jetBTagAlgo_puppi_tmp[npujet], jetBTagDefault_puppi_tmp[npujet],jetBTagPhysics_puppi_tmp[npujet], jetBTagNearest2_puppi_tmp[npujet], jetBTagNearest3_puppi_tmp[npujet], jetBTagHeaviest_puppi_tmp[npujet] ;
  float jetFlavourAlgo_puppi_tmp[npujet], jetFlavourDefault_puppi_tmp[npujet],jetFlavourPhysics_puppi_tmp[npujet], jetFlavourNearest2_puppi_tmp[npujet], jetFlavourNearest3_puppi_tmp[npujet], jetFlavourHeaviest_puppi_tmp[npujet];
  float jetptD_puppi_tmp[npujet], jetptDNe_puppi_tmp[npujet], jetptDCh_puppi_tmp[npujet];
  float jetnNeutral_puppi_tmp[npujet], jetnCharged_puppi_tmp[npujet], jetneuEMfrac_puppi_tmp[npujet], jetneuHadfrac_puppi_tmp[npujet];
  float jetbetaClassic_puppi_tmp[npujet],jetbetaClassicStar_puppi_tmp[npujet], jetbeta_puppi_tmp[npujet], jetbetaStar_puppi_tmp[npujet], jetconstituents_puppi_tmp[npujet], jetaxis2_puppi_tmp[npujet];
  float mjj_puppi_tmp, detajj_puppi_tmp;
  int njet_puppi_tmp, njetid_puppi_tmp, nbjet_puppi_tmp, hardbjpb_puppi_tmp, softbjpb_puppi_tmp;

  easyTree -> Branch("mjj_puppi",&mjj_puppi_tmp,"mjj_puppi/F");
  easyTree -> Branch("detajj_puppi",&detajj_puppi_tmp,"detajj_puppi/F");
  easyTree -> Branch("njet_puppi",&njet_puppi_tmp,"njet_puppi/I");
  easyTree -> Branch("nbjet_puppi",&nbjet_puppi_tmp,"nbjet_puppi/I");
  easyTree -> Branch("hardbjpb_puppi",&hardbjpb_puppi_tmp,"hardbjpb_puppi/I");
  easyTree -> Branch("softbjpb_puppi",&softbjpb_puppi_tmp,"softbjpb_puppi/I");
  easyTree -> Branch("njetid_puppi",&njetid_puppi_tmp,"njetid_puppi/I");

  for(int ijet = 0; ijet<npujet; ijet++){
    TString jetpt_puppiStr = "jetpt_puppi"; jetpt_puppiStr +=  (ijet+1);
    TString jeteta_puppiStr = "jeteta_puppi"; jeteta_puppiStr +=  (ijet+1);
    TString jetphi_puppiStr = "jetphi_puppi"; jetphi_puppiStr +=  (ijet+1);
    TString jetm_puppiStr = "jetmass_puppi"; jetm_puppiStr +=  (ijet+1);
    TString jetAreax_puppiStr = "jetAreaX_puppi"; jetAreax_puppiStr +=  (ijet+1);
    TString jetAreay_puppiStr = "jetAreaY_puppi"; jetAreay_puppiStr +=  (ijet+1);
    TString jetAreaz_puppiStr = "jetAreaZ_puppi"; jetAreaz_puppiStr +=  (ijet+1);
    TString jetAreat_puppiStr = "jetAreaT_puppi"; jetAreat_puppiStr +=  (ijet+1);
    TString jetBTagAlgo_puppiStr = "jetBTagAlgo_puppi"; jetBTagAlgo_puppiStr +=  (ijet+1);
    TString jetBTagDefault_puppiStr = "jetBTagDefault_puppi"; jetBTagDefault_puppiStr +=  (ijet+1);
    TString jetBTagPhysics_puppiStr = "jetBTagPhysics_puppi"; jetBTagPhysics_puppiStr +=  (ijet+1);
    TString jetBTagNearest2_puppiStr = "jetBTagNearest2_puppi"; jetBTagNearest2_puppiStr +=  (ijet+1);
    TString jetBTagNearest3_puppiStr = "jetBTagNearest3_puppi"; jetBTagNearest3_puppiStr +=  (ijet+1);
    TString jetBTagHeaviest_puppiStr = "jetBTagHeaviest_puppi"; jetBTagHeaviest_puppiStr +=  (ijet+1);
    TString jetFlavourAlgo_puppiStr = "jetFlavourAlgo_puppi"; jetFlavourAlgo_puppiStr +=  (ijet+1);
    TString jetFlavourDefault_puppiStr = "jetFlavourDefault_puppi"; jetFlavourDefault_puppiStr +=  (ijet+1);
    TString jetFlavourPhysics_puppiStr = "jetFlavourPhysics_puppi"; jetFlavourPhysics_puppiStr +=  (ijet+1);
    TString jetFlavourNearest2_puppiStr = "jetFlavourNearest2_puppi"; jetFlavourNearest2_puppiStr +=  (ijet+1);
    TString jetFlavourNearest3_puppiStr = "jetFlavourNearest3_puppi"; jetFlavourNearest3_puppiStr +=  (ijet+1);
    TString jetFlavourHeaviest_puppiStr = "jetFlavourHeaviest_puppi"; jetFlavourHeaviest_puppiStr +=  (ijet+1);
    TString jetptD_puppiStr = "jetptD_puppi"; jetptD_puppiStr +=  (ijet+1);
    TString jetptDNe_puppiStr = "jetptDNe_puppi"; jetptDNe_puppiStr +=  (ijet+1);
    TString jetptDCh_puppiStr = "jetptDCh_puppi"; jetptDCh_puppiStr +=  (ijet+1);
    TString jetnNeutral_puppiStr = "jetnNeutral_puppi"; jetnNeutral_puppiStr +=  (ijet+1);
    TString jetnCharged_puppiStr = "jetnCharged_puppi"; jetnCharged_puppiStr +=  (ijet+1);
    TString jetneuEMfrac_puppiStr = "jetneuEMfrac_puppi"; jetneuEMfrac_puppiStr +=  (ijet+1);
    TString jetneuHadfrac_puppiStr = "jetneuHadfrac_puppi"; jetneuHadfrac_puppiStr +=  (ijet+1);
    TString jetbetaClassic_puppiStr = "jetbetaClassic_puppiStr"; jetbetaClassic_puppiStr +=  (ijet+1);
    TString jetbetaClassicStar_puppiStr = "jetbetaClassicStar_puppiStr"; jetbetaClassicStar_puppiStr +=  (ijet+1);
    TString jetbeta_puppiStr = "jetbeta_puppiStr"; jetbeta_puppiStr +=  (ijet+1);
    TString jetbetaStar_puppiStr = "jetbetaStar_puppiStr"; jetbetaStar_puppiStr +=  (ijet+1);
    TString jetconstituents_puppiStr = "jetconstituents_puppiStr"; jetconstituents_puppiStr +=  (ijet+1);
    TString jetaxis2_puppiStr = "jetaxis2_puppiStr"; jetaxis2_puppiStr +=  (ijet+1);
    
    
    easyTree -> Branch(jetpt_puppiStr,&jetpt_puppi_tmp[ijet],jetpt_puppiStr+"/F");
    easyTree -> Branch(jeteta_puppiStr,&jeteta_puppi_tmp[ijet],jeteta_puppiStr+"/F");
    easyTree -> Branch(jetphi_puppiStr,&jetphi_puppi_tmp[ijet],jetphi_puppiStr+"/F");
    easyTree -> Branch(jetm_puppiStr,&jetmass_puppi_tmp[ijet],jetm_puppiStr+"/F");
    easyTree -> Branch(jetAreax_puppiStr,&jetAreaX_puppi_tmp[ijet],jetAreax_puppiStr+"/F");
    easyTree -> Branch(jetAreay_puppiStr,&jetAreaY_puppi_tmp[ijet],jetAreay_puppiStr+"/F");
    easyTree -> Branch(jetAreaz_puppiStr,&jetAreaZ_puppi_tmp[ijet],jetAreaz_puppiStr+"/F");
    easyTree -> Branch(jetAreat_puppiStr,&jetAreaT_puppi_tmp[ijet],jetAreat_puppiStr+"/F");
    easyTree -> Branch(jetBTagAlgo_puppiStr,&jetBTagAlgo_puppi_tmp[ijet],jetBTagAlgo_puppiStr+"/F");
    easyTree -> Branch(jetBTagDefault_puppiStr,&jetBTagDefault_puppi_tmp[ijet],jetBTagDefault_puppiStr+"/F");
    easyTree -> Branch(jetBTagPhysics_puppiStr,&jetBTagPhysics_puppi_tmp[ijet],jetBTagPhysics_puppiStr+"/F");
    easyTree -> Branch(jetBTagNearest2_puppiStr,&jetBTagNearest2_puppi_tmp[ijet],jetBTagNearest2_puppiStr+"/F");
    easyTree -> Branch(jetBTagNearest3_puppiStr,&jetBTagNearest3_puppi_tmp[ijet],jetBTagNearest3_puppiStr+"/F");
    easyTree -> Branch(jetBTagHeaviest_puppiStr,&jetBTagHeaviest_puppi_tmp[ijet],jetBTagHeaviest_puppiStr+"/F");
    easyTree -> Branch(jetFlavourAlgo_puppiStr,&jetFlavourAlgo_puppi_tmp[ijet],jetFlavourAlgo_puppiStr+"/F");
    easyTree -> Branch(jetFlavourDefault_puppiStr,&jetFlavourDefault_puppi_tmp[ijet],jetFlavourDefault_puppiStr+"/F");
    easyTree -> Branch(jetFlavourPhysics_puppiStr,&jetFlavourPhysics_puppi_tmp[ijet],jetFlavourPhysics_puppiStr+"/F");
    easyTree -> Branch(jetFlavourNearest2_puppiStr,&jetFlavourNearest2_puppi_tmp[ijet],jetFlavourNearest2_puppiStr+"/F");
    easyTree -> Branch(jetFlavourNearest3_puppiStr,&jetFlavourNearest3_puppi_tmp[ijet],jetFlavourNearest3_puppiStr+"/F");
    easyTree -> Branch(jetFlavourHeaviest_puppiStr,&jetFlavourHeaviest_puppi_tmp[ijet],jetFlavourHeaviest_puppiStr+"/F");
    easyTree -> Branch(jetptD_puppiStr,&jetptD_puppi_tmp[ijet],jetptD_puppiStr+"/F");
    easyTree -> Branch(jetptDNe_puppiStr,&jetptDNe_puppi_tmp[ijet],jetptDNe_puppiStr+"/F");
    easyTree -> Branch(jetptDCh_puppiStr,&jetptDCh_puppi_tmp[ijet],jetptDCh_puppiStr+"/F");
    easyTree -> Branch(jetnNeutral_puppiStr,&jetnNeutral_puppi_tmp[ijet],jetnNeutral_puppiStr+"/F");
    easyTree -> Branch(jetnCharged_puppiStr,&jetnCharged_puppi_tmp[ijet],jetnCharged_puppiStr+"/F");
    easyTree -> Branch(jetneuEMfrac_puppiStr,&jetneuEMfrac_puppi_tmp[ijet],jetneuEMfrac_puppiStr+"/F");
    easyTree -> Branch(jetneuHadfrac_puppiStr,&jetneuHadfrac_puppi_tmp[ijet],jetneuHadfrac_puppiStr+"/F");
    easyTree -> Branch(jetbetaClassic_puppiStr,&jetbetaClassic_puppi_tmp[ijet],jetbetaClassic_puppiStr+"/F");
    easyTree -> Branch(jetbetaClassicStar_puppiStr,&jetbetaClassicStar_puppi_tmp[ijet],jetbetaClassicStar_puppiStr+"/F");
    easyTree -> Branch(jetbeta_puppiStr,&jetbeta_puppi_tmp[ijet],jetbeta_puppiStr+"/F");
    easyTree -> Branch(jetbetaStar_puppiStr,&jetbetaStar_puppi_tmp[ijet],jetbetaStar_puppiStr+"/F");
    easyTree -> Branch(jetconstituents_puppiStr,&jetconstituents_puppi_tmp[ijet],jetconstituents_puppiStr+"/F");
    easyTree -> Branch(jetaxis2_puppiStr,&jetaxis2_puppi_tmp[ijet],jetaxis2_puppiStr+"/F");
  }


  //--------- LEPTON Branches

  int nlep=4;
  int nextra_tmp, sameflav_tmp, nlepton_tmp;
  int channel_tmp;	//0 mumu, 1 elel, 2 elmu, 3 muel
  float mll_tmp,  PTll_tmp, dPhill_tmp, dRll_tmp, dEtall_tmp, etall_tmp, yll_tmp;

  float pt_tmp[nlep], eta_tmp[nlep], phi_tmp[nlep], iso_tmp[nlep], isoDBeta_tmp[nlep], isoRhoCorr_tmp[nlep] ;
  float sumChargedHadron_tmp[nlep], sumNeutral_tmp[nlep], sumChargedPU_tmp[nlep], sumAllParticles_tmp[nlep];
  double  ch_tmp[nlep];

  easyTree -> Branch("mll",&mll_tmp,"mll/F");
  easyTree -> Branch("PTll",&PTll_tmp,"PTll/F");
  easyTree -> Branch("dPhill",&dPhill_tmp,"dPhill/F");
  easyTree -> Branch("dRll",&dRll_tmp,"dRll/F");
  easyTree -> Branch("dEtall",&dEtall_tmp,"dEtall/F");
  easyTree -> Branch("etall",&etall_tmp,"etall/F");
  easyTree -> Branch("yll",&yll_tmp,"yll/F");
  easyTree -> Branch("nextra",&nextra_tmp,"nextra/I");
  easyTree -> Branch("nlepton",&nlepton_tmp,"nlepton/I");
  easyTree -> Branch("sameflav",&sameflav_tmp,"sameflav/I");
  easyTree -> Branch("channel",&channel_tmp,"channel/I");

  for(int ilep =0; ilep<nlep; ilep++){
    TString ptStr = "pt"; ptStr += (ilep+1);
    TString etaStr = "eta"; etaStr += (ilep+1);
    TString phiStr = "phi"; phiStr += (ilep+1);
    TString chStr = "ch"; chStr += (ilep+1);
    TString isoStr = "iso"; isoStr += (ilep+1);
    TString isoDBetaStr = "isoDBeta"; isoDBetaStr += (ilep+1);
    TString isoRhoCorrStr = "isoRhoCorr"; isoRhoCorrStr += (ilep+1);
    TString sumChargedHadronStr = "sumChargedHadron"; sumChargedHadronStr += (ilep+1);
    TString sumNeutralStr = "sumNeutral"; sumNeutralStr += (ilep+1);
    TString sumChargedPUStr = "sumChargedPU"; sumChargedPUStr += (ilep+1);
    TString sumAllParticlesStr = "sumAllParticles"; sumAllParticlesStr += (ilep+1);
	    
    easyTree -> Branch(ptStr,&pt_tmp[ilep],ptStr+"/F");
    easyTree -> Branch(etaStr,&eta_tmp[ilep],etaStr+"/F");
    easyTree -> Branch(phiStr,&phi_tmp[ilep],phiStr+"/F");
    easyTree -> Branch(chStr,&ch_tmp[ilep],chStr+"/I");
    easyTree -> Branch(isoStr,&iso_tmp[ilep],isoStr+"/F");
    easyTree -> Branch(isoDBetaStr,&isoDBeta_tmp[ilep],isoDBetaStr+"/F");
    easyTree -> Branch(isoRhoCorrStr,&isoRhoCorr_tmp[ilep],isoRhoCorrStr+"/F");
    easyTree -> Branch(sumChargedHadronStr,&sumChargedHadron_tmp[ilep],sumChargedHadronStr+"/F");
    easyTree -> Branch(sumNeutralStr,&sumNeutral_tmp[ilep],sumNeutralStr+"/F");
    easyTree -> Branch(sumChargedPUStr,&sumChargedPU_tmp[ilep],sumChargedPUStr+"/F");
    easyTree -> Branch(sumAllParticlesStr,&sumAllParticles_tmp[ilep],sumAllParticlesStr+"/F");
  }


  //	MET
    
  float pfmet_tmp, pfmetphi_tmp;
  easyTree -> Branch("pfmet",&pfmet_tmp,"pfmet/F");
  easyTree -> Branch("pfmetphi",&pfmetphi_tmp,"pfmetphi/F");
	
  //	GENMET
    
  float metGenpt_tmp, metGenphi_tmp;
  easyTree -> Branch("metGenpt",&metGenpt_tmp,"metGenpt/F");
  easyTree -> Branch("metGenphi",&metGenphi_tmp,"metGenphi/F");
    
  //	Puppi MET
    
  float pfmet_puppi_tmp, pfmetphi_puppi_tmp;
  easyTree -> Branch("pfmet_puppi",&pfmet_puppi_tmp,"pfmet_puppi/F");
  easyTree -> Branch("pfmetphi_puppi",&pfmetphi_puppi_tmp,"pfmetphi_puppi/F");
  
  //  NPU
  
  float npu_tmp;
  easyTree->Branch("npu",&npu_tmp,"npu/F");
  
  //  Global RhoKt4
  
  float  globalRhokt4_tmp ;
  easyTree->Branch("globalRhokt4",&globalRhokt4_tmp,"globalRhokt4/F");
  
  //  Global RhoGridFastJet
  
  float  globalRhoGridFastJet_tmp ;
  easyTree->Branch("globalGridFastJet",&globalRhoGridFastJet_tmp,"globalGridFastJet/F");
 
  
  //  RhoKt4
  // 0: eta = 0 and eta = 2.5 range
  // 1: eta = 2.5 and eta = 4 range
  // 2: eta = 4 and eta = 5 range

  float Rhokt4_tmp[3] ;
  easyTree->Branch("Rhokt4_0",&Rhokt4_tmp[0],"Rhokt4_0/F");
  easyTree->Branch("Rhokt4_1",&Rhokt4_tmp[1],"Rhokt4_1/F");
  easyTree->Branch("Rhokt4_2",&Rhokt4_tmp[2],"Rhokt4_2/F");
  
  //  RhoGridFastJet
  // 0: eta = 0 and eta = 2.5 range
  // 1: eta = 2.5 and eta = 4 range
  // 2: eta = 4 and eta = 5 range
  float RhoGridFastJet_tmp[3] ;
  easyTree->Branch("RhoGridFastJet_0",&RhoGridFastJet_tmp[0],"RhoGridFastJet_0/F");
  easyTree->Branch("RhoGridFastJet_1",&RhoGridFastJet_tmp[1],"RhoGridFastJet_1/F");
  easyTree->Branch("RhoGridFastJet_2",&RhoGridFastJet_tmp[2],"RhoGridFastJet_2/F");
  
  //  PuppiRhoKt4
  // 0: eta = 0 and eta = 2.5 range
  // 1: eta = 2.5 and eta = 5 range
  float PuppiRhokt4_tmp[2] ;
  easyTree->Branch("PuppiRhokt4_0",&PuppiRhokt4_tmp[0],"PuppiRhokt4_0/F");
  easyTree->Branch("PuppiRhokt4_1",&PuppiRhokt4_tmp[1],"PuppiRhokt4_1/F");

  //  RhoGridFastJet
  // 0: eta = 0 and eta = 2.5 range
  // 1: eta = 2.5 and eta = 5 range
  float PuppiRhoGridFastJet_tmp[2] ;
  easyTree->Branch("PuppiRhoGridFastJet_0",&PuppiRhoGridFastJet_tmp[0],"PuppiRhoGridFastJet_0/F");
  easyTree->Branch("PuppiRhoGridFastJet_1",&PuppiRhoGridFastJet_tmp[1],"PuppiRhoGridFastJet_1/F");

  

    
    
  for(int iEvent = 0; iEvent <  numberOfEntries; iEvent++)
    {
      if (iEvent % 1000 == 0)
        {
	  cout << "iEvent = " << iEvent << endl;
        }
      delphesTree -> ReadEntry(iEvent);
        
      //------ LHE Particle Branches -----------------//
        
        
        
      vector< int> lhepartonID;
      vector< int> lhegluonID;
      vector< int> lheleptonID;
      vector< int> lheneutrinoID;
      vector< int> lhevbosonID;
        
      vector<LHEParticle*> lheParton;
      vector<LHEParticle*> lheGluon;
      vector<LHEParticle*> lheLepton;
      vector<LHEParticle*> lheNeutrino;
      vector<LHEParticle*> lheVBoson;
        
      for(int k =0; k<4; k++){
	leptonLHEpt_tmp[k]=-999;	leptonLHEeta_tmp[k]=-999;	leptonLHEphi_tmp[k]=-999;	 leptonLHEpid_tmp[k]=-999;	leptonLHEch_tmp[k]=-999;	leptonLHEm_tmp[k]=-999;
	neutrinoLHEpt_tmp[k]=-999; neutrinoLHEeta_tmp[k]=-999;  neutrinoLHEphi_tmp[k]=-999;  neutrinoLHEpid_tmp[k]=-999;
	jetLHEPartonpt_tmp[k]=-999; jetLHEPartoneta_tmp[k]=-999; jetLHEPartonphi_tmp[k]=-999; jetLHEPartonpid_tmp[k]=-999;
	jetLHEGluonpt_tmp[k]=-999;  jetLHEGluoneta_tmp[k]=-999;  jetLHEGluonphi_tmp[k]=-999;  jetLHEGluonpid_tmp[k]=-999;
	vbosonLHEpt_tmp[k]=-999;  vbosonLHEeta_tmp[k]=-999;  vbosonLHEphi_tmp[k]=-999;  vbosonLHEpid_tmp[k]=-999; vbosonLHEch_tmp[k]=-999;  vbosonLHEm_tmp[k]=-999;
      }
        
        
      int lhe_entries = branchLHEParticle->GetEntriesFast();
        
        
        
      for (int i = 0 ; i < lhe_entries  ; i++) {
	LHEParticle *lhepart = (LHEParticle*) branchLHEParticle->At(i);
	int type   =  lhepart-> PID;
	int status =  lhepart-> Status;
            
	// no incoming particle
	if(status == -1) continue;
            
            
	if (type < 6 && type > -6) {
	  lhepartonID.push_back(i);
	  lheParton.push_back(lhepart);
	}
            
	if (type == 21) {
	  lhegluonID.push_back(i);
	  lheGluon.push_back(lhepart);
	}
            
	if (type == 24 || type ==-24 || type == 23 ) {
	  lhevbosonID.push_back(i);
	  lheVBoson.push_back(lhepart);
	}
            
	if (type == 11 || type == 13 || type == 15 ||type == -11 || type == -13 || type == -15 ) {
	  lheleptonID.push_back(i);
	  lheLepton.push_back(lhepart);
	}
            
	if (type == 12 || type == 14 || type == 16 ||type == -12 || type == -14 || type == -16 ) {
	  lheneutrinoID.push_back(i);
	  lheNeutrino.push_back(lhepart);
	}
      } //lheloop
      
            
      // sorting in PT
            
      sort(lheParton.begin(), lheParton.end(),lheParticleDescendingPt());
      sort(lheLepton.begin(), lheLepton.end(),lheParticleDescendingPt());
      sort(lheNeutrino.begin(), lheNeutrino.end(),lheParticleDescendingPt());
      sort(lheGluon.begin(), lheGluon.end(),lheParticleDescendingPt());
      sort(lheVBoson.begin(), lheVBoson.end(),lheParticleDescendingPt());
            
      int jl = (lheleptonID.size()<4) ? lheleptonID.size():4;
      int jp = (lhepartonID.size()<4) ? lhepartonID.size():4;
      int jn = (lheneutrinoID.size()<4) ? lheneutrinoID.size():4;
      int jg = (lhegluonID.size()<2) ? lhegluonID.size():2;
      int jb = (lhevbosonID.size()<2) ? lhevbosonID.size():2;
            
      for(int j=0; j<jp; j++){
	jetLHEPartonpt_tmp[j] = lheParton.at(j)->PT;
	jetLHEPartoneta_tmp[j] = lheParton.at(j)->Eta;
	jetLHEPartonphi_tmp[j] = lheParton.at(j)->Phi;
	jetLHEPartonpid_tmp[j] = lheParton.at(j)->PID;
                
      }
            
      for(int j=0; j<jg; j++){
	jetLHEGluonpt_tmp[j] = lheGluon.at(j)->PT;
	jetLHEGluoneta_tmp[j] = lheGluon.at(j)->Eta;
	jetLHEGluonphi_tmp[j] = lheGluon.at(j)->Phi;
	jetLHEGluonpid_tmp[j] = lheGluon.at(j)->PID;
                
      }

      for(int j=0; j<jl; j++){
	leptonLHEpt_tmp[j] = lheLepton.at(j)->PT;
	leptonLHEeta_tmp[j] = lheLepton.at(j)->Eta;
	leptonLHEphi_tmp[j] = lheLepton.at(j)->Phi;
	leptonLHEpid_tmp[j] = lheLepton.at(j)->PID;
	leptonLHEch_tmp[j] = lheLepton.at(j)->Charge;
	leptonLHEm_tmp[j] = lheLepton.at(j)->Mass;
                
      }
	
	

            
      for(int j=0; j<jb; j++){
	vbosonLHEpt_tmp[j] = lheVBoson.at(j)->PT;
	vbosonLHEeta_tmp[j] = lheVBoson.at(j)->Eta;
	vbosonLHEphi_tmp[j] = lheVBoson.at(j)->Phi;
	vbosonLHEpid_tmp[j] = lheVBoson.at(j)->PID;
	vbosonLHEch_tmp[j] = lheVBoson.at(j)->Charge;
	vbosonLHEm_tmp[j] = lheVBoson.at(j)->Mass;
                
      }
			
      for(int j=0; j<jn; j++){
	neutrinoLHEpt_tmp[j] = lheNeutrino.at(j)->PT;
	neutrinoLHEeta_tmp[j] = lheNeutrino.at(j)->Eta;
	neutrinoLHEphi_tmp[j] = lheNeutrino.at(j)->Phi;
	neutrinoLHEpid_tmp[j] = lheNeutrino.at(j)->PID;
      }
            
      //--------- GenParticle filling

/*

      vector< int> leptonID;
      vector< int> neutrinoID;
      vector<GenParticle*> genLepton;
      vector<GenParticle*> genNeutrino;
	
      int gen_entries = branchGenParticle->GetEntries();
	
      for(int k =0; k<4;k++){
	
	leptonGenpt_tmp[k]=-99;
	leptonGeneta_tmp[k]=-99;
	leptonGenphi_tmp[k]=-99;
	leptonGenpid_tmp[k]=-99;
	neutrinoGenpt_tmp[k]=-99;
	neutrinoGeneta_tmp[k]=-99;
	neutrinoGenphi_tmp[k]=-99;
	neutrinoGenpid_tmp[k]=-99;
      }

	
	for (int i = 0 ; i < gen_entries  ; i++) {
	GenParticle *part = (GenParticle*) branchGenParticle->At(i);
	int type   =  part-> PID;
		

	if(part->Status!=1)continue;
			
	if (type == 11 || type == 13 || type == 15 ||type == -11 || type == -13 || type == -15 ) {
	  leptonID.push_back(i);
	  genLepton.push_back(part);
	}
			
	if (type == 12 || type == 14 || type == 16 ||type == -12 || type == -14 || type == -16 ) {
	  neutrinoID.push_back(i);
	  genNeutrino.push_back(part);
	}

      }


      sort(genLepton.begin(), genLepton.end(),genParticleDescendingPt());
      sort(genNeutrino.begin(), genNeutrino.end(),genParticleDescendingPt());
			
			

      int jgl = (leptonID.size()<4) ? leptonID.size():4;
      int jgn = (neutrinoID.size()<4) ? neutrinoID.size():4;
			
     
			
      for(int j=0; j<jgl; j++){
	leptonGenpt_tmp[j] = genLepton.at(j)->PT;
	leptonGeneta_tmp[j] = genLepton.at(j)->Eta;
	leptonGenphi_tmp[j] = genLepton.at(j)->Phi;
	leptonGenpid_tmp[j] = genLepton.at(j)->PID;

      }
			
      for(int j=0; j<jgn; j++){
	neutrinoGenpt_tmp[j] = genNeutrino.at(j)->PT;
	neutrinoGeneta_tmp[j] = genNeutrino.at(j)->Eta;
	neutrinoGenphi_tmp[j] = genNeutrino.at(j)->Phi;
	neutrinoGenpid_tmp[j] = genNeutrino.at(j)->PID;
      }
*/			
               
      //------ GEN JET Filling  -----------------//
        
      vector <Jet*> genJet;
		
      for(int k =0; k<4; k++){
	jetGenpt_tmp[k]=-999;
	jetGeneta_tmp[k]=-999;
	jetGenphi_tmp[k]=-999;
	jetGenm_tmp[k]=-999;
	jetGenAreaX_tmp[k]=-999;
	jetGenAreaY_tmp[k]=-999;
	jetGenAreaZ_tmp[k]=-999;
	jetGenAreaT_tmp[k]=-999;
      }
		
      int gjet_entries = branchGenJet->GetEntriesFast();
		
      for (int i = 0 ; i < gjet_entries  ; i++) {
	Jet *genjet = (Jet*) branchGenJet->At(i);
	genJet.push_back(genjet);
      }
						
			
      int njetsgen = (genJet.size()<4) ? genJet.size():4;
			
      for(int j=0; j<njetsgen; j++){
	jetGenpt_tmp[j] = genJet.at(j)->PT;
	jetGeneta_tmp[j] = genJet.at(j)->Eta;
	jetGenphi_tmp[j] = genJet.at(j)->Phi;
	jetGenm_tmp[j] = genJet.at(j)->Mass;
	jetGenAreaX_tmp[j] = genJet.at(j)->AreaX;
	jetGenAreaY_tmp[j] = genJet.at(j)->AreaY;
	jetGenAreaZ_tmp[j] = genJet.at(j)->AreaZ;
	jetGenAreaT_tmp[j] = genJet.at(j)->AreaT;
      }
            
			
      
        
        
      //------ TRACK JET Filling  -----------------//
        
      vector <Jet*> trackJet;
      vector <int> trackJetIndex;
		
    
		
      for(int k =0; k<8; k++){
	jetTrackpt_tmp[k]=-999;
	jetTracketa_tmp[k]=-999;
	jetTrackphi_tmp[k]=-999;
	jetTrackm_tmp[k]=-999;
	jetTrackAreaX_tmp[k]=-999;
	jetTrackAreaY_tmp[k]=-999;
	jetTrackAreaZ_tmp[k]=-999;
	jetTrackAreaT_tmp[k]=-999;
      }
		
      int tjet_entries = branchTrackJet->GetEntriesFast();
		
      for (int i = 0 ; i < tjet_entries  ; i++) {
	Jet *trackjet = (Jet*) branchTrackJet->At(i);
	trackJet.push_back(trackjet);
	trackJetIndex.push_back(i);
      }
			
						
      int njetstrack = (trackJetIndex.size()<8) ? trackJetIndex.size():8;
            
			
      for(int j=0; j<njetstrack; j++){
	jetTrackpt_tmp[j] = trackJet.at(j)->PT;
	jetTracketa_tmp[j] = trackJet.at(j)->Eta;
	jetTrackphi_tmp[j] = trackJet.at(j)->Phi;
	jetTrackm_tmp[j] = trackJet.at(j)->Mass;
	jetTrackAreaX_tmp[j] = trackJet.at(j)->AreaX;
	jetTrackAreaY_tmp[j] = trackJet.at(j)->AreaY;
	jetTrackAreaZ_tmp[j] = trackJet.at(j)->AreaZ;
	jetTrackAreaT_tmp[j] = trackJet.at(j)->AreaT;
                
      }
            
			
      
        
        
      //------  JET (JetPUID) Filling  -----------------//

      vector <Jet*> puidJet;
      vector <int> puidJetIndex;
        
      TLorentzVector jet1, jet2;
      int inbjet = 0, injetid = 0;
		
      mjj_tmp=-999.; detajj_tmp=-999.; 
      njet_tmp=-999; njetid_tmp=-999; nbjet_tmp=-999; hardbjpb_tmp=-999; softbjpb_tmp=-999;
        
      for(int k =0; k<8; k++){
	jetpt_tmp[k]=-999;   jeteta_tmp[k]=-999;  jetphi_tmp[k]=-999;    jetmass_tmp[k]=-999;
            
	jetAreaX_tmp[k]=-999;   jetAreaY_tmp[k]=-999;   jetAreaZ_tmp[k]=-999;    jetAreaT_tmp[k]=-999;
	jetBTagAlgo_tmp[k]=-999;  jetBTagDefault_tmp[k]=-999;  jetBTagPhysics_tmp[k]=-999;
	jetBTagNearest2_tmp[k]=-999;  jetBTagNearest3_tmp[k]=-999;  jetBTagHeaviest_tmp[k]=-999;  
	jetFlavourAlgo_tmp[k]=-999;   jetFlavourDefault_tmp[k]=-999; jetFlavourPhysics_tmp[k]=-999;
	jetFlavourNearest2_tmp[k]=-999;  jetFlavourNearest3_tmp[k]=-999;  jetFlavourHeaviest_tmp[k]=-999;  
            
	jetptD_tmp[k]=-999;  jetptDNe_tmp[k]=-999; jetptDCh_tmp[k]=-999;
	jetnNeutral_tmp[k]=-999;  jetnCharged_tmp[k]=-999;  jetneuEMfrac_tmp[k]=-999;  jetneuHadfrac_tmp[k]=-999;
	jetbetaClassic_tmp[k]=-999; jetbetaClassicStar_tmp[k]=-999;    jetbeta_tmp[k]=-999;  jetbetaStar_tmp[k]=-999;
	jetconstituents_tmp[k]=-999;  jetaxis2_tmp[k]=-999;
      }
        
      int jet_entries = branchJet->GetEntriesFast();
      njet_tmp = jet_entries;

		
      for (int i = 0 ; i < jet_entries  ; i++) {
	Jet *puidjet = (Jet*) branchJet->At(i);
	puidJet.push_back(puidjet);
	puidJetIndex.push_back(i);
			
	if (puidjet->PT > 30) injetid++;
	njetid_tmp = injetid;
			
	// using medium b tagging
	// 0 no btag, 1 loose, 1 medium, 3 tight

	if((puidjet->BTagAlgo == 2 ) &&  puidjet->PT >10 && puidjet->PT <30 )softbjpb_tmp = 1;
	if((puidjet->BTagAlgo == 2) &&  puidjet->PT >30 ){
	  hardbjpb_tmp = 1;
	  inbjet++;
	}
	nbjet_tmp = inbjet;

      }
			
      int njetspuid = (puidJetIndex.size()<8) ? puidJetIndex.size():8;
            
			
      for(int j=0; j<njetspuid; j++){
	jetpt_tmp[j] = puidJet.at(j)->PT;
	jeteta_tmp[j] = puidJet.at(j)->Eta;
	jetphi_tmp[j] = puidJet.at(j)->Phi;
	jetmass_tmp[j] = puidJet.at(j)->Mass;
	jetAreaX_tmp[j] = puidJet.at(j)->AreaX;
	jetAreaY_tmp[j] = puidJet.at(j)->AreaY;
	jetAreaZ_tmp[j] = puidJet.at(j)->AreaZ;
	jetAreaT_tmp[j] = puidJet.at(j)->AreaT;
	jetBTagAlgo_tmp[j] = puidJet.at(j)->BTagAlgo;
	jetBTagDefault_tmp[j] = puidJet.at(j)->BTagDefault;
	jetBTagPhysics_tmp[j] = puidJet.at(j)->BTagPhysics;
	jetBTagNearest2_tmp[j] = puidJet.at(j)->BTagNearest2;
	jetBTagNearest3_tmp[j] = puidJet.at(j)->BTagNearest3;
	jetBTagHeaviest_tmp[j] = puidJet.at(j)->BTagHeaviest;
	jetFlavourAlgo_tmp[j] = puidJet.at(j)->flavourAlgo;
	jetFlavourDefault_tmp[j] = puidJet.at(j)->flavourDefault;
	jetFlavourPhysics_tmp[j] = puidJet.at(j)->flavourPhysics;
	jetFlavourNearest2_tmp[j] = puidJet.at(j)->flavourNearest2;
	jetFlavourNearest3_tmp[j] = puidJet.at(j)->flavourNearest3;
	jetFlavourHeaviest_tmp[j] = puidJet.at(j)->flavourHeaviest;
	jetptD_tmp[j] = puidJet.at(j)->ptD;
	jetptDNe_tmp[j] = puidJet.at(j)->ptDNe;
	jetptDCh_tmp[j] = puidJet.at(j)->ptDCh;
	jetnNeutral_tmp[j] = puidJet.at(j)->nNeutral;
	jetnCharged_tmp[j] = puidJet.at(j)->nCharged;
	jetneuEMfrac_tmp[j] = puidJet.at(j)->neuEMfrac;
	jetneuHadfrac_tmp[j] = puidJet.at(j)->neuHadfrac;
	jetbetaClassic_tmp[j] = puidJet.at(j)->betaClassic;
	jetbetaClassicStar_tmp[j] = puidJet.at(j)->betaClassicStar;
	jetbeta_tmp[j] = puidJet.at(j)->beta;
	jetbetaStar_tmp[j] = puidJet.at(j)->betaStar;
	jetconstituents_tmp[j] = puidJet.at(j)->constituents;
      }
			 
      if (njetspuid >= 2){
	jet1.SetPtEtaPhiM(jetpt_tmp[0], jeteta_tmp[0],jetphi_tmp[0],jetmass_tmp[0]);
	jet2.SetPtEtaPhiM(jetpt_tmp[1], jeteta_tmp[1],jetphi_tmp[1],jetmass_tmp[1]);
	detajj_tmp = fabs(jeteta_tmp[0] -jeteta_tmp[1]);
	mjj_tmp = (jet1 + jet2).M(); 
      }
            
			
      
        
        
      //------  JET (Jetpuppi) Filling  -----------------//

      vector <Jet*> puppiJet;
      vector <int> puppiJetIndex;

      TLorentzVector jet1puppi, jet2puppi;
      int inbjetpuppi = 0, injetidpuppi = 0;

      mjj_puppi_tmp=-999.; detajj_puppi_tmp=-999.;
      njet_puppi_tmp=-999; njetid_puppi_tmp=-999; nbjet_puppi_tmp=-999; hardbjpb_puppi_tmp=-999; softbjpb_puppi_tmp=-999;

      for(int k =0; k<8; k++){
	jetpt_puppi_tmp[k]=-999;   jeteta_puppi_tmp[k]=-999;  jetphi_puppi_tmp[k]=-999;    jetmass_puppi_tmp[k]=-999;
    
	jetAreaX_puppi_tmp[k]=-999;   jetAreaY_puppi_tmp[k]=-999;   jetAreaZ_puppi_tmp[k]=-999;    jetAreaT_puppi_tmp[k]=-999;
	jetBTagAlgo_puppi_tmp[k]=-999;  jetBTagDefault_puppi_tmp[k]=-999;  jetBTagPhysics_puppi_tmp[k]=-999;
	jetBTagNearest2_puppi_tmp[k]=-999;  jetBTagNearest3_puppi_tmp[k]=-999;  jetBTagHeaviest_puppi_tmp[k]=-999;
	jetFlavourAlgo_puppi_tmp[k]=-999;   jetFlavourDefault_puppi_tmp[k]=-999; jetFlavourPhysics_puppi_tmp[k]=-999;
	jetFlavourNearest2_puppi_tmp[k]=-999;  jetFlavourNearest3_puppi_tmp[k]=-999;  jetFlavourHeaviest_puppi_tmp[k]=-999;
    
	jetptD_puppi_tmp[k]=-999;  jetptDNe_puppi_tmp[k]=-999; jetptDCh_puppi_tmp[k]=-999;
	jetnNeutral_puppi_tmp[k]=-999;  jetnCharged_puppi_tmp[k]=-999;  jetneuEMfrac_puppi_tmp[k]=-999;  jetneuHadfrac_puppi_tmp[k]=-999;
	jetbetaClassic_puppi_tmp[k]=-999; jetbetaClassicStar_puppi_tmp[k]=-999;    jetbeta_puppi_tmp[k]=-999;  jetbetaStar_puppi_tmp[k]=-999;
	jetconstituents_puppi_tmp[k]=-999;  jetaxis2_puppi_tmp[k]=-999;
      }

      int puppijet_entries = branchPuppiJet->GetEntriesFast();
      njet_puppi_tmp = puppijet_entries;


      for (int i = 0 ; i < puppijet_entries  ; i++) {
	Jet *puppijet = (Jet*) branchPuppiJet->At(i);
	puppiJet.push_back(puppijet);
	puppiJetIndex.push_back(i);
    
	if (puppijet->PT > 30) injetidpuppi++;
	njetid_puppi_tmp = injetidpuppi;
    
	// using medium b tagging
	// 0 no btag, 1 loose, 1 medium, 3 tight
    
	if((puppijet->BTagAlgo == 2 ) &&  puppijet->PT >10 && puppijet->PT <30 )softbjpb_puppi_tmp = 1;
        if((puppijet->BTagAlgo == 2) &&  puppijet->PT >30 ){
	  hardbjpb_puppi_tmp = 1;
	  inbjetpuppi++;
        }

      }
      nbjet_puppi_tmp = inbjetpuppi;
    
    
      int njetspuppi = (puppiJetIndex.size()<8) ? puppiJetIndex.size():8;
    
    
      for(int j=0; j<njetspuppi; j++){
	jetpt_puppi_tmp[j] = puppiJet.at(j)->PT;
	jeteta_puppi_tmp[j] = puppiJet.at(j)->Eta;
	jetphi_puppi_tmp[j] = puppiJet.at(j)->Phi;
	jetmass_puppi_tmp[j] = puppiJet.at(j)->Mass;
	jetAreaX_puppi_tmp[j] = puppiJet.at(j)->AreaX;
	jetAreaY_puppi_tmp[j] = puppiJet.at(j)->AreaY;
	jetAreaZ_puppi_tmp[j] = puppiJet.at(j)->AreaZ;
	jetAreaT_puppi_tmp[j] = puppiJet.at(j)->AreaT;
	jetBTagAlgo_puppi_tmp[j] = puppiJet.at(j)->BTagAlgo;
	jetBTagDefault_puppi_tmp[j] = puppiJet.at(j)->BTagDefault;
	jetBTagPhysics_puppi_tmp[j] = puppiJet.at(j)->BTagPhysics;
	jetBTagNearest2_puppi_tmp[j] = puppiJet.at(j)->BTagNearest2;
	jetBTagNearest3_puppi_tmp[j] = puppiJet.at(j)->BTagNearest3;
	jetBTagHeaviest_puppi_tmp[j] = puppiJet.at(j)->BTagHeaviest;
	jetFlavourAlgo_puppi_tmp[j] = puppiJet.at(j)->flavourAlgo;
	jetFlavourDefault_puppi_tmp[j] = puppiJet.at(j)->flavourDefault;
	jetFlavourPhysics_puppi_tmp[j] = puppiJet.at(j)->flavourPhysics;
	jetFlavourNearest2_puppi_tmp[j] = puppiJet.at(j)->flavourNearest2;
	jetFlavourNearest3_puppi_tmp[j] = puppiJet.at(j)->flavourNearest3;
	jetFlavourHeaviest_puppi_tmp[j] = puppiJet.at(j)->flavourHeaviest;
	jetptD_puppi_tmp[j] = puppiJet.at(j)->ptD;
	jetptDNe_puppi_tmp[j] = puppiJet.at(j)->ptDNe;
	jetptDCh_puppi_tmp[j] = puppiJet.at(j)->ptDCh;
	jetnNeutral_puppi_tmp[j] = puppiJet.at(j)->nNeutral;
	jetnCharged_puppi_tmp[j] = puppiJet.at(j)->nCharged;
	jetneuEMfrac_puppi_tmp[j] = puppiJet.at(j)->neuEMfrac;
	jetneuHadfrac_puppi_tmp[j] = puppiJet.at(j)->neuHadfrac;
	jetbetaClassic_puppi_tmp[j] = puppiJet.at(j)->betaClassic;
	jetbetaClassicStar_puppi_tmp[j] = puppiJet.at(j)->betaClassicStar;
	jetbeta_puppi_tmp[j] = puppiJet.at(j)->beta;
	jetbetaStar_puppi_tmp[j] = puppiJet.at(j)->betaStar;
	jetconstituents_puppi_tmp[j] = puppiJet.at(j)->constituents;
      }
    
      if (njetspuppi >= 2){
	jet1puppi.SetPtEtaPhiM(jetpt_puppi_tmp[0], jeteta_puppi_tmp[0],jetphi_puppi_tmp[0],jetmass_puppi_tmp[0]);
	jet2puppi.SetPtEtaPhiM(jetpt_puppi_tmp[1], jeteta_puppi_tmp[1],jetphi_puppi_tmp[1],jetmass_puppi_tmp[1]);
	detajj_puppi_tmp = fabs(jeteta_puppi_tmp[0] -jeteta_puppi_tmp[1]);
	mjj_puppi_tmp = (jet1puppi + jet2puppi).M();
      }
    
    
      


      //-------------LEPTON Filling-----------------//

      nextra_tmp=-999; sameflav_tmp=-999; nlepton_tmp=-999;
      channel_tmp=-999;	//0 mumu, 1 elel, 2 elmu, 3 muel
      mll_tmp = -999.;  PTll_tmp = -999. ; dPhill_tmp = -999. ; dRll_tmp = -999. ; dEtall_tmp = -999. ; etall_tmp = -999. ; yll_tmp = -999. ; 


      for(int k =0; k<4; k++){
	pt_tmp[k]=-999.; eta_tmp[k]=-999.;	phi_tmp[k]=-999.; iso_tmp[k]=-999;  
	isoDBeta_tmp[k]=-999;  isoRhoCorr_tmp[k]=-999;   sumChargedHadron_tmp[k]=-999; sumNeutral_tmp[k]=-999; 
	sumChargedPU_tmp[k]=-999;  sumAllParticles_tmp[k]=-999;
      }
	
	
      vector <Lepton> leptonvec;
      TLorentzVector lep1, lep2;

      //default isolation is calculated using DBeta Correction

      for(int i =0; i<branchEl->GetEntries();i++){
	Electron* elec = (Electron*) branchEl->At(i);
	Lepton l;
	l.lpt = elec->PT;
	l.type = 0;
	l.index = i;
	l.leta = elec->Eta;
	l.lphi =elec->Phi;
	l.lch = elec->Charge;
	l.liso =  elec->IsolationVarDBeta;
	l.lisoDBeta = elec->IsolationVarDBeta;
	l.lisoRhoCorr = elec->IsolationVarRhoCorr;
	l.lsumChargedHadron = elec->chargedHadronEnergy;
	l.lsumNeutral = elec->neutralEnergy;
	l.lsumChargedPU =  elec->chargedPUEnergy;
	l.lsumAllParticles = elec->allParticleEnergy;
	leptonvec.push_back(l);
      }

      for(int i =0; i<branchMu->GetEntries();i++){
	Muon* muo = (Muon*) branchMu->At(i);
	Lepton l;
	l.type = 1;
	l.index = i;
	l.lpt = muo->PT;
	l.leta = muo->Eta;
	l.lphi =muo->Phi;
	l.lch = muo->Charge;
	l.liso =  muo->IsolationVarDBeta;
	l.lisoDBeta = muo->IsolationVarDBeta;
	l.lisoRhoCorr = muo->IsolationVarRhoCorr;
	l.lsumChargedHadron = muo->chargedHadronEnergy;
	l.lsumNeutral = muo->neutralEnergy;
	l.lsumChargedPU =  muo->chargedPUEnergy;
	l.lsumAllParticles = muo->allParticleEnergy;
	leptonvec.push_back(l);
      }

	
      // sorting leptons in pt
      sort(leptonvec.begin(),leptonvec.end(),leptonDescendingPt());
	

      for(unsigned int i =0; i<leptonvec.size(); i++){
	pt_tmp[i] = leptonvec.at(i).lpt;
	eta_tmp[i] = leptonvec.at(i).leta;
	phi_tmp[i] = leptonvec.at(i).lphi;
	ch_tmp[i] = leptonvec.at(i).lch;
	iso_tmp[i] = leptonvec.at(i).liso;
	isoDBeta_tmp[i] = leptonvec.at(i).lisoDBeta;
	isoRhoCorr_tmp[i] = leptonvec.at(i).lisoRhoCorr;
	sumChargedHadron_tmp[i] = leptonvec.at(i).lsumChargedHadron;
	sumNeutral_tmp[i] = leptonvec.at(i).lsumNeutral;
	sumChargedPU_tmp[i] = leptonvec.at(i).lsumChargedPU;
	sumAllParticles_tmp[i] = leptonvec.at(i).lsumAllParticles;
      }
	

      nlepton_tmp = leptonvec.size();
	
      //getting nextra leptons
      if(nlepton_tmp >= 2) nextra_tmp = nlepton_tmp - 2;
	


      if(nlepton_tmp >= 2){
	lep1.SetPtEtaPhiM(pt_tmp[0], eta_tmp[0],phi_tmp[0],0);
	lep2.SetPtEtaPhiM(pt_tmp[1], eta_tmp[1],phi_tmp[1],0);
	mll_tmp = (lep1 + lep2).M(); 
	PTll_tmp = pt_tmp[0] + pt_tmp[1];
	dPhill_tmp = DeltaPhi(phi_tmp[0], phi_tmp[1]);
	dRll_tmp = DeltaR(eta_tmp[0], eta_tmp[1], phi_tmp[0], phi_tmp[1]);
	dEtall_tmp = fabs(eta_tmp[0] - eta_tmp[1]);
	etall_tmp = eta_tmp[0] + eta_tmp[1];
	yll_tmp = (lep1 + lep2).Rapidity(); 
	//channel info
	if(leptonvec.at(0).type == 0 && leptonvec.at(1).type ==1 ){
	  channel_tmp = 2;
	  sameflav_tmp =0;
	}
	if(leptonvec.at(0).type == 1 && leptonvec.at(1).type ==0 ){
	  channel_tmp = 3;
	  sameflav_tmp =0;
	}
	if(leptonvec.at(0).type == 0 && leptonvec.at(1).type ==0 ){
	  channel_tmp = 1;
	  sameflav_tmp =1;
	}
	if(leptonvec.at(0).type == 1 && leptonvec.at(1).type ==1 ){
	  channel_tmp = 0;
	  sameflav_tmp =1;
	}
			
      }
		
   
        
      //--------- MET BRANCHES
    	
      pfmet_tmp=-999; pfmetphi_tmp=-999;
      MissingET* met = (MissingET*) branchMET->At(0);
      pfmet_tmp = met->MET;
      pfmetphi_tmp = met->Phi;
		
      //--------- GEN MET BRANCHES
      metGenpt_tmp=-999;  metGenphi_tmp=-999;
      MissingET* genmet = (MissingET*) branchGenMET->At(0);
      metGenpt_tmp = genmet->MET;
      metGenphi_tmp = genmet->Phi;
        
      //--------- MET BRANCHES
    	
      pfmet_puppi_tmp=-999; pfmetphi_puppi_tmp=-999;
      MissingET* puppimet = (MissingET*) branchPuppiMET->At(0);
      pfmet_puppi_tmp = puppimet->MET;
      pfmetphi_puppi_tmp = puppimet->Phi;
      
      //  NPU  Filling
      npu_tmp=-999;
      ScalarHT* scht = (ScalarHT*)branchNPU->At(0);
      npu_tmp = scht->HT;
       
      //  Global RhoKt4 Filling
      globalRhokt4_tmp =-999.;
      Rho* grhokt4 = (Rho*)branchGlobalRhokt4->At(0);
      globalRhokt4_tmp = grhokt4->Rho; 
	
      //  Global RhoGridFastJet Filling
      globalRhoGridFastJet_tmp = -999.;
      Rho* grogfj = (Rho*)branchGlobalRhoGFJ->At(0);
      globalRhoGridFastJet_tmp = grogfj->Rho;
 
      //  RhoKt4 Filling

      int rk_entries = branchRhokt4->GetEntriesFast();

      for(int i=0; i<rk_entries; i++){
	Rho* rho = (Rho*)branchRhokt4->At(i);
	if(rho->Edges[0] == 0 && rho->Edges[1] == 2.5){
	  Rhokt4_tmp[0] = rho->Rho;
	}
	if(rho->Edges[0] == 2.5 && rho->Edges[1] == 4){
	  Rhokt4_tmp[1] = rho->Rho;
	}
	if(rho->Edges[0] == 4 && rho->Edges[1] == 5){
	  Rhokt4_tmp[2] = rho->Rho;
	}
      }

      int rgk_entries = branchRhoGFJ->GetEntriesFast();

      for(int i=0; i<rgk_entries; i++){
	Rho* rho = (Rho*)branchRhoGFJ->At(i);
	if(rho->Edges[0] == 0 && rho->Edges[1] == 2.5){
	  RhoGridFastJet_tmp[0] = rho->Rho;
	}
	if(rho->Edges[0] == 2.5 && rho->Edges[1] == 4){
	  RhoGridFastJet_tmp[1] = rho->Rho;
	}
	if(rho->Edges[0] == 4 && rho->Edges[1] == 5){
	  RhoGridFastJet_tmp[2] = rho->Rho;
	}
      }

      // Puppi Rhokt4 filling

      int rpk_entries = branchPuppiRhokt4->GetEntriesFast();

      for(int i=0; i<rpk_entries; i++){
	Rho* rho = (Rho*)branchPuppiRhokt4->At(i);
	if(rho->Edges[0] == 0 && rho->Edges[1] == 2.5){
	  PuppiRhokt4_tmp[0] = rho->Rho;
	}
	if(rho->Edges[0] == 2.5 && rho->Edges[1] == 5){
	  PuppiRhokt4_tmp[1] = rho->Rho;
	}
		
      }

      //Puppi GridFastJetRho filling

      int rpgk_entries = branchPuppiRhoGFJ->GetEntriesFast();

      for(int i=0; i<rpgk_entries; i++){
	Rho* rho = (Rho*)branchPuppiRhoGFJ->At(i);
	if(rho->Edges[0] == 0 && rho->Edges[1] == 2.5){
	  PuppiRhoGridFastJet_tmp[0] = rho->Rho;
	}
	if(rho->Edges[0] == 2.5 && rho->Edges[1] == 5){
	  PuppiRhoGridFastJet_tmp[1] = rho->Rho;
	}
		
      }

        
      easyTree -> Fill();
		
    }
    
	
	
  //easyTree -> Print("easyDelphes");
  outputFile -> Write();
  delete outputFile;
}










