
/*c++ -O2 -lm `root-config --cflags --glibs` -L /afs/cern.ch/work/d/depoyraz/VBS/TP/CMSSW_6_2_0_SLHC19/src/Delphes_Two/Delphes -I /afs/cern.ch/work/d/depoyraz/VBS/TP/CMSSW_6_2_0_SLHC19/src/Delphes_Two/Delphes -l Delphes -o DelphesDumper DelphesDumper.cpp
  ./DelphesDumper file.txt output.root
  file: list of the delphes trees
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
#include "TLorentzVector.h"
#include "TMath.h"
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
using namespace std;

//****************************************************************************************
struct eventType
{
  vector<long int> eventID;
  vector<int> eventChannel;
  // event channel: 0 DF, 1 ElEl, 2MuMu
};

//****************************************************************************************
// runs over all entries in the delphes tree 
// no inital cuts applied for the event selector



eventType event_preselector(ExRootTreeReader *delphesTree, TClonesArray* branchEl ,TClonesArray* branchMu)
{
  cout << endl << "################## EVENTS PRESELECTION ##################" << endl;
  eventType goodEvent;
  int iEvent=0;
  for(iEvent = 0; iEvent < delphesTree->GetEntries(); iEvent++)
    {
      if (iEvent % 10000 == 0)
	{
	  cout << "iEvent = " << iEvent << endl;
	}
      delphesTree -> ReadEntry(iEvent);
      // select DF events
      if(branchEl->GetEntriesFast() >= 1 && branchMu->GetEntriesFast() >= 1)
	{
	  goodEvent.eventID.push_back(iEvent);
	  goodEvent.eventChannel.push_back(0);
	}
      // select SF events
      if(branchMu->GetEntriesFast() >= 2 && branchEl->GetEntriesFast() == 0)
	{
	  goodEvent.eventID.push_back(iEvent);
	  goodEvent.eventChannel.push_back(2);
	}
	
      if(branchEl->GetEntriesFast() >= 2 && branchMu->GetEntriesFast() == 0)
	{
	  goodEvent.eventID.push_back(iEvent);
	  goodEvent.eventChannel.push_back(1);
	}
	
    }
  cout << "######### events from delphes: " << iEvent << endl;
  cout << "######### events after preselection cuts: "<<goodEvent.eventID.size()<<endl;
  return goodEvent;
}

//****************************************************************************************

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
  vector<string> inputFiles;
  TChain* delphesNtuples = new TChain("Delphes");
  ExRootTreeReader *delphesTree = new ExRootTreeReader(delphesNtuples);
  TFile* outputFile = TFile::Open(argv[2],"recreate");
  TTree* easyTree = new TTree("easyDelphes","easyDelphes");
  //----------------------------------------------------------------------------------------
  // reading input files
  if(argc < 3)
    {
      cout << "ERROR: not enough info provided" << endl;
      return 0;
    }
  ifstream inputList (argv[1], ios::in);
  string buffer;
  while(inputList >> buffer)
    {
      inputFiles.push_back(buffer);
      cout << "####### Input file #" << inputFiles.size() << ": " << inputFiles.back() << endl;
    }
  //--------- filling the TChain
  for (int iFiles = 0; iFiles < (int)inputFiles.size(); iFiles++)
    {
      delphesNtuples -> Add((inputFiles.at(iFiles)).c_str());
    }
  delphesNtuples -> BranchRef();
  //----------------------------------------------------------------------------------------
  //variable management
  //
  // -> all objects in the output tree are stored in decreasing pt order.

  //--------- getting objects from the delphes tree

  TClonesArray* branchEl = delphesTree->UseBranch("Electron");
  TClonesArray* branchMu = delphesTree->UseBranch("Muon");
  TClonesArray* branchJet = delphesTree->UseBranch("Jet");
  TClonesArray* branchGenParticle = delphesTree->UseBranch("Particle");
  TClonesArray* branchGenJet = delphesTree->UseBranch("GenJet");
  TClonesArray* branchMET = delphesTree->UseBranch("MissingET");
  TClonesArray* branchPuppiMET = delphesTree->UseBranch("PuppiMissingET");
  TClonesArray* branchGenMET = delphesTree->UseBranch("GenMissingET");
  TClonesArray* branchPuppiJet = delphesTree->UseBranch("PuppiJet");
    
  //--------- creating branches for the new (light) tree

  //--------- Gen Particles

  float jetGenPartonpt_tmp[3];
  float jetGenPartonpid_tmp[3];   
  float jetGenPartonphi_tmp[3];
  float jetGenPartoneta_tmp[3];
  float leptonGenpt_tmp[3];
  float leptonGenpid_tmp[3];
  float leptonGenphi_tmp[3];
  float leptonGeneta_tmp[3];
  float neutrinoGenpt_tmp[3];
  float neutrinoGenpid_tmp[3];
  float neutrinoGenphi_tmp[3];
  float neutrinoGeneta_tmp[3];
    
    
  easyTree -> Branch("jetGenPartonpt1",&jetGenPartonpt_tmp[0],"jetGenPartonpt1/F");
  easyTree -> Branch("jetGenPartonpt2",&jetGenPartonpt_tmp[1],"jetGenPartonpt1/F");
  easyTree -> Branch("jetGenPartonpt3",&jetGenPartonpt_tmp[2],"jetGenPartonpt3/F");
  easyTree -> Branch("jetGenPartonpt4",&jetGenPartonpt_tmp[3],"jetGenPartonpt4/F");
  easyTree -> Branch("jetGenPartonpid1",&jetGenPartonpid_tmp[0],"jetGenPartonpid1/F");
  easyTree -> Branch("jetGenPartonpid2",&jetGenPartonpid_tmp[1],"jetGenPartonpid2/F");
  easyTree -> Branch("jetGenPartonpid3",&jetGenPartonpid_tmp[2],"jetGenPartonpid3/F");
  easyTree -> Branch("jetGenPartonpid4",&jetGenPartonpid_tmp[3],"jetGenPartonpid4/F");
  easyTree -> Branch("jetGenPartonphi1",&jetGenPartonphi_tmp[0],"jetGenPartonphi1/F");
  easyTree -> Branch("jetGenPartonphi2",&jetGenPartonphi_tmp[1],"jetGenPartonphi2/F");    
  easyTree -> Branch("jetGenPartonphi3",&jetGenPartonphi_tmp[2],"jetGenPartonphi3/F");
  easyTree -> Branch("jetGenPartonphi4",&jetGenPartonphi_tmp[3],"jetGenPartonphi4/F");
  easyTree -> Branch("jetGenPartoneta1",&jetGenPartoneta_tmp[0],"jetGenPartoneta1/F");
  easyTree -> Branch("jetGenPartoneta2",&jetGenPartoneta_tmp[1],"jetGenPartoneta2/F");    
  easyTree -> Branch("jetGenPartoneta3",&jetGenPartoneta_tmp[2],"jetGenPartoneta3/F");
  easyTree -> Branch("jetGenPartoneta4",&jetGenPartoneta_tmp[3],"jetGenPartoneta4/F");
    
    
  easyTree -> Branch("leptonGenpt1",&leptonGenpt_tmp[0],"leptonGenpt1/F");
  easyTree -> Branch("leptonGenpt2",&leptonGenpt_tmp[1],"leptonGenpt2/F");
  easyTree -> Branch("leptonGenpt3",&leptonGenpt_tmp[2],"leptonGenpt3/F");
  easyTree -> Branch("leptonGenpt4",&leptonGenpt_tmp[3],"leptonGenpt4/F");
  easyTree -> Branch("leptonGenpid1",&leptonGenpid_tmp[0],"leptonGenpid1/F");
  easyTree -> Branch("leptonGenpid2",&leptonGenpid_tmp[1],"leptonGenpid2/F");
  easyTree -> Branch("leptonGenpid3",&leptonGenpid_tmp[2],"leptonGenpid3/F");
  easyTree -> Branch("leptonGenpid4",&leptonGenpid_tmp[3],"leptonGenpid4/F");
  easyTree -> Branch("leptonGenphi1",&leptonGenphi_tmp[0],"leptonGenphi1/F");
  easyTree -> Branch("leptonGenphi2",&leptonGenphi_tmp[1],"leptonGenphi2/F");
  easyTree -> Branch("leptonGenphi3",&leptonGenphi_tmp[2],"leptonGenphi3/F");
  easyTree -> Branch("leptonGenphi4",&leptonGenphi_tmp[3],"leptonGenphi4/F");
  easyTree -> Branch("leptonGeneta1",&leptonGeneta_tmp[0],"leptonGeneta1/F");
  easyTree -> Branch("leptonGeneta2",&leptonGeneta_tmp[1],"leptonGeneta2/F");
  easyTree -> Branch("leptonGeneta3",&leptonGeneta_tmp[2],"leptonGeneta3/F");
  easyTree -> Branch("leptonGeneta4",&leptonGeneta_tmp[3],"leptonGeneta4/F");
    
  easyTree -> Branch("neutrinoGenpt1",&neutrinoGenpt_tmp[0],"neutrinoGenpt1/F");
  easyTree -> Branch("neutrinoGenpt2",&neutrinoGenpt_tmp[1],"neutrinoGenpt2/F");
  easyTree -> Branch("neutrinoGenpt3",&neutrinoGenpt_tmp[2],"neutrinoGenpt3/F");
  easyTree -> Branch("neutrinoGenpt4",&neutrinoGenpt_tmp[3],"neutrinoGenpt4/F");
  easyTree -> Branch("neutrinoGenpid1",&neutrinoGenpid_tmp[0],"neutrinoGenpid1/F");
  easyTree -> Branch("neutrinoGenpid2",&neutrinoGenpid_tmp[1],"neutrinoGenpid2/F");
  easyTree -> Branch("neutrinoGenpid3",&neutrinoGenpid_tmp[2],"neutrinoGenpid3/F");
  easyTree -> Branch("neutrinoGenpid4",&neutrinoGenpid_tmp[3],"neutrinoGenpid4/F");
  easyTree -> Branch("neutrinoGenphi1",&neutrinoGenphi_tmp[0],"neutrinoGenphi1/F");
  easyTree -> Branch("neutrinoGenphi2",&neutrinoGenphi_tmp[1],"neutrinoGenphi2/F");
  easyTree -> Branch("neutrinoGenphi3",&neutrinoGenphi_tmp[2],"neutrinoGenphi3/F");
  easyTree -> Branch("neutrinoGenphi4",&neutrinoGenphi_tmp[3],"neutrinoGenphi4/F");
  easyTree -> Branch("neutrinoGeneta1",&neutrinoGeneta_tmp[0],"neutrinoGeneta1/F");
  easyTree -> Branch("neutrinoGeneta2",&neutrinoGeneta_tmp[1],"neutrinoGeneta2/F");
  easyTree -> Branch("neutrinoGeneta3",&neutrinoGeneta_tmp[2],"neutrinoGeneta3/F");
  easyTree -> Branch("neutrinoGeneta4",&neutrinoGeneta_tmp[3],"neutrinoGeneta4/F");
    
    
    
  //-------------------------------------------------------------------------------------------

  //--------- LHE




  //--------- Leptons
   

   
  int nextra_tmp=-1, sameflav_tmp=-1;;
  int channel_tmp=-1;	//0 mumu, 1 elel, 2 elmu, 3 muel
  float eta1_tmp=-99., eta2_tmp=-99., phi1_tmp=-99., phi2_tmp=-99., ch1_tmp=-99., ch2_tmp=-99., pt1_tmp=-99., pt2_tmp=-99., iso1_tmp=-99., iso2_tmp=-99.;
  int pdgid1_tmp = -99, pdgid2_tmp = -99;
  //float eta3_tmp=-99., eta4_tmp=-99., phi3_tmp=-99., phi4_tmp=-99., ch3_tmp=-99., ch4_tmp=-99., pt3_tmp=-99., pt4_tmp=-99., iso3_tmp=-99., iso4_tmp=-99.;
  //int pdgid3_tmp = -99, pdgid4_tmp = -99;
  float mll_tmp = -1.,  PTll_tmp = -1. , dPhill_tmp = -1. , dRll_tmp = -1. , dEtall_tmp = -1. , etall_tmp = -1. , yll_tmp = -1. ;
   

  vector<float> lepton_pt;
  TLorentzVector lep1, lep2;
    
  easyTree -> Branch("pt1",&pt1_tmp,"pt1/F");
  easyTree -> Branch("pt2",&pt2_tmp,"pt/F"); 
  easyTree -> Branch("eta1",&eta1_tmp,"eta1/F");
  easyTree -> Branch("eta2",&eta2_tmp,"eta2/F");
  easyTree -> Branch("phi1",&phi1_tmp,"phi1/F");
  easyTree -> Branch("phi2",&phi2_tmp,"phi2/F");
  easyTree -> Branch("ch1",&ch1_tmp,"ch1/F");
  easyTree -> Branch("ch2",&ch2_tmp,"ch2/F");
  easyTree -> Branch("iso1",&iso1_tmp,"iso1/F");
  easyTree -> Branch("iso2",&iso2_tmp,"iso2/F");
  easyTree -> Branch("pdgid1",&pdgid1_tmp,"pdgid1/I");
  easyTree -> Branch("pdgid2",&pdgid2_tmp,"pdgid2/I"); 
   
  easyTree -> Branch("mll",&mll_tmp,"mll/F");
  easyTree -> Branch("PTll",&PTll_tmp,"PTll/F");
  easyTree -> Branch("dPhill",&dPhill_tmp,"dPhill/F");
  easyTree -> Branch("dRll",&dRll_tmp,"dRll/F");  
  easyTree -> Branch("dEtall",&dEtall_tmp,"dEtall/F"); 
  easyTree -> Branch("etall",&etall_tmp,"etall/F"); 
  easyTree -> Branch("yll",&yll_tmp,"yll/F");     
  easyTree -> Branch("nextra",&nextra_tmp,"nextra/I");
  easyTree -> Branch("sameflav",&sameflav_tmp,"sameflav/I"); 
  easyTree -> Branch("channel",&channel_tmp,"channel/I");    
    
  /*    easyTree -> Branch("pt3",&pt3_tmp,"pt3/F");
	easyTree -> Branch("pt4",&pt4_tmp,"p4/F"); 
	easyTree -> Branch("eta3",&eta3_tmp,"eta3/F");
	easyTree -> Branch("eta4",&eta4_tmp,"eta4/F");
	easyTree -> Branch("phi3",&phi3_tmp,"phi3/F");
	easyTree -> Branch("phi4",&phi4_tmp,"phi4/F");
	easyTree -> Branch("ch3",&ch3_tmp,"ch3/F");
	easyTree -> Branch("ch4",&ch4_tmp,"ch4/F");
	easyTree -> Branch("iso3",&iso3_tmp,"iso3/F");
	easyTree -> Branch("iso4",&iso4_tmp,"iso4/F");
	easyTree -> Branch("pdgid3",&pdgid4_tmp,"pdgid4/I");
	easyTree -> Branch("pdgid4",&pdgid4_tmp,"pdgid4/I"); 
  */
    
  //----------------------------------------------------------------------------------------



  //---------  GENJets

  int ngenjet_tmp, ngenjetid_tmp, numbergenjet_tmp;

  float jetGenpt_tmp[3];
  float jetGenphi_tmp[3];
  float jetGeneta_tmp[3];
   
    
    
  easyTree -> Branch("jetGenpt1",&jetGenpt_tmp[0],"jetGenpt1/F");
  easyTree -> Branch("jetGenpt2",&jetGenpt_tmp[1],"jetGenpt1/F");
  easyTree -> Branch("jetGenpt3",&jetGenpt_tmp[2],"jetGenpt3/F");
  easyTree -> Branch("jetGenpt4",&jetGenpt_tmp[3],"jetGenpt4/F");
  easyTree -> Branch("jetGenphi1",&jetGenphi_tmp[0],"jetGenphi1/F");
  easyTree -> Branch("jetGenphi2",&jetGenphi_tmp[1],"jetGenphi2/F");    
  easyTree -> Branch("jetGenphi3",&jetGenphi_tmp[2],"jetGenphi3/F");
  easyTree -> Branch("jetGenphi4",&jetGenphi_tmp[3],"jetGenphi4/F");
  easyTree -> Branch("jetGeneta1",&jetGeneta_tmp[0],"jetGeneta1/F");
  easyTree -> Branch("jetGeneta2",&jetGeneta_tmp[1],"jetGeneta2/F");    
  easyTree -> Branch("jetGeneta3",&jetGeneta_tmp[2],"jetGeneta3/F");
  easyTree -> Branch("jetGeneta4",&jetGeneta_tmp[3],"jetGeneta4/F");


  //---------  Jets

  int njet_tmp, njetid_tmp, numberjet_tmp, nbjet_tmp, hardbjpb_tmp=-99, softbjpb_tmp=-99;
  float jeteta_tmp[8];
  float jetphi_tmp[8];
  float jetpt_tmp[8];
  float jetmass_tmp[8];
  float mjj_tmp=-99., detajj_tmp=-99.; 
	
  TLorentzVector jet1, jet2;
   
  easyTree -> Branch("njet",&njet_tmp,"njet/I");
  easyTree -> Branch("nbjet",&nbjet_tmp,"nbjet/I");
  easyTree -> Branch("hardbjpb",&hardbjpb_tmp,"hardbjpb/I");    
  easyTree -> Branch("softbjpb",&softbjpb_tmp,"softbjpb/I");    
  easyTree -> Branch("njetid",&njetid_tmp,"njetid/I");
    
  easyTree -> Branch("jeteta1",&jeteta_tmp[0],"jeteta1/F");
  easyTree -> Branch("jeteta2",&jeteta_tmp[1],"jeteta2/F");
  easyTree -> Branch("jeteta3",&jeteta_tmp[2],"jeteta3/F");
  easyTree -> Branch("jeteta4",&jeteta_tmp[3],"jeteta4/F");
  easyTree -> Branch("jeteta5",&jeteta_tmp[4],"jeteta5/F");
  easyTree -> Branch("jeteta6",&jeteta_tmp[5],"jeteta6/F");
  easyTree -> Branch("jeteta7",&jeteta_tmp[6],"jeteta7/F");
  easyTree -> Branch("jeteta8",&jeteta_tmp[7],"jeteta8/F");
    
  easyTree -> Branch("jetphi1",&jetphi_tmp[0],"jetphi1/F");
  easyTree -> Branch("jetphi2",&jetphi_tmp[1],"jetphi2/F");
  easyTree -> Branch("jetphi3",&jetphi_tmp[2],"jetphi3/F");
  easyTree -> Branch("jetphi4",&jetphi_tmp[3],"jetphi4/F");
  easyTree -> Branch("jetphi5",&jetphi_tmp[4],"jetphi5/F");
  easyTree -> Branch("jetphi6",&jetphi_tmp[5],"jetphi6/F");
  easyTree -> Branch("jetphi7",&jetphi_tmp[6],"jetphi7/F");
  easyTree -> Branch("jetphi8",&jetphi_tmp[7],"jetphi8/F");
    
  easyTree -> Branch("jetpt1",&jetpt_tmp[0],"jetpt1/F");
  easyTree -> Branch("jetpt2",&jetpt_tmp[1],"jetpt2/F");
  easyTree -> Branch("jetpt3",&jetpt_tmp[2],"jetpt3/F");
  easyTree -> Branch("jetpt4",&jetpt_tmp[3],"jetpt4/F");
  easyTree -> Branch("jetpt5",&jetpt_tmp[4],"jetpt5/F");
  easyTree -> Branch("jetpt6",&jetpt_tmp[5],"jetpt6/F");
  easyTree -> Branch("jetpt7",&jetpt_tmp[6],"jetpt7/F");
  easyTree -> Branch("jetpt8",&jetpt_tmp[7],"jetpt8/F");
    
  easyTree -> Branch("jetmass1",&jetmass_tmp[0],"jetmass1/F");
  easyTree -> Branch("jetmass2",&jetmass_tmp[1],"jetmass2/F");
  easyTree -> Branch("jetmass3",&jetmass_tmp[2],"jetmass3/F");
  easyTree -> Branch("jetmass4",&jetmass_tmp[3],"jetmass4/F");
  easyTree -> Branch("jetmass5",&jetmass_tmp[4],"jetmass5/F");
  easyTree -> Branch("jetmass6",&jetmass_tmp[5],"jetmass6/F");
  easyTree -> Branch("jetmass7",&jetmass_tmp[6],"jetmass7/F");
  easyTree -> Branch("jetmass8",&jetmass_tmp[7],"jetmass8/F");
    
  // to do ?
  //  easyTree -> Branch("dphilljet",&dphilljet_tmp,"dphilljet/F");
  //  easyTree -> Branch("dphilljetjet",&dphilljetjet_tmp,"dphilljetjet/F");
  
  easyTree -> Branch("mjj",&mjj_tmp,"mjj/F");
  easyTree -> Branch("detajj",&detajj_tmp,"detajj/F");
    
    
    
    
    
  //-------------------------------------------------------------------------------------------

  // Puppi Jets
  int njet_puppi_tmp, njetid_puppi_tmp, numberjet_puppi_tmp, nbjet_puppi_tmp, hardbjpb_puppi_tmp=-99, softbjpb_puppi_tmp=-99;
  float jeteta_puppi_tmp[8];
  float jetphi_puppi_tmp[8];
  float jetpt_puppi_tmp[8];
  float jetmass_puppi_tmp[8];
  float mjj_puppi_tmp=-99., detajj_puppi_tmp=-99.; 
	
  TLorentzVector jet1puppi, jet2puppi;

  easyTree -> Branch("njet_puppi",&njet_puppi_tmp,"njet/I");
  easyTree -> Branch("nbjet_puppi",&nbjet_puppi_tmp,"nbjet/I");
  easyTree -> Branch("hardbjpb_puppi",&hardbjpb_puppi_tmp,"hardbjpb/I");    
  easyTree -> Branch("softbjpb_puppi",&softbjpb_puppi_tmp,"softbjpb/I");    
  easyTree -> Branch("njetid_puppi",&njetid_puppi_tmp,"njetid/I");
    
  easyTree -> Branch("jeteta1_puppi",&jeteta_puppi_tmp[0],"jeteta1/F");
  easyTree -> Branch("jeteta2_puppi",&jeteta_puppi_tmp[1],"jeteta2/F");
  easyTree -> Branch("jeteta3_puppi",&jeteta_puppi_tmp[2],"jeteta3/F");
  easyTree -> Branch("jeteta4_puppi",&jeteta_puppi_tmp[3],"jeteta4/F");
  easyTree -> Branch("jeteta5_puppi",&jeteta_puppi_tmp[4],"jeteta5/F");
  easyTree -> Branch("jeteta6_puppi",&jeteta_puppi_tmp[5],"jeteta6/F");
  easyTree -> Branch("jeteta7_puppi",&jeteta_puppi_tmp[6],"jeteta7/F");
  easyTree -> Branch("jeteta8_puppi",&jeteta_puppi_tmp[7],"jeteta8/F");
    
  easyTree -> Branch("jetphi1_puppi",&jetphi_puppi_tmp[0],"jetphi1/F");
  easyTree -> Branch("jetphi2_puppi",&jetphi_puppi_tmp[1],"jetphi2/F");
  easyTree -> Branch("jetphi3_puppi",&jetphi_puppi_tmp[2],"jetphi3/F");
  easyTree -> Branch("jetphi4_puppi",&jetphi_puppi_tmp[3],"jetphi4/F");
  easyTree -> Branch("jetphi5_puppi",&jetphi_puppi_tmp[4],"jetphi5/F");
  easyTree -> Branch("jetphi6_puppi",&jetphi_puppi_tmp[5],"jetphi6/F");
  easyTree -> Branch("jetphi7_puppi",&jetphi_puppi_tmp[6],"jetphi7/F");
  easyTree -> Branch("jetphi8_puppi",&jetphi_puppi_tmp[7],"jetphi8/F");
    
  easyTree -> Branch("jetpt1_puppi",&jetpt_puppi_tmp[0],"jetpt1/F");
  easyTree -> Branch("jetpt2_puppi",&jetpt_puppi_tmp[1],"jetpt2/F");
  easyTree -> Branch("jetpt3_puppi",&jetpt_puppi_tmp[2],"jetpt3/F");
  easyTree -> Branch("jetpt4_puppi",&jetpt_puppi_tmp[3],"jetpt4/F");
  easyTree -> Branch("jetpt5_puppi",&jetpt_puppi_tmp[4],"jetpt5/F");
  easyTree -> Branch("jetpt6_puppi",&jetpt_puppi_tmp[5],"jetpt6/F");
  easyTree -> Branch("jetpt7_puppi",&jetpt_puppi_tmp[6],"jetpt7/F");
  easyTree -> Branch("jetpt8_puppi",&jetpt_puppi_tmp[7],"jetpt8/F");
    
  easyTree -> Branch("jetmass1_puppi",&jetmass_puppi_tmp[0],"jetmass1/F");
  easyTree -> Branch("jetmass2_puppi",&jetmass_puppi_tmp[1],"jetmass2/F");
  easyTree -> Branch("jetmass3_puppi",&jetmass_puppi_tmp[2],"jetmass3/F");
  easyTree -> Branch("jetmass4_puppi",&jetmass_puppi_tmp[3],"jetmass4/F");
  easyTree -> Branch("jetmass5_puppi",&jetmass_puppi_tmp[4],"jetmass5/F");
  easyTree -> Branch("jetmass6_puppi",&jetmass_puppi_tmp[5],"jetmass6/F");
  easyTree -> Branch("jetmass7_puppi",&jetmass_puppi_tmp[6],"jetmass7/F");
  easyTree -> Branch("jetmass8_puppi",&jetmass_puppi_tmp[7],"jetmass8/F");
    
  // to do ?
  //  easyTree -> Branch("dphilljet_puppi",&dphilljet_puppi_tmp,"dphilljet/F");
  //  easyTree -> Branch("dphilljetjet_puppi",&dphilljetjet_puppi_tmp,"dphilljetjet/F");
  
  easyTree -> Branch("mjj_puppi",&mjj_puppi_tmp,"mjj/F");
  easyTree -> Branch("detajj_puppi",&detajj_puppi_tmp,"detajj/F");
    


  //----------------------------------------------------------------------------------------
   
  //	MET

  float pfmet_tmp=0,pfmetphi_tmp=0;
  easyTree -> Branch("pfmet",&pfmet_tmp,"pfmet/F");
  easyTree -> Branch("pfmetphi",&pfmetphi_tmp,"pfmetphi/F");
	
  //	GENMET

  float metGenpt_tmp=0,metGeneta_tmp=0, metGenphi_tmp=0;
  easyTree -> Branch("metGenpt",&metGenpt_tmp,"metGenpt/F");
  easyTree -> Branch("metGenphi",&metGenphi_tmp,"metGenphi/F");
    
  //	Puppi MET

  float pfmet_puppi_tmp=0,pfmetphi_puppi_tmp=0;
  easyTree -> Branch("pfmet_puppi",&pfmet_puppi_tmp,"pfmet_puppi/F");
  easyTree -> Branch("pfmetphi_puppi",&pfmetphi_puppi_tmp,"pfmetphi_puppi/F");


	
	
  //filling the new (plain) tree
  eventType goodEvent = event_preselector(delphesTree,branchEl,branchMu);
  //  eventType goodEvent = event_preselector(delphesTree,branchEl,branchMu,branchMET);
  cout << endl << "################# TREE CREATION STARTED #################" << endl;
  for(int iEvent = 0; iEvent < (long int)goodEvent.eventID.size(); iEvent++)
    {
      if (iEvent % 1000 == 0)
	{
	  cout << "iEvent = " << iEvent << endl;
	}
      delphesTree -> ReadEntry(goodEvent.eventID.at(iEvent));
	
      //--------- lepton branches
      int nlep =0;
      //--------- DF
      if(goodEvent.eventChannel.at(iEvent) == 0)
	{
	  Electron* el = (Electron*) branchEl->At(0);
	  lepton_pt.push_back(el->PT);
	  Muon* mu = (Muon*) branchMu->At(0);
	  lepton_pt.push_back(mu->PT);
	    
	  if(branchEl->GetEntries() >1) nlep++;
	  if(branchMu->GetEntries() >1) nlep++;
	  nextra_tmp = nlep;
		
	  sort(lepton_pt.begin(), lepton_pt.end());
	  pt1_tmp = lepton_pt[1];
	  pt2_tmp = lepton_pt[0];
	   
	  if(pt1_tmp == el->PT){
	    eta1_tmp = el->Eta;
	    phi1_tmp = el->Phi;
	    ch1_tmp = el->Charge;
	    iso1_tmp = el->IsolationVar;
	    eta2_tmp = mu->Eta;
	    phi2_tmp = mu->Phi;
	    ch2_tmp = mu->Charge;
    	    iso2_tmp = mu->IsolationVar;
	    channel_tmp = 2.;	//channel =2 elmu
	    lep1.SetPtEtaPhiM(pt1_tmp, eta1_tmp,phi1_tmp,0);
	    lep2.SetPtEtaPhiM(pt2_tmp, eta2_tmp,phi2_tmp,0);
	    //lepton mass are set to zero
	    mll_tmp = (lep1 + lep2).M(); 
	    PTll_tmp = pt1_tmp + pt2_tmp;
	    dPhill_tmp = DeltaPhi(phi1_tmp, phi2_tmp);
	    dRll_tmp = DeltaR(eta1_tmp, eta2_tmp, phi1_tmp, phi2_tmp);
	    dEtall_tmp = fabs(eta1_tmp - eta2_tmp);
	    etall_tmp = eta1_tmp + eta2_tmp;
	    yll_tmp = (lep1 + lep2).Rapidity(); 
	    sameflav_tmp =0;
	  }
	    
	  if(pt1_tmp == mu->PT){
	    eta1_tmp = mu->Eta;
	    phi1_tmp = mu->Phi;
	    ch1_tmp = mu->Charge;
	    iso1_tmp = mu->IsolationVar;
	    eta2_tmp = el->Eta;
	    phi2_tmp = el->Phi;
	    ch2_tmp = el->Charge;
	    iso2_tmp = el->IsolationVar;
	    channel_tmp = 3.; //channel =3 muel
	    lep1.SetPtEtaPhiM(pt1_tmp, eta1_tmp,phi1_tmp,0);
	    lep2.SetPtEtaPhiM(pt2_tmp, eta2_tmp,phi2_tmp,0);
	    mll_tmp = (lep1 + lep2).M(); 
	    PTll_tmp = pt1_tmp + pt2_tmp;
	    dPhill_tmp = DeltaPhi(phi1_tmp, phi2_tmp);
	    dRll_tmp = DeltaR(eta1_tmp, eta2_tmp, phi1_tmp, phi2_tmp);
	    dEtall_tmp = fabs(eta1_tmp - eta2_tmp);
	    etall_tmp = eta1_tmp + eta2_tmp;
	    yll_tmp = (lep1 + lep2).Rapidity(); 
	    sameflav_tmp =0;
	  	  
	  }
	    
	   
	}
      lepton_pt.clear();
		
      //--------- SF- ElEl
      if(goodEvent.eventChannel.at(iEvent) == 1)
	{
	  Electron* el1 = (Electron*) branchEl->At(0);
	  pt1_tmp = el1->PT;
	  eta1_tmp = el1->Eta;
	  phi1_tmp = el1->Phi;
	  ch1_tmp = el1->Charge;
	  iso1_tmp = el1->IsolationVar;
	    
	  Electron* el2 = (Electron*) branchEl->At(1);
	  pt2_tmp = el2->PT;
	  eta2_tmp = el2->Eta;
	  phi2_tmp = el2->Phi;
	  ch2_tmp = el2->Charge;
	  iso2_tmp = el2->IsolationVar;
	  channel_tmp = 1.; //channel =1 elel
	    
	  if(branchEl->GetEntries() >2) nlep = branchEl->GetEntries() -2;
	  nextra_tmp = nlep;
	  lep1.SetPtEtaPhiM(pt1_tmp, eta1_tmp,phi1_tmp,0);
	  lep2.SetPtEtaPhiM(pt2_tmp, eta2_tmp,phi2_tmp,0);
	  mll_tmp = (lep1 + lep2).M(); 
	  PTll_tmp = pt1_tmp + pt2_tmp;
	  dPhill_tmp = DeltaPhi(phi1_tmp, phi2_tmp);
	  dRll_tmp = DeltaR(eta1_tmp, eta2_tmp, phi1_tmp, phi2_tmp);
	  dEtall_tmp = fabs(eta1_tmp - eta2_tmp);
	  etall_tmp = eta1_tmp + eta2_tmp;
	  yll_tmp = (lep1 + lep2).Rapidity(); 
	  sameflav_tmp =1;

	}
	
      if(goodEvent.eventChannel.at(iEvent) == 2)
	{
	  //--------- SF- MuMu
	  Muon* mu1 = (Muon*) branchMu->At(0);
	  Muon* mu2 = (Muon*) branchMu->At(1);
	  pt1_tmp = mu1->PT;
	  eta1_tmp = mu1->Eta;
	  phi1_tmp = mu1->Phi;
	  ch1_tmp = mu1->Charge;
	  iso1_tmp =mu1->IsolationVar;

	  pt2_tmp = mu2->PT;
	  eta2_tmp = mu2->Eta;
	  phi2_tmp = mu2->Phi;
	  ch2_tmp = mu2->Charge;
	  iso2_tmp =mu2->IsolationVar;
	  channel_tmp = 0.;	//channel =0 mumu
	    
	  if(branchMu->GetEntries() >2) nlep = branchMu->GetEntries() -2;
	  nextra_tmp = nlep;
	  lep1.SetPtEtaPhiM(pt1_tmp, eta1_tmp,phi1_tmp,0);
	  lep2.SetPtEtaPhiM(pt2_tmp, eta2_tmp,phi2_tmp,0);
	  mll_tmp = (lep1 + lep2).M(); 
	  PTll_tmp = pt1_tmp + pt2_tmp;
	  dPhill_tmp = DeltaPhi(phi1_tmp, phi2_tmp);
	  dRll_tmp = DeltaR(eta1_tmp, eta2_tmp, phi1_tmp, phi2_tmp);
	  dEtall_tmp = fabs(eta1_tmp - eta2_tmp);
	  etall_tmp = eta1_tmp + eta2_tmp;
	  yll_tmp = (lep1 + lep2).Rapidity(); 
	  sameflav_tmp =1;

	}
	
	
      //--------- GenParticle Branches


      vector< int> partonID;
      vector< int> leptonID;
      vector< int> neutrinoID;
    
      vector<GenParticle*> genParton;
      vector<GenParticle*> genLepton;
      vector<GenParticle*> genNeutrino;
	
      struct genParticleDescendingPt 
      {
	bool operator() (GenParticle* a, GenParticle* b) 
	{     
	  return a->PT > b->PT;
	}
      };
    

      int gen_entries = branchGenParticle->GetEntries();
	
      for(int k =0; k<4;k++){
	jetGenPartonpt_tmp[k]=-99;
	jetGenPartoneta_tmp[k]=-99;
	jetGenPartonphi_tmp[k]=-99;
	jetGenPartonpid_tmp[k]=-99;
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
		

	if ((type < 6 && type > -6) || type == 21) {
	  partonID.push_back(i);
	  genParton.push_back(part);
	}
			
	if (type == 11 || type == 13 || type == 15 ||type == -11 || type == -13 || type == -15 ) {
	  leptonID.push_back(i);
	  genLepton.push_back(part);
	}
			
	if (type == 12 || type == 14 || type == 16 ||type == -12 || type == -14 || type == -16 ) {
	  neutrinoID.push_back(i);
	  genNeutrino.push_back(part);
	}

      }

      sort(genParton.begin(), genParton.end(),genParticleDescendingPt());
      sort(genLepton.begin(), genLepton.end(),genParticleDescendingPt());
      sort(genNeutrino.begin(), genNeutrino.end(),genParticleDescendingPt());
			
			
      int jp = (partonID.size()<4) ? partonID.size():4;
      int jl = (leptonID.size()<4) ? leptonID.size():4;
      int jn = (neutrinoID.size()<4) ? neutrinoID.size():4;
			
      for(int j=0; j<jp; j++){
	jetGenPartonpt_tmp[j] = genParton.at(j)->PT;
	jetGenPartoneta_tmp[j] = genParton.at(j)->Eta;
	jetGenPartonphi_tmp[j] = genParton.at(j)->Phi;
	jetGenPartonpid_tmp[j] = genParton.at(j)->PID;
      }
			
      for(int j=0; j<jl; j++){
	leptonGenpt_tmp[j] = genLepton.at(j)->PT;
	leptonGeneta_tmp[j] = genLepton.at(j)->Eta;
	leptonGenphi_tmp[j] = genLepton.at(j)->Phi;
	leptonGenpid_tmp[j] = genLepton.at(j)->PID;

      }
			
      for(int j=0; j<jn; j++){
	neutrinoGenpt_tmp[j] = genNeutrino.at(j)->PT;
	neutrinoGeneta_tmp[j] = genNeutrino.at(j)->Eta;
	neutrinoGenphi_tmp[j] = genNeutrino.at(j)->Phi;
	neutrinoGenpid_tmp[j] = genNeutrino.at(j)->PID;
      }
			
			
      //--------- Gen Jet branches
	
      numbergenjet_tmp = branchGenJet->GetEntriesFast();
      ngenjet_tmp=0;
      ngenjetid_tmp=0;
	
    
      for(int k =0; k<4;k++){
    	jetGenpt_tmp[k]=-99;	
    	jetGeneta_tmp[k]=-99;	
    	jetGenphi_tmp[k]=-99;	
      }
   	
      ngenjet_tmp = numbergenjet_tmp;
	
      int Ngenjets = (numbergenjet_tmp < 4) ? numbergenjet_tmp : 4;
      // taking maximally first 4 jets	
      for(int i=0; i<Ngenjets; i++){ 
	
	Jet* genjet[njetid_tmp];
    	
	for(int ij =0; ij<ngenjetid_tmp;ij++){

	  genjet[ij] = (Jet*) branchGenJet->At(ij);
	  jetGeneta_tmp[ij] = genjet[ij]->Eta;
	  jetGenpt_tmp[ij] = genjet[ij]->PT;
	  jetGenphi_tmp[ij] = genjet[ij]->Phi;
	}
      }
		
	
      //--------- jet branches
	
      numberjet_tmp = branchJet->GetEntriesFast();
      int nbjetcounter=0;
      njet_tmp=0;
      njetid_tmp=0;
      nbjet_tmp=0;
    
      for(int k =0; k<8;k++){
    	jeteta_tmp[k]=-99;	
    	jetpt_tmp[k]=-99;	
    	jetphi_tmp[k]=-99;	
    	jetmass_tmp[k]=-99;	
      }
   	
      njet_tmp = numberjet_tmp;
	
      int Njets = (numberjet_tmp < 8) ? numberjet_tmp : 8;
      // taking maximally first 8 jets	
      for(int i=0; i<Njets; i++){ 
	
	Jet* jet0 = (Jet*) branchJet->At(i);
		
	// using medium b tagging

	if((jet0->BTag & (1 << 1)) &&  jet0->PT >10 && jet0->PT <30 )softbjpb_tmp = 1;
	if((jet0->BTag & (1 << 1)) &&  jet0->PT >30 ){
	  hardbjpb_tmp = 1;
	  nbjetcounter++;
	}
	nbjet_tmp = nbjetcounter;
	// bjets are counted using hard b tag discriminator

		
	if(jet0->PT < 30) continue;
	njetid_tmp++;
	    
	    
    	Jet* jet[njetid_tmp];
    	
	for(int ij =0; ij<njetid_tmp;ij++){

	  jet[ij] = (Jet*) branchJet->At(ij);
	  jeteta_tmp[ij] = jet[ij]->Eta;
	  jetpt_tmp[ij] = jet[ij]->PT;
	  jetphi_tmp[ij] = jet[ij]->Phi;
	  jetmass_tmp[ij] = jet[ij]->Mass;

	  if(njetid_tmp ==2){
	    jet1.SetPtEtaPhiM(jetpt_tmp[0], jeteta_tmp[0],jetphi_tmp[0],0);
	    jet2.SetPtEtaPhiM(jetpt_tmp[1], jeteta_tmp[1],jetphi_tmp[1],0);
	    if (jetpt_tmp[1]==-99 || jeteta_tmp[1]==-99 || jetphi_tmp[1]==-99)  continue;
	    detajj_tmp = fabs(jeteta_tmp[0] -jeteta_tmp[1]);
	    mjj_tmp = (jet1 + jet2).M(); 
	  }
	}
      }
    
      //--------- MET BRANCHES
    	
      MissingET* met = (MissingET*) branchMET->At(0);
      pfmet_tmp = met->MET;
      pfmetphi_tmp = met->Phi;
		
      //--------- GEN MET BRANCHES
      MissingET* genmet = (MissingET*) branchGenMET->At(0);
      metGenpt_tmp = genmet->MET;
      metGenphi_tmp = genmet->Phi;

      //--------- MET BRANCHES
    	
      MissingET* puppimet = (MissingET*) branchPuppiMET->At(0);
      pfmet_puppi_tmp = puppimet->MET;
      pfmetphi_puppi_tmp = puppimet->Phi;
		
		
      //--------- LHEBbranches

	
		
		
      //--------- Puppi Jet Branches
	
      numberjet_puppi_tmp = branchPuppiJet->GetEntriesFast();
      int nbjetcounter_puppi=0;
      njet_puppi_tmp=0;
      njetid_puppi_tmp=0;
      nbjet_puppi_tmp=0;
    
      for(int k =0; k<8;k++){
    	jeteta_puppi_tmp[k]=-99;	
    	jetpt_puppi_tmp[k]=-99;	
    	jetphi_puppi_tmp[k]=-99;	
    	jetmass_puppi_tmp[k]=-99;	
      }
   	
      njet_puppi_tmp = numberjet_puppi_tmp;
	

      int Njets_puppi = (numberjet_puppi_tmp < 8) ? numberjet_puppi_tmp : 8;
      // taking maximally first 8 jets	
      for(int i=0; i<Njets_puppi; i++){ 
	
	Jet* jet0 = (Jet*) branchPuppiJet->At(i);
		
	if((jet0->BTag & (1 << 1)) &&  jet0->PT >30 )hardbjpb_puppi_tmp = 1;
	if((jet0->BTag & (1 << 1)) &&  jet0->PT >10 && jet0->PT <30 )softbjpb_puppi_tmp = 1;
	// not sure about this definition
		
	if(jet0->PT < 30) continue;
	njetid_puppi_tmp++;
	    
	if(jet0->BTag & (1 << 1)) {nbjetcounter_puppi++;}
	nbjet_puppi_tmp = nbjet_puppi_tmp;
	//no b puppi jets?
	    
    	
    	Jet* jet[njetid_puppi_tmp];
    	
	for(int ij =0; ij<njetid_puppi_tmp;ij++){
	    
	  jet[ij] = (Jet*) branchJet->At(ij);
	  jeteta_puppi_tmp[ij] = jet[ij]->Eta;
	  jetpt_puppi_tmp[ij] = jet[ij]->PT;
	  jetphi_puppi_tmp[ij] = jet[ij]->Phi;
	  jetmass_puppi_tmp[ij] = jet[ij]->Mass;
    	    
	  if(njetid_puppi_tmp ==2){
	    jet1.SetPtEtaPhiM(jetpt_puppi_tmp[0], jeteta_puppi_tmp[0],jetphi_puppi_tmp[0],0);
	    jet2.SetPtEtaPhiM(jetpt_puppi_tmp[1], jeteta_puppi_tmp[1],jetphi_puppi_tmp[1],0);
	    if (jetpt_puppi_tmp[1]==-99 || jeteta_puppi_tmp[1]==-99 || jetphi_puppi_tmp[1]==-99)  continue;
	    detajj_puppi_tmp = fabs(jeteta_puppi_tmp[0] -jeteta_puppi_tmp[1]);
	    mjj_puppi_tmp = (jet1 + jet2).M(); 

	  }
			
		
	}
			
			
	        
      }
    

		
      easyTree -> Fill();

	
    }
	

  easyTree -> Print("easyDelphes");
  outputFile -> Write();
  delete outputFile;
}


	





