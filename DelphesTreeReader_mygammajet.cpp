/****************************************************************************************
- this program reads the output from delphes and store the interesting parameters in
a new tree, it also apllies some preselection cuts
- before compiling ---> source ../Decay/setup_slc6.sh
- compile with ---> c++ -O2 -lm `root-config --cflags --glibs` -L /afs/cern.ch/user/s/spigazzi/work/DelphesStuff/delphes_code/ -I /afs/cern.ch/user/s/spigazzi/work/DelphesStuff/delphes_code/ -lDelphes -o delphesTreeReader delphesTreeReader.cpp
*****************************************************************************************/

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
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
#include "TRefArray.h"
using namespace std;

//****************************************************************************************
struct eventType
{
    vector<long int> eventID;
    vector<int> lepFlavor;
};

//****************************************************************************************
// runs over all entries in the delphes tree and applies some preselection cuts
// -> semileptonic final state.

/*
eventType event_preselector(ExRootTreeReader *delphesTree, TClonesArray* branchEl ,TClonesArray* branchMu, TClonesArray* branchMET)
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
	if(branchEl->GetEntriesFast() == 1 && branchMu->GetEntriesFast() == 0)
	{
	    Electron* el = (Electron*) branchEl->At(0);
//MissingET* met = (MissingET*) branchMET->At(0);
	    if(el->PT > 27)// && met->MET > 65)
	    {
		goodEvent.eventID.push_back(iEvent);
		goodEvent.lepFlavor.push_back(11);
	    }
	}
	if(branchMu->GetEntriesFast() == 1 && branchEl->GetEntriesFast() == 0)
	{
	    Muon* mu = (Muon*) branchMu->At(0);
//MissingET* met = (MissingET*) branchMET->At(0);
	    if(mu->PT > 24)// && met->MET > 50)
	    {
		goodEvent.eventID.push_back(iEvent);
		goodEvent.lepFlavor.push_back(13);
	    }
	}
    }
    cout << "######### events from delphes: " << iEvent << endl;
    cout << "######### events after preselection cuts: " << goodEvent.eventID.size() << endl;
    return goodEvent;
}*/
//****************************************************************************************
//delta phi
float DeltaPhi (float phi1, float phi2)
{
float delta_phi = TMath::Abs(phi1 - phi2);
if (delta_phi > 2*TMath::Pi())
delta_phi -= 2*TMath::Pi();
if (delta_phi > TMath::Pi() && delta_phi < 2*TMath::Pi())
delta_phi = 2*TMath::Pi() - delta_phi;
return delta_phi;
}
//****************************************************************************************
//delta R
float DeltaR (float eta1, float eta2, float phi1, float phi2)
{
float d_phi = DeltaPhi(phi1, phi2);
float d_eta = TMath::Abs(eta1 - eta2);
return TMath::Sqrt(pow(d_eta,2)+pow(d_phi,2));
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
// -> jets: jets of both type in the output tree are stored in decreasing pt order.
// -> leptons: the isolation variable is defined as pt_sum(all particle in a cone)/pt_lep
//--------- getting objects from the delphes tree

    //TClonesArray* branch_PvJetPUID = delphesTree->UseBranch("PvJetPUID");
    //TClonesArray* branch_PvJetPUID = delphesTree->UseBranch("mergedjets");
    TClonesArray* branch_PvJetPUID = delphesTree->UseBranch("PvJet");
    //TClonesArray* branch_PvJetPUID = delphesTree->UseBranch("noTimeJet");
    //TClonesArray* branch_JetPUID = delphesTree->UseBranch("JetPUID");
    TClonesArray* branch_JetPUID = delphesTree->UseBranch("Jet");
    TClonesArray* branch_GenJet = delphesTree->UseBranch("GenJet");
    TClonesArray* branch_LHEParticles = delphesTree->UseBranch("LHEParticles");
    
//--------- creating branches for the new (light) tree


    vector<float> JetPUID_Pt ,JetPUID_E,JetPUID_m,JetPUID_Eta,JetPUID_Phi, vbs_cnd_pt, vbs_sec_cnd_pt;
    easyTree -> Branch("JetPUID_Pt","vector<float>",&JetPUID_Pt);
    easyTree -> Branch("JetPUID_Eta","vector<float>",&JetPUID_Eta);
    easyTree -> Branch("JetPUID_Phi","vector<float>",&JetPUID_Phi);
    easyTree -> Branch("JetPUID_m","vector<float>",&JetPUID_m);
    easyTree -> Branch("JetPUID_E","vector<float>",&JetPUID_E);
    
//----------Jet with time cut
    vector<float> PvJetPUID_Pt ,PvJetPUID_E,PvJetPUID_m,PvJetPUID_Eta,PvJetPUID_Phi, Pv_vbs_cnd_pt, Pv_vbs_sec_cnd_pt;
    easyTree -> Branch("PvJetPUID_Pt","vector<float>",&PvJetPUID_Pt);
    easyTree -> Branch("PvJetPUID_Eta","vector<float>",&PvJetPUID_Eta);
    easyTree -> Branch("PvJetPUID_Phi","vector<float>",&PvJetPUID_Phi);
    easyTree -> Branch("PvJetPUID_m","vector<float>",&PvJetPUID_m);
    easyTree -> Branch("PvJetPUID_E","vector<float>",&PvJetPUID_E);
//--------- GenJet
    vector<float> GenJet_Pt ,GenJet_E,GenJet_m,GenJet_Eta,GenJet_Phi, gen_vbs_cnd_pt, gen_vbs_sec_cnd_pt;
    easyTree -> Branch("GenJet_Pt","vector<float>",&GenJet_Pt);
    easyTree -> Branch("GenJet_Eta","vector<float>",&GenJet_Eta);
    easyTree -> Branch("GenJet_Phi","vector<float>",&GenJet_Phi);
    easyTree -> Branch("GenJet_m","vector<float>",&GenJet_m);
    //easyTree -> Branch("GenJet_E","vector<float>",&GenJet_E);
//---------- LHE Particles
    vector<float> part_Pt ,part_E,part_m,part_Eta,part_Phi, part_status;
    vector<int> part_ID;
    easyTree -> Branch("part_Pt","vector<float>",&part_Pt);
    easyTree -> Branch("part_Eta","vector<float>",&part_Eta);
    easyTree -> Branch("part_Phi","vector<float>",&part_Phi);
    easyTree -> Branch("part_m","vector<float>",&part_status);
    easyTree -> Branch("part_PID","vector<int>",&part_ID);
//-----------------------------------------------------
float pv_e_frac, e_frac, d_r, d_r2, e_gen, reco_frac, dR1, dR2, dR3, dR4, pv_dR1, pv_dR2, pv_dR3, pv_dR4; 
vector<float> vbf_etaphi, pv_vbf_etaphi, gen_vbf_etaphi, dR, pv_dR, dr_tmpvector, reco_dR, efrac_jet, pvefrac_jet, ejet_vect, pvejet_vect, genejet_vect;
    easyTree -> Branch("PvEfrac",&pv_e_frac,"PvEfrac/F");    
    easyTree -> Branch("Efrac", &e_frac, "Efrac/F");
    easyTree -> Branch("deltaR",&d_r, "deltaR/F");
   // easyTree -> Branch("deltaR", &d_r2, "deltaR/F");
   easyTree -> Branch("e_gen",&e_gen,"e_gen/F");
   easyTree -> Branch("reco_frac",&reco_frac,"reco_frac/F");
   //easyTree -> Branch("vbf_eta","vector<float",&vbf_eta);
   //easyTree -> Branch("pv_vbf_eta","vector<float>",&pv_vbf_eta);
   //easyTree -> Branch("gen_vbf_eta","vector<float>",&gen_vbf_eta);
   easyTree -> Branch("dR","vector<float>",&dR);
   easyTree -> Branch("pv_dR","vector<float>",&pv_dR);
   easyTree -> Branch("dr_tmpvector","vector<float>",&dr_tmpvector);
   easyTree -> Branch("reco_dR","vector<float>",&gen_vbf_etaphi);
   easyTree -> Branch("dR1",&pv_dR1,"dR1/F");
   easyTree -> Branch("dR2",&pv_dR2,"dR2/F");
   easyTree -> Branch("dR3",&pv_dR1,"dR3/F");
   easyTree -> Branch("dR4",&pv_dR1,"dR4/F");
   easyTree -> Branch("efrac_jet","vector<float>",&efrac_jet);
   easyTree -> Branch("pvefrac_jet","vector<float>",&pvefrac_jet);
   easyTree -> Branch("ejet_vect","vector<float>",&ejet_vect);
   easyTree -> Branch("pvejet_vect","vector<float>",&pvejet_vect);
   easyTree -> Branch("genejet_vect","vector<float>",&genejet_vect);
//----------------------------------------------------------------------------------------
//filling the new (plain) tree
   // eventType goodEvent = event_preselector(delphesTree,branchEl,branchMu,branchMET);
    cout << endl << "################# TREE CREATION STARTED #################" << endl;
    //for(int iEvent = 0; iEvent < 1000000; iEvent++)
    int iEvent = 0;
    for(iEvent = 0; iEvent < delphesTree->GetEntries() ; iEvent++)
    //for(iEvent = 0; iEvent < 249; iEvent++)
    {
	if (iEvent % 1000 == 0)
	{
	    cout << "iEvent = " << iEvent << endl;
	}
	//delphesTree -> ReadEntry(goodEvent.eventID.at(iEvent));
	delphesTree -> ReadEntry(iEvent);
//--------- Jet old school
float Pt_tmp = 0., sum_pt = 0., Pt_tmp2 = 0., sum_pt_tmp = 0., pt1 = 0., pt2 = 0., eta_tmp, eta_tmp2, phi_tmp, phi_tmp2, m_tmp, m_tmp2, E_sum=0., pt_tmp, pt_sum=0., eta1, eta2, phi1, phi2, m1, m2, pt_tmp2;        
int jetID, jetID2;
TLorentzVector temp, temp2, sum;
   int JetPUID_entries = branch_JetPUID->GetEntriesFast(), JetPUID_i = 0;
 //  if(JetPUID_entries != 0) {
//   TLorentzVector temp, temp2, sum;
   while(JetPUID_i < JetPUID_entries)   {
        Jet* JetPUID = (Jet*) branch_JetPUID->At(JetPUID_i);
        JetPUID_Pt.push_back(JetPUID->PT);
        JetPUID_m.push_back(JetPUID->Mass);
        JetPUID_Eta.push_back(JetPUID->Eta);
        JetPUID_Phi.push_back(JetPUID->Phi);
        ++JetPUID_i;
        }

/*float Pt_tmp = 0., sum_pt = 0., Pt_tmp2 = 0., sum_pt_tmp = 0., pt1 = 0., pt2 = 0., eta_tmp, eta_tmp2, phi_tmp, phi_tmp2, m_tmp, m_tmp2, E_sum=0., pt_tmp, pt_sum=0., eta1, eta2, phi1, phi2, m1, m2, pt_tmp2;        
int jetID, jetID2;  */
/*for(unsigned it = 0; it != JetPUID_i; it++)    {
    //float Pt_tmp = *it;
    pt_tmp = JetPUID_Pt.at(it);
    eta_tmp = JetPUID_Eta.at(it);
    phi_tmp = JetPUID_Phi.at(it);
    m_tmp = JetPUID_m.at(it);
    //temp.SetPtEtaPhiM(Pt_tmp,eta_tmp,phi_tmp,m_tmp);
    temp.SetPtEtaPhiM(pt_tmp,eta_tmp,phi_tmp,m_tmp);
for(unsigned it2 =  1; it2 != JetPUID_i; it2 ++) {
    pt_tmp2 = JetPUID_Pt.at(it2);
    eta_tmp2 = JetPUID_Eta.at(it2);
    phi_tmp2 = JetPUID_Phi.at(it2);
    m_tmp2 = JetPUID_m.at(it2);
    temp2.SetPtEtaPhiM(pt_tmp2,eta_tmp2,phi_tmp2,m_tmp2);
    sum = temp + temp2;
  /*  if(sum.Pt() > pt_sum)   {
        pt_sum = sum.Pt();
        jetID = it;
        jetID2 = it2;  
          }
        }
        }

*/


float tmp = 0.;       
 for(unsigned it = 0; it != JetPUID_i; it++)    {
    pt_tmp = JetPUID_Pt.at(it);
    eta_tmp = JetPUID_Eta.at(it);
    phi_tmp = JetPUID_Phi.at(it);
    m_tmp = JetPUID_m.at(it);
    if(pt_tmp > tmp) {
        temp.SetPtEtaPhiM(pt_tmp, eta_tmp, phi_tmp, m_tmp);
        tmp = pt_tmp;}
        }


float tmp2 = 0;
for(unsigned it2 =  1; it2 != JetPUID_i; it2 ++) {
    pt_tmp2 =  JetPUID_Pt.at(it2);
    eta_tmp2 = JetPUID_Eta.at(it2);
    phi_tmp2 = JetPUID_Phi.at(it2);
    m_tmp2 = JetPUID_m.at(it2);
    if(pt_tmp2 > tmp2)   {
        temp2.SetPtEtaPhiM(pt_tmp2, eta_tmp2, phi_tmp2, m_tmp2);
        tmp2 = pt_tmp2;}
        }
        
/*float tmp3 = 0;
TLorentzVector temp3;
for(unsigned it3 =  2; it3 != JetPUID_i; it3 ++) {
    pt_tmp2 =  JetPUID_Pt.at(it3);
    eta_tmp2 = JetPUID_Eta.at(it3);
    phi_tmp2 = JetPUID_Phi.at(it3);
    m_tmp2 = JetPUID_m.at(it3);
    if(pt_tmp2 > tmp3)   {
        temp3.SetPtEtaPhiM(pt_tmp2, eta_tmp2, phi_tmp2, m_tmp2);
        tmp3 = pt_tmp2;}
        }*/
       
 /*   if(sum.E() > pt_sum)   {
         pt_sum = sum.E();
         pt1 = pt_tmp;
         pt2 = pt_tmp2;
         eta1 = eta_tmp;
         eta2 = eta_tmp2;
         phi1 = phi_tmp;
         phi2 = phi_tmp2;
         m1 = m_tmp;
         m2 = m_tmp2;
         }
        }
        }
        //float Pt1 = PvJetPUID_Eta.at(jet_id), m=9.,t=5.,y=7.;
        //float Pt1 , m=9.,t=5.,y=7.;
        //cout<<PvJetPUID_Pt.size()<<endl;
        temp.SetPtEtaPhiM(pt1, eta1, phi1, m1);
        temp2.SetPtEtaPhiM(pt2, eta2, phi2, m2); */
sum = temp ;//+ temp2;
E_sum = sum.Pt();
/*vbf_etaphi.push_back(eta1);
vbf_etaphi.push_back(eta2); 
vbf_etaphi.push_back(phi1);
vbf_etaphi.push_back(phi2);       */       
vbf_etaphi.push_back(temp.Eta());
//vbf_etaphi.push_back(temp2.Eta()); 
vbf_etaphi.push_back(temp.Phi());
//vbf_etaphi.push_back(temp2.Phi());  
JetPUID_E.push_back(E_sum);
JetPUID_i = 0;
pt_sum = 0.;  


//--------- Jet with time cut
//PvJetPUID_Pt ,PvJetPUID_E,PvJetPUID_m,PvJetPUID_Eta,PvJetPUID_Phi;
   int PvJetPUID_entries = branch_PvJetPUID->GetEntriesFast(), PvJetPUID_i = 0;
 if(PvJetPUID_entries!=0)    {  
   while(PvJetPUID_i < PvJetPUID_entries)   {
        Jet* PvJetPUID = (Jet*) branch_PvJetPUID->At(PvJetPUID_i);
        PvJetPUID_Pt.push_back(PvJetPUID->PT);
        PvJetPUID_m.push_back(PvJetPUID->Mass);
        PvJetPUID_Eta.push_back(PvJetPUID->Eta);
        PvJetPUID_Phi.push_back(PvJetPUID->Phi);
        ++PvJetPUID_i;
        }
        
        //cout<<"check1:entries:"<<"\t"<<PvJetPUID_entries<<endl;
TLorentzVector pvtemp, pvtemp2; 
float pv_eta_tmp, pv_phi_tmp, pv_m_tmp, pv_E_sum=0., pv_pt_tmp, pv_pt_sum= 0., pv_eta, pv_phi, pv_pt, pv_m;    
float pv_pt1, pv_pt2, pv_eta1, pv_eta2, pv_phi1, pv_phi2, pv_m1, pv_m2, pv_pt_tmp2, pv_eta_tmp2, pv_phi_tmp2, pv_m_tmp2;
int jet_id=0, jet_id2=1;    
/*for(unsigned it = 0; it != PvJetPUID_i; it++)    {
    pv_pt_tmp = PvJetPUID_Pt.at(it);
    pv_eta_tmp = PvJetPUID_Eta.at(it);
    pv_phi_tmp = PvJetPUID_Phi.at(it);
    pv_m_tmp = PvJetPUID_m.at(it);
    pvtemp.SetPtEtaPhiM(pv_pt_tmp, pv_eta_tmp, pv_phi_tmp, pv_m_tmp);
for(unsigned it2 =  1; it2 != PvJetPUID_i; it2 ++) {
    pv_pt_tmp2 =  PvJetPUID_Pt.at(it2);
    pv_eta_tmp2 = PvJetPUID_Eta.at(it2);
    pv_phi_tmp2 = PvJetPUID_Phi.at(it2);
    pv_m_tmp2 = PvJetPUID_m.at(it2);
    pvtemp2.SetPtEtaPhiM(pv_pt_tmp2, pv_eta_tmp2, pv_phi_tmp2, pv_m_tmp2);
    sum = pvtemp + pvtemp2;
    if(sum.E() > pv_pt_sum)   {
        pv_pt_sum = sum.E();
         jet_id = it;
         jet_id2 = it2;
         pv_pt1 = pv_pt_tmp;
         pv_pt2 = pv_pt_tmp2;
         pv_eta1 = pv_eta_tmp;
         pv_eta2 = pv_eta_tmp2;
         pv_phi1 = pv_phi_tmp;
         pv_phi2 = pv_phi_tmp2;
         pv_m1 = pv_m_tmp;
         pv_m2 = pv_m_tmp2;
         }
        }
        } */
/*if(PvJetPUID_entries > 2)//{
         pv_pt1 = PvJetPUID_Pt.at(0);
         pv_pt2 = PvJetPUID_Pt.at(1);
         pv_eta1 = PvJetPUID_Eta.at(0);
         pv_eta2 = PvJetPUID_Eta.at(1);
         pv_phi1 = PvJetPUID_Phi.at(0);
         pv_phi2 = PvJetPUID_Phi.at(1);
         pv_m1 = PvJetPUID_m.at(0);
         pv_m2 = PvJetPUID_m.at(1);*/
 float pvtmp = 0.;       
 for(unsigned it = 0; it != PvJetPUID_i; it++)    {
    pv_pt_tmp = PvJetPUID_Pt.at(it);
    pv_eta_tmp = PvJetPUID_Eta.at(it);
    pv_phi_tmp = PvJetPUID_Phi.at(it);
    pv_m_tmp = PvJetPUID_m.at(it);
    if(pv_pt_tmp > pvtmp) {
        pvtemp.SetPtEtaPhiM(pv_pt_tmp, pv_eta_tmp, pv_phi_tmp, pv_m_tmp);
        pvtmp = pv_pt_tmp;}
        }
        //cout<<pvtemp.Pt()<<endl;

float pvtmp2 = 0;
for(unsigned it2 =  1; it2 != PvJetPUID_i; it2 ++) {
    pv_pt_tmp2 =  PvJetPUID_Pt.at(it2);
    pv_eta_tmp2 = PvJetPUID_Eta.at(it2);
    pv_phi_tmp2 = PvJetPUID_Phi.at(it2);
    pv_m_tmp2 = PvJetPUID_m.at(it2);
    if(pv_pt_tmp2 > pvtmp2)   {
        pvtemp2.SetPtEtaPhiM(pv_pt_tmp2, pv_eta_tmp2, pv_phi_tmp2, pv_m_tmp2);
        pvtmp2 = pv_pt_tmp2;}
        }
        
  /*  sum = pvtemp + pvtemp2;
    if(sum.E() > pv_pt_sum)   {
        pv_pt_sum = sum.E();
         jet_id = it;
         jet_id2 = it2;
         pv_pt1 = pv_pt_tmp;
         pv_pt2 = pv_pt_tmp2;
         pv_eta1 = pv_eta_tmp;
         pv_eta2 = pv_eta_tmp2;
         pv_phi1 = pv_phi_tmp;
         pv_phi2 = pv_phi_tmp2;
         pv_m1 = pv_m_tmp;
         pv_m2 = pv_m_tmp2;
         }
        }
        }   */
        
        
        
   //     pvtemp.SetPtEtaPhiM(pv_pt1, pv_eta1, pv_phi1, pv_m1);
     //   pvtemp2.SetPtEtaPhiM(pv_pt2, pv_eta2, pv_phi2, pv_m2);
pvtemp.SetPtEtaPhiM(PvJetPUID_Pt.at(0),PvJetPUID_Eta.at(0), PvJetPUID_Phi.at(0), PvJetPUID_m.at(0));
//pvtemp2.SetPtEtaPhiM(PvJetPUID_Pt.at(1),PvJetPUID_Eta.at(1), PvJetPUID_Phi.at(1), PvJetPUID_m.at(1));
sum =  pvtemp;// + pvtemp2;
//pv_E_sum = sum.E();
pv_E_sum = sum.Pt();
//cout<<pv_E_sum<<endl;
PvJetPUID_E.push_back(pv_E_sum);
/*pv_vbf_etaphi.push_back(pv_eta1);
pv_vbf_etaphi.push_back(pv_eta2);  
pv_vbf_etaphi.push_back(pv_phi1);
pv_vbf_etaphi.push_back(pv_phi2);*/
pv_vbf_etaphi.push_back(pvtemp.Eta());
//pv_vbf_etaphi.push_back(pvtemp2.Eta());  
pv_vbf_etaphi.push_back(pvtemp.Phi());
//pv_vbf_etaphi.push_back(pvtemp2.Phi());
PvJetPUID_i = 0;
pv_pt_sum = 0.;
//}
//cout<<"check2"<<endl;
//--------- Gen Jet
GenJet_Pt ,GenJet_E,GenJet_m,GenJet_Eta,GenJet_Phi;
   int GenJet_entries = branch_GenJet->GetEntriesFast(), GenJet_i = 0;
   
   while(GenJet_i < GenJet_entries)   {
        Jet* GenJet = (Jet*) branch_GenJet->At(GenJet_i);
        GenJet_Pt.push_back(GenJet->PT);        
        GenJet_m.push_back(GenJet->Mass);
        GenJet_Eta.push_back(GenJet->Eta);
        GenJet_Phi.push_back(GenJet->Phi);
        ++GenJet_i;
        }         
        
 /*       
float gen_eta_tmp, gen_phi_tmp, gen_m_tmp, gen_E_sum=0., gen_pt_tmp, gen_pt_sum= 0.;
float gen_pt_tmp2, gen_phi_tmp2, gen_eta_tmp2, gen_m_tmp2;
float gentmp = 0.;   
TLorentzVector gentemp, gentemp2;    
 for(unsigned it = 0; it != GenJet_i; it++)    {
    gen_pt_tmp = GenJet_Pt.at(it);
    gen_eta_tmp = GenJet_Eta.at(it);
    gen_phi_tmp = GenJet_Phi.at(it);
    gen_m_tmp = GenJetID_m.at(it);
    if(gen_pt_tmp > gentmp) {
        gentemp.SetPtEtaPhiM(gen_pt_tmp, gen_eta_tmp, gen_phi_tmp, gen_m_tmp);
        gentmp = gen_pt_tmp;}
        }
        //cout<<pvtemp.Pt()<<endl;

float gentmp2 = 0;
for(unsigned it2 =  1; it2 != GenJet_i; it2 ++) {
    gen_pt_tmp2 =  GenJet_Pt.at(it2);
    gen_eta_tmp2 = GenJet_Eta.at(it2);
    gen_phi_tmp2 = GenJet_Phi.at(it2);
    gen_m_tmp2 = GenJet_m.at(it2);
    if(gen_pt_tmp2 > gentmp2)   {
        gentemp2.SetPtEtaPhiM(gen_pt_tmp2, gen_eta_tmp2, gen_phi_tmp2, gen_m_tmp2);
        gentmp2 = gen_pt_tmp2;}
        }        
    */    
//----------LHEParticles
part_Pt ,part_E,part_m, part_Eta, part_Phi;
   int part_entries = branch_LHEParticles->GetEntriesFast(), part_i = 0;
   while(part_i < part_entries)   {
        GenParticle* part = (GenParticle*) branch_LHEParticles->At(part_i);
        //part_Pt.push_back(part->PT);
        part_Pt.push_back(part->Eta);                
        part_m.push_back(part->Mass);
        //part_Eta.push_back(part->Eta);
        //part_Phi.push_back(part->Phi);
        part_Eta.push_back(part->Phi);
        part_Phi.push_back(part->Rapidity);
        part_ID.push_back(part->PID);
        part_status.push_back(part->Status);
        ++part_i;
        }   
        
TLorentzVector genjet1, genjet2;
/*int j = part_ID.size(), part_j, flag = 0, IDjet=0;
float d_r = 1000., gen_pt, gen_phi, gen_eta, gen_m, d_r2=1000.;
for(part_j = 0; part_j < j; part_j++)    {
    if(fabs(part_ID.at(part_j)) < 7 && part_status.at(part_j) == 1)  {
    if(flag == 0)   {
        float eta_part = part_Eta.at(part_j);
        float phi_part = part_Phi.at(part_j);
        for(int k = 0; k < GenJet_Eta.size(); k++)   {
            float gen_eta_tmp = GenJet_Eta.at(k);
            float gen_phi_tmp = GenJet_Phi.at(k);
            float gen_pt_tmp = GenJet_Pt.at(k);
            float gen_m_tmp = GenJet_m.at(k);
            float d_r_tmp = DeltaR(eta_part, gen_eta_tmp, phi_part, gen_phi_tmp);
                if(d_r_tmp < d_r)   {
                     d_r = d_r_tmp;
                     gen_pt = gen_pt_tmp;
                     gen_eta = gen_eta_tmp;
                     gen_phi = gen_phi_tmp;
                     gen_m = gen_m_tmp;                     
                    }
                }
                flag = 1;
                genjet1.SetPtEtaPhiM(gen_pt, gen_eta, gen_phi, gen_m);  
                }
                else {
                float eta_part2 = part_Eta.at(part_j);
                float phi_part2 = part_Phi.at(part_j);
        for(int k = 0; k < GenJet_Eta.size(); k++)   {
            float gen_eta_tmp2 = GenJet_Eta.at(k);
            float gen_phi_tmp2 = GenJet_Phi.at(k);
            float gen_pt_tmp2 = GenJet_Pt.at(k);
            float gen_m_tmp2 = GenJet_m.at(k);
            float d_r_tmp2 = DeltaR(eta_part2, gen_eta_tmp2, phi_part2, gen_phi_tmp2);
                if(d_r_tmp2 < d_r2)   {
                    d_r2 = d_r_tmp2;
                     IDjet = k;
                     gen_pt = gen_pt_tmp2;
                     gen_eta = gen_eta_tmp2;
                     gen_phi = gen_phi_tmp2;
                     gen_m = gen_m_tmp2;                     
                    }
                    
                }
                genjet2.SetPtEtaPhiM(gen_pt, gen_eta, gen_phi, gen_m);
                }
                    }
                    }*/

genjet1.SetPtEtaPhiM(GenJet_Pt.at(0), GenJet_Eta.at(0), GenJet_Phi.at(0), GenJet_m.at(0));
//genjet2.SetPtEtaPhiM(GenJet_Pt.at(1), GenJet_Eta.at(1), GenJet_Phi.at(1), GenJet_m.at(1));
                    
TLorentzVector sum_jet = genjet1 ;//+ genjet2;
//TLorentzVector sum_jet = gentemp + gentemp2;
     //e_gen = sum_jet.E();

     e_gen = sum_jet.Pt();
     genejet_vect.push_back(e_gen);
     pv_e_frac = pv_E_sum/e_gen;
     e_frac = E_sum/e_gen;
     reco_frac = E_sum/pv_E_sum;
     float genjet1_eta = genjet1.Eta();
     //float genjet2_eta = genjet2.Eta();
     float genjet1_phi = genjet1.Phi();
     //float genjet2_phi = genjet2.Phi(); 
     gen_vbf_etaphi.push_back(genjet1_eta);
     //gen_vbf_etaphi.push_back(genjet2_eta);
     gen_vbf_etaphi.push_back(genjet1_phi);
     //gen_vbf_etaphi.push_back(genjet2_phi);
        
     float min_tmp=9000., dR_tmp=90000.;
     int jflag = 0;
     //straight
    // float dR_tmp1 = DeltaR(gen_vbf_etaphi.at(0),vbf_etaphi.at(0),gen_vbf_etaphi.at(2), vbf_etaphi.at(2));
     //float dR_tmp2 = DeltaR(gen_vbf_etaphi.at(1),vbf_etaphi.at(1),gen_vbf_etaphi.at(3), vbf_etaphi.at(3));
     //crossed
     //float dR_tmp3 = DeltaR(gen_vbf_etaphi.at(0),vbf_etaphi.at(1),gen_vbf_etaphi.at(2), vbf_etaphi.at(3));
     //float dR_tmp4 = DeltaR(gen_vbf_etaphi.at(1),vbf_etaphi.at(0),gen_vbf_etaphi.at(3), vbf_etaphi.at(2));  
   
      dR1 = DeltaR(gen_vbf_etaphi.at(0),vbf_etaphi.at(0),gen_vbf_etaphi.at(1), vbf_etaphi.at(1));
      //dR2 = DeltaR(gen_vbf_etaphi.at(1),vbf_etaphi.at(1),gen_vbf_etaphi.at(3), vbf_etaphi.at(3));
     //crossed
      //dR3 = DeltaR(gen_vbf_etaphi.at(0),vbf_etaphi.at(1),gen_vbf_etaphi.at(2), vbf_etaphi.at(3));
      //dR4 = DeltaR(gen_vbf_etaphi.at(1),vbf_etaphi.at(0),gen_vbf_etaphi.at(3), vbf_etaphi.at(2));
        dR.push_back(dR1);
        float tmp_efrac_jet1 = 9000., ejet=99999., tmp_efrac_jet2 = 9000.;
        if (dR1 < 100000.05)  {
             tmp_efrac_jet1 = (-genjet1.Pt() + temp.Pt())/genjet1.Pt();
             ejet = temp.Pt();}
         
   /*   float min_tmp1 = TMath::Min(dR1,dR2);
      float min_tmp2 = TMath::Min(dR3, min_tmp1);
      float min_tmp3 = TMath::Min(dR4, min_tmp2);
      dR_tmp = min_tmp3;
      float tmp_efrac_jet1 = 9000., ejet=99999., tmp_efrac_jet2 = 9000.;
      if(dR_tmp == dR1) {
        if(dR_tmp < 0.05 && dR2 < 0.05){
        jflag = 1,
        tmp_efrac_jet1 = (genjet1.Pt() - temp.Pt())/genjet1.Pt();   
        tmp_efrac_jet2 = (genjet2.Pt() -temp2.Pt())/genjet2.Pt();
        dR.push_back(dR1);
        dR.push_back(dR2);
        ejet = temp.Pt();}}
      else if(dR_tmp == dR2)   {
        if(dR_tmp < 0.05 && dR1 < 0.05){
        jflag = 2;
       tmp_efrac_jet1 = (genjet2.Pt() -temp2.Pt())/genjet2.Pt();
       tmp_efrac_jet2 = (genjet1.Pt() - temp.Pt())/genjet1.Pt();
       ejet = temp2.Pt();
       dR.push_back(dR2);
       dR.push_back(dR1);
       }}
      //else if ((dR_tmp == dR3) || (dR_tmp == dR4))  {
      else if (dR_tmp == dR3) {
      if(dR_tmp < 0.05 && dR4 < 0.05){
      jflag = 3;
        tmp_efrac_jet1 = (genjet1.Pt()-temp2.Pt())/genjet1.Pt();
        tmp_efrac_jet2 = (genjet2.Pt() - temp.Pt())/genjet2.Pt();
        ejet = temp2.Pt();
        dR.push_back(dR3);
        dR.push_back(dR4);
        }}
      else if (dR_tmp == dR4)   {
        if(dR_tmp < 0.05 && dR3 < 0.05)    {
        jflag = 4;
            tmp_efrac_jet1 = (genjet2.Pt() - temp.Pt())/genjet2.Pt();
            tmp_efrac_jet2 = (genjet1.Pt()-temp2.Pt())/genjet1.Pt();
            dR.push_back(dR4);
            dR.push_back(dR3);
            }}  
        */
      if(jflag==0 || jflag ==2 || jflag==3||jflag==4)    {
      efrac_jet.push_back(tmp_efrac_jet1);
      //efrac_jet.push_back(tmp_efrac_jet2);
      ejet_vect.push_back(ejet);
       }      
        
           pv_dR1 = DeltaR(gen_vbf_etaphi.at(0),pv_vbf_etaphi.at(0),gen_vbf_etaphi.at(0), pv_vbf_etaphi.at(0));
           //pv_dR1 = 100.;
          // pv_dR2 = DeltaR(gen_vbf_etaphi.at(1),pv_vbf_etaphi.at(1),gen_vbf_etaphi.at(3), pv_vbf_etaphi.at(3));
     //crossed
           //pv_dR3 = DeltaR(gen_vbf_etaphi.at(0),pv_vbf_etaphi.at(1),gen_vbf_etaphi.at(2), pv_vbf_etaphi.at(3));
           //pv_dR4 = DeltaR(gen_vbf_etaphi.at(1),pv_vbf_etaphi.at(0),gen_vbf_etaphi.at(3), pv_vbf_etaphi.at(2));
           cout<<"check1"<<endl;
           pv_dR.push_back(pv_dR1);
           
        float pv_tmp_efrac_jet = 9000., pv_tmp_efrac_jet2 = 9000.;   
        float gen_eta = 99., vbs_eta = 99.;  
     
     float pv_ejet=99999., gen_ejet = 90000.;
     int jetflag = 0;
     
     if ( pv_dR1 < 1111110.05){
        pv_tmp_efrac_jet = (-genjet1.Pt()+pvtemp.Pt())/genjet1.Pt();
        //pv_tmp_efrac_jet = 9.;
        pv_ejet = pvtemp.Pt();}
        //pv_ejet = 3.;}
     
   /*  float pv_min_tmp1 = TMath::Min(pv_dR1,pv_dR2);
     float pv_min_tmp2 = TMath::Min(pv_dR3, pv_min_tmp1);
     float pv_min_tmp3 = TMath::Min(pv_dR4, pv_min_tmp2);
     float pv_dR_tmp = pv_min_tmp3;
     float gen_eta = 99., vbs_eta = 99.;  
     float pv_tmp_efrac_jet = 9000., pv_tmp_efrac_jet2 = 9000.;
     float pv_ejet=99999., gen_ejet = 90000.;
     int jetflag = 0;
     if(pv_dR_tmp == pv_dR1)    {
        if(pv_dR_tmp < 0.05 && pv_dR2 < 0.05){
        jetflag = 1;
        //pv_tmp_efrac_jet = pvtemp.Pt()/genjet1.Pt();
        //pv_tmp_efrac_jet2 = pvtemp2.Pt()/genjet2.Pt();
        pv_tmp_efrac_jet = (-genjet1.Pt()+pvtemp.Pt())/genjet1.Pt();
        pv_tmp_efrac_jet2 = (-genjet2.Pt()+pvtemp2.Pt())/genjet2.Pt();
        pv_dR.push_back(pv_dR1);
        pv_dR.push_back(pv_dR2);
        pv_ejet = pvtemp.Pt();
        pv_ejet = pvtemp2.Pt();
        gen_ejet = genjet1.Pt();
        gen_eta = gen_vbf_etaphi.at(0);
        vbs_eta = pv_vbf_etaphi.at(0);
        }
        }
     else if(pv_dR_tmp == pv_dR2)   {
        if(pv_dR_tmp < 0.05 && pv_dR1 < 0.05){
        jetflag = 2;
        //pv_tmp_efrac_jet = pvtemp.Pt()/genjet1.Pt();
        //pv_tmp_efrac_jet2 = pvtemp2.Pt()/genjet2.Pt();
         pv_tmp_efrac_jet = (-genjet2.Pt()+pvtemp2.Pt())/genjet2.Pt();
         pv_tmp_efrac_jet2 = (-genjet1.Pt()+pvtemp.Pt())/genjet1.Pt();
         pv_dR.push_back(pv_dR2);
         pv_dR.push_back(pv_dR1);
         pv_ejet = pvtemp2.Pt();
         pv_ejet = pvtemp.Pt();
         gen_ejet = genjet2.Pt();
         gen_eta = gen_vbf_etaphi.at(1);
         vbs_eta = pv_vbf_etaphi.at(1);
         }}
     else if ((pv_dR_tmp == pv_dR3)){// || (pv_dR_tmp == pv_dR4))   {
        if(pv_dR_tmp < 0.05 && pv_dR4 < 0.05){
        jetflag = 3;
        //pv_tmp_efrac_jet = pvtemp2.Pt()/genjet1.Pt();
        //pv_tmp_efrac_jet2 = pvtemp.Pt()/genjet2.Pt();
        pv_tmp_efrac_jet = (-genjet1.Pt()+pvtemp2.Pt())/genjet1.Pt();
        pv_tmp_efrac_jet2 = (-genjet2.Pt()+pvtemp.Pt())/genjet2.Pt();
        pv_dR.push_back(pv_dR3);
        pv_dR.push_back(pv_dR4);
        pv_ejet = pvtemp2.Pt();
        pv_ejet = pvtemp.Pt();
        gen_ejet = genjet1.Pt();
        gen_eta = gen_vbf_etaphi.at(0);
        vbs_eta = pv_vbf_etaphi.at(1);
        }}
         else if (pv_dR_tmp == pv_dR4)  {
        if(pv_dR_tmp < 0.05 && pv_dR3 < 0.05){
        jetflag = 4;
        //pv_tmp_efrac_jet = pvtemp.Pt()/genjet2.Pt();
        //pv_tmp_efrac_jet2 = pvtemp2.Pt()/genjet1.Pt();
        pv_tmp_efrac_jet = (-genjet2.Pt()+pvtemp.Pt())/genjet2.Pt();
        pv_tmp_efrac_jet2 = (-genjet1.Pt()+pvtemp2.Pt())/genjet1.Pt();
        pv_dR.push_back(pv_dR3);
        pv_dR.push_back(pv_dR4);
        pv_ejet = pvtemp.Pt();
        pv_ejet = pvtemp2.Pt();
        gen_ejet = genjet2.Pt();
        gen_eta = gen_vbf_etaphi.at(1);   
        vbs_eta = pv_vbf_etaphi.at(0);
        }}*/
        //cout<<jetflag<<endl;
      if(jetflag== 0  || jetflag ==2 || jetflag==3||jetflag==4)    {
      if(pv_ejet > 0){
      pvefrac_jet.push_back(pv_tmp_efrac_jet);
      //pvefrac_jet.push_back(pv_tmp_efrac_jet2); 
      pvejet_vect.push_back(pv_ejet);}}
      //genejet_vect.push_back(gen_ejet);  }}
 /*    float reco_dR_tmp1 = DeltaR(vbf_etaphi.at(0),pv_vbf_etaphi.at(0),vbf_etaphi.at(2), pv_vbf_etaphi.at(2));
     float reco_dR_tmp2 = DeltaR(vbf_etaphi.at(1),pv_vbf_etaphi.at(1),vbf_etaphi.at(3), pv_vbf_etaphi.at(3));
     float reco_dR_tmp3 = DeltaR(vbf_etaphi.at(0),pv_vbf_etaphi.at(1),vbf_etaphi.at(2), pv_vbf_etaphi.at(3));
     float reco_dR_tmp;
             
     if(reco_dR_tmp1 < reco_dR_tmp2)  {
         //min_tmp == dR_tmp1;
         if(reco_dR_tmp1 < reco_dR_tmp3)  reco_dR_tmp =reco_dR_tmp1;    
         else   reco_dR_tmp = reco_dR_tmp3;
         }
    else    {
        if(reco_dR_tmp2 < reco_dR_tmp3)   reco_dR_tmp = reco_dR_tmp2;    
        else reco_dR_tmp = reco_dR_tmp3;
        }
     
    */
      //dR.push_back(dR_tmp);
      //dR.push_back(vbs_eta);
      //pv_dR.push_back(pv_dR_tmp);
      //reco_dR.push_back(reco_dR_tmp);
      reco_dR.push_back(gen_eta);
             
//--------- fill the output tree
	easyTree -> Fill();
//--------- Cleaning all tmp variables for the next entry
	JetPUID_Pt.erase(JetPUID_Pt.begin(),JetPUID_Pt.end());
	JetPUID_m.erase(JetPUID_m.begin(),JetPUID_m.end());
	JetPUID_Eta.erase(JetPUID_Eta.begin(),JetPUID_Eta.end());
	JetPUID_Phi.erase(JetPUID_Phi.begin(),JetPUID_Phi.end());
	JetPUID_E.erase(JetPUID_E.begin(),JetPUID_E.end());
	PvJetPUID_Pt.erase(PvJetPUID_Pt.begin(),PvJetPUID_Pt.end());
	PvJetPUID_m.erase(PvJetPUID_m.begin(),PvJetPUID_m.end());
	PvJetPUID_Eta.erase(PvJetPUID_Eta.begin(),PvJetPUID_Eta.end());
	PvJetPUID_Phi.erase(PvJetPUID_Phi.begin(),PvJetPUID_Phi.end());
	PvJetPUID_E.erase(PvJetPUID_E.begin(),PvJetPUID_E.end());
	GenJet_Pt.erase(GenJet_Pt.begin(),GenJet_Pt.end());
	GenJet_m.erase(GenJet_m.begin(),GenJet_m.end());
	GenJet_Eta.erase(GenJet_Eta.begin(),GenJet_Eta.end());
	GenJet_Phi.erase(GenJet_Phi.begin(),GenJet_Phi.end());
	GenJet_E.erase(GenJet_E.begin(),GenJet_E.end());
	part_Pt.erase(part_Pt.begin(),part_Pt.end());
	part_Eta.erase(part_Eta.begin(), part_Eta.end());
	part_Phi.erase(part_Phi.begin(), part_Phi.end());
	part_m.erase(part_m.begin(), part_m.end());
	part_ID.erase(part_ID.begin(),part_ID.end());
	part_status.erase(part_status.begin(), part_status.end());
	//vbf_eta.erase(vbf_eta.begin(), vbf_eta.end());
	//pv_vbf_eta.erase(pv_vbf_eta.begin(), pv_vbf_eta.end());
	//gen_vbf_eta.erase(gen_vbf_eta.begin(), gen_vbf_eta.end());
	vbf_etaphi.erase(vbf_etaphi.begin(), vbf_etaphi.end());
	pv_vbf_etaphi.erase(pv_vbf_etaphi.begin(), pv_vbf_etaphi.end());
	gen_vbf_etaphi.erase(gen_vbf_etaphi.begin(), gen_vbf_etaphi.end());
	dR.erase(dR.begin(),dR.end());
	pv_dR.erase(pv_dR.begin(), pv_dR.end());
	reco_dR.erase(reco_dR.begin(), reco_dR.end());
	dr_tmpvector.erase(dr_tmpvector.begin(),dr_tmpvector.end());
	efrac_jet.erase(efrac_jet.begin(),efrac_jet.end());
	pvefrac_jet.erase(pvefrac_jet.begin(),pvefrac_jet.end());
	ejet_vect.erase(ejet_vect.begin(), ejet_vect.end());
	pvejet_vect.erase(pvejet_vect.begin(), pvejet_vect.end());
	genejet_vect.erase(genejet_vect.begin(), genejet_vect.end());
	pv_e_frac = 0.;
	e_frac = 0.;
	pv_E_sum = 0.;
	E_sum = 0.;
	e_gen = 0.;
	//flag = 0;
	//dR1 = 0.;
	//dR2 = 0.;
	//dR3 = 0.;
    //dR4 = 0.;
	//d_r = 0.;
	//d_r2 = 0.;
	}
	cout<<iEvent<<endl;
    }
    easyTree -> Print("easyDelphes");
    outputFile -> Write();
    delete outputFile;
}
