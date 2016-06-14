/** \class Calorimeter
 *  Fills calorimeter towers, performs calorimeter resolution smearing,
 *  preselects towers hit by photons and creates energy flow objects.
 *  \author P. Demin - UCL, Louvain-la-Neuve
*/

#include "modules/ECalorimeter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TF2.h"
#include "TProfile.h"
#include "TH2F.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

static const double c_light     = 0.299792458 ; // mm/ps , or m/ns
static const double inv_c_light = 3.335640952 ; // ps/mm or ns/m

using namespace std;

//------------------------------------------------------------------------------

inline bool isMatching (TLorentzVector & part, vector<TLorentzVector> & set) 
{
  for (unsigned int i = 0 ; i < set.size () ; ++i)
    {
      if (part.DeltaR (set.at (i)) < 0.1) return true ;
    }
  return false ;
}

//------------------------------------------------------------------------------

ECalorimeter::ECalorimeter() :
  fECalResolutionFormula(0), fHCalResolutionFormula(0),
  fItParticleInputArray(0), fItTrackInputArray(0),
  fTowerTrackArray(0), fItTowerTrackArray(0),
  fItLHEPartonInputArray (0), fItVertexInputArray(0)
  {
  
  //fDelayBarrel and fDelayEndcap by Pietro
  //fDelayBarrel = new TF1 ("fDelayBarrel", "pol5", 0, 5) ;
  //fDelayEndcap = new TF1 ("fDelayEndcap", "pol5", 0, 5) ;
  
    
  //fDelayBarrel = new TF2 ("fDelayBarrel", "x/y", 0, 100, 0, 100);    
  //fDelayEndcap = new TF2 ("fDelayEndcap", "x/y", 0, 100, 0, 100);
  
  //corr_timeTower = new TF1 ( "corr_timeTower", "5 + TMath::Abs(x)", -1.479,1.479);
  //corr_timeTower = new TF1 ( "corr_timeTower", "35 * x^2", -4,4);
  //corr_timeTower2 = new TF1 ( "corr_timeTower2", "35 * x^2 - 5 * x", -4,4);
  
  fECalResolutionFormula = new DelphesFormula;
  fHCalResolutionFormula = new DelphesFormula;

  fTowerTrackArray = new TObjArray;
  fItTowerTrackArray = fTowerTrackArray->MakeIterator();
  
  fVertexArray = new TObjArray;
  fItVertexArray = fVertexArray->MakeIterator();
  
  fTimeVertexArray = new TObjArray;
  fItTimeVertexArray = fTimeVertexArray->MakeIterator();
  
  fCalTrackArray = new TObjArray;
  fItCalTrackArray = fCalTrackArray->MakeIterator();
  
  fRecoVertexArray = new TObjArray;
  fItRecoVertexArray = fRecoVertexArray->MakeIterator();
  
  fEFlowTowerOutputArray = new TObjArray;
  fItEFlowTowerOutputArray = fEFlowTowerOutputArray->MakeIterator();
  
  fPvEFlowTowerOutputArray = new TObjArray;
  fItPvEFlowTowerOutputArray = fPvEFlowTowerOutputArray->MakeIterator();
  
  fPvEFlowTrackOutputArray = new TObjArray;
  fItPvEFlowTrackOutputArray = fPvEFlowTrackOutputArray->MakeIterator();
  
  fPvEFlowTowerCandidateArray = new TObjArray;
  fItPvEFlowTowerCandidateArray = fPvEFlowTowerCandidateArray->MakeIterator();
  
  fnotimeTowerOutputArray = new TObjArray;
  fItnotimeTowerOutputArray = fnotimeTowerOutputArray->MakeIterator();

  fnotimeTrackOutputArray = new TObjArray;
  fItnotimeTrackOutputArray = fnotimeTrackOutputArray->MakeIterator();

}

//------------------------------------------------------------------------------

ECalorimeter::~ECalorimeter(){

  delete fDelayBarrel ;
  delete fDelayEndcap ;
  delete corr_timeTower;
  if(fECalResolutionFormula)        delete fECalResolutionFormula;
  if(fHCalResolutionFormula)        delete fHCalResolutionFormula;
//  if(fTowerTrackArray)              delete fTowerTrackArray;
  if(fItTowerTrackArray)            delete fItTowerTrackArray;
  if(fItLHEPartonInputArray)        delete fItLHEPartonInputArray;
  if(fItVertexInputArray)           delete fItVertexInputArray;
  if(fItTimeVertexArray)            delete fItTimeVertexArray;
  if(fItVertexingTrackInputArray)   delete fItVertexingTrackInputArray;
  //if(fCalTrackArray)                delete fCalTrackArray;
  //if(fItEFlowTowerOutputArray)      delete fItEFlowTowerOutputArray;
  //if(fItEFlowTrackOutputArray)      delete fItEFlowTrackOutputArray;
  if(fRecoVertexArray)              delete fRecoVertexArray;
  //if(fItPvEFlowTowerOutputArray)    delete fItPvEFlowTowerOutputArray;
  //if(fItPvEFlowTrackOutputArray)    delete fItPvEFlowTrackOutputArray;
  //if(fItPvEFlowTowerCandidateArray) delete fItPvEFlowTowerCandidateArray;
  //if(fItPvEFlowTrackOutputArray)    delete fItPvEFlowTrackOutputArray;
  //time_vs_eta->Reset();
}

//------------------------------------------------------------------------------

void ECalorimeter::Init(){

  ExRootConfParam param, paramEtaBins, paramPhiBins, paramFractions;
  Long_t i, j, k, size, sizeEtaBins, sizePhiBins;
  Double_t ecalFraction, hcalFraction;
  TBinMap::iterator itEtaBin;
  set< Double_t >::iterator itPhiBin;
  vector< Double_t > *phiBins;

  // read eta and phi bins
  param = GetParam("EtaPhiBins");
  size = param.GetSize();

  fBinMap.clear();
  fEtaBins.clear();
  fPhiBins.clear();

  for(i = 0; i < size/2; ++i){

    paramEtaBins = param[i*2];
    sizeEtaBins = paramEtaBins.GetSize();
    paramPhiBins = param[i*2 + 1];
    sizePhiBins = paramPhiBins.GetSize();

    for(j = 0; j < sizeEtaBins; ++j){
      for(k = 0; k < sizePhiBins; ++k){
        fBinMap[paramEtaBins[j].GetDouble()].insert(paramPhiBins[k].GetDouble());
      }
    }
  }

  // for better performance we transform map of sets to parallel vectors:
  // vector< double > and vector< vector< double >* >
  for(itEtaBin = fBinMap.begin(); itEtaBin != fBinMap.end(); ++itEtaBin){
    fEtaBins.push_back(itEtaBin->first);
    phiBins = new vector< double >(itEtaBin->second.size());
    fPhiBins.push_back(phiBins);
    phiBins->clear();
    for(itPhiBin = itEtaBin->second.begin(); itPhiBin != itEtaBin->second.end(); ++itPhiBin){
      phiBins->push_back(*itPhiBin);
    }
 }

  // read energy fractions for different particles
  param = GetParam("EnergyFraction");
  size = param.GetSize();

  // set default energy fractions values
  fFractionMap.clear();
  //fFractionMap[0] = make_pair(0.0, 1.0);
  fFractionMap[0] = 1.0;

  //consider only a fraction between ecal and hcal

  for(i = 0; i < size/2; ++i){
    paramFractions = param[i*2 + 1];
    ecalFraction = paramFractions[0].GetDouble();
    fFractionMap[param[i*2].GetInt()] = ecalFraction;
    //hcalFraction = paramFractions[1].GetDouble();
    //fFractionMap[param[i*2].GetInt()] = make_pair(ecalFraction, hcalFraction);
  }

  // read resolution formulas
  fECalResolutionFormula->Compile(GetString("ECalResolutionFormula", "0"));
  //fHCalResolutionFormula->Compile(GetString("HCalResolutionFormula", "0"));

  // import array with output from other modules
  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "ParticlePropagator/particles"));
  fItParticleInputArray = fParticleInputArray->MakeIterator();

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "ParticlePropagator/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  //PG for matching to LHE particles, for timing studies
  fLHEPartonInputArray = ImportArray(GetString("LHEPartonInputArray", "Delphes/LHEParticles"));
  fItLHEPartonInputArray = fLHEPartonInputArray->MakeIterator();
    
  // import smeared vertex array 
  fVertexInputArray = ImportArray(GetString("VertexInputArray", "ModifyBeamSpot/vertices"));
  fItVertexInputArray = fVertexInputArray->MakeIterator();
  
  //import tracks from vertexing
  fVertexingTrackInputArray = ImportArray(GetString("VertexingTrackInputArray", "Vertexing/vertexingtracks"));
  fItVertexingTrackInputArray = fVertexingTrackInputArray->MakeIterator();
  
  fPVInputArray = ImportArray (GetString("PVInputArray","ModifyBeamspot/PV"));
  fPVItInputArray = fPVInputArray->MakeIterator();

  // create output arrays
  // calo tower output
  fTowerOutputArray  = ExportArray(GetString("TowerOutputArray", "ecalTowers"));
  fItTowerOutputArray = fTowerOutputArray->MakeIterator();
  // photons 
  fPhotonOutputArray = ExportArray(GetString("PhotonOutputArray", "photons"));
  // track for particle flow
  fEFlowTrackOutputArray = ExportArray(GetString("EFlowTrackOutputArray", "ecalEflowTracks"));
  fPvEFlowTrackOutputArray = ExportArray(GetString("PvEFlowTrackOutputArray", "ecalPvEflowTracks"));
  fItPvEFlowTrackOutputArray = fEFlowTrackOutputArray->MakeIterator();
  // tower for particle flow
  fEFlowTowerOutputArray = ExportArray(GetString("EFlowTowerOutputArray", "ecalEflowTowers"));
  //fEFlowTowerOutputArray = ExportArray(GetString("PvEFlowTowerCandidateArray", "ecalEflowTowers"));
  fPvEFlowTowerOutputArray = ExportArray(GetString("PvEFlowTowerOutputArray", "ecalPvEflowTowers"));
  //tower with no assigned time
  fnotimeTowerOutputArray = ExportArray(GetString("notimeTowerOutputArray", "notimeEcalEflowTowers"));
  fnotimeTrackOutputArray = ExportArray(GetString("notimeTracksOutputArray", "notimeEcalEflowTracks"));

  // For timing
  // So far this flag needs to be false
  // since the curved extrapolation not supported
  // if this value is true, electron timing is taken from the track,
  // otherwise is taken from the particle collection
  fElectronsFromTrack  = false;
  
  //PG FIXME where does this come from?
  // suggested from A. Bornheim, reasonable according to him    
  
  fTimingEMin = GetDouble ("TimingEMin", 0.) ;
  fTimeCut = GetDouble ("TimeCut", 30.);

  //simple outputs during running
  fDebugOutputCollector.addVariable3D ("DR_tmp") ;
  fDebugOutputCollector.addVariable3D ("m_partType:eta:E") ;
  fDebugOutputCollector.addVariable3D ("Nm_partType:ID:IsPU") ;
  fDebugOutputCollector.addVariable3D ("PUindex:trackID:vertexID") ;
  fDebugOutputCollector.addVariable3D ("m_eta:time:pt") ;
  fDebugOutputCollector.addVariable3D ("m_eta:time_nc:pt") ;
  fDebugOutputCollector.addVariable4D ("Nm_eta:time:pt:size") ;
  fDebugOutputCollector.addVariable4D ("eta:time:pt:PID") ;

  fDebugOutputCollector.addVariable3D ("tower_eta:time:pt") ;
  fDebugOutputCollector.addVariable4D ("tower_eta:time:E:res") ;
  fEventCounter = 0 ;
  
  time_vs_eta = new TProfile ("time_vs_eta","timespread_vs_eta",200, -2.,2.,3800.,12000.);
//    time_vs_eta = new TH2F ( "time_vs_eta","time_vs_eta", 200, -3,3,200,3000,12000);
}

//------------------------------------------------------------------------------

void ECalorimeter::Finish(){

  std::string outfile = GetString ("simpleOutputFileName", "simpleOutput_ECal.root") ;
  fDebugOutputCollector.save (outfile) ;
//time_vs_eta->Draw();
  if(fItParticleInputArray)  delete fItParticleInputArray;
  if(fItTrackInputArray)     delete fItTrackInputArray;
  if(fItVertexInputArray)    delete fItVertexInputArray;
  if(fItVertexingTrackInputArray) delete fItVertexingTrackInputArray;
  if(fItRecoVertexArray)     delete fItRecoVertexArray;
//  if(fCalTrackArray)         delete fCalTrackArray;
//time_vs_eta->Reset();
  vector< vector< Double_t >* >::iterator itPhiBin;
  for (itPhiBin = fPhiBins.begin(); itPhiBin != fPhiBins.end(); ++itPhiBin){
    delete *itPhiBin;
  
  }
}

//------------------------------------------------------------------------------

void ECalorimeter::Process(){

  Candidate *particle, *track, *vertices, *vertex, *timevertex;
  TLorentzVector position, momentum;
  Short_t etaBin, phiBin, flags;
  Int_t number;
  Long64_t towerHit, towerEtaPhi, hitEtaPhi;
  Double_t ecalFraction, hcalFraction;
  Double_t ecalEnergy, hcalEnergy, totalEnergy;
  Int_t pdgCode;

  TFractionMap::iterator itFractionMap;

  vector< Double_t >::iterator itEtaBin;
  vector< Double_t >::iterator itPhiBin;
  vector< Double_t > *phiBins;

  vector<Long64_t>::iterator itTowerHits;

  DelphesFactory *factory = GetFactory();
  fTowerHits.clear();
  fTowerECalFractions.clear();
  fTowerHCalFractions.clear();
  fTrackECalFractions.clear();
  fTrackHCalFractions.clear();
  fParticlePDGId.clear();
  fTrackPDGId.clear();
  
  gen_energy = 0.;
  pu_gen_energy = 0.;
  
  // get the eta, phi, pt of the LHE partons that could generate
  // a Delphes jet
  vector<TLorentzVector> LHEParticles ;
  fItLHEPartonInputArray->Reset () ;
  while (Candidate * LHEparticle = static_cast<Candidate*> (fItLHEPartonInputArray->Next ()))
  
    { 
      if (LHEparticle->Status != 1) continue ;
      if (fabs (LHEparticle->PID) == 11 ||   // electron
//          1) 
          fabs (LHEparticle->PID) < 7   ||   // quarks
          fabs (LHEparticle->PID) == 21 ||   // gluon
          fabs (LHEparticle->PID) == 22)     // photon
        LHEParticles.push_back (LHEparticle->Momentum) ;
    }



fItVertexInputArray->Reset();
number=-1;
while((vertices = static_cast<Candidate*>(fItVertexInputArray->Next()))){

//vertex = static_cast<Candidate*>(fVertexInputArray->At(number));
    //while(fItVertexInputArrary->Next()){
//    fDebugOutputCollector.fillContainer ("PUindex", vertices->Position.Z()) ;
    number++;


}
  // loop over all tracks to get the deposited energy due to the
  // charged hadrons and electrons. 
  // THis energy is added to the tower as additional information,
  // that will be used in FinalizeTower to mimic the particle flow
  
  fItTrackInputArray->Reset();
  //fItVertexingTrackInputArray->Reset();
  number = -1;

  //while((track = static_cast<Candidate*>(fItVertexingTrackInputArray->Next()))){
  while((track = static_cast<Candidate*>(fItTrackInputArray->Next()))){

    const TLorentzVector &trackPosition = track->Position;
    ++number;  
   
    pdgCode = TMath::Abs(track->PID);

    itFractionMap = fFractionMap.find(pdgCode);
    if(itFractionMap == fFractionMap.end()){
      itFractionMap = fFractionMap.find(0);
    }

    ecalFraction = itFractionMap->second;
    fTrackECalFractions.push_back(ecalFraction);
    fTrackPDGId.push_back(pdgCode);
    /*ecalFraction = itFractionMap->second.first;
    hcalFraction = itFractionMap->second.second;

    fTrackECalFractions.push_back(ecalFraction);
    fTrackHCalFractions.push_back(hcalFraction);
    fTrackPDGId.push_back(pdgCode); */
 
    // find eta bin [1, fEtaBins.size - 1]
    itEtaBin = lower_bound(fEtaBins.begin(), fEtaBins.end(), trackPosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    etaBin = distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    phiBins = fPhiBins[etaBin];

    // find phi bin [1, phiBins.size - 1]
    itPhiBin = lower_bound(phiBins->begin(), phiBins->end(), trackPosition.Phi());
    if(itPhiBin == phiBins->begin() || itPhiBin == phiBins->end()) continue;
    phiBin = distance(phiBins->begin(), itPhiBin);

    flags = 1; // charged tracks flags equal one add to the tower hits
 
    // make tower hit:
    //    {16-bits for eta bin number, 
    //     16-bits for phi bin number, 
    //     8-bits for flags, 
    //     24-bits for track number}
    towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);

    fTowerHits.push_back(towerHit);
  } // loop over all tracks

  // loop over all particles of the event,
  // to get the energy of all the particles that hit the calorimeter
  fItParticleInputArray->Reset();
  number = -1;
  int ph_counter = 0;
  double eta1,eta2;
  while((particle = static_cast<Candidate*>(fItParticleInputArray->Next()))){
   
   
  if ( particle->PID == 22 && ph_counter == 0 && particle->IsPU == 0) {
   eta1 = particle->Position.Eta();
  ++ph_counter;}
  if ( particle->PID == 22 && ph_counter == 1 && particle->IsPU == 0)  {
   eta2 = particle->Position.Eta();
  ++ph_counter;
  }
  if (ph_counter == 2 && TMath::Abs(eta1 - eta2) < 0.8 ){
fDebugOutputCollector.fillContainer3D ("m_eta:time:pt", particle->Momentum.Pt(), eta1, ph_counter);
++ph_counter;
}


    const TLorentzVector &particlePosition = particle->Position;
    ++number;
    pdgCode = TMath::Abs(particle->PID);
    itFractionMap = fFractionMap.find(pdgCode); // find the particle in the fraction map

    if(itFractionMap == fFractionMap.end()){
      itFractionMap = fFractionMap.find(0);
    }
    
    //fill ecal tower fraction vectors
    ecalFraction = itFractionMap->second;
    fTowerECalFractions.push_back(ecalFraction);
    fParticlePDGId.push_back(pdgCode);

    //ecalFraction = itFractionMap->second.first; // take ecal fraction
    //hcalFraction = itFractionMap->second.second; // take Hcal fraction
  
    // fill tower fraction vectors
    /*fTowerECalFractions.push_back(ecalFraction);
    fTowerHCalFractions.push_back(hcalFraction);
    fParticlePDGId.push_back(pdgCode);  */

    //if(ecalFraction < 1.0E-9 && hcalFraction < 1.0E-9) continue;
    if(ecalFraction < 1.0E-9) continue;

    // find eta bin [1, fEtaBins.size - 1] --> find seed position
    itEtaBin = lower_bound(fEtaBins.begin(), fEtaBins.end(), particlePosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    etaBin = distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    phiBins = fPhiBins[etaBin];

    // find phi bin [1, phiBins.size - 1]
    itPhiBin = lower_bound(phiBins->begin(), phiBins->end(), particlePosition.Phi());
    if(itPhiBin == phiBins->begin() || itPhiBin == phiBins->end()) continue;
    phiBin = distance(phiBins->begin(), itPhiBin);

    flags  = 0;
    flags |= (pdgCode == 11 || pdgCode == 22) << 1; // flag for electrons and photons

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for particle number}
    towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);
    fTowerHits.push_back(towerHit); // fill tower hit vectors
  } // loop over all particles of the event

  // all hits are sorted first by eta bin number, then by phi bin number,
  // then by flags and then by particle or track number
  sort(fTowerHits.begin(),fTowerHits.end()); // sort the hit

  // loop over all hits (tracks and calo deposits)
  towerEtaPhi = 0;
  fTower = 0;
  
  // loop over the hits
  for(itTowerHits = fTowerHits.begin(); itTowerHits != fTowerHits.end(); ++itTowerHits) { 
    towerHit  = (*itTowerHits);
    flags     = (towerHit >> 24) & 0x00000000000000FFLL; // take the flag
    number    = (towerHit) & 0x0000000000FFFFFFLL;       // take the number
    hitEtaPhi = towerHit >> 32;
    if(towerEtaPhi != hitEtaPhi){  
      // first time hit, no  more than once since we have sorted tower 
      // as a function of eta and phi hits

      // switch to next tower
      towerEtaPhi = hitEtaPhi;
      // finalize previous tower
      FinalizeTower();

      // create new tower using the calorimeter information
      fTower = factory->NewCandidate();
      // store which type of particle it belongs to
      phiBin = (towerHit >> 32) & 0x000000000000FFFFLL;
      etaBin = (towerHit >> 48) & 0x000000000000FFFFLL;
      // phi bins for given eta bin
      phiBins = fPhiBins[etaBin];

      // calculate eta and phi of the tower's center
      fTowerEta = 0.5*(fEtaBins[etaBin - 1] + fEtaBins[etaBin]);
      fTowerPhi = 0.5*((*phiBins)[phiBin - 1] + (*phiBins)[phiBin]);
      

      // take the edges of the tower
      fTowerEdges[0] = fEtaBins[etaBin - 1];
      fTowerEdges[1] = fEtaBins[etaBin];
      fTowerEdges[2] = (*phiBins)[phiBin - 1];
      fTowerEdges[3] = (*phiBins)[phiBin];

      //smear eta and phi of tower's center
      //fTowerEta = gRandom->Uniform(fTowerEdges[0], fTowerEdges[1]);
      //fTowerPhi = gRandom->Uniform(fTowerEdges[2], fTowerEdges[3]);
      //fTowerEta = 0.5*(fEtaBins[etaBin - 1] + fEtaBins[etaBin]);
      //fTowerPhi = 0.5*((*phiBins)[phiBin - 1] + (*phiBins)[phiBin]);
      
      fTowerECalEnergy = 0.0;
      //fTowerHCalEnergy = 0.0;

      fTrackECalEnergy = 0.0;
      //fTrackHCalEnergy = 0.0;

      fTowerTrackHits = 0;
      fTowerPhotonHits = 0;
//fTower->IsPU = 10;
      fTowerTrackArray->Clear();
    }  // first time hit

    // check for track hits
    if(flags & 1){  
      ++fTowerTrackHits;
      //track = static_cast<Candidate*>(fTrackInputArray->At(number));
      track = static_cast<Candidate*>(fVertexingTrackInputArray->At(number));
      momentum = track->Momentum;
      
      Candidate *prt = static_cast<Candidate*>(track->GetCandidates()->Last());
        const TLorentzVector &ini = prt->Position;
      
      bool dbg_scz = false;
      if (dbg_scz) {
        cout << "   Calorimeter input track has x y z t " << track->Position.X() 
             << " " << track->Position.Y() << " " << track->Position.Z() << " " << track->Position.T() 
             << endl;
        //Candidate *prt = static_cast<Candidate*>(track->GetCandidates()->Last());
        //const TLorentzVector &ini = prt->Position;
        cout << "                and parent has x y z t " << ini.X() << " " << ini.Y() 
             << " " << ini.Z() << " " << ini.T();
      }

      ecalEnergy = momentum.E() * fTrackECalFractions[number];
      //hcalEnergy = momentum.E() * fTrackHCalFractions[number];

      fTrackECalEnergy += ecalEnergy;
      //fTrackHCalEnergy += hcalEnergy;

      //totalEnergy = ecalEnergy + hcalEnergy ;
      totalEnergy = ecalEnergy;
      
      //float time = track->Position.T () ;
      float TrackPositionZ, TrackPositionX, TrackPositionY;
      TrackPositionZ= track->Position.Z ();
      TrackPositionX= track->Position.X ();
      TrackPositionY= track->Position.Y ();
      float TrackMomentumX, TrackMomentumY, TrackMomentumZ;
      TrackMomentumX= track->Momentum.X ();
      TrackMomentumY= track->Momentum.Y ();
      TrackMomentumZ= track->Momentum.Z();     
           
      //double distance = TMath::Sqrt(TMath::Power(TrackPositionZ,2)+TMath::Power(TrackPositionX,2)+TMath::Power(TrackPositionY,2));
      //double velocity = (TMath::Sqrt(TMath::Power(TrackMomentumZ,2)+TMath::Power(TrackMomentumX,2)+TMath::Power(TrackMomentumY,2))/ totalEnergy);// * c_light ;
     
      if ( totalEnergy > fTimingEMin && fTower && fElectronsFromTrack) {
        
        float delay = 0. ;
        float TOF = 0. ;        
        float feta = fabs (track->Position.Eta ()) ;        
        double distance = TMath::Sqrt(TMath::Power(TrackPositionX,2)+TMath::Power(TrackPositionY,2));
        double velocity = (TMath::Sqrt(TMath::Power(track->Momentum.Z(),2)+TMath::Power(track->Momentum.X(),2)+TMath::Power(track->Momentum.Y(),2))/track->Momentum.E()) ;
                   
                   float Pt = track->Momentum.Pt();
                   float R = (Pt/(0.3*3.8));//Bz; m
                   //secant with respect to 0.0.0
                   float s = (TMath::Sqrt(TMath::Power(track->Position.X(),2)+TMath::Power(track->Position.Y(),2)));
                   float alpha = TMath::ASin(s/(1000*2*R));
                   float arc = 2 * alpha * R * 1000;
                   
                   float z_distance = track->Position.Z();
                   TOF = (TMath::Sqrt((TMath::Power(track->Position.Z(),2) + TMath::Power(arc,2)))) *inv_c_light;
        
        
        float time = track->Position.T () * inv_c_light ;
        //cout<<time<<"\t"<<TOF<<endl;
        fTower->ecal_E_t.push_back(
            //std::make_pair<float,float>(totalEnergy, time - TOF));
            std::make_pair<float,float>(totalEnergy, (double)time));
            fTower->ecal_E_eta.push_back(std::make_pair<float,float>(totalEnergy, track->Position.Eta()));
      
            
      /*if(particle->IsPU == 0 && pdgCode != 11) {
        fTower->IsPU = 0;
        gen_energy += totalEnergy;}
       if(particle->IsPU == 1 && fTower->IsPU != 0) fTower->IsPU = 1; */      
      }     
      fTowerTrackArray->Add(track);

      continue; // go to the next hit
    } // check for track hits

    // check for photon and electron hits in current tower
    if(flags & 2) ++fTowerPhotonHits;

    particle = static_cast<Candidate*>(fParticleInputArray->At(number));
    momentum = particle->Momentum;

    // fill current tower
    ecalEnergy = momentum.E() * fTowerECalFractions[number];
    //hcalEnergy = momentum.E() * fTowerHCalFractions[number];

    fTowerECalEnergy += ecalEnergy;
    //fTowerHCalEnergy += hcalEnergy;

    //totalEnergy = ecalEnergy + hcalEnergy ;
    totalEnergy = ecalEnergy;
    
    if ( (totalEnergy > fTimingEMin && fTower) &&
         (abs(particle->PID) != 11 || !fElectronsFromTrack) ) {
/*if( particle->PID == 22){
fDebugOutputCollector.fillContainer3D ("PUindex:trackID:vertexID", particle->Position.X(), particle->Position.Y(), TMath::Sqrt(TMath::Power(particle->Position.X(),2)+TMath::Power(particle->Position.Y(),2)));
 fDebugOutputCollector.fillContainer3D ("Nm_partType:ID:IsPU", particle->Position.T()*inv_c_light,particle->Position.Z(), TMath::Sqrt(TMath::Power(particle->Position.X(),2)+TMath::Power(particle->Position.Y(),2)) ) ;
 }*/
        // the time in the delphes propagation is in mm
        float time = particle->Position.T () *inv_c_light;
        //float time_theta;
        float delay = 0. ;
        double TOF = 0.;
        //float feta = fabs (particle->Position.Eta ()) ;
        float feta = particle->Position.Eta () ;
        double ParticlePositionZ, ParticlePositionX, ParticlePositionY;
        ParticlePositionZ = particle->Position.Z ();
        ParticlePositionX = particle->Position.X ();        
        ParticlePositionY = particle->Position.Y ();   
        double ParticleMomentumZ, ParticleMomentumX, ParticleMomentumY;
        ParticleMomentumZ = particle-> Momentum.Z();     
        ParticleMomentumX = particle-> Momentum.X();        
        ParticleMomentumY = particle-> Momentum.Y();    
                       
        double distance = TMath::Sqrt(TMath::Power(ParticlePositionZ,2)+TMath::Power(ParticlePositionX,2)+TMath::Power(ParticlePositionY,2));;
        //double velocity = (TMath::Sqrt(TMath::Power(ParticleMomentumZ,2)+TMath::Power(ParticleMomentumX,2)+TMath::Power(ParticleMomentumY,2))/ totalEnergy) ;//* c_light ;
         
         //double distance = TMath::Sqrt(TMath::Power(ParticlePositionX,2)+TMath::Power(ParticlePositionY,2));;
         //double velocity = (TMath::Sqrt(TMath::Power(particle->Momentum.Z(),2)+TMath::Power(particle->Momentum.X(),2)+TMath::Power(particle->Momentum.Y(),2))/track->Momentum.E()) ;
                   
                   float Pt = particle->Momentum.Pt();
                   float R = (Pt/(0.3*3.8));//Bz; m
                   //secant with respect to 0.0.0
                   float s = (TMath::Sqrt(TMath::Power(particle->Position.X(),2)+TMath::Power(particle->Position.Y(),2)));
                   float alpha = TMath::ASin(s/(1000*2*R));
                   float arc = 2 * alpha * R * 1000;
                   
                   float z_distance = track->Position.Z();
                   //TOF = (TMath::Sqrt((TMath::Power(particle->Position.Z(),2) + TMath::Power(arc,2)))) *inv_c_light;       
                   TOF = distance*inv_c_light;
        
      /*  if (feta < 1.6) TOF = fDelayBarrel->Eval(distance,velocity);
        
        else                  TOF = fDelayEndcap->Eval(distance,velocity);  */
            
                  float time_corr = particle->Position.T()*inv_c_light - TOF;
                  //particle->Position.SetT(time_corr);
                  fDebugOutputCollector.fillContainer3D ("PUindex:trackID:vertexID", TOF, particle->Position.Y(), time_corr);
        //time = time_corr;
        //cout<<time<<"\t"<<TOF<<endl;
        // the timing is corrected in the map,
        // since I assume that the knowledge of the position of the detID
        // is enough to determine the correction.
        // The PV z position could be used to introduce a correction
        // to the delay of the order of cz
        fTower->ecal_E_t.push_back(
          //std::make_pair<float,float> (totalEnergy, time - delay);
          std::make_pair<float,float> (totalEnergy, (double)time));   
          
           fTower->ecal_E_eta.push_back(std::make_pair<float,float>(totalEnergy, particle->Position.Eta())); 
        //  std::make_pair<float,float> (totalEnergy, time - TOF));
//         double Rmin = 100. ; 
//         for (unsigned int i = 0 ; i < LHEParticles.size () ; ++i)
//           {
//             if (particle->Position.DeltaR (LHEParticles.at (i)) < Rmin) 
//               Rmin = particle->Position.DeltaR (LHEParticles.at (i)) ;
//           }
//         fDebugOutputCollector.fillContainer ("DR", Rmin) ;
       fTower->IsPU = 1;
       
       if(particle->IsPU == 0) {
       fTower->IsPU = 0;
       gen_energy += totalEnergy;}
       if(particle->IsPU == 1 && fTower->IsPU != 0) {
       fTower->IsPU = 1;
       pu_gen_energy += totalEnergy;
       }
       //fTower->Charge = 0;
       //cout<<"charge"<<"\t"<<particle->Charge<<endl;
       if(particle->Charge != 0 )    fTower->Charge == particle->Charge;
       //if(particle->PID != 22)  fTower->PID == 1000;
       fTower->PID = particle->PID;
        
        

       // fDebugOutputCollector.fillContainer ("PUindex", particle->IsPU) ;



//fDebugOutputCollector.fillContainer3D ("m_eta:time:pt", particle->PID, particle->Eta, time);
      
//}
        if (isMatching (particle->Position, LHEParticles))  
         {               
               
                        
            //fDebugOutputCollector.fillContainer ("m_partType", particle->PID) ;
            //fDebugOutputCollector.fillContainer3D ("m_eta:time:pt", 
            //feta, (time - TOF) , particle->Momentum.Pt ()) ; //feta
         //   fDebugOutputCollector.fillContainer3D ("m_eta:time_nc:pt", 
           // feta, time * inv_c_light, particle->Momentum.Pt ()) ;    //feta
            
          }
            
            
          //}
        else
          {
           // fDebugOutputCollector.fillContainer ("Nm_partType", particle->PID) ;
            //fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt", 
              // feta, (time - TOF) , particle->Momentum.Pt ()) ;
          }
    }

    fTower->PID = fParticlePDGId.at(number);
    fTower->AddCandidate(particle);
  } // loop over the hits


//tracks from vertexing

//set vertices time
/*
int n = 1;
float vertex_time = 0;
int sum_track_time = 0;
int vertexID;
int trackID ;
float TOF;


fItVertexInputArray->Reset();
//fVertexArray->Clear();
    while((vertex = static_cast<Candidate*>(fItVertexInputArray->Next())))  {
       vertexID = vertex->VertexID_gen;
       //cout<<"vertex:"<<vertexID<<endl;
       fItVertexingTrackInputArray->Reset();
        while((track = static_cast<Candidate*>(fItVertexingTrackInputArray->Next())))  {
        pdgCode = TMath::Abs(track->PID);
        itFractionMap = fFractionMap.find(pdgCode);
        ecalFraction = itFractionMap->second.first;
       // if(track->Momentum.E()*ecalFraction>2.)    {
        
            if(track->Momentum.Pt()>1. && fabs(track->Position.Eta())<4 && track->Momentum.P()>1.) {
                  
                  trackID = track->VertexID_gen; 
                
                //cout<<"track:"<<trackID<<endl; 
                    if(trackID == vertexID)   {
                    //++n;
                    double distance = TMath::Sqrt(TMath::Power(track->Position.Z(),2)+TMath::Power(track->Position.X(),2)+TMath::Power(track->Position.Y(),2));
                    double velocity = (TMath::Sqrt(TMath::Power(track->Momentum.Z(),2)+TMath::Power(track->Momentum.X(),2)+TMath::Power(track->Momentum.Y(),2))/track->Momentum.E()) ;
                   
                   float Pt = track->Momentum.Pt();
                   float R = (Pt/(0.3*3.8));//Bz; m
                   //secant with respect to 0.0.0
                   float s = (TMath::Sqrt(TMath::Power(track->Position.X(),2)+TMath::Power(track->Position.Y(),2)));
                   float alpha = TMath::ASin(s/(1000*2*R));
                   float arc = 2 * alpha * R * 1000;
                   
                   float z_distance = track->Position.Z();
                   TOF = (TMath::Sqrt((TMath::Power(track->Position.Z(),2) + TMath::Power(arc,2)))) *inv_c_light;
                   
                   //cout<<Pt<<"\t"<<track->Momentum.Pz()<<"\t"<<track->Momentum.P()<<"\t"<<track->Position.Eta()<<"\t"<<R<<"\t"<<alpha<<"\t"<<TOF - track->Position.T()*inv_c_light<<endl;
                                    
                   //cout<<track->Position.T()*inv_c_light<<"\t"<<TOF<<"\t"<<track->Position.T()*inv_c_light - TOF<<"\t"<<vertexID<<endl;                   
                   float track_time = (track->Position.T() * inv_c_light - TOF);
                   if(track_time > 160.) { 
                    track_time = 0.;
                    n = n - 1;
                    }
                   // cout<<track->Position.T()*inv_c_light<<"\t"<<TOF<<"\t"<<track_time<<"\t"<<vertexID<<endl;                   
                  
                   sum_track_time += track_time;
                   ++n;
                   
                   
                   fDebugOutputCollector.fillContainer3D ("m_eta:time_nc:pt", 
                   track->Position.T()*inv_c_light , track_time, TOF ) ;  
                   
                   /*fDebugOutputCollector.fillContainer3D ("tower_eta:time:pt",
                   track_time*1.E-7,
                   n,
                   Pt ); */
                   /*
                  // n++;
                    } 
                                
         }
           //}
         }   
            vertex_time = sum_track_time/n;
            //cout<<vertex_time<<endl;
            //cout<<sum_track_time<<"\t"<<n<<endl;
            //fDebugOutputCollector.fillContainer3D ("m_eta:time_nc:pt", 
             //vertexID , vertex_time, trackID ) ;  
             //vertex->Position.SetT(vertex_time);
             //fDebugOutputCollector.fillContainer3D ("m_eta:time_nc:pt", 
             //vertex->Position.T() , vertex_time, trackID ) ;  

             //fVertexOutputArray->Add(vertex);
          /*   factory = GetFactory();
             timevertex = factory->NewCandidate();
             timevertex = static_cast<Candidate*>(vertex->Clone());
             timevertex->Position.SetT(vertex_time);
             cout<<"checka"<<endl;
             fVertexOutputArray->Add(timevertex);
             fDebugOutputCollector.fillContainer3D ("m_eta:time_nc:pt", 
             timevertex->Position.T() , vertex_time, trackID ) ;  
             cout<<"check2"<<endl; 
           */  
           /*  
            timevertex = vertex ;
            timevertex = static_cast<Candidate*> (vertex->Clone()) ;
            //vertexcandidate = static_cast<Candidate*> (vertexcandidate->Next ()) ;
            timevertex->Position.SetT(vertex_time) ;
            //vertexcandidate->AddCandidate (vertexmother) ;
   // cout<<"aa"<<endl;
            fTimeVertexArray->Add (timevertex) ;
            /*fDebugOutputCollector.fillContainer3D ("tower_eta:time:pt",
            timevertex->Position.T(),
            n,
            timevertex->VertexID_gen
    ) ; 
        
         n = 1;
         sum_track_time = 0;
        }
        */
       
         
  


  // finalize last tower
  FinalizeTower();
  PileUpSubtractor();
  ++fEventCounter ;

}


//------------------------------------------------------------------------------

void ECalorimeter::FinalizeTower(){

  Candidate *track, *tower;
  Double_t energy, pt, eta, phi;
  Double_t ecalEnergy, hcalEnergy;
  Double_t ecalSigma, hcalSigma;


  if(!fTower) return;

  // ECAL resolution
  ecalSigma  = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerECalEnergy);
  ecalEnergy = LogNormal(fTowerECalEnergy, ecalSigma);
    
  // HCAL resolution    
  /*
  hcalSigma  = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerHCalEnergy);
  hcalEnergy = LogNormal(fTowerHCalEnergy, hcalSigma);  */
    
  //energy     = ecalEnergy + hcalEnergy;
  energy = ecalEnergy;

  //smeared point of arrival of particles at tower
  //eta = gRandom->Uniform(fTowerEdges[0], fTowerEdges[1]);
  //phi = gRandom->Uniform(fTowerEdges[2], fTowerEdges[3]);
  
  eta = fTowerEta;
  phi = fTowerPhi;

  pt = energy / TMath::CosH(eta);

  // Time calculation for tower
  fTower->nTimes = 0;
  float tow_sumT = 0;
  float tow_sumW = 0;
   
  //double PVT; 
    // per particle timing cleaning
/*    fItVertexInputArray->Reset();
Candidate *pv = static_cast<Candidate*>(fItVertexInputArray->Next());
  PVT = pv->Position.T(); */

  float mostE_deposit = -100., Edeposit = -50., towerTime = 100000., Tdeposit;;
  for (unsigned int i = 0 ; i < fTower->ecal_E_t.size() ; i++) { 
    float w = TMath::Sqrt(fTower->ecal_E_t[i].first);
     Edeposit = fTower->ecal_E_t[i].first;
     Tdeposit = fTower->ecal_E_t[i].second;
    
     if (Edeposit > mostE_deposit) {
        mostE_deposit = Edeposit;
        towerTime = fTower->ecal_E_t[i].second;
        }
          //fDebugOutputCollector.fillContainer ("m_partType", fTower->Position.Z()) ;
    //cout<<fTower->ecal_E_t[i].first<<"\t"<<fTower->ecal_E_t[i].second<<endl;
    tow_sumT += w*fTower->ecal_E_t[i].second;
    tow_sumW += w;
    //tow_sumW += 1;
    
    fTower->nTimes++;
  }
    float mostE_eta;
    
  for (unsigned int j = 0; j < fTower->ecal_E_eta.size() ; j ++)    {
        if (fTower->ecal_E_eta[j].first == mostE_deposit) mostE_eta = fTower->ecal_E_eta[j].second;}
//        mostE_eta = fTower->ecal_E_t.find(mostE_deposit)->second;
        //fTower->Position.SetEta(mostE_eta);
mostE_deposit = 0.;
  //if (tow_sumW > 0.) {
 // eta = mostE_eta;
  
  if (tow_sumW > fTimingEMin) {
  float dt_tower = gRandom->Gaus(0.0, 15.);
  
    //fTower->Position.SetPtEtaPhiE(1.0, eta, phi,((tow_sumT/tow_sumW)));// + dt_tower));
    //fTower->Position.SetPtEtaPhiE(1.29, eta, phi,(towerTime) + dt_tower);
    fTower->Position.SetPtEtaPhiE(1.29, mostE_eta, phi,(towerTime + dt_tower));
    //fDebugOutputCollector.fillContainer3D ("m_partType:eta:E", fTower->Position.T(),fTower->Position.Eta(), fTower->Momentum.E()) ;
    //fTower->Position.SetT(tow_sumT/tow_sumW);
    //fDebugOutputCollector.fillContainer3D ("m_eta:time:pt", fTower->Position.T(), fTower->Momentum.E(), fTower->Momentum.Pt());
   /* fDebugOutputCollector.fillContainer3D ("tower_eta:time:E",
     fabs (fTower->Position.Eta ()),
     (fTower->Position.T ()) ,
     fTower->Momentum.E ()
    ) ; */
    
  } else {

    //fTower->Position.SetPtEtaPhiE(1.0,eta,phi,999999.);
    fTower->Position.SetPtEtaPhiE(1.29,eta,phi,999999.);
  }
  
//pt = energy / TMath::CosH(eta);

  //  fTower->Position.SetPtEtaPhiE(1.0, eta, phi, 0.);
  //fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
  fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
  fTower->Eem = ecalEnergy;
 // fTower->Ehad = hcalEnergy;
fDebugOutputCollector.fillContainer3D ("m_partType:eta:E", fTower->Position.T(),fTower->Position.Eta(), fTower->Momentum.E()) ; 
//time_vs_eta = new TProfile ("time_vs_eta","timespread_vs_eta",200, -2.,2.,3800.,12000.);
/*cout<<"hh"<<endl;
//TH2F *time_vs_eta = new TH2F ( "time_eta","time_eta", 200, -3,3,200,3000,12000);
time_vs_eta->Fill(fTower->Position.Eta(),fTower->Position.T());
TProfile * prof = time_vs_eta->ProfileY();
prof->Fit(corr_timeTower);

cout<<"tprofile"<<"\t"<<time_vs_eta->GetEntries()<<endl;*/
//time_vs_eta->Fit(corr_timeTower);
//time_eta->Delete();
//time_vs_eta->Reset();
//fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt", fTower->Position.Phi(), fTower->Position.T()*inv_c_light , fTower->Position.Eta()) ;

/*  fDebugOutputCollector.fillContainer3D ("tower_eta:time:E",
     fabs (fTower->Position.Eta ()),
     (fTower->Position.T ()) ,
     fTower->Momentum.E ()
    ) ; */

  fTower->Edges[0] = fTowerEdges[0];
  fTower->Edges[1] = fTowerEdges[1];
  fTower->Edges[2] = fTowerEdges[2];
  fTower->Edges[3] = fTowerEdges[3];

  // fill calorimeter towers and photon candidates
  if(energy > 0.0){
    if(fTowerPhotonHits > 0 && fTowerTrackHits == 0){      
      fPhotonOutputArray->Add(fTower);
    }
    
/*    fDebugOutputCollector.fillContainer3D ("tower_eta:time:E",
     fabs (fTower->Position.Eta ()),
     (fTower->Position.T ()) ,
     fTower->Momentum.E ()
    ) ; */
    
    fTowerOutputArray->Add(fTower);
  }

  // save all the tracks as energy flow tracks
  fItTowerTrackArray->Reset();
  while((track = static_cast<Candidate*>(fItTowerTrackArray->Next()))){
    fEFlowTrackOutputArray->Add(track);
  }

  ecalEnergy -= fTrackECalEnergy;
  if(ecalEnergy < 0.0) ecalEnergy = 0.0;

 // hcalEnergy -= fTrackHCalEnergy;
  //if(hcalEnergy < 0.0) hcalEnergy = 0.0;

//  energy = ecalEnergy + hcalEnergy;

    energy = ecalEnergy;
  // save ECAL and/or HCAL energy excess as an energy flow tower
  if(energy > 0.0){
    // create new tower
    tower = static_cast<Candidate*>(fTower->Clone());
    pt = energy / TMath::CosH(eta);
    tower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
    tower->Eem = ecalEnergy;
    //tower->Ehad = hcalEnergy;
    fEFlowTowerOutputArray->Add(tower);
   
  }
  
/*  
  
  // Time calculation for tower
  fTower->nTimes = 0;
  float tow_sumT = 0;
  float tow_sumW = 0;
   
  double PVT; 
    // per particle timing cleaning
    fItVertexInputArray->Reset();
Candidate *pv = static_cast<Candidate*>(fItVertexInputArray->Next());
  PVT = pv->Position.T();

  float mostE_deposit = -100., Edeposit = -50., towerTime = 100000., Tdeposit;;
  for (unsigned int i = 0 ; i < fTower->ecal_E_t.size() ; i++) { 
    float w = TMath::Sqrt(fTower->ecal_E_t[i].first);
     Edeposit = fTower->ecal_E_t[i].first;
     Tdeposit = fTower->ecal_E_t[i].second;
     /*if ( TMath::Abs(Tdeposit - Edeposit) > fTimeCut)   {
        energy -= Edeposit;
        Edeposit = 0.;
        } */
 /*    if (Edeposit > mostE_deposit) {
        mostE_deposit = Edeposit;
        towerTime = fTower->ecal_E_t[i].second;
        }
          //fDebugOutputCollector.fillContainer ("m_partType", fTower->Position.Z()) ;
    //cout<<fTower->ecal_E_t[i].first<<"\t"<<fTower->ecal_E_t[i].second<<endl;
    tow_sumT += w*fTower->ecal_E_t[i].second;
    tow_sumW += w;
    //tow_sumW += 1;
    
    fTower->nTimes++;
  }
    float mostE_eta;
    
  for (unsigned int j = 0; j < fTower->ecal_E_eta.size() ; j ++)    {
        if (fTower->ecal_E_eta[j].first == mostE_deposit) mostE_eta = fTower->ecal_E_eta[j].second;}
//        mostE_eta = fTower->ecal_E_t.find(mostE_deposit)->second;
        //fTower->Position.SetEta(mostE_eta);
mostE_deposit = 0.;
  //if (tow_sumW > 0.) {
  if (tow_sumW > fTimingEMin) {
  float dt_tower = gRandom->Gaus(0.0, 10.);
    //fTower->Position.SetPtEtaPhiE(1.0, eta, phi,((tow_sumT/tow_sumW)));// + dt_tower));
    //fTower->Position.SetPtEtaPhiE(1.29, eta, phi,(towerTime) + dt_tower);
    fTower->Position.SetPtEtaPhiE(1.29, mostE_eta, phi,(towerTime + dt_tower));
    //fDebugOutputCollector.fillContainer3D ("m_partType:eta:E", fTower->Position.T(),fTower->Position.Eta(), fTower->Momentum.E()) ;
    //fTower->Position.SetT(tow_sumT/tow_sumW);
    //fDebugOutputCollector.fillContainer3D ("m_eta:time:pt", fTower->Position.T(), fTower->Momentum.E(), fTower->Momentum.Pt());
   /* fDebugOutputCollector.fillContainer3D ("tower_eta:time:E",
     fabs (fTower->Position.Eta ()),
     (fTower->Position.T ()) ,
     fTower->Momentum.E ()
    ) ; */
    
/*  } else {

    //fTower->Position.SetPtEtaPhiE(1.0,eta,phi,999999.);
    fTower->Position.SetPtEtaPhiE(1.29,eta,phi,999999.);
  }
  
pt = energy / TMath::CosH(eta);

  //  fTower->Position.SetPtEtaPhiE(1.0, eta, phi, 0.);
  fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
  fTower->Eem = ecalEnergy;
  
*/

  
}
//------------------------------------------------------------------------------
void ECalorimeter::PileUpSubtractor()    {
Candidate *vertex, *tower, *track, *cal_track, *gen_vertex, *reco_vertex, *pvTower, *eflowtrack, *pvTrack, *notime_track;
double_t vertex_time, tower_time, lepton_energy = 0., ttime, corr_time = 8000.;
float_t vertex_energy, tmp_t, TOF;
int leptons=0;
//link tracks to towers
//fItTowerTrackArray->Reset();
    int numero = 0;
fItVertexingTrackInputArray = fVertexingTrackInputArray->MakeIterator();
fItVertexingTrackInputArray->Reset();
while((track = static_cast<Candidate*>(fItVertexingTrackInputArray->Next())))    {

        if(track->IsPU == 0. && (fabs(track->PID) == 11  || fabs(track->PID)==13))   {
         lepton_energy += track->Momentum.E(); 
         ++leptons; }
 //cal_track = static_cast<Candidate*> (track->Clone ()) ;    
    float track_eta = track->Position.Eta();
    float track_phi = track->Position.Phi();

//fDebugOutputCollector.fillContainer2D ("Nm_partType:ID",track->Position.Eta(),track->Position.T() ) ;
    //fItTowerOutputArray->Reset();
    //while((tower = static_cast<Candidate*>(fItTowerOutputArray->Next()))){
            
            fItEFlowTowerOutputArray->Reset();
    while((tower = static_cast<Candidate*>(fItEFlowTowerOutputArray->Next()))){
            
        //float tower_eta = tower->Position.Eta();
        //float tower_phi = tower->Position.Phi();
        if(track_eta < tower->Edges[1] && track_eta > tower->Edges[0] && track_phi < tower->Edges[3] && track_phi > tower->Edges[2]  && track->Momentum.Pt()> 2. && track->Momentum.P()>5. && tower->Position.T()<999099.){// && tower->ecal_E_t.size()== 1 ) {
        ++numero;
      // if(track->Momentum.Pt()> 2. && track->Momentum.P()>10.){
            //tower->VertexID_gen = track->VertexID_gen;
            cal_track = static_cast<Candidate*> (track->Clone ()) ;    
            cal_track->Position.SetT(tower->Position.T()); // in ps
            //fDebugOutputCollector.fillContainer2D ("Nm_partType:ID", cal_track->Momentum.Pt(),cal_track->Position.T() ) ;
            //fDebugOutputCollector.fillContainer4D ("tower_eta:time:E:res", cal_track->Position.T(),cal_track->IsPU, tower->Position.T(), tower->Position.T()*inv_c_light); 
            fCalTrackArray->Add (cal_track) ;
            //if ( track->Charge != 0)    cout<<"pid tower:   "<<tower->Charge<<endl;
            }         
    }        
}
cout<<numero<<endl;
 double gen_pv_z, vtx_t;
fItVertexInputArray->Reset();
/*while((vertex = static_cast<Candidate*>(fItVertexInputArray->Next()))) {
    if (vertex->IsPU == 0)  {   
        gen_pv_z = vertex->Position.Z();
        vtx_t = vertex->Position.T()*inv_c_light;}
        }
  */       
 fPVItInputArray->Reset();
Candidate *pv_track = static_cast<Candidate*>(fPVItInputArray->Next());
  if (pv_track) {
    vtx_t = pv_track->Position.T()*inv_c_light;
    //gen_pv_z = pv_track->Position.Z();
    //gen_pv_z += gRandom->Gaus(0.,0.065);
    }
 
 fItCalTrackArray = fCalTrackArray->MakeIterator();
    fItCalTrackArray->Reset();
    while((track = static_cast<Candidate*>(fItCalTrackArray->Next())))  { 
        //fDebugOutputCollector.fillContainer2D ("Nm_partType:ID", track->Position.Z(),track->Position.T() ) ;
        
        fItVertexInputArray->Reset();
        while((vertex = static_cast<Candidate*>(fItVertexInputArray->Next())))  {
//        vtx_t = vertex->Position.T() * inv_c_light;
//        gen_pv_z = vertex->Position.Z();
            if(track->VertexID_gen == vertex->VertexID_gen) {
           // if(track->IsPU == 0)       
                   float R = (track->Momentum.Pt()/(0.3*3.8));//Bz; m
                   //secant with respect to 0.0.0
                   float s = (TMath::Sqrt(TMath::Power(track->Position.X(),2)+TMath::Power(track->Position.Y(),2)));// + TMath::Power(track->Position.Z()-pv_z,2)));
                   float alpha = TMath::ASin(s/(1000*2*R));
                   float arc = 2 * alpha * R * 1000;
                   double eff_velocity = TMath::Sqrt((TMath::Power(track->Momentum.Py(),2) + TMath::Power(track->Momentum.Pz(),2) + TMath::Power(track->Momentum.Px(),2))) / tower->Momentum.E();
                        double inv_v = TMath::Power(eff_velocity, -1);
                   float z_distance = track->Position.Z() - vertex->Position.Z();
                   //float z_distance = track->Position.Z() - gen_pv_z;
                   TOF = (TMath::Sqrt((TMath::Power(z_distance,2) + TMath::Power(arc,2)))) ; 
                         
                    corr_time = track->Position.T() - TOF*inv_c_light;
                   ttime = track->Position.T();
                   //track->Position.SetT(corr_time);
                  //fDebugOutputCollector.fillContainer4D ("tower_eta:time:E:res", track->Position.T(),track->IsPU, corr_time - vtx_t, track->Momentum.Pt());  
                   }        
        }       
        track->Position.SetT(corr_time);
 fDebugOutputCollector.fillContainer3D ("Nm_partType:ID:IsPU", track->Position.T() - vtx_t,track->IsPU, track->Position.T() ) ;
//fDebugOutputCollector.fillContainer4D ("tower_eta:time:E:res", track->Position.T(),track->IsPU, corr_time-vtx_t, track->Momentum.Pt()); 
}
 

//set vertices time and energy
double track_molt = 0.;
vertex_time = 0.;
float sum_track_time = 0.;
float track_time = 0.;
vertex_energy = 0.;
float track_energy = 0.;
int vertexnumber = 0., j;
double sumpt2 = 0.;
double track_pt2, time_resolution,timing_pt2, gen_time, track_eta, track_phi;
Candidate *tower_candidate, *tower_candidate2, *tower_candidate3, *notime_tower;

fItVertexInputArray->Reset();
while((vertex = static_cast<Candidate*>(fItVertexInputArray->Next())))  {
    int vertexID = vertex->VertexID_gen;
    gen_time = vertex->Position.T()*inv_c_light;
    ++vertexnumber;
    fItCalTrackArray = fCalTrackArray->MakeIterator();
    fItCalTrackArray->Reset();
    while((track = static_cast<Candidate*>(fItCalTrackArray->Next())))  { 
        int trackID = track->VertexID_gen;      
       // fDebugOutputCollector.fillContainer2D ("Nm_partType:ID",track->Momentum.Pt(),track->Position.T() ) ;
        if(trackID == vertexID && track->Momentum.Pt() > 4. && track->Momentum.P()>15.){ // && fabs(track->Position.T())<3. && fabs(track->Position.T()>1.479)) 
        //fDebugOutputCollector.fillContainer2D ("Nm_partType:ID",track->Momentum.Pt(),track->Position.T() ) ;
        //++track_molt;
   // fDebugOutputCollector.fillContainer ("Nm_partType",gen_time ) ;
            track_time = track->Position.T();
            track_energy = track->Momentum.E();
            track_pt2 = track->Momentum.Pt();
            vertex_energy += track_energy;
                //cout<<"inside:    "<<track_time <<endl;// - gen_time<<endl;       
            if (TMath::Abs(track_time - gen_time)< 4000.)   {
                sum_track_time += track_time;// * TMath::Power(track_pt2,4) ;
                timing_pt2 += TMath::Power(track->Momentum.Pt(),4);
                
                ++track_molt;   
           
                sumpt2 += track_pt2;}
               
               }//fDebugOutputCollector.fillContainer2D ("Nm_partType:ID",TOF,track->Position.T() ) ;
   //}
    } 
    vertex_time = (sum_track_time)/(track_molt);//(timing_pt2);
             
    if(track_molt > 0. && TMath::Abs(vertex_time)<1700.) { // && TMath::Abs(vertex_time)<640.)    {
    reco_vertex = static_cast<Candidate*>(vertex->Clone());
    //vertex_time = (sum_track_time)/(timing_pt2);
    reco_vertex->Position.SetT(vertex_time);
    reco_vertex->Momentum.SetE(vertex_energy);
    //reco_vertex->Momentum.SetPerp(sumpt2);
//    reco_vertex->
    reco_vertex->Momentum.SetPtEtaPhiE(sumpt2,0.,0.,vertex_energy);
    //if(TMath::Abs(reco_vertex->Position.T())<2000.)  {
        //cout<<"inside"<<endl;       
    fRecoVertexArray->Add(reco_vertex);
    //fDebugOutputCollector.fillContainer3D ("Nm_partType:ID:IsPU",reco_vertex->Position.T(),track_molt, reco_vertex->IsPU ) ;
    //}
    }
    sum_track_time = 0.;
    vertex_energy = 0.;
    sumpt2 = 0.;
    track_molt =0.;
    timing_pt2=0.;
    }
    
fItRecoVertexArray = fRecoVertexArray->MakeIterator();        
fItRecoVertexArray->Reset();    
while((reco_vertex = static_cast<Candidate*>(fItRecoVertexArray->Next())))    {
    //cout<<reco_vertex->IsPU<<endl;
    //fDebugOutputCollector.fillContainer ("m_partType", reco_vertex->Position.T()) ;
    //wrap = 0;
fItVertexInputArray = fVertexInputArray->MakeIterator();
    fItVertexInputArray->Reset();
    while((gen_vertex = static_cast<Candidate*>(fItVertexInputArray->Next())))  {
       // fDebugOutputCollector.fillContainer ("DR",gen_vertex->Position.T() );
        if(reco_vertex->VertexID_gen == gen_vertex->VertexID_gen)    {
        
        fDebugOutputCollector.fillContainer4D ("tower_eta:time:E:res",
        reco_vertex->Position.T(),
        //reco_vertex->Momentum.E(),
        reco_vertex->IsPU,
        reco_vertex->Position.T() - gen_vertex->Position.T()*inv_c_light,
        //reco_vertex->VertexID_gen,
        //reco_vertex->Momentum.E()
        gen_vertex->IsPU
        //vertex_energy
        ) ;
        }

    }

    

    
            
}  
//fRecoVertexArray->Clear(); 
//cout<<vertexnumber<<endl;
 

//chose primary vertex 
double tmp_eem = 0., tmp_pt2 = -90.;
double vertex_eem = 0. , vertex_pt2 = 0.;
int pvID = 1, pvIsPU;
double pv_t = 19999.;
double pv_energy = 0.;
double energy_left = 0.;
float gen_pv_time;
//fDebugOutputCollector.fillContainer ("DR",tmp_pt2 ); 
fItRecoVertexArray->Reset();
while((vertex = static_cast<Candidate*>(fItRecoVertexArray->Next())))  {
           
//fDebugOutputCollector.fillContainer2D ("Nm_partType:ID",vertex->Position.T(), vertex->IsPU ) ;
    //vertex_eem = vertex->Momentum.E();
    vertex_pt2 = TMath::Power(vertex->Momentum.Pt(),2);
    //cout<<vertex->Momentum.Pt()<<endl;
    //vertex_eem = vertex->sumPtSquare;
    //cout<<"v_energy:"<<"\t"<<vertex_eem<<endl;
 /*   if(vertex_eem > tmp_eem)    {
        tmp_eem = vertex_eem;
        pvID = vertex->VertexID_gen;
        pvIsPU = vertex->IsPU;
      //  cout<<"pvispu:"<<"\t"<<pvIsPU<<endl;
        pv_t = vertex->Position.T();
        
        } */
       //fDebugOutputCollector.fillContainer ("DR",vertex_pt2 ); 
       //fDebugOutputCollector.fillContainer2D ("DR_tmp",vertex_pt2 , tmp_pt2);
           if(vertex_pt2 > tmp_pt2)    {
          //     cout<<"eo"<<endl;
        tmp_pt2 = vertex_pt2;
        pvID = vertex->VertexID_gen;
        //pvID = 10;
        pvIsPU = vertex->IsPU;
      //  cout<<"pvispu:"<<"\t"<<pvIsPU<<endl;
        pv_t = vertex->Position.T();
        pv_z = vertex->Position.Z();
        
        }
        
    } 
    
float gen_z, gen_x, gen_y;
int flag_pv = 0;  
    
fItRecoVertexArray = fRecoVertexArray->MakeIterator();               
fItRecoVertexArray->Reset();
while((vertex = static_cast<Candidate*>(fItRecoVertexArray->Next())))  {  
  if(vertex->VertexID_gen == pvID)  {
    vertex->Position.SetT(pv_t);
    }
    //}
    /*float gen_z, gen_x, gen_y;
    int flag_pv = 0;*/
fItVertexInputArray = fVertexInputArray->MakeIterator();
fItVertexInputArray->Reset();
    while((gen_vertex = static_cast<Candidate*>(fItVertexInputArray->Next())))  {
    //if(pvID == gen_vertex->VertexID_gen && flag_pv == 0) {
    if(gen_vertex->VertexID_gen == vertex->VertexID_gen && flag_pv == 0 ) { 
    cout<<"inside"<<endl;       
    flag_pv = 1;
     gen_z = gen_vertex->Position.Z();
     gen_y = gen_vertex->Position.Y();
     gen_x = gen_vertex->Position.X();
     gen_pv_time = gen_vertex->Position.T()*inv_c_light;
    fDebugOutputCollector.fillContainer4D ("eta:time:pt:PID",
                  //tmp_pt2,
                  vertex->Position.T() - gen_vertex->Position.T()*inv_c_light,
                  gen_vertex->IsPU,
                  gen_vertex->VertexID_gen,
                  gen_vertex->Position.T()*inv_c_light);
                  /*pv_t,
                  pvIsPU,
                  tmp_pt2); */
                //  fDebugOutputCollector.fillContainer4D ("tower_eta:time:E:res", vertex->Position.T()-gen_pv_time,vertex->IsPU, gen_vertex->IsPU, vertex->Momentum.Pt()/gen_vertex->Momentum.Pt()); 
}
}
}
flag_pv = 0;
double err_vertex_time,std;
int n = 0;

/*fVertexInputArray->MakeIterator();
fItVertexInputArray->Reset();
while((gen_vertex = static_cast<Candidate*>(fItVertexInputArray->Next())))  {
if ( gen_vertex->IsPU == 0) {
         gen_pv_time = gen_vertex->Position.T()*inv_c_light;
        gen_z = gen_vertex->Position.Z(); }
fDebugOutputCollector.fillContainer4D ("eta:time:pt:PID",
                  //tmp_pt2,
                  gen_vertex->Position.T()*inv_c_light,
                  gen_pv_time,
                  gen_vertex->IsPU,
                  gen_vertex->Position.Z()); 
                  }
 */
                  
/*fItCalTrackArray->Reset();
while((track = static_cast<Candidate*>(fItCalTrackArray->Next())))  { 
    if(track->VertexID_gen == pvID ){//&& fabs(track->Position.T())<4000) {
        std = TMath::Power((track->Position.T()-pv_t),2);
        ++n;}
        }
        err_vertex_time = std/n;
//fDebugOutputCollector.fillContainer2D ("DR_tmp",n , track->Position.T());}
fRecoVertexArray->Clear();       
              */

              
double neutralTOF, neutr_corrTime;
int nflag = 10;
//pv_z = 0;
pv_z = gen_z;
fItEFlowTowerOutputArray = fEFlowTowerOutputArray->MakeIterator();
fItEFlowTowerOutputArray->Reset();
//fItPhotonOutputArray = fPhotonOutputArray->MakeIterator();
//fItPhotonOutputArray->Reset();
    while((tower = static_cast<Candidate*>(fItEFlowTowerOutputArray->Next()))){
    //while((tower = static_cast<Candidate*>(fItPhotonOutputArray->Next()))){
    
    
    
    
    
//     fDebugOutputCollector.fillContainer3D ("tower_eta:time:pt", tower->Position.T(), tower->Momentum.E(), tower->Momentum.Pt());
 //   fDebugOutputCollector.fillContainer ("m_partType",tower->IsPU) ;
         if (tower->Position.T()==999999.) {
         //cout<<"notimetower"<<endl,
            notime_tower = static_cast<Candidate*>(tower->Clone());
            fnotimeTowerOutputArray->Add(notime_tower);
            }
            
    //fDebugOutputCollector.fillContainer3D ("DR_tmp",tower->Position.Eta() , tower->Position.Phi(), tower->Position.T());
   // fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt", 
     //          tower->Position.Phi(), tower->Position.T()*inv_c_light , tower->Position.Eta()) ;
    
//    while((tower = static_cast<Candidate*>(fItTowerOutputArray->Next()))){
  /*  fItCalTrackArray = fCalTrackArray->MakeIterator();    
    fItCalTrackArray->Reset();
    while((track = static_cast<Candidate*>(fItCalTrackArray->Next())))    {
        float cal_track_eta = track->Position.Eta();
        float cal_track_phi = track->Position.Phi();
// fDebugOutputCollector.fillContainer ("m_partType", track->Position.T()) ;
//    fItEFlowTowerOutputArray->Reset();
  //  while((tower = static_cast<Candidate*>(fItEFlowTowerOutputArray->Next())))

        if(cal_track_eta < tower->Edges[1] && cal_track_eta > tower->Edges[0] && cal_track_phi < tower->Edges[3] && cal_track_phi > tower->Edges[2]){// && cal_track->Momentum.Pt()> 2. && cal_track->Momentum.Pt()>2.)   
//        nflag = 1;
        if(tower->Position.T()<99999.){
        nflag = 1;
        //fDebugOutputCollector.fillContainer3D ("DR_tmp",tower->Position.Z() , tower->Position.X(), tower->Position.T());//tower
        tower_candidate = static_cast<Candidate*>(tower->Clone());
        tower_candidate->Position.SetT(track->Position.T());
        //fDebugOutputCollector.fillContainer3D ("DR_tmp",tower->Position.T() , tower->Position.Eta(), nflag);//tower->Position.T());
        //tower->Position.SetT(60.);
        fDebugOutputCollector.fillContainer3D ("DR_tmp",        tower_candidate->Charge, tower_candidate->Position.T() , tower_candidate->Position.Eta()) ;
        fPvEFlowTowerCandidateArray->Add(tower_candidate);
        }
        
        }
        /*else if (tower->Position.T()==999999.) {
            notime_tower = static_cast<Candidate*>(tower->Clone());
            fnotimeTowerOutputArray->Add(notime_tower);
            } */
  //  } 
    
   //     else    
//fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt",        nflag, tower->Position.T()*inv_c_light , tower->Position.Eta()) ;
   if(nflag !=100000000) {
   //fDebugOutputCollector.fillContainer3D ("DR_tmp",tower->Position.Eta() , nflag, tower->Position.T());

   //fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt",        tower->Position.Phi(), tower->Position.T()*inv_c_light , tower->Position.Eta()) ;
                       //double hTower = 1290. * TMath::Tan ( TMath::Pi()/2. - 2*TMath::ATan(exp(- tower->Position.Eta())));
                        if(TMath::Abs(tower->Position.Eta())<1.6)  {
                        //if(hTower>0)  {
                        nflag = 2; 
                       // fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt",        tower->Position.Phi(), tower->Position.T()*inv_c_light , tower->Position.Eta()) ;
                        //if(tower->Position.Eta()> -1.65 && tower->Position.Eta() < 1.65 )  {
                        //neutralTOF = 1. * inv_c_light * TMath::Sqrt(TMath::Power(1290,2)+TMath::Power(1290* TMath::Tan(0.5*TMath::Pi() - (2*TMath::ATan(exp(- tower->Position.Eta())))),2));
                        double neutrTowerT = tower->Position.T();
                        //double xTower = 1290 * cos(2*TMath::ATan(exp(- tower->Position.Eta())));
                        double xTower = 1.29 * 1000. * TMath::Cos ( tower->Position.Phi());
                        //double yTower = 1.29*1000 * sin(2*TMath::ATan(exp(- tower->Position.Eta())));
                        //double hTower = 1290. * TMath::Tan ( TMath::Pi()/2. - 2*TMath::ATan(exp(- tower->Position.Eta())));
                        double yTower = 1290. * TMath::Sin ( tower->Position.Phi());
                        //double hTower = 1.29 * 1000 * tan(TMath::Pi() - (2*TMath::ATan(exp(- tower->Position.Eta()))));
                        //double hTower = 1290. * TMath::Cos ( TMath::ATan(exp(- tower->Position.Eta()))) / TMath::Sin ( TMath::ATan(exp(- tower->Position.Eta())));
                        double hTower = 1290. * TMath::Tan ( TMath::Pi()/2. - 2*TMath::ATan(exp(- tower->Position.Eta())));
                        //double hTower = tower->Position.Z()*1000*1.3;
                        //double eff_velocity = TMath::Sqrt(TMath::Power(tower->Momentum.Py(),2)+ TMath::Power( tower->Momentum.Pz(),2) + TMath::Power(tower->Momentum.Px(),2)) / tower->Momentum.E();
                        //double inv_v = TMath::Power(eff_velocity, -1);
                         //neutralTOF = 10E2 * TMath::Sqrt(TMath::Power(tower->Position.Z(),2)+TMath::Power(tower->Position.X(),2)+TMath::Power(tower->Position.Y(),2));
   //                      fDebugOutputCollector.fillContainer4D ("Nm_eta:time:pt:size",       tower->ecal_E_t.size(), hTower - pv_z  , tower->Position.Eta(), hTower) ;
  // pv_z = 0.;
                          //neutralTOF = TMath::Sqrt(TMath::Power(xTower,2)+TMath::Power(yTower,2)+TMath::Power((hTower -gen_z),2))*inv_c_light;                       
                        //pv_z =0;
                        //cout<<double(tower->Momentum.P()/tower->Momentum.E())<<endl;
                        //hTower += gRandom->Gaus(0.,1.);
                        double hit_pos = TMath::Sqrt(TMath::Power(xTower,2)+TMath::Power(yTower,2)+TMath::Power((hTower),2));
                        hit_pos += gRandom->Gaus(0.,1.);
                        neutralTOF = TMath::Sqrt(TMath::Power(xTower - gen_x ,2)+TMath::Power(yTower - gen_y,2)+TMath::Power((hTower - gen_z),2))*inv_c_light;
                        double vert_pos  = TMath::Sqrt(TMath::Power( gen_x ,2)+TMath::Power( gen_y,2)+TMath::Power(( pv_z),2));
                        //neutralTOF = TMath::Sqrt(TMath::Power(hit_pos,2) - TMath::Power(vert_pos,2)) * inv_c_light;
                        //neutralTOF = abs(1.29 * inv_c_light * TMath::Power ( TMath::Sin(2*TMath::ATan(exp(- tower->Position.Eta()))),-1));
                        //neutralTOF = TMath::Sqrt(1.29*1.29+TMath::Power((hTower - pv_z),2))*inv_c_light;
                        //cout<<"c:"<<"\t"<<inv_c_light<<"v"<<"\t"<<eff_velocity<<endl;
                        neutr_corrTime = neutrTowerT - neutralTOF;
                        
                       //  fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt",       neutr_corrTime, neutrTowerT , neutralTOF) ;
                        //if(tower->Position.T()<999999.)  {
                        tower_candidate = static_cast<Candidate*>(tower->Clone());                      
                       tower_candidate->Position.SetT(neutrTowerT - neutralTOF) ;
                       fPvEFlowTowerCandidateArray->Add(tower_candidate);
                       //TProfile *hprof  = new TProfile("hprof","Profile of pz versus px",100,-3,3,-100,100);
                      //fDebugOutputCollector.fillContainer4D ("Nm_eta:time:pt:size",      neutrTowerT, tower_candidate->Position.Eta() , neutralTOF, hTower) ;
                        time_vs_eta->Fill(neutrTowerT - neutralTOF, tower_candidate->Position.Eta());
                      // fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt",        tower_candidate->Position.Phi(), tower_candidate->Position.T()*inv_c_light , tower_candidate->Position.Eta()) ;
                      // fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt",        tower->Position.Phi(), tower->Position.T()*inv_c_light , tower->Position.Eta()) ;
                      //if(tower_candidate->ecal_E_t.size() == 1){
                      //fDebugOutputCollector.fillContainer4D ("Nm_eta:time:pt:size",      gen_pv_time, tower_candidate->IsPU , neutralTOF, tower_candidate->Position.T() - gen_pv_time) ;}
                        //tower->Position.SetT(neutrTowerT - neutralTOF);
                        //fDebugOutputCollector.fillContainer ("m_partType", tower->Position.T()) ;
                      //  fDebugOutputCollector.fillContainer3D ("DR_tmp",tower_candidate->Position.T() , neutr_corrTime, neutrTowerT);//tower->Position.T());
                        //}
                    }
                    //TMath::Abs(tower->Position.Eta())< 3. &&
                        else  if ( TMath::Abs(tower->Position.Eta())>1.6 && nflag !=2) {
                        //else  if ( (hTower<3000 || hTower == 3000) && nflag !=2) {
                        float effR;
                        //else  if ( tower->Position.Eta()<-1.479 && nflag != 2 ) { //|| tower->Position.Eta()<-1.479) && nflag !=2) {
                      //else if ( nflag != 1 && nflag!=2) {
                      //float effR = 3000. - pv_z ;
                      if ( tower->Position.Eta() > 0. || tower->Position.Eta()==0.) {
                       effR = 3000. - pv_z ;
                      }
                      else  {
                                 effR = -3000. - pv_z ;
                                }
                      
                      double zTower = 3000.;
                      double xTower = 3000. * TMath::Tan ( (2*TMath::ATan(exp(- tower->Position.Eta()))));
                      double yTower = xTower * TMath::Tan ( tower->Position.Phi());
                      //xTower += gRandom->Gaus(0.,1.);
                      //yTower += gRandom->Gaus(0.,1.);
                      //cout<<"PV_Z:"<<"\t"<<pv_z<<endl;
                            //neutralTOF = 1000. * inv_c_light*TMath::Sqrt(9.+ 9 * TMath::Power(TMath::Tan((2*TMath::ATan(exp(- tower->Position.Eta())))),2));
                            //neutralTOF = 1. * inv_c_light*TMath::Sqrt(TMath::Power(effR,2) + TMath::Sqrt(TMath::Power(effR,2)*TMath::Power(TMath::Tan((2*TMath::ATan(exp(- tower->Position.Eta())))),2)) + TMath::Power(effR * TMath::Tan((2*TMath::ATan(exp(- tower->Position.Eta()))) )* TMath::Tan(tower->Position.Phi()) ,2  ));
                            //neutralTOF = TMath::Sqrt(TMath::Power(xTower,2)+TMath::Power(yTower,2)+TMath::Power((zkTower - pv_z),2))*inv_c_light;
                            //fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt",        tower->Position.Phi(), tower->Position.T() , tower->Position.Eta()) ;
                            //neutralTOF = (3000. - pv_z) * inv_c_light;
                            neutralTOF = TMath::Abs(effR *    inv_c_light * TMath::Power ( TMath::Cos((2*TMath::ATan(exp(- tower->Position.Eta())))),-1));
                            //neutralTOF = (3000.) *    TMath::Power(tower->Momentum.Pz()/tower->Momentum.M(),-1) * TMath::Power ( TMath::Cos((2*TMath::ATan(exp(- tower->Position.Eta())))),-1);
                        //tower->Position.SetT(150.);
                        double neutrTowerT = tower->Position.T();  //in ps
                        double neutr_corrTime = neutrTowerT - neutralTOF;
                         //fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt",       neutr_corrTime, neutrTowerT  , tower->Momentum.E()) ;
                        //tower->Position.SetT(neutrTowerT - neutralTOF) ;
                        //if(tower->Position.T()<9999999.)  {
                        //fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt",        neutr_corrTime, tower->Position.T() , tower->Position.Eta()) ;
                        tower_candidate = static_cast<Candidate*>(tower->Clone());
                       //tower_candidate->Position.SetXYZT(xTower,yTower,zTower * tower->Position.Eta()/abs(tower->Position.Eta()),neutrTowerT - neutralTOF) ;
                       tower_candidate->Position.SetT(neutrTowerT - neutralTOF) ;
                       //if(tower_candidate->ecal_E_t.size() == 1){
                       //fDebugOutputCollector.fillContainer4D ("Nm_eta:time:pt:size",        tower_candidate->Position.Eta(), tower_candidate->Position.T() - gen_pv_time, tower_candidate->IsPU, effR) ;}
                       
                       fPvEFlowTowerCandidateArray->Add(tower_candidate);
                         //fDebugOutputCollector.fillContainer3D ("Nm_partType:ID:IsPU", tower->Position.T(),neutr_corrTime , nflag) ;
                       //  fDebugOutputCollector.fillContainer3D ("DR_tmp",tower_candidate->Position.T() , neutr_corrTime, neutrTowerT);//tower->Position.T());

                               //fDebugOutputCollector.fillContainer3D ("DR_tmp",neutralTOF , neutr_corrTime, neutrTowerT);//tower->Position.T());
                                //fDebugOutputCollector.fillContainer3D ("DR_tmp",tower->Position.Z() , tower->Position.X(), neutrTowerT);//tower->Position.T());
                                //fDebugOutputCollector.fillContainer4D ("Nm_eta:time:pt:size",      neutrTowerT, tower_candidate->Position.T() , tower_candidate->IsPU, neutralTOF) ;
                                //}
                                }
                        neutr_corrTime=0;
                        neutralTOF =0;
                        //nflag = 10;
                        }
                        nflag = 10;
                    }
                       // fDebugOutputCollector.fillContainer3D ("DR_tmp",tower->Position.T() , neutralTOF, nflag);//tower->Position.T());
                        //nflag = 10;
            
//}
//fDebugOutputCollector.fillContainer2D ("DR_tmp",tower->Position.Eta() , tower->Position.T());
//}



/*fItPvEFlowTowerCandidateArray = fPvEFlowTowerCandidateArray->MakeIterator();
fItPvEFlowTowerCandidateArray->Reset();
    while((tower = static_cast<Candidate*>(fItPvEFlowTowerCandidateArray->Next()))){
        fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt", 
               tower->Position.Phi(), tower->Position.T()*inv_c_light , tower->Position.Eta()) ;
               }
 */

fItnotimeTowerOutputArray = fnotimeTowerOutputArray->MakeIterator();
fItnotimeTowerOutputArray->Reset();
while((tower = static_cast<Candidate*>(fItnotimeTowerOutputArray->Next())))  {
    fItVertexingTrackInputArray = fVertexingTrackInputArray->MakeIterator();
    fItVertexingTrackInputArray->Reset();
while((track = static_cast<Candidate*>(fItVertexingTrackInputArray->Next())))    {
     if(track_eta < tower->Edges[1] && track_eta > tower->Edges[0] && track_phi < tower->Edges[3] && track_phi > tower->Edges[2])   {
        notime_track = static_cast<Candidate*>(track->Clone());
        fnotimeTrackOutputArray->Add(notime_track);
        }
    }
}



double neutr_energy = 0. , deposit_time=0.,pileup_energy = 0., PVT ;
int k = 0;
fItPvEFlowTowerCandidateArray = fPvEFlowTowerCandidateArray->MakeIterator();
fItPvEFlowTowerCandidateArray->Reset();

fPVItInputArray->Reset();
Candidate *pv = static_cast<Candidate*>(fPVItInputArray->Next());
  if (pv) PVT = pv->Position.T()*inv_c_light;
PVT += gRandom->Gaus(0.,30.);
gen_pv_time = PVT;
while((tower = static_cast<Candidate*>(fItPvEFlowTowerCandidateArray->Next())))  {
fDebugOutputCollector.fillContainer4D ("Nm_eta:time:pt:size",       tower->ecal_E_t.size(), tower->Position.T() - PVT , tower->Position.Eta() , tower->IsPU) ;
//fDebugOutputCollector.fillContainer3D ("DR_tmp",gen_pv_time , tower->IsPU, tower->Position.T());               
//fDebugOutputCollector.fillContainer3D ("DR_tmp",tower->Position.Eta() , tower->Position.Phi(), tower->Position.T());
//fItEFlowTowerOutputArray = fEFlowTowerOutputArray->MakeIterator();
//fItEFlowTowerOutputArray->Reset();
//while((tower = static_cast<Candidate*>(fItEFlowTowerOutputArray->Next())))  {
//fDebugOutputCollector.fillContainer3D ("DR_tmp",tower->Position.T() , tower->Position.Eta(), nflag);//tower->Position.T());
//            double neutralTOF = TMath::Sqrt(TMath::Power(tower->Position.Z(),2)+TMath::Power(tower->Position.X(),2)+TMath::Power(tower->Position.Y(),2));
            
//fDebugOutputCollector.fillContainer2D ("DR_tmp",tower->Position.Eta() , tower->Position.T());
 double corrTimeCut;
              //if(fabs(tower->Position.T() - pv_t) < 10.) {
              //if(TMath::Abs(tower->Position.T() - pv_t) < fTimeCut) {
             /*   if ( TMath::Abs(tower->Position.Eta()) < 0.5) fTimeCut = 0.9*fTimeCut;
                if ( TMath::Abs(tower->Position.Eta()) < 1. &&  TMath::Abs(tower->Position.Eta()>0.5)) fTimeCut = 8*fTimeCut;                             
                if ( TMath::Abs(tower->Position.Eta()) < 1.479 &&  TMath::Abs(tower->Position.Eta()>1.)) fTimeCut = 95.; 
                if ( TMath::Abs(tower->Position.Eta()) > 1.479) fTimeCut = 1.75*fTimeCut; */
                /*if(TMath::Abs(tower->Position.Eta())<1.59 || TMath::Abs(tower->Position.Eta()) == 1.59) {
                corrTimeCut = 1.8 * fTimeCut * corr_timeTower->Eval(tower->Position.Eta()) + 20.;
                }
                else if ( tower->Position.Eta()>1.59 )
                {
                    corrTimeCut = 2 * corr_timeTower2->Eval(TMath::Abs(3. - tower->Position.Eta())) + 20.;
                    }
                else if ( tower->Position.Eta()<-1.59)
                {
                    corrTimeCut = 2 * corr_timeTower2->Eval(TMath::Abs(-3. - tower->Position.Eta())) + 20.;
               }     
                   // if ( tower->ecal_E_t.size() == 1){
                //fDebugOutputCollector.fillContainer4D ("Nm_eta:time:pt:size",      tower->Position.T(), tower->Position.Eta() , (tower->Position.T() - gen_pv_time), tower->IsPU) ;}
//                if ( TMath::Abs(tower->Position.Eta()) > 1.56)  corrTimeCut = 90;*/

                corrTimeCut = fTimeCut;
                //if(tower->ecal_E_t.size() > 1)  corrTimeCut = 300.;
                                //if ( TMath::Abs(tower->Position.Eta()) > 1.6)  {
                                  //  corrTimeCut = 60;}
              if(TMath::Abs(tower->Position.T() - gen_pv_time) < corrTimeCut) {
              //fDebugOutputCollector.fillContainer3D ("DR_tmp",tower->IsPU, tower->Momentum.E(), tower->Position.T());
              //if((tower->Position.T() - pv_t) < 10. && (tower->Position.T() - pv_t) > -10.) {
            //  fDebugOutputCollector.fillContainer3D ("DR_tmp",tower->Position.Eta() , tower->Momentum.E(), tower->Position.T());
              //if(fabs(tower->Position.T() - pv_t) < err_vertex_time) {
              //cout<<"deltaT"<<endl;
              pvTower = static_cast<Candidate*>(tower->Clone());
              fPvEFlowTowerOutputArray->Add(pvTower);
              }
 }             
fItPvEFlowTowerOutputArray = fPvEFlowTowerOutputArray->MakeIterator();                
fItPvEFlowTowerOutputArray->Reset();                
//cout<<"here"<<endl;
while((pvTower = static_cast<Candidate*>(fItPvEFlowTowerOutputArray->Next())))  {
 //fDebugOutputCollector.fillContainer3D ("DR_tmp",gen_pv_time , pvTower->IsPU, pvTower->Position.T());               
    if(pvTower->IsPU == 0){
    neutr_energy += pvTower->Momentum.E();  }
    else if(pvTower->IsPU == 1) {
        pileup_energy += pvTower->Momentum.E();}
 //fDebugOutputCollector.fillContainer3D ("DR_tmp",pvTower->Position.Eta() , pvTower->IsPU, pvTower->Position.T());               
                }
double chs_neutr_energy, chs_pileupenergy;
fItEFlowTowerOutputArray = fEFlowTowerOutputArray->MakeIterator();                
fItEFlowTowerOutputArray->Reset();                         
while((pvTower = static_cast<Candidate*>(fItEFlowTowerOutputArray->Next())))  {
if(pvTower->IsPU == 0){
    chs_neutr_energy += pvTower->Momentum.E();  
                }
            else if (pvTower->IsPU == 1){
                chs_pileupenergy += pvTower->Momentum.E();
                }     
                }
//fDebugOutputCollector.fillContainer ("m_partType", neutr_energy/chs_neutr_energy) ;
fDebugOutputCollector.fillContainer3D ("DR_tmp",neutr_energy/chs_neutr_energy , pileup_energy/chs_pileupenergy, (neutr_energy + pileup_energy)/chs_neutr_energy);               
//fDebugOutputCollector.fillContainer3D ("DR_tmp",neutr_energy/gen_energy , chs_neutr_energy/gen_energy, (neutr_energy + pileup_energy)/(gen_energy + pu_gen_energy));
//fItCalTrackArray->Reset();                                
//while((track = static_cast<Candidate*>(fItCalTrackArray->Next()))) {
fItVertexingTrackInputArray->Reset();
while((track = static_cast<Candidate*>(fItVertexingTrackInputArray->Next())))   {
//cout<<"track"<<endl;
    //if(fabs(track->Position.T() - pv_t) < 30.)  {
    if(track->VertexID_gen == pvID)  {
    
        pvTrack = static_cast<Candidate*>(track->Clone());
         //fDebugOutputCollector.fillContainer ("m_partType", pvTrack->VertexID_gen) ;
        fPvEFlowTrackOutputArray->Add(pvTrack);
        }
    }
//fItPvEFlowTrackOutputArray = fPvEFlowTrackOutputArray->MakeIterator();
//fItPvEFlowTrackOutputArray->Reset();    
  //  while((track = static_cast<Candidate*>(fItPvEFlowTrackOutputArray->Next())))    {
          // fDebugOutputCollector.fillContainer ("m_partType", track->VertexID_gen) ;
      //      }
fCalTrackArray->Clear();
fRecoVertexArray->Clear();
fPvEFlowTowerCandidateArray->Clear();
++k;
}               
                
                
                
                
                
                
              
//------------------------------------------------------------------------------

Double_t ECalorimeter::LogNormal(Double_t mean, Double_t sigma)
{
  Double_t a, b;

  if(mean > 0.0)
  {
    b = TMath::Sqrt(TMath::Log((1.0 + (sigma*sigma)/(mean*mean))));
    a = TMath::Log(mean) - 0.5*b*b;

    return TMath::Exp(a + b*gRandom->Gaus(0, 1));
  }
  else
  {
    return 0.0;
  }
}
