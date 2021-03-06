/** \class Calorimeter
 *  Fills calorimeter towers, performs calorimeter resolution smearing,
 *  preselects towers hit by photons and creates energy flow objects.
 *  \author P. Demin - UCL, Louvain-la-Neuve
*/

#include "modules/HCalorimeter.h"

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

HCalorimeter::HCalorimeter() :
  fECalResolutionFormula(0), fHCalResolutionFormula(0),
  fItParticleInputArray(0), fItTrackInputArray(0),
  fTowerTrackArray(0), fItTowerTrackArray(0),
  fItLHEPartonInputArray (0), fItVertexInputArray(0)
  {
  
  //fDelayBarrel and fDelayEndcap by Pietro
  //fDelayBarrel = new TF1 ("fDelayBarrel", "pol5", 0, 5) ;
  //fDelayEndcap = new TF1 ("fDelayEndcap", "pol5", 0, 5) ;
  
    
  fDelayBarrel = new TF2 ("fDelayBarrel", "x/y", 0, 100, 0, 100);    
  fDelayEndcap = new TF2 ("fDelayEndcap", "x/y", 0, 100, 0, 100);
  
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
}

//------------------------------------------------------------------------------

HCalorimeter::~HCalorimeter(){

  delete fDelayBarrel ;
  delete fDelayEndcap ;
  if(fECalResolutionFormula)        delete fECalResolutionFormula;
  if(fHCalResolutionFormula)        delete fHCalResolutionFormula;
  if(fTowerTrackArray)              delete fTowerTrackArray;
  if(fItTowerTrackArray)            delete fItTowerTrackArray;
  if(fItLHEPartonInputArray)        delete fItLHEPartonInputArray;
  if(fItVertexInputArray)           delete fItVertexInputArray;
  if(fTimeVertexArray)              delete fTimeVertexArray;
  if(fItVertexingTrackInputArray)   delete fItVertexingTrackInputArray;
  if(fCalTrackArray)                delete fCalTrackArray;
}

//------------------------------------------------------------------------------

void HCalorimeter::Init(){

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
    hcalFraction = paramFractions[0].GetDouble();
    fFractionMap[param[i*2].GetInt()] = hcalFraction;
    //hcalFraction = paramFractions[1].GetDouble();
    //fFractionMap[param[i*2].GetInt()] = make_pair(ecalFraction, hcalFraction);
  }

  // read resolution formulas
  //fECalResolutionFormula->Compile(GetString("ECalResolutionFormula", "0"));
  fHCalResolutionFormula->Compile(GetString("HCalResolutionFormula", "0"));

  // import array with output from other modules
  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "ParticlePropagator/particles"));
  fItParticleInputArray = fParticleInputArray->MakeIterator();

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "ParticlePropagator/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  //PG for matching to LHE particles, for timing studies
  fLHEPartonInputArray = ImportArray(GetString("LHEPartonInputArray", "Delphes/LHEParticles"));
  fItLHEPartonInputArray = fLHEPartonInputArray->MakeIterator();
    
  // import smeared vertex array 
  fVertexInputArray = ImportArray(GetString("VertexInputArray", "PileUpMerger/vertices"));
  fItVertexInputArray = fVertexInputArray->MakeIterator();
  
  //import tracks from vertexing
  fVertexingTrackInputArray = ImportArray(GetString("VertexingTrackInputArray", "Vertexing/vertexingtracks"));
  fItVertexingTrackInputArray = fVertexingTrackInputArray->MakeIterator();
  

  // create output arrays
  // calo tower output
  fTowerOutputArray  = ExportArray(GetString("TowerOutputArray", "hcalTowers"));
  fItTowerOutputArray = fTowerOutputArray->MakeIterator();
  // photons 
  //fPhotonOutputArray = ExportArray(GetString("neutralHadronsOutputArray", "neutralHadrons"));
  fPhotonOutputArray = ExportArray(GetString("photonOutputArray", "photons"));
  // track for particle flow
  fEFlowTrackOutputArray = ExportArray(GetString("EFlowTrackOutputArray", "hcalEflowTracks"));
  fItEFlowTrackOutputArray = fEFlowTrackOutputArray->MakeIterator();
  // tower for particle flow
  fEFlowTowerOutputArray = ExportArray(GetString("EFlowTowerOutputArray", "hcalEflowTowers"));

  // For timing
  // So far this flag needs to be false
  // since the curved extrapolation not supported
  // if this value is true, electron timing is taken from the track,
  // otherwise is taken from the particle collection
  fChargedFromTrack  = false;
  
  //PG FIXME where does this come from?
  // suggested from A. Bornheim, reasonable according to him    
  
  fTimingEMin = GetDouble ("TimingEMin", 0.) ;

  //simple outputs during running
//   fDebugOutputCollector.addVariable ("DR") ;
//   fDebugOutputCollector.addVariable ("m_partType") ;
//   fDebugOutputCollector.addVariable ("Nm_partType") ;
  fDebugOutputCollector.addVariable ("PUindex") ;
  fDebugOutputCollector.addVariable3D ("m_eta:time:pt") ;
  fDebugOutputCollector.addVariable3D ("m_eta:time_nc:pt") ;
  fDebugOutputCollector.addVariable3D ("Nm_eta:time:pt") ;
  fDebugOutputCollector.addVariable4D ("eta:time:pt:PID") ;

  fDebugOutputCollector.addVariable3D ("tower_eta:time:pt") ;
  fDebugOutputCollector.addVariable3D ("tower_eta:time:E") ;
  fEventCounter = 0 ;

}

//------------------------------------------------------------------------------

void HCalorimeter::Finish(){

  std::string outfile = GetString ("simpleOutputFileName", "simpleOutput_HCa.root") ;
  fDebugOutputCollector.save (outfile) ;

  if(fItParticleInputArray)  delete fItParticleInputArray;
  if(fItTrackInputArray)     delete fItTrackInputArray;
  if(fItVertexInputArray)    delete fItVertexInputArray;
  if(fItVertexingTrackInputArray) delete fItVertexingTrackInputArray;

  vector< vector< Double_t >* >::iterator itPhiBin;
  for (itPhiBin = fPhiBins.begin(); itPhiBin != fPhiBins.end(); ++itPhiBin){
    delete *itPhiBin;
  }
}

//------------------------------------------------------------------------------

void HCalorimeter::Process(){

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

  // get the eta, phi, pt of the LHE partons that could generate
  // a Delphes jet
  vector<TLorentzVector> LHEParticles ;
  fItLHEPartonInputArray->Reset () ;
  
  while (Candidate * LHEparticle = static_cast<Candidate*> (fItLHEPartonInputArray->Next ()))
  
    { 
      if (LHEparticle->Status != 1) continue ;
      if (fabs (LHEparticle->PID) != 11 ||   // electron
//          1) 
          //fabs (LHEparticle->PID) < 7   ||   // quarks
          fabs (LHEparticle->PID) != 21 ||   // gluon
          fabs (LHEparticle->PID) != 22)     // photon
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

    hcalFraction = itFractionMap->second;
    fTrackHCalFractions.push_back(hcalFraction);
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
  while((particle = static_cast<Candidate*>(fItParticleInputArray->Next()))){
  

    const TLorentzVector &particlePosition = particle->Position;
    ++number;
    pdgCode = TMath::Abs(particle->PID);
    itFractionMap = fFractionMap.find(pdgCode); // find the particle in the fraction map

    if(itFractionMap == fFractionMap.end()){
      itFractionMap = fFractionMap.find(0);
    }
    
    //fill ecal tower fraction vectors
    hcalFraction = itFractionMap->second;
    fTowerHCalFractions.push_back(hcalFraction);
    fParticlePDGId.push_back(pdgCode);

    //ecalFraction = itFractionMap->second.first; // take ecal fraction
    //hcalFraction = itFractionMap->second.second; // take Hcal fraction
  
    // fill tower fraction vectors
    /*fTowerECalFractions.push_back(ecalFraction);
    fTowerHCalFractions.push_back(hcalFraction);
    fParticlePDGId.push_back(pdgCode);  */

    //if(ecalFraction < 1.0E-9 && hcalFraction < 1.0E-9) continue;
    if( hcalFraction < 1.0E-9) continue;

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

      //fTowerECalEnergy = 0.0;
      fTowerHCalEnergy = 0.0;

      //fTrackECalEnergy = 0.0;
      fTrackHCalEnergy = 0.0;

      fTowerTrackHits = 0;
      fTowerPhotonHits = 0;

      fTowerTrackArray->Clear();
    }  // first time hit

    // check for track hits
    if(flags & 1){

      ++fTowerTrackHits;
      track = static_cast<Candidate*>(fTrackInputArray->At(number));

      //track = static_cast<Candidate*>(fVertexingTrackInputArray->At(number));

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

      //ecalEnergy = momentum.E() * fTrackECalFractions[number];
      hcalEnergy = momentum.E() * fTrackHCalFractions[number];

      //fTrackECalEnergy += ecalEnergy;
      fTrackHCalEnergy += hcalEnergy;

      //totalEnergy = ecalEnergy + hcalEnergy ;
      totalEnergy = hcalEnergy;
      
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
     
      if ( totalEnergy > fTimingEMin && fTower && fChargedFromTrack) {
        
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
            std::make_pair<float,float>(totalEnergy, time - TOF));
         // std::make_pair<float,float>(totalEnergy, (double)time));
      }
            
      fTowerTrackArray->Add(track);

      continue; // go to the next hit
    } // check for track hits

    // check for photon and electron hits in current tower
    if(flags & 2) ++fTowerPhotonHits;

    particle = static_cast<Candidate*>(fParticleInputArray->At(number));
    momentum = particle->Momentum;

    // fill current tower
    hcalEnergy = momentum.E() * fTowerHCalFractions[number];
    //hcalEnergy = momentum.E() * fTowerHCalFractions[number];

    //fTowerECalEnergy += ecalEnergy;
    fTowerHCalEnergy += hcalEnergy;

    //totalEnergy = ecalEnergy + hcalEnergy ;
    totalEnergy = hcalEnergy;

    if ( (totalEnergy > fTimingEMin && fTower) && !fChargedFromTrack)   {
         //(abs(particle->PID) != 11 ||  !fChargedFromTrack) ) {

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
                       
        //double distance = TMath::Sqrt(TMath::Power(ParticlePositionZ,2)+TMath::Power(ParticlePositionX,2)+TMath::Power(ParticlePositionY,2));;
        //double velocity = (TMath::Sqrt(TMath::Power(ParticleMomentumZ,2)+TMath::Power(ParticleMomentumX,2)+TMath::Power(ParticleMomentumY,2))/ totalEnergy) ;//* c_light ;
         
         double distance = TMath::Sqrt(TMath::Power(ParticlePositionX,2)+TMath::Power(ParticlePositionY,2));;
         double velocity = (TMath::Sqrt(TMath::Power(particle->Momentum.Z(),2)+TMath::Power(particle->Momentum.X(),2)+TMath::Power(particle->Momentum.Y(),2))/track->Momentum.E()) ;
                   
                   float Pt = particle->Momentum.Pt();
                   float R = (Pt/(0.3*3.8));//Bz; m
                   //secant with respect to 0.0.0
                   float s = (TMath::Sqrt(TMath::Power(particle->Position.X(),2)+TMath::Power(particle->Position.Y(),2)));
                   float alpha = TMath::ASin(s/(1000*2*R));
                   float arc = 2 * alpha * R * 1000;
                   
                   float z_distance = track->Position.Z();
                   TOF = (TMath::Sqrt((TMath::Power(particle->Position.Z(),2) + TMath::Power(arc,2)))) *inv_c_light;       
        
      /*  if (feta < 1.6) TOF = fDelayBarrel->Eval(distance,velocity);
        
        else                  TOF = fDelayEndcap->Eval(distance,velocity);  */
        
        
        //cout<<time<<"\t"<<TOF<<endl;
        // the timing is corrected in the map,
        // since I assume that the knowledge of the position of the detID
        // is enough to determine the correction.
        // The PV z position could be used to introduce a correction
        // to the delay of the order of cz
        fTower->ecal_E_t.push_back(
          //std::make_pair<float,float> (totalEnergy, time - delay);
          //std::make_pair<float,float> (totalEnergy, (double)time));     
            std::make_pair<float,float> (totalEnergy, time - TOF));
//         double Rmin = 100. ; 
//         for (unsigned int i = 0 ; i < LHEParticles.size () ; ++i)
//           {
//             if (particle->Position.DeltaR (LHEParticles.at (i)) < Rmin) 
//               Rmin = particle->Position.DeltaR (LHEParticles.at (i)) ;
//           }
//         fDebugOutputCollector.fillContainer ("DR", Rmin) ;

       
        
        
        

       // fDebugOutputCollector.fillContainer ("PUindex", particle->IsPU) ;




//}
           /* if (isMatching (particle->Position, LHEParticles))  

         {               
               
            //fDebugOutputCollector.fillContainer ("m_partType", particle->PID) ;
            fDebugOutputCollector.fillContainer3D ("m_eta:time:pt", 
            feta, (time - TOF) , particle->Momentum.Pt ()) ; //feta
            fDebugOutputCollector.fillContainer3D ("m_eta:time_nc:pt", 
            feta, time * inv_c_light, particle->Momentum.Pt ()) ;    //feta
            
          }
            
            
          //}
          else
            {
            fDebugOutputCollector.fillContainer ("Nm_partType", particle->PID) ;
            fDebugOutputCollector.fillContainer3D ("Nm_eta:time:pt", 
               feta, (time - TOF) , particle->Momentum.Pt ()) ;
          }
   
*/
    }

    fTower->PID = fParticlePDGId.at(number);
    fTower->AddCandidate(particle);
              //cout<<"hcal"<<endl;
  
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
  //PileUpSubtractor();
  ++fEventCounter ;

}


//------------------------------------------------------------------------------

void HCalorimeter::FinalizeTower(){

  Candidate *track, *tower;
  Double_t energy, pt, eta, phi;
  Double_t ecalEnergy, hcalEnergy;
  Double_t ecalSigma, hcalSigma;

  if(!fTower) return;

  // ECAL resolution
 /* ecalSigma  = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerECalEnergy);
  ecalEnergy = LogNormal(fTowerECalEnergy, ecalSigma);  */

  // HCAL resolution    

  hcalSigma  = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerHCalEnergy);
  hcalEnergy = LogNormal(fTowerHCalEnergy, hcalSigma);  
    
  //energy     = ecalEnergy + hcalEnergy;
  energy = hcalEnergy;

  eta = fTowerEta;
  phi = fTowerPhi;

  pt = energy / TMath::CosH(eta);

  // Time calculation for tower
  fTower->nTimes = 0;
  float tow_sumT = 0;
  float tow_sumW = 0;

  for (unsigned int i = 0 ; i < fTower->ecal_E_t.size() ; i++) { 
    float w = TMath::Sqrt(fTower->ecal_E_t[i].first);
    //cout<<fTower->ecal_E_t[i].first<<"\t"<<fTower->ecal_E_t[i].second<<endl;
    tow_sumT += w*fTower->ecal_E_t[i].second;
    tow_sumW += w;
    
    fTower->nTimes++;
  }
  
  if (tow_sumW > 0.) {
    fTower->Position.SetPtEtaPhiE(1.0, eta, phi,tow_sumT/tow_sumW);
    //fTower->Position.SetT(tow_sumT/tow_sumW);
    
   /* fDebugOutputCollector.fillContainer3D ("tower_eta:time:E",
     fabs (fTower->Position.Eta ()),
     (fTower->Position.T ()) ,
     fTower->Momentum.E ()
    ) ; */
    
  } else {
    fTower->Position.SetPtEtaPhiE(1.0,eta,phi,999999.);
  }

  //  fTower->Position.SetPtEtaPhiE(1.0, eta, phi, 0.);
  fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
  //fTower->Eem = ecalEnergy;
  fTower->Ehad = hcalEnergy;


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
    
   /* fDebugOutputCollector.fillContainer3D ("tower_eta:time:E",
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

  hcalEnergy -= fTrackHCalEnergy;
  //if(ecalEnergy < 0.0) ecalEnergy = 0.0;
  
 // hcalEnergy -= fTrackHCalEnergy;
  if(hcalEnergy < 0.0) hcalEnergy = 0.0;

//  energy = ecalEnergy + hcalEnergy;

    energy = hcalEnergy;
  // save ECAL and/or HCAL energy excess as an energy flow tower
  if(energy > 0.0){

    // create new tower
    tower = static_cast<Candidate*>(fTower->Clone());
    pt = energy / TMath::CosH(eta);
    tower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
    //tower->Eem = ecalEnergy;
    tower->Ehad = hcalEnergy;
    fEFlowTowerOutputArray->Add(tower);

    
  }
  

  

}
//------------------------------------------------------------------------------
void HCalorimeter::PileUpSubtractor()    {
Candidate *vertex, *tower, *track;
double_t vertex_time, tower_time;
float_t vertex_energy, tmp_t, TOF;

//link tracks to towers
//fItTowerTrackArray->Reset();

fItVertexingTrackInputArray->Reset();
while((track = static_cast<Candidate*>(fItVertexingTrackInputArray->Next())))    {
//fItEFlowTrackOutputArray->Reset();
//while((track = static_cast<Candidate*>(fItEFlowTrackOutputArray->Next())))  {
//cout<<"pu"<<endl;
    float track_eta = track->Position.Eta();
    float track_phi = track->Position.Phi();
       
       
       
     /*  float Pt = track->Momentum.Pt();
                   float R = (Pt/(0.3*3.8));//Bz; m
                   //secant with respect to 0.0.0
                   float s = (TMath::Sqrt(TMath::Power(track->Position.X(),2)+TMath::Power(track->Position.Y(),2)));
                   float alpha = TMath::ASin(s/(1000*2*R));
                   float arc = 2 * alpha * R * 1000;
                   
                   float z_distance = track->Position.Z();
                   TOF = (TMath::Sqrt((TMath::Power(track->Position.Z(),2) + TMath::Power(arc,2)))) *inv_c_light;
        
        float track_T = track->Position.T() - TOF;
        float time = track->Position.T () * inv_c_light ;
       track->Position.SetT(track_T);
      */ 
        
    /* fDebugOutputCollector.fillContainer3D ("tower_eta:time:pt",
     track_eta,
     track_phi,
     track->Position.T()
    ) ; */
    
    fItTowerOutputArray->Reset();
    while((tower = static_cast<Candidate*>(fItTowerOutputArray->Next()))){
    
        
       
        
    
        //float tower_eta = tower->Position.Eta();
        //float tower_phi = tower->Position.Phi();
        if(track_eta < tower->Edges[1] && track_eta > tower->Edges[0] && track_phi < tower->Edges[3] && track_phi > tower->Edges[2])    {
            track->Position.SetT(tower->Position.T());
            track->Momentum.SetE(tower->Ehad);
            
             fCalTrackArray->Add(track);
             
   /*           fDebugOutputCollector.fillContainer3D ("tower_eta:time:pt",
     track->Momentum.E(),
     tower->Position.T(),
     track->Position.T()
    ) ;    */
                            
            //cout<<track->Position.T()*inv_c_light<<endl;
            }
            
        /*    fDebugOutputCollector.fillContainer3D ("tower_eta:time:pt",
     track->Position.Phi(),
     tower->Position.T(),
     track->Position.T()
    ) ; */
            
        }
        
      /*  
        fDebugOutputCollector.fillContainer3D ("tower_eta:time:E",
      fabs(fTower->Position.Eta()),
     (fTower->Position.T ()) * inv_c_light,
     fTower->Momentum.E ()
    ) ;  */
}
              
//set vertices time and energy
int n = 1;
vertex_time = 0;
float sum_track_time = 0;
float track_time = 0;
vertex_energy = 0;
float track_energy = 0;

fItVertexInputArray->Reset();
while((vertex = static_cast<Candidate*>(fItVertexInputArray->Next())))  {
    int vertexID = vertex->VertexID_gen;

fItCalTrackArray->Reset();
while((track = static_cast<Candidate*>(fItCalTrackArray->Next())))  {

        int trackID = track->VertexID_gen;       
        if(trackID == vertexID) {
 
            track_time = track->Position.T();
            track_energy = track->Momentum.E();
            if(track_time > 600.) {
 
                track_time = 0.;
                n = n - 1.;
                track_energy = 0.;
                
                }   
                sum_track_time += track_time;
                vertex_energy += track_energy;
            ++n;
            }
  //  }
    } 
    vertex_time = sum_track_time/n;
    //cout<<vertex_time<<"\t"<<vertex_energy<<"\t"<<vertex->IsPU<<"\t"<<vertex->VertexID_gen<<endl;
    vertex->Position.SetT(vertex_time);
    vertex->Momentum.SetE(vertex_energy);
  //  if(vertex->VertexID_gen == 0 && vertex->IsPU == 1 ) goto LABEL;
        //vertex_time = 99999.;
        //vertex_energy = 999999;
        //}
 /*   fDebugOutputCollector.fillContainer3D ("tower_eta:time:E",
        vertex->Position.T(),
        vertex->Momentum.E(),
        vertex->IsPU
        ) ;  */
    //    cout<<vertex_time<<"\t"<<vertex_energy<<"\t"<<vertex->IsPU<<"\t"<<vertex->VertexID_gen<<endl;
       
    //LABEL:n = 1;
    sum_track_time = 0;
    vertex_energy = 0;
    //}
    //cout<<vertex->Position.T()<<endl;

    
            
}






//chose primary vertex 
/*
fItVertexArray = fVertexArray->MakeIterator();
//fItEFlowTowerOutputArray = fEFlowTowerOutputArray->MakeIterator();
TIterator *fItTowerOutputArray = fEFlowTowerOutputArray->MakeIterator();

float tmp_energy = 0;
fItTimeVertexArray->Reset();
while((vertex = static_cast<Candidate*>(fItTimeVertexArray->Next())))   {
    vertex_time = vertex->Position.T();
    //vertex->Momentum.SetE(0.);
    //fItEFlowTowerOutputArray->Reset();
    fItTowerOutputArray->Reset();
   // while((tower = static_cast<Candidate*>(fItEFlowTowerOutputArray->Next())))   {
    while((tower = static_cast<Candidate*>(fItTowerOutputArray->Next())))   {
        tower_time = fTower->Position.T()*inv_c_light;
        //cout<<vertex_time<<"\t"<<tower_time<<endl;
        /*if(vertex_time - tower_time > 30.)  {
        tower->Eem = 0.;
        }
        else vertex_energy += tower->Eem;
        if(fabs(vertex_time - tower_time)<100.)  {
            vertex_energy += tower->Eem;
            //cout<<vertex_energy<<endl;
            }
        }
        vertex->Momentum.SetE(vertex_energy);
       /* fDebugOutputCollector.fillContainer3D ("tower_eta:time:pt",
        vertex_time,
        vertex_energy,
        vertex->VertexID_gen
        ) ; 
        }  
//tmp_t = 0;
/*
fItTimeVertexArray->Reset();

while((vertex = static_cast<Candidate*>(fItTimeVertexArray->Next())))   {
        float v_energy = vertex->Momentum.E();
        //cout<<v_energy<<endl;
    //if(tmp_energy < vertex->Momentum.E())   {
    if(tmp_energy < vertex_energy)   {
    //cout<<vertex_energy<<endl;
        //tmp_energy = vertex->Momentum.E();
        tmp_energy = vertex_energy;
        tmp_t = vertex->Position.T();
       // v_energy = tmp_energy;
       // v_t = tmp_t;
        }
    //cout<<tmp_t<<"\t"<<tmp_energy<<endl;
    
}   //cout<<tmp_t<</*"\t"<<tmp_energy<<endl;
   // fDebugOutputCollector.fillContainer ("PUindex", tmp_t) ;
//    cout<<tmp_t<<endl;
/*fDebugOutputCollector.fillContainer3D ("tower_eta:time:pt",
        tmp_energy,
        tmp_t,
        tmp_energy
        );     
    
    //cout<<tmp_energy<<"\t"<<tmp_t<<endl;
    
//set to 0 energy of out of time towers
    //fItEFlowTowerOutputArray->Reset();
    fItTowerOutputArray->Reset();
    int j = 0;
    //while((tower = static_cast<Candidate*>(fItEFlowTowerOutputArray->Next())))  {
    while((tower = static_cast<Candidate*>(fItTowerOutputArray->Next())))  {
        //cout<<tower->Position.T()*inv_c_light<<"\t"<<tmp_t<<endl;
        //if(fabs(tower->Position.Eta())>1.6)   {
        //cout<<tower->Position.T()<<endl;
        //}
       
        if(fabs(fabs(tower->Position.T() * inv_c_light) - fabs(tmp_t)) <30.)  {
        ++j;
        //tower->Eem = 0.;
        //}   else { 
        //cout<<tower->Eem<<endl;
        fDebugOutputCollector.fillContainer3D ("tower_eta:time:pt",
        tower->Eem,
        tower->Position.T() * inv_c_light,
        tmp_t       
        ) ;
        } 
        /*if((fabs(tower->Position.T() * inv_c_light) - fabs(tmp_t)) > 30.)   {
        tower->Momentum.SetE(999999.);
        tower->Position.SetT(999999.);
        }
        fDebugOutputCollector.fillContainer3D ("tower_eta:time:pt",
        tower->Eem,
        tower->Position.T() * inv_c_light,
        tmp_t       
        ) ;
    //cout<<tower->Eem<<endl;
    }
    //cout<<fTowerOutputArray->GetEntriesFast()<<"\t"<<j<<endl;*/
//}

}

//------------------------------------------------------------------------------

Double_t HCalorimeter::LogNormal(Double_t mean, Double_t sigma)
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
