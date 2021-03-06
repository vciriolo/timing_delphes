#ifndef HCalorimeter_h
#define HCalorimeter_h

/** \class Calorimeter
 *  Fills calorimeter towers, performs calorimeter resolution smearing,
 *  preselects towers hit by photons and creates energy flow objects.
 *  $Date: 2013-09-03 17:54:56 +0200 (Tue, 03 Sep 2013) $
 *  \author P. Demin - UCL, Louvain-la-Neuve
 */

#include "classes/DelphesModule.h"
#include "modules/simpleVariableCollector.h"

#include <map>
#include <set>
#include <vector>

#include "TF1.h"

class TObjArray;
class DelphesFormula;
class Candidate;

class HCalorimeter: public DelphesModule {
public:

  HCalorimeter();
  ~HCalorimeter();

  void Init();
  void Process();
  //void PileUpSubtractor();
  void Finish();

private:

  //typedef std::map< Long64_t, std::pair< Double_t, Double_t > > TFractionMap; //!
  typedef std::map< Long64_t, Double_t > TFractionMap; //!
  typedef std::map< Double_t, std::set< Double_t > > TBinMap; //!

  Candidate *fTower; // candidate
  Double_t fTowerEta, fTowerPhi, fTowerEdges[4];
  Double_t fTowerECalEnergy, fTowerHCalEnergy;
  Double_t fTrackECalEnergy, fTrackHCalEnergy;
  Int_t    fTowerTrackHits, fTowerPhotonHits;

  Double_t fTimingEMin;

  TFractionMap fFractionMap; //!
  
  TBinMap      fBinMap; //!

  std::vector < Double_t > fEtaBins;
  std::vector < std::vector < Double_t >* > fPhiBins;

  std::vector < Long64_t > fTowerHits;

  std::vector < Double_t > fTowerECalFractions;
  std::vector < Double_t > fTowerHCalFractions;

  std::vector < Double_t > fTrackECalFractions;
  std::vector < Double_t > fTrackHCalFractions;

  std::vector<int> fParticlePDGId ;
  std::vector<int> fTrackPDGId ;


  DelphesFormula *fECalResolutionFormula; //!
  DelphesFormula *fHCalResolutionFormula; //!

  TIterator *fItParticleInputArray; //!
  TIterator *fItTrackInputArray; //!
  TIterator *fItVertexInputArray;
  TIterator *fItVertexingTrackInputArray;

  const TObjArray *fParticleInputArray; //!
  const TObjArray *fTrackInputArray; //!
  const TObjArray *fVertexInputArray;
  const TObjArray *fVertexingTrackInputArray;

  // HCAL+ECAL towers with all the energy accounted
  TObjArray *fTowerOutputArray; //!
  TIterator *fItTowerOutputArray;
  // photons
  TObjArray *fPhotonOutputArray; //!

  // tracks, taken in input and replicated in output
  TObjArray *fEFlowTrackOutputArray; //!
  TIterator *fItEFlowTrackOutputArray;
  // PF towers: energy deposits cleaned from tracks
  TObjArray *fEFlowTowerOutputArray; //!
  TIterator *fItEFlowTowerOutputArray;

  TObjArray *fTowerTrackArray; //!
  TIterator *fItTowerTrackArray; //!

  TObjArray *fVertexArray;
  TIterator *fItVertexArray;
  
  TObjArray *fTimeVertexArray;
  TIterator *fItTimeVertexArray;
  
  TObjArray *fCalTrackArray;
  TIterator *fItCalTrackArray;

  void FinalizeTower();
  Double_t LogNormal(Double_t mean, Double_t sigma);

  void PileUpSubtractor();
  bool fChargedFromTrack; // for timing

  const TObjArray *fLHEPartonInputArray; //! 
  TIterator *fItLHEPartonInputArray; //!  

  TF1 * fDelayBarrel ;
  TF1 * fDelayEndcap ;

  simpleVariableCollector fDebugOutputCollector ;    
  int fEventCounter ;                                                                                                                                             

  ClassDef(HCalorimeter, 1)
};

#endif
