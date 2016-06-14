#ifndef FastJetFinder_h
#define FastJetFinder_h

/** \class FastJetFinder
 *
 *  Finds jets using FastJet library.
 *  $Date: 2013-11-04 11:59:27 +0100 (Mon, 04 Nov 2013) $
 *  $Revision: 1315 $
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"
#include "modules/simpleVariableCollector.h"

#include <map>

class TObjArray;
class TIterator;

namespace fastjet {
  class JetDefinition;
  class AreaDefinition;
  class Selector;
}

class FastJetFinder: public DelphesModule {

public:

  FastJetFinder();
  ~FastJetFinder();

  void Init();
  void Process();
  void Finish();

private:

  void *fPlugin; //!
  fastjet::JetDefinition *fDefinition; //!

  // For genjets mostly
  Int_t fKeepPileUp;
  Int_t fJetAlgorithm;
  Int_t fMaxIterations;
  Int_t fMaxPairSize;
  Int_t fIratch;

  Double_t fParameterR;
  Double_t fJetPTMin;
  Double_t fConeRadius;
  Double_t fSeedThreshold;
  Double_t fConeAreaFraction;
  Double_t fAdjacencyCut;
  Double_t fOverlapThreshold;

  // --- FastJet Area method --------
  fastjet::AreaDefinition *fAreaDefinition;
  Int_t  fAreaAlgorithm;
  Bool_t fComputeRho;
  Bool_t fComputeRhoGrid;
  Bool_t fComputeRhoGridParticles;

  // -- ghost based areas --
  Double_t fGhostEtaMax;
  Int_t    fRepeat;
  Double_t fGhostArea;
  Double_t fGridScatter;
  Double_t fPtScatter;
  Double_t fMeanGhostPt;

  // -- voronoi areas --
  Double_t fEffectiveRfact;

  Double_t fParticlePTMin;

  std::map< Double_t, Double_t > fEtaRangeMap; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!
  TObjArray *fRhoOutputArray; //!

  simpleVariableCollector fDebugOutputCollector ;    
  int fEventCounter ;                                                                                                                                             

  ClassDef(FastJetFinder, 1)

};

#endif
