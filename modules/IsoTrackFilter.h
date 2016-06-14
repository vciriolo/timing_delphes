#ifndef IsoTrackFilter_h
#define IsoTrackFilter_h

/** \class IsoTrackFilter
 *
 *  Sums transverse momenta of IsoTrack objects (tracks) within a DeltaR 
 *  cone around a candidate and calculates fraction of this sum
 *  to the candidate's transverse momentum. outputs candidates that have
 *  the transverse momenta fraction within (PTRatioMin, PTRatioMax].
 *
 *  $Date: 2013-11-04 13:14:33 +0100 (Mon, 04 Nov 2013) $
 *  $Revision: 1317 $
 *
 *
 *  \author Z. Wu - Baylor
 *
 */

#include "classes/DelphesModule.h"

class TObjArray;

class ExRootFilter;
class IsoTrackClassifier;

class IsoTrackFilter: public DelphesModule
{
public:

  IsoTrackFilter();
  ~IsoTrackFilter();

  void Init();
  void Process();
  void Finish();
  void IsoTrackSelector(std::vector< TIterator * >& fAllInputList, ExRootFilter* fCandidateFilter, bool IsEM);

private:

  Double_t fDeltaRMax;

  Double_t fPTRatioMax;

  Double_t fPTSumMax;
  
  Double_t fIsoTrackPTMin;

  Bool_t fUsePTSum;



  ExRootFilter *fEleFilter;
  ExRootFilter *fMuonFilter;
  ExRootFilter *fHADFilter;


  TIterator *fItEleInputArray; //!
  TIterator *fItMuonInputArray; //!
  TIterator *fItHADInputArray; //!

  const TObjArray *fEleInputArray; //!
  const TObjArray *fMuonInputArray; //!
  const TObjArray *fHADInputArray; //!

  // Isolation using all input list
  std::vector< TIterator * > fAllInputList; //!

  IsoTrackClassifier *fClassifier; //!

  /// Output Array
  TObjArray *fOutputArray; //!

  ClassDef(IsoTrackFilter, 1)
};

#endif
