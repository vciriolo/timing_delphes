
/** \class IsoTrackFilter
 *
 *  Sums transverse momenta of Track objects within a DeltaR cone around a
 *  candidate and calculates fraction of this sum to the candidate's
 *  transverse momentum. outputs candidates that have the transverse momenta
 *  fraction within (PTRatioMin, PTRatioMax].
 *
 *  $Date: 2013-11-04 13:14:33 +0100 (Mon, 04 Nov 2013) $
 *  $Revision: 1317 $
 *
 *
 *  \author  Z. Wu - Baylor
 *
 */

#include "modules/IsoTrackFilter.h"

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

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

class IsoTrackClassifier : public ExRootClassifier
{
public:

  IsoTrackClassifier() {}

  Int_t GetCategory(TObject *object);

  Double_t fPTMin;

};

//------------------------------------------------------------------------------

Int_t IsoTrackClassifier::GetCategory(TObject *object)
{
  Candidate *track = static_cast<Candidate*>(object);
  const TLorentzVector &momentum = track->Momentum;

  if(momentum.Pt() < fPTMin) return -1;

  return 0;
}

//------------------------------------------------------------------------------

IsoTrackFilter::IsoTrackFilter():
  fEleFilter(0),fMuonFilter(0),fHADFilter(0),
  fItEleInputArray(0), fItMuonInputArray(0), fItHADInputArray(0),
  fEleInputArray(0), fMuonInputArray(0), fHADInputArray(0),
  fClassifier(0), fOutputArray(0)

{
  fClassifier = new IsoTrackClassifier;
}

//------------------------------------------------------------------------------

IsoTrackFilter::~IsoTrackFilter()
{
}

//------------------------------------------------------------------------------

void IsoTrackFilter::Init()
{

  fDeltaRMax = GetDouble("DeltaRMax", 0.5);

  fPTRatioMax = GetDouble("PTRatioMax", 0.1);

  fPTSumMax = GetDouble("PTSumMax", 5.0);

  fUsePTSum = GetBool("UsePTSum", false);


  fIsoTrackPTMin = GetDouble("IsoTrackPTMin", 10);

  fClassifier->fPTMin = GetDouble("PTMin", 0.5);

  // import input array(s)
  
  // Electron Input Array
  fEleInputArray = ImportArray(GetString("ElectronInputArray", "ElectronEnergySmearing/electrons"));
  fEleFilter = new ExRootFilter(fEleInputArray);

  // Muon Input Array
  fMuonInputArray = ImportArray(GetString("MuonInputArray", "MuonMomentumSmearing/muons"));
  fMuonFilter = new ExRootFilter(fMuonInputArray);

  // Charged Hadron Input Array
  fHADInputArray = ImportArray(GetString("HADInputArray", "ChargedHadronMomentumSmearing/chargedHadrons"));
  fHADFilter = new ExRootFilter(fHADInputArray);

  // All tracks used for isolation calculation
  fAllInputList.push_back(fEleInputArray->MakeIterator());
  fAllInputList.push_back(fMuonInputArray->MakeIterator());
  fAllInputList.push_back(fHADInputArray->MakeIterator());

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "IsoTrack"));
}

//------------------------------------------------------------------------------

void IsoTrackFilter::Finish()
{
  // delete Filter
  if(fEleFilter) delete fEleFilter;
  if(fMuonFilter) delete fMuonFilter;
  if(fHADFilter) delete fHADFilter;

  if (fClassifier) delete fClassifier;
}

//------------------------------------------------------------------------------

void IsoTrackFilter::Process()
{
  IsoTrackSelector(fAllInputList, fEleFilter, true);
  IsoTrackSelector(fAllInputList, fMuonFilter, true);
  IsoTrackSelector(fAllInputList, fHADFilter, false);
}

//------------------------------------------------------------------------------
//
void IsoTrackFilter::IsoTrackSelector(std::vector< TIterator * >& fAllInputList, ExRootFilter*  fCandidateFilter, bool IsEM)
{
  Candidate *candidate, *IsoTrack;
  TObjArray *CanTrackArray;

  // Prepare the candidate array
  fCandidateFilter->Reset();
  CanTrackArray = fCandidateFilter->GetSubArray(fClassifier, 0);

  if(fAllInputList.size() == 0) return;

  Double_t sum, ratio;
  Int_t counter;

  // loop over all input candidates
  TIter itCanTrackArray(CanTrackArray);
  itCanTrackArray.Reset();
  while((candidate = static_cast<Candidate*>(itCanTrackArray.Next())))
  {
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    sum = 0.0;
    counter = 0;

    // loop over all input tracks
    for(std::vector< TIterator * >::iterator itInputList = fAllInputList.begin(); itInputList != fAllInputList.end(); ++itInputList)
    {
      TIterator *iterator = *itInputList;
      // loop over all candidates
      iterator->Reset();

      while((IsoTrack = static_cast<Candidate*>(iterator->Next())))
      {
        if (fClassifier->GetCategory(IsoTrack) == -1) continue;

        const TLorentzVector &IsoTrackMomentum = IsoTrack->Momentum;

        if(candidateMomentum.DeltaR(IsoTrackMomentum) <= fDeltaRMax &&
            !candidate->Overlaps(IsoTrack))
        {
          sum += IsoTrackMomentum.Pt();
          ++counter;
        }
      }
    }

    ratio = sum/candidateMomentum.Pt();

    candidate->TrackIsolationVar = ratio;

    if (IsEM) candidate->IsEMCand = 1;

    if((fUsePTSum && sum > fPTSumMax) || ratio > fPTRatioMax) continue;

    if (candidate->Momentum.Pt() > fIsoTrackPTMin )
    {
      fOutputArray->Add(candidate);
    }
  }
}

