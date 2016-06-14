//------------------------------------------------------------------------------

#ifndef GenBeamSpotFilter_h
#define GenBeamSpotFilter_h

/** \class Efficiency
 *
 *  Removes all generated particles except electrons, muons, taus,
 *  and particles with status == 3.
 *
 *  $Date: 2013-04-07 00:12:34 +0200 (Sun, 07 Apr 2013) $
 *  $Revision: 1079 $
 *
 *
 *  \author J. Hirschauer - FNAL
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;

class GenBeamSpotFilter: public DelphesModule
{
public:

  GenBeamSpotFilter();
  ~GenBeamSpotFilter();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fPTMin; //!

  Float_t fPassedOne;

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(GenBeamSpotFilter, 1)
};

#endif
