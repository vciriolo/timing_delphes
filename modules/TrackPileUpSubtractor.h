#ifndef TrackPileUpSubtractor_h
#define TrackPileUpSubtractor_h

/** \class TrackPileUpSubtractor
 *
 *  Subtract pile-up contribution from tracks.
 *
 *  $Date: 2013-03-24 15:08:05 +0100 (Sun, 24 Mar 2013) $
 *  $Revision: 1069 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"
#include "modules/simpleVariableCollector.h"
#include <map>

class TIterator;
class TObjArray;

class TrackPileUpSubtractor: public DelphesModule
{
public:

  TrackPileUpSubtractor();
  ~TrackPileUpSubtractor();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fZVertexResolution;
  Double_t fPTMin ; 

  std::map< TIterator *, TObjArray * > fInputMap; //!

  TIterator *fPVItInputArray; //!
  const TObjArray *fPVInputArray; //!
  simpleVariableCollector fDebugOutputCollector ;    

  ClassDef(TrackPileUpSubtractor, 1)
};

#endif
