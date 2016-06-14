#ifndef ModifyBeamSpot_h
#define ModifyBeamSpot_h

/** \class ModifyBeamSpot
 *
 *  \author S. Zenz
 *
 */

#include "classes/DelphesModule.h"
#include "modules/simpleVariableCollector.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class ModifyBeamSpot: public DelphesModule
{
public:

  ModifyBeamSpot();
  ~ModifyBeamSpot();

  void Init();
  void Process();
  void Finish();

private:

  DelphesFormula *fFormula; //!

  TIterator *fItInputArray; //!
  TIterator * fItVertexInputArray;
  const TObjArray *fInputArray; //!
  const TObjArray *fVertexInputArray;
  TObjArray *fOutputArray; //!
  TObjArray *fVertexOutputArray;

  Double_t fOutputBSX, fOutputBSY, fOutputBSZ, fOutputBST ;
  Double_t fZVertexSpread;
  Double_t fTVertexSpread;
  Double_t currentZ, currentT;

  // Store Z of PV
  TObjArray *fPVOutputArray; //!

  simpleVariableCollector fDebugOutputCollector ;    
  int fEventCounter ;                                                                                                                                             

  ClassDef(ModifyBeamSpot, 1)
};

#endif
