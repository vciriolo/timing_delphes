#ifndef RhoMerger_h
#define RhoMerger_h


#include "classes/DelphesModule.h"

#include <vector>

class TIterator;
class TObjArray;
class DelphesFormula;

class RhoMerger: public DelphesModule
{
public:

  RhoMerger();
  ~RhoMerger();

  void Init();
  void Process();
  void Finish();

private:

  std::vector< TIterator * > fInputList; //!

  TObjArray *fOutputArray; //!
  TObjArray *fECalRhoInputArray;
  TObjArray *fHCalRhoInputArray;
  TIterator *fItECalRhoInputArray;
  TIterator *fItHCalRhoInputArray;

  ClassDef(RhoMerger, 1)
};

#endif
