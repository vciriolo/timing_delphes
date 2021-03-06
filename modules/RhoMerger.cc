#include "modules/RhoMerger.h"

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

//---------------------------//
RhoMerger::RhoMerger()
{
}

//------------------------------------------------------------------------------

RhoMerger::~RhoMerger()
{
}

//------------------------------------------------------------------------------

void RhoMerger::Init()  {

/*ExRootConfParam param = GetParam("RhoInputArray");
  Long_t i, size;
  const TObjArray *array;
  TIterator *iterator;

  size = param.GetSize();

  for(i = 0; i < size; ++i){
    array = ImportArray(param[i].GetString());
    iterator = array->MakeIterator();

    fInputList.push_back(iterator);
  } */
  
  fECalRhoInputArray   = ImportArray(GetString("EcalRhoInputArray", "PvRhoKt4/rho"));
  fItECalRhoInputArray = fECalRhoInputArray->MakeIterator();

fHCalRhoInputArray   = ImportArray(GetString("HCalRhoInputArray", "2PvRhoKt4/2rho"));
  fItHCalRhoInputArray = fHCalRhoInputArray->MakeIterator();
// create output arrays

  fOutputArray = ExportArray(GetString("OutputArray", "rho"));
}

//------------------------------------------------------------------------------

void RhoMerger::Finish()
{
/*
  vector< TIterator * >::iterator itInputList;
  TIterator *iterator;

  for(itInputList = fInputList.begin(); itInputList != fInputList.end(); ++itInputList)
  {
    iterator = *itInputList;
    if(iterator) delete iterator;
  } */
  
  if(fItECalRhoInputArray) delete fItECalRhoInputArray;
  if(fItHCalRhoInputArray) delete fItHCalRhoInputArray;
}

//------------------------------------------------------------------------------
void RhoMerger::Process()   {
    Candidate *candidate, *candidate2, *rhocandidate;
    TLorentzVector rho;
    Double_t tmpRhoPT, tmpRhoE, tmpRho = 0.;  
    vector< TIterator * >::iterator itInputList;
    TIterator *iterator;

    DelphesFactory *factory = GetFactory();
  
  /*  momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    sumPT = 0;
    sumE = 0;   */

  // loop over all input arrays
  /*for(itInputList = fInputList.begin(); itInputList != fInputList.end(); ++itInputList)
  {
    iterator = *itInputList;
    iterator->Reset();
    while((candidate = static_cast<Candidate*>(iterator->Next())))  {
        tmpRhoPT += candidate->Momentum.Pt();
        tmpRhoE += candidate->Momentum.E();
        cout<<tmpRhoPT<<endl;
        }
    }   */
    int ecal_size = fECalRhoInputArray->GetEntriesFast();
    for (int i = 0; i != ecal_size; i++)   {
    //while((candidate = static_cast<Candidate*>(fItECalRhoInputArray->Next())))  {
        candidate = static_cast<Candidate*>(fECalRhoInputArray->At(i));
        tmpRhoPT += candidate->Momentum.Pt();
        tmpRhoE += candidate->Momentum.E();
        candidate2 = static_cast<Candidate*>(fHCalRhoInputArray->At(i));
        tmpRhoPT += candidate2->Momentum.Pt();
        tmpRhoE += candidate2->Momentum.E();
        rhocandidate = factory->NewCandidate();
  
  rhocandidate->Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  rhocandidate->Momentum.SetPtEtaPhiE(tmpRhoPT,0.0,0.0,tmpRhoE);
  fOutputArray->Add(rhocandidate);
  tmpRhoPT = 0.;
  tmpRhoE = 0.;
       }
    
  /*  for (unsigned i = fHCalRhoInputArray.begin(); i != fHCalRhoInputArray.end(); i++)   {
        candidate = fHCalRhoInputArray.at(i);
        tmpRhoPT += candidate->Momentum.Pt();
        tmpRhoE += candidate->Momentum.E();
    }   */
        
        
    
   /* candidate = factory->NewCandidate();
  
  candidate->Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  candidate->Momentum.SetPtEtaPhiE(tmpRhoPT,0.0,0.0,tmpRhoE);
  fOutputArray->Add(candidate); */
}
    
    






