/** \class ModifyBeamSpot
 *  \author S. Zenz
*/

#include "modules/ModifyBeamSpot.h"

//#include "CLHEP/Units/GlobalSystemOfUnits.h"
//#include "CLHEP/Units/GlobalPhysicalConstants.h"

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

// static const double mm  = 1.;
// static const double m = 1000.*mm;
// static const double ns  = 1.;
// static const double s = 1.e+9 *ns;
// static const double c_light   = 2.99792458e+8 * m/s;

static const double c_light     = 0.299792458 ; // mm/ps , or m/ns
static const double inv_c_light = 3.335640952 ; // ps/mm or ns/m


using namespace std;

//------------------------------------------------------------------------------

ModifyBeamSpot::ModifyBeamSpot() :
  fFormula(0), fItInputArray(0), fItVertexInputArray(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

ModifyBeamSpot::~ModifyBeamSpot()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void ModifyBeamSpot::Init()
{
  // read resolution formula

  fZVertexSpread = GetDouble ("ZVertexSpread", 53)  ;
  // mm in the cfg file, mm in the code
  fTVertexSpread = GetDouble ("TVertexSpread", 160) ;
  // ps in the cfg file, ps in the code

  fOutputBSX = GetDouble ("OutputBSX", 0.) ;
  fOutputBSY = GetDouble ("OutputBSY", 0.) ;
  fOutputBSZ = GetDouble ("OutputBSZ", 0.) ;
  fOutputBST = GetDouble ("OutputBST", 0.) ;

  // import input array
  
  //fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fVertexInputArray = ImportArray(GetString("VertexInputArray", "PileUpMerger/vertices"));
  fInputArray = ImportArray(GetString("InputArray", "PileUpMerger/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();
  fItVertexInputArray = fVertexInputArray->MakeIterator();
  

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
  fPVOutputArray = ExportArray(GetString("PVOutputArray", "PV"));
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));
  fDebugOutputCollector.addVariable ("PUindex") ;
  
}

//------------------------------------------------------------------------------

void ModifyBeamSpot::Finish()
{  
  std::string outfile = GetString ("simpleOutputFileName", "simpleOutput_Ca.root") ;
  fDebugOutputCollector.save (outfile) ;

  if(fItInputArray)       delete fItInputArray;
  if(fItVertexInputArray) delete fItVertexInputArray;
}

//------------------------------------------------------------------------------

void ModifyBeamSpot::Process()
{
  Candidate *candidate, *mother, *vertexcandidate, *vertexmother;
  // Average position of primary particles
  Double_t PVX = 0., PVY = 0., PVZ = 0., PVT = 0.; 
  int VertexID;
  
  
  PVX = fOutputBSX ;
  PVY = fOutputBSY ;
  PVZ = gRandom->Gaus (0., fZVertexSpread) + fOutputBSZ ;
  PVT = gRandom->Gaus (0., fTVertexSpread * c_light) + fOutputBST * c_light;

  fItInputArray->Reset () ;

  //fItVertexInputArray->Reset();
  // loop over particles in the event
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;

    double X = candidatePosition.X () ;
    double Y = candidatePosition.Y () ;
    double Z = candidatePosition.Z () ;
    double T = candidatePosition.T () ;
    //double dz = gRandom->Gaus(0.0, 4.0);
    // LHE files are centered in (0,0,0,0)
    // while PU has been generated with smearing already, therefore I smear only 
    // the hard scattering.
    if (candidate->IsPU == 0)
      {
      //cout<<"ifPU:"<<"\t"<<candidate->IsPU<<endl;
        X += PVX ;
        Y += PVY ;
        Z += PVZ ;
        T += PVT ;
        //fDebugOutputCollector.fillContainer ("PUindex", PVT) ;
    }       

    //fDebugOutputCollector.fillContainer ("PUindex", candidate->IsPU) ;

    mother = candidate ;
    candidate = static_cast<Candidate*> (candidate->Clone ()) ;
    candidate->Position.SetXYZT (X, Y, Z, T) ;
    candidate->AddCandidate (mother) ;
    //cout<<candidate->VertexID_gen<<endl;
    fOutputArray->Add (candidate) ;
//cout<<"while1"<<endl;
//fDebugOutputCollector.fillContainer ("PUindex", candidate->Position.T()*inv_c_light) ;
  } // loop over particles in the event
//cout<<"fuoriwhile1"<<endl;
  // Store the PV "beam spot"
  DelphesFactory *factory = GetFactory () ;
  candidate               = factory->NewCandidate () ;
  candidate->Position.SetXYZT (PVX, PVY, PVZ, PVT) ;
  fPVOutputArray->Add (candidate) ;
  
  
    fItVertexInputArray->Reset();
    while((vertexcandidate = static_cast<Candidate*>(fItVertexInputArray->Next()))) {
//        if(vertexcandidate->IsPU == 0)  {
        const TLorentzVector &vertexcandidatePosition = vertexcandidate->Position;
//cout<<"vertexwhile"<<endl;
    double X = vertexcandidatePosition.X () ;
    double Y = vertexcandidatePosition.Y () ;
    double Z = vertexcandidatePosition.Z () ;
    double T = vertexcandidatePosition.T () ;
        if(vertexcandidate->IsPU == 0)  {
  //      cout<<"vertexif"<<endl;
            X += PVX;
            Y += PVY;
            Z += PVZ;
            T += PVT;
            fDebugOutputCollector.fillContainer ("PUindex", PVT) ;
            }//
            vertexmother = vertexcandidate ;
            vertexcandidate = static_cast<Candidate*> (vertexcandidate->Clone ()) ;
            //vertexcandidate = static_cast<Candidate*> (vertexcandidate->Next ()) ;
            vertexcandidate->Position.SetXYZT (X, Y, Z,T) ;
            vertexcandidate->AddCandidate (vertexmother) ;
            //cout<<"vertexcandidate:"<<"\t"<<vertexcandidate->Position.T()*inv_c_light<<endl;
            fVertexOutputArray->Add (vertexcandidate) ;
            //if(vertexcandidate->IsPU == 0)  fDebugOutputCollector.fillContainer ("PUindex", vertexcandidate->Position.T()*inv_c_light) ;
       // }
        }
/*
  vertexmother = vertexcandidate;
  vertexcandidate = 
  factory = GetFactory();
  vertexcandidate = factory->NewCandidate();
  vertexcandidate->Position.SetXYZT ( PVX, PVY, PVZ, PVT);
  vertexcandidate->VertexID_gen = VertexID;
  fVertexOutputArray->Add(vertexcandidate);*/
  //VertexID++;

  ++fEventCounter ;   
/*
  while((vertexcandidate = static_cast<Candidate*>(fItVertexInputArray->Next())))
  {
    int VertexID;
    const TLorentzVector &vertexcandidatePosition = vertexcandidate->Position;

    double X = vertexcandidatePosition.X () ;
    double Y = vertexcandidatePosition.Y () ;
    double Z = vertexcandidatePosition.Z () ;
    double T = vertexcandidatePosition.T () ;
    //double dz = gRandom->Gaus(0.0, 4.0);
    // LHE files are centered in (0,0,0,0)
    // while PU has been generated with smearing already, therefore I smear only 
    // the hard scattering.
    if (vertexcandidate->IsPU == 0)
      {
        //float dz = gRandom->Gaus(0.0, 4.0);
        X += PVX ;
        Y += PVY ;
        Z += PVZ ;
        T += PVT ;
       }

    vertexmother = vertexcandidate ;
    vertexcandidate = static_cast<Candidate*> (vertexcandidate->Clone ()) ;
    vertexcandidate->Position.SetXYZT (X, Y, Z, T) ;
    vertexcandidate->AddCandidate (vertexmother) ;

    fVertexOutputArray->Add (vertexcandidate) ;
*/

    /*factory = GetFactory();
    vertexcandidate = factory->NewCandidate();
    vertexcandidate->Position.SetXYZT ( X, Y, Z, T);
    vertexcandidate->VertexID_gen = VertexID;
    fVertexOutputArray->Add(vertexcandidate);
    VertexID++;
  cout<<"while2"<<endl;
  
        */                                                                                                                                    
//}

//cout<<"fuoriwhile2"<<endl;
}

//------------------------------------------------------------------------------
