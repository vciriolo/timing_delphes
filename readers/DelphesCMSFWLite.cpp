#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <memory>

#include <map>
#include <vector>

#include <stdlib.h>
#include <signal.h>
#include <stdio.h>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

#include "modules/Delphes.h"
#include "classes/DelphesStream.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LesHouches.h"


using namespace std;

//---------------------------------------------------------------------------

void ConvertInput(fwlite::Event &event, Long64_t eventCounter, ExRootTreeBranch *branch, DelphesFactory *factory, TObjArray *allParticleOutputArray, TObjArray *stableParticleOutputArray, TObjArray *partonOutputArray)
{
    LHEFEvent *lheEvt;
    lheEvt = static_cast<LHEFEvent *>(branch->NewEntry());

    fwlite::Handle<LHEEventProduct> lheEvtInfo;

    lheEvtInfo.getByLabel(event, "source");

    if (lheEvtInfo.isValid()) {
      lheEvt->Number = eventCounter;
      lheEvt->Weight = lheEvtInfo->originalXWGTUP();
      lheEvt->ProcessID = ((lhef::HEPEUP)lheEvtInfo->hepeup()).IDPRUP;
      lheEvt->ScalePDF = ((lhef::HEPEUP)lheEvtInfo->hepeup()).IDPRUP;
      lheEvt->AlphaQED = ((lhef::HEPEUP)lheEvtInfo->hepeup()).SCALUP;
      lheEvt->AlphaQCD = ((lhef::HEPEUP)lheEvtInfo->hepeup()).AQCDUP;
    }
    else {
      lheEvt->Number = 0;
      lheEvt->Weight = 1;
      lheEvt->ProcessID = 0;
      lheEvt->ScalePDF = 0;
      lheEvt->AlphaQED = 0;
      lheEvt->AlphaQCD = 0;
    }

    fwlite::Handle< vector< reco::GenParticle > > handleParticle;
    vector< reco::GenParticle >::const_iterator itParticle;

    vector< const reco::Candidate * > vectorCandidate;
    vector< const reco::Candidate * >::iterator itCandidate;

    handleParticle.getByLabel(event, "genParticles");

    Candidate *candidate;
    TDatabasePDG *pdg;
    TParticlePDG *pdgParticle;
    Int_t pdgCode;

    Int_t pid, status;
    Double_t px, py, pz, e, mass;
    Double_t x, y, z;

    pdg = TDatabasePDG::Instance();

    for(itParticle = handleParticle->begin(); itParticle != handleParticle->end(); ++itParticle)
    {
        vectorCandidate.push_back(&*itParticle);
    }

    for(itParticle = handleParticle->begin(); itParticle != handleParticle->end(); ++itParticle)
    {
        const reco::GenParticle &particle = *itParticle;

        pid = particle.pdgId();
        status = particle.status();
        px = particle.px(); py = particle.py(); pz = particle.pz(); e = particle.energy(); mass = particle.mass();
        x = particle.vx(); y = particle.vy(); z = particle.vz();

        candidate = factory->NewCandidate();

        candidate->PID = pid;
        pdgCode = TMath::Abs(candidate->PID);

        candidate->Status = status;

        //M1
        if(particle.mother()){
            itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), particle.mother());
            if(itCandidate != vectorCandidate.end()) candidate->M1 = distance(vectorCandidate.begin(), itCandidate);
        }

        //D1
        itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), particle.daughter(0));
        if(itCandidate != vectorCandidate.end()) candidate->D1 = distance(vectorCandidate.begin(), itCandidate);

        //D2
        if(particle.numberOfDaughters() > 1)
            itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), particle.daughter(1));
        else
            itCandidate = find(vectorCandidate.begin(), vectorCandidate.end(), particle.daughter(particle.numberOfDaughters() - 1));

        if(itCandidate != vectorCandidate.end()) candidate->D2 = distance(vectorCandidate.begin(), itCandidate);

        pdgParticle = pdg->GetParticle(pid);
        candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
        candidate->Mass = mass;

        candidate->Momentum.SetPxPyPzE(px, py, pz, e);

        candidate->Position.SetXYZT(x, y, z, 0.0);

        allParticleOutputArray->Add(candidate);

        if(!pdgParticle) continue;

        if(status == 1)
        {
            stableParticleOutputArray->Add(candidate);
        }
        else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
        {
            partonOutputArray->Add(candidate);
        }
    }
}

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
    interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "DelphesCMSFWLite";
  stringstream message;
  TFile *inputFile = 0;
  TFile *outputFile = 0;
  TStopwatch eventStopWatch;
  ExRootTreeWriter *treeWriter = 0;
  ExRootTreeBranch *branchEvent = 0;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;
  TObjArray *allParticleOutputArray = 0, *stableParticleOutputArray = 0, *partonOutputArray = 0;
  Int_t i;
  Int_t maxEvents, skipEvents;
  Long64_t eventCounter, numberOfEvents;
  
  if(argc < 4)
    {
      cout << " Usage: " << appName << " config_file" << " output_file" << " input_file(s)" << endl;
      cout << " config_file - configuration file in Tcl format," << endl;
      cout << " output_file - output file in ROOT format," << endl;
      cout << " input_file(s) - input file(s) in ROOT format." << endl;
      return 1;
    }
  
  signal(SIGINT, SignalHandler);
  
  gROOT->SetBatch();
  
  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);
  
  AutoLibraryLoader::enable();
  
  try
    {
      outputFile = TFile::Open(argv[2], "CREATE");
      
      if(outputFile == NULL)
        {
	  message << "can't open " << argv[2] << endl;
	  throw runtime_error(message.str());
        }
      
      treeWriter = new ExRootTreeWriter(outputFile, "Delphes");
      
      branchEvent = treeWriter->NewBranch("Event", LHEFEvent::Class());
      
      confReader = new ExRootConfReader;
      confReader->ReadFile(argv[1]);
      
      maxEvents = confReader->GetInt("::MaxEvents", 0);
      skipEvents = confReader->GetInt("::SkipEvents", 0);
      
      if (maxEvents<0) {
	throw runtime_error("MaxEvents must be zero or positive");
      }
      if (skipEvents<0){ 
	throw runtime_error("SkipEvents must be zero or positive");
      }
      
      modularDelphes = new Delphes("Delphes");
      modularDelphes->SetConfReader(confReader);
      modularDelphes->SetTreeWriter(treeWriter);

      factory = modularDelphes->GetFactory();
      allParticleOutputArray = modularDelphes->ExportArray("allParticles");
      stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
      partonOutputArray = modularDelphes->ExportArray("partons");
      
      modularDelphes->InitTask();
      
      int totEventCounter = 0;
      
      for(i = 3; i < argc && !interrupted; ++i)
	{
	  cout << "** Reading " << argv[i] << endl;
	  
	  inputFile = TFile::Open(argv[i]);
	  
	  if(inputFile == NULL)
	    {
	      message << "can't open " << argv[i] << endl;
	      throw runtime_error(message.str());
            }
	  
            fwlite::Event event(inputFile);
	    
            numberOfEvents = event.size();
	    
            if(numberOfEvents <= 0) continue;
	    
            ExRootProgressBar progressBar(-1);
	    
            // Loop over all objects
            eventCounter = 0;
            modularDelphes->Clear();
            treeWriter->Clear();
	    
            for(event.toBegin(); !event.atEnd() && !interrupted; ++event) 
	      {
		
		if (eventCounter<skipEvents) 
		  {
		    ++eventCounter;
		    continue;
		  }
		else if (totEventCounter>maxEvents-1 && maxEvents>0) 
		  {
		    break;
		  }
		else 
		  {
		    ConvertInput(event, eventCounter, branchEvent, factory, allParticleOutputArray, stableParticleOutputArray, partonOutputArray);
		    modularDelphes->ProcessTask();
		    
		    treeWriter->Fill();
		    
		    modularDelphes->Clear();
		    treeWriter->Clear();
		    
		    progressBar.Update(totEventCounter, totEventCounter);
		    ++eventCounter;
		    ++totEventCounter;
		  }
	      }

	    progressBar.Update(totEventCounter, totEventCounter, kTRUE);
	    progressBar.Finish();

            inputFile->Close();
        }
      
      modularDelphes->FinishTask();
      treeWriter->Write();
      
      cout << "** Exiting..." << endl;
      
      delete modularDelphes;
      delete confReader;
      delete treeWriter;
      delete outputFile;
      
      return 0;
    }
  catch(runtime_error &e)
    {
      if(treeWriter) delete treeWriter;
      if(outputFile) delete outputFile;
      cerr << "** ERROR: " << e.what() << endl;
      return 1;
    }
}
