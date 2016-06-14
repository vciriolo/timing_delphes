/****************************************************************************************************************************
        -this program create a sample on minimum bias events for SLHC-14TeV
        -reference for the Tune:pp option used can be found here: https://cds.cern.ch/record/1697700/files/GEN-14-001-pas.pdf 
        -source minbias_setup_slc6, then:

        -compile with ---> gcc -I$HEPMC/include/ -L$HEPMC/lib -I$PYTHIA8DATA/../include -L$PYTHIA8DATA/../lib -I$LHAPDF/include -L$LHAPDF/lib 
                               -lLHAPDF -lHepMC -lpythia8tohepmc -lpythia8 -o genMinBias_14TeV genMinBias_14TeV.cpp
****************************************************************************************************************************/

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "Pythia8Plugins/LHAPDF5.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/IO_AsciiParticles.h"

using namespace Pythia8;  
using namespace std;

//********************************************************************************************
int main(int argc, char** argv) {

    if(argc < 4){
        cout << "missing argument!" << endl
             << "-----------------" << endl
	     << "usage: " << endl
             << "genMinBias_14TeV events_number output_file tune" << endl
             << "-----------------" << endl << endl;
        return 0;
    }

    int nEvents    = atoi(argv[1]);
    string outFile = argv[2];
    int tunePythia = atoi(argv[3]);

    //-----Pythia setup-----
    Pythia pythia;     

    //---general setup
    pythia.readString("Random:seed = 0");
    pythia.readString("Random:setSeed = on"); 

    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eA  = 7000.");
    pythia.readString("Beams:eB  = 7000.");
    pythia.readString("Beams:eCM = 14000.");

    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("SoftQCD:singleDiffractive = on");
    pythia.readString("SoftQCD:doubleDiffractive = on");

    pythia.readString("HadronLevel:Hadronize = on");


    pythia.readString("Check:epTolErr = 0.01");
    pythia.readString("SLHA:keepSM    = on");
    pythia.readString("SLHA:minMassSM = 1000.");

    //---CMS tuned min bias
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("ParticleDecays:tau0Max = 10");
    pythia.readString("ParticleDecays:allowPhotonRadiation = on");

    if( tunePythia == 5){

     pythia.readString("MultipartonInteractions:pT0Ref = 2.1006");
     pythia.readString("MultipartonInteractions:ecmPow = 0.21057");
     pythia.readString("MultipartonInteractions:expPow = 1.6089");
     pythia.readString("MultipartonInteractions:a1     = 0.00");

     pythia.readString("Tune:pp  = 5");  // cms tune   
     pythia.readString("Tune:ee  = 3");

     pythia.readString("PDF:pSet = LHAPDF5:cteq6ll.LHpdf");

    }
    else if (tunePythia == 14){

      pythia.readString("MultipartonInteractions:pT0Ref=2.4024");
      pythia.readString("MultipartonInteractions:ecmPow=0.25208");
      pythia.readString("MultipartonInteractions:expPow=1.6");
      pythia.readString("Tune:pp = 14");
      pythia.readString("Tune:ee = 7");

    }
    else if (tunePythia == 15){

     pythia.readString("MultipartonInteractions:pT0Ref = 2.1006");
     pythia.readString("MultipartonInteractions:ecmPow = 0.21057");
     pythia.readString("MultipartonInteractions:expPow = 1.6089");
     pythia.readString("MultipartonInteractions:a1     = 0.00");

     pythia.readString("Tune:pp  = 15");  // cms tune   
     pythia.readString("Tune:ee  = 3");

     pythia.readString("PDF:pSet = LHAPDF5:cteq6ll.LHpdf");

    }

    //---Pythia initialization
    pythia.init();

    HepMC::IO_GenEvent hepmc_file_out(outFile, std::ios::out);
    HepMC::Pythia8ToHepMC ToHepMC;
	
    for (int iEvent = 0; iEvent < nEvents; ++iEvent){
	if(iEvent%(nEvents/10) == 0) {
	    cout << "Events:  " << iEvent << endl;
	}
	if(!pythia.next()){
	    cout << "CRASH! ---> skip" << endl;
	    continue;
	}
	// construct new HepMC event setting units.	
	HepMC::GenEvent* hepmcevt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
	// fill the event including PDF infos
	ToHepMC.fill_next_event( pythia, hepmcevt );
	hepmc_file_out << hepmcevt;
        delete hepmcevt;
    }    

    // Done.                           
    return 0;
}
