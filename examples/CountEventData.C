/*
root -l examples/Example1.C\(\"delphes_output.root\"\)

or just:

./CountEventData delphes_output.root

N.B. you must touch CountEventData.cpp in order to convince Make that there is anything to do to build this!
*/

//------------------------------------------------------------------------------

void CountEventData(const char *inputFile, const char *outputFile)
{
  //  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  bool verbose = true;
  
  TFile *out = new TFile(outputFile,"RECREATE");
  TTree *tree = new TTree("tree","tree");

  int nj = 0; // 1,2,3,4
  int ng = 0; // 9,21
  int nb = 0; // 5                                                                                                                                                             
  int nt = 0; // 6                                                                                                                                                             
  int ngam = 0; // 22                                                                                                                                                          
  int nZ = 0; // 23                                                                                                                                                           
  int nW = 0; // 24                                                                                                                                                            
  int nH = 0; // 25                                                                                                                                                            
  int ne = 0; //11                                                                                                                                                             
  int nmu = 0; // 13                                                                                                                                                           
  int ntau = 0; // 15                                                                                                                                                          
  int nv = 0;  // 12 14 16                                                                                                                                                     
  int py_nj = 0; // 1,2,3,4
  int py_ng = 0; // 9,21
  int py_nb = 0; // 5
  int py_nt = 0; // 6
  int py_ngam = 0; // 22
  int py_nZ = 0; // 23
  int py_nW = 0; // 24
  int py_nH = 0; // 25
  int py_ne = 0; //11
  int py_nmu = 0; // 13
  int py_ntau = 0; // 15
  int py_nv = 0;  // 12 14 16

  float weight = -1;

  tree->Branch("nj",&nj,"nj/I");
  tree->Branch("ng",&ng,"ng/I");
  tree->Branch("nb",&nb,"nb/I");
  tree->Branch("nt",&nt,"nt/I");
  tree->Branch("ngam",&ngam,"ngam/I");
  tree->Branch("nZ",&nZ,"nZ/I");
  tree->Branch("nW",&nW,"nW/I");
  tree->Branch("nH",&nH,"nH/I");
  tree->Branch("ne",&ne,"ne/I");
  tree->Branch("nmu",&nmu,"nmu/I");
  tree->Branch("ntau",&ntau,"ntau/I");
  tree->Branch("nv",&nv,"nv/I");
  tree->Branch("py_nj",&py_nj,"py_nj/I");
  tree->Branch("py_ng",&py_ng,"py_ng/I");
  tree->Branch("py_nb",&py_nb,"py_nb/I");
  tree->Branch("py_nt",&py_nt,"py_nt/I");
  tree->Branch("py_ngam",&py_ngam,"py_ngam/I");
  tree->Branch("py_nZ",&py_nZ,"py_nZ/I");
  tree->Branch("py_nW",&py_nW,"py_nW/I");
  tree->Branch("py_nH",&py_nH,"py_nH/I");
  tree->Branch("py_ne",&py_ne,"py_ne/I");
  tree->Branch("py_nmu",&py_nmu,"py_nmu/I");
  tree->Branch("py_ntau",&py_ntau,"py_ntau/I");
  tree->Branch("py_nv",&py_nv,"py_nv/I");

  tree->Branch("weight",&weight,"weight/F");

  float minphopt = 9999.;
  float minphoe = 9999.;

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    if (verbose||entry%5000==0) cout << "Event " << entry << " / " << numberOfEntries << endl;
    
    nj = 0;
    ng = 0;
    nb = 0;
    nt = 0;
    ngam = 0;
    nZ = 0;
    nW = 0;
    nH = 0;
    ne = 0;
    nmu = 0;
    ntau = 0;
    nv = 0;  
    py_nj = 0;
    py_ng = 0;
    py_nb = 0;
    py_nt = 0;
    py_ngam = 0;
    py_nZ = 0;
    py_nW = 0;
    py_nH = 0;
    py_ne = 0;
    py_nmu = 0;
    py_ntau = 0;
    py_nv = 0;


    for (int i = 0 ; i < branchGenParticle->GetEntries() ; i++ ) {
      GenParticle *part = (GenParticle*) branchGenParticle->At(i);
      if (verbose) {
	cout << "     Status code 3 generator particle PID Pt Eta Phi Z T (at origin) "  << part->PID << " "
	     << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Z << " " << part->T << endl;
      }
      signed int apid = abs(part->PID);
      if (fabs(part->Z) < 1e-20) {
	if (apid >= 1 && apid <= 4) nj++;
	if (apid == 9 || apid == 21) ng++;
	if (apid == 5) nb++;
	if (apid == 6) nt++;
	if (apid == 11) ne++;
	if (apid == 13) nmu++;
	if (apid == 15) ntau++;
	if (apid == 12 || apid == 14 || apid == 16) nv++;
	if (apid == 22) {
	  ngam++;
	  if (part->PT < minphopt) { 
	    minphopt = part->PT;
	    cout << " New minphopt " << minphopt << endl;
	  }
          if (part->E < minphoe) {
            minphoe = part->E;
            cout << " New minphoe " << minphoe << endl;
          }

	}
	if (apid == 23) nZ++;
	if (apid == 24) nW++;
	if (apid == 25) nH++;
      } else {
        if (apid >= 1 && apid <= 4) py_nj++;
        if (apid == 9 || apid == 21) py_ng++;
        if (apid == 5) py_nb++;
        if (apid == 6) py_nt++;
        if (apid == 11) py_ne++;
        if (apid == 13) py_nmu++;
        if (apid == 15) py_ntau++;
        if (apid == 12 || apid == 14 || apid == 16) py_nv++;
        if (apid == 22) py_ngam++;
        if (apid == 23) py_nZ++;
        if (apid == 24) py_nW++;
	if (apid == 25) py_nH++;
      }
    }

    weight = -1.;

    if (branchEvent) {
      LHEFEvent *ev = (LHEFEvent*) branchEvent->At(0);
      weight = ev->Weight;
      if (verbose) {
	cout << "Event weight is: " << weight << endl;
      }
    }

    cout << "Listing   nj ng nb nt ne nmu ntau nv ngam nZ nW nH weight" << endl;
    cout << "MadGraph: " << nj << "  " << ng << "  " << nb << "  " << nt << "  " << ne << "  " << nmu << "   " << ntau << "    " << nv << "  " 
	 << ngam << "    " << nZ << "  " << nW << "  " << nH << "  " << weight << endl;
    cout << "Pythia:   " << py_nj << "  " << py_ng << "  " << py_nb << "  " << py_nt << "  " << py_ne << "  " << py_nmu << "   " << py_ntau << "    " << py_nv << "  "
	 << py_ngam << "    " << py_nZ << "  " << py_nW << "  " << py_nH << "  " << weight << endl;
    tree->Fill();

  } // per event

  cout << "minphopt: " << minphopt << endl;
  cout << "minphoe: " << minphoe << endl;


  tree->Write();
  out->Close();
      
}
