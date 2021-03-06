#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TTreeReader.h" 
#include "TProfile.h"
 
void time_correction()
{
   // Variables used to store the data
  TProfile *prof = new TProfile ("time_vs_eta","timespread_vs_eta",200, -2.,2.,3800.,12000.);
    corr_timeTower = new TF1 ( "corr_timeTower", "22.4 * [0] + [1]*TMath::Abs(x)", -1.479,1.479);
 // open the file
   TFile *f = TFile::Open("simpleOutput_ECal.root");
   if (f == 0) {
      // if we cannot open the file, print an error message and return immediatly
      printf("Error: cannot open !\n");
      return;
   }
 
   // Create tyhe tree reader and its data containers
   TTreeReader myReader("Nm_eta_time_pt_size", f);
 
TTreeReaderValue<Double_t> time(myReader, "fir");
TTreeReaderValue<Double_t> eta(myReader, "sec");

time_vs_eta->Fill (time, eta);
time_vs_eta->Fit("corr_timeTower");
time_vs_eta->Draw();       
   
}
