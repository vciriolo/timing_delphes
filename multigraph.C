
#include "TCanvas.h"
#include "math.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TROOT.h"
#include "TStyle.h"
#include <iostream>
#include "TLegend.h"

    void multigraph()
    {
       //TCanvas *c = new TCanvas("c","Equivalent MIPS vs X0 for 0.5 cm^2 pad",600, 400);
       TCanvas *c = new TCanvas();
	//c->SetGrid();
	//c->SetLogy();
       TMultiGraph * mg = new TMultiGraph("mg","Equivalent MIPS vs X0, 0.5 cm^2 pad");

       /*const Int_t size = 10;
              
       double x[size];
       double y1[size];
       double y2[size];
       double y3[size];

       for ( int i = 0; i <  size ; ++i ) {
          x[i] = i;
          y1[i] = size - i;
          y2[i] = size - 0.5 * i;
          y3[i] = size - 0.6 * i;
       }
	*/
       TGraphErrors * gr1 = new TGraphErrors("ptres_pt_chs.txt");
       gr1->SetName("gr1");
       gr1->SetTitle("chs");
       gr1->SetMarkerStyle(21);
       gr1->SetMarkerColor(kRed);
       gr1->SetDrawOption("ACP");
       gr1->SetLineColor(kRed);
       gr1->SetLineWidth(2);
       gr1->SetFillStyle(0);
       gr1->GetXaxis()->SetRange(65,100);
       
	//gr1->Draw("alp");
	//mg->Add(gr1);
	gr1->SetTitle("lead");

      TGraphErrors * gr2 = new TGraphErrors("ptres_pt.txt");
       gr2->SetName("gr2");
       gr2->SetTitle("putrid");
       gr2->SetMarkerStyle(21);
       gr2->SetMarkerColor(4);
       gr2->SetDrawOption("ACP");
       gr2->SetLineColor(kBlue);
       gr2->SetLineWidth(2);
       gr2->SetFillStyle(0);
       gr2->GetXaxis()->SetRange(65,100);       
	//gr2->Draw("same");
	//mg->Add(gr2);
       /*TGraph * gr3 = new TGraph( size, x, y3 );
       gr3->SetName("gr3");
       gr3->SetTitle("graph 3");
       gr3->SetMarkerStyle(23);
       gr3->SetLineColor(4);
       gr3->SetLineWidth(4);
       gr3->SetFillStyle(0);
*/
       mg->Add( gr1 );
       mg->Add( gr2 );
	
  //     gr3->Draw("ALP");
       leg = new TLegend(0.1,0.7,0.48,0.9);
       leg->AddEntry(gr1,"CHS","lp");	
       leg->AddEntry(gr2,"CHS + Time Cut","lp");
	   leg->SetFillStyle(0);
	   leg->SetBorderSize(0);
	mg->Draw("alp");
	leg->Draw("same");
//       c->BuildLegend();

       //c->Print("multigraphleg.gif");
    }


