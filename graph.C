void graph()

{
	c1 = new TCanvas("c1","molt",200,10,700,500);
	//gr1 = new TGraphErrors("vbs/njets_tcut.txt");
	gr1 = new TGraph("vbs/pt_res_tcut.txt");  
	//c1->SetLogy();
       gr1->SetMarkerColor(kBlue);
       gr1->SetMarkerStyle(21);
       gr1->Draw("APC");
	gr1->SetTitle("Pt vs TCut value ");
	
  gr1->GetXaxis()->SetTitle("TCut (ps)");
  //gr1->GetYaxis()->SetTitle("(Pt_reco - Pt_gen)/Pt_reco");
  gr1->GetYaxis()->SetTitle("(#sigma_{Chs} - #sigma_{TCut})/#sigma_{Chs}");
  //gr1->GetYaxis()->SetRange(9,13);

//	return c1;
}
