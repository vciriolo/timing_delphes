TCanvas *c1 = new TCanvas ()
easyDelphes->Draw("pvefrac_jet>>h(200,-1,1)","TMath::Abs(reco_dR)<0.5","")
easyDelphes->Draw("efrac_jet>>h1(200,-1,1)","TMath::Abs(reco_dR)<0.5","same")
TCanvas *c2 = new TCanvas ()
easyDelphes->Draw("pvefrac_jet>>h2(200,-1,1)","TMath::Abs(reco_dR)>0.5 && TMath::Abs(reco_dR)<1","")
easyDelphes->Draw("efrac_jet>>h3(200,-1,1)","TMath::Abs(reco_dR)>0.5 && TMath::Abs(reco_dR)<1","same")
TCanvas *c3 = new TCanvas ()
easyDelphes->Draw("pvefrac_jet>>h4(200,-1,1)","TMath::Abs(reco_dR)>1. && TMath::Abs(reco_dR)<1.479","")
easyDelphes->Draw("efrac_jet>>h5(200,-1,1)","TMath::Abs(reco_dR)>1. && TMath::Abs(reco_dR)<1.479","same")
TCanvas *c4 = new TCanvas ()
easyDelphes->Draw("pvefrac_jet>>h6(200,-1,1)","TMath::Abs(reco_dR)>1.479","")
easyDelphes->Draw("efrac_jet>>h7(200,-1,1)","TMath::Abs(reco_dR)>1.479","same")
