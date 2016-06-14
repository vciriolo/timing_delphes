/** \class Vertexing
*
* Merges multiple input arrays into one output array
* and sums transverse momenta of all input objects.
*
* $Date: 2013-04-26 12:39:14 +0200 (Fri, 26 Apr 2013) $
* $Revision: 1099 $
*
*
* \author P. Demin - UCL, Louvain-la-Neuve
*
*/
#include "modules/Vertexing.h"
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
#include <vector>
using namespace std;
//------------------------------------------------------------------------------
Vertexing::Vertexing():
fpt_min(0)
{
}
//------------------------------------------------------------------------------
Vertexing::~Vertexing()
{
}
//------------------------------------------------------------------------------
void Vertexing::Init()
{
// import arrays with output from other modules
//fParticleInputArray = ImportArray(GetString("ParticleInputArray, stableParticle"));
//fItParticleInputArray = fParticleInputArray->MakeIterator();
fTrackInputArray = ImportArray(GetString("TrackInputArray", "tracks"));
//fTrackInputArray = ImportArray(GetString("TrackInputArray", "stableParticles"));
fVertexInputArray = ImportArray(GetString("VertexInputArray", "vertices"));
//fGenVertexInputArray = ImportArray(GetString("GenVertexInputArray", "stableParticles"));
fItVertexInputArray = fVertexInputArray->MakeIterator();
//fItGenVertexInputArray = fGenVertexInputArray->MakeIterator();
fItTrackInputArray = fTrackInputArray->MakeIterator();
//output array
fVertexingTrackOutputArray = ExportArray(GetString("VertexingTrackOutputArray", "vertexingtracks"));
//fItVertexingTrackOutputArray = fVertexingTrackOutputArray->MakeIterator();
fInitialTrackOutput = ExportArray(GetString("InitialTrackOutput", "initialtrack"));
fDebugOutputCollector.addVariable4D ("imp_par:Pt:P:PID") ;
fDebugOutputCollector.addVariable ("z_at_vert") ;
fDebugOutputCollector.addVariable3D ("candidate_track:track_molt:vertex_time") ;
fDebugOutputCollector.addVariable2D ("candidate_track_cut:track_molt_cut") ;
fpt_min = GetDouble("pt_min",0.1);

}
//------------------------------------------------------------------------------
void Vertexing::Finish()
{
if(fTrackInputArray)      delete fTrackInputArray;
if(fVertexInputArray)     delete fVertexInputArray;
//if(fGenVertexInputArray)  delete fGenVertexInputArray;
std::string outfile = GetString ("simpleOutputFileName", "simpleOutput_VE.root") ;
  fDebugOutputCollector.save (outfile) ;

}
//------------------------------------------------------------------------------
void Vertexing::Process()
{
Candidate *candidate_track, *candidate_vertex, *vertex, *particle , *gen_vertex, *track, *smeared_vertex, *reco_vertex, *initial_track;
TIter track_iter(fTrackInputArray);
TIter vertex_iter(fVertexInputArray);
vector<double> sumPt;

int track_id=0;
int vertex_id=0;
const double_t c_light = 0.299792458E8;  //mm/ps
const double inv_c_light = 3.335640952 ; // ps/mm or ns/m
float imp_par=0.;
float imp_par_tmp=5.;
float z_vertex=0.;
double lepton_energy =0.;
int leptons = 0;
fItTrackInputArray->Reset();
while((candidate_track = static_cast<Candidate*>(fItTrackInputArray->Next())))  {
  /*      
    if(candidate_track->IsPU==0 && ((candidate_track->PID) == 11 || candidate_track->PID == -11) ){//|| candidate_track->PID == 13 || candidate_track->PID == -13){// || fabs(candidate_track->PID)==13)){
        ++leptons;
                lepton_energy += candidate_track->Momentum.E();}
    if(candidate_track->IsPU == 0 && (candidate_track->PID == 13 || candidate_track->PID == -13)){  
    lepton_energy += candidate_track->Momentum.E();
    cout<<"hhhh"<<endl;
    ++leptons;}                
                */
   //fDebugOutputCollector.fillContainer ("z_at_vert", candidate_track->Position.T()*inv_c_light);
   track_id = candidate_track->VertexID_gen;
   //track position in calorimeter
   const TLorentzVector &trackposition = candidate_track->Position;
   float x_track = trackposition.X();
   float y_track = trackposition.Y();
   float z_track = trackposition.Z();
   float t_track = trackposition.T();
   //fDebugOutputCollector.fillContainer ("imp_par", candidate_track->Position.Z()*inv_c_light);  
   vertex_iter.Reset();
   while((vertex = static_cast<Candidate*>(vertex_iter.Next()))) {
    //cout<<"vertexT:"<<"\t"<<vertex->Position.T()<<endl;
    vertex_id = vertex->VertexID_gen;
    //vertex position
    smeared_vertex = static_cast<Candidate*>(vertex->Clone());
    float dz_vertex = gRandom->Gaus (0.0, 0.014);
    const TLorentzVector &nonsmeared = vertex->Position;
    float x = nonsmeared.X();
    float y = nonsmeared.Y();
    float z = nonsmeared.Z();
    float t = nonsmeared.T();
    //z += dz_vertex;
    //smear vertex
    smeared_vertex->Position.SetXYZT(x, y, z , t);
    
    int flag = 0;
    if(track_id == vertex_id)   {
        flag = 1;
        //track position = non smeared vertex position and then smear track position
        const TLorentzVector &initialPosition = vertex->Position;
        float x_vert = initialPosition.X();
        float y_vert = initialPosition.Y();
        float z_vert = initialPosition.Z();
        //z_vert *= 1.0E3;
        //float t_vert = initialPosition.T();
        //t_vert *= inv_c_light;
        //smear track z
        float dz_track= gRandom -> Gaus (0.0, 0.065);
        z_vert += dz_track;
        candidate_track->Position.SetX(x_vert);
        candidate_track->Position.SetY(y_vert);
        candidate_track->Position.SetZ(z_vert); 
//        initial_track = static_cast<Candidate*>((candidate_track->Clone()));
  //      fInitialTrackOutput->Add(initial_track);
        //cout<<initial_track->Position.Z()<<endl;
        //fDebugOutputCollector.fillContainer ("candidate_track", z_vert);
    }
    }
    initial_track = static_cast<Candidate*>((candidate_track->Clone()));
        fInitialTrackOutput->Add(initial_track);           
        //impact parameter      
        fItVertexInputArray->Reset();
        while (fItVertexInputArray->Next())   {
            imp_par = smeared_vertex->Position.Z() - candidate_track->Position.Z();
            //fDebugOutputCollector.fillContainer ("imp_par", imp_par);
                if(fabs(imp_par) < fabs(imp_par_tmp))  {
                    imp_par_tmp = imp_par;
                    candidate_track->VertexID_gen = smeared_vertex->VertexID_gen;
                    }         
                //fDebugOutputCollector.fillContainer ("imp_par", imp_par_tmp);  
                //fDebugOutputCollector.fillContainer ("z_at_vert", smeared_vertex->VertexID_gen);   
                //fVertexingTrackOutputArray->Add(candidate_track);  
                }
                //after calculating impact parameter, set track position at position in calorimeter
                /*
                candidate_track->Position.SetX(x_track);
                candidate_track->Position.SetY(y_track);
                candidate_track->Position.SetZ(z_track);
                candidate_track->Position.SetT(t_track);
                fDebugOutputCollector.fillContainer ("candidate_track", candidate_track->Position.Z());
                //fDebugOutputCollector.fillContainer ("imp_par", candidate_track->Position.Z());  
                fVertexingTrackOutputArray->Add(candidate_track);  */
            //}
            
        //}
        candidate_track->Position.SetX(x_track);
                candidate_track->Position.SetY(y_track);
                candidate_track->Position.SetZ(z_track);
                candidate_track->Position.SetT(t_track);
                //fDebugOutputCollector.fillContainer ("z_at_vert", candidate_track->Position.T()*inv_c_light); 
                //candidate_track->Position.Z());
                //fDebugOutputCollector.fillContainer ("imp_par", candidate_track->Position.Z());  
                //cout<<candidate_track->Position.T()<<endl;
                fVertexingTrackOutputArray->Add(candidate_track);  
    }
//cout<<lepton_energy<<"\t"<<leptons<<endl;
/*
double vertex_time, sum_weight = 0., weight = 0., track_time = 0.;
int track_molt = 0;
int vertices = 0;
int counter = 0;
double x_reco, y_reco, z_reco, z_at_cal, ct, TOF, t_gen, time_resolution;
TIterator *fItVertexingTrackOutputArray;
fItVertexingTrackOutputArray = fVertexingTrackOutputArray->MakeIterator();
fItVertexInputArray->Reset();
while((reco_vertex = static_cast<Candidate*>(fItVertexInputArray->Next()))) {
    ++ vertices;
    //fDebugOutputCollector.fillContainer ("z_at_vert", vertices);
    t_gen = reco_vertex->Position.T()*inv_c_light;
    //cout<<"\t"<<"nihuh"<<"\t"<<reco_vertex->Position.T()<<endl;
    int reco_vertexID = reco_vertex->VertexID_gen;
    //fDebugOutputCollector.fillContainer ("z_at_vert", reco_vertex->Position.T());
    fItVertexingTrackOutputArray->Reset();    
    while((track = static_cast<Candidate*>(fItVertexingTrackOutputArray->Next())))  {
    //fDebugOutputCollector.fillContainer ("z_at_vert", track->Position.T()*inv_c_light);
    if(track->Momentum.Pt()>0.5){ // && track->Momentum.P()>20.){
    //fDebugOutputCollector.fillContainer ("imp_par", track->Position.T()*inv_c_light);  
        if(reco_vertex->VertexID_gen == track->VertexID_gen) {
            
           //TOF
            z_at_cal = track->Position.Z();
            double dt = gRandom->Gaus(0.0, 30.);
            float Pt = track->Momentum.Pt();
                   float R = (Pt/(0.3*3.8));//Bz; m
                   //secant with respect to 0.0.0
                   float s = (TMath::Sqrt(TMath::Power(track->Position.X(),2)+TMath::Power(track->Position.Y(),2)));
                   float alpha = TMath::ASin(s/(1000*2*R));
                   float arc = 2 * alpha * R * 1000;
                   
                   float z_distance = z_at_cal - reco_vertex->Position.Z();
                   //TOF = (TMath::Sqrt((TMath::Power(z_at_cal,2) + TMath::Power(arc,2)))) *inv_c_light;
                   TOF = (TMath::Sqrt((TMath::Power(z_distance,2) + TMath::Power(arc,2)))) *inv_c_light;    
                   track_time = track->Position.T()*inv_c_light;
                   track_time += gRandom->Gaus(0.,30.);
                   ct = track_time - TOF;          
                   //ct = ((track->Position.T())*inv_c_light ) - TOF;
            //track->Position.SetT(ct);
            
            fDebugOutputCollector.fillContainer4D ("imp_par:Pt:P:PID", ct - t_gen,track->PID, track->Momentum.Pt(), track->Momentum.P());  
            
            if(ct<48000. && track->Momentum.Pt() > 2. && track->Momentum.P() > 5.){
                weight = TMath::Power(track->Momentum.Pt(),2);
                //weight =1.;
                sum_weight += weight;
            //vertex_time += track->Position.T();
            //vertex_time += ct * weight;
            vertex_time += track->Position.T();
            
            ++track_molt;
            //fDebugOutputCollector.fillContainer4D ("imp_par:Pt:P:PID", ct,track->PID, track->Momentum.Pt(), track->Momentum.P());  
            //track->Position.SetT(ct);
            //vertex_time += track->Position.T();
            }}
        }
        }   
    if(track_molt==0) ++counter;
    if(track_molt != 0){
    reco_vertex->Position.SetT(vertex_time/track_molt);
    //reco_vertex->Position.SetT(vertex_time/sum_weight);
    time_resolution = (reco_vertex->Position.T()-t_gen);///t_gen;
    //fDebugOutputCollector.fillContainer ("z_at_vert", t_gen);
    //if (track_molt != 0){
    fDebugOutputCollector.fillContainer3D ("candidate_track:track_molt:vertex_time", time_resolution, track_molt, reco_vertex->IsPU); }
    //fDebugOutputCollector.fillContainer2D ("imp_par:Pt", time_resolution, track->Momentum.Pt());  
    track_molt = 0;
    vertex_time = 0.;
    sum_weight = 0.;
    weight = 0.;

    //cout<<"vertex_time:"<<"\t"<<reco_vertex->Position.T()<<"\t"<<"time_resolution:"<<"\t"<<time_resolution<<"\t"<<t_gen<<endl;
}
    //cout<<vertices<<endl;

        int n = 0;
 //       float vertex_time = 0;
 vertex_time = 0.;
        fItVertexInputArray->Reset();
        while((vertex = static_cast<Candidate*>(fItVertexInputArray->Next())))  {
      
            fItTrackInputArray->Reset();
             while((candidate_track = static_cast<Candidate*>(fItTrackInputArray->Next())))  {
                //int n = 0;
                //float vertex_time = 0;
                if(candidate_track->VertexID_gen == vertex->VertexID_gen)   {
                    ++n;
                    float track_time = candidate_track->Position.T();
                    //cout<<track_time<<endl;
                    track_time *= inv_c_light;
                    vertex_time = track_time;
                    vertex_time += vertex_time;
                    //cout<<vertex_time<<endl;
                    }
        }
            vertex_time = vertex_time/n;
          // fDebugOutputCollector.fillContainer ("z_at_vert", vertex_time);       
            n = 0;  
            vertex_time = 0;
        }
         
        
    /*    double eff =1.- (double(counter)/double(vertices));
         fDebugOutputCollector.fillContainer ("z_at_vert", eff);                                                             
    cout<<"vert"<<""<<counter<<"\t"<<vertices<<"\t"<<eff<<endl; */
   
     
    





/*
for(int i=0; i<fVertexInputArray->GetEntries(); i++)
{
sumPt.push_back(0);
}
// loop over all input tracks arrays
track_iter.Reset();
vertex_iter.Reset();
while((candidate_track = static_cast<Candidate*>(track_iter.Next())))
{
int vertex_id_tmp = candidate_track->VertexID_gen;
Double_t pt_tmp = (candidate_track->Momentum).Pt();
if( abs(pt_tmp) > fpt_min )
{
sumPt.at(vertex_id_tmp) = sumPt.at(vertex_id_tmp) + pt_tmp*pt_tmp;
}
}
while((candidate_vertex = static_cast<Candidate*>(vertex_iter.Next())))
{
//fDebugOutputCollector.fillContainer ("z_at_vert", candidate_vertex->Position.T()*inv_c_light);
int vertex_id_tmp = candidate_vertex->VertexID_gen;
candidate_vertex->sumPtSquare = sumPt.at(vertex_id_tmp);
}

double tmp_sumpt2 = 0.;
int PV = 0;
double pvcounter = 0.;
double pveff;
double vert_counter = 0.;
vertex_iter.Reset();
while((vertex = static_cast<Candidate*>(vertex_iter.Next())))   {
    ++vert_counter;
    if(vertex->sumPtSquare > tmp_sumpt2)    {
        tmp_sumpt2 = vertex->sumPtSquare;
        PV = vertex->IsPU;
        if(PV != 0) {
            PV =1;
            ++pvcounter;}
        }
    }
    //cout<<"IsPV:"<<"\t"<<PV<<endl;
    fDebugOutputCollector.fillContainer ("z_at_vert", (1.- pvcounter/vert_counter));
 */   
}
//------------------------------------------------------------------------------

