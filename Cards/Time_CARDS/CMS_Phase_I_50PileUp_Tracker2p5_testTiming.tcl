#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  PileUpMerger 

  ModifyBeamSpot

  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronEnergySmearing
  MuonMomentumSmearing

  TrackMerger
  Vertexing
  
  ECalorimeter
  HCalorimeter
  
  TrackMergerWithMuon
  TrackPileUpSubtractor
  
  
  
  
  
  
  EFlowMerger
  PvEFlowMerger
  notimeEFlowMerger

  RhoKt4
  PvRhoKt4

  noTimeTowersRhoKt4
  2PvRhoKt4

  FastJetFinder
  TrackJetFinder
  
  PvFastJetFinder
  PvTrackJetFinder
  

  NeutrinoFilter
  GenJetFinderNoNu

  RhoMerger

  JetPileUpSubtractor
  JetFlavourAssociation
  
   
  
  PvJetPileUpSubtractor
  

  BTagging

  
  
  PileUpJetID
  PvPileUpJetID
  

  PhotonEfficiency
  PhotonIsolation 

  ElectronEfficiency 
  ElectronIsolation 
 
  MuonEfficiency
  MuonIsolation  

  GenMissingET
  MissingET

  GenScalarHT
  ScalarHT
  
  TreeWriter 
}

### remove some modules

###kenzo
# GlobalRhoKt4
 # GlobalRhoGridFastJet
  
#  PvGlobalRhoKt4
 # PvGlobalRhoGridFastJet
#RhoGridFastJet
#PvRhoGridFastJet

# GenBeamSpotFilter
# StatusPid
# JetPileUpSubtractorGrid
# JetPileUpSubtractor4VArea
# PuppiJetPileUpSubtractorGrid
# PuppiJetPileUpSubtractor4VArea

#### remove the module which do the filter of jet constituent
# ConstituentFilter
# PuppiConstituentFilter

#### removing all puppi modules for the time being

#  RunPUPPI
#  PuppiRhoKt4
#  PuppiRhoGrid
#  PuppiJetFinder

#  PuppiJetPileUpSubtractor
#  PuppiJetFlavourAssociation

#  PuppiBTagging

#  PuppiPileUpJetID
#  PuppiMissingET
#  PuppiScalarHT

#  TreeWriter


#################
# PileUp Merger #
#################

module PileUpMerger PileUpMerger {
 ## inputs are status 1 HEPMC particles --> real final state one
 set InputArray  Delphes/stableParticles
 ## output array is called stable particles
 set OutputArray stableParticles
 ## store NPU
 set NPUOutputArray NPU
 ## output vertices array
 set VertexOutputArray vertices
 # Get rid of beam spot from http://red-gridftp11.unl.edu/Snowmass/MinBias100K_14TeV.pileup ...
# set InputBSX 2.44
# set InputBSY 3.39
 # replace it with beam spot from CMSSW files  
# set OutputBSX 0.24
# set OutputBSY 0.39  

#in mm
 set OutputBSX 0.
 set OutputBSY 0.
 set OutputBSZ 0.

#in ps
 set TimeDelay 0.
 
  # pre-generated minbias input file
 #set PileUpFile /afs/cern.ch/user/g/govoni/work/TP/Delphes/testSample/MB_1.mb
 set PileUpFile /afs/cern.ch/user/g/govoni/work/public/4Vincenzo/MB_1.mb
 #set PileUpFile  /afs/cern.ch/work/s/spigazzi/public/PU_14TeV/MinBias_14TeV_100k_TunePP15.pileup
 #average expected pile up (poissonian generation)
 set MeanPileUp 0
 # spread in the beam direction in mm (assumed gaussian) ; 
 set ZVertexSpread 0.
 # time spread in the events happening in ps (assumed gaussian) ; 
 # baseline = 160 ps, crab-kiss = 100 ps
 set TVertexSpread 160.
 
 # set TVertexSpread 100
 #set simpleOutputFileName /afs/cern.ch/user/g/govoni/work/TP/Delphes/testSample/simpleOutput_PU.root
 set simpleOutputFileName /afs/cern.ch/work/l/lbrianza/VINCENZO/CMSSW_6_2_0_SLHC23/src/Delphes/simpleOutput_PU.root
}

##################
# ModifyBeamSpot #
##################


module ModifyBeamSpot ModifyBeamSpot {
  # in ps
  set TVertexSpread 160.
  # in mm
  set ZVertexSpread 0.

  set OutputBSX 0.
  set OutputBSY 0.  
  set OutputBSZ 0.
  set OutputBST 0.

  ##set InputArray    ParticlePropagator/stableParticles 
  set InputArray    PileUpMerger/stableParticles
  set VertexInputArray PileUpMerger/vertices
  set OutputArray   stableParticles
  set VertexOutputArray vertices
  set PVOutputArray PV  

  #set simpleOutputFileName /afs/cern.ch/user/g/govoni/work/TP/Delphes/testSample/simpleOutput_MB.root
  set simpleOutputFileName /afs/cern.ch/work/l/lbrianza/VINCENZO/CMSSW_6_2_0_SLHC23/src/Delphes/simpleOutput_MB.root
}


#####################################################################################
# Propagate particles in cylinder and divide charged particles in different classes #
#####################################################################################

module ParticlePropagator ParticlePropagator {
  ## take particles after beam spot smearing
  set InputArray ModifyBeamSpot/stableParticles
  ## produce independent output collection: all particles, only charged hadrons, electrons and muons
  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons
  ## radius of the magnetic field coverage, in m
  set Radius 1.29
  ## half-length of the magnetic field coverage, in m
  set HalfLength 3.00
  ## magnetic field
  set Bz 3.8
}

###############################################################################################################
# StatusPidFilter: this module removes all generated particles except electrons, muons, taus, and status == 3 #
###############################################################################################################

module StatusPidFilter StatusPid {
    ## take the particles from Pythia8 not adding pile-up
    set InputArray  Delphes/allParticles
    set OutputArray filteredParticles
    set PTMin 0.35
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  ## particles after propagation
  set InputArray  ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons
  # tracking efficiency formula for charged hadrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.85) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.97) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.85) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.90) + \
                         (abs(eta) > 2.5)                                                  * (0.00)
  }
}


#####################################
# Electron tracking efficiency - ID
####################################

module Efficiency ElectronTrackingEfficiency {
  set InputArray  ParticlePropagator/electrons
  set OutputArray electrons
  # tracking efficiency formula for electrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.85) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e2) * (0.97) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e2)                * (0.99) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.85) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e2) * (0.90) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e2)                * (0.95) + \
                         (abs(eta) > 2.5)                                                  * (0.00)
  }
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons
  # tracking efficiency formula for muons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.998) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.9998) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.98) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.98) + \
                         (abs(eta) > 2.5)                                                  * (0.00)
  }
}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  ## hadrons after having applied the tracking efficiency
  set InputArray  ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons
  # resolution formula for charged hadrons
  set ResolutionFormula {                  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.013) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.02) + \
                                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta)>2.5)*(0.00)
  }
}

#################################
# Energy resolution for electrons
#################################

module EnergySmearing ElectronEnergySmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons
  # set ResolutionFormula {resolution formula as a function of eta and energy}
  set ResolutionFormula {  (abs(eta) <= 2.5) * (energy > 0.1   && energy <= 2.5e1) * (energy*0.025) + \
                           (abs(eta) <= 2.5) * (energy > 2.5e1)                    * (energy*0.035) + \
                           (abs(eta) > 2.5 && abs(eta) <= 3.0)                     * (energy*0.035) + \
                           (abs(eta) > 3.0 && abs(eta) <= 5.0)                     * (energy*0.07)
  }                                                  
}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons
  # resolution formula for muons
  set ResolutionFormula {  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.012) + \
                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.015) + \
                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.03) + \
                           (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                           (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.015) + \
                           (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.025) + \
                           (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.03)
  }
}

#################################################################
# Track merger : merge two collection of object in a output one #
#################################################################

module Merger TrackMerger {
  ## take smeared charged hadron and electrons as only tracks 
  ## to take into account in the calorimeter simulation
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray MuonMomentumSmearing/muons
  add InputArray ElectronEnergySmearing/electrons
  set OutputArray tracks
}

##############################################################
# Vertex Reco
##############################################################
module Vertexing Vertexing {
set TrackInputArray TrackMerger/tracks
##set TrackInputArray ModifyBeamSpot/stableParticles
##set VertexInputArray PileUpMerger/vertices
set VertexInputArray ModifyBeamSpot/vertices
##set GenVertexInputArray ModifyBeamSpot/stableParticles

set VertexingOutputArray vertexingtracks
set InitialTrackOutput  initialtrack
set pt_min 0.0
#set pt_min 0.5
}


##############################################################
# ECalorimeter : emulate electromagnetic calorimeter answer making caloTowers #
##############################################################

module ECalorimeter ECalorimeter {
  ## particle from the propagation without any efficiency or smearing (for neutrals)
  set ParticleInputArray    ParticlePropagator/stableParticles
  ## track after smearing and efficiency: used for charged particles
  ##set VerticesInputArray    PileUpMerger/vertices
  set VerticesInputArray    ModifyBeamSpot/vertices
  set TrackInputArray       TrackMerger/tracks
  set LHEPartonInputArray   Delphes/LHEParticles
  set VertexingTrackArray   Vertexing/vertexingtracks
  set PVInputArray          ModifyBeamSpot/PV

  ## output collections
  set TowerOutputArray        ecalTowers
  set PhotonOutputArray       photons
  set PvEFlowTrackOutputArray ecalPvEflowTracks
  set PvEFlowTowerOutputArray ecalPvEflowTowers
  set EFlowTowerOutputArray   ecalEflowTowers
  set EFlowTrackOutputArray   ecalEflowTracks
  set notimeTowerOutputArray  notimeEcalEflowTowers
  set notimeTrackOutputArray  notimeEcalEflowTracks
                              
  
  # flat ntuples output
  #set simpleOutputFileName /afs/cern.ch/user/g/govoni/work/TP/Delphes/testSample/simpleOutput_Ca.root
  set simpleOutputFileName /afs/cern.ch/work/l/lbrianza/VINCENZO/CMSSW_6_2_0_SLHC23/src/Delphes/simpleOutput_ECal.root

  # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
  #   timing of the calo deposits
  # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
  
  # Minimum threshold for a particle deposit in a tower
  # to enter in the computation of the timing associated to a tower
  # in GeV.
  set TimingEMin  0.
  set TimeCut 100.

  # profile of the delay for a particle traveling from the interaction vertex
  # to the calorimetry, determined from LHE particles in Delphes
  # originating from (0,0,0)
  # up to eta = 1.6, delay in millimeters
  set DelayBarrel_0 1287.51    
  set DelayBarrel_1 5.55276    
  set DelayBarrel_2 635.669    
  set DelayBarrel_3 19.2171    
  set DelayBarrel_4 30.0051    
  set DelayBarrel_5 12.1216    

  # from eta = 1.6, delay in millimeters
  set DelayEndcap_0  6.41347e+03
  set DelayEndcap_1 -4.49583e+03
  set DelayEndcap_2  2.40551e+03
  set DelayEndcap_3 -6.47156e+02
  set DelayEndcap_4  8.70166e+01
  set DelayEndcap_5 -4.66147e+00

  # up to eta = 1.6, delay in picoseconds
#  set DelayBarrel_0   4294.67
#  set DelayBarrel_1     18.52
#  set DelayBarrel_2   2120.36
#  set DelayBarrel_3     64.10
#  set DelayBarrel_4    100.09
#  set DelayBarrel_5     40.43

  # from eta = 1.6, delay in picoseconds
#  set DelayEndcap_0  21393.03
#  set DelayEndcap_1 -14996.47
#  set DelayEndcap_2   8023.92
#  set DelayEndcap_3  -2158.68
#  set DelayEndcap_4    290.26
#  set DelayEndcap_5    -15.55


  # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
  #   calorimetry granularity
  # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  set pi [expr {acos(-1)}]

  ## Granularity for |eta| < 1.65
  set PhiBins {}
  ##for {set i -36} {$i <= 36} {incr i} 
  for {set i -180} {$i <= 180} {incr i} {
    ##add PhiBins [expr {$i * $pi/36.0}]
    add PhiBins [expr {$i * $pi/180.0}]
  }
  ##foreach eta -1.566 -1.479 -1.392 -1.305 -1.218 -1.131 -1.044 -0.957 -0.87 -0.783 -0.696 -0.609 -0.522 -0.435 -0.348 -0.261 -0.174 -0.087 0 0.087 0.174 0.261 0.348 0.435 0.522 0.609 0.696 0.783 0.87 0.957 1.044 1.131 1.218 1.305 1.392 1.479 1.566 1.653
 
foreach eta {-1.479 1.461 -1.444 -1.426 -1.409 -1.392 -1.374 -1.356 -1.339 -1.321 -1.305 -1.286 -1.269 -1.251 -1.235 -1.218 -1.200 -1.183 -1.165 -1.148 -1.131 -1.113 -1.096 -1.078 -1.061 -1.044 -1.026 -1.009 -0.991 -0.973  -0.957 -0.939 -0.922 -0.904  -0.887  -0.87 -0.852 -0.835 -0.817 -0.8 -0.783 -0.765 -0.748 -0.730 -0.713 -0.696 -0.678 -0.661 -0.643 -0.626 -0.609 -0.591 -0.574 -0.556 -0.539 -0.522 -0.504 -0.487 -0.469 -0.452  -0.435 -0.417 -0.4 -0.382 -0.365 -0.348 -0.33 -0.313 -0.295 -0.278 -0.261 -0.243 -0.226 -0.208 -0.191 -0.174 -0.156 -0.139 -0.121 -0.103 -0.087 -0.069 -0.052 -0.034 -0.017 0 0.017 0.034 0.052 0.069 0.087 0.103 0.121 0.139 0.156 0.174 0.191 0.208 0.226 0.243 0.261 0.278 0.295 0.313 0.33 0.348 0.365 0.382 0.4 0.417 0.435 0.452 0.469 0.487 0.504 0.522 0.539 0.556 0.574 0.591 0.609 0.626 0.643 0.661 0.678 0.696 0.713 0.730 0.748 0.765 0.783 0.8 0.817 0.835 0.852 0.87 0.887 0.904 0.922 0.939 0.957 0.973 0.991 1.009 1.026 1.044 1.061 1.078 1.096 1.113 1.131 1.148 1.165 1.183 1.2 1.218 1.235 1.251 1.269 1.286 1.305 1.321 1.339 1.356 1.374 1.392 1.409 1.426 1.444 1.461 1.479}   {
    add EtaPhiBins $eta $PhiBins
  }

  ## Granularity for 1.65 < |eta| < 4.5
  ## Granularity for 1.479 < |eta| < 3   ecal endcap coverage
  set PhiBins {}
  ##for {set i -18} {$i <= 18} {incr i} 
  for {set i -90} {$i <= 90} {incr i} {
    ##add PhiBins [expr {$i * $pi/18.0}]
    add PhiBins [expr {$i * $pi/90.0}]
  }
    foreach eta {-3 -2.95 -2.933 -2.883 -2.833 -2.783 -2.733 -2.683 -2.639 -2.596 -2.553 -2.51 -2.46 -2.41 -2.36 -2.31 -2.267 -2.224 -2.182 -2.145 -2.108 -2.071 -2.034 -2.002 -1.97 -1.938 -1.921 -1.904 -1.886 -1.868 -1.851 -1.833 -1.816 -1.798 -1.781 -1.763 -1.746 -1.728 -1.711 -1.693 -1.676  
   -1.658 -1.641 -1.623 -1.606 -1.588 -1.571 -1.553 -1.536 -1.518 -1.501 -1.483  1.483 1.501 1.518 1.536 1.553 1.571 1.588 1.606 1.623 1.641 1.658 1.676 1.693 1.711 1.728 1.746 1.763 1.781 1.798 1.816 1.833 1.851 1.868 1.886 1.904 1.921 1.938 1.97 2.002 2.034 2.071 2.108 2.145 2.182 2.224 2.267 2.31 2.36 2.41 2.46 2.51 2.553 2.596 2.639 2.683 2.733 2.783 2.833 2.883 2.933 2.95 3}    {
   add EtaPhiBins $eta $PhiBins
  }

  ## Granularity for 4.5 < |eta| < 5
  ##set PhiBins {}
  ##for {set i -9} {$i <= 9} {incr i} {
  ##  add PhiBins [expr {$i * $pi/9.0}]
  ##} 
  ##foreach eta {-5 -4.7 -4.525 4.7 5} {
  ##  add EtaPhiBins $eta $PhiBins
  ##}

  ### energy deposition for each particle type
  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
  add EnergyFraction {0} {0.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11}  {1.0}
  add EnergyFraction {22}  {1.0}
  add EnergyFraction {111} {1.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.3}
  add EnergyFraction {3122} {0.3}

  # set ECalResolutionFormula {resolution formula as a function of eta and energy}
  set ECalResolutionFormula {                  (abs(eta) <= 3.0) * sqrt(energy^2*0.005^2 + energy*0.027^2 + 0.15^2) + \
                             (abs(eta) > 3.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.08^2 + energy*1.97^2)
  }

}

##############################################################
# HCalorimeter : emulate hadronic calorimeter answer making caloTowers #
##############################################################

module HCalorimeter HCalorimeter {
  ## particle from the propagation without any efficiency or smearing (for neutrals)
  set ParticleInputArray    ParticlePropagator/stableParticles
  ## track after smearing and efficiency: used for charged particles
  ##set VerticesInputArray    PileUpMerger/vertices
  set VerticesInputArray    ModifyBeamSpot/vertices
  set TrackInputArray       TrackMerger/tracks
  ##set TrackInputArray       ECalorimeter/ecalEflowTracks
  set LHEPartonInputArray   Delphes/LHEParticles
  set VertexingTrackArray   Vertexing/vertexingtracks

  ## output collections
  set TowerOutputArray      hcalTowers
  set PhotonOutputArray     photons
  set EFlowTrackOutputArray hcalEflowTracks
  set EFlowTowerOutputArray hcalEflowTowers
  
  # flat ntuples output
  #set simpleOutputFileName /afs/cern.ch/user/g/govoni/work/TP/Delphes/testSample/simpleOutput_Ca.root
  set simpleOutputFileName /afs/cern.ch/work/l/lbrianza/VINCENZO/CMSSW_6_2_0_SLHC23/src/Delphes/simpleOutput_HCa.root

  # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
  #   timing of the calo deposits
  # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
  
  # Minimum threshold for a particle deposit in a tower
  # to enter in the computation of the timing associated to a tower
  # in GeV.
  set TimingEMin 0.0

  # profile of the delay for a particle traveling from the interaction vertex
  # to the calorimetry, determined from LHE particles in Delphes
  # originating from (0,0,0)
  # up to eta = 1.6, delay in millimeters
  set DelayBarrel_0 1287.51    
  set DelayBarrel_1 5.55276    
  set DelayBarrel_2 635.669    
  set DelayBarrel_3 19.2171    
  set DelayBarrel_4 30.0051    
  set DelayBarrel_5 12.1216    

  # from eta = 1.6, delay in millimeters
  set DelayEndcap_0  6.41347e+03
  set DelayEndcap_1 -4.49583e+03
  set DelayEndcap_2  2.40551e+03
  set DelayEndcap_3 -6.47156e+02
  set DelayEndcap_4  8.70166e+01
  set DelayEndcap_5 -4.66147e+00

  # up to eta = 1.6, delay in picoseconds
#  set DelayBarrel_0   4294.67
#  set DelayBarrel_1     18.52
#  set DelayBarrel_2   2120.36
#  set DelayBarrel_3     64.10
#  set DelayBarrel_4    100.09
#  set DelayBarrel_5     40.43

  # from eta = 1.6, delay in picoseconds
#  set DelayEndcap_0  21393.03
#  set DelayEndcap_1 -14996.47
#  set DelayEndcap_2   8023.92
#  set DelayEndcap_3  -2158.68
#  set DelayEndcap_4    290.26
#  set DelayEndcap_5    -15.55


  # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
  #   calorimetry granularity
  # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  set pi [expr {acos(-1)}]

  ## Granularity for |eta| < 1.65
  set PhiBins {}
  for {set i -36} {$i <= 36} {incr i} {
    add PhiBins [expr {$i * $pi/36.0}]
  }
  foreach eta {-1.566 -1.479 -1.392 -1.305 -1.218 -1.131 -1.044 -0.957 -0.87 -0.783 -0.696 -0.609 -0.522 -0.435 -0.348 -0.261 -0.174 -0.087 0 0.087 0.174 0.261 0.348 0.435 0.522 0.609 0.696 0.783 0.87 0.957 1.044 1.131 1.218 1.305 1.392 1.479 1.566 1.653} {
    add EtaPhiBins $eta $PhiBins
  }

  ## Granularity for 1.65 < |eta| < 4.5
  set PhiBins {}
  for {set i -18} {$i <= 18} {incr i} {
    add PhiBins [expr {$i * $pi/18.0}]
  }
  foreach eta {-4.35 -4.175 -4 -3.825 -3.65 -3.475 -3.3 -3.125 -2.95 -2.868 -2.65 -2.5 -2.322 -2.172 -2.043 -1.93 -1.83 -1.74 -1.653 1.74 1.83 1.93 2.043 2.172 2.322 2.5 2.65 2.868 2.95 3.125 3.3 3.475 3.65 3.825 4 4.175 4.35 4.525} {
    add EtaPhiBins $eta $PhiBins
  }

  ## Granularity for 4.5 < |eta| < 5
  set PhiBins {}
  for {set i -9} {$i <= 9} {incr i} {
    add PhiBins [expr {$i * $pi/9.0}]
  }
  foreach eta {-5 -4.7 -4.525 4.7 5} {
    add EtaPhiBins $eta $PhiBins
  }

  ### energy deposition for each particle type
  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
  add EnergyFraction {0} {1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11}  {0.0}
  add EnergyFraction {22}  {0.0}
  add EnergyFraction {111} {0.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.7}
  add EnergyFraction {3122} {0.7}

  # set HCalResolutionFormula {resolution formula as a function of eta and energy}
  set HCalResolutionFormula {                  (abs(eta) <= 1.7) * sqrt(energy^2*0.0302^2 + energy*0.5205^2 + 1.59^2) + \
                             (abs(eta) > 1.7 && abs(eta) <= 3.2) * sqrt(energy^2*0.050^2 + energy*0.706^2) + \
                             (abs(eta) > 3.0 && abs(eta) <= 4.9) * sqrt(energy^2*0.05^2 + energy*1.00^2)
  }


}

####################################################################
## Track pile-up subtractor: apply CHS on top of track collection ##
####################################################################

module TrackPileUpSubtractor TrackPileUpSubtractor {
  ## take tracks from calorimeter module, smeared electrons and smeared muon. Take the PV from the ModifyBeamSpot --> pileup subtraction or CHS
  ## Is not useful to run this module on top of NoPU collections
  #add InputArray ECalorimeter/ecalEflowTracks    ecalEflowTracks
  ##add InputArray HCalorimeter/hcalEflowTracks    hcalEflowTracks
  add InputArray ElectronEnergySmearing/electrons electrons
  add InputArray MuonMomentumSmearing/muons muons
  add InputArray Vertexing/initialtrack initialtrack
  
  set PVInputArray  ModifyBeamSpot/PV 
  # assume perfect pile-up subtraction for tracks with |z| > fZVertexResolution in mm
  #set ZVertexResolution 0.0002
  #set ZVertexResolution 0.02
  set ZVertexResolution 0.2

}

####################################################################
## PvTrack pile-up subtractor: apply CHS on top of track collection ##
####################################################################

module TrackPileUpSubtractor PvTrackPileUpSubtractor {
  ## take tracks from calorimeter module, smeared electrons and smeared muon. Take the PV from the ModifyBeamSpot --> pileup subtraction or CHS
  ## Is not useful to run this module on top of NoPU collections
  add InputArray ECalorimeter/ecalPvEflowTracks    ecalPvEflowTracks
  #add InputArray ECalorimeter/ecalEflowTracks    ecalEflowTracks
  ##add InputArray/TrackMerger tracks
  ##add InputArray HCalorimeter/hcalEflowTracks    hcalEflowTracks
  #add InputArray ElectronEnergySmearing/electrons electrons
  #add InputArray MuonMomentumSmearing/muons muons
  #add InputArray Vertexing/initialtrack initialtrack
  
  set PVInputArray  ModifyBeamSpot/PV 
  # assume perfect pile-up subtraction for tracks with |z| > fZVertexResolution in mm
  set ZVertexResolution 0.02

}


####################################################################
## notimeTrack pile-up subtractor: apply CHS on top of track collection ##
####################################################################

module TrackPileUpSubtractor notimeTrackPileUpSubtractor {
  ## take tracks from calorimeter module, smeared electrons and smeared muon. Take the PV from the ModifyBeamSpot --> pileup subtraction or CHS
  ## Is not useful to run this module on top of NoPU collections
  #add InputArray ECalorimeter/ecalPvEflowTracks    ecalPvEflowTracks
  add InputArray ECalorimeter/notimeEcalEflowTracks    notimeEcalEflowTracks
  ##add InputArray/TrackMerger tracks
  ##add InputArray HCalorimeter/hcalEflowTracks    hcalEflowTracks
  #add InputArray ElectronEnergySmearing/electrons electrons
  #add InputArray MuonMomentumSmearing/muons muons
  add InputArray Vertexing/initialtrack initialtrack
  
  set PVInputArray  ModifyBeamSpot/PV 
  # assume perfect pile-up subtraction for tracks with |z| > fZVertexResolution in m
  set ZVertexResolution 0.00002

}


########################
## Energy flow merger ##
########################

module Merger EFlowMerger {
  ## charged particles after having applied CHS
  #add InputArray TrackPileUpSubtractor/ecalEflowTracks
  add InputArray TrackPileUpSubtractor/initialtrack
  #kenzo hcaleflowtracka
  ##add InputArray TrackPileUpSubtractor/hcalEflowTracks
  ## calorimeter towers to get also photons and neutral hadrons for the whole detector coverage
  add InputArray ECalorimeter/ecalEflowTowers
  add InputArray HCalorimeter/hcalEflowTowers
  ##add InputArray MuonMomentumSmearing/muons

  set OutputArray eflow
}

###########################
## Pv Energy flow merger ##
###########################

module Merger PvEFlowMerger {
  ## charged particles after having applied CHS
  #add InputArray PvTrackPileUpSubtractor/ecalPvEflowTracks
  #add InputArray PvTrackPileUpSubtractor/ecalEflowTracks
  add InputArray TrackPileUpSubtractor/initialtrack
  #add InputArray PvTrackPileUpSubtractor/tracks    
  #add InputArray ECalorimeter/ecalPvEflowTracks
  # add kenzo
  #add InputArray HCalorimeter/hcalEflowTracks
  ##add InputArray PvTrackPileUpSubtractor/hcalEflowTracks
  ## calorimeter towers to get also photons and neutral hadrons for the whole detector coverage
  add InputArray ECalorimeter/ecalPvEflowTowers
  add InputArray HCalorimeter/hcalEflowTowers
  ##add InputArray MuonMomentumSmearing/muons
  #add InputArray ECalorimeter/notimeEcalEflowTowers
  set OutputArray Pveflow
}

################################
## No time towers flow merger ##
################################


module Merger notimeEFlowMerger {
  ## charged particles after having applied CHS
  #add InputArray TrackPileUpSubtractor/ecalEflowTracks
  #add InputArray TrackPileUpSubtractor/initialtrack
  #add InputArray TrackMerger/tracks
  #add InputArray notimeTrackPileUpSubtractor/notimeEcalEflowTracks
  #kenzo hcaleflowtracka
  ##add InputArray TrackPileUpSubtractor/hcalEflowTracks
  ## calorimeter towers to get also photons and neutral hadrons for the whole detector coverage
  add InputArray ECalorimeter/notimeEcalEflowTowers
  add InputArray HCalorimeter/hcalEflowTowers
  ##add InputArray MuonMomentumSmearing/muons

  set OutputArray eflow
}



##############################################
## Calculate Rho using KtClustering or grid ##
##############################################

module FastJetFinder GlobalRhoKt4 {
  ## take as input the particle flow particles
  set InputArray EFlowMerger/eflow
  ## output name
  set RhoOutputArray rho
  ## compute rho clustering the event  
  set ComputeRho     true
  ## not compute rho using the grid
  set ComputeRhoGrid false
  ## area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1 
  ## jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 4
  ## Clustering and Ghost parameter
  set ParameterR  0.4
  set GhostEtaMax 5.0
  set RhoEtaMax   5.0
  ## eta bins for rho evaluation
  add RhoEtaRange 0.0 5.0  
  set JetPTMin 0.0
}

module FastJetFinder GlobalRhoGridFastJet {
  ## take as input the particle flow particles
  set InputArray EFlowMerger/eflow
  ## compute rho clustering the event  
  set ComputeRho     false
  ## not compute rho using the grid
  set ComputeRhoGrid true
  set RhoOutputArray rho  
  ## area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1  
  ## jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 4
  ## Clustering and Ghost parameter
  set ParameterR  0.4
  set GhostEtaMax 5.0
  set RhoEtaMax   5.0
  ## eta bins for rho evaluation
  add RhoEtaRange 0.0 5.0  
  set JetPTMin 0.0
}

####################################################
## Calculate Rho using KtClustering or grid for Pv##
####################################################

module FastJetFinder PvGlobalRhoKt4 {
  ## take as input the particle flow particles
  set InputArray PvEFlowMerger/Pveflow
  ## output name
  set RhoOutputArray rho
  ## compute rho clustering the event  
  set ComputeRho     true
  ## not compute rho using the grid
  set ComputeRhoGrid false
  ## area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1 
  ## jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 4
  ## Clustering and Ghost parameter
  set ParameterR  0.4
  set GhostEtaMax 5.0
  set RhoEtaMax   5.0
  ## eta bins for rho evaluation
  add RhoEtaRange 0.0 5.0  
  set JetPTMin 0.0
}

module FastJetFinder PvGlobalRhoGridFastJet {
  ## take as input the particle flow particles
  set InputArray PvEFlowMerger/Pveflow
  ## compute rho clustering the event  
  set ComputeRho     false
  ## not compute rho using the grid
  set ComputeRhoGrid true
  set RhoOutputArray rho  
  ## area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1  
  ## jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 4
  ## Clustering and Ghost parameter
  set ParameterR  0.4
  set GhostEtaMax 5.0
  set RhoEtaMax   5.0
  ## eta bins for rho evaluation
  add RhoEtaRange 0.0 5.0  
  set JetPTMin 0.0
}



#######################################
## Rho pile-up in different eta bins ##
#######################################

module FastJetFinder RhoKt4 {
  # input particles
  set InputArray EFlowMerger/eflow 
  # output name
  set RhoOutputArray rho
  ## compute rho clustering the event  
  set ComputeRho     true
  ## not compute rho using the grid
  set ComputeRhoGrid false
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 4
  set ParameterR   0.4
  set GhostEtaMax  5.0
  set RhoEtaMax    5.0
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0
  set JetPTMin 0.0
}

module FastJetFinder RhoGridFastJet {
  ## take as input the particle flow particles
  set InputArray EFlowMerger/eflow
  ## output name
  set RhoOutputArray rho
  ## compute rho clustering the event  
  set ComputeRho     false
  ## not compute rho using the grid
  set ComputeRhoGrid true
  ## area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1  
  ## Clustering and Ghost parameter
  set GhostEtaMax 5.0
  set RhoEtaMax   5.0
  ## eta bins for rho evaluation
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0  
  set JetPTMin 0.0
}

#######################################
## Rho pile-up in different eta bins ##
#######################################

module FastJetFinder PvRhoKt4 {
  # input particles
  set InputArray PvEFlowMerger/Pveflow
  #set InputArray ECalorimeter/ecalPvEflowTowers
  #set InputArray HCalorimeter/hcalEflowTowers
  #set InputArray notimeEFlowMerger/eflow
  # output name
  set RhoOutputArray rho
  ## compute rho clustering the event  
  set ComputeRho     true
  #set ComputeRho     false
  ## not compute rho using the grid
  set ComputeRhoGrid false
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 4
  set ParameterR   0.4
  set GhostEtaMax  5.0
  set RhoEtaMax    5.0
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0
  set JetPTMin 0.0
}

#######################################
## 2 Rho pile-up in different eta bins ##
#######################################

module FastJetFinder 2PvRhoKt4 {
  # input particles
  #set InputArray PvEFlowMerger/Pveflow
  set InputArray ECalorimeter/ecalPvEflowTowers
  #set InputArray HCalorimeter/hcalEflowTowers    
  # output name
  set RhoOutputArray 2rho
  ## compute rho clustering the event  
  set ComputeRho     true
  #set ComputeRho     false
  ## not compute rho using the grid
  set ComputeRhoGrid false
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 4
  set ParameterR   0.4
  set GhostEtaMax  5.0
  set RhoEtaMax    5.0
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0
  set JetPTMin 0.0
}



module FastJetFinder PvRhoGridFastJet {
  ## take as input the particle flow particles
  set InputArray PvEFlowMerger/Pveflow
  ## output name
  set RhoOutputArray rho
  ## compute rho clustering the event  
  set ComputeRho     false
  ## not compute rho using the grid
  set ComputeRhoGrid true
  ## area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1  
  ## Clustering and Ghost parameter
  set GhostEtaMax 5.0
  set RhoEtaMax   5.0
  ## eta bins for rho evaluation
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0  
  set JetPTMin 0.0
}

###########################################################
## Rho computation for Ecal towers with no assigned time ##
###########################################################


module FastJetFinder noTimeTowersRhoKt4 {
  # input particles
  set InputArray notimeEFlowMerger/eflow
  #set InputArray ECalorimeter/notimeEcalEflowTowers  
  #set InputArray notimeEFlowMerger/eflow
  # output name
  set RhoOutputArray rho
  ## compute rho clustering the event  
  set ComputeRho     true
  ## not compute rho using the grid
  set ComputeRhoGrid false
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 4
  set ParameterR   0.4
  set GhostEtaMax  5.0
  set RhoEtaMax    5.0
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0
  set JetPTMin 0.0
}





################
## Jet finder ##
################

module FastJetFinder FastJetFinder {
  set InputArray   EFlowMerger/eflow
  set OutputArray  jets
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm  6
  set ParameterR    0.4
  set JetPTMin      20.0
  #set simpleOutputFileName /afs/cern.ch/user/g/govoni/work/TP/Delphes/testSample/simpleOutput_FJ.root
  set simpleOutputFileName /afs/cern.ch/work/l/lbrianza/VINCENZO/CMSSW_6_2_0_SLHC23/src/Delphes/simpleOutput_FJ.root
}

################
## Jet finder ##
################

module FastJetFinder PvFastJetFinder {
  set InputArray   PvEFlowMerger/Pveflow
  set OutputArray  jets
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm  6
  set ParameterR    0.4
  set JetPTMin      20.0
  #set simpleOutputFileName /afs/cern.ch/user/g/govoni/work/TP/Delphes/testSample/simpleOutput_FJ.root
  set simpleOutputFileName /afs/cern.ch/work/l/lbrianza/VINCENZO/CMSSW_6_2_0_SLHC23/src/Delphes/simpleOutput_FJ_pv.root
}

##################################
## Jet finder in no time towers ##
##################################

module FastJetFinder noTimeTowersFastJetFinder {
  #set InputArray   PvEFlowMerger/Pveflow
  #set InputArray   ECalorimeter/notimeEcalEflowTowers
  set InputArray   notimeEFlowMerger/eflow
  set OutputArray  jets
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm  6
  set ParameterR    0.4
  set JetPTMin      20.0
  #set simpleOutputFileName /afs/cern.ch/user/g/govoni/work/TP/Delphes/testSample/simpleOutput_FJ.root
  set simpleOutputFileName /afs/cern.ch/work/l/lbrianza/VINCENZO/CMSSW_6_2_0_SLHC23/src/Delphes/simpleOutput_FJ.root
}







####################
### Track jets ###FJ
####################

module FastJetFinder TrackJetFinder {
  #set InputArray  TrackPileUpSubtractor/ecalEflowTracks
  set InputArray  Vertexing/initialtrack
  set OutputArray jets
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR   0.4
  set JetPTMin      1.0
  set ParticlePTMin 0.35
}

#######################
### Pv Track jets ###FJ
#######################

module FastJetFinder PvTrackJetFinder {
  #set InputArray  TrackPileUpSubtractor/ecalEflowTracks
  set InputArray  Vertexing/initialtrack
  set OutputArray jets
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR   0.4
  set JetPTMin      1.0
  set ParticlePTMin 0.35
}

#######################
# MC truth jet finder #
#######################

module FastJetFinder GenJetFinder { 
  ## generator level particle, no smearing and no efficiecny
  set InputArray Delphes/stableParticles  
  set OutputArray jets
  set AreaAlgorithm 1
  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR   0.4
  set JetPTMin     20.0
  #set JetPTMin     15.0
}


############################################
## Jet Pile-Up Subtraction: L1 correction ##
############################################

module JetPileUpSubtractor JetPileUpSubtractor { ## make the rho correction 
  ## input jets
  set JetInputArray FastJetFinder/jets
  ## take the Rho from cluster the event with kt jets (decide to use this or the median grid or the safeAreaSbtraction)
  set RhoInputArray RhoKt4/rho
  ## output jets
  set OutputArray jets
  set doSafe4VAreaSubtraction false
  set JetPTMin 20.0
}

module JetPileUpSubtractor PvJetPileUpSubtractor { ## make the rho correction 
  ## input jets
  set JetInputArray PvFastJetFinder/jets
  #set InputArray PvEFlowMerger/Pveflow
  #set JetInputArray JetMerger/jets
  ## take the Rho from cluster the event with kt jets (decide to use this or the median grid or the safeAreaSbtraction)
  set RhoInputArray PvRhoKt4/rho
  #set RhoInputArray RhoMerger/rho
  #set RhoInputArray RhoKt4/rho
  #set RhoInputArray noTimeTowersRhoKt4/rho
  ## output jets
  set OutputArray jets
  set doSafe4VAreaSubtraction false
  set JetPTMin 20.0
}

## rho correction for no time towers

module JetPileUpSubtractor noTimeTowersJetPileUpSubtractor { ## make the rho correction 
  ## input jets
  #set JetInputArray noTimeTowersFastJetFinder/jets
  set JetInputArray PvJetPileUpSubtractor/jets
  ## take the Rho from cluster the event with kt jets (decide to use this or the median grid or the safeAreaSbtraction)
  #set RhoInputArray noTimeTowersRhoKt4/rho
  set RhoInputArray 2PvRhoKt4/2rho
  ## output jets
  set OutputArray jets
  set doSafe4VAreaSubtraction false
  set JetPTMin 20.0
}





module JetPileUpSubtractor JetPileUpSubtractorGrid { ## make the rho correction 
  ## input jets
  set JetInputArray FastJetFinder/jets
  ## take the Rho from cluster the event with kt jets (decide to use this or the median grid or the safeareasbtraction)
  set RhoInputArray RhoGridFastJet/rho
  ## output jets
  set OutputArray jets
  set doSafe4VAreaSubtraction false
  set JetPTMin 20.0
}

module JetPileUpSubtractor JetPileUpSubtractor4VArea { ## make the rho correction using safe 4V subtraction
  ## input jets
  set JetInputArray FastJetFinder/jets
  ## not used when doSafe4VAreaSubtraction is true
  set RhoInputArray RhoGridFastJet/rho
  ## output jets
  set OutputArray jets
  set JetPTMin 20.0

  ## options for 4V safe subtracion
  set doSafe4VAreaSubtraction true
  ## use this info only if doSafe4VAreaSubtraction is set to true
  set InputArray EFlowMerger/eflow
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithmRho 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithmRho 4
  set ParameterRRho   0.4
  set GhostEtaMaxRho  5.0
  set RhoEtaMaxRho    5.0

  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                       
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                
  set JetAlgorithm 6
  set ParameterR   0.4


  ## eta bins for rho evaluation
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0

}

module JetFlavourAssociation  JetFlavourAssociation {

  set PartonInputArray    Delphes/partons
  set ParticleInputArray  Delphes/allParticles
  set LHEPartonInputArray Delphes/LHEParticles
  set JetInputArray       JetPileUpSubtractor/jets
  set DeltaR        0.4
  set PartonPTMin   0.5
  set PartonEtaMax  2.5
 
}


module Merger JetMerger {
  ## take smeared charged hadron and electrons as only tracks 
  ## to take into account in the calorimeter simulation
  add InputArray PvJetPileUpSubtractor/jets
  add InputArray noTimeTowersJetPileUpSubtractor/jets
  #add InputArray noTimeTowersFastJetFinder/jets
  #add InputArray PvFastJetFinder/jets
  set OutputArray jets
}

module RhoMerger RhoMerger {
  ## take smeared charged hadron and electrons as only tracks 
  ## to take into account in the calorimeter simulation
  add InputArray PvRhoKt4/rho
  add InputArray 2PvRhoKt4/2rho
  #add InputArray noTimeTowersFastJetFinder/jets
  #add InputArray PvFastJetFinder/jets
  set OutputArray rho
}

################################################################################
### Neutrino Filter on generated particles  of status 1 without any smearing ###
################################################################################
module NeutrinoFilter NeutrinoFilter {
  set InputArray  Delphes/stableParticles  
  set OutputArray stableParticles
}

### make GenJets without neutrino
module FastJetFinder GenJetFinderNoNu {
  ## input particles
  set InputArray NeutrinoFilter/stableParticles
  ## output name
  set OutputArray jets
  set AreaAlgorithm 1  
  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4  
  set JetPTMin   20.0
}

### -sum of all particles after filtering neutrinos
module Merger GenMissingET {
  add InputArray NeutrinoFilter/stableParticles
  set MomentumOutputArray momentum
}

###########################
### Run the puppi code  ###
###########################
module Merger TrackMergerWithMuon {
  ## take smeared charged hadron and electrons as only tracks 
  ## to take into account in the calorimeter simulation
  add InputArray ECalorimeter/ecalEflowTracks    
  add InputArray MuonMomentumSmearing/muons 
  set OutputArray tracks
}


module RunPUPPI RunPUPPI {
  ## input information
  set TrackInputArray   TrackMergerWithMuon/tracks
  set NeutralInputArray Calorimeter/eflowTowers
  set PVInputArray      ModifyBeamSpot/PV
  ## min puppi weight and use dZ vertex option
  set MinPuppiWeight    0.05
  set UseExp            false
  ## define puppi algorithm parameters (more than one for the same eta region is possible) 
  add EtaMinBin           0.    2.5    2.5    3.0   3.0      
  add EtaMaxBin           2.5   3.0    3.0    10.0  10.0
  add PtMinBin            0.    0.5    0.5    0.5   0.5    
  add ConeSizeBin         0.25  0.25   0.25   0.25  0.25
  add RMSPtMinBin         0.1   0.5    0.5    0.5   0.5
  add RMSScaleFactorBin   1.0   1.0    1.0    1.0   1.0
  add NeutralMinEBin      0.2   1.0    1.0    1.5   1.5
  add NeutralPtSlope      0.02  0.02   0.02   0.02  0.02
  add ApplyCHS            true  true   true   true  true
  add UseCharged          true  false  false  false false
  add ApplyLowPUCorr      true  true   true   true  true
  add MetricId            5     5      1      5     1
   ## output name
  set OutputArray PuppiParticles
  set OutputArrayTracks   puppiTracks
  set OutputArrayNeutrals puppiNeutrals
} 

#####################
## Make puppi jets ##
#####################

module FastJetFinder PuppiJetFinder {
  set InputArray RunPUPPI/PuppiParticles
  set OutputArray jets
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR   0.4
  set JetPTMin     10.
}

#####################################
## Maker rho corrections for Puppi ##
#####################################

module FastJetFinder PuppiRhoKt4 {
  set InputArray RunPUPPI/PuppiParticles
  ## compute rho clustering the event  
  set ComputeRho     true
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 4
  set ParameterR   0.4
  set GhostEtaMax  5.0
  set RhoEtaMax    5.0  
  add RhoEtaRange  0.0 2.5
  add RhoEtaRange  2.5 4.0
  add RhoEtaRange  4.0 5.0
  set JetPTMin     0.0
}


module FastJetFinder PuppiRhoGridFastJet {
  ## take as input the particle flow particles
  set InputArray RunPUPPI/PuppiParticles
  ## compute rho clustering the event  
  set ComputeRho     false
  ## not compute rho using the grid
  set ComputeRhoGrid true
  set RhoOutputArray rho
  ## area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 1 
  ## Clustering and Ghost parameter
  set GhostEtaMax 5.0
  set RhoEtaMax   5.0
  ## eta bins for rho evaluation
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 4.0
  add RhoEtaRange 4.0 5.0
  set JetPTMin 0.0
}

module FastJetFinder PuppiRhoGrid {
  ## take as input the particle flow particles
  set InputArray RunPUPPI/PuppiParticles
  ## compute rho clustering the event  
  set ComputeRho     false
  ## not compute rho using the grid
  set ComputeRhoGridParticles true
  set RhoOutputArray rho
  ## eta bins for rho evaluation
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 5.0
  set JetPTMin 0.0
}

########################
## Correct puppi jets ##
########################

module JetPileUpSubtractor PuppiJetPileUpSubtractor { ## make the rho correction 
  set JetInputArray PuppiJetFinder/jets
  ## take the Rho from cluster the event with kt jets (decide to use this or the median grid or the safeAreaSbtraction)
  set RhoInputArray PuppiRhoGrid/rho
  set OutputArray jets
  set doSafe4VAreaSubtraction false
  set JetPTMin 20.0
}



module JetPileUpSubtractor PuppiJetPileUpSubtractorGrid { ## make the rho correction 
  set JetInputArray PuppiJetFinder/jets
  ## take the Rho from cluster the event with kt jets (decide to use this or the median grid or the safeareasbtraction)
  set RhoInputArray PuppiRhoGridFastJet/rho
  set OutputArray jets
  set doSafe4VAreaSubtraction false
  set JetPTMin 20.0
}

module JetPileUpSubtractor PuppiJetPileUpSubtractor4VArea { ## make the rho correction 
  set JetInputArray PuppiJetFinder/jets
  ## not used when doSafe4VAreaSubtraction is true
  set RhoInputArray PuppiRhoGridFastJet/rho
  set OutputArray   jets
  set doSafe4VAreaSubtraction true
  set JetPTMin      20.0
  ## use this info only if doSafe4VAreaSubtraction is set to true
  set InputArray    RunPUPPI/PuppiParticles
  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithmRho 1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithmRho 4
  set ParameterRRho   0.4
  set GhostEtaMaxRho  5.0
  set RhoEtaMaxRho    5.0

  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area                                       
  set AreaAlgorithm   1
  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt                                                                                                
  set JetAlgorithm    6
  set ParameterR      0.4

  ## eta bins for rho evaluation
  add RhoEtaRange 0.0 2.5
  add RhoEtaRange 2.5 3.0
  add RhoEtaRange 3.0 10.0
}

#############################
## Jet Flavour Association ##
#############################

module JetFlavourAssociation  PuppiJetFlavourAssociation {

  set PartonInputArray    Delphes/partons
  set ParticleInputArray  Delphes/allParticles
  set LHEPartonInputArray Delphes/LHEParticles
  set JetInputArray       PuppiJetPileUpSubtractor/jets

  set DeltaR       0.4
  set PartonPTMin  0.5
  set PartonEtaMax 2.5
 
}

#####################
# Missing ET merger #
#####################

module Merger MissingET {
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}


###############
## PUPPI MET ##
###############

module Merger PuppiMissingET {
  add InputArray RunPUPPI/PuppiParticles
  set MomentumOutputArray momentum
}

#####################
# Photon efficiency #
#####################

module Efficiency PhotonEfficiency {
  ## input particles
  set InputArray ECalorimeter/photons
  ## output particles
  set OutputArray photons
  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for photons
  set EfficiencyFormula {                                      (pt <= 10.0) * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 10.0)  * (0.9635) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 10.0)  * (0.9624) + \
                         (abs(eta) > 2.5)                                   * (0.00)
  }
}

####################
# Photon isolation #
####################

module Isolation PhotonIsolation {
  # particle for which calculate the isolation
  set CandidateInputArray        PhotonEfficiency/photons 
  # neutral and charged particles for the whole event (no CHS applied)
  set NeutralIsolationInputArray ECalorimeter/ecalEflowTowers
  set ChargedIsolationInputArray TrackMergerWithMuon/tracks
  # select a rho for the isolation
  set RhoInputArray RhoKt4/rho
  # output array
  set OutputArray photons
  # isolation cone
  set DeltaRMax 0.3
  # minimum pT  
  set PTMin     0.5
  # iso ratio to cut
  set PTRatioMax 9999.
}


#######################
# Electron efficiency #
#######################

module Efficiency ElectronEfficiency {
  set InputArray TrackPileUpSubtractor/electrons
  set OutputArray electrons
  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for electrons
  set EfficiencyFormula {    (pt <= 4.0)  * (0.00) + \
                             (abs(eta) <= 1.45 ) * (pt >  4.0 && pt <= 6.0)   * (0.50) + \
                             (abs(eta) <= 1.45 ) * (pt >  6.0 && pt <= 8.0)   * (0.70) + \
                             (abs(eta) <= 1.45 ) * (pt >  8.0 && pt <= 10.0)  * (0.85) + \
                             (abs(eta) <= 1.45 ) * (pt > 10.0 && pt <= 30.0)  * (0.94) + \                                                      
                             (abs(eta) <= 1.45 ) * (pt > 30.0 && pt <= 50.0)  * (0.97) + \                          
                             (abs(eta) <= 1.45 ) * (pt > 50.0 && pt <= 70.0)  * (0.98) + \          
                             (abs(eta) <= 1.45 ) * (pt > 70.0 )  * (1.0) + \   
                             (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt >  4.0 && pt <= 10.0) * (0.35) + \
                             (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 10.0 && pt <= 30.0) * (0.40) + \   
                             (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 30.0 && pt <= 70.0) * (0.45) + \                                 
                             (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 70.0 )  * (0.45) + \    
                             (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt >  4.0 && pt <= 10.0)  * (0.75) + \
                             (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 10.0 && pt <= 30.0)  * (0.85) + \                                                      
                             (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 30.0 && pt <= 50.0)  * (0.95) + \                          
                             (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 50.0 && pt <= 70.0)  * (0.95) + \          
                             (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 70.0 )  * (1.0) + \   
                             (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt >  4.0 && pt <= 10.0)  * (0.65) + \
                             (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt > 10.0 && pt <= 30.0)  * (0.75) + \                                                      
                             (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt > 30.0 && pt <= 50.0)  * (0.85) + \                          
                             (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt > 50.0 && pt <= 70.0)  * (0.85) + \          
                             (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt > 70.0 )  * (0.85) + \                                                                       
     	                     (abs(eta) > 2.5)                              * (0.00)
    }
}

######################
# Electron isolation #
######################

module Isolation ElectronIsolation {
  set CandidateInputArray        ElectronEfficiency/electrons
  set NeutralIsolationInputArray ECalorimeter/ecalEflowTowers
  set ChargedIsolationInputArray TrackMergerWithMuon/tracks
  set RhoInputArray RhoKt4/rho
  set OutputArray electrons
  set DeltaRMax 0.3
  set PTMin 0.5
  set PTRatioMax 9999.
}

###################
# Muon efficiency #
###################

module Efficiency MuonEfficiency {
  set InputArray TrackPileUpSubtractor/muons
  set OutputArray muons
  # set EfficiencyFormula {efficiency as a function of eta and pt}
  # efficiency formula for muons
  set EfficiencyFormula {                                  (pt <= 2.0)  * (0.00) + \  
                         (abs(eta) <= 2.40) * (pt >  2.0 && pt <= 3.0)  * (0.51) + \
                         (abs(eta) <= 2.40) * (pt >  3.0 && pt <= 4.0)  * (0.85) + \ 
                         (abs(eta) <= 2.40) * (pt >  4.0 && pt <= 11.0) * (0.93) + \               
                         (abs(eta) <= 2.40) * (pt >  11. && pt <= 50.)  * (0.96) + \   
                         (abs(eta) <= 2.40) * (pt >  50. && pt <= 70.)  * (0.98) + \                      
                         (abs(eta) <= 2.40) * (pt > 70.0 )  * (1.00) + \   
 	                 (abs(eta) > 2.40)  * (0.00)
  }
}

##################
# Muon isolation #
##################

module Isolation MuonIsolation {
  set CandidateInputArray MuonEfficiency/muons
  set NeutralIsolationInputArray ECalorimeter/ecalEflowTowers
  set ChargedIsolationInputArray TrackMergerWithMuon/tracks 
  set RhoInputArray RhoKt4/rho
  set OutputArray muons
  set DeltaRMax 0.3
  set PTMin 0.5
  set PTRatioMax 9999.
}


#############
# b-tagging #
#############
module BTagging BTagging {

  set JetInputArray JetPileUpSubtractor/jets

  add EfficiencyFormulaLoose {0} {(pt <= 20.0) * (0.000) + \
			          (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.096) + \
 			          (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.112) + \
				  (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.0812) + \
				  (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.0908) + \
				  (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.0792) + \
				  (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.0858) + \
				  (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.0917) + \
				  (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.084) + \
				  (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.0906) + \
				  (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.0989) + \
				  (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.1022) + \
				  (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.1035) + \
				  (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.1095) + \
				  (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.1201) + \
				  (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.1348) + \
				  (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.1482) + \
				  (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.1629) + \
				  (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.1775) + \
				  (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.2002) + \
				  (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.1977) + \
				  (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.2084) + \
				  (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.2195) + \
				  (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.2424) + \
				  (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.2909) + \
				  (abs(eta) <= 1.8) * (pt > 2000.0) * (0.3457) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.074765) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.100053) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.071492) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.084796) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.076927) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.08424) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.093118) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.084629) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.092977) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.10206) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.102344) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.098435) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.105507) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.112841) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.126329) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.140759) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.153193) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.107869) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.119527) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.08688) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.089324) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.097172) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.000) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.000) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.000) + \
				  (abs(eta) > 2.4) * (0.000)
  }

  add EfficiencyFormulaLoose {4} { (pt <= 20.0) * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.387) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.448) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.408) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.427) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.408) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.425) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.426) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.4) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.415) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.416) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.405) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.387) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.39) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.389) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.389) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.381) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.381) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.367) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.369) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.326) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.335) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.326) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.341) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.403) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0) * (0.47) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.2497) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.31891) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.273) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.28445) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.28036) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.28453) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.29495) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.27024) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.28912) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.29048) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.27507) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.23327) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.2493) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.2416) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.26652) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.24852) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.26927) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.19302) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.19433) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.17523) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.14981) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.16666) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.) + \
				       (abs(eta) > 2.4) * (0.000)
  }

   # efficiency formula for b-jets
  add EfficiencyFormulaLoose {5} { (pt <= 20.0) * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.75) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.827) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.837) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.853) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.855) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.862) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.868) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.865) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.863) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.857) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.851) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.838) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.831) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.817) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.796) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.772) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.76) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.743) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.703) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.638) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.605) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.572) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.541) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.567) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0) * (0.603) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.6063) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.7188) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.72) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.7365) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.7462) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.7454) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.7415) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.727) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.7112) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.7112) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.6754) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.6359) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.6348) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.6115) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.5585) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.5608) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.5208) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.456) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.4524) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.388) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.3928) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.3823) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
				       (abs(eta) > 2.4) * (0.000)
  }
 
  add EfficiencyFormulaMedium {0} { (pt <= 20.0) * (0.000) + \
				    (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.00654) + \
				    (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.00921) + \
					(abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.00573) + \
					(abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.00694) + \
					(abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.0062) + \
					(abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.00708) + \
					(abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.00779) + \
					(abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.00693) + \
					(abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.00777) + \
					(abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.00862) + \
					(abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.01038) + \
					(abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.01189) + \
					(abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.01279) + \
					(abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.01452) + \
					(abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.01696) + \
					(abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.01958) + \
					(abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.02253) + \
					(abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.01787) + \
					(abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.02154) + \
					(abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.01839) + \
					(abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.01987) + \
					(abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.02351) + \
					(abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.02937) + \
					(abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.04001) + \
					(abs(eta) <= 1.8) * (pt > 2000.0) * (0.0542) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.004236) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.006653) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.005512) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.00754) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.005813) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.006439) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.008063) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.00647) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.007583) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.008543) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.01034) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.011253) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.012945) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.014474) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.017361) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.020912) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.023139) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.010756) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.012569) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.006046) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.006428) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.00887) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.000) + \
					(abs(eta) > 2.4) * (0.000)
  }

 # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormulaMedium {4} { (pt <= 20.0) * (0.000) + \
					(abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.1102) + \
					(abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.1344) + \
					(abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.1025) + \
					(abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.1025) + \
					(abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.1063) + \
					(abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.1087) + \
					(abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.1124) + \
					(abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.109) + \
					(abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.111) + \
					(abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.1091) + \
					(abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.1087) + \
					(abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.1091) + \
					(abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.1107) + \
					(abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.1061) + \
					(abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.1017) + \
					(abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.1) + \
					(abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.0966) + \
					(abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.0697) + \
					(abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.0679) + \
					(abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.0503) + \
					(abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.0514) + \
					(abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.0481) + \
					(abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.0667) + \
					(abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.0861) + \
					(abs(eta) <= 1.8) * (pt > 2000.0) * (0.092) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.03331) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.04361) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.03863) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.04287) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.04431) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.04452) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.04339) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.0436) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.0456) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.05138) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.04794) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.04004) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.04713) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.04515) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.05314) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.05143) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.05936) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.02357) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.03222) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.01523) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.02621) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.01709) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.) + \
					(abs(eta) > 2.4) * (0.000)
  }

 # efficiency formula for b-jets
  add EfficiencyFormulaMedium {5} { (pt <= 20.0) * (0.000) + \
					(abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.536) + \
					(abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.6439) + \
					(abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.6504) + \
					(abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.6716) + \
					(abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.6841) + \
					(abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.6896) + \
					(abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.6916) + \
					(abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.6882) + \
					(abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.6838) + \
					(abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.6715) + \
					(abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.6554) + \
					(abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.6366) + \
					(abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.6192) + \
					(abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.595) + \
					(abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.5551) + \
					(abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.5138) + \
					(abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.4884) + \
					(abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.4009) + \
					(abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.3459) + \
					(abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.2523) + \
					(abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.2404) + \
					(abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.2198) + \
					(abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.2263) + \
					(abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.2614) + \
					(abs(eta) <= 1.8) * (pt > 2000.0) * (0.3194) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.3254) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.4339) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.4499) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.4716) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.4766) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.4788) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.4863) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.4891) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.462) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.4583) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.4247) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.3775) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.3734) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.3348) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.2939) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.285) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.2421) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.1565) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.1522) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.1231) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.1607) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.1323) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
					(abs(eta) > 2.4) * (0.000)
  }

  add EfficiencyFormulaTight {0} { (pt <= 20.0) * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.000164) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.000235) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.000266) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.000329) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.000309) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.000309) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.000546) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.000499) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.000642) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.000742) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.000928) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.001323) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.001392) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.00154) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.002094) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.002427) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.002927) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.001854) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.002355) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.002297) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.002433) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.002706) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.003602) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.004987) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0) * (0.007414) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.000573) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.000574) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.000521) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.000786) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.000539) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.000673) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.000934) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.000781) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.000949) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.000977) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.001168) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.000879) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.000812) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.001215) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.001679) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.001893) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.002723) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.003555) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.003881) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.006046) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.005563) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.007611) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.000) + \
				       (abs(eta) > 2.4) * (0.000)
  }

 # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormulaTight {4} { (pt <= 20.0) * (0.000) + \
				  (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.00531) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.00567) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.0064) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.00673) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.00766) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.00729) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.00674) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.00824) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.00888) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.00919) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.01021) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.01041) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.01027) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.00999) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.01047) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.01014) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.01021) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.00601) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.0054) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.00487) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.00519) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.00469) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.00651) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.01299) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0) * (0.00897) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.01014) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.01288) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.01392) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.01533) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.01508) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.01579) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.01106) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.01346) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.01315) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.01156) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.0082) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.00439) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.00744) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.00685) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.00755) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.00706) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.00428) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.01031) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.00976) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.01523) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.00374) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.00854) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
				       (abs(eta) > 2.4) * (0.000)
  }

  # efficiency formula for b-jets
  add EfficiencyFormulaTight {5} { (pt <= 20.0) * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.2426) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.327) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.3559) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.3704) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.3824) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.3844) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.3848) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.3862) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.3778) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.3622) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.3299) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.2889) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.2815) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.253) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.221) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.1963) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.1739) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.0992) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.0788) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.0581) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.0534) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.0521) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.0626) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.0826) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0) * (0.1022) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.1562) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.2499) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.2956) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.3128) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.3147) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.3222) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.304) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.3051) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.2657) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.2578) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.2087) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.1634) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.1651) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.1353) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.109) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.0799) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.0699) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.054) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.0718) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.0746) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.0803) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.0882) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
				       (abs(eta) > 2.4) * (0.000)
  }

}


module BTagging PuppiBTagging {

  set JetInputArray PuppiJetPileUpSubtractor/jets

  add EfficiencyFormulaLoose {0} {(pt <= 20.0) * (0.000) + \
			          (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.096) + \
 			          (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.112) + \
				  (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.0812) + \
				  (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.0908) + \
				  (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.0792) + \
				  (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.0858) + \
				  (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.0917) + \
				  (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.084) + \
				  (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.0906) + \
				  (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.0989) + \
				  (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.1022) + \
				  (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.1035) + \
				  (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.1095) + \
				  (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.1201) + \
				  (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.1348) + \
				  (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.1482) + \
				  (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.1629) + \
				  (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.1775) + \
				  (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.2002) + \
				  (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.1977) + \
				  (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.2084) + \
				  (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.2195) + \
				  (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.2424) + \
				  (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.2909) + \
				  (abs(eta) <= 1.8) * (pt > 2000.0) * (0.3457) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.074765) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.100053) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.071492) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.084796) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.076927) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.08424) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.093118) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.084629) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.092977) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.10206) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.102344) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.098435) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.105507) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.112841) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.126329) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.140759) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.153193) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.107869) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.119527) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.08688) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.089324) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.097172) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.000) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.000) + \
				  (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.000) + \
				  (abs(eta) > 2.4) * (0.000)
  }

  add EfficiencyFormulaLoose {4} { (pt <= 20.0) * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.387) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.448) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.408) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.427) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.408) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.425) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.426) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.4) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.415) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.416) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.405) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.387) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.39) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.389) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.389) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.381) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.381) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.367) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.369) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.326) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.335) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.326) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.341) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.403) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0) * (0.47) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.2497) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.31891) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.273) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.28445) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.28036) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.28453) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.29495) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.27024) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.28912) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.29048) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.27507) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.23327) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.2493) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.2416) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.26652) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.24852) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.26927) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.19302) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.19433) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.17523) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.14981) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.16666) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.) + \
				       (abs(eta) > 2.4) * (0.000)
  }

   # efficiency formula for b-jets
  add EfficiencyFormulaLoose {5} { (pt <= 20.0) * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.75) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.827) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.837) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.853) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.855) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.862) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.868) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.865) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.863) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.857) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.851) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.838) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.831) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.817) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.796) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.772) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.76) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.743) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.703) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.638) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.605) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.572) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.541) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.567) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0) * (0.603) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.6063) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.7188) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.72) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.7365) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.7462) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.7454) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.7415) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.727) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.7112) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.7112) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.6754) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.6359) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.6348) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.6115) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.5585) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.5608) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.5208) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.456) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.4524) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.388) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.3928) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.3823) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
				       (abs(eta) > 2.4) * (0.000)
  }
 
  add EfficiencyFormulaMedium {0} { (pt <= 20.0) * (0.000) + \
				    (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.00654) + \
				    (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.00921) + \
					(abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.00573) + \
					(abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.00694) + \
					(abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.0062) + \
					(abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.00708) + \
					(abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.00779) + \
					(abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.00693) + \
					(abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.00777) + \
					(abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.00862) + \
					(abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.01038) + \
					(abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.01189) + \
					(abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.01279) + \
					(abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.01452) + \
					(abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.01696) + \
					(abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.01958) + \
					(abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.02253) + \
					(abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.01787) + \
					(abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.02154) + \
					(abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.01839) + \
					(abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.01987) + \
					(abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.02351) + \
					(abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.02937) + \
					(abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.04001) + \
					(abs(eta) <= 1.8) * (pt > 2000.0) * (0.0542) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.004236) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.006653) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.005512) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.00754) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.005813) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.006439) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.008063) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.00647) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.007583) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.008543) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.01034) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.011253) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.012945) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.014474) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.017361) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.020912) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.023139) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.010756) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.012569) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.006046) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.006428) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.00887) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.000) + \
					(abs(eta) > 2.4) * (0.000)
  }

 # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormulaMedium {4} { (pt <= 20.0) * (0.000) + \
					(abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.1102) + \
					(abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.1344) + \
					(abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.1025) + \
					(abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.1025) + \
					(abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.1063) + \
					(abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.1087) + \
					(abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.1124) + \
					(abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.109) + \
					(abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.111) + \
					(abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.1091) + \
					(abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.1087) + \
					(abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.1091) + \
					(abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.1107) + \
					(abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.1061) + \
					(abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.1017) + \
					(abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.1) + \
					(abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.0966) + \
					(abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.0697) + \
					(abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.0679) + \
					(abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.0503) + \
					(abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.0514) + \
					(abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.0481) + \
					(abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.0667) + \
					(abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.0861) + \
					(abs(eta) <= 1.8) * (pt > 2000.0) * (0.092) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.03331) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.04361) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.03863) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.04287) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.04431) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.04452) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.04339) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.0436) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.0456) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.05138) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.04794) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.04004) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.04713) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.04515) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.05314) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.05143) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.05936) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.02357) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.03222) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.01523) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.02621) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.01709) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.) + \
					(abs(eta) > 2.4) * (0.000)
  }

 # efficiency formula for b-jets
  add EfficiencyFormulaMedium {5} { (pt <= 20.0) * (0.000) + \
					(abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.536) + \
					(abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.6439) + \
					(abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.6504) + \
					(abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.6716) + \
					(abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.6841) + \
					(abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.6896) + \
					(abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.6916) + \
					(abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.6882) + \
					(abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.6838) + \
					(abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.6715) + \
					(abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.6554) + \
					(abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.6366) + \
					(abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.6192) + \
					(abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.595) + \
					(abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.5551) + \
					(abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.5138) + \
					(abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.4884) + \
					(abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.4009) + \
					(abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.3459) + \
					(abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.2523) + \
					(abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.2404) + \
					(abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.2198) + \
					(abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.2263) + \
					(abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.2614) + \
					(abs(eta) <= 1.8) * (pt > 2000.0) * (0.3194) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.3254) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.4339) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.4499) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.4716) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.4766) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.4788) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.4863) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.4891) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.462) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.4583) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.4247) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.3775) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.3734) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.3348) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.2939) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.285) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.2421) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.1565) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.1522) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.1231) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.1607) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.1323) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
					(abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
					(abs(eta) > 2.4) * (0.000)
  }

  add EfficiencyFormulaTight {0} { (pt <= 20.0) * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.000164) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.000235) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.000266) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.000329) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.000309) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.000309) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.000546) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.000499) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.000642) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.000742) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.000928) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.001323) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.001392) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.00154) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.002094) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.002427) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.002927) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.001854) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.002355) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.002297) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.002433) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.002706) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.003602) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.004987) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0) * (0.007414) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.000573) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.000574) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.000521) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.000786) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.000539) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.000673) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.000934) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.000781) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.000949) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.000977) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.001168) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.000879) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.000812) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.001215) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.001679) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.001893) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.002723) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.003555) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.003881) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.006046) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.005563) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.007611) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.000) + \
				       (abs(eta) > 2.4) * (0.000)
  }

 # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormulaTight {4} { (pt <= 20.0) * (0.000) + \
				  (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.00531) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.00567) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.0064) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.00673) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.00766) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.00729) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.00674) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.00824) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.00888) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.00919) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.01021) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.01041) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.01027) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.00999) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.01047) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.01014) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.01021) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.00601) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.0054) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.00487) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.00519) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.00469) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.00651) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.01299) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0) * (0.00897) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.01014) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.01288) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.01392) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.01533) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.01508) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.01579) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.01106) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.01346) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.01315) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.01156) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.0082) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.00439) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.00744) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.00685) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.00755) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.00706) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.00428) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.01031) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.00976) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.01523) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.00374) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.00854) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
				       (abs(eta) > 2.4) * (0.000)
  }

  # efficiency formula for b-jets
  add EfficiencyFormulaTight {5} { (pt <= 20.0) * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30) * (0.2426) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40) * (0.327) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50) * (0.3559) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60) * (0.3704) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70) * (0.3824) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80) * (0.3844) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90) * (0.3848) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100) * (0.3862) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120) * (0.3778) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140) * (0.3622) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160) * (0.3299) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180) * (0.2889) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200) * (0.2815) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250) * (0.253) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300) * (0.221) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350) * (0.1963) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400) * (0.1739) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500) * (0.0992) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600) * (0.0788) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700) * (0.0581) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800) * (0.0534) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000) * (0.0521) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.0626) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.0826) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0) * (0.1022) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0) * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30) * (0.1562) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40) * (0.2499) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50) * (0.2956) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60) * (0.3128) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70) * (0.3147) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80) * (0.3222) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90) * (0.304) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100) * (0.3051) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120) * (0.2657) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140) * (0.2578) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160) * (0.2087) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180) * (0.1634) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200) * (0.1651) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250) * (0.1353) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300) * (0.109) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350) * (0.0799) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400) * (0.0699) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500) * (0.054) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600) * (0.0718) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 700) * (0.0746) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 700.0 && pt <= 800) * (0.0803) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 800.0 && pt <= 1000) * (0.0882) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1000.0 && pt <= 1400) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 1400.0 && pt <= 2000) * (0.0) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 2000.0) * (0.0) + \
				       (abs(eta) > 2.4) * (0.000)
  }

}


###################
## PileUp Jet ID ##
###################

module PileUpJetID PileUpJetID {
  set JetInputArray            JetPileUpSubtractor/jets
  set TrackInputArray          TrackMergerWithMuon/tracks
  set NeutralInputArray        ECalorimeter/ecalEflowTowers
  set NeutralHadronsInputArray HCalorimeter/hcalEflowTowers
  set PVInputArray             ModifyBeamSpot/PV

  set OutputArray   jets

  set ParameterR      0.4  
  set JetPTMin        20.0
  set UseConstituents 0
 
  add Cones  0.1 0.2 0.3 0.4 0.5 0.6 0.7

  add Pt010_Tight_betaStar   0.15 0.15 999. 999.
  add Pt1020_Tight_betaStar  0.15 0.15 999. 999.
  add Pt2030_Tight_betaStar  0.15 0.15 999. 999.
  add Pt3050_Tight_betaStar  0.15 0.15 999. 999.

  add Pt010_Tight_RMS  0.06 0.07 0.04 0.05 
  add Pt1020_Tight_RMS 0.06 0.07 0.04 0.05
  add Pt2030_Tight_RMS 0.05 0.07 0.03 0.045
  add Pt3050_Tight_RMS 0.05 0.06 0.03 0.04

  add Pt010_Medium_betaStar  0.2 0.3 999. 999.
  add Pt1020_Medium_betaStar 0.2 0.3 999. 999.
  add Pt2030_Medium_betaStar 0.2 0.3 999. 999.
  add Pt3050_Medium_betaStar 0.2 0.3 999. 999.

  add Pt010_Medium_RMS    0.06 0.03 0.03 0.04
  add Pt1020_Medium_RMS   0.06 0.03 0.03 0.04
  add Pt2030_Medium_RMS   0.06 0.03 0.03 0.04
  add Pt3050_Medium_RMS   0.06 0.03 0.03 0.04

  add Pt010_Loose_betaStar  0.2 0.3 999. 999
  add Pt1020_Loose_betaStar 0.2 0.3 999. 999
  add Pt2030_Loose_betaStar 0.2 0.3 999. 999
  add Pt3050_Loose_betaStar 0.2 0.3 999. 999

  add Pt010_Loose_RMS   0.06 0.05 0.05 0.07
  add Pt1020_Loose_RMS  0.06 0.05 0.05 0.07
  add Pt2030_Loose_RMS  0.06 0.05 0.05 0.055
  add Pt3050_Loose_RMS  0.06 0.05 0.05 0.055
  
}


module PileUpJetID PvPileUpJetID {
  set JetInputArray            PvJetPileUpSubtractor/jets
  set TrackInputArray          TrackMergerWithMuon/tracks
  set NeutralInputArray        ECalorimeter/ecalPvEflowTowers
  set NeutralHadronsInputArray HCalorimeter/hcalEflowTowers
  set PVInputArray             ModifyBeamSpot/PV

  set OutputArray   jets

  set ParameterR      0.4  
  set JetPTMin        20.0
  set UseConstituents 0
 
  add Cones  0.1 0.2 0.3 0.4 0.5 0.6 0.7

  add Pt010_Tight_betaStar   0.15 0.15 999. 999.
  add Pt1020_Tight_betaStar  0.15 0.15 999. 999.
  add Pt2030_Tight_betaStar  0.15 0.15 999. 999.
  add Pt3050_Tight_betaStar  0.15 0.15 999. 999.

  add Pt010_Tight_RMS  0.06 0.07 0.04 0.05 
  add Pt1020_Tight_RMS 0.06 0.07 0.04 0.05
  add Pt2030_Tight_RMS 0.05 0.07 0.03 0.045
  add Pt3050_Tight_RMS 0.05 0.06 0.03 0.04

  add Pt010_Medium_betaStar  0.2 0.3 999. 999.
  add Pt1020_Medium_betaStar 0.2 0.3 999. 999.
  add Pt2030_Medium_betaStar 0.2 0.3 999. 999.
  add Pt3050_Medium_betaStar 0.2 0.3 999. 999.

  add Pt010_Medium_RMS    0.06 0.03 0.03 0.04
  add Pt1020_Medium_RMS   0.06 0.03 0.03 0.04
  add Pt2030_Medium_RMS   0.06 0.03 0.03 0.04
  add Pt3050_Medium_RMS   0.06 0.03 0.03 0.04

  add Pt010_Loose_betaStar  0.2 0.3 999. 999
  add Pt1020_Loose_betaStar 0.2 0.3 999. 999
  add Pt2030_Loose_betaStar 0.2 0.3 999. 999
  add Pt3050_Loose_betaStar 0.2 0.3 999. 999

  add Pt010_Loose_RMS   0.06 0.05 0.05 0.07
  add Pt1020_Loose_RMS  0.06 0.05 0.05 0.07
  add Pt2030_Loose_RMS  0.06 0.05 0.05 0.055
  add Pt3050_Loose_RMS  0.06 0.05 0.05 0.055
  
}

#################################################
## pile up jet id for different rho correction ##
#################################################


module PileUpJetID noTimeTowersPileUpJetID {
  #set JetInputArray            PvJetPileUpSubtractor/jets
  #set JetInputArray            noTimeTowersJetPileUpSubtractor/jets
  set JetInputArray            JetMerger/jets
  set TrackInputArray          TrackMergerWithMuon/tracks
  set NeutralInputArray        ECalorimeter/ecalPvEflowTowers
  set NeutralHadronsInputArray HCalorimeter/hcalEflowTowers
  set PVInputArray             ModifyBeamSpot/PV

  set OutputArray   jets

  set ParameterR      0.4
  set JetPTMin        20.0
  set UseConstituents 0
 
  add Cones  0.1 0.2 0.3 0.4 0.5 0.6 0.7

  add Pt010_Tight_betaStar   0.15 0.15 999. 999.
  add Pt1020_Tight_betaStar  0.15 0.15 999. 999.
  add Pt2030_Tight_betaStar  0.15 0.15 999. 999.
  add Pt3050_Tight_betaStar  0.15 0.15 999. 999.

  add Pt010_Tight_RMS  0.06 0.07 0.04 0.05 
  add Pt1020_Tight_RMS 0.06 0.07 0.04 0.05
  add Pt2030_Tight_RMS 0.05 0.07 0.03 0.045
  add Pt3050_Tight_RMS 0.05 0.06 0.03 0.04

  add Pt010_Medium_betaStar  0.2 0.3 999. 999.
  add Pt1020_Medium_betaStar 0.2 0.3 999. 999.
  add Pt2030_Medium_betaStar 0.2 0.3 999. 999.
  add Pt3050_Medium_betaStar 0.2 0.3 999. 999.

  add Pt010_Medium_RMS    0.06 0.03 0.03 0.04
  add Pt1020_Medium_RMS   0.06 0.03 0.03 0.04
  add Pt2030_Medium_RMS   0.06 0.03 0.03 0.04
  add Pt3050_Medium_RMS   0.06 0.03 0.03 0.04

  add Pt010_Loose_betaStar  0.2 0.3 999. 999
  add Pt1020_Loose_betaStar 0.2 0.3 999. 999
  add Pt2030_Loose_betaStar 0.2 0.3 999. 999
  add Pt3050_Loose_betaStar 0.2 0.3 999. 999

  add Pt010_Loose_RMS   0.06 0.05 0.05 0.07
  add Pt1020_Loose_RMS  0.06 0.05 0.05 0.07
  add Pt2030_Loose_RMS  0.06 0.05 0.05 0.055
  add Pt3050_Loose_RMS  0.06 0.05 0.05 0.055
  
}








module PileUpJetID PuppiPileUpJetID {
  set JetInputArray     PuppiJetPileUpSubtractor/jets
  set TrackInputArray   RunPUPPI/puppiTracks
  set NeutralInputArray RunPUPPI/puppiNeutrals
  set PVInputArray      ModifyBeamSpot/PV

  set OutputArray   jets

  set ParameterR      0.4  
  set JetPTMin        20.0
  set UseConstituents 0
 
  add Cones  0.1 0.2 0.3 0.4 0.5 0.6 0.7

  add Pt010_Tight_betaStar   0.15 0.15 999. 999.
  add Pt1020_Tight_betaStar  0.15 0.15 999. 999.
  add Pt2030_Tight_betaStar  0.15 0.15 999. 999.
  add Pt3050_Tight_betaStar  0.15 0.15 999. 999.

  add Pt010_Tight_RMS  0.06 0.07 0.04 0.04 
  add Pt1020_Tight_RMS 0.06 0.07 0.04 0.04
  add Pt2030_Tight_RMS 0.05 0.07 0.03 0.04
  add Pt3050_Tight_RMS 0.05 0.06 0.03 0.04

  add Pt010_Medium_betaStar  0.2 0.3 999. 999.
  add Pt1020_Medium_betaStar 0.2 0.3 999. 999.
  add Pt2030_Medium_betaStar 0.2 0.3 999. 999.
  add Pt3050_Medium_betaStar 0.2 0.3 999. 999.

  add Pt010_Medium_RMS    0.06 0.03 0.03 0.05
  add Pt1020_Medium_RMS   0.06 0.03 0.03 0.05
  add Pt2030_Medium_RMS   0.06 0.03 0.03 0.045
  add Pt3050_Medium_RMS   0.06 0.03 0.03 0.04

  add Pt010_Loose_betaStar  0.2 0.3 999. 999
  add Pt1020_Loose_betaStar 0.2 0.3 999. 999
  add Pt2030_Loose_betaStar 0.2 0.3 999. 999
  add Pt3050_Loose_betaStar 0.2 0.3 999. 999

  add Pt010_Loose_RMS   0.06 0.05 0.05 0.07
  add Pt1020_Loose_RMS  0.06 0.05 0.05 0.07
  add Pt2030_Loose_RMS  0.06 0.05 0.05 0.055
  add Pt3050_Loose_RMS  0.06 0.05 0.05 0.055
  
}


module Merger PileUpJetIDMissingET {
  add InputArray TrackPileUpSubtractor/eflowTracks
  add InputArray MuonMomentumSmearing/muons
  add InputArray PileUpJetID/eflowTowers
  set MomentumOutputArray momentum
}

########################
## Constituent filter ##
########################

module ConstituentFilter ConstituentFilter {

  set ConEMin 0.

  add JetInputArray GenJetFinderNoNu/jets
  add JetInputArray PileUpJetID/jets

  add ConstituentInputArray Delphes/stableParticles stableParticles
  add ConstituentInputArray TrackPileUpSubtractor/eflowTracks eflowTracks
  add ConstituentInputArray Calorimeter/eflowTowers eflowTowers
  add ConstituentInputArray MuonMomentumSmearing/muons muons

} 


module ConstituentFilter ConstituentFilterPUPPI {

  set ConEMin 0.

  add JetInputArray GenJetFinderNoNu/jets
  add JetInputArray PileUpJetIDPUPPI/jets

  add ConstituentInputArray Delphes/stableParticles stableParticles
  add ConstituentInputArray TrackPileUpSubtractor/eflowTracks eflowTracks
  add ConstituentInputArray Calorimeter/eflowTowers eflowTowers
  add ConstituentInputArray MuonMomentumSmearing/muons muons

} 

##############
## Scalr HT ##
##############

module Merger ScalarHT {
 # add InputArray InputArray
 add InputArray  EFlowMerger/eflow
 set EnergyOutputArray energy
}

module Merger GenScalarHT {
 # add InputArray InputArray
 add InputArray NeutrinoFilter/stableParticles  
 set EnergyOutputArray energy
}

module Merger PuppiScalarHT {
 # add InputArray InputArray
 add InputArray RunPUPPI/PuppiParticles
 set EnergyOutputArray energy
}

 
##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
  ## branch notation : <particle collection> <branch name> <type of object in classes/DelphesClass.h<

  ## input status 1 particle from Pythia8
  #add Branch Delphes/stableParticles GenParticles GenParticle
  ## input partons
  #add Branch Delphes/partons GenParton GenParticle 
  ## LHE particles
  add Branch Delphes/LHEParticles LHEParticles LHEParticle

  ## NPU after Pileup Merging
  add Branch PileUpMerger/NPU NPU ScalarHT

  ## gen particles after vertex smearing 
  #add Branch GenBeamSpotFilter/beamSpotParticles GenBeamSpotParticles GenParticle

  ## particle after B field propagation   
  #add Branch ParticlePropagator/stableParticles particlePropagator GenParticle
  #add Branch ParticlePropagator/electrons       electronPropagator GenParticle
  #add Branch ParticlePropagator/muons           muonPropagator GenParticle 

  ## after Pt filter: all delphes particles, not only status 1 
  #add Branch StatusPid/filteredParticles GenParticles GenParticle
 
  ## track collection after: charged hadrons smearing and track eff, electron smearing and track eff
  #add Branch TrackMerger/tracks trackCollectionNoMU Track

  ## output of the calorimeter simulation
  #add Branch Calorimeter/towers caloTowers Tower
  #add Branch Calorimeter/photons RawPhotons Photon
  #add Branch Calorimeter/eflowTracks eflowTracks Track
  #add Branch Calorimeter/eflowTracks eflowTowers Tower
  
  ## tracks after CHS
  #add Branch TrackPileUpSubtractor/eflowTracks trackCollectionCHS Track

  ## eflow output
  #add Branch EFlowMerger/eflow eflowCandidates Track

  ## Rho values
  ##add Branch GlobalRhoKt4/rho GlobalRhoKt4 Rho
  ##add Branch GlobalRhoGridFastJet/rho GlobalRhoGridFastJet Rho
  add Branch RhoKt4/rho RhoKt4 Rho
  ##add Branch RhoGridFastJet/rho RhoGridFastJet Rho

  ##Rho values for pv
  ##add Branch PvGlobalRhoKt4/rho PvGlobalRhoKt4 Rho
  ##add Branch PvGlobalRhoGridFastJet/rho PvGlobalRhoGridFastJet Rho
  add Branch PvRhoKt4/rho PvRhoKt4 Rho
  add Branch RhoMerger/rho RhoMerger Rho
  ##add Branch PvRhoGridFastJet/rho PvRhoGridFastJet Rho
  
  ## Standard Jets 
  #add Branch FastJetFinder/jets RawJet Jet
  #add Branch noTimeTowersFastJetFinder/jets notimeRawJet Jet
  #add Branch GenJetFinder/jets GenJetWithNu Jet
  add Branch GenJetFinderNoNu/jets GenJet Jet 
  add Branch JetPileUpSubtractor/jets Jet Jet
  add Branch PvJetPileUpSubtractor/jets PvJet Jet
  #add Branch noTimeTowersJetPileUpSubtractor/jets noTimeJet Jet
  #add Branch TrackJetFinder/jets TrackJet Jet si 
  #add Branch PvTrackJetFinder/jets PvTrackJet Jet 
  #add Branch JetPileUpSubtractorGrid/jets Jet Jet
  #add Branch JetPileUpSubtractor4VArea/jets Jet4VArea Jet
  add Branch PileUpJetID/jets JetPUID Jet 
  #add Branch PvPileUpJetID/jets PvJetPUID Jet 
  #add Branch noTimeTowersPileUpJetID/jets PvJetPUID Jet
  add Branch PvFastJetFinder/jets Pvjets Jet 
  #add Branch JetMerger/jets mergedjets Jet
  
  ## PUPPI
  #add Branch RunPUPPI/PuppiParticles puppiParticles GenParticle
  #add Branch PuppiRhoKt4/rho         PuppiRhoKt4 Rho
  #add Branch PuppiRhoGrid/rho PuppiRhoGrid Rho
  #add Branch PuppiJetFinder/jets     RawPuppiJet Jet
  #add Branch PuppiJetPileUpSubtractor/jets PuppiJet Jet
  #add Branch PuppiJetPileUpSubtractorGrid/jets PuppiJetGrid Jet
  #add Branch PuppiJetPileUpSubtractor4VArea/jets PuppiJet4VArea Jet
  #add Branch PuppiPileUpJetID/jets PuppiJetPUID Jet

  ## MET
  ##add Branch GenMissingET/momentum GenMissingET MissingET
  ##add Branch MissingET/momentum MissingET MissingET
  #add Branch PuppiMissingET/momentum PuppiMissingET MissingET

  ## HT
 # add Branch ScalarHT/energy HT ScalarHT si
#  add Branch GenScalarHT/energy GenHT ScalarHT si
  #add Branch PuppiScalarHT/energy PuppiHT ScalarHT

  ## photons and leptons
  #add Branch ElectronIsolation/electrons Electron Electron si
 # add Branch PhotonIsolation/photons Photon Photon si 
#  add Branch MuonIsolation/muons Muon Muon si

  set fOffsetFromModifyBeamSpot 10 
}


#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

#module UniqueObjectFinder UniqueObjectFinderGJ {
#   add InputArray PhotonIsolation/photons photons
#   add InputArray JetPileUpSubtractor/jets jets
#}

#module UniqueObjectFinder UniqueObjectFinderEJ {
#   add InputArray ElectronIsolation/electrons electrons
#   add InputArray UniqueObjectFinderGJ/jets jets
#}

#module UniqueObjectFinder UniqueObjectFinderMJ {
#   add InputArray MuonIsolation/muons muons
#   add InputArray UniqueObjectFinderEJ/jets jets
#}


