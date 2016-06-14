
#
# Makefile for ExRootAnalysis
#
# Author: P. Demin - UCL, Louvain-la-Neuve
#
# multi-platform configuration is taken from ROOT (root/test/Makefile.arch)
#

include doc/Makefile.arch

ifeq ($(ARCH),macosx64)
UNDEFOPT = dynamic_lookup
endif

SrcSuf = cc

LHAPDF=$(LHAPATH)/../../..
HEPMC=$(LHAPATH)/../../../../../hepmc/2.06.07-cms

shell export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(PYTHIA8DATA)/../../../lib
shell export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(LHAPDF)/lib
shell export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(HEPMC)/lib

FASTJET_BASE_TMP =$(shell scram tool info fastjet | grep FASTJET_BASE)
FASTJET_BASE     =$(subst FASTJET_BASE=,,$(FASTJET_BASE_TMP))
FASTJET_LIB_TMP  =$(shell scram tool info fastjet | grep LIB | grep -v LIBDIR)
FASTJET_LIB      =$(subst LIB=,,$(FASTJET_LIB_TMP))
FASTJET_LIBDIR_TMP =$(shell scram tool info fastjet | grep LIBDIR)
FASTJET_LIBDIR   =$(subst LIBDIR=,,$(FASTJET_LIBDIR_TMP))
FASTJET_INCLUDE_TMP  =$(shell scram tool info fastjet | grep INCLUDE)
FASTJET_INCLUDE  =$(subst INCLUDE=,,$(FASTJET_INCLUDE_TMP))

CXXFLAGS += $(ROOTCFLAGS) -Wno-write-strings -D_FILE_OFFSET_BITS=64 -DDROP_CGAL -I. -Iexternal -Iexternal/tcl -I$(FASTJET_INCLUDE)  -Iexternal/LHEActions

DELPHES_LIBS = $(shell $(RC) --libs) -lEG $(SYSLIBS)
DISPLAY_LIBS = $(shell $(RC) --evelibs) $(SYSLIBS)

ifneq ($(CMSSW_FWLITE_INCLUDE_PATH),)
HAS_CMSSW = true
CXXFLAGS += -std=c++11 -I$(subst :, -I,$(CMSSW_FWLITE_INCLUDE_PATH)) 
DELPHES_LIBS += -L$(subst include,lib,$(subst :, -L,$(CMSSW_FWLITE_INCLUDE_PATH))) -lGenVector
ifneq ($(CMSSW_RELEASE_BASE),)
CXXFLAGS += -I$(CMSSW_RELEASE_BASE)/src
endif
ifneq ($(LD_LIBRARY_PATH),)
DELPHES_LIBS += -L$(subst include,lib,$(subst :, -L,$(LD_LIBRARY_PATH))) -lGenVector -lGenVector -lfastjet -lfastjetcontribfragile -lfastjetplugins -lfastjettools -lsiscone -lsiscone_spherical
endif
DELPHES_LIBS += -lFWCoreFWLite -lDataFormatsFWLite -lDataFormatsPatCandidates -lDataFormatsLuminosity -lCommonToolsUtils -lMathCore -lDataFormatsMath  -lGenVector -lGenVector -lfastjet -lfastjetcontribfragile -lfastjetplugins -lfastjettools -lsiscone -lsiscone_spherical
endif

ifneq ($(PROMC),)
HAS_PROMC = true
CXXFLAGS += -I$(PROMC)/include
DELPHES_LIBS += -L$(PROMC)/lib -lprotoc -lprotobuf -lprotobuf-lite -lcbook -lz
endif

ifneq ($(PYTHIA8),)
HAS_PYTHIA8 = true
CXXFLAGS += -I$(PYTHIA8)/include
DELPHES_LIBS += -L$(PYTHIA8)/lib -lpythia8 -lLHAPDF -lgfortran -lz
else
ifneq ($(PYTHIA8DATA),)
HAS_PYTHIA8 = true
CXXFLAGS += -I$(PYTHIA8DATA)/../../../include
DELPHES_LIBS += -L$(PYTHIA8DATA)/../../../lib -lpythia8 -lLHAPDF -lgfortran -lz
endif
endif

###

DELPHES    = libDelphes.$(DllSuf)
DELPHESLIB = libDelphes.lib

DISPLAY    = libDelphesDisplay.$(DllSuf)
DISPLAYLIB = libDelphesDisplay.lib

VERSION = $(shell cat VERSION)
DISTDIR = Delphes-$(VERSION)
DISTTAR = $(DISTDIR).tar.gz

all:


hepmc2pileup$(ExeSuf): \
	tmp/converters/hepmc2pileup.$(ObjSuf)

tmp/converters/hepmc2pileup.$(ObjSuf): \
	converters/hepmc2pileup.cpp \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesHepMCReader.h \
	classes/DelphesPileUpWriter.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
lhco2root$(ExeSuf): \
	tmp/converters/lhco2root.$(ObjSuf)

tmp/converters/lhco2root.$(ObjSuf): \
	converters/lhco2root.cpp \
	modules/Delphes.h \
	classes/DelphesStream.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
pileup2root$(ExeSuf): \
	tmp/converters/pileup2root.$(ObjSuf)

tmp/converters/pileup2root.$(ObjSuf): \
	converters/pileup2root.cpp \
	classes/DelphesStream.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesPileUpReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
root2lhco$(ExeSuf): \
	tmp/converters/root2lhco.$(ObjSuf)

tmp/converters/root2lhco.$(ObjSuf): \
	converters/root2lhco.cpp \
	classes/DelphesClasses.h \
	external/ExRootAnalysis/ExRootTreeReader.h \
	external/ExRootAnalysis/ExRootProgressBar.h
root2pileup$(ExeSuf): \
	tmp/converters/root2pileup.$(ObjSuf)

tmp/converters/root2pileup.$(ObjSuf): \
	converters/root2pileup.cpp \
	classes/DelphesClasses.h \
	classes/DelphesPileUpWriter.h \
	external/ExRootAnalysis/ExRootTreeReader.h \
	external/ExRootAnalysis/ExRootProgressBar.h
stdhep2pileup$(ExeSuf): \
	tmp/converters/stdhep2pileup.$(ObjSuf)

tmp/converters/stdhep2pileup.$(ObjSuf): \
	converters/stdhep2pileup.cpp \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesSTDHEPReader.h \
	classes/DelphesPileUpWriter.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
CountEventData$(ExeSuf): \
	tmp/examples/CountEventData.$(ObjSuf)

tmp/examples/CountEventData.$(ObjSuf): \
	examples/CountEventData.cpp \
	classes/DelphesClasses.h \
	external/ExRootAnalysis/ExRootTreeReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootUtilities.h
Example1$(ExeSuf): \
	tmp/examples/Example1.$(ObjSuf)

tmp/examples/Example1.$(ObjSuf): \
	examples/Example1.cpp \
	classes/DelphesClasses.h \
	external/ExRootAnalysis/ExRootTreeReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootUtilities.h
GeneralExample$(ExeSuf): \
	tmp/examples/GeneralExample.$(ObjSuf)

tmp/examples/GeneralExample.$(ObjSuf): \
	examples/GeneralExample.cpp \
	classes/DelphesClasses.h \
	external/ExRootAnalysis/ExRootTreeReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootUtilities.h
JetExample$(ExeSuf): \
	tmp/examples/JetExample.$(ObjSuf)

tmp/examples/JetExample.$(ObjSuf): \
	examples/JetExample.cpp \
	classes/DelphesClasses.h \
	external/ExRootAnalysis/ExRootTreeReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootUtilities.h
DelphesComparison$(ExeSuf): \
	tmp/external/DelphesAnalysisCodes/DelphesComparison.$(ObjSuf)

tmp/external/DelphesAnalysisCodes/DelphesComparison.$(ObjSuf): \
	external/DelphesAnalysisCodes/DelphesComparison.cpp \
	classes/DelphesClasses.h
DelphesDumper$(ExeSuf): \
	tmp/external/DelphesAnalysisCodes/DelphesDumper.$(ObjSuf)

tmp/external/DelphesAnalysisCodes/DelphesDumper.$(ObjSuf): \
	external/DelphesAnalysisCodes/DelphesDumper.cpp \
	classes/DelphesClasses.h \
	external/ExRootAnalysis/ExRootTreeReader.h
JETResponseStudies$(ExeSuf): \
	tmp/external/DelphesAnalysisCodes/JETResponseStudies.$(ObjSuf)

tmp/external/DelphesAnalysisCodes/JETResponseStudies.$(ObjSuf): \
	external/DelphesAnalysisCodes/JETResponseStudies.cpp \
	classes/DelphesClasses.h
METResponseStudies$(ExeSuf): \
	tmp/external/DelphesAnalysisCodes/METResponseStudies.$(ObjSuf)

tmp/external/DelphesAnalysisCodes/METResponseStudies.$(ObjSuf): \
	external/DelphesAnalysisCodes/METResponseStudies.cpp \
	classes/DelphesClasses.h
EXECUTABLE +=  \
	hepmc2pileup$(ExeSuf) \
	lhco2root$(ExeSuf) \
	pileup2root$(ExeSuf) \
	root2lhco$(ExeSuf) \
	root2pileup$(ExeSuf) \
	stdhep2pileup$(ExeSuf) \
	CountEventData$(ExeSuf) \
	Example1$(ExeSuf) \
	GeneralExample$(ExeSuf) \
	JetExample$(ExeSuf) \
	DelphesComparison$(ExeSuf) \
	DelphesDumper$(ExeSuf) \
	JETResponseStudies$(ExeSuf) \
	METResponseStudies$(ExeSuf)

EXECUTABLE_OBJ +=  \
	tmp/converters/hepmc2pileup.$(ObjSuf) \
	tmp/converters/lhco2root.$(ObjSuf) \
	tmp/converters/pileup2root.$(ObjSuf) \
	tmp/converters/root2lhco.$(ObjSuf) \
	tmp/converters/root2pileup.$(ObjSuf) \
	tmp/converters/stdhep2pileup.$(ObjSuf) \
	tmp/examples/CountEventData.$(ObjSuf) \
	tmp/examples/Example1.$(ObjSuf) \
	tmp/examples/GeneralExample.$(ObjSuf) \
	tmp/examples/JetExample.$(ObjSuf) \
	tmp/external/DelphesAnalysisCodes/DelphesComparison.$(ObjSuf) \
	tmp/external/DelphesAnalysisCodes/DelphesDumper.$(ObjSuf) \
	tmp/external/DelphesAnalysisCodes/JETResponseStudies.$(ObjSuf) \
	tmp/external/DelphesAnalysisCodes/METResponseStudies.$(ObjSuf)

DelphesHepMC$(ExeSuf): \
	tmp/readers/DelphesHepMC.$(ObjSuf)

tmp/readers/DelphesHepMC.$(ObjSuf): \
	readers/DelphesHepMC.cpp \
	modules/Delphes.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesHepMCReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
DelphesLHEF$(ExeSuf): \
	tmp/readers/DelphesLHEF.$(ObjSuf)

tmp/readers/DelphesLHEF.$(ObjSuf): \
	readers/DelphesLHEF.cpp \
	modules/Delphes.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesLHEFReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
DelphesSTDHEP$(ExeSuf): \
	tmp/readers/DelphesSTDHEP.$(ObjSuf)

tmp/readers/DelphesSTDHEP.$(ObjSuf): \
	readers/DelphesSTDHEP.cpp \
	modules/Delphes.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesSTDHEPReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
EXECUTABLE +=  \
	DelphesHepMC$(ExeSuf) \
	DelphesLHEF$(ExeSuf) \
	DelphesSTDHEP$(ExeSuf)

EXECUTABLE_OBJ +=  \
	tmp/readers/DelphesHepMC.$(ObjSuf) \
	tmp/readers/DelphesLHEF.$(ObjSuf) \
	tmp/readers/DelphesSTDHEP.$(ObjSuf)

ifeq ($(HAS_CMSSW),true)
DelphesCMSFWLite$(ExeSuf): \
	tmp/readers/DelphesCMSFWLite.$(ObjSuf)

tmp/readers/DelphesCMSFWLite.$(ObjSuf): \
	readers/DelphesCMSFWLite.cpp \
	modules/Delphes.h \
	classes/DelphesStream.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
EXECUTABLE +=  \
	DelphesCMSFWLite$(ExeSuf)

EXECUTABLE_OBJ +=  \
	tmp/readers/DelphesCMSFWLite.$(ObjSuf)

endif

ifeq ($(HAS_PROMC),true)
DelphesProMC$(ExeSuf): \
	tmp/readers/DelphesProMC.$(ObjSuf)

tmp/readers/DelphesProMC.$(ObjSuf): \
	readers/DelphesProMC.cpp \
	modules/Delphes.h \
	classes/DelphesStream.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h \
	external/ProMC/ProMCBook.h
EXECUTABLE +=  \
	DelphesProMC$(ExeSuf)

EXECUTABLE_OBJ +=  \
	tmp/readers/DelphesProMC.$(ObjSuf)

tmp/external/ProMC/ProMC.pb.$(ObjSuf): \
	external/ProMC/ProMC.pb.$(SrcSuf)
tmp/external/ProMC/ProMCBook.$(ObjSuf): \
	external/ProMC/ProMCBook.$(SrcSuf)
tmp/external/ProMC/ProMCDescription.pb.$(ObjSuf): \
	external/ProMC/ProMCDescription.pb.$(SrcSuf)
tmp/external/ProMC/ProMCHeader.pb.$(ObjSuf): \
	external/ProMC/ProMCHeader.pb.$(SrcSuf)
tmp/external/ProMC/ProMCStat.pb.$(ObjSuf): \
	external/ProMC/ProMCStat.pb.$(SrcSuf)
DELPHES_OBJ +=  \
	tmp/external/ProMC/ProMC.pb.$(ObjSuf) \
	tmp/external/ProMC/ProMCBook.$(ObjSuf) \
	tmp/external/ProMC/ProMCDescription.pb.$(ObjSuf) \
	tmp/external/ProMC/ProMCHeader.pb.$(ObjSuf) \
	tmp/external/ProMC/ProMCStat.pb.$(ObjSuf)

ifeq ($(HAS_PYTHIA8),true)
DELPHES_OBJ +=  \
	
endif

endif

ifeq ($(HAS_PYTHIA8),true)
DelphesPythia8$(ExeSuf): \
	tmp/readers/DelphesPythia8.$(ObjSuf)

tmp/readers/DelphesPythia8.$(ObjSuf): \
	readers/DelphesPythia8.cpp \
	modules/Delphes.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h \
	external/LHEActions/LHEF.h
EXECUTABLE +=  \
	DelphesPythia8$(ExeSuf)

EXECUTABLE_OBJ +=  \
	tmp/readers/DelphesPythia8.$(ObjSuf)

tmp/modules/Pythia8Dict.$(SrcSuf): \
	modules/Pythia8LinkDef.h \
	modules/PileUpMergerPythia8.h
DELPHES_DICT +=  \
	tmp/modules/Pythia8Dict.$(SrcSuf)

DELPHES_DICT_OBJ +=  \
	tmp/modules/Pythia8Dict.$(ObjSuf)

endif

tmp/classes/ClassesDict.$(SrcSuf): \
	classes/ClassesLinkDef.h \
	classes/DelphesModule.h \
	classes/DelphesFactory.h \
	classes/SortableObject.h \
	classes/DelphesClasses.h
tmp/modules/ModulesDict.$(SrcSuf): \
	modules/ModulesLinkDef.h \
	modules/Delphes.h \
	modules/FastJetFinder.h \
	modules/ParticlePropagator.h \
	modules/Efficiency.h \
	modules/EnergySmearing.h \
	modules/MomentumSmearing.h \
	modules/Isolation.h \
	modules/IsoTrackFilter.h \
	modules/EnergyScale.h \
	modules/UniqueObjectFinder.h \
	modules/BTagging.h \
	modules/TauTagging.h \
	modules/TreeWriter.h \
	modules/Merger.h \
	modules/LeptonDressing.h \
	modules/PileUpMerger.h \
	modules/JetPileUpSubtractor.h \
	modules/TrackPileUpSubtractor.h \
	modules/ConstituentFilter.h \
	modules/StatusPidFilter.h \
	modules/Cloner.h \
	modules/Weighter.h \
	modules/ExampleModule.h \
	modules/JetFlavourAssociation.h \
	modules/PileUpJetID.h \
	modules/ModifyBeamSpot.h \
	modules/GenBeamSpotFilter.h \
	modules/RunPUPPI.h \
	modules/NeutrinoFilter.h \
	modules/FakeLepton.h \
	modules/Vertexing.h \
	modules/ECalorimeter.h \
	modules/HCalorimeter.h \
	modules/RhoMerger.h
tmp/external/ExRootAnalysis/ExRootAnalysisDict.$(SrcSuf): \
	external/ExRootAnalysis/ExRootAnalysisLinkDef.h \
	external/ExRootAnalysis/ExRootTreeReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootUtilities.h \
	external/ExRootAnalysis/ExRootClassifier.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootProgressBar.h \
	external/ExRootAnalysis/ExRootConfReader.h \
	external/ExRootAnalysis/ExRootTask.h
DELPHES_DICT +=  \
	tmp/classes/ClassesDict.$(SrcSuf) \
	tmp/modules/ModulesDict.$(SrcSuf) \
	tmp/external/ExRootAnalysis/ExRootAnalysisDict.$(SrcSuf)

DELPHES_DICT_OBJ +=  \
	tmp/classes/ClassesDict.$(ObjSuf) \
	tmp/modules/ModulesDict.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootAnalysisDict.$(ObjSuf)

DELPHES_DICT_v2 +=  \
	tmp/classes/ClassesDict.$(SrcSuf) \
	tmp/modules/ModulesDict.$(SrcSuf) \
	tmp/external/ExRootAnalysis/ExRootAnalysisDict.$(SrcSuf)

DELPHES_DICT_v2_OBJ +=  \
	tmp/classes/ClassesDict.$(ObjSuf) \
	tmp/modules/ModulesDict.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootAnalysisDict.$(ObjSuf)

tmp/display/DisplayDict.$(SrcSuf): \
	display/DisplayLinkDef.h \
	display/DelphesDisplay.h \
	display/DelphesCaloData.h
DISPLAY_DICT +=  \
	tmp/display/DisplayDict.$(SrcSuf)

DISPLAY_DICT_OBJ +=  \
	tmp/display/DisplayDict.$(ObjSuf)

tmp/classes/DelphesClasses.$(ObjSuf): \
	classes/DelphesClasses.$(SrcSuf) \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/SortableObject.h
tmp/classes/DelphesFactory.$(ObjSuf): \
	classes/DelphesFactory.$(SrcSuf) \
	classes/DelphesFactory.h \
	classes/DelphesClasses.h \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/classes/DelphesFormula.$(ObjSuf): \
	classes/DelphesFormula.$(SrcSuf) \
	classes/DelphesFormula.h
tmp/classes/DelphesHepMCReader.$(ObjSuf): \
	classes/DelphesHepMCReader.$(SrcSuf) \
	classes/DelphesHepMCReader.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesStream.h \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/classes/DelphesLHEFReader.$(ObjSuf): \
	classes/DelphesLHEFReader.$(SrcSuf) \
	classes/DelphesLHEFReader.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesStream.h \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/classes/DelphesModule.$(ObjSuf): \
	classes/DelphesModule.$(SrcSuf) \
	classes/DelphesModule.h \
	classes/DelphesFactory.h \
	external/ExRootAnalysis/ExRootTreeReader.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootResult.h
tmp/classes/DelphesPileUpReader.$(ObjSuf): \
	classes/DelphesPileUpReader.$(SrcSuf) \
	classes/DelphesPileUpReader.h
tmp/classes/DelphesPileUpWriter.$(ObjSuf): \
	classes/DelphesPileUpWriter.$(SrcSuf) \
	classes/DelphesPileUpWriter.h
tmp/classes/DelphesSTDHEPReader.$(ObjSuf): \
	classes/DelphesSTDHEPReader.$(SrcSuf) \
	classes/DelphesSTDHEPReader.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/classes/DelphesStream.$(ObjSuf): \
	classes/DelphesStream.$(SrcSuf) \
	classes/DelphesStream.h
tmp/modules/BTagging.$(ObjSuf): \
	modules/BTagging.$(SrcSuf) \
	modules/BTagging.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/Cloner.$(ObjSuf): \
	modules/Cloner.$(SrcSuf) \
	modules/Cloner.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/ConstituentFilter.$(ObjSuf): \
	modules/ConstituentFilter.$(SrcSuf) \
	modules/ConstituentFilter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/Delphes.$(ObjSuf): \
	modules/Delphes.$(SrcSuf) \
	modules/Delphes.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h \
	external/ExRootAnalysis/ExRootConfReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h
tmp/modules/Efficiency.$(ObjSuf): \
	modules/Efficiency.$(SrcSuf) \
	modules/Efficiency.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/EnergyScale.$(ObjSuf): \
	modules/EnergyScale.$(SrcSuf) \
	modules/EnergyScale.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/EnergySmearing.$(ObjSuf): \
	modules/EnergySmearing.$(SrcSuf) \
	modules/EnergySmearing.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/ExampleModule.$(ObjSuf): \
	modules/ExampleModule.$(SrcSuf) \
	modules/ExampleModule.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/FakeLepton.$(ObjSuf): \
	modules/FakeLepton.$(SrcSuf) \
	modules/FakeLepton.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/FastJetFinder.$(ObjSuf): \
	modules/FastJetFinder.$(SrcSuf) \
	modules/FastJetFinder.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/GenBeamSpotFilter.$(ObjSuf): \
	modules/GenBeamSpotFilter.$(SrcSuf) \
	modules/GenBeamSpotFilter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/IsoTrackFilter.$(ObjSuf): \
	modules/IsoTrackFilter.$(SrcSuf) \
	modules/IsoTrackFilter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/Isolation.$(ObjSuf): \
	modules/Isolation.$(SrcSuf) \
	modules/Isolation.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/JetFlavourAssociation.$(ObjSuf): \
	modules/JetFlavourAssociation.$(SrcSuf) \
	modules/JetFlavourAssociation.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/JetPileUpSubtractor.$(ObjSuf): \
	modules/JetPileUpSubtractor.$(SrcSuf) \
	modules/JetPileUpSubtractor.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/LeptonDressing.$(ObjSuf): \
	modules/LeptonDressing.$(SrcSuf) \
	modules/LeptonDressing.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/Merger.$(ObjSuf): \
	modules/Merger.$(SrcSuf) \
	modules/Merger.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/ModifyBeamSpot.$(ObjSuf): \
	modules/ModifyBeamSpot.$(SrcSuf) \
	modules/ModifyBeamSpot.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/MomentumSmearing.$(ObjSuf): \
	modules/MomentumSmearing.$(SrcSuf) \
	modules/MomentumSmearing.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/NeutrinoFilter.$(ObjSuf): \
	modules/NeutrinoFilter.$(SrcSuf) \
	modules/NeutrinoFilter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/ParticlePropagator.$(ObjSuf): \
	modules/ParticlePropagator.$(SrcSuf) \
	modules/ParticlePropagator.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/PileUpJetID.$(ObjSuf): \
	modules/PileUpJetID.$(SrcSuf) \
	modules/PileUpJetID.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/PileUpMerger.$(ObjSuf): \
	modules/PileUpMerger.$(SrcSuf) \
	modules/PileUpMerger.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	classes/DelphesPileUpReader.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/PileUpMergerPythia8.$(ObjSuf): \
	modules/PileUpMergerPythia8.$(SrcSuf) \
	modules/PileUpMergerPythia8.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	classes/DelphesPileUpReader.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/RunPUPPI.$(ObjSuf): \
	modules/RunPUPPI.$(SrcSuf) \
	modules/RunPUPPI.h \
	external/PUPPI/puppiCleanContainer.hh \
	external/PUPPI/RecoObj.hh \
	external/PUPPI/puppiParticle.hh \
	external/PUPPI/puppiAlgoBin.hh \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h
tmp/modules/StatusPidFilter.$(ObjSuf): \
	modules/StatusPidFilter.$(SrcSuf) \
	modules/StatusPidFilter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/TauTagging.$(ObjSuf): \
	modules/TauTagging.$(SrcSuf) \
	modules/TauTagging.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/TrackPileUpSubtractor.$(ObjSuf): \
	modules/TrackPileUpSubtractor.$(SrcSuf) \
	modules/TrackPileUpSubtractor.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/TreeWriter.$(ObjSuf): \
	modules/TreeWriter.$(SrcSuf) \
	modules/TreeWriter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/modules/UniqueObjectFinder.$(ObjSuf): \
	modules/UniqueObjectFinder.$(SrcSuf) \
	modules/UniqueObjectFinder.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/Weighter.$(ObjSuf): \
	modules/Weighter.$(SrcSuf) \
	modules/Weighter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/HCalorimeter.$(ObjSuf): \
	modules/HCalorimeter.$(SrcSuf) \
	modules/HCalorimeter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/Vertexing.$(ObjSuf): \
	modules/Vertexing.$(SrcSuf) \
	modules/Vertexing.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/RhoMerger.$(ObjSuf): \
	modules/RhoMerger.$(SrcSuf) \
	modules/RhoMerger.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/ECalorimeter.$(ObjSuf): \
	modules/ECalorimeter.$(SrcSuf) \
	modules/ECalorimeter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/external/ExRootAnalysis/ExRootConfReader.$(ObjSuf): \
	external/ExRootAnalysis/ExRootConfReader.$(SrcSuf) \
	external/ExRootAnalysis/ExRootConfReader.h \
	external/tcl/tcl.h
tmp/external/ExRootAnalysis/ExRootFilter.$(ObjSuf): \
	external/ExRootAnalysis/ExRootFilter.$(SrcSuf) \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/external/ExRootAnalysis/ExRootProgressBar.$(ObjSuf): \
	external/ExRootAnalysis/ExRootProgressBar.$(SrcSuf) \
	external/ExRootAnalysis/ExRootProgressBar.h
tmp/external/ExRootAnalysis/ExRootResult.$(ObjSuf): \
	external/ExRootAnalysis/ExRootResult.$(SrcSuf) \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootUtilities.h
tmp/external/ExRootAnalysis/ExRootTask.$(ObjSuf): \
	external/ExRootAnalysis/ExRootTask.$(SrcSuf) \
	external/ExRootAnalysis/ExRootTask.h \
	external/ExRootAnalysis/ExRootConfReader.h
tmp/external/ExRootAnalysis/ExRootTreeBranch.$(ObjSuf): \
	external/ExRootAnalysis/ExRootTreeBranch.$(SrcSuf) \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/external/ExRootAnalysis/ExRootTreeReader.$(ObjSuf): \
	external/ExRootAnalysis/ExRootTreeReader.$(SrcSuf) \
	external/ExRootAnalysis/ExRootTreeReader.h
tmp/external/ExRootAnalysis/ExRootTreeWriter.$(ObjSuf): \
	external/ExRootAnalysis/ExRootTreeWriter.$(SrcSuf) \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/external/ExRootAnalysis/ExRootUtilities.$(ObjSuf): \
	external/ExRootAnalysis/ExRootUtilities.$(SrcSuf) \
	external/ExRootAnalysis/ExRootUtilities.h
tmp/external/PUPPI/puppiCleanContainer.$(ObjSuf): \
	external/PUPPI/puppiCleanContainer.$(SrcSuf)
DELPHES_OBJ +=  \
	tmp/classes/DelphesClasses.$(ObjSuf) \
	tmp/classes/DelphesFactory.$(ObjSuf) \
	tmp/classes/DelphesFormula.$(ObjSuf) \
	tmp/classes/DelphesHepMCReader.$(ObjSuf) \
	tmp/classes/DelphesLHEFReader.$(ObjSuf) \
	tmp/classes/DelphesModule.$(ObjSuf) \
	tmp/classes/DelphesPileUpReader.$(ObjSuf) \
	tmp/classes/DelphesPileUpWriter.$(ObjSuf) \
	tmp/classes/DelphesSTDHEPReader.$(ObjSuf) \
	tmp/classes/DelphesStream.$(ObjSuf) \
	tmp/modules/BTagging.$(ObjSuf) \
	tmp/modules/Cloner.$(ObjSuf) \
	tmp/modules/ConstituentFilter.$(ObjSuf) \
	tmp/modules/Delphes.$(ObjSuf) \
	tmp/modules/Efficiency.$(ObjSuf) \
	tmp/modules/EnergyScale.$(ObjSuf) \
	tmp/modules/EnergySmearing.$(ObjSuf) \
	tmp/modules/ExampleModule.$(ObjSuf) \
	tmp/modules/FakeLepton.$(ObjSuf) \
	tmp/modules/FastJetFinder.$(ObjSuf) \
	tmp/modules/GenBeamSpotFilter.$(ObjSuf) \
	tmp/modules/IsoTrackFilter.$(ObjSuf) \
	tmp/modules/Isolation.$(ObjSuf) \
	tmp/modules/JetFlavourAssociation.$(ObjSuf) \
	tmp/modules/JetPileUpSubtractor.$(ObjSuf) \
	tmp/modules/LeptonDressing.$(ObjSuf) \
	tmp/modules/Merger.$(ObjSuf) \
	tmp/modules/ModifyBeamSpot.$(ObjSuf) \
	tmp/modules/MomentumSmearing.$(ObjSuf) \
	tmp/modules/NeutrinoFilter.$(ObjSuf) \
	tmp/modules/ParticlePropagator.$(ObjSuf) \
	tmp/modules/PileUpJetID.$(ObjSuf) \
	tmp/modules/PileUpMerger.$(ObjSuf) \
	tmp/modules/RunPUPPI.$(ObjSuf) \
	tmp/modules/StatusPidFilter.$(ObjSuf) \
	tmp/modules/TauTagging.$(ObjSuf) \
	tmp/modules/TrackPileUpSubtractor.$(ObjSuf) \
	tmp/modules/TreeWriter.$(ObjSuf) \
	tmp/modules/UniqueObjectFinder.$(ObjSuf) \
	tmp/modules/Weighter.$(ObjSuf) \
	tmp/modules/HCalorimeter.$(ObjSuf) \
	tmp/modules/Vertexing.$(ObjSuf) \
	tmp/modules/RhoMerger.$(ObjSuf) \
	tmp/modules/ECalorimeter.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootConfReader.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootFilter.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootProgressBar.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootResult.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootTask.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootTreeBranch.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootTreeReader.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootTreeWriter.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootUtilities.$(ObjSuf) \
	tmp/external/PUPPI/puppiCleanContainer.$(ObjSuf)

ifeq ($(HAS_PYTHIA8),true)
DELPHES_OBJ +=  \
	tmp/modules/PileUpMergerPythia8.$(ObjSuf)
endif

tmp/display/DelphesCaloData.$(ObjSuf): \
	display/DelphesCaloData.$(SrcSuf) \
	display/DelphesCaloData.h
tmp/display/DelphesDisplay.$(ObjSuf): \
	display/DelphesDisplay.$(SrcSuf) \
	display/DelphesDisplay.h
DISPLAY_OBJ +=  \
	tmp/display/DelphesCaloData.$(ObjSuf) \
	tmp/display/DelphesDisplay.$(ObjSuf)

ifeq ($(HAS_PYTHIA8),true)
DISPLAY_OBJ +=  \
	
endif

tmp/external/tcl/panic.$(ObjSuf): \
	external/tcl/panic.c
tmp/external/tcl/tclAlloc.$(ObjSuf): \
	external/tcl/tclAlloc.c
tmp/external/tcl/tclAsync.$(ObjSuf): \
	external/tcl/tclAsync.c
tmp/external/tcl/tclBasic.$(ObjSuf): \
	external/tcl/tclBasic.c
tmp/external/tcl/tclCkalloc.$(ObjSuf): \
	external/tcl/tclCkalloc.c
tmp/external/tcl/tclCmdAH.$(ObjSuf): \
	external/tcl/tclCmdAH.c
tmp/external/tcl/tclCmdIL.$(ObjSuf): \
	external/tcl/tclCmdIL.c
tmp/external/tcl/tclCmdMZ.$(ObjSuf): \
	external/tcl/tclCmdMZ.c
tmp/external/tcl/tclCompExpr.$(ObjSuf): \
	external/tcl/tclCompExpr.c
tmp/external/tcl/tclCompile.$(ObjSuf): \
	external/tcl/tclCompile.c
tmp/external/tcl/tclExecute.$(ObjSuf): \
	external/tcl/tclExecute.c
tmp/external/tcl/tclGet.$(ObjSuf): \
	external/tcl/tclGet.c
tmp/external/tcl/tclHash.$(ObjSuf): \
	external/tcl/tclHash.c
tmp/external/tcl/tclHistory.$(ObjSuf): \
	external/tcl/tclHistory.c
tmp/external/tcl/tclIndexObj.$(ObjSuf): \
	external/tcl/tclIndexObj.c
tmp/external/tcl/tclLink.$(ObjSuf): \
	external/tcl/tclLink.c
tmp/external/tcl/tclListObj.$(ObjSuf): \
	external/tcl/tclListObj.c
tmp/external/tcl/tclNamesp.$(ObjSuf): \
	external/tcl/tclNamesp.c
tmp/external/tcl/tclObj.$(ObjSuf): \
	external/tcl/tclObj.c
tmp/external/tcl/tclParse.$(ObjSuf): \
	external/tcl/tclParse.c
tmp/external/tcl/tclPosixStr.$(ObjSuf): \
	external/tcl/tclPosixStr.c
tmp/external/tcl/tclPreserve.$(ObjSuf): \
	external/tcl/tclPreserve.c
tmp/external/tcl/tclProc.$(ObjSuf): \
	external/tcl/tclProc.c
tmp/external/tcl/tclResolve.$(ObjSuf): \
	external/tcl/tclResolve.c
tmp/external/tcl/tclStringObj.$(ObjSuf): \
	external/tcl/tclStringObj.c
tmp/external/tcl/tclUtil.$(ObjSuf): \
	external/tcl/tclUtil.c
tmp/external/tcl/tclVar.$(ObjSuf): \
	external/tcl/tclVar.c
TCL_OBJ +=  \
	tmp/external/tcl/panic.$(ObjSuf) \
	tmp/external/tcl/tclAlloc.$(ObjSuf) \
	tmp/external/tcl/tclAsync.$(ObjSuf) \
	tmp/external/tcl/tclBasic.$(ObjSuf) \
	tmp/external/tcl/tclCkalloc.$(ObjSuf) \
	tmp/external/tcl/tclCmdAH.$(ObjSuf) \
	tmp/external/tcl/tclCmdIL.$(ObjSuf) \
	tmp/external/tcl/tclCmdMZ.$(ObjSuf) \
	tmp/external/tcl/tclCompExpr.$(ObjSuf) \
	tmp/external/tcl/tclCompile.$(ObjSuf) \
	tmp/external/tcl/tclExecute.$(ObjSuf) \
	tmp/external/tcl/tclGet.$(ObjSuf) \
	tmp/external/tcl/tclHash.$(ObjSuf) \
	tmp/external/tcl/tclHistory.$(ObjSuf) \
	tmp/external/tcl/tclIndexObj.$(ObjSuf) \
	tmp/external/tcl/tclLink.$(ObjSuf) \
	tmp/external/tcl/tclListObj.$(ObjSuf) \
	tmp/external/tcl/tclNamesp.$(ObjSuf) \
	tmp/external/tcl/tclObj.$(ObjSuf) \
	tmp/external/tcl/tclParse.$(ObjSuf) \
	tmp/external/tcl/tclPosixStr.$(ObjSuf) \
	tmp/external/tcl/tclPreserve.$(ObjSuf) \
	tmp/external/tcl/tclProc.$(ObjSuf) \
	tmp/external/tcl/tclResolve.$(ObjSuf) \
	tmp/external/tcl/tclStringObj.$(ObjSuf) \
	tmp/external/tcl/tclUtil.$(ObjSuf) \
	tmp/external/tcl/tclVar.$(ObjSuf)

modules/ModifyBeamSpot.h: \
	classes/DelphesModule.h \
	modules/simpleVariableCollector.h
	@touch $@

modules/EnergySmearing.h: \
	classes/DelphesModule.h
	@touch $@

modules/LeptonDressing.h: \
	classes/DelphesModule.h
	@touch $@

modules/NeutrinoFilter.h: \
	classes/DelphesModule.h
	@touch $@

modules/ConstituentFilter.h: \
	classes/DelphesModule.h
	@touch $@

classes/DelphesModule.h: \
	external/ExRootAnalysis/ExRootTask.h
	@touch $@

modules/Isolation.h: \
	classes/DelphesModule.h
	@touch $@

modules/EnergyScale.h: \
	classes/DelphesModule.h
	@touch $@

modules/Merger.h: \
	classes/DelphesModule.h
	@touch $@

modules/ExampleModule.h: \
	classes/DelphesModule.h
	@touch $@

modules/JetPileUpSubtractor.h: \
	classes/DelphesModule.h
	@touch $@

modules/Vertexing.h: \
	modules/simpleVariableCollector.h \
	classes/DelphesModule.h
	@touch $@

modules/Efficiency.h: \
	classes/DelphesModule.h
	@touch $@

modules/TrackPileUpSubtractor.h: \
	classes/DelphesModule.h \
	modules/simpleVariableCollector.h
	@touch $@

modules/IsoTrackFilter.h: \
	classes/DelphesModule.h
	@touch $@

modules/HCalorimeter.h: \
	classes/DelphesModule.h \
	modules/simpleVariableCollector.h
	@touch $@

modules/FakeLepton.h: \
	classes/DelphesModule.h
	@touch $@

modules/PileUpMerger.h: \
	classes/DelphesModule.h \
	modules/simpleVariableCollector.h
	@touch $@

modules/RhoMerger.h: \
	classes/DelphesModule.h
	@touch $@

modules/Cloner.h: \
	classes/DelphesModule.h
	@touch $@

modules/RunPUPPI.h: \
	classes/DelphesModule.h
	@touch $@

modules/ECalorimeter.h: \
	classes/DelphesModule.h \
	modules/simpleVariableCollector.h
	@touch $@

modules/PileUpJetID.h: \
	classes/DelphesModule.h \
	modules/simpleVariableCollector.h
	@touch $@

modules/MomentumSmearing.h: \
	classes/DelphesModule.h
	@touch $@

modules/TauTagging.h: \
	classes/DelphesModule.h
	@touch $@

modules/Delphes.h: \
	classes/DelphesModule.h
	@touch $@

modules/UniqueObjectFinder.h: \
	classes/DelphesModule.h
	@touch $@

modules/JetFlavourAssociation.h: \
	classes/DelphesModule.h \
	classes/DelphesClasses.h
	@touch $@

modules/PileUpMergerPythia8.h: \
	classes/DelphesModule.h
	@touch $@

modules/ParticlePropagator.h: \
	classes/DelphesModule.h \
	modules/simpleVariableCollector.h
	@touch $@

external/PUPPI/puppiCleanContainer.hh: \
	external/PUPPI/RecoObj.hh \
	external/PUPPI/puppiParticle.hh \
	external/PUPPI/puppiAlgoBin.hh
	@touch $@

modules/BTagging.h: \
	classes/DelphesModule.h \
	classes/DelphesClasses.h
	@touch $@

modules/GenBeamSpotFilter.h: \
	classes/DelphesModule.h
	@touch $@

modules/Weighter.h: \
	classes/DelphesModule.h
	@touch $@

external/ExRootAnalysis/ExRootTask.h: \
	external/ExRootAnalysis/ExRootConfReader.h
	@touch $@

modules/TreeWriter.h: \
	classes/DelphesModule.h
	@touch $@

modules/StatusPidFilter.h: \
	classes/DelphesModule.h
	@touch $@

classes/DelphesClasses.h: \
	classes/SortableObject.h
	@touch $@

modules/FastJetFinder.h: \
	classes/DelphesModule.h \
	modules/simpleVariableCollector.h
	@touch $@



###

all: $(DELPHES) $(EXECUTABLE) genMinBias_14TeV countEvents

display: $(DISPLAY)

$(DELPHES): $(DELPHES_DICT_OBJ) $(DELPHES_OBJ) $(TCL_OBJ)
	@mkdir -p $(@D)
	@echo ">> Building $@"
ifeq ($(ARCH),aix5)
	@$(MAKESHARED) $(OutPutOpt) $@ $(DELPHES_LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	@$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@ $(DELPHES_LIBS)
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
	@ln -sf $@ $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
	@bindexplib $* $^ > $*.def
	@lib -nologo -MACHINE:IX86 $^ -def:$*.def $(OutPutOpt)$(DELPHESLIB)
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(DELPHES_LIBS) $(OutPutOpt)$@
	@$(MT_DLL)
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(DELPHES_LIBS)
	@$(MT_DLL)
endif
endif
endif

$(DISPLAY): $(DELPHES_DICT_OBJ) $(DISPLAY_DICT_OBJ) $(DELPHES_OBJ) $(DISPLAY_OBJ) $(TCL_OBJ)
	@mkdir -p $(@D)
	@echo ">> Building $@"
ifeq ($(ARCH),aix5)
	@$(MAKESHARED) $(OutPutOpt) $@ $(DISPLAY_LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	@$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@ $(DISPLAY_LIBS)
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
	@ln -sf $@ $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
	@bindexplib $* $^ > $*.def
	@lib -nologo -MACHINE:IX86 $^ -def:$*.def $(OutPutOpt)$(DISPLAYLIB)
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(DISPLAY_LIBS) $(OutPutOpt)$@
	@$(MT_DLL)
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(DISPLAY_LIBS)
	@$(MT_DLL)
endif
endif
endif

clean:
	@rm -f $(DELPHES_DICT_OBJ) $(DISPLAY_DICT_OBJ) $(DELPHES_OBJ) $(DISPLAY_OBJ) $(TCL_OBJ) core
	@rm -rf tmp

distclean: clean
	@rm -f $(DELPHES) $(DELPHESLIB) $(DISPLAY) $(DISPLAYLIB) $(EXECUTABLE)

dist:
	@echo ">> Building $(DISTTAR)"
	@mkdir -p $(DISTDIR)
	@cp -a CREDITS README VERSION Makefile configure classes converters display doc examples external modules python readers $(DISTDIR)
	@find $(DISTDIR) -depth -name .\* -exec rm -rf {} \;
	@tar -czf $(DISTTAR) $(DISTDIR)
	@rm -rf $(DISTDIR)

###

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

%Dict.$(SrcSuf):
	@mkdir -p $(@D)
	@echo ">> Generating $@"
	@rootcint -f $@ -c -Iexternal $<
	@echo "#define private public" > $@.arch
	@echo "#define protected public" >> $@.arch
	@mv $@ $@.base
	@cat $@.arch $< $@.base > $@
	@rm $@.arch $@.base

$(DELPHES_OBJ): tmp/%.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(DISPLAY_OBJ): tmp/%.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(DELPHES_DICT_OBJ): %.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(DISPLAY_DICT_OBJ): %.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(TCL_OBJ): tmp/%.$(ObjSuf): %.c
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@gcc $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(EXECUTABLE_OBJ): tmp/%.$(ObjSuf): %.cpp
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(EXECUTABLE): %$(ExeSuf): $(DELPHES_DICT_OBJ) $(DELPHES_OBJ) $(TCL_OBJ)
	@echo ">> Building $@"
	@$(LD) $(LDFLAGS) $^ $(DELPHES_LIBS) $(OutPutOpt)$@

genMinBias_14TeV: external/MinBiasProduction/genMinBias_14TeV.cpp
	@echo ">> Compiling $<"
	@$(CXX) -o $@ $< -I$(HEPMC)/include -L$(HEPMC)/lib -I$(PYTHIA8DATA)/../../../include -L$(PYTHIA8DATA)/../../../lib -I$(LHAPDF)/include -L$(LHAPDF)/lib -lHepMC  -lpythia8 -lLHAPDF -lgfortran -lpythia8lhapdf5

countEvents: external/LHEActions/countEvents.cpp
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -o $@ $< 



###


