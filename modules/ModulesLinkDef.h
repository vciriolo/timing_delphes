
/** \class
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  $Date: 2013-07-04 13:19:15 +0200 (Thu, 04 Jul 2013) $
 *  $Revision: 1177 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Delphes.h"

#include "modules/FastJetFinder.h"
#include "modules/ParticlePropagator.h"
#include "modules/Efficiency.h"
#include "modules/EnergySmearing.h"
#include "modules/MomentumSmearing.h"
//#include "modules/Calorimeter.h"
#include "modules/Isolation.h"
#include "modules/IsoTrackFilter.h"
#include "modules/EnergyScale.h"
#include "modules/UniqueObjectFinder.h"
#include "modules/BTagging.h"
#include "modules/TauTagging.h"
#include "modules/TreeWriter.h"
#include "modules/Merger.h"
#include "modules/LeptonDressing.h"
#include "modules/PileUpMerger.h"
#include "modules/JetPileUpSubtractor.h"
#include "modules/TrackPileUpSubtractor.h"
#include "modules/ConstituentFilter.h"
#include "modules/StatusPidFilter.h"
#include "modules/Cloner.h"
#include "modules/Weighter.h"
#include "modules/ExampleModule.h"
#include "modules/JetFlavourAssociation.h"
//#include "modules/CleansedJetFinder.h"
#include "modules/PileUpJetID.h"
#include "modules/ModifyBeamSpot.h"
#include "modules/GenBeamSpotFilter.h"
#include "modules/RunPUPPI.h"
#include "modules/NeutrinoFilter.h"
#include "modules/FakeLepton.h"
#include "modules/Vertexing.h"
#include "modules/ECalorimeter.h"
#include "modules/HCalorimeter.h"
#include "modules/RhoMerger.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class Delphes+;

#pragma link C++ class FastJetFinder+;
#pragma link C++ class ParticlePropagator+;
#pragma link C++ class Efficiency+;
#pragma link C++ class EnergySmearing+;
#pragma link C++ class MomentumSmearing+;
//#pragma link C++ class Calorimeter+;
#pragma link C++ class Isolation+;
#pragma link C++ class IsoTrackFilter+;
#pragma link C++ class EnergyScale+;
#pragma link C++ class UniqueObjectFinder+;
#pragma link C++ class BTagging+;
#pragma link C++ class TauTagging+;
#pragma link C++ class TreeWriter+;
#pragma link C++ class Merger+;
#pragma link C++ class LeptonDressing+;
#pragma link C++ class PileUpMerger+;
#pragma link C++ class JetPileUpSubtractor+;
#pragma link C++ class TrackPileUpSubtractor+;
#pragma link C++ class JetFlavourAssociation+;
#pragma link C++ class ConstituentFilter+;
#pragma link C++ class StatusPidFilter+;
#pragma link C++ class Cloner+;
#pragma link C++ class Weighter+;
#pragma link C++ class ExampleModule+;
#pragma link C++ class FakeLepton+;
#pragma link C++ class std::vector<TLorentzVector>+;

//#pragma link C++ class CleansedJetFinder+;
#pragma link C++ class PileUpJetID+;
#pragma link C++ class ModifyBeamSpot+;
#pragma link C++ class GenBeamSpotFilter+;
#pragma link C++ class RunPUPPI+;
#pragma link C++ class NeutrinoFilter+;
#pragma link C++ class Vertexing+;
#pragma link C++ class ECalorimeter+;
#pragma link C++ class HCalorimeter+;
#pragma link C++ class RhoMerger+;

#endif
