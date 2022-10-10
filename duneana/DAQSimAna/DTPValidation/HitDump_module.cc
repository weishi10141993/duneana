////////////////////////////////////////////////////////////////////////
// Class:       HitDump
// Module Type: analyzer
// File:        HitDump_module.cc
//
// Dump the three-view reob::Hit objects to a txt file for quick validation 
// Can pass either software-generated or firmware tps. 
////////////////////////////////////////////////////////////////////////


#include <fstream>

// Framework includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

class HitDump : public art::EDAnalyzer {

public:

  explicit HitDump(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  HitDump(HitDump const &) = delete;
  HitDump(HitDump &&) = delete;
  HitDump & operator = (HitDump const &) = delete;
  HitDump & operator = (HitDump &&) = delete;

  void analyze(art::Event const & evt) override;

private:

  std::string m_HitLabel;

  //files where we dump our TPs
  std::string m_outputFilename;
  std::ofstream m_outputFile; 
 
  // --- Declare our services
  //  art::ServiceHandle<geo::Geometry> geo;
  
};

//......................................................
HitDump::HitDump(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p), 
  m_HitLabel(p.get<std::string>("HitLabel")),
  m_outputFilename(p.get<std::string>("OutputFile")),
  m_outputFile(m_outputFilename)
{
  
}

//......................................................
void HitDump::analyze(art::Event const & evt)
{


  // --- Lift out the reco hits:
  auto reco_hits = evt.getValidHandle<std::vector<recob::Hit> >(m_HitLabel);
 
  for(size_t hit = 0; hit < reco_hits->size(); ++hit) {
    recob::Hit const& ThisHit = reco_hits->at(hit);   // current hit 

    m_outputFile << ThisHit.StartTick()  << " " << ThisHit.EndTick() - ThisHit.StartTick() 
		 << " " << ThisHit.PeakTime()  << " " << ThisHit.Channel()
		 << " " << ThisHit.SummedADC() << " "<< ThisHit.PeakAmplitude() << std::endl; 
  } // Loop over reco_hits.
} // Analyze HitDump.


//......................................................
DEFINE_ART_MODULE(HitDump)
