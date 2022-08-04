////////////////////////////////////////////////////////////////////////
// Class:       DAQTPFinderDump
// Module Type: analyzer
// File:        DAQTPFinderDump_module.cc
//
// Add noise to the waveforms and dump the three-view hit information 
// to a txt file in a format suitable for TriggerAlgs
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

class DAQTPFinderDump : public art::EDAnalyzer {

public:

  explicit DAQTPFinderDump(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  DAQTPFinderDump(DAQTPFinderDump const &) = delete;
  DAQTPFinderDump(DAQTPFinderDump &&) = delete;
  DAQTPFinderDump & operator = (DAQTPFinderDump const &) = delete;
  DAQTPFinderDump & operator = (DAQTPFinderDump &&) = delete;

  void analyze(art::Event const & evt) override;

private:

  std::string m_HitLabel;
  std::string m_RawDigitLabel;

  //files where we dump our TPs
  std::string m_outputFilename;
  std::ofstream m_outputFile; 
  std::string m_outputFilename1;
  std::ofstream m_outputFile1; 
  bool m_VDCBgeo;

  // --- Declare our services
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  
};

//......................................................
DAQTPFinderDump::DAQTPFinderDump(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p), 
  m_HitLabel(p.get<std::string>("HitLabel")),
  m_RawDigitLabel(p.get<std::string>("RawDigitLabel")),
  m_outputFilename(p.get<std::string>("OutputFile")),
  m_outputFile(m_outputFilename),
  m_outputFilename1(p.get<std::string>("OutputFile1")),
  m_outputFile1(m_outputFilename1),
  m_VDCBgeo(p.get<bool>("VDCBgeo",false))
{
  
}

//......................................................
void DAQTPFinderDump::analyze(art::Event const & evt)
{
  // --- Lift out the reco hits:
  auto reco_hits = evt.getValidHandle<std::vector<recob::Hit> >(m_HitLabel);
 
  //auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

  int choffset = 0; //add a channel offset if dealing with the VDCB geometry
  int rate = 25; // change sampling frequency from 2 MHz to 50 MHz for dunedaq 
  int increment = 4492; //readout timewindow
  if (m_VDCBgeo){ 
    choffset=1600;
    increment = 6000; 
  }

  //like in DAQ, using "timestamps" instead of event numbers to tell different events apart 
  int timeoffset =  (evt.event() - 1) * increment;
  for(size_t hit = 0; hit < reco_hits->size(); ++hit) {
    recob::Hit const& ThisHit = reco_hits->at(hit);   // current hit 

    //dump the collection hits to a file 
    //if (ThisHit.View() == 2){
    if (ThisHit.LocalIndex() == 0){
      m_outputFile << (ThisHit.StartTick() + timeoffset)*rate << " " << (ThisHit.EndTick() - ThisHit.StartTick())*rate 
		   << " " << (ThisHit.PeakTime() + timeoffset)*rate << " " << ThisHit.Channel() + choffset 
		   << " " << ThisHit.SummedADC() << " "<< ThisHit.PeakAmplitude() << " "
		   << ThisHit.WireID().TPC  << " " << "1" << std::endl; 
    }
    
    //dump the induction hits to a file
    //if (ThisHit.View() == 0 || ThisHit.View() == 1){
    if (ThisHit.LocalIndex() == 1){      
      m_outputFile1 << (ThisHit.StartTick() + timeoffset)*rate << " " << (ThisHit.EndTick() - ThisHit.StartTick())*rate 
		    << " " << (ThisHit.PeakTime() + timeoffset)*rate << " " << ThisHit.Channel() + choffset 
		    << " " << ThisHit.SummedADC() << " "<< ThisHit.PeakAmplitude() << " "
		    << ThisHit.WireID().TPC  << " " << "1" << std::endl; 
    }
  } // Loop over reco_hits.
} // Analyze DAQTPFinderDump.


//......................................................
DEFINE_ART_MODULE(DAQTPFinderDump)
