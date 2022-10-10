/////////////////////////////////////////////////////////////////////////
// Class:       TaggedHitDump
// Module Type: analyzer
// File:        TaggedHitDump_module.cc
// Author:      Klaudia Wawrowska
//
// Add noise to the waveforms and dump the three-view hit information 
// to a txt file in a format suitable for TriggerAlgs
// tag each hit to determine whether it was a signal/noise/bgd hit. 
////////////////////////////////////////////////////////////////////////


#include <fstream>
#include <iterator> 


#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
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
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"


// kUnknown == noise 
enum PType{ kUnknown=0, kGen, kAPA, kCPA, kAr39, kNeut, kKryp, kPol, kRdn, kAr42};

class TaggedHitDump : public art::EDAnalyzer {

public:

  explicit TaggedHitDump(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  TaggedHitDump(TaggedHitDump const &) = delete;
  TaggedHitDump(TaggedHitDump &&) = delete;
  TaggedHitDump & operator = (TaggedHitDump const &) = delete;
  TaggedHitDump & operator = (TaggedHitDump &&) = delete;

  void analyze(art::Event const & evt) override;

private:

  //--- any custom functions 
  static void BuildMCParticleHitMaps( const art::Event& evt, 
				      const lar_pandora::HitVector& hitVector,
				      const lar_pandora::SimChannelVector& simChannelVector, 
				      lar_pandora::HitsToTrackIDEs& hitsToTrackIDEs,
				      bool VDCBgeo); // need larger hit window if working with the messed up VDCB geometry
  //--- producer labels
  //generator used for signal events 
  std::string m_GenLabel; // typically single particle gun or some neutrino generator
  std::string m_GeantLabel; //g4
  std::string m_SimChanLabel; 
  std::string m_HitLabel;
  std::string m_RawDigitLabel;
  //--- background producer labels
  std::string m_Ar39Gen; 
  std::string m_Ar42Gen; 
  std::string m_APAGen; //Co60 in APA frames
  std::string m_NeutronGen; // neutrons from concrete 
  std::string m_CPAGen; // K40 from CPA
  std::string m_Kr85Gen; //Kr in LAr 
  std::string m_Rn222Gen; //Rn in LAr 
  std::string m_Po210Gen; //polonium-210 - alpha emitter

  //--- files where we dump our collection and induction TPs
  std::string m_outputFilename;
  std::ofstream m_outputFile; 
  std::string m_outputFilename1;
  std::ofstream m_outputFile1; 

  //--- VD coldbox geo mods
  bool m_VDCBgeo;
  

  // --- Declare our services
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  
};

//......................................................
TaggedHitDump::TaggedHitDump(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p), 
  m_GenLabel(      p.get<std::string>("GenLabel")),
  m_GeantLabel(    p.get<std::string>("GEANT4Label")),
  m_SimChanLabel(  p.get<std::string>("SimChanLabel")),
  m_HitLabel(      p.get<std::string>("HitLabel")),
  m_RawDigitLabel( p.get<std::string>("RawDigitLabel")),
  m_Ar39Gen(       p.get<std::string>("Argon39Label")),
  m_Ar42Gen(       p.get<std::string>("Argon42Label")),
  m_APAGen(        p.get<std::string>("APALabel")),
  m_NeutronGen(    p.get<std::string>("NeutronLabel")),
  m_CPAGen(        p.get<std::string>("CPALabel")),
  m_Kr85Gen(       p.get<std::string>("KryptonLabel")),
  m_Rn222Gen(      p.get<std::string>("RadonLabel")),
  m_Po210Gen(      p.get<std::string>("PoloniumLabel")),
  m_outputFilename(p.get<std::string>("OutputFile")),
  m_outputFile(m_outputFilename),
  m_outputFilename1(p.get<std::string>("OutputFile1")),
  m_outputFile1(m_outputFilename1),
  m_VDCBgeo(       p.get<bool>("VDCBgeo",false))
{
  
}

// build mapping between Hits and TrackIDEs
void TaggedHitDump::BuildMCParticleHitMaps(const art::Event& evt, 
						const lar_pandora::HitVector& hitVector,
						const lar_pandora::SimChannelVector& simChannelVector, 
						lar_pandora::HitsToTrackIDEs& hitsToTrackIDEs,
						bool VDCBgeo)
{
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  lar_pandora::SimChannelMap simChannelMap; 
  for (lar_pandora::SimChannelVector::const_iterator iter = simChannelVector.begin(), 
	                                             iterEnd = simChannelVector.end();
       iter != iterEnd;  ++iter) {
    const art::Ptr<sim::SimChannel> simChannel = *iter; 
    simChannelMap.insert(lar_pandora::SimChannelMap::value_type(simChannel->Channel(), simChannel));
  }

  for (lar_pandora::HitVector::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end();
       iter != iterEnd; ++iter) {
    const art::Ptr<recob::Hit> hit = *iter;

    lar_pandora::SimChannelMap::const_iterator sIter = simChannelMap.find(hit->Channel());
    if (simChannelMap.end() == sIter) continue; //Hit has no truth information, move on 
    
    //Convert TDCtick (int) to TDC (unsigned int) before passing to simChannel 
    const raw::TDCtick_t start_tick(clock_data.TPCTick2TDC( hit->PeakTimeMinusRMS() ));
    const raw::TDCtick_t end_tick(clock_data.TPCTick2TDC( hit->PeakTimePlusRMS() ));
    const unsigned int start_tdc((start_tick < 0) ? 0 : start_tick);
    const unsigned int end_tdc(end_tick);
    if (start_tdc > end_tdc) continue; // hit undershoots the readout window 

    unsigned int t_delay = 20; 
    if (VDCBgeo) t_delay = 300; //for some reason the IDE matching in VDCB behaves badly so need bigger window

    const art::Ptr<sim::SimChannel> simChannel = sIter->second;
    const lar_pandora::TrackIDEVector trackCollection(simChannel->TrackIDEs(start_tdc-t_delay, end_tdc+t_delay));
    //increase hit time window as induction hits are often mistagged due to increased TOT span
 
    if (trackCollection.empty()) continue; 


    for (unsigned int iTrack = 0, iTrackEnd = trackCollection.size(); iTrack < iTrackEnd;  ++iTrack){
      const sim::TrackIDE trackIDE = trackCollection.at(iTrack);
      hitsToTrackIDEs[hit].push_back(trackIDE);
    }
  }
}
//......................................................
void TaggedHitDump::analyze(art::Event const & evt)
{

  std::map< int, PType> TrkIdToPType;

  // --- First want to associate all the g4 tracks to their respective generators:

  auto GenTrue = evt.getValidHandle< std::vector<simb::MCTruth> >(m_GenLabel);
  if (GenTrue.isValid()){
    art::FindManyP<simb::MCParticle> GenParts( GenTrue, evt, m_GeantLabel);
    for (size_t i = 0; i < GenTrue->size(); i ++){
      for (size_t j = 0; j < GenParts.at(i).size(); j++){
	simb::MCParticle ThisParticle = (*GenParts.at(i).at(j));
	TrkIdToPType.insert( std::pair< int, PType >( ThisParticle.TrackId(), kGen ));
      }
    }
  }

  auto Ar39True = evt.getValidHandle< std::vector<simb::MCTruth> >(m_Ar39Gen);
  if (Ar39True.isValid()){
    art::FindManyP<simb::MCParticle> Ar39Parts(Ar39True, evt, m_GeantLabel);
    for (size_t i = 0; i < Ar39True->size(); i++){
      for (size_t j = 0; j < Ar39Parts.at(i).size(); j++){
	simb::MCParticle ThisParticle = (*Ar39Parts.at(i).at(j));
	TrkIdToPType.insert( std::pair< int, PType>( ThisParticle.TrackId(), kAr39));
      }
    }
  }

  auto Ar42True = evt.getValidHandle< std::vector<simb::MCTruth> >(m_Ar42Gen);
  if (Ar42True.isValid()){
    art::FindManyP<simb::MCParticle> Ar42Parts(Ar42True, evt, m_GeantLabel);
    for (size_t i = 0; i < Ar42True->size(); i++){
      for (size_t j = 0; j < Ar42Parts.at(i).size(); j++){
	simb::MCParticle ThisParticle = (*Ar42Parts.at(i).at(j));
	TrkIdToPType.insert( std::pair< int, PType>( ThisParticle.TrackId(), kAr42));
      }
    }
  }
  
  auto APATrue = evt.getValidHandle< std::vector<simb::MCTruth> >(m_APAGen);
  if (APATrue.isValid()){
    art::FindManyP<simb::MCParticle> APAParts(APATrue, evt, m_GeantLabel);
    for (size_t i = 0; i < APATrue->size(); i++){
      for (size_t j = 0; j < APAParts.at(i).size(); j++){
	simb::MCParticle ThisParticle = (*APAParts.at(i).at(j));
	TrkIdToPType.insert( std::pair< int, PType>( ThisParticle.TrackId(), kAPA));
      }
    }
  }
  
  auto CPATrue = evt.getValidHandle< std::vector<simb::MCTruth> >(m_CPAGen);
  if (CPATrue.isValid()){
    art::FindManyP<simb::MCParticle> CPAParts(CPATrue, evt, m_GeantLabel);
    for (size_t i = 0; i < CPATrue->size(); i++){
      for (size_t j = 0; j < CPAParts.at(i).size(); j++){
	simb::MCParticle ThisParticle = (*CPAParts.at(i).at(j));	
	TrkIdToPType.insert( std::pair< int, PType>( ThisParticle.TrackId(), kCPA));
      }   
    }
  }
  
  auto NeutTrue = evt.getValidHandle< std::vector<simb::MCTruth> >(m_NeutronGen);
  if (NeutTrue.isValid()){
    art::FindManyP<simb::MCParticle> NeutParts(NeutTrue, evt, m_GeantLabel);
    for (size_t i = 0; i < NeutTrue->size(); i++){
      for (size_t j = 0; j < NeutParts.at(i).size(); j++){
	simb::MCParticle ThisParticle = (*NeutParts.at(i).at(j));
	TrkIdToPType.insert( std::pair< int, PType>( ThisParticle.TrackId(), kNeut));
      }   
    }
  }

  auto KrypTrue = evt.getValidHandle< std::vector<simb::MCTruth> >(m_Kr85Gen);
  if (KrypTrue.isValid()){
    art::FindManyP<simb::MCParticle> KrypParts(KrypTrue, evt, m_GeantLabel);
    for (size_t i = 0; i < KrypTrue->size(); i++){
      for (size_t j = 0; j < KrypParts.at(i).size(); j++){
	simb::MCParticle ThisParticle = (*KrypParts.at(i).at(j));
	TrkIdToPType.insert( std::pair< int, PType>( ThisParticle.TrackId(), kKryp));
      }    
    }
  }

  auto PolTrue = evt.getValidHandle< std::vector<simb::MCTruth> >(m_Po210Gen);
  if (PolTrue.isValid()){
    art::FindManyP<simb::MCParticle> PolParts(PolTrue, evt, m_GeantLabel);
    for (size_t i = 0; i < PolTrue->size(); i++){
      for (size_t j = 0; j < PolParts.at(i).size(); j++){
	simb::MCParticle ThisParticle = (*PolParts.at(i).at(j));
	TrkIdToPType.insert( std::pair< int, PType>( ThisParticle.TrackId(),kPol));
      }   
    }
  }

  auto RdnTrue = evt.getValidHandle< std::vector<simb::MCTruth> >(m_Rn222Gen);
  if (RdnTrue.isValid()){
    art::FindManyP<simb::MCParticle> RdnParts(RdnTrue, evt, m_GeantLabel);
    for (size_t i = 0; i < RdnTrue->size(); i++){
      for (size_t j = 0; j < RdnParts.at(i).size(); j++){
	simb::MCParticle ThisParticle = (*RdnParts.at(i).at(j));
	TrkIdToPType.insert( std::pair< int,PType>(ThisParticle.TrackId(), kRdn));
      }   
    }
  }
  
  //Check that everything works as expected
  /*
  for (auto it : TrkIdToPType){
    m_outputFile1 << evt.event() << " " << it.first << " " << it.second << std::endl;
  }
  */

  //Get the reco hits
  std::vector< art::Ptr<recob::Hit> > recoHits;
  auto HitHandle = evt.getValidHandle <std::vector<recob::Hit> >(m_HitLabel);
    if (HitHandle){    art::fill_ptr_vector(recoHits, HitHandle);   }

  if (HitHandle){

    //Get the SimChannels
    bool areSimChannelsValid; 
    lar_pandora::SimChannelVector simChannels;
    lar_pandora::LArPandoraHelper::CollectSimChannels( evt, m_SimChanLabel, simChannels, areSimChannelsValid) ; 

    //Make association maps between hits and trackIDEs
    lar_pandora::HitsToTrackIDEs hitsToTrackIDEs;
    
    BuildMCParticleHitMaps( evt, recoHits, simChannels, hitsToTrackIDEs, m_VDCBgeo);

    for (unsigned int hit = 0; hit < recoHits.size(); ++hit) {
       
      //Which GenParticle (if any) produced this hit?
      const art::Ptr<recob::Hit> ThisHit = recoHits.at(hit);
      PType ThisPType = kUnknown; 
      
      lar_pandora::HitsToTrackIDEs::const_iterator uIter = hitsToTrackIDEs.find( ThisHit ); 
      if (hitsToTrackIDEs.end() != uIter){
	const lar_pandora::TrackIDEVector& trackCollection = uIter->second;
	
	int bestTrackID = -1;
	float bestEnergyFrac(0.f);
	for (lar_pandora::TrackIDEVector::const_iterator 
	       iter2 = trackCollection.begin(),
	       iterEnd2 = trackCollection.end();
	     iter2 != iterEnd2; 
	     ++iter2){
	  const sim::TrackIDE& trackIDE = *iter2; 
	  const int trackID( std::abs( trackIDE.trackID ));
	  const float energyFrac(trackIDE.energyFrac); 
	  if (energyFrac > bestEnergyFrac){
	    bestEnergyFrac = energyFrac; 
	    bestTrackID = trackID; 
	  }
	}
	
	if (bestTrackID >= 0){
	  for (auto& x: TrkIdToPType){
	    if (x.first == bestTrackID) ThisPType = x.second; 
	  }
	}
      }
    

      //Output TP information to a file 
      int increment = 4492; //readout timewindow
      int choffset = 0;
      int rate =1;  
      int timeoffset = 0; 
      
      if (m_VDCBgeo){ 
	choffset=1600;  //add a channel offset if dealing with the VDCB geometry
	increment = 6000; 
	rate = 25; //change sampling frequency 2 MHz-> 50 MHz for DUNEDAQ
	timeoffset =  (evt.event() - 1) * increment; //using DAQ-like "timestamps" instead of event numbers to tell different events apart
      
      }
      
    
      //--- dump the collection hits to a file  
      //if (ThisHit->View() == 2){
      if (ThisHit->LocalIndex() == 0){
	m_outputFile << (ThisHit->StartTick() + timeoffset)*rate << " " << (ThisHit->EndTick() - ThisHit->StartTick())*rate 
		     << " " << (ThisHit->PeakTime() + timeoffset)*rate << " " << ThisHit->Channel() + choffset 
		     << " " << ThisHit->SummedADC() << " "<< ThisHit->PeakAmplitude() << " "
		     << ThisHit->WireID().TPC  << " " << "1" << " " << ThisPType <<std::endl; 
      }
      
      //--- dump the induction hits to a file
      //if (ThisHit->View() == 0 || ThisHit->View() == 1){
      if (ThisHit->LocalIndex() == 1){      
	m_outputFile1 << (ThisHit->StartTick() + timeoffset)*rate << " " << (ThisHit->EndTick() - ThisHit->StartTick())*rate 
		      << " " << (ThisHit->PeakTime() + timeoffset)*rate << " " << ThisHit->Channel() + choffset 
		      << " " << ThisHit->SummedADC() << " "<< ThisHit->PeakAmplitude() << " "
		      << ThisHit->WireID().TPC  << " " << "1" << " " << ThisPType <<  std::endl; 
      }

    }// Loop over recoHits.
  } // if HitHandleValide
}// Analyze TaggedHitDump.
 

//......................................................
DEFINE_ART_MODULE(TaggedHitDump)
