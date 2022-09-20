////////////////////////////////////////////////////////////////////////
// Class:      ROIAna 
// Plugin Type: analyzer (art v3_00_00)
// File:        ROIAna_module.cc
// Written by Tejin Cai
// Reach out for questions/issues/bugs
////////////////////////////////////////////////////////////////////////
#ifndef WIREANA_H
#define WIREANA_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "lardataobj/RawData/RDTimeStamp.h"


#include "cetlib/search_path.h"
#include "cetlib/filesystem.h"

#include "art_root_io/TFileService.h"

// ROOT includes
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TROOT.h"
#include "TLorentzVector.h"

//Others
#define DEFAULT_VALUE -99999

namespace roiana {

  class ROIAna;

}

class roiana::ROIAna : public art::EDAnalyzer {
public:
  explicit ROIAna(fhicl::ParameterSet const& pset);
  ROIAna(ROIAna const&) = delete;
  ROIAna(ROIAna&&) = delete;
  ROIAna& operator=(ROIAna const&) = delete;
  ROIAna& operator=(ROIAna&&) = delete;
  virtual ~ROIAna() noexcept {};

  /////////////////////////////////////////////
  // Required functions.
  void analyze(art::Event const& evt) override;

  /////////////////////////////////////////////
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  /////////////////////////////////////////////
  std::map<int,std::vector<std::pair<int,int>>>
    ProcessROI(std::vector<art::Ptr<raw::RawDigit>>& rawList);

  bool ispeak(float x) { return (x != 0); }
  bool isZero(float x) { return x == 0 ; }
private:
  ////////////////////////////////////////////
  // Labels
  const art::InputTag fWireProducerLabel; 
  const art::InputTag fRawProducerLabel;
  const art::InputTag fSimChannelLabel;
  const art::InputTag fSimulationProducerLabel;
  ////////////////////////////////////////////
  // Log Control
  int fLogLevel;


  /////////////////////////////////////////////
  // Geometry Options && Tool options
  int fNPlanes;
  int fNChanPerApa;
  unsigned int fNTicksPerWire;

  int fROI_Peak;
  int fROI_Range;
  int fROI_CH;

  std::map<raw::ChannelID_t, std::pair<art::Ptr<raw::RawDigit>, art::Ptr<sim::SimChannel>>> 
    ch_w_sc;
  std::map<int, std::string> 
    trkid_to_label_map;

  /////////////////////////////////////////////
  // Geometry services
  const geo::GeometryCore* geo; //= lar::providerFrom<geo::Geometry>();

  int run;
  int subrun;
  int event;
  int MC;

  
  /////////////////////////////////////////////
  // Private Functions
  template<class T>
  void SortWirePtrByChannel( std::vector<art::Ptr<T>> &vec, bool increasing );

  
  /////////////////////////////////////////////
  // Truth Filter
  art::ServiceHandle<cheat::ParticleInventoryService> PIS;
  void TruthFilter();
  /////////////////////////////////////////////
  // Declare output data
  TTree *fTree;
  std::string fTreeName;

  TH1F* TrueEnergyDeposited;
  TH1F* TrueEnergyDepositedInROI; //Filled in TagROITruth
  TH1F* TrueChargeDeposited;
  TH1F* TrueChargeDepositedInROI; //Filled in TagROITruth
  float fECMin; //minimum energy and charge to begin accumulation
  float fHistEnergyMax;
  float fHistChargeMax;

  enum SType{
    kSAll, kSNeutrino, kSRad, kSEnd
  };


  std::vector<std::string> partTypes{"","_electron","_proton","_neutron","_photon","_other"};

  enum PType {
    kAll, kElectron, kProton, kNeutron, kPhoton, kNuc, kPEnd
  };
};





#endif