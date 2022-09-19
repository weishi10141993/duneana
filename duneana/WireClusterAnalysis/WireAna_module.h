////////////////////////////////////////////////////////////////////////
// Class:      WireAna 
// Plugin Type: analyzer (art v3_00_00)
// File:        WireAna_module.cc
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

#include "c2numpy.h"
#include "hep_hpc/hdf5/File.hpp"
#include "hep_hpc/hdf5/Column.hpp"
#include "hep_hpc/hdf5/Ntuple.hpp"
#include "hep_hpc/hdf5/PropertyList.hpp"
#include "hep_hpc/hdf5/errorHandling.hpp"
#include "hep_hpc/hdf5/make_column.hpp"
#include "hep_hpc/hdf5/make_ntuple.hpp"


// ROOT includes
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "WireAna_Utils.h"

//Others
#define DEFAULT_VALUE -99999
#define PDGDECIMAL 10000000

namespace wireana {

  //using recob::SpacePoint;
  using namespace hep_hpc::hdf5;
  using wire_nt_t = Ntuple<float>;
  using evt_nt_t = Ntuple<unsigned int, double>;
  using imager_t = Ntuple<
    int, //event id , 5,run,subrun,event, isMC,clusterID  
    float, //neutrino P, 4
    float, //int vtx, 4
    float, //part P, 4
    std::string,  //label, scalar
    std::string,  //label, scalar
    std::string,  //label, scalar
    int, //nPart, scalar
    int, //tpcid
    int, //view
    int, //pdg ,5 sorted by energy
    int, //trkid ,5 sorted by energy
    float, //charge, 5
    float, //energy, 5
    int, //n: e,gamma, p,n, meson, alpha, 6
    short, //image data (width,height,min_channel,min_time),4
    Column<double,3>, //image 
    Column<short,3>, //image pdg mask 
    Column<short,3> > ; //image pid mask 

  using apaview_xuv_t = Ntuple<
    int, //event id , 4,run,subrun,event, isMC,
    int, //apaid
    int, //neutrinotpcid
    std::string, //intType
    float, //neutrino P, 4
    float, //int vtx, 4
    float, //part P, 4
    int, //n: nu, e,gamma, p,n, meson, alpha, 7
    Column<float,2>, //image 1
//    Column<float,2>, //image 1 sig
//    Column<float,2>, //image 1 energy
//    Column<float,2>, //image 1 charge
    Column<int,2>, //image pdg mask 
//    Column<int,2>, //image pid mask 
    Column<float,2>, //image 2
//    Column<float,2>, //image 2 sig
//    Column<float,2>, //image 2 energy
//    Column<float,2>, //image 2 charge
    Column<int,2>, //image pdg mask 
//    Column<int,2>,//image pid mask 
    Column<float,2>, //image 3
//    Column<float,2>, //image 3 sig
//    Column<float,2>, //image 3 energy
//    Column<float,2>, //image 3 charge
    Column<int,2>, //image pdg mask 
//    Column<int,2>,//image pid mask 
    Column<float,2>,  //image 4
//    Column<float,2>, //image 4 sig
//    Column<float,2>, //image 4 energy
//    Column<float,2>, //image 4 charge
    Column<int,2> //image pdg mask 
//    Column<int,2>  //image pid mask 
    >;


  struct DataBlock_Truth
  {
    //Truth information
    //Incident Particle
    int truth_intType; //ES,CC,Rad: 0,1,2
    double truth_nu_momentum_x;
    double truth_nu_momentum_y;
    double truth_nu_momentum_z;
    double truth_nu_momentum_e;
    double truth_nu_vtx_x;
    double truth_nu_vtx_y;
    double truth_nu_vtx_z;
    double truth_nu_vtx_t;
    double truth_nu_PDG;
    // // Electron Info
    double truth_lep_momentum_x;
    double truth_lep_momentum_y;
    double truth_lep_momentum_z;
    double truth_lep_momentum_e;
    double truth_lep_vtx_x;
    double truth_lep_vtx_y;
    double truth_lep_vtx_z;
    double truth_lep_vtx_t;
    double truth_lep_PDG;
    // // Hadron Info
    double truth_had_momentum_x;
    double truth_had_momentum_y;
    double truth_had_momentum_z;
    double truth_had_momentum_e;
    double truth_had_vtx_x;
    double truth_had_vtx_y;
    double truth_had_vtx_z;
    double truth_had_vtx_t;
    double truth_had_PDG;

    // //ROI Info
    int nROIs;
    std::vector<int> roi_has_truth;
    std::vector<int> lead_pdg;
    std::vector<double> lead_pdg_total_energy;
    std::vector<double> lead_pdg_total_charge;
    std::vector<double> total_energy;
    std::vector<double> total_charge;
    std::vector<double> hit_fraction;
  };


  class WireAna;

}

class wireana::WireAna : public art::EDAnalyzer {
public:
  explicit WireAna(fhicl::ParameterSet const& pset);
  WireAna(WireAna const&) = delete;
  WireAna(WireAna&&) = delete;
  WireAna& operator=(WireAna const&) = delete;
  WireAna& operator=(WireAna&&) = delete;
  virtual ~WireAna() noexcept {};

  /////////////////////////////////////////////
  // Required functions.
  void analyze(art::Event const& evt) override;

  /////////////////////////////////////////////
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  /////////////////////////////////////////////
  void reset();

private:
  
  /////////////////////////////////////////////
  // Geometry Options && Tool options
  int fNPlanes;
  unsigned int fNChanPerApa, fNChanPerApaX, fNChanPerApaUV;
  unsigned int fNTicksPerWire;
  int fChannelDistance;
  int fTickDistance;
  int fMinClusterSize;

  unsigned int fMinPts;
  float fEps;
  float fDrift;
  float fPitch;
  float fDeltaMetric;

  float fMatcherDist;
  float fMatcherScore;

  unsigned int image_channel_width;
  unsigned int image_tick_width;
  unsigned int image_rebin_tick;
  unsigned int image_size;

  float fHistEnergyMax;
  float fHistChargeMax;

  /////////////////////////////////////////////
  // Backtracker services
  art::ServiceHandle<cheat::ParticleInventoryService> PIS;
 
  /////////////////////////////////////////////
  // Geometry services
  const geo::GeometryCore* geo; //= lar::providerFrom<geo::Geometry>();

  /////////////////////////////////////////////
  // config
  int fLogLevel;
  bool fDoAssns;
  bool fMakeCluster;
  bool fTagTruth;



  /////////////////////////////////////////////
  // Wire Filtering/Clustering Functions
  template<class T>
  void SortWirePtrByChannel( std::vector<art::Ptr<T>> &vec, bool increasing );
  std::vector<wireana::wirecluster> BuildInitialClusters( std::vector<art::Ptr<recob::Wire>> &vec, int dW, int dTick );

  int GetAPAID( int channel );
  int GetAPAViewID( int channel );
  int GetAPAWirePos( int channel );

  void BuildPlaneViewROIMap(  std::vector<art::Ptr<recob::Wire>> &wires );
  void BuildInitialROIClusters();
  std::vector<art::Ptr<recob::Wire>> FilterWires(std::vector<art::Ptr<recob::Wire>> &vec, int dC1, int dT1, int dCn, int dTn );
  bool HasHit( const art::Ptr<recob::Wire> &wire, int minTick );


  void PrintClusters( std::vector<wirecluster> &clusters );
  void PrintROIs( const std::vector<roi> &ROIs);

  void WriteNumPy( matchedroicluster& cluster, std::vector<art::Ptr<recob::Wire>>& wirelist );

  /////////////////////////////////////////////
  // Cluster Matching

  /////////////////////////////////////////////
  // Declare output data
  TTree *fTree;

  std::map<std::string, TH2F* > TH2FMap;

  std::vector<std::string> selTypes{"","_neutrino","_rad"};
  enum SType{
    kSAll, kSNeutrino, kSRad, kSEnd
  };
  std::vector<std::string> partTypes{"","_electron","_proton","_neutron","_photon","_other"};
  enum PType {
    kAll, kElectron, kProton, kNeutron, kPhoton, kNuc, kPEnd
  };
  //selTypes=std::vector<std::string>({"","_neutrino","_rad"});
  //partTypes=std::vector<std::string>({"","_electron","_proton","_neutron","_photon","_other"});

  void FillHistogram(std::vector<float> energies, std::vector<float>charges ,SType s, bool inROI);
  //TH2F* TrueEnergyChargeDeposited;
  //TH2F* TrueEnergyChargeDepositedInROI; //Filled in TagROITruth

  //TH2F* TrueEnergyChargeDeposited_electron;
  //TH2F* TrueEnergyChargeDepositedInROI_electron; //Filled in TagROITruth
  //TH2F* TrueEnergyChargeDeposited_proton;
  //TH2F* TrueEnergyChargeDepositedInROI_proton; //Filled in TagROITruth
  //TH2F* TrueEnergyChargeDeposited_neutron;
  //TH2F* TrueEnergyChargeDepositedInROI_neutron; //Filled in TagROITruth
  //TH2F* TrueEnergyChargeDeposited_photon;
  //TH2F* TrueEnergyChargeDepositedInROI_photon; //Filled in TagROITruth
  //TH2F* TrueEnergyChargeDeposited_other;
  //TH2F* TrueEnergyChargeDepositedInROI_other; //Filled in TagROITruth


  TH1F* TrueEnergyDeposited;
  TH1F* TrueEnergyDepositedInROI; //Filled in TagROITruth

  TH1F* TrueChargeDeposited;
  TH1F* TrueChargeDepositedInROI; //Filled in TagROITruth
  float fECMin; //minimum energy and charge to begin accumulation




  /////////////////////////////////////////////
  //Module Labels and Settting
  const art::InputTag fWireProducerLabel; 
  const art::InputTag fSimChannelLabel;
  const art::InputTag fSimulationProducerLabel;


  art::InputTag fRawProducerLabel; 
  art::InputTag fNoiselessRawProducerLabel; 

  std::string fDumpFileName;
  std::string fHDF5DumpFileName;
  std::string fHDF5DumpFileNameAPA;
  int fDumpMaxRow;
  int fDumpNClusters;

  std::string fTreeName;

  /////////////////////////////////////////////
  // Event level information
  int run;
  int subrun;
  int event;
  int MC;
  int signal_Apa;
 
  ////////////////////////////////////////////
  // Truth Operation
  // declare truth branch
  wireana::DataBlock_Truth truth_data;

  void DeclareTruthBranches(TTree*t, DataBlock_Truth &blk);
  void FillTruthBranches(art::Event const& evt, TTree*t, DataBlock_Truth &blk);
  void ResetTruthBranches(DataBlock_Truth &blk);

  // Truth tagging
  void FillTrackIDtoLabelMap( art::Event const& evt );
  void TagAllROIClusterTruth( const detinfo::DetectorClocksData &clock );
  void TagROIClusterTruth( wireana::roicluster &cluster, const detinfo::DetectorClocksData &clock );
  void TagAllROITruth( const detinfo::DetectorClocksData &clock );
  void TagROITruth( wireana::roi &roi, const detinfo::DetectorClocksData &clock );

  ////////////////////////////////////////////
  // Internal Data Structure
  std::map<raw::ChannelID_t, std::pair<art::Ptr<recob::Wire>, art::Ptr<sim::SimChannel>>> 
    ch_w_sc;
  std::map<int, std::string> 
    trkid_to_label_map;
  PlaneViewROIMap plane_view_roi_map; //type is std::map<int, std::map< geo::View_t, std::vector<roi> > >
  PlaneViewROIClusterMap plane_view_roicluster_map; //type is std::map<int, std::map< geo::View_t, std::vector<roi> > >

  ////////////////////////////////////////////
  // Image Forming:
  std::vector<double> GetArrayFromWire( std::vector<art::Ptr<recob::Wire>> &wirelist, wireana::roicluster &cluster, int channel_width, int new_tickwidth );
  std::vector<double> CombineTicks( const std::vector<double> &input, int channel_width, int nticks);
  std::vector<double> ScaleArray( const std::vector<double> &input, double min, double max );
  template <class T> 
    std::vector<T> CombineArrays( std::vector<T> &u, std::vector<T> &v, std::vector<T> &z); 
  int CalculateIndex( int c, int t, int c_width, int t_width );

  std::pair<std::vector<short>,std::vector<short>>
  GetMaskFromWire( std::vector<art::Ptr<recob::Wire>> &wirelist, wireana::roicluster &cluster, int channel_width, int new_tickwidth );

  std::pair<int,int> CalculateCT( int index, int c_width );

  std::map<int, std::string> fViewMap;

  c2numpy_writer npywriter;

  /////////////////////////////////////////////
  // HDF5 Create
  void CreateHDF5DataSet( art::Event const & evt, wireana::matchedroicluster &mcluster,  std::vector<art::Ptr<recob::Wire>>& wirelist, std::vector<art::Ptr<sim::SimChannel>>& chlist) ;


  std::map<int, std::vector<std::vector<float>>>
    GetRawMap( std::vector<art::Ptr<raw::RawDigit>>& rawList );

  void GetMaskMap( std::vector<art::Ptr<sim::SimChannel>>& chlist );

  void CreateHDF5ApaView( 
      art::Event const & evt, 
      std::vector<art::Ptr<raw::RawDigit>>& rawList,
      std::vector<art::Ptr<raw::RawDigit>>& noiselessRawList,
      std::vector<art::Ptr<sim::SimChannel>>& chlist);
  std::map<int, std::vector<std::vector<float>>> plane_raw_map;
  std::map<int, std::vector<std::vector<float>>> plane_rawsig_map;
  std::map<int, std::vector<std::vector<float>>> plane_energy_map;
  std::map<int, std::vector<std::vector<float>>> plane_charge_map;
  std::map<int, std::vector<std::vector<int>>> plane_pdg_map;
  std::map<int, std::vector<std::vector<int>>> plane_trkid_map;

  File hdffile;
  File hdffileAPA;
  imager_t *image_tuple;
  apaview_xuv_t *apa_tuple_xuv;
  //apaview_t *apa_tuple_x;
  //apaview_t *apa_tuple_u;
  //apaview_t *apa_tuple_v;
  //imager_t image_tuple;

  //unsigned int -- event id,4
  //string --> generator label

  int fApaNColums;
  int fApaNRows;

  int fBuffSize;

  /////////////////////////////////////////////
  // PDG parser codes
  int PDGReducer( int pdg ) { return pdg%PDGDECIMAL; }
  int PDGEncoder( int pdg, int gencode ){ return pdg%PDGDECIMAL + gencode*PDGDECIMAL; }

  std::map<std::string, int> TruthCode;

};





#endif
