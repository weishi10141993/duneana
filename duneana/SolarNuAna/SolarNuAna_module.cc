////////////////////////////////////////////////////////////////////////
// Class:       SolarNuAna
// Module Type: analyzer
//
// Marcos Vínius, Marco Dias, built on DAQSimAna_module.cc
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// 
//  Erros: 1057, 1145, 1167, 1351, 1503, 1596, 1774, 1969 ;  NOVO: 1045  
//             
//  /cvmfs/larsoft.opensciencegrid.org/products/larsim/v09_22_06/include/larsim/MCCheater/BackTrackerService.h
// 
////////////////////////////////////////////////////////////////////////



// C++ includes

// ROOT includes
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TGeoMatrix.h"

// Framework includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"         
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"                        // ...
#include "larsim/MCCheater/PhotonBackTrackerService.h"
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

#include "TMinuit.h"

const int nMaxHits = 8000;
const int nMaxOpFlash = 1000;
//const int nMaxDigs = 4492; // unused

enum PType{ kUnknown, kMarl, kAPA, kCPA, kAr39, kAr42, kNeut, kKryp, kPlon, kRdon };

class SolarNuAna : public art::EDAnalyzer {

public:

  explicit SolarNuAna(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  SolarNuAna(SolarNuAna const &) = delete;
  SolarNuAna(SolarNuAna &&) = delete;
  SolarNuAna & operator = (SolarNuAna const &) = delete;
  SolarNuAna & operator = (SolarNuAna &&) = delete;

  // The main guts...
  void analyze(art::Event const & evt) override;

  void reconfigure(fhicl::ParameterSet const & p);

  void beginJob() override;

private:

  // --- Some of our own functions.
  int ViewMatchCluster(std::vector< std::vector<recob::Hit> > testc,
                       double t, double z, double a0, int xsign,
                       double &curY, double & curZ, int minN=3);
  void ResetVariables();
  void  FillMyMaps  ( std::map< int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle< std::vector<simb::MCTruth> > Hand );
  PType WhichParType( int TrID );
  bool  InMyMap     ( int TrID, std::map< int, simb::MCParticle> ParMap );
  void  CalcAdjHits ( std::vector< recob::Hit > MyVec,
                      std::vector< std::vector<recob::Hit> >& clusters,
                      TH1I* MyHist, TH1F* MyADCIntHist, bool HeavDebug );


  // --- Our fcl parameter labels for the modules that made the data products
  std::string fRawDigitLabel;
  std::string fHitLabel;
  std::string fOpHitLabel;
  std::string fOpFlashLabel;

  std::string fGEANTLabel;
  std::string fMARLLabel; std::map< int, simb::MCParticle > MarlParts;
  std::string fAPALabel;  std::map< int, simb::MCParticle > APAParts;
  std::string fCPALabel;  std::map< int, simb::MCParticle > CPAParts;
  std::string fAr39Label; std::map< int, simb::MCParticle > Ar39Parts;
  std::string fAr42Label; std::map< int, simb::MCParticle > Ar42Parts;
  std::string fNeutLabel; std::map< int, simb::MCParticle > NeutParts;
  std::string fKrypLabel; std::map< int, simb::MCParticle > KrypParts;
  std::string fPlonLabel; std::map< int, simb::MCParticle > PlonParts;
  std::string fRdonLabel; std::map< int, simb::MCParticle > RdonParts;

  // --- Other variables
  //int nADC; // no longer used

  // --- Our TTree, and its associated variables.
  TTree* fDAQSimTree;
  TTree* fHitTree;
  TTree* fOpHitTree;
  TTree* fTrkHitTree;
  TTree* fMarlHitTree;
  TTree* fNearHitTree;
  TTree* fDAQClusterTree;

  TTree* fNeutronTruthTree;
  TTree* fAlphaTruthTree;
  // General event info.
  int Run;
  int SubRun;
  int Event;
  // Raw digits
  //int NTotDigs; // unused

  // Truth information about our generated neutrino event
  float NuEnergy;
  float EEnergy;
  float NeutPhotEnergy;
  float NeutEnergy;
  float HiEEnergy;
  float PartX;
  float PartY;
  float PartZ;
  float MomX;
  float MomY;
  float MomZ;
  float TrueX;
  float TrueY;
  float TrueZ;

  float BiggestT0;
  float BiggestT1;
  float BiggestT2;
  float delayedPE;

  // OpFlash variables
  int NOpFlash;
  float OpFlashMarlPur[nMaxOpFlash];
  float OpFlashPE[nMaxOpFlash];
  float OpFlashT[nMaxOpFlash];
  float OpFlashDeltaT[nMaxOpFlash];
  float OpFlashNHit[nMaxOpFlash];
  float OpFlashY[nMaxOpFlash];
  float OpFlashdY[nMaxOpFlash];
  float OpFlashZ[nMaxOpFlash];
  float OpFlashdZ[nMaxOpFlash];
  float OpFlashFrame[nMaxOpFlash];

  // The reconstructed hits
  int   NTotHits;
  int   NColHits;
  int   NIndHits;

  int   HitView[nMaxHits]; ///< View i.e Coll, U, V
  int   HitSize[nMaxHits]; ///< Time width (ticks) Start - End time
  int   HitTPC [nMaxHits]; ///< The TPC which the hit occurs in
  int   HitChan[nMaxHits]; ///< The channel which the hit occurs on
  float HitTime[nMaxHits]; ///< The time of the hit (ticks)
  float HitRMS [nMaxHits]; ///< The RMS of the hit
  float HitSADC[nMaxHits]; ///< The summed ADC of the hit
  float HitInt [nMaxHits]; ///< The ADC integral of the hit
  float HitPeak[nMaxHits]; ///< The peak ADC value of the hit
  float HitFrame[nMaxHits];///< The frame hit was recorded in
  int   GenType[nMaxHits]; ///< The generator which generated the particle responsible for the hit
  int   TotGen_Marl;
  int   TotGen_APA;
  int   TotGen_CPA;
  int   TotGen_Ar39;
  int   TotGen_Neut;
  int   TotGen_Kryp;
  int   TotGen_Plon;
  int   TotGen_Rdon;

  int trkIdx;
  float trkLen;

  float hitPur, hitT, hitZ, hitCharge;
  float fracE, fracGa, fracPr;

  int clustIdx;
  float NhitPur, NhitT, NhitZ, NhitCharge;

  // DAQ cluster variables
  int view;
  float dzdy;
  float charge;
  float nhit;
  float nhit50;
  float clustY;
  float clustZ;
  float frame;
  float clustT;
  float maxHit;
  int nN;
  float ntotCharge;

  int nNeut;
  float neutCaptureTime;

  std::vector<float> nNHit, ndZ, ndT, ndY, nR, nCharge, neutE, rdonE;
  std::vector<int> nPur;

  float m1T;
  float m1Z;
  float m1Y;
  float m1dzdy;
  float m1Charge;
  float m1NHit;
  float m150NHit;
  float m1Pur;
  float m2T;
  float m2Z;
  float m2Y;
  float m2dzdy;
  float m2Charge;
  float m2NHit;
  float m250NHit;
  float m2Pur;
  float fitY;
  float OpT;
  float OpdT;
  float OpY;
  float OpZ;
  float OpPE;
  float OpHiFrac;

  int OpHitIdx;
  float OpHitT;
  float OpHitY;
  float OpHitZ;
  float OpHitPE;

  float OpPur;
  float marlPur;
  float apaPur;
  float cpaPur;
  float Ar39Pur;
  float Ar42Pur;
  float neutPur;
  float krypPur;
  float plonPur;
  float rdonPur;

  float marlEDep;

  float edirx;
  float ediry;
  float edirz;
  float nudirx;
  float nudiry;
  float nudirz;
  float trkdirx;
  float trkdiry;
  float trkdirz;
  float trkl;
  float trknhit;
  float trkcos;
  float trkcharge;
  float enuCos;
  float etrkCos;
  float enuFlippedCos;
  float etrkFlippedCos;



  // histograms to fill about Collection plane hits
  TH1I* hAdjHits;
  TH1I* hAdjHits_Marl;
  TH1I* hAdjHits_APA;
  TH1I* hAdjHits_CPA;
  TH1I* hAdjHits_Ar39;
  TH1I* hAdjHits_Neut;
  TH1I* hAdjHits_Kryp;
  TH1I* hAdjHits_Plon;
  TH1I* hAdjHits_Rdon;
  TH1I* hAdjHits_Oth;

  // histograms saving the total ADC integral for each group of adjacent hits
  TH1F* hAdjHitsADCInt;
  TH1F* hAdjHitsADCInt_Marl;
  TH1F* hAdjHitsADCInt_APA;
  TH1F* hAdjHitsADCInt_CPA;
  TH1F* hAdjHitsADCInt_Ar39;
  TH1F* hAdjHitsADCInt_Neut;
  TH1F* hAdjHitsADCInt_Kryp;
  TH1F* hAdjHitsADCInt_Plon;
  TH1F* hAdjHitsADCInt_Rdon;
  TH1F* hAdjHitsADCInt_Oth;

  // --- Declare our services
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::PhotonBackTrackerService> pbt;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

};

//......................................................
SolarNuAna::SolarNuAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  this->reconfigure(p);
}

//......................................................
void SolarNuAna::reconfigure(fhicl::ParameterSet const & p)
{
  fRawDigitLabel = p.get<std::string> ("RawDigitLabel");
  fHitLabel      = p.get<std::string> ("HitLabel");
  fOpFlashLabel  = p.get<std::string>("OpFlashLabel");
  fOpHitLabel    = p.get<std::string>("OpHitLabel");

  fGEANTLabel = p.get<std::string> ("GEANT4Label");
  fMARLLabel = p.get<std::string> ("MARLEYLabel");
  fAPALabel  = p.get<std::string> ("APALabel");
  fCPALabel  = p.get<std::string> ("CPALabel");
  fAr39Label = p.get<std::string> ("Argon39Label");
  fAr42Label = p.get<std::string> ("Argon42Label");
  fNeutLabel = p.get<std::string> ("NeutronLabel");
  fKrypLabel = p.get<std::string> ("KryptonLabel");
  fPlonLabel = p.get<std::string> ("PoloniumLabel");
  fRdonLabel = p.get<std::string> ("RadonLabel");

} // Reconfigure

//......................................................
void SolarNuAna::ResetVariables()
{
  // Clear my MCParticle maps.
  MarlParts.clear(); APAParts .clear(); CPAParts .clear(); Ar39Parts.clear();
  NeutParts.clear(); KrypParts.clear(); PlonParts.clear(); RdonParts.clear();
  Ar42Parts.clear();

  // General event info.
  Run = SubRun = Event = -1;

  // Set Number of GenParts to 0
  TotGen_Marl = TotGen_APA  = TotGen_CPA  = TotGen_Ar39 = 0;
  TotGen_Neut = TotGen_Kryp = TotGen_Plon = TotGen_Rdon = 0;

  // Summary of OpFlash's in event
  NOpFlash = 0;
  /*
  for (int hh=0; hh<nMaxHits; ++hh) {
    OpFlashPE[hh] = OpFlashT[hh] = OpFlashMarlPur[hh] = 0;
    OpFlashDeltaT[hh] = OpFlashNHit[hh] = 0;
    OpFlashY[hh] = OpFlashdY[hh] = 0;
    OpFlashZ[hh] = OpFlashdZ[hh] = 0;
    OpFlashFrame[hh] = 0;
  }
  */

  // Reconstructed hits
  NTotHits = NColHits = NIndHits = 0;

  /*
  for (int hh=0; hh<nMaxHits; ++hh) {
    HitView[hh] = HitSize[hh] = HitChan[hh] = GenType[hh] = 0;
    HitTime[hh] = HitRMS [hh] = HitSADC[hh] = 0;
    HitInt [hh] = HitPeak[hh] = HitFrame[hh] = HitTPC [hh] = 0;
  }
  */
} // ResetVariables

//......................................................
void SolarNuAna::beginJob()
{

  // --- Make our handle to the TFileService
  art::ServiceHandle<art::TFileService> tfs;

  fTrkHitTree = tfs->make<TTree>("TrkHitTree","DAQ simulation analysis tree");
  fTrkHitTree->Branch("Event", &Event);
  fTrkHitTree->Branch("SubRun", &SubRun);
  fTrkHitTree->Branch("Run", &Run);
  fTrkHitTree->Branch("trkIdx", &trkIdx);
  fTrkHitTree->Branch("trkLen", &trkLen);
  fTrkHitTree->Branch("hitZ", &hitZ);
  fTrkHitTree->Branch("hitT", &hitT);
  fTrkHitTree->Branch("hitCharge", &hitCharge);

  fMarlHitTree = tfs->make<TTree>("MarlHitTree","");
  fMarlHitTree->Branch("Event", &Event);
  fMarlHitTree->Branch("SubRun", &SubRun);
  fMarlHitTree->Branch("Run", &Run);
  fMarlHitTree->Branch("TrueX", &TrueX);
  fMarlHitTree->Branch("TrueY", &TrueY);
  fMarlHitTree->Branch("TrueZ", &TrueZ);
  fMarlHitTree->Branch("hitZ", &hitZ);
  fMarlHitTree->Branch("hitT", &hitT);
  fMarlHitTree->Branch("hitCharge", &hitCharge);

  fNearHitTree = tfs->make<TTree>("NearHitTree","");
  fNearHitTree->Branch("Event", &Event);
  fNearHitTree->Branch("SubRun", &SubRun);
  fNearHitTree->Branch("Run", &Run);
  fNearHitTree->Branch("hitZ", &NhitZ);
  fNearHitTree->Branch("hitT", &NhitT);
  fNearHitTree->Branch("hitCharge", &NhitCharge);
  fNearHitTree->Branch("hitPur", &NhitPur);
  fNearHitTree->Branch("clustIdx", &clustIdx);

  fOpHitTree = tfs->make<TTree>("OpHitTree","DAQ simulation analysis tree");
  fOpHitTree->Branch("Event", &Event);
  fOpHitTree->Branch("SubRun", &SubRun);
  fOpHitTree->Branch("Run", &Run);
  fOpHitTree->Branch("TrueX", &TrueX);
  fOpHitTree->Branch("TrueY", &TrueY);
  fOpHitTree->Branch("TrueZ", &TrueZ);
  fOpHitTree->Branch("OpHitIdx", &OpHitIdx);
  fOpHitTree->Branch("OpHitT", &OpHitT);
  fOpHitTree->Branch("OpHitPE", &OpHitPE);
  fOpHitTree->Branch("OpHitY", &OpHitY);
  fOpHitTree->Branch("OpHitZ", &OpHitZ);


  fHitTree = tfs->make<TTree>("HitTree","DAQ simulation analysis tree");
  fHitTree->Branch("Event", &Event);
  fHitTree->Branch("SubRun", &SubRun);
  fHitTree->Branch("Run", &Run);
  fHitTree->Branch("hitPur", &hitPur);
  fHitTree->Branch("hitT", &hitT);
  fHitTree->Branch("hitZ", &hitZ);
  fHitTree->Branch("hitCharge", &hitCharge);
  fHitTree->Branch("fracE", &fracE);
  fHitTree->Branch("fracGa",&fracGa);
  fHitTree->Branch("fracPr",&fracPr);

  // --- Our TTree
  fDAQSimTree = tfs->make<TTree>("DAQSimTree","DAQ simulation analysis tree");
  // General event information...
  fDAQSimTree -> Branch( "Run"   , &Run   , "Run/I"    );
  fDAQSimTree -> Branch( "SubRun", &SubRun, "SubRun/I" );
  fDAQSimTree -> Branch( "Event" , &Event , "Event/I"  );
  // Neutrino truth information
  fDAQSimTree -> Branch( "NuEnergy"  , &NuEnergy  , "NuEnergy/F"  );
  fDAQSimTree -> Branch( "EEnergy"  , &EEnergy  , "EEnergy/F"  );
  fDAQSimTree -> Branch( "HiEEnergy"     , &HiEEnergy     , "HiEEnergy/F" );
  fDAQSimTree -> Branch( "enuCos"     , &enuCos     , "enuCos/F"     );
  fDAQSimTree -> Branch( "TrueX"     , &TrueX     , "TrueX/F"     );
  fDAQSimTree -> Branch( "TrueY"     , &TrueY     , "TrueY/F"     );
  fDAQSimTree -> Branch( "TrueZ"     , &TrueZ     , "TrueZ/F"     );

  fDAQSimTree -> Branch( "BiggestT0" , &BiggestT0 , "BiggestT0/F" );
  fDAQSimTree -> Branch( "BiggestT1" , &BiggestT1 , "BiggestT1/F" );
  fDAQSimTree -> Branch( "BiggestT2" , &BiggestT2 , "BiggestT2/F" );


  fDAQClusterTree = tfs->make<TTree>("DAQClusterTree","DAQ Cluster Tree");

  fDAQClusterTree -> Branch( "Run"   , &Run   , "Run/I"    );
  fDAQClusterTree -> Branch( "SubRun", &SubRun, "SubRun/I" );
  fDAQClusterTree -> Branch( "Event" , &Event , "Event/I"  );
  fDAQClusterTree->Branch("view",   &view,   "view/I");
  fDAQClusterTree->Branch("dzdy",   &dzdy,   "dzdy/F");
  fDAQClusterTree->Branch("charge", &charge, "charge/F");
  fDAQClusterTree->Branch("nhit",   &nhit,   "nhit/F");
  fDAQClusterTree->Branch("nhit50",   &nhit50,   "nhit50/F");
  fDAQClusterTree->Branch("clustY", &clustY, "clustY/F");
  fDAQClusterTree->Branch("clustZ", &clustZ, "clustZ/F");
  fDAQClusterTree->Branch("frame",  &frame,  "frame/F");
  fDAQClusterTree->Branch("clustT",  &clustT,  "clustT/F");
  fDAQClusterTree->Branch("maxHit",   &maxHit,   "maxHit/F");

  fDAQClusterTree->Branch("neutE",&neutE);
  fDAQClusterTree->Branch("rdonE",&rdonE);

  fDAQClusterTree->Branch("nN",&nN,"nN/I");
  fDAQClusterTree->Branch("nPur",&nPur);
  fDAQClusterTree->Branch("nNHit",&nNHit);
  fDAQClusterTree->Branch("nCharge",&nCharge);
  fDAQClusterTree->Branch("ndZ",&ndZ);
  fDAQClusterTree->Branch("ndT",&ndT);
  fDAQClusterTree->Branch("ndY",&ndY);
  fDAQClusterTree->Branch("nR",&nR);
  fDAQClusterTree->Branch("ntotCharge",&ntotCharge);


  fDAQClusterTree->Branch("m1T",     &m1T,     "m1T/F");
  fDAQClusterTree->Branch("m1Y",     &m1Y,     "m1Y/F");
  fDAQClusterTree->Branch("m1Z",     &m1Z,     "m1Z/F");
  fDAQClusterTree->Branch("m1dzdy",  &m1dzdy,  "m1dzdy/F");
  fDAQClusterTree->Branch("m1Charge",&m1Charge,"m1Charge/F");
  fDAQClusterTree->Branch("m1NHit",&m1NHit,"m1NHit/F");
  fDAQClusterTree->Branch("m150NHit",&m150NHit,"m150NHit/F");
  fDAQClusterTree->Branch("m1Pur",   &m1Pur,   "m1Pur/F");
  fDAQClusterTree->Branch("m2T",     &m2T,     "m2T/F");
  fDAQClusterTree->Branch("m2Y",     &m2Y,     "m2Y/F");
  fDAQClusterTree->Branch("m2Z",     &m2Z,     "m2Z/F");
  fDAQClusterTree->Branch("m2dzdy",  &m2dzdy,  "m2dzdy/F");
  fDAQClusterTree->Branch("m2Charge",&m2Charge,"m2Charge/F");
  fDAQClusterTree->Branch("m2NHit",&m2NHit,"m2NHit/F");
  fDAQClusterTree->Branch("m250NHit",&m250NHit,"m250NHit/F");
  fDAQClusterTree->Branch("m2Pur",   &m2Pur,   "m2Pur/F");
  fDAQClusterTree->Branch("fitY",    &fitY,    "fitY/F");

  fDAQClusterTree->Branch("OpT",    &OpT,    "OpT/F");
  fDAQClusterTree->Branch("OpdT",   &OpdT,   "OpdT/F");
  fDAQClusterTree->Branch("OpY",    &OpY,    "OpY/F");
  fDAQClusterTree->Branch("OpZ",    &OpZ,    "OpZ/F");
  fDAQClusterTree->Branch("OpPE",   &OpPE,   "OpPE/F");
  fDAQClusterTree->Branch("OpHiFrac",   &OpHiFrac,   "OpHiFrac/F");
  fDAQClusterTree->Branch("OpPur",  &OpPur,  "OpPur/F");
  fDAQClusterTree->Branch("delayedPE",   &delayedPE,   "delayedPE/F");

  fDAQClusterTree->Branch("TrueX",   &TrueX,   "TrueX/F");
  fDAQClusterTree->Branch("TrueY",   &TrueY,   "TrueY/F");
  fDAQClusterTree->Branch("TrueZ",   &TrueZ,   "TrueZ/F");
  fDAQClusterTree->Branch("PartX",   &PartX,   "PartX/F");
  fDAQClusterTree->Branch("PartY",   &PartY,   "PartY/F");
  fDAQClusterTree->Branch("PartZ",   &PartZ,   "PartZ/F");
  fDAQClusterTree->Branch("MomX",   &MomX,   "MomX/F");
  fDAQClusterTree->Branch("MomY",   &MomY,   "MomY/F");
  fDAQClusterTree->Branch("MomZ",   &MomZ,   "MomZ/F");
  fDAQClusterTree->Branch("NuEnergy",&NuEnergy,"NuEnergy/F");
  fDAQClusterTree->Branch("EEnergy",&EEnergy,"EEnergy/F");
  fDAQClusterTree->Branch("NeutEnergy",&NeutEnergy,"NeutEnergy/F");
  fDAQClusterTree->Branch("NeutPhotEnergy",&NeutPhotEnergy,"NeutPhotEnergy/F");
  fDAQClusterTree->Branch("HiEEnergy",&HiEEnergy,"HiEEnergy/F");

  fDAQClusterTree->Branch("nudirx", &nudirx, "nudirx/F");
  fDAQClusterTree->Branch("nudiry", &nudiry, "nudiry/F");
  fDAQClusterTree->Branch("nudirz", &nudirz, "nudirz/F");
  fDAQClusterTree->Branch("edirx", &edirx, "edirx/F");
  fDAQClusterTree->Branch("ediry", &ediry, "ediry/F");
  fDAQClusterTree->Branch("edirz", &edirz, "edirz/F");

  fDAQClusterTree->Branch("marlPur", &marlPur, "marlpur/F");
  fDAQClusterTree->Branch("apaPur",  &apaPur,  "apapur/F");
  fDAQClusterTree->Branch("cpaPur",  &cpaPur,  "cpapur/F");
  fDAQClusterTree->Branch("Ar39Pur", &Ar39Pur, "Ar39pur/F");
  fDAQClusterTree->Branch("Ar42Pur", &Ar42Pur, "Ar42pur/F");
  fDAQClusterTree->Branch("neutPur", &neutPur, "neutpur/F");
  fDAQClusterTree->Branch("krypPur", &krypPur, "kryppur/F");
  fDAQClusterTree->Branch("plonPur", &plonPur, "plonpur/F");
  fDAQClusterTree->Branch("rdonPur", &rdonPur, "rdonpur/F");

  fDAQClusterTree->Branch("trkdirx", &trkdirx, "trkdirx/F");
  fDAQClusterTree->Branch("trkdiry", &trkdiry, "trkdiry/F");
  fDAQClusterTree->Branch("trkdirz", &trkdirz, "trkdirz/F");
  fDAQClusterTree->Branch("trkl", &trkl, "trkl/F");
  fDAQClusterTree->Branch("trknhit", &trknhit, "trknhit/F");
  fDAQClusterTree->Branch("trkcharge", &trkcharge, "trkcharge/F");
  fDAQClusterTree->Branch("trkcos", &trkcos, "trkcos/F");

  fDAQClusterTree->Branch("nNeut", &nNeut, "nNeut/I");
  fDAQClusterTree->Branch("neutCaptureTime", &neutCaptureTime, "neutCaptureTime/F");
  fDAQClusterTree->Branch("marlEDep", &marlEDep, "marlEDep/F");

  fDAQClusterTree->Branch("enuCos", &enuCos, "enuCos/F");
  fDAQClusterTree->Branch("etrkCos", &etrkCos, "etrkCos/F");
  fDAQClusterTree->Branch("etrkFlippedCos", &etrkFlippedCos, "etrkFlippedCos/F");




  // --- Our Histograms...
  hAdjHits = tfs->make<TH1I>("hAdjHits", "Number of adjacent collection plane hits; Number of adjacent collection plane hits; Number of events"  , 21, -0.5, 20.5 );
  hAdjHits_Marl = tfs->make<TH1I>("hAdjHits_Marl", "Number of adjacent collection plane hits for MARLEY; Number of adjacent collection plane hits; Number of events"  , 21, -0.5, 20.5 );
  hAdjHits_APA  = tfs->make<TH1I>("hAdjHits_APA" , "Number of adjacent collection plane hits for APAs; Number of adjacent collection plane hits; Number of events"    , 21, -0.5, 20.5 );
  hAdjHits_CPA  = tfs->make<TH1I>("hAdjHits_CPA" , "Number of adjacent collection plane hits for CPAs; Number of adjacent collection plane hits; Number of events"    , 21, -0.5, 20.5 );
  hAdjHits_Ar39 = tfs->make<TH1I>("hAdjHits_Ar39", "Number of adjacent collection plane hits for Argon39; Number of adjacent collection plane hits; Number of events" , 21, -0.5, 20.5 );
  hAdjHits_Neut = tfs->make<TH1I>("hAdjHits_Neut", "Number of adjacent collection plane hits for Neutrons; Number of adjacent collection plane hits; Number of events", 21, -0.5, 20.5 );
  hAdjHits_Kryp = tfs->make<TH1I>("hAdjHits_Kryp", "Number of adjacent collection plane hits for Krypton; Number of adjacent collection plane hits; Number of events" , 21, -0.5, 20.5 );
  hAdjHits_Plon = tfs->make<TH1I>("hAdjHits_Plon", "Number of adjacent collection plane hits for Polonium; Number of adjacent collection plane hits; Number of events", 21, -0.5, 20.5 );
  hAdjHits_Rdon = tfs->make<TH1I>("hAdjHits_Rdon", "Number of adjacent collection plane hits for Radon; Number of adjacent collection plane hits; Number of events"   , 21, -0.5, 20.5 );
  hAdjHits_Oth  = tfs->make<TH1I>("hAdjHits_Oth" , "Number of adjacent collection plane hits for Others; Number of adjacent collection plane hits; Number of events"  , 21, -0.5, 20.5 );


  hAdjHitsADCInt = tfs->make<TH1F>("hAdjHitsADCInt", "Total summed ADC Integrals for clusters; Total summed ADC Integrals for clusters; Number of events"  , 1000, 0, 10000 );
  hAdjHitsADCInt_Marl = tfs->make<TH1F>("hAdjHitsADCInt_Marl", "Total summed ADC Integrals for clusters for MARLEY; Total summed ADC Integrals for clusters; Number of events"  , 1000, 0, 10000 );
  hAdjHitsADCInt_APA  = tfs->make<TH1F>("hAdjHitsADCInt_APA" , "Total summed ADC Integrals for clusters for APAs; Total summed ADC Integrals for clusters; Number of events"    , 1000, 0, 10000 );
  hAdjHitsADCInt_CPA  = tfs->make<TH1F>("hAdjHitsADCInt_CPA" , "Total summed ADC Integrals for clusters for CPAs; Total summed ADC Integrals for clusters; Number of events"    , 1000, 0, 10000 );
  hAdjHitsADCInt_Ar39 = tfs->make<TH1F>("hAdjHitsADCInt_Ar39", "Total summed ADC Integrals for clusters for Argon39; Total summed ADC Integrals for clusters; Number of events" , 1000, 0, 10000 );
  hAdjHitsADCInt_Neut = tfs->make<TH1F>("hAdjHitsADCInt_Neut", "Total summed ADC Integrals for clusters for Neutrons; Total summed ADC Integrals for clusters; Number of events", 1000, 0, 10000 );
  hAdjHitsADCInt_Kryp = tfs->make<TH1F>("hAdjHitsADCInt_Kryp", "Total summed ADC Integrals for clusters for Krypton; Total summed ADC Integrals for clusters; Number of events" , 1000, 0, 10000 );
  hAdjHitsADCInt_Plon = tfs->make<TH1F>("hAdjHitsADCInt_Plon", "Total summed ADC Integrals for clusters for Polonium; Total summed ADC Integrals for clusters; Number of events", 1000, 0, 10000 );
  hAdjHitsADCInt_Rdon = tfs->make<TH1F>("hAdjHitsADCInt_Rdon", "Total summed ADC Integrals for clusters for Radon; Total summed ADC Integrals for clusters; Number of events"   , 1000, 0, 10000 );
  hAdjHitsADCInt_Oth  = tfs->make<TH1F>("hAdjHitsADCInt_Oth" , "Total summed ADC Integrals for clusters for Others; Total summed ADC Integrals for clusters; Number of events"  , 1000, 0, 10000 );

} // BeginJob

//......................................................
void SolarNuAna::analyze(art::Event const & evt)
{

  // --- We want to reset all of our TTree variables...
  ResetVariables();

  // --- Set all of my general event information...
  Run    = evt.run();
  SubRun = evt.subRun();
  Event  = evt.event();

  NuEnergy = -1;
  EEnergy = -1;
  TrueX = -1e3;
  TrueY = -1e3;
  TrueZ = -1e3;

  // --- Lift out the MARLEY particles.

  // Needed for marley running
  auto MarlTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fMARLLabel);
  if (MarlTrue->size() != 0){
    if (MarlTrue->at(0).NParticles() > 1){
      const simb::MCParticle& part = MarlTrue->at(0).GetParticle(0);

      NuEnergy = part.E();
      const simb::MCParticle& lep = MarlTrue->at(0).GetParticle(1);
      TrueX = lep.Vx();
      TrueY = lep.Vy();
      TrueZ = lep.Vz();
      EEnergy = lep.E();

      std::cout << "Parts and leps " << part.P() <<"   " << lep.P() << std::endl;
      std::cout << "Parts and leps " << MarlTrue->at(0).NParticles() << std::endl;

      nudirx = part.Px()/part.P();
      nudiry = part.Py()/part.P();
      nudirz = part.Pz()/part.P();

      std::cout << nudirx << "  " << nudiry << "  "<< nudirz << std::endl;
      //std::cout << edirx << "  " << ediry << "  "<< edirz << std::endl;

    }
  }

  // Find OpFlashes associated with the event
  art::Handle< std::vector< recob::OpFlash > > FlashHandle;
  std::vector<art::Ptr<recob::OpFlash> > flashlist;
  if (evt.getByLabel(fOpFlashLabel, FlashHandle)) {
    art::fill_ptr_vector(flashlist, FlashHandle);
    //std::sort(flashlist.begin(), flashlist.end(), recob::OpFlashPtrSortByPE);
  }
  // Grab assns with OpHits to get match to neutrino purity
  art::FindManyP< recob::OpHit > OpAssns(flashlist, evt, fOpFlashLabel);
  std::cout << "total number of flashes we saw: " << flashlist.size() << std::endl;
  NOpFlash = flashlist.size();

  std::vector<double> fracHiPE;

  for ( int i = 0; i < int(flashlist.size()); i++ ){
    recob::OpFlash TheFlash = *flashlist[i];
    std::vector< art::Ptr< recob::OpHit > > matchedHits = OpAssns.at(i);

    double totPE = 0;
    double hiPE = 0;
    for (int j = 0; j < int(matchedHits.size()); j++){
      recob::OpHit ohit = *matchedHits[j];
      totPE += ohit.PE();
      OpHitIdx = i;
      OpHitPE = ohit.PE();
      OpHitT = ohit.PeakTimeAbs();

      if (OpHitPE > hiPE) hiPE = OpHitPE;

      int chan = ohit.OpChannel();
      double xyz[3];
      geo->OpDetGeoFromOpChannel(chan).GetCenter(xyz);

      OpHitY = xyz[1];
      OpHitZ = xyz[2];

      fOpHitTree->Fill();
      //if (abs(OpHitT)<10) fOpHitTree->Fill();

    }

    fracHiPE.push_back(hiPE/totPE);

  }


  // --- Lift out the TPC raw digits:
  //auto rawdigits = evt.getValidHandle<std::vector<raw::RawDigit> >(fRawDigitLabel);

  // Lift out reconstructed tracks, and pull associations for hits
  art::Handle< std::vector< recob::Track > > trkHandle;
  std::vector<art::Ptr<recob::Track> > trklist;
  if (evt.getByLabel("pmtracktc", trkHandle)) {
    art::fill_ptr_vector(trklist, trkHandle);
    //std::sort(flashlist.begin(), flashlist.end(), recob::OpFlashPtrSortByPE);
  }
  int nHitBiggestTrk = -1;
  int bigTrkIdx = -1;
  for (int i = 0; i < int(trklist.size()); i++){
    trkIdx = i;
    art::Ptr<recob::Track> tr(trkHandle,i);
    if (int(tr->NPoints())>nHitBiggestTrk){
      bigTrkIdx = i;
      nHitBiggestTrk = tr->NPoints();
    }
  }
  trkdirx = -5;
  trkdiry = -5;
  trkdirz = -5;
  if (bigTrkIdx>=0){
    art::Ptr<recob::Track> tr(trkHandle,bigTrkIdx);
    std::vector< art::Ptr<recob::Track> > trks = {tr};
    recob::Track::Vector_t trkdir = tr->StartDirection();

    //std::vector<art::Ptr<recob::Hit> > trkhitlist;
    art::FindManyP< recob::Hit > HitAssns(trks, evt, "pmtracktc");
    for ( int i = 0; i < int(trks.size()); i++ ){
      trkIdx = i;

      recob::Track TheFlash = *trks[i];
      std::vector< art::Ptr< recob::Hit > > maedHits = HitAssns.at(i);

      trkLen = TheFlash.Length();

      for (int j = 0; j < int(maedHits.size()); j++){
        recob::Hit chit = *maedHits[j];
        if (chit.View() != 2) continue;
        std::cout << "Track Hit " << chit.Integral() << "  " << chit.Channel() << std::endl;

        hitT = 0.5*(chit.EndTick() + chit.StartTick());
        hitCharge = chit.Integral();
        const geo::WireGeo* wire = geo->GeometryCore::WirePtr(chit.WireID());
        double hXYZ[3];
        wire->GetCenter(hXYZ);
        hitZ = hXYZ[2];

        fTrkHitTree->Fill();
      }

      trknhit = maedHits.size();

    }


    std::cout << "trk list " << trknhit << std::endl;
    std::cout << "trk list " << trknhit << std::endl;
    std::cout << "trk list " << trknhit << std::endl;
    std::cout << "trk list " << trknhit << std::endl;
    std::cout << "trk list " << trknhit << std::endl;
    std::cout << "trk list " << trknhit << std::endl;
    std::cout << "trk list " << trknhit << std::endl;

    std::cout << "Found biggest trk with " << tr->NPoints() << std::endl;
    trkdirx = trkdir.X();
    trkdiry = trkdir.Y();
    trkdirz = trkdir.Z();

    trkl = tr->Length();
  }
  // Grab assns with OpHits to get match to neutrino purity
  //art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trkHandle, evt, "pmtracktc");




  art::FindManyP<simb::MCParticle> MarlAssn(MarlTrue,evt,fGEANTLabel);
  FillMyMaps( MarlParts, MarlAssn, MarlTrue );
  TotGen_Marl = MarlParts.size();
  std::cout << "--- The size of MarleyParts is " << MarlParts.size() << std::endl;

  std::set< int > marl_trackids;
  for ( std::map<int,simb::MCParticle>::iterator iter = MarlParts.begin(); iter != MarlParts.end(); iter++ )
    marl_trackids.insert( iter->first );// Contains a list of Marly TrIDs
  std::set< int > apa_trackids;
  for ( std::map<int,simb::MCParticle>::iterator iter = APAParts.begin(); iter != APAParts.end(); iter++ )
    apa_trackids.insert( iter->first );// Contains a list of APA TrIDs
  std::set< int > cpa_trackids;
  for ( std::map<int,simb::MCParticle>::iterator iter = CPAParts.begin(); iter != CPAParts.end(); iter++ )
    cpa_trackids.insert( iter->first );// Contains a list of CPA TrIDs
  std::set< int > Ar39_trackids;
  for ( std::map<int,simb::MCParticle>::iterator iter = Ar39Parts.begin(); iter != Ar39Parts.end(); iter++ )
    Ar39_trackids.insert( iter->first );// Contains a list of Ar39 TrIDs
  std::set< int > neut_trackids;
  for ( std::map<int,simb::MCParticle>::iterator iter = NeutParts.begin(); iter != NeutParts.end(); iter++ )
    neut_trackids.insert( iter->first );// Contains a list of Neutron TrIDs
  std::set< int > kryp_trackids;
  for ( std::map<int,simb::MCParticle>::iterator iter = KrypParts.begin(); iter != KrypParts.end(); iter++ )
    kryp_trackids.insert( iter->first );// Contains a list of Krypton TrIDs
  std::set< int > plon_trackids;
  for ( std::map<int,simb::MCParticle>::iterator iter = PlonParts.begin(); iter != PlonParts.end(); iter++ )
    plon_trackids.insert( iter->first );// Contains a list of Plon TrIDs
  std::set< int > rdon_trackids;
  for ( std::map<int,simb::MCParticle>::iterator iter = RdonParts.begin(); iter != RdonParts.end(); iter++ )
    rdon_trackids.insert( iter->first );// Contains a list of Radon TrIDs



  std::set< int > signal_trackids;
  // Needed for marley running
  std::map< int, simb::MCParticle >::iterator iter;
  std::cout << "We're starting the marley bits " << std::endl;
  nNeut = 0;
  neutCaptureTime = -1e3;
  int neutTrkId = 0;
  for ( size_t i = 0; i < MarlAssn.size(); i++) {
    auto parts = MarlAssn.at(i);
    for (auto part = parts.begin(); part != parts.end(); part++) {
      std::cout << "Marley particle " << (*part)->PdgCode() << "  " << (*part)->E() << "  " << (*part)->T() << "  " << (*part)->TrackId() << "  " << (*part)->Mother() << std::endl;
      signal_trackids.emplace((*part)->TrackId());
      if ((*part)->PdgCode()==2112) nNeut++;
      if ((*part)->PdgCode()==2112) neutTrkId = (*part)->TrackId();
    }

    if (nNeut>0){
      for (auto part = parts.begin(); part != parts.end(); part++) {
        if ((*part)->Mother()==neutTrkId) neutCaptureTime = (*part)->T()*1e-3;
      }
    }

  }

  // --- Finally, get a list of all of my particles in one chunk.
  const sim::ParticleList& PartList = pi_serv->ParticleList();
  std::cout << "There are a total of " << PartList.size() << " MCParticles in the event " << std::endl;


  TVector3 primP(0,0,0);
  TVector3 primPE(0,0,0);
  for ( size_t i = 0; i < MarlAssn.size(); i++) {
    auto parts = MarlAssn.at(i);
    for (auto part = parts.begin(); part != parts.end(); part++) {
      //std::cout << "Marley Part " << (*part)->PdgCode() <<"  " << (*part)->Process() << "  " << (*part)->E() << "  " <<  (*part)->Vx() << "  " << (*part)->Vy() << "  " << (*part)->Vz() <<"  " << (*part)->T() <<  std::endl;
      //std::cout << "Initial KE " << (*part)->E()-(*part)->Mass() <<"  " << (*part)->EndE()-(*part)->Mass() << std::endl;
      if ((*part)->Process()=="primary"){
        TVector3 curPrimP((*part)->Px(),(*part)->Py(),(*part)->Pz());
        primP += curPrimP;
        if ((*part)->PdgCode()==11) primPE = curPrimP;
      }
      /*
      const simb::MCTrajectory traj = (*part)->Trajectory();
      for (int pt = 0; pt < int((*part)->NumberTrajectoryPoints()); pt++){
        std::cout << traj.X(pt) << "  " << traj.Y(pt) <<"  " << traj.Z(pt) << "  " << traj.E(pt)-(*part)->Mass() << std::endl;
      }
      */
    }
  }
  enuCos =  primPE.Dot(primP)/sqrt(primP.Mag2()*primPE.Mag2());


  // --- Lift out the APA particles.
  auto APATrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fAPALabel);
  art::FindManyP<simb::MCParticle> APAAssn(APATrue,evt,fGEANTLabel);
  FillMyMaps( APAParts, APAAssn, APATrue );
  TotGen_APA = APAParts.size();
  std::cout << "--- The size of APAParts is " << APAParts.size() << std::endl;

  // --- Lift out the CPA particles.
  auto CPATrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fCPALabel);
  art::FindManyP<simb::MCParticle> CPAAssn(CPATrue,evt,fGEANTLabel);
  FillMyMaps( CPAParts, CPAAssn, CPATrue );
  TotGen_CPA = CPAParts.size();
  std::cout << "--- The size of CPAParts is " << CPAParts.size() << std::endl;

  // --- Lift out the Ar39 particles.
  auto Ar39True = evt.getValidHandle<std::vector<simb::MCTruth> >(fAr39Label);
  art::FindManyP<simb::MCParticle> Ar39Assn(Ar39True,evt,fGEANTLabel);
  FillMyMaps( Ar39Parts, Ar39Assn, Ar39True );
  TotGen_Ar39 = Ar39Parts.size();
  std::cout << "--- The size of Ar39Parts is " << Ar39Parts.size() << std::endl;

  // --- Lift out the Ar42 particles.
  auto Ar42True = evt.getValidHandle<std::vector<simb::MCTruth> >(fAr42Label);
  art::FindManyP<simb::MCParticle> Ar42Assn(Ar42True,evt,fGEANTLabel);
  FillMyMaps( Ar42Parts, Ar42Assn, Ar42True );
  //TotGen_Ar42 = Ar42Parts.size();
  std::cout << "--- The size of Ar42Parts is " << Ar42Parts.size() << std::endl;

  // --- Lift out the Neut particles.
  auto NeutTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fNeutLabel);
  art::FindManyP<simb::MCParticle> NeutAssn(NeutTrue,evt,fGEANTLabel);
  FillMyMaps( NeutParts, NeutAssn, NeutTrue );
  TotGen_Neut = NeutParts.size();
  std::cout << "--- The size of NeutParts is " << NeutParts.size() << std::endl;
  NeutEnergy = 0;
  NeutPhotEnergy = 0;
  neutE.clear();
  rdonE.clear();

  /*
  for ( size_t i = 0; i < NeutAssn.size(); i++) {
    auto parts = NeutAssn.at(i);
    for (auto part = parts.begin(); part != parts.end(); part++) {
      //if ((*part)->Process()!="primary") continue;
      //std::cout << (*part)->Mother() << "  " << (*part)->Process() << std::endl;
      std::cout << "Found a neutron Particle " << (*part)->PdgCode() << "  " << (*part)->Process() << "  " << (*part)->E() << "  " << (*part)->Vx() << "   " << (*part)->Vy() << " " << (*part)->Vz() << std::endl;

      if ((*part)->Process()=="nCapture" && (*part)->PdgCode()==22){
        neutE.push_back((*part)->E());
      }

      //if ((*part)->Process()=="primary" || (*part)->Process()=="nCapture"){
      //  const simb::MCTrajectory traj = (*part)->Trajectory();
      //  for (int pt = 0; pt < int((*part)->NumberTrajectoryPoints()); pt++){
      //    std::cout << traj.X(pt) << "  " << traj.Y(pt) <<"  " << traj.Z(pt) << "  " << traj.E(pt)-(*part)->Mass() << std::endl;
      //  }
      // }

      if ((*part)->Process() != "primary" && (*part)->Mother()>=0){
        sim::ParticleList::const_iterator mom_it = PartList.find((*part)->Mother());
        const simb::MCParticle *mom = mom_it->second;
        if (mom->Process()=="primary"){
          std::cout << "Found a daughter " << (*part)->PdgCode() << "  " << (*part)->E() << std::endl;
          NeutEnergy += (*part)->E()-(*part)->Mass();
          if ((*part)->PdgCode()==22) NeutPhotEnergy += (*part)->E();
        }
      }
    }
  }
  */

  /*
  for ( size_t i = 0; i < NeutAssn.size(); i++) {
    auto parts = NeutAssn.at(i);
    for (auto part = parts.begin(); part != parts.end(); part++) {
      std::cout << "NEutron Part " << (*part)->PdgCode() <<"  " << (*part)->Process() << "  " << (*part)->E() << "  " <<  (*part)->Vx() << "  " << (*part)->Vy() << "  " << (*part)->Vz() <<"  " << (*part)->T() <<  std::endl;
      std::cout << "Initial KE " << (*part)->E()-(*part)->Mass() <<"  " << (*part)->EndE()-(*part)->Mass() << std::endl;
    }
  }
      */


  // --- Lift out the Kryp particles.
  auto KrypTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fKrypLabel);
  art::FindManyP<simb::MCParticle> KrypAssn(KrypTrue,evt,fGEANTLabel);
  FillMyMaps( KrypParts, KrypAssn, KrypTrue );
  TotGen_Kryp = KrypParts.size();
  std::cout << "--- The size of KrypParts is " << KrypParts.size() << std::endl;
  /*
  for ( size_t i = 0; i < KrypAssn.size(); i++) {
    auto parts = KrypAssn.at(i);
    for (auto part = parts.begin(); part != parts.end(); part++) {
      std::cout << "Kryp Part " << (*part)->PdgCode() <<"  " << (*part)->E() << "  " << (*part)->Vx() << "  " << (*part)->Vy() << "  " << (*part)->Vz() << std::endl;
    }
  }
  */

  // --- Lift out the Plon particles.
  auto PlonTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fPlonLabel);
  art::FindManyP<simb::MCParticle> PlonAssn(PlonTrue,evt,fGEANTLabel);
  FillMyMaps( PlonParts, PlonAssn, PlonTrue );
  TotGen_Plon = PlonParts.size();
  std::cout << "--- The size of PlonParts is " << PlonParts.size() << std::endl;

  // --- Lift out the Rdon particles.
  auto RdonTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fRdonLabel);
  art::FindManyP<simb::MCParticle> RdonAssn(RdonTrue,evt,fGEANTLabel);
  FillMyMaps( RdonParts, RdonAssn, RdonTrue );
  TotGen_Rdon = RdonParts.size();
  std::cout << "--- The size of RdonParts is " << RdonParts.size() << std::endl;
  for ( size_t i = 0; i < RdonAssn.size(); i++) {
    auto parts = RdonAssn.at(i);
    for (auto part = parts.begin(); part != parts.end(); part++){
      if ((*part)->Process() != "alphaInelastic") continue;

      if ((*part)->PdgCode()==22){
        rdonE.push_back((*part)->E());
      }

      //std::cout << "Radon Part " << (*part)->PdgCode() <<"  " << (*part)->Process() << "  "<< (*part)->E() << "  " << (*part)->Vx() << "  " << (*part)->Vy() << "  " << (*part)->Vz() << std::endl;
    }
  }

  /*
  for ( int i = 0; i < int(flashlist.size()); i++ ){
    recob::OpFlash TheFlash = *flashlist[i];
    std::vector< art::Ptr< recob::OpHit > > matchedHits = OpAssns.at(i);

    // Calculate the flash purity, only for the Marley events
    double purity = pbt->OpHitCollectionPurity(signal_trackids, matchedHits);
    OpFlashMarlPur[i] = purity;
    OpFlashPE[i] = TheFlash.TotalPE();
    OpFlashT[i] = TheFlash.Time();
    OpFlashDeltaT[i] = TheFlash.TimeWidth();
    OpFlashNHit[i] = matchedHits.size();

    for (int j = 0; j < int(matchedHits.size()); j++){
      std::vector< sim::TrackSDP > sdps =
        pbt->OpHitToEveTrackSDPs(matchedHits[j]);
      //std::cout << "Op Hit from flash with " << matchedHits.size() << " hits and " << sdps.size() << std::endl;
      for (int k = 0; k < int(sdps.size()); k++)
        std::cout << k << "  " << sdps[k].trackID << "  " << sdps[k].energyFrac << "  " << sdps[k].energy << std::endl;
    }

  }
*/

  /*
  std::vector< recob::Hit > ColHits_Marl;
  std::vector< recob::Hit > ColHits_CPA;
  std::vector< recob::Hit > ColHits_APA;
  std::vector< recob::Hit > ColHits_Ar39;
  std::vector< recob::Hit > ColHits_Neut;
  std::vector< recob::Hit > ColHits_Kryp;
  std::vector< recob::Hit > ColHits_Plon;
  std::vector< recob::Hit > ColHits_Rdon;
  std::vector< recob::Hit > ColHits_Oth;
  */

  std::vector<recob::Hit> ColHits0;
  std::vector<recob::Hit> ColHits1;
  std::vector<recob::Hit> ColHits2;



  // --- Lift out the reco hits:
  auto reco_hits = evt.getValidHandle<std::vector<recob::Hit> >(fHitLabel);

  //*
  // --- Loop over the reconstructed hits to determine the "size" of each hit 
  NTotHits = reco_hits->size();
  int LoopHits = std::min( NTotHits, nMaxHits );
  std::cout << "---- There are " << NTotHits << " hits in the event, but array is of size " << nMaxHits << ", so looping over first " << LoopHits << " hits." << std::endl;
  int totIdx = 0;
  marlEDep = 0;
  
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);   // NOVO
  
  for(int hit = 0; hit < NTotHits; ++hit) {
    totIdx++;
    //for(int hit = 0; hit < LoopHits; ++hit) {
    // --- Let access this particular hit.
    recob::Hit const& ThisHit = reco_hits->at(hit);
    //if (ThisHit.View() != 2) continue;

    // --- Lets figure out which particle contributed the most charge to this hit...
    int MainTrID    = -1;
    double TopEFrac = -DBL_MAX;
    std::vector< sim::TrackIDE > ThisHitIDE = bt_serv->HitToTrackIDEs(clockData,ThisHit);             // Erro
    for (size_t ideL=0; ideL < ThisHitIDE.size(); ++ideL) {
      if ( ThisHitIDE[ideL].energyFrac > TopEFrac ) {
	TopEFrac = ThisHitIDE[ideL].energyFrac;
	MainTrID = ThisHitIDE[ideL].trackID;
      }
    }
    // --- Lets figure out how that particle was generated...
    PType ThisPType = WhichParType( MainTrID );
    if (ThisPType == 1 && ThisHit.View()==2){
      std::cout << "Marley reco hit  " << ThisHit.SummedADC() << "  " << ThisHit.View() << "  " << ThisHit.Channel() << std::endl;
      marlEDep+=ThisHit.Integral();

      hitT = 0.5*(ThisHit.EndTick() + ThisHit.StartTick());
      hitCharge = ThisHit.Integral();
      const geo::WireGeo* wire = geo->GeometryCore::WirePtr(ThisHit.WireID());
      double hXYZ[3];
      wire->GetCenter(hXYZ);
      hitZ = hXYZ[2];

      fMarlHitTree->Fill();
    }

    /*
    // --- Write out some information about this hit....
    std::cout << "Looking at hit on channel " << ThisHit.Channel() << " corresponding to TPC " << ThisHit.WireID().TPC << ", wire " << ThisHit.WireID().Wire << ", plane " << ThisHit.WireID().Plane << ".\n"
	      << "\tIt was at time " << ThisHit.PeakTime() << ", with amplitude " << ThisHit.PeakAmplitude() << ", it was caused by " << ThisHitIDE.size() << " particles, the main one being"
	      << " TrackID " << MainTrID << " which was generated by " << ThisPType
	      << std::endl;
    //*/
    // --- Check which view this hit is on...
    if(ThisHit.View() == geo::kU || ThisHit.View() == geo::kV) {
      ++NIndHits;
    } else { // If not induction then must be collection.
      ++NColHits;
    }

    // --- Now fill in all of the hit level variables.
    /*
    HitView[hit] = ThisHit.View();
    HitSize[hit] = ThisHit.EndTick() - ThisHit.StartTick();
    HitTPC [hit] = ThisHit.WireID().TPC;
    HitChan[hit] = ThisHit.Channel();
    HitTime[hit] = ThisHit.PeakTime();
    HitRMS [hit] = ThisHit.RMS();
    HitSADC[hit] = ThisHit.SummedADC();
    HitInt [hit] = ThisHit.Integral();
    HitPeak[hit] = ThisHit.PeakAmplitude();
    // HitFrame[hit]= ThisHit.Frame();
    GenType[hit] = ThisPType;
    */


    hitCharge = ThisHit.Integral();
    //if (hitCharge < 40) continue;

    // --- I want to fill a vector of coll plane hits, for each of the different kinds of generator.
    if (ThisHit.View() == 2) {
      /*
      if (ThisPType == 0)      ColHits_Oth .push_back( ThisHit );
      else if (ThisPType == 1) ColHits_Marl.push_back( ThisHit );
      else if (ThisPType == 2) ColHits_APA .push_back( ThisHit );
      else if (ThisPType == 3) ColHits_CPA .push_back( ThisHit );
      else if (ThisPType == 4) ColHits_Ar39.push_back( ThisHit );
      else if (ThisPType == 5) ColHits_Neut.push_back( ThisHit );
      else if (ThisPType == 6) ColHits_Kryp.push_back( ThisHit );
      else if (ThisPType == 7) ColHits_Plon.push_back( ThisHit );
      else if (ThisPType == 8) ColHits_Rdon.push_back( ThisHit );
      */
      ColHits2.push_back( ThisHit ); // Finally add every col plane hit to here


      hitT = 0.5*(ThisHit.EndTick() + ThisHit.StartTick());
      //hitT = ThisHit.PeakTime();
      const geo::WireGeo* wire = geo->GeometryCore::WirePtr(ThisHit.WireID());
      double hXYZ[3];
      wire->GetCenter(hXYZ);
      hitZ = hXYZ[2];

      fracE = 0;
      fracGa= 0;
      fracPr= 0;

      //hitPur = ThisPType==6 ? 1 : 0; // neut
      //hitPur = ThisPType==9 ? 1 : 0; // radon
      hitPur = ThisPType==1 ? 1 : 0; // marley
      if (hitPur==0) continue;

      const std::vector<sim::TrackIDE> ides = bt_serv->HitToTrackIDEs(clockData,ThisHit);           // Erro

      if (ThisPType == -1000){
        //std::cout << "Starting a hit " << ThisHit.Integral() << "  " << ThisPType << std::endl;
        for (sim::TrackIDE ide : ides){
          int tID = ide.trackID;
          sim::ParticleList::const_iterator part_it = PartList.find(tID);
          if (part_it == PartList.end()) continue;
          const simb::MCParticle *part = part_it->second;
          if (abs(part->PdgCode())==11) fracE += ide.energyFrac;
          if (part->PdgCode()==22)      fracGa+= ide.energyFrac;
          if (part->PdgCode()==2212)    fracPr+= ide.energyFrac;
          //std::cout << part->PdgCode() << "  " << ide.energyFrac << "  " << tID << "  " << part->E() << "  " << part->Vx() << "  " << part->Vy() << "  " << part->Vz() << "  " << part->T() << std::endl;
          //std::cout << part->Mother() << "  " << part->Process() << std::endl;
          if (part->Process() != "primary" && part->Mother()>=0){
            sim::ParticleList::const_iterator mom_it = PartList.find(part->Mother());
            const simb::MCParticle *mom = mom_it->second;
            if (mom->PdgCode()==0) std::cout <<"Mother  " <<  mom->PdgCode() << "  " << mom->Process() << "  " << mom->E() << "  " << mom->T() << std::endl;
          }
        }
      }

      const std::vector<int> trks = bt_serv->HitToTrackIds(clockData,ThisHit);     // Erro
      if (trks.size()==0) continue;

      //const simb::MCParticle *part = pi_serv->TrackIdToParticle_P(trks[0]);
      //std::cout << part->PdgCode() <<"  " << hitCharge << "  " << hitZ << "   " << hitT <<  " " << part->E() << "   " << part->Vx() << "  " << part->Vy() << " " << part->Vz() << "  " << part->T() << "  " << part->Process() << std::endl;

      fHitTree->Fill();

      //std::cout << ThisPType << "  " << hitCharge << std::endl;
    }


    if (ThisHit.View() == 0)
      ColHits0.push_back( ThisHit ); // Finally add every col plane hit to here
    if (ThisHit.View() == 1)
      ColHits1.push_back( ThisHit ); // Finally add every col plane hit to here

    /*
    if (ThisPType == 1){
      const std::vector<sim::IDE> hits = bt_serv->HitToAvgSimIDEs(ThisHit);
      if (hits.size() == 0) continue;
      std::vector<geo::WireID> wires =
        geo->GeometryCore::ChannelToWire(ThisHit.Channel());
      for (geo::WireID wireid : wires){
        const geo::WireGeo* wire = geo->GeometryCore::WirePtr(wireid);
        double hXYZ[3];
        wire->GetCenter(hXYZ);
        std::cout << "Wire center x,y,z:     "<< hXYZ[0] << ", " << hXYZ[1] << ", " << hXYZ[2] << "   " << std::endl;
        wire->GetStart(hXYZ);
        std::cout << "Wire start x,y,z:      "<< hXYZ[0] << ", " << hXYZ[1] << ", " << hXYZ[2] << "   " << std::endl;
        wire->GetEnd(hXYZ);
        std::cout << "Wire end x,y,z:        "<< hXYZ[0] << ", " << hXYZ[1] << ", " << hXYZ[2] << "   "  << std::endl;
        std::cout << "BackTracked hit x,y,z: " << hits[0].x << ",  " << hits[0].y << ",  " << hits[0].z << std::endl;
      }
    }
    */


  } // Loop over reco_hits.

  std::vector<recob::OpFlash> FlashCands;


  for ( int i = 0; i < int(flashlist.size()); i++ ){
    recob::OpFlash TheFlash = *flashlist[i];

    FlashCands.push_back(TheFlash);
  }
/*
    std::cout << "Flash " << i << "  " << TheFlash.Time() << "  " << TheFlash.TotalPE() << std::endl;
    std::vector< art::Ptr< recob::OpHit > > matchedHits = OpAssns.at(i);


    // Calculate the flash purity, only for the Marley events
    double purity = pbt->OpHitCollectionPurity(signal_trackids, matchedHits);
    OpFlashMarlPur[i] = purity;
    OpFlashPE[i] = TheFlash.TotalPE();
    OpFlashT[i] = TheFlash.Time();
    OpFlashDeltaT[i] = TheFlash.TimeWidth();
    OpFlashNHit[i] = matchedHits.size();
    OpFlashY[i] = TheFlash.YCenter();
    OpFlashdY[i] = TheFlash.YWidth();
    OpFlashZ[i] = TheFlash.ZCenter();
    OpFlashdZ[i] = TheFlash.ZWidth();


    //if (purity > 0.5) std::cout << "We have a  Marley OpFlash " << TheFlash.TotalPE() << "  " << TheFlash.Time() << "  " << TheFlash.YCenter()-TrueY << "  " << TheFlash.ZCenter()-TrueZ << std::endl;

    for (int j = 0; j < int(matchedHits.size()); j++){
      std::vector< sim::TrackSDP > sdps =
        pbt->OpHitToEveTrackSDPs(matchedHits[j]);
      std::cout << "Op Hit from flash with " << matchedHits.size() << " hits and " << sdps.size() << " " <<  matchedHits[j]->Frame() << std::endl;
      for (int k = 0; k < int(sdps.size()); k++){
        std::cout << k << "  " << sdps[k].trackID << "  " << sdps[k].energyFra  << "  " << sdps[k].energy << std::endl;
        PType ThisPType = WhichParType( sdps[k].trackID );
        std::cout << "and this is a " << ThisPType << " particle" << std::endl;
      }
}

    */
  //std::cout << "We have " << FlashCands.size() << " flash candidates" << std::endl;
  //std::cout << "FlashCands size " << FlashCands.size() << std::endl;



  // --- Now calculate all of the hits...
  std::vector< std::vector<recob::Hit> > Clusters0;
  CalcAdjHits( ColHits0, Clusters0,
               hAdjHits, hAdjHitsADCInt, false );
  std::vector< std::vector<recob::Hit> > Clusters1;
  CalcAdjHits( ColHits1, Clusters1,
               hAdjHits, hAdjHitsADCInt, false );
  std::vector< std::vector<recob::Hit> > Clusters2;
  CalcAdjHits( ColHits2, Clusters2,
               hAdjHits, hAdjHitsADCInt, false );

  //for (std::vector< std::vector<recob::Hit> > Clusters : {Clusters2}){
  //std::cout << "Vec sizes " << Clusters0.size() << "  " << ColHits0.size() << "  " << Clusters1.size() << "  " << ColHits1.size() << "  " << Clusters2.size() << "  " << ColHits2.size() << std::endl;
  std::vector<double> BiggestT;
  std::vector< std::vector< std::vector<recob::Hit> > > Clusterss = {Clusters0,Clusters1,Clusters2};
  for (int idx = 0; idx < 3; idx++){
    std::vector< std::vector<recob::Hit> > Clusters = Clusterss[idx];
    BiggestT.push_back(-1e4);
    double BiggestCharge = -1e4;

    for (int i = 0; i < int(Clusters.size()); i++){


      nPur.clear();
      nNHit.clear();
      ndZ.clear();
      ndT.clear();
      ndY.clear();
      nR.clear();
      nCharge.clear();

      //nhit = Clusters[i].size();
      nhit = 0;
      nhit50 = 0;
      std::vector<int> curChans;
      std::vector<int> curPeaks;
      for (recob::Hit hit : Clusters[i]){
        bool contributes = true;
        for (int i = 0; i < int(curChans.size()); i++){
          if (int(hit.Channel()) == int(curChans[i]) &&
              int(hit.PeakTime()) == int(curPeaks[i]))
            contributes = false;
        }
        if (!contributes) continue;
        curChans.push_back(hit.Channel());
        curPeaks.push_back( int(hit.PeakTime()) );
        nhit++;
        if (hit.Integral()>50) nhit50++;
      }
      //std::cout << idx << "  " << i << "  " << nhit << "  " << Clusters[i].size() << std::endl;

      if (nhit < 3) continue;
      clustY = 0;
      clustZ = 0;
      clustT = 0;
      charge = 0;
      marlPur = 0;
      apaPur = 0;
      cpaPur = 0;
      Ar39Pur = 0;
      Ar42Pur = 0;
      neutPur = 0;
      krypPur = 0;
      plonPur = 0;
      rdonPur = 0;
      OpPE = -1;
      OpHiFrac = -1;
      OpZ = -1e4;
      OpY = -1e4;
      OpT = -1e4;
      OpdT = -1e4;
      OpPur = 0;
      HiEEnergy = 0;

      edirx = -5;
      ediry = -5;
      edirz = -5;

      int ourTPC = Clusters[i][0].WireID().TPC;

      int xsign = 0;

      maxHit = 0;
      for (recob::Hit hit : Clusters[i]){
        view = hit.View(); // same for all hits in current view
        //frame = hit.Frame();
        charge += hit.Integral();

        if (hit.Integral()>maxHit) maxHit = hit.Integral();

        const geo::WireGeo* wire = geo->GeometryCore::WirePtr(hit.WireID());
        double dyds = wire->Direction()[1];
        double dzds = wire->Direction()[2];
        double hXYZ[3];
        wire->GetCenter(hXYZ);
        clustY += hit.Integral() * hXYZ[1];
        clustZ += hit.Integral() * hXYZ[2];
        clustT += hit.Integral() * hit.PeakTime();
        dzdy = dzds/dyds; // same for all hits in current view
        const std::vector<int> trks = bt_serv->HitToTrackIds(clockData,hit);             // Erro
        if (trks.size()==0) continue;


        sim::ParticleList::const_iterator part_it = PartList.find(trks[0]);
        if (part_it == PartList.end()) continue;
        const simb::MCParticle *part = part_it->second;
        if (abs(part->PdgCode())==11 && part->E()>HiEEnergy){
          edirx = part->Px()/part->P();
          ediry = part->Py()/part->P();
          edirz = part->Pz()/part->P();
          HiEEnergy = part->E();
        }

        if (part->Process() != "primary"){
          int momID = part->Mother();
          sim::ParticleList::const_iterator mother_it = PartList.find(momID);
          if (mother_it != PartList.end()){
            const simb::MCParticle *mom = mother_it->second;
            MomX = mom->Vx();
            MomY = mom->Vy();
            MomZ = mom->Vz();
          }
        }
        else{
          MomX = part->Vx();
          MomY = part->Vy();
          MomZ = part->Vz();
        }
        PartX = part->Vx();
        PartY = part->Vy();
        PartZ = part->Vz();


        PType ThisPType = WhichParType( trks[0] );
        if (ThisPType == 1)
          marlPur += hit.Integral();
        if (ThisPType == 2)
          apaPur += hit.Integral();
        if (ThisPType == 3)
          cpaPur += hit.Integral();
        if (ThisPType == 4)
          Ar39Pur += hit.Integral();
        if (ThisPType == 5)
          Ar42Pur += hit.Integral();
        if (ThisPType == 6)
          neutPur += hit.Integral();
        if (ThisPType == 7)
          krypPur += hit.Integral();
        if (ThisPType == 8)
          plonPur += hit.Integral();
        if (ThisPType == 9)
          rdonPur += hit.Integral();

        if (view==2){
          const geo::WireGeo* wire = geo->GeometryCore::WirePtr(hit.WireID());
          double hXYZ[3];
          wire->GetCenter(hXYZ);
          xsign = hXYZ[0] > 0 ? 1 : -1;
        }

      }
      clustY /= charge;
      clustZ /= charge;
      clustT /= charge;
      marlPur /= charge; apaPur  /= charge; cpaPur  /= charge; Ar39Pur/=charge;
      Ar42Pur /= charge;
      neutPur /= charge; krypPur /= charge; plonPur /= charge; rdonPur/=charge;

      if (charge > BiggestCharge){
        BiggestT.back() = clustT;
        BiggestCharge = charge;
      }




      // We'll need to find the nearest cluster to current
      /*
      double minD = 1e6;
      for (int j = 0; j < int(Clusters.size()); j++){
        if (j == i) continue;
        double curPE = 0;
        double curZ = 0;
        double curT = 0;
        double curPur = 0;
        for (recob::Hit hit : Clusters[j]){
          curPE += hit.Integral();
          const geo::WireGeo* wire = geo->GeometryCore::WirePtr(hit.WireID());
          double hXYZ[3];
          wire->GetCenter(hXYZ);
          curZ += hit.Integral() * hXYZ[2];
          curT += hit.Integral() * hit.PeakTime();

          const std::vector<int> trks = bt_serv->HitToTrackIds(clockData,hit);
          if (trks.size()==0) continue;
          PType ThisPType = WhichParType( trks[0] );
          if (ThisPType == 1) curPur += hit.Integral();
        }
        curZ /= curPE;
        curT /= curPE;
        curPur /= curPE;

        if (std::abs(curT - clustT) > 200) continue;

        double curD = sqrt(pow(curZ-clustZ,2));
        if (curD < minD){
          nCharge = curPE;
          ndZ = curZ - clustZ;
          ndT = curT - clustT;
          nPur = curPur;
          nNHit = Clusters[j].size();
          minD = curD;
        }
      }

      */
      //if (marlPur > 0) std::cout << "Marl pur " << view << "  " << charge << "  " << nhit << "  " << clustT << "  " << std::endl;


      //std::cout << "Finding match for " << NuEnergy << std::endl;
      if (idx == 2){

        for (recob::Hit rhit : Clusters[i])
          std::cout << rhit.View() << "  " << rhit.Integral() << "  " << rhit.PeakTime() << "  " << rhit.Channel() << std::endl;



        m1Charge = -1e3;
        m1Pur = 0;
        m1dzdy = -1e3;
        m1NHit = 0;
        m150NHit = 0;
        m1T = -1e4;
        m1Z = -1e4;
        m1Y = -1e4;
        int ourTPC0 = 0;
        int m1Idx = -1;
        for (int tIdx = 0; tIdx < int(Clusterss[0].size()); tIdx++){
          std::vector<recob::Hit> testc = Clusterss[0][tIdx];
          if (int(testc.size())<2) continue;

          double mCharge = 0;
          double mClustT = 0;
          double mClustY = 0;
          double mClustZ = 0;
          double mPur = 0;
          int hiIdx = 0;
          for (int hIdx = 0; hIdx < int(testc.size()); hIdx++){
            mCharge += testc[hIdx].Integral();
            mClustT += testc[hIdx].Integral() * testc[hIdx].PeakTime();

            const std::vector<int> trks = bt_serv->HitToTrackIds(clockData,testc[hIdx]);      // Erro
            if (trks.size()==0) continue;
            PType ThisPType = WhichParType( trks[0] );
            if (ThisPType == 1)
              mPur += testc[hIdx].Integral();

            if (testc[hIdx].Integral() > testc[hiIdx].Integral())
              hiIdx = hIdx;
          }
          mClustT /= mCharge;
          if (abs(mClustT - clustT)>100) continue;

          for (int hIdx = hiIdx; hIdx <= hiIdx; hIdx++){
            recob::Hit hit = testc[hIdx];
            int wireID = -1;
            std::vector<geo::WireID> wires =
              geo->GeometryCore::ChannelToWire(hit.Channel());
            for (int i = 0; i < int(wires.size()); i++){
              const geo::WireGeo* wire = geo->GeometryCore::WirePtr(wires[i]);

              double strtXYZ[3];
              wire->GetStart(strtXYZ);

              double stopXYZ[3];
              wire->GetEnd(stopXYZ);

              if (strtXYZ[0] * xsign > 0 &&
                  ( (strtXYZ[2] < clustZ && stopXYZ[2] > clustZ) ||
                    (strtXYZ[2] > clustZ && stopXYZ[2] < clustZ)) )
                wireID = i;
            }
            if (wireID == -1) continue;
            const geo::WireGeo* wire=geo->GeometryCore::WirePtr(wires[wireID]);

            double dyds = wire->Direction()[1];
            double dzds = wire->Direction()[2];
            double hXYZ[3];
            wire->GetCenter(hXYZ);

            mClustY = hXYZ[1];
            mClustZ = hXYZ[2];
            //mClustY += hit.Integral() * hXYZ[1];
            //mClustZ += hit.Integral() * hXYZ[2];
            m1dzdy = dzds/dyds; // same for all hits in current view
            const std::vector<int> trks = bt_serv->HitToTrackIds(clockData,hit);
            if (trks.size()==0) continue;
            PType ThisPType = WhichParType( trks[0] );
            if (ThisPType == 1)
              mPur += hit.Integral();

            wire->GetStart(hXYZ);

          }
          //mClustT /= mCharge;
          if (clustT - mClustT < 25 && clustT - mClustT > -5){
            m1Idx = tIdx;
            m1NHit = int(testc.size());
            for (recob::Hit hit : testc) if (hit.Integral()>50) m150NHit++;
            m1T = mClustT;
            m1Y = mClustY;///mCharge;
            m1Z = mClustZ;///mCharge;
            m1Charge = mCharge;
            m1Pur = mPur/mCharge;
            ourTPC0 = testc[0].WireID().TPC;
          }
        }


        m2Charge = -1e3;
        m2Pur = 0;
        m2dzdy = -1e3;
        m2NHit = 0;
        m250NHit =0;
        m2T = -1e4;
        m2Z = -1e4;
        m2Y = -1e4;
        int ourTPC1 = 0;
        int m2Idx = 0;
        if (m2Idx==-1) std::cout <<"Here " << std::endl;
        for (int tIdx = 0; tIdx < int(Clusterss[1].size()); tIdx++){
          std::vector<recob::Hit> testc = Clusterss[1][tIdx];
          if (int(testc.size())<2) continue;
          if (testc.size() == 0) continue;
          double mCharge = 0;
          double mClustT = 0;
          double mClustY = 0;
          double mClustZ = 0;
          double mPur = 0;
          int hiIdx = 0;
          for (int hIdx = 0; hIdx < int(testc.size()); hIdx++){
            mCharge += testc[hIdx].Integral();
            mClustT += testc[hIdx].Integral() * testc[hIdx].PeakTime();

            const std::vector<int> trks = bt_serv->HitToTrackIds(clockData,testc[hIdx]);            // Erro
            if (trks.size()==0) continue;
            PType ThisPType = WhichParType( trks[0] );
            if (ThisPType == 1)
              mPur += testc[hIdx].Integral();

            if (testc[hIdx].Integral() > testc[hiIdx].Integral())
              hiIdx = hIdx;
          }
          mClustT /= mCharge;
          if (abs(mClustT - clustT)>100) continue;

          for (int hIdx = hiIdx; hIdx <= hiIdx; hIdx++){
            recob::Hit hit = testc[hIdx];
            int wireID = -1;
            std::vector<geo::WireID> wires =
              geo->GeometryCore::ChannelToWire(hit.Channel());
            for (int i = 0; i < int(wires.size()); i++){
              const geo::WireGeo* wire = geo->GeometryCore::WirePtr(wires[i]);

              double strtXYZ[3];
              wire->GetStart(strtXYZ);

              double stopXYZ[3];
              wire->GetEnd(stopXYZ);

              if (strtXYZ[0] * xsign > 0 &&
                  ( (strtXYZ[2] < clustZ && stopXYZ[2] > clustZ) ||
                    (strtXYZ[2] > clustZ && stopXYZ[2] < clustZ)) )
                wireID = i;
            }
            if (wireID == -1) continue;
            const geo::WireGeo* wire=geo->GeometryCore::WirePtr(wires[wireID]);

            double dyds = wire->Direction()[1];
            double dzds = wire->Direction()[2];
            double hXYZ[3];
            wire->GetCenter(hXYZ);

            mClustY = hXYZ[1];
            mClustZ = hXYZ[2];
            //mClustY += hit.Integral() * hXYZ[1];
            //mClustZ += hit.Integral() * hXYZ[2];
            m2dzdy = dzds/dyds; // same for all hits in current view

            wire->GetStart(hXYZ);

          }
          //mClustT /= mCharge;
          if (clustT - mClustT < 20 && clustT - mClustT > -10){
            for (recob::Hit hit : testc) if (hit.Integral()>50) m250NHit++;
            m2Idx = tIdx;
            m2T = mClustT;
            m2NHit = int(testc.size());
            m2Y = mClustY;///mCharge;
            m2Z = mClustZ;///mCharge;
            m2Charge = mCharge;
            m2Pur = mPur/mCharge;
            ourTPC1 = testc[0].WireID().TPC;
          }
        }

        std::vector< std::vector<recob::Hit> > cands0;
        for (int tIdx = 0; tIdx < int(Clusterss[0].size()); tIdx++){
          if (tIdx == m1Idx) continue;
          if (Clusterss[0][tIdx].size() == 0) continue;
          int curTPC = Clusterss[0][tIdx][0].WireID().TPC;
          if (curTPC != ourTPC0) continue;
          cands0.push_back(Clusterss[0][tIdx]);
        }
        std::vector< std::vector<recob::Hit> > cands1;
        for (int tIdx = 0; tIdx < int(Clusterss[1].size()); tIdx++){
          if (tIdx == m1Idx) continue;
          if (Clusterss[1][tIdx].size() == 0) continue;
          int curTPC = Clusterss[1][tIdx][0].WireID().TPC;
          if (curTPC != ourTPC1) continue;
          cands1.push_back(Clusterss[1][tIdx]);
        }
        //std::cout << "Size In each " << cands0.size() << "  " << cands1.size() << std::endl;
        //std::cout << m1Idx <<"  " << m2Idx << "  " << std::endl;

        double fitY02 = -1e4;
        double fitZ02 = -1e4;
        if (fitZ02==0) std::cout << "Here " << std::endl;
        if (m1T > -1e3){
          fitY02 = m1Y + (clustZ - m1Z)/m1dzdy;
          fitZ02 = clustZ;
        }

        double fitY12 = -1e4;
        double fitZ12 = -1e4;
        if (fitZ12==0) std::cout << "Here " << std::endl;
        if (m2T > -1e3){
          fitY12 = m2Y + (clustZ - m2Z)/m2dzdy;
          fitZ12 = clustZ;
        }
        //std::cout << m1Pur << "  " << m2Pur << std::endl;
        //std::cout << m1Charge <<  "  " << m2Charge << std::endl;
        //std::cout << m1T <<  "  " << m2T << std::endl;
        //std::cout << "Fits " << fitY02 << "  " << fitZ02 << " " << fitY12 << "  " << fitZ12 << std::endl;
        //if (m1Charge > m2Charge)
          //std::cout << TrueY << "  " << TrueZ <<"   " <<  fitY02 <<"  " << fitZ02 << "  " << std::endl;
        //if (m1Charge < m2Charge)
        //std::cout << TrueY << "  " << TrueZ <<"   " <<  fitY12 <<"  " << fitZ12 << "  " << std::endl;

        if (m1Charge > 0 && m2Charge > 0){
          fitY = 0.5*(fitY12+fitY02);
        }
        else if (m1Charge > 0){
          fitY = fitY02;
        }
        else if (m2Charge > 0){
          fitY = fitY12;
        }
        else{
          fitY = -1e4;
        }



        ntotCharge = 0;
        nN = 0;

        for (int j = 0; j < int(Clusters.size()); j++){
          if (i == j) continue;
          if (Clusters[j].size() == 0) continue;

          recob::Hit tHit = Clusters[j][0];
          const geo::WireGeo* wire = geo->GeometryCore::WirePtr(tHit.WireID());
          double hXYZ[3];
          wire->GetCenter(hXYZ);
          double curZ = hXYZ[2];
          double curT = tHit.PeakTime();

          double dist = pow((curT-clustT)*0.08,2) + pow(curZ-clustZ,2);

          dist = sqrt(dist);
          if (dist < 300){
            double curCharge = 0;
            for (recob::Hit ttHit : Clusters[j]) curCharge += ttHit.Integral();

            int curTPC = Clusters[j][0].WireID().TPC;
            if (curTPC != ourTPC) continue;

            double curY1 = -20000;
            double curY2 = -20000;

            ViewMatchCluster(cands0,curT,curZ,
                             10,xsign,curY1,curZ,1);
            ViewMatchCluster(cands1,curT,curZ,
                             5, xsign,curY2,curZ,1);

            double curY = curY1;
            int n0Idx = ViewMatchCluster(cands0,curT,curZ,
                                         10,xsign,curY1,curZ,1);
            int n1Idx = ViewMatchCluster(cands1,curT,curZ,
                                         5, xsign,curY2,curZ,1);
            if (curY1==-20000) curY = curY2;
            if (curY1>-2000 && curY2>-2000){
              curY = 0.5*(curY1+curY2);
              if (abs(curY1-curY2)>12){
                double cur0Charge = 0;
                for (recob::Hit curhit : cands0[n0Idx])
                  cur0Charge += curhit.Integral();
                double cur1Charge = 0;
                for (recob::Hit curhit : cands1[n1Idx])
                  cur1Charge += curhit.Integral();
                curY = cur0Charge>cur1Charge ? curY1 : curY2;
              }
            }
            dist = sqrt(dist*dist + (curY-fitY)*(curY-fitY));
            if (dist>200) continue;
            ntotCharge += curCharge;

            std::cout << "Foudnd a match " << dist <<  "  " << curCharge << std::endl;
            //std::cout << n0Idx << "  " << n1Idx << "  " << curY1 << "  " << curY2 << std::endl;

            bool isMarl = false;
            const std::vector<int> trks = bt_serv->HitToTrackIds(clockData,Clusters[j][0]);                       // Erro
            if (trks.size()>0){
              PType curPType = WhichParType(trks[0]);
              isMarl = int(curPType)==1 || int(curPType)==0;
              //std::cout << WhichParType(trks[0]) << std::endl;;
              sim::ParticleList::const_iterator cur_it=PartList.find(trks[0]);
              if (cur_it != PartList.end()){
                const simb::MCParticle *part = cur_it->second;
                if (part->PdgCode()==-1) std::cout <<"Part " << part->Vx() <<"  " << part->Vy() <<"  " << part->Vz() << std::endl;
              }
            }
            nPur.push_back(isMarl);
            ndY.push_back(curY-fitY);
            nR.push_back(dist);
            ndZ.push_back(curZ-clustZ);
            ndT.push_back(curT-clustT);
            nCharge.push_back(curCharge);
            nNHit.push_back(Clusters[j].size());


          }
        }
        nN = nPur.size();


        // Check to see if we have a matched OpFlash for this cluster
        // If there are multiple, take the largest PE cluster
        int Opidx = -1;
        std::cout << "Start flash match " << FlashCands.size() <<  std::endl;
        for (int i = 0; i < int(FlashCands.size()); i++){

          double r = sqrt(pow(FlashCands[i].YCenter()-fitY,2)+pow(FlashCands[i].ZCenter()-clustZ,2));

          if (FlashCands[i].TotalPE()>15)
            std::cout << fracHiPE[i] << "  " << FlashCands[i].ZCenter()-clustZ << "  " << FlashCands[i].YCenter()-fitY << "  " << clustT-2*FlashCands[i].Time() << " " <<FlashCands[i].Time() << "  " << FlashCands[i].TotalPE() << std::endl;

          /*
          if ((FlashCands[i].TotalPE()>200||fracHiPE[i]<0.4)){
            std::cout << "Here 1" << std::endl;
            if (FlashCands[i].TotalPE()>=10){
              std::cout << "Here 2" << std::endl;
              if (std::abs(FlashCands[i].ZCenter()-clustZ) < 300){
                std::cout << "Here 3" << std::endl;
                if (std::abs(FlashCands[i].YCenter()-fitY) < 200){
                  std::cout << "Here 4" << std::endl;
                  if (clustT-2*FlashCands[i].Time()<4400 && clustT-2*FlashCands[i].Time()>=0){
                    std::cout << Opidx << std::endl;
                  }
                }
              }
            }
          }
          */


          double dt = clustT - 2 * FlashCands[i].Time();
          //std::abs(FlashCands[i].ZCenter()-clustZ) < 300 &&
          if ((FlashCands[i].TotalPE()>200 || fracHiPE[i]<0.4) &&
              FlashCands[i].TotalPE() >= 10 &&
              r<250 && dt>0 && dt<4400 &&
              FlashCands[i].TotalPE()>200*exp(-9e-4*dt)){
            if (Opidx == -1)
              Opidx = i;
            else{
              double r0 = sqrt(pow(FlashCands[Opidx].YCenter()-fitY,2)+pow(FlashCands[Opidx].ZCenter()-clustZ,2));

              if (FlashCands[i].TotalPE()/r>FlashCands[Opidx].TotalPE()/r0)
                Opidx = i;
            }
          }
        }
        std::cout << Opidx << std::endl;
        if (Opidx > -1){
          std::vector< art::Ptr<recob::OpHit> > matchedHits=OpAssns.at(Opidx);
          OpHiFrac = 0;
          for (int j = 0; j < int(matchedHits.size()); j++){
            recob::OpHit ohit = *matchedHits[j];
            if (ohit.PE() > OpHiFrac) OpHiFrac = ohit.PE();
          }
          OpHiFrac /= FlashCands[Opidx].TotalPE();
          OpPE = FlashCands[Opidx].TotalPE();
          OpZ = FlashCands[Opidx].ZCenter();
          OpY = FlashCands[Opidx].YCenter();
          OpT = FlashCands[Opidx].Time();
          OpdT= FlashCands[Opidx].TimeWidth();
          OpPur = pbt->OpHitCollectionPurity(signal_trackids, matchedHits);
          //std::cout << "OpPur " << OpPur << std::endl;
          //std::cout << "OpPur " << OpPur << std::endl;
          //std::cout << "OpPur " << OpPur << std::endl;
          //std::cout << "OpPur " << OpPur << std::endl;
          //std::cout << "OpPur " << OpPur << std::endl;
          std::cout << "OpPur " << OpPur << std::endl;
        }
        delayedPE = 0;
        if (OpPE > 20){
          art::Handle< std::vector< recob::OpHit > > HitHandle;
          std::vector<art::Ptr<recob::OpHit> > hitlist;
          if (evt.getByLabel(fOpHitLabel, HitHandle)) {
            art::fill_ptr_vector(hitlist, HitHandle);
          }
          //std::cout << "Fill delayed " << hitlist.size() << std::endl;
          for (int i = 0; i < int(hitlist.size()); i++){
            art::Ptr<recob::OpHit> ohit(HitHandle,i);
            int channel = ohit->OpChannel();
            double xyz[3];
            geo->OpDetGeoFromOpChannel(channel).GetCenter(xyz);;
            double ohy = xyz[1];
            double ohz = xyz[2];
            double r = sqrt(pow(ohy-OpY,2)+pow(ohz-OpZ,2));
            double pe= ohit->PE();
            double dt= ohit->PeakTime()-OpT;
            if (r < 500 && dt < 5 && dt > 0.5)
              delayedPE += pe;

          }
        }



      }

      std::cout << "Filling DAQCluster  " << OpPE  << "  " << OpZ << "  " << delayedPE << std::endl;



      /*
      etrkCos = 0;
      etrkFlippedCos = 0;
      if (idx == 2){
        double hiTrkCharge = -1;
        for(int trk = 0; trk < int(trklist.size()); ++trk) {

          art::Ptr<recob::Track> track(trkHandle,trk);
          auto curRHits = fmthm.at(trk);
          auto curTHM = fmthm.data(trk);
          auto trkDir = track->DirectionAtPoint(curTHM[0]->Index());
          double curCharge = 0;
          double trkZ = 0;
          double trkT = 0;
          for (int tIdx = 0; tIdx < int(curTHM.size()); tIdx++){
            auto xyz2 = track->LocationAtPoint(curTHM[tIdx]->Index());
            std::cout << curRHits[tIdx]->Integral() << "  " << curRHits[tIdx]->View() << "  " << curRHits[tIdx]->Channel() <<"   " << curRHits[tIdx]->PeakTime() <<"  " << xyz2.X() << "  " << xyz2.Y() <<"   " << xyz2.Z() << std::endl;
            curCharge += curRHits[tIdx]->Integral();
            trkZ += xyz2.Z() * curRHits[tIdx]->Integral();
            trkT += curRHits[tIdx]->PeakTime() * curRHits[tIdx]->Integral();
          }
          trkZ /= curCharge;  trkT /= curCharge;
          if (curCharge > hiTrkCharge && abs(clustZ-trkZ)<5&&abs(clustT-trkT)<60){
            hiTrkCharge = curCharge;
            trkcharge = curCharge;
            trknhit = int(curRHits.size());
            trkl = track->Length();
            TVector3 fooTrk(trkDir.X(),trkDir.Y(),trkDir.Z());
            trkcos = fooTrk.Z()/sqrt(fooTrk.Mag2());
            TVector3 foobarTrk(-trkDir.X(),-trkDir.Y(),-trkDir.Z());
            if (fooTrk.Z()>0) foobarTrk = fooTrk;
            etrkCos =  primPE.Dot(fooTrk)/sqrt(fooTrk.Mag2()*primPE.Mag2());
            etrkFlippedCos = primPE.Dot(foobarTrk)/sqrt(foobarTrk.Mag2()*primPE.Mag2());
          }
          //std::cout << "Track Direction " << trkDir.X() << " " << trkDir.Y() <<"   " << trkDir.Z() << "  " << trkDir.Dot(primPE)/sqrt(primPE.Mag2()*trkDir.Mag2()) << "  " << primPE.Dot(primP)/sqrt(primP.Mag2()*primPE.Mag2()) << std::endl;


          //std::cout << track->Start().Z() <<"   " << track->End().Z() << std::endl;
          std::cout << "We found a reco'd track, at " << clustZ-trkZ << "  " << clustT-trkT << "  " << etrkCos << std::endl;
        }
      }

      */


      std::cout << "Attemtp near hit tree " << clustT << "  " << clustZ << std::endl;
      if (OpPE>1){
        for(int hit = 0; hit < NTotHits; ++hit) {
          totIdx++;
          //for(int hit = 0; hit < LoopHits; ++hit) {
          // --- Let access this particular hit.
          recob::Hit const& ThisHit = reco_hits->at(hit);
          if (ThisHit.View()!=2) continue;

          const geo::WireGeo* wire = geo->GeometryCore::WirePtr(ThisHit.WireID());
          double hXYZ[3];
          wire->GetCenter(hXYZ);
          double curZ = hXYZ[2];
          double curT = ThisHit.PeakTime();

          double dist = pow((curT-clustT)*0.08,2) + pow(curZ-clustZ,2);


          dist = sqrt(dist);

          //if (dist<200) std::cout << curT << " " << curZ << "  " << dist << std::endl;

          if (dist < 300){

            bool isMarl = false;
            const std::vector<int> trks = bt_serv->HitToTrackIds(clockData,ThisHit);             // Erro
            if (trks.size()>0){
              PType curPType = WhichParType(trks[0]);
              isMarl = int(curPType)==1 || int(curPType)==0;
              //std::cout << WhichParType(trks[0]) << std::endl;;
              sim::ParticleList::const_iterator cur_it=PartList.find(trks[0]);
              if (cur_it != PartList.end()){
                const simb::MCParticle *part = cur_it->second;
                if (part->PdgCode()==-1) std::cout <<"Part " << part->Vx() <<"  " << part->Vy() <<"  " << part->Vz() << std::endl;
              }
            }

            NhitPur = isMarl ? 1.0 : 0.0;
            NhitZ = curZ-clustZ;
            NhitT = curT-clustT;
            NhitCharge = ThisHit.Integral();
            clustIdx = i;
            //if (dist<200) fNearHitTree->Fill();

          }

        }
      }


      fDAQClusterTree->Fill();
    }
  }

  BiggestT0 = BiggestT[0];
  BiggestT1 = BiggestT[1];
  BiggestT2 = BiggestT[2];

  //std::cout <<" We reconstructed " << trklist.size() <<"   " << std::endl;





  /*
  // Find OpFlashes associated with the event
  art::Handle< std::vector< recob::Cluster > > ClustHandle;
  std::vector<art::Ptr<recob::Cluster> > clustlist;
  if (evt.getByLabel("blurredcluster", ClustHandle)) {
    art::fill_ptr_vector(clustlist, ClustHandle);
    //std::sort(flashlist.begin(), flashlist.end(), recob::OpFlashPtrSortByPE);
  }
  // Grab assns with OpHits to get match to neutrino purity
  art::FindManyP< recob::Hit > ClustHitAssns(clustlist, evt, "blurredcluster");


  for (int i = 0; i < int(clustlist.size()); i++){
    std::vector< art::Ptr<recob::Hit> > matchedHits = ClustHitAssns.at(i);
    std::cout << "Cluster  " << i << std::endl;
    std::cout << matchedHits.size() << "  " << std::endl;
    for (int j = 0; j < int(matchedHits.size()); j++){
      std::cout << j << "  " << matchedHits[j]->View() << " " << matchedHits[j]->Integral() <<"   " << matchedHits[j]->PeakTime() << "  " << matchedHits[j]->Channel() << std::endl;
    }
  }

  auto clust_list = evt.getValidHandle<std::vector<recob::Cluster> >("blurredcluster");
  if (MarlTrue->size() != 0){
    if (MarlTrue->at(0).NParticles() != 0){
      const simb::MCParticle& part = MarlTrue->at(0).GetParticle(0);
      std::cout << "start " << part.E() << std::endl;

      NuEnergy = part.E();
      //const simb::MCParticle& lep = MarlTrue->at(0).GetParticle(1);
      TrueX = part.Vx();
      TrueY = part.Vy();
      TrueZ = part.Vz();

    }
  }
  */


  /*
  */



  /*
  std::cerr << "\nAnd now for APA hits..." << std::endl;
  CalcAdjHits( ColHits_Marl, hAdjHits_Marl, hAdjHitsADCInt_Marl, true );
  std::cerr << "\nAnd now for APA hits..." << std::endl;
  CalcAdjHits( ColHits_APA , hAdjHits_APA , hAdjHitsADCInt_APA, false  );
  std::cerr << "\nAnd now for CPA hits..." << std::endl;
  CalcAdjHits( ColHits_CPA , hAdjHits_CPA , hAdjHitsADCInt_CPA, false  );
  std::cerr << "\nAnd now for Ar39 hits..." << std::endl;
  CalcAdjHits( ColHits_Ar39, hAdjHits_Ar39, hAdjHitsADCInt_Ar39, false );
  std::cerr << "\nAnd now for Neuton hits..." << std::endl;
  CalcAdjHits( ColHits_Neut, hAdjHits_Neut, hAdjHitsADCInt_Neut, false );
  std::cerr << "\nAnd now for Krypton hits..." << std::endl;
  CalcAdjHits( ColHits_Kryp, hAdjHits_Kryp, hAdjHitsADCInt_Kryp, false );
  std::cerr << "\nAnd now for Polonium hits..." << std::endl;
  CalcAdjHits( ColHits_Plon, hAdjHits_Plon, hAdjHitsADCInt_Plon, false );
  std::cerr << "\nAnd now for Radon hits..." << std::endl;
  CalcAdjHits( ColHits_Rdon, hAdjHits_Rdon, hAdjHitsADCInt_Rdon, false );
  std::cerr << "\nAnd now for Other hits..." << std::endl;
  CalcAdjHits( ColHits_Oth , hAdjHits_Oth , hAdjHitsADCInt_Oth, false  );
  */

  //*/
  // --- Now loop through the particle list.
  /*
  std::cout << "\n\nNow to loop through the truth information." << std::endl;
  for ( sim::ParticleList::const_iterator ipar = PartList.begin(); ipar!=PartList.end(); ++ipar) {
    // --- Grab this particle.
    simb::MCParticle *particle = ipar->second;
    // Let's just write out what our primary particles are...
    if (particle->Process() != "primary") continue; // Can also check that particle->Mother() != 0.
    std::cout << "-- Particle with TrackID " << particle->TrackId() << ", which was a " << particle->PdgCode() << " was a primary and had initial energy " << particle->E()
	      << ", " << particle->NumberTrajectoryPoints() << " trajectory points, and " << particle->NumberDaughters() << " daughters, and Process - " << particle->Process()
	      << std::endl;
  }
/  //*/
  /*
  std::vector<short> uADCs;
  for (unsigned int dig=0; dig<rawdigits->size(); ++dig){
    // --- Lets access this particular RawDigit
    raw::RawDigit ThisDig = rawdigits->at(dig);
    nADC=0;
    // --- Uncompress the ADC vector.
    if (dig==0){
      std::cout << "uADCs.size() before Uncompress = " << uADCs.size() << std::endl;
      raw::Uncompress(ThisDig.ADCs(), uADCs, ThisDig.Compression());
      std::cout << "uADCs.size() after Uncompress = " << uADCs.size() << "\n\n";
    }
    // --- Print some stuff about the first RawDigit
    if (dig==0){
      std::cout << "\nLooking at rawdigit["<<dig<<"]. It was on channel " << ThisDig.Channel() << ". "
		<< "It had " << ThisDig.Samples() << " samples. "                                 // The readout length for 1x2x6 is 4492 ticks
		<< "There were a total of " << ThisDig.NADC() << " ADCs saved with compression " // This is the readout length with compression
		<< "level " << ThisDig.Compression()
		<< std::endl;
    }

  } // Loop over RawDigits.
  */
  // --- Finally, fill our TTree once per event.
  std::cout << "Pre fill " << std::endl;
  fDAQSimTree -> Fill();
  std::cout << "Post fill " << std::endl;





} // Analyze SolarNuAna.


//......................................................
void SolarNuAna::FillMyMaps( std::map< int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle< std::vector<simb::MCTruth> > Hand )
{
  for ( size_t L1=0; L1 < Hand->size(); ++L1 ) {
    for ( size_t L2=0; L2 < Assn.at(L1).size(); ++L2 ) {
      const simb::MCParticle ThisPar = (*Assn.at(L1).at(L2));
      MyMap[ThisPar.TrackId()] = ThisPar;
    }
  }
  return;
}

//......................................................
PType SolarNuAna::WhichParType( int TrID )
{
  // Check if Ar42
  if ( InMyMap( TrID, Ar42Parts ) ) {
    return kAr42;
  }
  else if ( InMyMap( TrID, Ar39Parts ) ) {
    return kAr39;
  // Check if MARLEY
  } else  if ( InMyMap( TrID, MarlParts ) ) {
    return kMarl;
  // Check if APA
  } else if ( InMyMap( TrID, APAParts  ) ) {
    return kAPA;
    // Check if CPA
  } else if ( InMyMap( TrID, CPAParts  ) ) {
    return kCPA;
    // Check if Neut
  } else if ( InMyMap( TrID, NeutParts ) ) {
    return kNeut;
    // Check if Kryp
  } else if ( InMyMap( TrID, KrypParts ) ) {
    return kKryp;
    // Check if Plon
  } else if ( InMyMap( TrID, PlonParts ) ) {
    return kPlon;
    // Check if Rdon
  } else if ( InMyMap( TrID, RdonParts ) ) {
    return kRdon;
  }
  // If get this far then who knows???
  return kUnknown;
}

//......................................................
bool SolarNuAna::InMyMap( int TrID, std::map< int, simb::MCParticle> ParMap )
{
  std::map<int, simb::MCParticle>::iterator ParIt;
  ParIt = ParMap.find( TrID );
  if ( ParIt != ParMap.end() ) {
    //std::cout << "In Map " << (ParMap.at(TrID)).second->TrackID() << "  " << (ParMap.at(TrID)).second->PdgCode() < std::endl;
    return true;
  } else
    return false;
}

//......................................................
void SolarNuAna::CalcAdjHits( std::vector< recob::Hit > MyVec,
                             std::vector< std::vector<recob::Hit> >& clusters,
                             TH1I* MyHist, TH1F* ADCIntHist, bool HeavDebug ) {
  const double TimeRange = 10;
  const int    ChanRange = 2;
  unsigned int FilledHits = 0;
  unsigned int NumOriHits = MyVec.size();

  //std::vector< std::vector< recob::Hit > > HitClusters;

  while( NumOriHits != FilledHits ) {
    if (HeavDebug) std::cerr << "\nStart of my while loop" << std::endl;
    std::vector< recob::Hit > AdjHitVec;
    AdjHitVec.push_back ( MyVec[0] );
    MyVec.erase( MyVec.begin()+0 );
    int LastSize = 0;
    int NewSize  = AdjHitVec.size();
    while ( LastSize != NewSize ) {
      std::vector<int> AddNow;
      for (size_t aL=0; aL < AdjHitVec.size(); ++aL) {
	for (size_t nL=0; nL < MyVec.size(); ++nL) {
	  if (HeavDebug) {
	    std::cerr << "\t\tLooping though AdjVec " << aL << " and  MyVec " << nL
		      << " AdjHitVec - " << AdjHitVec[aL].Channel() << " & " << AdjHitVec[aL].PeakTime()
		      << " MVec - " << MyVec[nL].Channel() << " & " << MyVec[nL].PeakTime()
		      << " Channel " << abs( (int)AdjHitVec[aL].Channel()  - (int)MyVec[nL].Channel()  )  << " bool " << (bool)(abs( (int)AdjHitVec[aL].Channel() - (int)MyVec[nL].Channel()  ) <= ChanRange)
		      << " Time " << abs( AdjHitVec[aL].PeakTime() - MyVec[nL].PeakTime() ) << " bool " << (bool)(abs( (double)AdjHitVec[aL].PeakTime() - (double)MyVec[nL].PeakTime() ) <= TimeRange)
		      << std::endl;
	  }
	  if ( abs( (int)AdjHitVec[aL].Channel()  - (int)MyVec[nL].Channel()  ) <= ChanRange &&
	       abs( (double)AdjHitVec[aL].PeakTime() - (double)MyVec[nL].PeakTime() ) <= TimeRange ) {
	    if (HeavDebug) std::cerr << "\t\t\tFound a new thing!!!" << std::endl;
	    // --- Check that this element isn't already in AddNow.
	    bool AlreadyPres = false;
	    for (size_t zz=0; zz<AddNow.size(); ++zz) {
	      if (AddNow[zz] == (int)nL) AlreadyPres = true;
	    }
	    if (!AlreadyPres)
	      AddNow.push_back( nL );
	  } // If this hit is within the window around one of my other hits.
	} // Loop through my vector of colleciton plane hits.
      } // Loop through AdjHitVec
      // --- Now loop through AddNow and remove from Marley whilst adding to AdjHitVec
      for (size_t aa=0; aa<AddNow.size(); ++aa) {
	if (HeavDebug) {
	  std::cerr << "\tRemoving element " << AddNow.size()-1-aa << " from MyVec ===> "
		    << MyVec[ AddNow[AddNow.size()-1-aa] ].Channel() << " & " << MyVec[ AddNow[AddNow.size()-1-aa] ].PeakTime()
		    << std::endl;
	}
	AdjHitVec.push_back ( MyVec[ AddNow[AddNow.size()-1-aa] ] );
	MyVec.erase( MyVec.begin() + AddNow[AddNow.size()-1-aa] );
      }
      LastSize = NewSize;
      NewSize  = AdjHitVec.size();
      if (HeavDebug) {
	std::cerr << "\t---After that pass, AddNow was size " << AddNow.size() << " ==> LastSize is " << LastSize << ", and NewSize is " << NewSize
		  << "\nLets see what is in AdjHitVec...." << std::endl;
	for (size_t aL=0; aL < AdjHitVec.size(); ++aL) {
	  std::cout << "\tElement " << aL << " is ===> " << AdjHitVec[aL].Channel() << " & " << AdjHitVec[aL].PeakTime() << std::endl;
	}
      }

  } // while ( LastSize != NewSize )
    int NumAdjColHits = AdjHitVec.size();
    float SummedADCInt = 0;
    for ( recob::Hit hit : AdjHitVec)
      SummedADCInt += hit.Integral();
    if (HeavDebug) std::cerr << "After that loop, I had " << NumAdjColHits << " adjacent collection plane hits." << std::endl;
    MyHist -> Fill( NumAdjColHits );
    ADCIntHist -> Fill( SummedADCInt );
    FilledHits += NumAdjColHits;

    if (AdjHitVec.size() > 0)
      clusters.push_back(AdjHitVec);

  }

  if (HeavDebug){
    std::vector<double> avgChannel;
    std::vector<double> avgTick;
    std::vector<double> summedADCInt;
    for (std::vector< recob::Hit > hits : clusters){
      double adcInt = 0;
      double channel = 0;
      double tick = 0;
      for (recob::Hit hit : hits){
        tick += hit.Integral()*hit.PeakTime();
        channel += hit.Integral()*hit.Channel();
        adcInt += hit.Integral();
      }
      tick /= adcInt;
      channel /= adcInt;
      summedADCInt.push_back(adcInt);
      avgTick.push_back(tick);
      avgChannel.push_back(channel);
    }

    for (int i = 0; i < int(avgTick.size()-1); i++){
      for (int j = i+1; j < int(avgTick.size()); j++){
        std::cout << avgChannel[i] << " " << avgChannel[j] << "  " << std::abs(avgChannel[i]-avgChannel[j]) << std::endl;
        std::cout << avgTick[i] << " " << avgTick[j] << "  " << std::abs(avgTick[i]-avgTick[j]) << std::endl;
        std::cout << summedADCInt[i] << " " << summedADCInt[j] << std::endl;
      }
    }

  }


  //MarryNearClusters(AdjHitVec);

  return;
}

int SolarNuAna::ViewMatchCluster(std::vector< std::vector<recob::Hit> > testcs,
                               double t, double z,
                               double a0, int xsign,
                               double &curY, double &curZ,int minN)
{
  int idx = -1;

  //for (std::vector<recob::Hit> testc : Clusterss[1]){
  double mdzdy = 0;
  double minA0dist = 1e5;
  for (int tIdx = 0; tIdx < int(testcs.size()); tIdx++){
    std::vector<recob::Hit> testc = testcs[tIdx];
    if (int(testc.size())<minN) continue;
    if (testc.size() == 0) continue;
    //double mCharge = 0;
    double mClustT = 0;
    double mClustY = 0;
    double mClustZ = 0;
    //double mPur = 0;
    int hiIdx = 0;
    /*
    for (int hIdx = 0; hIdx < int(testc.size()); hIdx++){
      mCharge += testc[hIdx].Integral();
      mClustT += testc[hIdx].Integral() * testc[hIdx].PeakTime();

      const std::vector<int> trks = bt_serv->HitToTrackIds(clockData,testc[hIdx]);
      if (trks.size()==0) continue;
      PType ThisPType = WhichParType( trks[0] );
      if (ThisPType == 1)
        mPur += testc[hIdx].Integral();

      if (testc[hIdx].Integral() > testc[hiIdx].Integral())
        hiIdx = hIdx;
    }
    */
    //mClustT /= mCharge;
    //if (abs(mClustT - t)>100) continue;

    for (int hIdx = 0; hIdx <= int(testc.size()); hIdx++){
      recob::Hit hit = testc[hIdx];
      mClustT = hit.PeakTime();

      int wireID = -1;
      std::vector<geo::WireID> wires =
        geo->GeometryCore::ChannelToWire(hit.Channel());
      for (int i = 0; i < int(wires.size()); i++){
        const geo::WireGeo* wire = geo->GeometryCore::WirePtr(wires[i]);

        double strtXYZ[3];
        wire->GetStart(strtXYZ);

        double stopXYZ[3];
        wire->GetEnd(stopXYZ);

        if (strtXYZ[0] * xsign > 0 &&
            ( (strtXYZ[2] < z && stopXYZ[2] > z) ||
              (strtXYZ[2] > z && stopXYZ[2] < z)) )
          wireID = i;
      }
      if (wireID == -1) continue;
      const geo::WireGeo* wire=geo->GeometryCore::WirePtr(wires[wireID]);

      double dyds = wire->Direction()[1];
      double dzds = wire->Direction()[2];
      double hXYZ[3];
      wire->GetCenter(hXYZ);

      mClustY = hXYZ[1];
      mClustZ = hXYZ[2];
      //mClustY += hit.Integral() * hXYZ[1];
      //mClustZ += hit.Integral() * hXYZ[2];
      mdzdy = dzds/dyds; // same for all hits in current view

      //mClustT /= mCharge;
      if (t - mClustT < a0+15 && t - mClustT > a0-15 && abs(t-mClustT-a0)<minA0dist){
        minA0dist = abs(t-mClustT-a0);
        idx = tIdx;
        curY = mClustY + (z - mClustZ)/mdzdy;
        std::cout << mClustY << " " << z << " " << mClustZ << " " << mdzdy << " " << hiIdx << "  " << t-mClustT-a0 << std::endl;
      }
      if (mClustY < -1e20 || mClustZ < -1e20) continue;
    }
  }



  return idx;
}



//void MarryNearClusters(AdjHitVec);


//......................................................
DEFINE_ART_MODULE(SolarNuAna)
