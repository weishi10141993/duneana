#include <functional>   // std::greater

//ROOT includes
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TMath.h"
#include "TGeoMatrix.h" // TGeoHMatrix

// LArSoft libraries

//LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SupernovaTruth.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"

//ART includes
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
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

enum PType{kUnknown=0, kMarl, kAr39_Lar, kAr42_Lar, kK42_Lar, kKr85_Lar, kRn222_Lar, kPo218_Lar, kPb214_Lar, kBi214_Lar, kPb210_Lar, kK40_CPA, kAr42_CPA, kU238_CPA, kPo218_CPA, kPb214_CPA, kBi214_CPA, kPb210_CPA, kBi210_CPA, kCo60_APA, kU238_APA, kRn222_PDS, kNeut, kCavGam};
std::map<PType, std::string> PTypeString{
  {kUnknown  ,"Unknown"  },
  {kMarl     ,"Marl"     },
  {kAr39_Lar ,"Ar39_Lar" },
  {kAr42_Lar ,"Ar42_Lar" },
  {kK42_Lar  ,"K42_Lar"  },
  {kKr85_Lar ,"Kr85_Lar" },
  {kRn222_Lar,"Rn222_Lar"},
  {kPo218_Lar,"Po218_Lar"},
  {kPb214_Lar,"Pb214_Lar"},
  {kBi214_Lar,"Bi214_Lar"},
  {kPb210_Lar,"Pb210_Lar"},
  {kK40_CPA  ,"K40_CPA"  },
  {kAr42_CPA ,"Ar42_CPA" },
  {kU238_CPA ,"U238_CPA" },
  {kPo218_CPA,"Po218_CPA"},
  {kPb214_CPA,"Pb214_CPA"},
  {kBi214_CPA,"Bi214_CPA"},
  {kPb210_CPA,"Pb210_CPA"},
  {kBi210_CPA,"Bi210_CPA"},
  {kCo60_APA ,"Co60_APA" },
  {kU238_APA ,"U238_APA" },
  {kRn222_PDS,"Rn222_PDS"},
  {kNeut     ,"Neut"     },
  {kCavGam   ,"CavGam"   }
};


class SNAna : public art::EDAnalyzer {

public:
  explicit SNAna(fhicl::ParameterSet const & p);

  SNAna(SNAna const &) = delete;
  SNAna(SNAna &&) = delete;
  SNAna & operator = (SNAna const &) = delete;
  SNAna & operator = (SNAna &&) = delete;

  void analyze(art::Event const & evt) override;
  void reconfigure(fhicl::ParameterSet const & p);
  void beginJob() override;
  void endJob() override;

private:

  void ResetVariables();
  void FillMyMaps  ( std::map< int, simb::MCParticle> &MyMap,
                     art::FindManyP<simb::MCParticle> Assn,
                     art::Handle< std::vector<simb::MCTruth> > Hand,
                     std::map<int, int>* indexMap=nullptr);
  PType WhichParType( int TrID );
  bool  InMyMap     ( int TrID, std::map< int, simb::MCParticle> ParMap );
  void FillTruth(const art::FindManyP<simb::MCParticle> Assn,
                 const art::Handle<std::vector<simb::MCTruth>>& Hand,
                 const PType type) ;
  void SaveNeighbourADC(int channel,
                        art::Handle< std::vector<raw::RawDigit> >rawDigitsVecHandle,
                        std::set<int> badChannels,
                        recob::Hit const& hits);

  void SaveIDEs(art::Event const & evt);
  std::string fname;
  
  int firstCatch;
  int secondCatch;
  int thirdCatch;

  std::string fRawDigitLabel;
  std::string fHitLabel;
  std::string fCalDataModuleLabel;
  std::string fOpHitModuleLabel;

  std::string fGEANTLabel;
  std::string fMARLLabel      ; std::map< int, simb::MCParticle > MarlParts;
  std::string fAr39_LarLabel  ; std::map< int, simb::MCParticle > Ar39_LarParts;
  std::string fAr42_LarLabel  ; std::map< int, simb::MCParticle > Ar42_LarParts;
  std::string fK42_LarLabel   ; std::map< int, simb::MCParticle > K42_LarParts;
  std::string fKr85_LarLabel  ; std::map< int, simb::MCParticle > Kr85_LarParts;
  std::string fRn222_LarLabel ; std::map< int, simb::MCParticle > Rn222_LarParts;
  std::string fPo218_LarLabel ; std::map< int, simb::MCParticle > Po218_LarParts;
  std::string fPb214_LarLabel ; std::map< int, simb::MCParticle > Pb214_LarParts;
  std::string fBi214_LarLabel ; std::map< int, simb::MCParticle > Bi214_LarParts;
  std::string fPb210_LarLabel ; std::map< int, simb::MCParticle > Pb210_LarParts;

  std::string fK40_CPALabel   ; std::map< int, simb::MCParticle > K40_CPAParts;
  std::string fAr42_CPALabel  ; std::map< int, simb::MCParticle > Ar42_CPAParts;
  std::string fU238_CPALabel  ; std::map< int, simb::MCParticle > U238_CPAParts;
  std::string fPo218_CPALabel ; std::map< int, simb::MCParticle > Po218_CPAParts;
  std::string fPb214_CPALabel ; std::map< int, simb::MCParticle > Pb214_CPAParts;
  std::string fBi214_CPALabel ; std::map< int, simb::MCParticle > Bi214_CPAParts;
  std::string fPb210_CPALabel ; std::map< int, simb::MCParticle > Pb210_CPAParts;
  std::string fBi210_CPALabel ; std::map< int, simb::MCParticle > Bi210_CPAParts;

  std::string fCo60_APALabel  ; std::map< int, simb::MCParticle > Co60_APAParts;
  std::string fU238_APALabel  ; std::map< int, simb::MCParticle > U238_APAParts;

  std::string fRn222_PDSLabel ; std::map< int, simb::MCParticle > Rn222_PDSParts;

  std::string fNeutLabel      ; std::map< int, simb::MCParticle > NeutParts;
  std::string fCavGamLabel    ; std::map< int, simb::MCParticle > CavGamParts;

  std::map<int, const simb::MCParticle*> truthmap;
  // Which MARLEY interaction (if any) caused this true track ID?
  std::map<int, int> trkIDToMarleyIndex;

  // Mapping from track ID to particle type, for use in WhichParType()
  std::map<int, PType> trkIDToPType;

  bool fSaveNeighbourADCs;
  bool fSaveIDEs;
  bool fSaveTruth;
  bool fSaveTPC;
  bool fSavePDS;

  TTree* fSNAnaTree;

  int Run;
  int SubRun;
  int Event;

  int NTotHit   ;
  int NColHit   ;
  int NIndHit   ;
  int NHitNoBT  ;

  std::vector<int>                  Hit_View                 ;
  std::vector<int>                  Hit_Size                 ;
  std::vector<int>                  Hit_TPC                  ;
  std::vector<int>                  Hit_Chan                 ;
  std::vector<double>               Hit_X_start              ;
  std::vector<double>               Hit_Y_start              ;
  std::vector<double>               Hit_Z_start              ;
  std::vector<double>               Hit_X_end                ;
  std::vector<double>               Hit_Y_end                ;
  std::vector<double>               Hit_Z_end                ;
  std::vector<float>                Hit_Time                 ;
  std::vector<float>                Hit_RMS                  ;
  std::vector<float>                Hit_SADC                 ;
  std::vector<float>                Hit_Int                  ;
  std::vector<float>                Hit_Peak                 ;
  std::vector<int>                  Hit_True_GenType         ;
  std::vector<int>                  Hit_True_MainTrID        ;
  std::vector<int>                  Hit_True_TrackID         ;
  std::vector<float>                Hit_True_EvEnergy        ;
  std::vector<int>                  Hit_True_MarleyIndex     ;
  std::vector<float>                Hit_True_X               ;
  std::vector<float>                Hit_True_Y               ;
  std::vector<float>                Hit_True_Z               ;
  std::vector<float>                Hit_True_Energy          ;
  std::vector<float>                Hit_True_nElec           ;
  std::vector<int>                  Hit_True_nIDEs           ;
  std::vector<int>                  Hit_AdjM5SADC            ;
  std::vector<int>                  Hit_AdjM2SADC            ;
  std::vector<int>                  Hit_AdjM1SADC            ;
  std::vector<int>                  Hit_AdjP1SADC            ;
  std::vector<int>                  Hit_AdjP2SADC            ;
  std::vector<int>                  Hit_AdjP5SADC            ;
  std::vector<int>                  Hit_AdjM5Chan            ;
  std::vector<int>                  Hit_AdjM2Chan            ;
  std::vector<int>                  Hit_AdjM1Chan            ;
  std::vector<int>                  Hit_AdjP1Chan            ;
  std::vector<int>                  Hit_AdjP2Chan            ;
  std::vector<int>                  Hit_AdjP5Chan            ;

  std::vector<int>                  PDS_OpHit_OpChannel      ;
  std::vector<double>               PDS_OpHit_X              ;
  std::vector<double>               PDS_OpHit_Y              ;
  std::vector<double>               PDS_OpHit_Z              ;
  std::vector<double>               PDS_OpHit_PeakTimeAbs    ;
  std::vector<double>               PDS_OpHit_PeakTime       ;
  std::vector<unsigned short>       PDS_OpHit_Frame          ;
  std::vector<double>               PDS_OpHit_Width          ;
  std::vector<double>               PDS_OpHit_Area           ;
  std::vector<double>               PDS_OpHit_Amplitude      ;
  std::vector<double>               PDS_OpHit_PE             ;
  std::vector<double>               PDS_OpHit_FastToTotal    ;
  std::vector<int>                  PDS_OpHit_True_GenType   ;
  std::vector<int>                  PDS_OpHit_True_Index     ;
  std::vector<double>               PDS_OpHit_True_Energy    ;
  std::vector<int>                  PDS_OpHit_True_TrackID   ;
  std::vector<int>                  PDS_OpHit_True_GenTypeAll;
  std::vector<double>               PDS_OpHit_True_EnergyAll ;
  std::vector<int>                  PDS_OpHit_True_TrackIDAll;
  std::vector<int>                  PDS_OpHit_True_IndexAll  ;

  std::vector<int>                  True_VertexChan          ;
  std::vector<int>                  True_Nu_Type             ;
  std::vector<int>                  True_Nu_Lep_Type         ;
  std::vector<int>                  True_Mode                ;
  std::vector<int>                  True_CCNC                ;
  std::vector<int>                  True_HitNucleon          ;
  std::vector<int>                  True_Target              ;
  std::vector<int>                  True_MarlSample          ;
  std::vector<float>                True_MarlTime            ;
  std::vector<float>                True_MarlWeight          ;
  std::vector<float>                True_ENu                 ;
  std::vector<float>                True_ENu_Lep             ;
  std::vector<float>                True_VertX               ;
  std::vector<float>                True_VertY               ;
  std::vector<float>                True_VertZ               ;
  std::vector<float>                True_VertexT             ;
  std::vector<float>                True_Px                  ;
  std::vector<float>                True_Py                  ;
  std::vector<float>                True_Pz                  ;
  std::vector<float>                True_Dirx                ;
  std::vector<float>                True_Diry                ;
  std::vector<float>                True_Dirz                ;
  std::vector<float>                True_Time                ;

  std::vector<int>                  True_Bck_Mode            ;   
  std::vector<int>                  True_Bck_PDG             ;
  std::vector<int>                  True_Bck_ID              ;
  std::vector<std::string>          True_Bck_Process         ; // str because why not
  std::vector<std::string>          True_Bck_EndProcess      ;
  std::vector<int>                  True_Bck_Mother          ;
  std::vector<double>               True_Bck_P               ;
  std::vector<double>               True_Bck_VertX           ;
  std::vector<double>               True_Bck_VertY           ;
  std::vector<double>               True_Bck_VertZ           ;
  std::vector<double>               True_Bck_Time            ;
  std::vector<double>               True_Bck_Energy          ;
  std::vector<double>               True_Bck_EndX            ;
  std::vector<double>               True_Bck_EndY            ;
  std::vector<double>               True_Bck_EndZ            ;
  std::vector<double>               True_Bck_EndT            ;
  std::vector<double>               True_Bck_EndE            ;   

  int   NTotIDEs;
  std::vector<int>                  IDEChannel               ;
  std::vector<int>                  IDEStartTime             ;
  std::vector<int>                  IDEEndTime               ;
  std::vector<float>                IDEEnergy                ;
  std::vector<float>                IDEElectrons             ;
  std::vector<int>                  IDEParticle              ;

  int   TotGen_Marl     ;
  int   TotGen_Ar39_Lar ;
  int   TotGen_Ar42_Lar ;
  int   TotGen_K42_Lar  ;
  int   TotGen_Kr85_Lar ;
  int   TotGen_Rn222_Lar;
  int   TotGen_Po218_Lar;
  int   TotGen_Pb214_Lar;
  int   TotGen_Bi214_Lar;
  int   TotGen_Pb210_Lar;
  int   TotGen_K40_CPA  ;
  int   TotGen_Ar42_CPA ;
  int   TotGen_U238_CPA ;
  int   TotGen_Po218_CPA;
  int   TotGen_Pb214_CPA;
  int   TotGen_Bi214_CPA;
  int   TotGen_Pb210_CPA;
  int   TotGen_Bi210_CPA;
  int   TotGen_Co60_APA ;
  int   TotGen_U238_APA ;
  int   TotGen_Rn222_PDS;
  int   TotGen_Neut     ;
  int   TotGen_CavGam   ;



  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::PhotonBackTrackerService> pbt_serv;

};

SNAna::SNAna(fhicl::ParameterSet const & p):EDAnalyzer(p),
                                            fname("SNAna_module")
{
  this->reconfigure(p);
}


void SNAna::reconfigure(fhicl::ParameterSet const & p)
{
  
  fRawDigitLabel      = p.get<std::string>("RawDigitLabel"     );
  fHitLabel           = p.get<std::string>("HitLabel"          );
  fCalDataModuleLabel = p.get<std::string>("CalDataModuleLabel");
  fOpHitModuleLabel   = p.get<std::string>("OpHitModuleLabel"  );

  fGEANTLabel     = p.get<std::string> ("GEANT4Label"      );
  fMARLLabel      = p.get<std::string> ("MARLEYLabel"      );
  fAr39_LarLabel  = p.get<std::string> ("Argon39_Lar_Label");
  fAr42_LarLabel  = p.get<std::string> ("Argon42_Lar_Label");
  fK42_LarLabel   = p.get<std::string> ("Potas42_Lar_Label");
  fKr85_LarLabel  = p.get<std::string> ("Krypton_Lar_Label");
  fRn222_LarLabel = p.get<std::string> ("Rado222_Lar_Label");
  fPo218_LarLabel = p.get<std::string> ("Po218_Lar_Label");
  fPb214_LarLabel = p.get<std::string> ("Pb214_Lar_Label");
  fBi214_LarLabel = p.get<std::string> ("Bi214_Lar_Label");
  fPb210_LarLabel = p.get<std::string> ("Pb210_Lar_Label");
  fK40_CPALabel   = p.get<std::string> ("Potas40_CPA_Label");
  fAr42_CPALabel  = p.get<std::string> ("Argon42_CPA_Label");
  fU238_CPALabel  = p.get<std::string> ("U238Cha_CPA_Label");
  fPo218_CPALabel = p.get<std::string> ("Po218_CPA_Label");
  fPb214_CPALabel = p.get<std::string> ("Pb214_CPA_Label");
  fBi214_CPALabel = p.get<std::string> ("Bi214_CPA_Label");
  fPb210_CPALabel = p.get<std::string> ("Pb210_CPA_Label");
  fBi210_CPALabel = p.get<std::string> ("Bi210_CPA_Label");
  fCo60_APALabel  = p.get<std::string> ("Cobal60_APA_Label");
  fU238_APALabel  = p.get<std::string> ("U238Cha_APA_Label");
  fRn222_PDSLabel = p.get<std::string> ("Rado222_PDS_Label");
  fNeutLabel      = p.get<std::string> ("NeutronLabel"     );
  fCavGamLabel    = p.get<std::string> ("CavernGammasLabel");


  fSaveNeighbourADCs = p.get<bool> ("SaveNeighbourADCs",0);
  fSaveTruth = p.get<bool>("SaveTruth",0);
  fSaveIDEs  = p.get<bool>("SaveIDEs",0);
  fSaveTPC   = p.get<bool>("SaveTPC",1);
  fSavePDS   = p.get<bool>("SavePDS",1);

  mf::LogInfo(fname) << "Reconfigured " << this->processName() << " with "
                     << " SaveNeighbourADCs: " << std::boolalpha << fSaveNeighbourADCs
                     << " SaveTruth: " << std::boolalpha << fSaveTruth
                     << " SaveIDEs: " << std::boolalpha << fSaveIDEs << std::endl;
}


void SNAna::ResetVariables()
{
  trkIDToPType.clear();

  MarlParts     .clear(); Ar39_LarParts .clear(); Ar42_LarParts .clear(); K42_LarParts  .clear();
  Kr85_LarParts .clear(); Rn222_LarParts.clear(); Po218_LarParts.clear(); Pb214_LarParts.clear();
  Bi214_LarParts.clear(); Pb210_LarParts.clear(); K40_CPAParts  .clear(); Ar42_CPAParts .clear();
  U238_CPAParts .clear(); Po218_CPAParts.clear(); Pb214_CPAParts.clear(); Bi214_CPAParts.clear();
  Pb210_CPAParts.clear(); Bi210_CPAParts.clear(); Co60_APAParts .clear(); U238_APAParts .clear();
  Rn222_PDSParts.clear(); NeutParts     .clear(); CavGamParts   .clear();

  Run = SubRun = Event = -1;

  TotGen_Marl      = TotGen_Ar39_Lar  = TotGen_Ar42_Lar  = TotGen_K42_Lar   = 0;
  TotGen_Kr85_Lar  = TotGen_Rn222_Lar = TotGen_Po218_Lar = TotGen_Pb214_Lar = 0;
  TotGen_Bi214_Lar = TotGen_Pb210_Lar = TotGen_K40_CPA   = TotGen_Ar42_CPA  = 0;
  TotGen_U238_CPA  = TotGen_Po218_CPA = TotGen_Pb214_CPA = TotGen_Bi214_CPA = 0;
  TotGen_Pb210_CPA = TotGen_Bi210_CPA = TotGen_Co60_APA  = TotGen_U238_APA  = 0;
  TotGen_Rn222_PDS = TotGen_Neut      = TotGen_CavGam    = 0;

  NTotHit    = 0;
  NColHit    = 0;
  NIndHit    = 0;
  NHitNoBT   = 0;

  Hit_View                 .clear();
  Hit_Size                 .clear();
  Hit_TPC                  .clear();
  Hit_Chan                 .clear();
  Hit_X_start              .clear();
  Hit_Y_start              .clear();
  Hit_Z_start              .clear();
  Hit_X_end                .clear();
  Hit_Y_end                .clear();
  Hit_Z_end                .clear();
  Hit_Time                 .clear();
  Hit_RMS                  .clear();
  Hit_SADC                 .clear();
  Hit_Int                  .clear();
  Hit_Peak                 .clear();
  Hit_True_GenType         .clear();
  Hit_True_MainTrID        .clear();
  Hit_True_TrackID         .clear();
  Hit_True_EvEnergy        .clear();
  Hit_True_MarleyIndex     .clear();
  Hit_True_X               .clear();
  Hit_True_Y               .clear();
  Hit_True_Z               .clear();
  Hit_True_Energy          .clear();
  Hit_True_nElec           .clear();
  Hit_True_nIDEs           .clear();

  Hit_AdjM5SADC            .clear();
  Hit_AdjM2SADC            .clear();
  Hit_AdjM1SADC            .clear();
  Hit_AdjP1SADC            .clear();
  Hit_AdjP2SADC            .clear();
  Hit_AdjP5SADC            .clear();
  Hit_AdjM5Chan            .clear();
  Hit_AdjM2Chan            .clear();
  Hit_AdjM1Chan            .clear();
  Hit_AdjP1Chan            .clear();
  Hit_AdjP2Chan            .clear();
  Hit_AdjP5Chan            .clear();

  PDS_OpHit_OpChannel      .clear();
  PDS_OpHit_X              .clear();
  PDS_OpHit_Y              .clear();
  PDS_OpHit_Z              .clear();
  PDS_OpHit_PeakTimeAbs    .clear();
  PDS_OpHit_PeakTime       .clear();
  PDS_OpHit_Frame          .clear();
  PDS_OpHit_Width          .clear();
  PDS_OpHit_Area           .clear();
  PDS_OpHit_Amplitude      .clear();
  PDS_OpHit_PE             .clear();
  PDS_OpHit_FastToTotal    .clear();
  PDS_OpHit_True_GenType   .clear();
  PDS_OpHit_True_Index     .clear();
  PDS_OpHit_True_Energy    .clear();
  PDS_OpHit_True_TrackID   .clear();
  PDS_OpHit_True_GenTypeAll.clear();
  PDS_OpHit_True_EnergyAll .clear();
  PDS_OpHit_True_TrackIDAll.clear();
  PDS_OpHit_True_IndexAll  .clear();

  True_VertexChan          .clear();
  True_Nu_Type             .clear();
  True_Nu_Lep_Type         .clear();
  True_Mode                .clear();
  True_CCNC                .clear();
  True_HitNucleon          .clear();
  True_Target              .clear();
  True_MarlSample          .clear();
  True_MarlTime            .clear();
  True_MarlWeight          .clear();
  True_ENu                 .clear();
  True_ENu_Lep             .clear();
  True_VertX               .clear();
  True_VertY               .clear();
  True_VertZ               .clear();
  True_VertexT             .clear();
  True_Px                  .clear();
  True_Py                  .clear();
  True_Pz                  .clear();
  True_Dirx                .clear();
  True_Diry                .clear();
  True_Dirz                .clear();
  True_Time                .clear();

  True_Bck_Mode            .clear();
  True_Bck_PDG             .clear();
  True_Bck_ID              .clear();
  True_Bck_Process         .clear();
  True_Bck_EndProcess      .clear();
  True_Bck_Mother          .clear();
  True_Bck_P               .clear();
  True_Bck_VertX           .clear();
  True_Bck_VertY           .clear();
  True_Bck_VertZ           .clear();
  True_Bck_Time            .clear();
  True_Bck_Energy          .clear();
  True_Bck_EndX            .clear();
  True_Bck_EndY            .clear();
  True_Bck_EndZ            .clear();
  True_Bck_EndT            .clear();
  True_Bck_EndE            .clear();

  // IDEs
  NTotIDEs=0;
  IDEChannel               .clear();
  IDEStartTime             .clear();
  IDEEndTime               .clear();
  IDEEnergy                .clear();
  IDEElectrons             .clear();
  IDEParticle              .clear();

}


void SNAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;


  fSNAnaTree = tfs->make<TTree>("SNSimTree","SN simulation analysis tree");

  fSNAnaTree->Branch("Run"       , &Run       , "Run/I"       );
  fSNAnaTree->Branch("SubRun"    , &SubRun    , "SubRun/I"    );
  fSNAnaTree->Branch("Event"     , &Event     , "Event/I"     );
  fSNAnaTree->Branch("NTotHit"   , &NTotHit   , "NTotHits/I"  );
  fSNAnaTree->Branch("NColHit"   , &NColHit   , "NColHits/I"  );
  fSNAnaTree->Branch("NIndHit"   , &NIndHit   , "NIndHits/I"  );
  fSNAnaTree->Branch("NHitNoBT"  , &NHitNoBT  , "NHitNoBT/I"  );

  if (fSaveTPC) {
    fSNAnaTree->Branch("Hit_View"                 , &Hit_View                 );
    fSNAnaTree->Branch("Hit_Size"                 , &Hit_Size                 );
    fSNAnaTree->Branch("Hit_TPC"                  , &Hit_TPC                  );
    fSNAnaTree->Branch("Hit_Chan"                 , &Hit_Chan                 );
    fSNAnaTree->Branch("Hit_X_start"              , &Hit_X_start              );
    fSNAnaTree->Branch("Hit_Y_start"              , &Hit_Y_start              );
    fSNAnaTree->Branch("Hit_Z_start"              , &Hit_Z_start              );
    fSNAnaTree->Branch("Hit_X_end"                , &Hit_X_end                );
    fSNAnaTree->Branch("Hit_Y_end"                , &Hit_Y_end                );
    fSNAnaTree->Branch("Hit_Z_end"                , &Hit_Z_end                );
    fSNAnaTree->Branch("Hit_Time"                 , &Hit_Time                 );
    fSNAnaTree->Branch("Hit_RMS"                  , &Hit_RMS                  );
    fSNAnaTree->Branch("Hit_SADC"                 , &Hit_SADC                 );
    fSNAnaTree->Branch("Hit_Int"                  , &Hit_Int                  );
    fSNAnaTree->Branch("Hit_Peak"                 , &Hit_Peak                 );
    fSNAnaTree->Branch("Hit_True_GenType"         , &Hit_True_GenType         );
    fSNAnaTree->Branch("Hit_True_MainTrID"        , &Hit_True_MainTrID        );
    fSNAnaTree->Branch("Hit_True_TrackID"         , &Hit_True_TrackID         );
    fSNAnaTree->Branch("Hit_True_EvEnergy"        , &Hit_True_EvEnergy        );
    fSNAnaTree->Branch("Hit_True_MarleyIndex"     , &Hit_True_MarleyIndex     );
    fSNAnaTree->Branch("Hit_True_X"               , &Hit_True_X               );
    fSNAnaTree->Branch("Hit_True_Y"               , &Hit_True_Y               );
    fSNAnaTree->Branch("Hit_True_Z"               , &Hit_True_Z               );
    fSNAnaTree->Branch("Hit_True_Energy"          , &Hit_True_Energy          );
    fSNAnaTree->Branch("Hit_True_nElec"           , &Hit_True_nElec           );
    fSNAnaTree->Branch("Hit_True_nIDEs"           , &Hit_True_nIDEs           );
  }

  if (fSaveNeighbourADCs) {
    fSNAnaTree->Branch("Hit_AdjM5SADC"            , &Hit_AdjM5SADC            );
    fSNAnaTree->Branch("Hit_AdjM2SADC"            , &Hit_AdjM2SADC            );
    fSNAnaTree->Branch("Hit_AdjM1SADC"            , &Hit_AdjM1SADC            );
    fSNAnaTree->Branch("Hit_AdjP1SADC"            , &Hit_AdjP1SADC            );
    fSNAnaTree->Branch("Hit_AdjP2SADC"            , &Hit_AdjP2SADC            );
    fSNAnaTree->Branch("Hit_AdjP5SADC"            , &Hit_AdjP5SADC            );
    fSNAnaTree->Branch("Hit_AdjM5Chan"            , &Hit_AdjM5Chan            );
    fSNAnaTree->Branch("Hit_AdjM2Chan"            , &Hit_AdjM2Chan            );
    fSNAnaTree->Branch("Hit_AdjM1Chan"            , &Hit_AdjM1Chan            );
    fSNAnaTree->Branch("Hit_AdjP1Chan"            , &Hit_AdjP1Chan            );
    fSNAnaTree->Branch("Hit_AdjP2Chan"            , &Hit_AdjP2Chan            );
    fSNAnaTree->Branch("Hit_AdjP5Chan"            , &Hit_AdjP5Chan            );
  }

  if (fSavePDS) {
    fSNAnaTree->Branch("PDS_OpHit_OpChannel"      , &PDS_OpHit_OpChannel      );
    fSNAnaTree->Branch("PDS_OpHit_X"              , &PDS_OpHit_X              );
    fSNAnaTree->Branch("PDS_OpHit_Y"              , &PDS_OpHit_Y              );
    fSNAnaTree->Branch("PDS_OpHit_Z"              , &PDS_OpHit_Z              );
    fSNAnaTree->Branch("PDS_OpHit_PeakTimeAbs"    , &PDS_OpHit_PeakTimeAbs    );
    fSNAnaTree->Branch("PDS_OpHit_PeakTime"       , &PDS_OpHit_PeakTime       );
    fSNAnaTree->Branch("PDS_OpHit_Frame"          , &PDS_OpHit_Frame          );
    fSNAnaTree->Branch("PDS_OpHit_Width"          , &PDS_OpHit_Width          );
    fSNAnaTree->Branch("PDS_OpHit_Area"           , &PDS_OpHit_Area           );
    fSNAnaTree->Branch("PDS_OpHit_Amplitude"      , &PDS_OpHit_Amplitude      );
    fSNAnaTree->Branch("PDS_OpHit_PE"             , &PDS_OpHit_PE             );
    fSNAnaTree->Branch("PDS_OpHit_FastToTotal"    , &PDS_OpHit_FastToTotal    );
    fSNAnaTree->Branch("PDS_OpHit_True_GenType"   , &PDS_OpHit_True_GenType   );
    fSNAnaTree->Branch("PDS_OpHit_True_Energy"    , &PDS_OpHit_True_Energy    );
    fSNAnaTree->Branch("PDS_OpHit_True_TrackID"   , &PDS_OpHit_True_TrackID   );
    fSNAnaTree->Branch("PDS_OpHit_True_GenTypeAll", &PDS_OpHit_True_GenTypeAll);
    fSNAnaTree->Branch("PDS_OpHit_True_EnergyAll" , &PDS_OpHit_True_EnergyAll );
    fSNAnaTree->Branch("PDS_OpHit_True_TrackIDAll", &PDS_OpHit_True_TrackIDAll);
    fSNAnaTree->Branch("PDS_OpHit_True_IndexAll"  , &PDS_OpHit_True_IndexAll  );
  }

  fSNAnaTree->Branch("True_VertexChan"          , &True_VertexChan          );
  fSNAnaTree->Branch("True_Nu_Type"             , &True_Nu_Type             );
  fSNAnaTree->Branch("True_Nu_Lep_Type"         , &True_Nu_Lep_Type         );
  fSNAnaTree->Branch("True_Mode"                , &True_Mode                );
  fSNAnaTree->Branch("True_CCNC"                , &True_CCNC                );
  fSNAnaTree->Branch("True_HitNucleon"          , &True_HitNucleon          );
  fSNAnaTree->Branch("True_Target"              , &True_Target              );
  fSNAnaTree->Branch("True_MarlSample"          , &True_MarlSample          );
  fSNAnaTree->Branch("True_MarlTime"            , &True_MarlTime            );
  fSNAnaTree->Branch("True_MarlWeight"          , &True_MarlWeight          );
  fSNAnaTree->Branch("True_ENu"                 , &True_ENu                 );
  fSNAnaTree->Branch("True_ENu_Lep"             , &True_ENu_Lep             );
  fSNAnaTree->Branch("True_VertX"               , &True_VertX               );
  fSNAnaTree->Branch("True_VertY"               , &True_VertY               );
  fSNAnaTree->Branch("True_VertZ"               , &True_VertZ               );
  fSNAnaTree->Branch("True_VertexT"             , &True_VertexT             );
  fSNAnaTree->Branch("True_Px"                  , &True_Px                  );
  fSNAnaTree->Branch("True_Py"                  , &True_Py                  );
  fSNAnaTree->Branch("True_Pz"                  , &True_Pz                  );
  fSNAnaTree->Branch("True_Dirx"                , &True_Dirx                );
  fSNAnaTree->Branch("True_Diry"                , &True_Diry                );
  fSNAnaTree->Branch("True_Dirz"                , &True_Dirz                );
  fSNAnaTree->Branch("True_Time"                , &True_Time                );


  fSNAnaTree->Branch("True_Bck_Mode"            , &True_Bck_Mode            );
  fSNAnaTree->Branch("True_Bck_PDG"             , &True_Bck_PDG             );
  fSNAnaTree->Branch("True_Bck_ID"              , &True_Bck_ID              );
  fSNAnaTree->Branch("True_Bck_Process"         , &True_Bck_Process         );
  fSNAnaTree->Branch("True_Bck_EndProcess"      , &True_Bck_EndProcess      );
  fSNAnaTree->Branch("True_Bck_Mother"          , &True_Bck_Mother          );
  fSNAnaTree->Branch("True_Bck_P"               , &True_Bck_P               );
  fSNAnaTree->Branch("True_Bck_VertX"           , &True_Bck_VertX           );
  fSNAnaTree->Branch("True_Bck_VertY"           , &True_Bck_VertY           );
  fSNAnaTree->Branch("True_Bck_VertZ"           , &True_Bck_VertZ           );
  fSNAnaTree->Branch("True_Bck_Time"            , &True_Bck_Time            );
  fSNAnaTree->Branch("True_Bck_Energy"          , &True_Bck_Energy          );
  fSNAnaTree->Branch("True_Bck_EndX"            , &True_Bck_EndX            );
  fSNAnaTree->Branch("True_Bck_EndY"            , &True_Bck_EndY            );
  fSNAnaTree->Branch("True_Bck_EndZ"            , &True_Bck_EndZ            );
  fSNAnaTree->Branch("True_Bck_EndT"            , &True_Bck_EndT            );
  fSNAnaTree->Branch("True_Bck_EndE"            , &True_Bck_EndE            );  

  // IDEs
  if(fSaveIDEs) {
    fSNAnaTree->Branch("NTotIDEs"                 , &NTotIDEs  , "NTotIDEs/I" );
    fSNAnaTree->Branch("IDEChannel"               , &IDEChannel               );
    fSNAnaTree->Branch("IDEStartTime"             , &IDEStartTime             );
    fSNAnaTree->Branch("IDEEndTime"               , &IDEEndTime               );
    fSNAnaTree->Branch("IDEEnergy"                , &IDEEnergy                );
    fSNAnaTree->Branch("IDEElectrons"             , &IDEElectrons             );
    fSNAnaTree->Branch("IDEParticle"              , &IDEParticle              );
  }

  fSNAnaTree->Branch("TotGen_Marl"     , &TotGen_Marl     , "TotGen_Marl/I"     );
  fSNAnaTree->Branch("TotGen_Ar39_Lar" , &TotGen_Ar39_Lar , "TotGen_Ar39_Lar/I" );
  fSNAnaTree->Branch("TotGen_Ar42_Lar" , &TotGen_Ar42_Lar , "TotGen_Ar42_Lar/I" );
  fSNAnaTree->Branch("TotGen_K42_Lar"  , &TotGen_K42_Lar  , "TotGen_K42_Lar/I"  );
  fSNAnaTree->Branch("TotGen_Kr85_Lar" , &TotGen_Kr85_Lar , "TotGen_Kr85_Lar/I" );
  fSNAnaTree->Branch("TotGen_Rn222_Lar", &TotGen_Rn222_Lar, "TotGen_Rn222_Lar/I");
  fSNAnaTree->Branch("TotGen_Po218_Lar", &TotGen_Po218_Lar, "TotGen_Po218_Lar/I");
  fSNAnaTree->Branch("TotGen_Pb214_Lar", &TotGen_Pb214_Lar, "TotGen_Pb214_Lar/I");
  fSNAnaTree->Branch("TotGen_Bi214_Lar", &TotGen_Bi214_Lar, "TotGen_Bi214_Lar/I");
  fSNAnaTree->Branch("TotGen_Pb210_Lar", &TotGen_Pb210_Lar, "TotGen_Pb210_Lar/I");
  fSNAnaTree->Branch("TotGen_K40_CPA"  , &TotGen_K40_CPA  , "TotGen_K40_CPA/I"  );
  fSNAnaTree->Branch("TotGen_Ar42_CPA" , &TotGen_Ar42_CPA , "TotGen_Ar42_CPA/I" );
  fSNAnaTree->Branch("TotGen_U238_CPA" , &TotGen_U238_CPA , "TotGen_U238_CPA/I" );
  fSNAnaTree->Branch("TotGen_Po218_CPA", &TotGen_Po218_CPA, "TotGen_Po218_CPA/I");
  fSNAnaTree->Branch("TotGen_Pb214_CPA", &TotGen_Pb214_CPA, "TotGen_Pb214_CPA/I");
  fSNAnaTree->Branch("TotGen_Bi214_CPA", &TotGen_Bi214_CPA, "TotGen_Bi214_CPA/I");
  fSNAnaTree->Branch("TotGen_Pb210_CPA", &TotGen_Pb210_CPA, "TotGen_Pb210_CPA/I");
  fSNAnaTree->Branch("TotGen_Bi210_CPA", &TotGen_Bi210_CPA, "TotGen_Bi210_CPA/I");
  fSNAnaTree->Branch("TotGen_Co60_APA" , &TotGen_Co60_APA , "TotGen_Co60_APA/I" );
  fSNAnaTree->Branch("TotGen_U238_APA" , &TotGen_U238_APA , "TotGen_U238_APA/I" );
  fSNAnaTree->Branch("TotGen_Rn222_PDS", &TotGen_Rn222_PDS, "TotGen_Rn222_PDS/I");
  fSNAnaTree->Branch("TotGen_Neut"     , &TotGen_Neut     , "TotGen_Neut/I"     );
  fSNAnaTree->Branch("TotGen_CavGam"   , &TotGen_CavGam   , "TotGen_CavGam/I"   );

}


void SNAna::analyze(art::Event const & evt)
{
  ResetVariables();

  Run    = evt.run();
  SubRun = evt.subRun();
  Event  = evt.event();

  //GET INFORMATION ABOUT THE DETECTOR'S GEOMETRY.
  auto const* geo = lar::providerFrom<geo::Geometry>();

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

  //LIFT OUT THE MARLEY PARTICLES.
  //auto MarlTrue = evt.getValidHandle<std::vector<simb::MCTruth> >(fMARLLabel);
  auto MarlTrue = evt.getHandle< std::vector<simb::MCTruth> >(fMARLLabel);
  if (MarlTrue) {
    art::FindManyP<simb::MCParticle> MarlAssn(MarlTrue,evt,fGEANTLabel);
    FillMyMaps( MarlParts, MarlAssn, MarlTrue, &trkIDToMarleyIndex );
    TotGen_Marl = MarlParts.size();
    double Px_(0), Py_(0), Pz_(0), Pnorm(1);

    for(size_t i = 0; i < MarlTrue->size(); i++)
    {
      True_Nu_Type    .push_back(MarlTrue->at(i).GetNeutrino().Nu().PdgCode());
      True_Nu_Lep_Type.push_back(MarlTrue->at(i).GetNeutrino().Lepton().PdgCode());
      True_Mode       .push_back(MarlTrue->at(i).GetNeutrino().Mode());
      True_CCNC       .push_back(MarlTrue->at(i).GetNeutrino().CCNC());
      True_Target     .push_back(MarlTrue->at(i).GetNeutrino().Target());
      True_HitNucleon .push_back(MarlTrue->at(i).GetNeutrino().HitNuc());
      True_VertX      .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Vx());
      True_VertY      .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Vy());
      True_VertZ      .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Vz());
      True_ENu_Lep    .push_back(MarlTrue->at(i).GetNeutrino().Lepton().E());
      True_ENu        .push_back(MarlTrue->at(i).GetNeutrino().Nu().E());
      True_Px         .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Px());
      True_Py         .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Py());
      True_Pz         .push_back(MarlTrue->at(i).GetNeutrino().Lepton().Pz());
      Pnorm = std::sqrt(Px_*Px_+Py_*Py_+Pz_*Pz_);
      double Px = Px_/Pnorm;
      double Py = Py_/Pnorm;
      double Pz = Pz_/Pnorm;
      True_Dirx.push_back(Px);
      True_Diry.push_back(Py);
      True_Dirz.push_back(Pz);
      True_Time.push_back(MarlTrue->at(i).GetNeutrino().Lepton().T());

      try {
        art::FindManyP<sim::SupernovaTruth> SNTruth(MarlTrue, evt, fMARLLabel);
        for (size_t j = 0; j < SNTruth.at(i).size(); j++) {
          const sim::SupernovaTruth ThisTr = (*SNTruth.at(i).at(j));
          True_MarlTime  .push_back(ThisTr.SupernovaTime);
          True_MarlWeight.push_back(ThisTr.Weight);
          True_MarlSample.push_back(ThisTr.SamplingMode);
        }
      } catch (...) {
        mf::LogDebug(fname) << "Didnt find SN truth (few things wont available):\n"
                            << " - MarlTime\n"
                            << " - MarlWeight\n"
                            << " - MarlSample\n";
      }
    }

    for(size_t i=0; i<True_VertX.size(); ++i) {
      bool caught = false;
      double Vertex[3] = {True_VertX[i],
                          True_VertY[i],
                          True_VertZ[i]};
      geo::WireID WireID;
      geo::PlaneID Plane(geo->FindTPCAtPosition(Vertex),geo::kZ);
      try
      {
        WireID = geo->NearestWireID(Vertex, Plane);
      }
      catch(...)
      {
        caught = true;
      }
      if(caught==true)
      {
        True_VertexChan.push_back(-1);
      }
      else
      {
        True_VertexChan.push_back(geo->PlaneWireToChannel(WireID));
      }

      //CM/MICROSECOND.
      double drift_velocity = detProp.DriftVelocity(detProp.Efield(),detProp.Temperature());
      //CM/TICK
      drift_velocity = drift_velocity*0.5;
      True_VertexT.push_back(True_VertX.back()/drift_velocity);
    }
    if (fSaveTruth)
      FillTruth(MarlAssn, MarlTrue, kMarl);
  }

  auto Ar39_LarTrue = evt.getHandle< std::vector<simb::MCTruth> >(fAr39_LarLabel);
  if (Ar39_LarTrue) {
    art::FindManyP<simb::MCParticle> Ar39_LarAssn(Ar39_LarTrue,evt,fGEANTLabel);
    FillMyMaps( Ar39_LarParts, Ar39_LarAssn, Ar39_LarTrue );
    TotGen_Ar39_Lar = Ar39_LarParts.size();
    if(fSaveTruth) FillTruth(Ar39_LarAssn , Ar39_LarTrue , kAr39_Lar );

  }

  auto Ar42_LarTrue = evt.getHandle< std::vector<simb::MCTruth> >(fAr42_LarLabel);
  if(Ar42_LarTrue) {
    art::FindManyP<simb::MCParticle> Ar42_LarAssn(Ar42_LarTrue,evt,fGEANTLabel);
    FillMyMaps( Ar42_LarParts, Ar42_LarAssn, Ar42_LarTrue );
    TotGen_Ar42_Lar = Ar42_LarParts.size();
    if(fSaveTruth) FillTruth(Ar42_LarAssn , Ar42_LarTrue , kAr42_Lar );
  }

  auto K42_LarTrue = evt.getHandle< std::vector<simb::MCTruth> > (fK42_LarLabel);
  if(K42_LarTrue) {
    art::FindManyP<simb::MCParticle> K42_LarAssn(K42_LarTrue,evt,fGEANTLabel);
    FillMyMaps( K42_LarParts, K42_LarAssn, K42_LarTrue );
    TotGen_K42_Lar = K42_LarParts.size();
    if(fSaveTruth) FillTruth(K42_LarAssn , K42_LarTrue , kK42_Lar );
  }

  auto Kr85_LarTrue = evt.getHandle< std::vector<simb::MCTruth> > (fKr85_LarLabel);
  if(Kr85_LarTrue) {
    art::FindManyP<simb::MCParticle> Kr85_LarAssn(Kr85_LarTrue,evt,fGEANTLabel);
    FillMyMaps( Kr85_LarParts, Kr85_LarAssn, Kr85_LarTrue );
    TotGen_Kr85_Lar = Kr85_LarParts.size();
    if(fSaveTruth) FillTruth(Kr85_LarAssn , Kr85_LarTrue , kKr85_Lar );
  }

  auto Rn222_LarTrue = evt.getHandle< std::vector<simb::MCTruth> > (fRn222_LarLabel);
  if(Rn222_LarTrue) {
    art::FindManyP<simb::MCParticle> Rn222_LarAssn(Rn222_LarTrue,evt,fGEANTLabel);
    FillMyMaps( Rn222_LarParts, Rn222_LarAssn, Rn222_LarTrue );
    TotGen_Rn222_Lar = Rn222_LarParts.size();
    if(fSaveTruth) FillTruth(Rn222_LarAssn , Rn222_LarTrue , kRn222_Lar );
  }

  auto Po218_LarTrue = evt.getHandle< std::vector<simb::MCTruth> > (fPo218_LarLabel);
  if(Po218_LarTrue) {
    art::FindManyP<simb::MCParticle> Po218_LarAssn(Po218_LarTrue,evt,fGEANTLabel);
    FillMyMaps( Po218_LarParts, Po218_LarAssn, Po218_LarTrue );
    TotGen_Po218_Lar = Po218_LarParts.size();
    if(fSaveTruth) FillTruth(Po218_LarAssn , Po218_LarTrue , kPo218_Lar );
  }

  auto Pb214_LarTrue = evt.getHandle< std::vector<simb::MCTruth> > (fPb214_LarLabel);
  if(Pb214_LarTrue) {
    art::FindManyP<simb::MCParticle> Pb214_LarAssn(Pb214_LarTrue,evt,fGEANTLabel);
    FillMyMaps( Pb214_LarParts, Pb214_LarAssn, Pb214_LarTrue );
    TotGen_Pb214_Lar = Pb214_LarParts.size();
    if(fSaveTruth) FillTruth(Pb214_LarAssn , Pb214_LarTrue , kPb214_Lar );
  }

  auto Bi214_LarTrue = evt.getHandle< std::vector<simb::MCTruth> > (fBi214_LarLabel);
  if(Bi214_LarTrue) {
    art::FindManyP<simb::MCParticle> Bi214_LarAssn(Bi214_LarTrue,evt,fGEANTLabel);
    FillMyMaps( Bi214_LarParts, Bi214_LarAssn, Bi214_LarTrue );
    TotGen_Bi214_Lar = Bi214_LarParts.size();
    if(fSaveTruth) FillTruth(Bi214_LarAssn , Bi214_LarTrue , kBi214_Lar );
  }

  auto Pb210_LarTrue = evt.getHandle< std::vector<simb::MCTruth> > (fPb210_LarLabel);
  if(Pb210_LarTrue) {
    art::FindManyP<simb::MCParticle> Pb210_LarAssn(Pb210_LarTrue,evt,fGEANTLabel);
    FillMyMaps( Pb210_LarParts, Pb210_LarAssn, Pb210_LarTrue );
    TotGen_Pb210_Lar = Pb210_LarParts.size();
    if(fSaveTruth) FillTruth(Pb210_LarAssn , Pb210_LarTrue , kPb210_Lar );
  }

  auto K40_CPATrue = evt.getHandle< std::vector<simb::MCTruth> > (fK40_CPALabel);
  if(K40_CPATrue) {
    art::FindManyP<simb::MCParticle> K40_CPAAssn(K40_CPATrue,evt,fGEANTLabel);
    FillMyMaps( K40_CPAParts, K40_CPAAssn, K40_CPATrue );
    TotGen_K40_CPA = K40_CPAParts.size();
    if(fSaveTruth) FillTruth(K40_CPAAssn , K40_CPATrue , kK40_CPA );
  }

  auto Ar42_CPATrue = evt.getHandle< std::vector<simb::MCTruth> > (fAr42_CPALabel);
  if(Ar42_CPATrue) {
    art::FindManyP<simb::MCParticle> Ar42_CPAAssn(Ar42_CPATrue,evt,fGEANTLabel);
    FillMyMaps( Ar42_CPAParts, Ar42_CPAAssn, Ar42_CPATrue );
    TotGen_Ar42_CPA = Ar42_CPAParts.size();
    if(fSaveTruth) FillTruth(Ar42_CPAAssn , Ar42_CPATrue , kAr42_CPA );
  }

  auto U238_CPATrue = evt.getHandle< std::vector<simb::MCTruth> > (fU238_CPALabel);
  if(U238_CPATrue) {
    art::FindManyP<simb::MCParticle> U238_CPAAssn(U238_CPATrue,evt,fGEANTLabel);
    FillMyMaps( U238_CPAParts, U238_CPAAssn, U238_CPATrue );
    TotGen_U238_CPA = U238_CPAParts.size();
    if(fSaveTruth) FillTruth(U238_CPAAssn , U238_CPATrue , kU238_CPA );
  }

  auto Po218_CPATrue = evt.getHandle< std::vector<simb::MCTruth> > (fPo218_CPALabel);
  if(Po218_CPATrue) {
    art::FindManyP<simb::MCParticle> Po218_CPAAssn(Po218_CPATrue,evt,fGEANTLabel);
    FillMyMaps( Po218_CPAParts, Po218_CPAAssn, Po218_CPATrue );
    TotGen_Po218_CPA = Po218_CPAParts.size();
    if(fSaveTruth) FillTruth(Po218_CPAAssn , Po218_CPATrue , kPo218_CPA );
  }

  auto Pb214_CPATrue = evt.getHandle< std::vector<simb::MCTruth> > (fPb214_CPALabel);
  if(Pb214_CPATrue) {
    art::FindManyP<simb::MCParticle> Pb214_CPAAssn(Pb214_CPATrue,evt,fGEANTLabel);
    FillMyMaps( Pb214_CPAParts, Pb214_CPAAssn, Pb214_CPATrue );
    TotGen_Pb214_CPA = Pb214_CPAParts.size();
    if(fSaveTruth) FillTruth(Pb214_CPAAssn , Pb214_CPATrue , kPb214_CPA );
  }

  auto Bi214_CPATrue = evt.getHandle< std::vector<simb::MCTruth> > (fBi214_CPALabel);
  if(Bi214_CPATrue) {
    art::FindManyP<simb::MCParticle> Bi214_CPAAssn(Bi214_CPATrue,evt,fGEANTLabel);
    FillMyMaps( Bi214_CPAParts, Bi214_CPAAssn, Bi214_CPATrue );
    TotGen_Bi214_CPA = Bi214_CPAParts.size();
    if(fSaveTruth) FillTruth(Bi214_CPAAssn , Bi214_CPATrue , kBi214_CPA );
  }

  auto Pb210_CPATrue = evt.getHandle< std::vector<simb::MCTruth> > (fPb210_CPALabel);
  if(Pb210_CPATrue) {
    art::FindManyP<simb::MCParticle> Pb210_CPAAssn(Pb210_CPATrue,evt,fGEANTLabel);
    FillMyMaps( Pb210_CPAParts, Pb210_CPAAssn, Pb210_CPATrue );
    TotGen_Pb210_CPA = Pb210_CPAParts.size();
    if(fSaveTruth) FillTruth(Pb210_CPAAssn , Pb210_CPATrue , kPb210_CPA );
  }

  auto Bi210_CPATrue = evt.getHandle< std::vector<simb::MCTruth> > (fBi210_CPALabel);
  if(Bi210_CPATrue) {
    art::FindManyP<simb::MCParticle> Bi210_CPAAssn(Bi210_CPATrue,evt,fGEANTLabel);
    FillMyMaps( Bi210_CPAParts, Bi210_CPAAssn, Bi210_CPATrue );
    TotGen_Bi210_CPA = Bi210_CPAParts.size();
    if(fSaveTruth) FillTruth(Bi210_CPAAssn , Bi210_CPATrue , kBi210_CPA );
  }

  auto Co60_APATrue = evt.getHandle< std::vector<simb::MCTruth> > (fCo60_APALabel);
  if(Co60_APATrue) {
    art::FindManyP<simb::MCParticle> Co60_APAAssn(Co60_APATrue,evt,fGEANTLabel);
    FillMyMaps( Co60_APAParts, Co60_APAAssn, Co60_APATrue );
    TotGen_Co60_APA = Co60_APAParts.size();
    if(fSaveTruth) FillTruth(Co60_APAAssn , Co60_APATrue , kCo60_APA );
  }

  auto U238_APATrue = evt.getHandle< std::vector<simb::MCTruth> > (fU238_APALabel);
  if(U238_APATrue) {
    art::FindManyP<simb::MCParticle> U238_APAAssn(U238_APATrue,evt,fGEANTLabel);
    FillMyMaps( U238_APAParts, U238_APAAssn, U238_APATrue );
    TotGen_U238_APA = U238_APAParts.size();
    if(fSaveTruth) FillTruth(U238_APAAssn , U238_APATrue , kU238_APA );
  }

  auto Rn222_PDSTrue = evt.getHandle< std::vector<simb::MCTruth> > (fRn222_PDSLabel);
  if(Rn222_PDSTrue) {
    art::FindManyP<simb::MCParticle> Rn222_PDSAssn(Rn222_PDSTrue,evt,fGEANTLabel);
    FillMyMaps( Rn222_PDSParts, Rn222_PDSAssn, Rn222_PDSTrue );
    TotGen_Rn222_PDS = Rn222_PDSParts.size();
    if(fSaveTruth) FillTruth(Rn222_PDSAssn , Rn222_PDSTrue , kRn222_PDS );
  }

  auto NeutTrue = evt.getHandle< std::vector<simb::MCTruth> > (fNeutLabel);
  if(NeutTrue) {
    art::FindManyP<simb::MCParticle> NeutAssn(NeutTrue,evt,fGEANTLabel);
    FillMyMaps( NeutParts, NeutAssn, NeutTrue );
    TotGen_Neut = NeutParts.size();
    if(fSaveTruth) FillTruth(NeutAssn , NeutTrue , kNeut );
  }

  auto CavGamTrue = evt.getHandle< std::vector<simb::MCTruth> > (fCavGamLabel);
  if(CavGamTrue) {
    art::FindManyP<simb::MCParticle> CavGamAssn(CavGamTrue,evt,fGEANTLabel);
    FillMyMaps( CavGamParts, CavGamAssn, CavGamTrue );
    TotGen_CavGam = CavGamParts.size();
    if(fSaveTruth) FillTruth(CavGamAssn , CavGamTrue , kCavGam );
  }

  std::vector<simb::MCParticle> allTruthParts;
  for(auto& it: Ar39_LarParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Ar42_LarParts)
    allTruthParts.push_back(it.second);
  for(auto& it: K42_LarParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Kr85_LarParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Rn222_LarParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Po218_LarParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Pb214_LarParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Bi214_LarParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Pb210_LarParts)
    allTruthParts.push_back(it.second);
  for(auto& it: K40_CPAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Ar42_CPAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: U238_CPAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Po218_CPAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Pb214_CPAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Bi214_CPAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Pb210_CPAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Bi210_CPAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Co60_APAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: U238_APAParts)
    allTruthParts.push_back(it.second);
  for(auto& it: Rn222_PDSParts)
    allTruthParts.push_back(it.second);
  for(auto& it: NeutParts)
    allTruthParts.push_back(it.second);
  for(auto& it: CavGamParts)
    allTruthParts.push_back(it.second);

  std::map<PType, std::map< int, simb::MCParticle >&> PTypeToMap{
    {kMarl,      MarlParts     },
    {kAr39_Lar,  Ar39_LarParts },
    {kAr42_Lar,  Ar42_LarParts },
    {kK42_Lar,   K42_LarParts  },
    {kKr85_Lar,  Kr85_LarParts },
    {kRn222_Lar, Rn222_LarParts},
    {kPo218_Lar, Po218_LarParts},
    {kPb214_Lar, Pb214_LarParts},
    {kBi214_Lar, Bi214_LarParts},
    {kPb210_Lar, Pb210_LarParts},
    {kK40_CPA,   K40_CPAParts  },
    {kAr42_CPA,  Ar42_CPAParts },
    {kU238_CPA,  U238_CPAParts },
    {kPo218_CPA, Po218_CPAParts},
    {kPb214_CPA, Pb214_CPAParts},
    {kBi214_CPA, Bi214_CPAParts},
    {kPb210_CPA, Pb210_CPAParts},
    {kBi210_CPA, Bi210_CPAParts},
    {kCo60_APA,  Co60_APAParts },
    {kU238_APA,  U238_APAParts },
    {kRn222_PDS, Rn222_PDSParts},
    {kNeut,      NeutParts     },
    {kCavGam,    CavGamParts   }
  };

  std::map<PType,std::set<int>> PTypeToTrackID;

  for(auto const& it : PTypeToMap){
    const PType p=it.first;
    auto const& m=it.second;
    for(auto const& it2 : m){
      trkIDToPType.insert(std::make_pair(it2.first, p));
      PTypeToTrackID[p].insert(it2.first);
    }
  }

  mf::LogInfo(fname) << "THE EVENTS NUMBER IS: " << Event << std::endl;

  if (fSaveTPC) {
    auto reco_hits = evt.getHandle< std::vector<recob::Hit> >(fHitLabel);
    auto rawDigitsVecHandle = evt.getHandle< std::vector<raw::RawDigit> >(fRawDigitLabel);
    if ( reco_hits && rawDigitsVecHandle ) {
      std::vector< recob::Hit > ColHits_Oth      ;
      std::vector< recob::Hit > ColHits_Marl     ;
      std::vector< recob::Hit > ColHits_Ar39_Lar ;
      std::vector< recob::Hit > ColHits_Ar42_Lar ;
      std::vector< recob::Hit > ColHits_K42_Lar  ;
      std::vector< recob::Hit > ColHits_Kr85_Lar ;
      std::vector< recob::Hit > ColHits_Rn222_Lar;
      std::vector< recob::Hit > ColHits_Po218_Lar;
      std::vector< recob::Hit > ColHits_Pb214_Lar;
      std::vector< recob::Hit > ColHits_Bi214_Lar;
      std::vector< recob::Hit > ColHits_Pb210_Lar;
      std::vector< recob::Hit > ColHits_K40_CPA  ;
      std::vector< recob::Hit > ColHits_Ar42_CPA ;
      std::vector< recob::Hit > ColHits_U238_CPA ;
      std::vector< recob::Hit > ColHits_Po218_CPA;
      std::vector< recob::Hit > ColHits_Pb214_CPA;
      std::vector< recob::Hit > ColHits_Bi214_CPA;
      std::vector< recob::Hit > ColHits_Pb210_CPA;
      std::vector< recob::Hit > ColHits_Bi210_CPA;
      std::vector< recob::Hit > ColHits_Co60_APA ;
      std::vector< recob::Hit > ColHits_U238_APA ;
      std::vector< recob::Hit > ColHits_Rn222_PDS;
      std::vector< recob::Hit > ColHits_Neut     ;
      std::vector< recob::Hit > ColHits_CavGam   ;

      NTotHit = reco_hits->size();
      int colHitCount(0);
      int LoopHits = NTotHit;

      for(int hit = 0; hit < LoopHits; ++hit) {
        recob::Hit const& ThisHit = reco_hits->at(hit);
        if(ThisHit.View() == geo::kU || ThisHit.View() == geo::kV) {
          ++NIndHit;
        } else {
          ++NColHit;
        }
      }


      // Fill a set of the channels that we don't want to consider in the
      // adjacent-channel loop below. Do it here once so we don't have to
      // re-do it for every single hit
      std::set<int> badChannels;
      for(size_t i=0; i<rawDigitsVecHandle->size(); ++i) {
        int rawWireChannel=(*rawDigitsVecHandle)[i].Channel();
        std::vector<geo::WireID> adjacentwire = geo->ChannelToWire(rawWireChannel);

        if (adjacentwire.size() < 1 || adjacentwire[0].Plane == geo::kU ||
            adjacentwire[0].Plane == geo::kV){
          badChannels.insert(rawWireChannel);
        }
      }
      // std::cout << "Inserted " << badChannels.size() << " out of " << rawDigitsVecHandle->size() << " channels into set" << std::endl;

      auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
      for(int hit = 0; hit < LoopHits; ++hit) {
        recob::Hit const& ThisHit = reco_hits->at(hit);

        if (ThisHit.View() == 2) {
          std::vector<sim::TrackIDE> ThisHitIDE;
          //GETTING HOLD OF THE SIM::IDEs.

          std::vector<const sim::IDE*> ThisSimIDE;
          try {
            // HitToTrackIDEs opens a specific window around the hit. I want a
            // wider one, because the filtering can delay the hit. So this bit
            // is a copy of HitToTrackIDEs from the backtracker, with some
            // modification
            const double start = ThisHit.PeakTime()-20;
            const double end   = ThisHit.PeakTime()+ThisHit.RMS()+20;
            ThisHitIDE = bt_serv->ChannelToTrackIDEs(clockData, ThisHit.Channel(), start, end);

            // ThisHitIDE = bt_serv->HitToTrackIDEs(clockData,  ThisHit );
          } catch(...){
            // std::cout << "FIRST CATCH" << std::endl;
            firstCatch++;
            try {
              ThisSimIDE = bt_serv->HitToSimIDEs_Ps(clockData, ThisHit);
            } catch(...) {
               // std::cout << "SECOND CATCH" << std::endl;
              secondCatch++;
              // continue;
            }
            // continue;
          }

          // Get the simIDEs.
          try {
            ThisSimIDE = bt_serv->HitToSimIDEs_Ps(clockData, ThisHit);
          } catch(...) {
            // std::cout << "THIRD CATCH" << std::endl;
            thirdCatch++;
            // continue;
          }

          Hit_View.push_back(ThisHit.View());
          Hit_Size.push_back(ThisHit.EndTick() - ThisHit.StartTick());
          Hit_TPC .push_back(ThisHit.WireID().TPC);
          int channel = ThisHit.Channel();
          Hit_Chan.push_back(channel);

          if(fSaveNeighbourADCs)
            SaveNeighbourADC(channel,rawDigitsVecHandle, badChannels, ThisHit);

          double wire_start[3] = {0,0,0};
          double wire_end[3] = {0,0,0};
          auto& wgeo = geo->WireIDToWireGeo(ThisHit.WireID());
          wgeo.GetStart(wire_start);
          wgeo.GetEnd(wire_end);
          Hit_X_start.push_back(wire_start[0]);
          Hit_Y_start.push_back(wire_start[1]);
          Hit_Z_start.push_back(wire_start[2]);
          Hit_X_end  .push_back(wire_end[0]);
          Hit_Y_end  .push_back(wire_end[1]);
          Hit_Z_end  .push_back(wire_end[2]);
          Hit_Time   .push_back(ThisHit.PeakTime());
          Hit_RMS    .push_back(ThisHit.RMS());
          Hit_SADC   .push_back(ThisHit.SummedADC());
          Hit_Int    .push_back(ThisHit.Integral());
          Hit_Peak   .push_back(ThisHit.PeakAmplitude());
          Hit_True_nIDEs.push_back(ThisHitIDE.size());

          if(ThisHitIDE.size()==0)
            NHitNoBT++;

          //WHICH PARTICLE CONTRIBUTED MOST TO THE HIT.
          double TopEFrac = -DBL_MAX;

          Hit_True_EvEnergy.push_back(0);
          Hit_True_MainTrID.push_back(-1);
          for (size_t ideL=0; ideL < ThisHitIDE.size(); ++ideL) {
            Hit_True_TrackID.push_back(ThisHitIDE[ideL].trackID);
            for (size_t ipart=0; ipart<allTruthParts.size(); ++ipart) {

              if (allTruthParts[ipart].TrackId() == ThisHitIDE[ideL].trackID) {
                Hit_True_EvEnergy.at(colHitCount) += allTruthParts[ipart].E();
              }
            }
            if (ThisHitIDE[ideL].energyFrac > TopEFrac) {
              TopEFrac = ThisHitIDE[ideL].energyFrac;
              Hit_True_MainTrID.at(colHitCount) = ThisHitIDE[ideL].trackID;

            }
          }

          PType ThisPType = WhichParType(Hit_True_MainTrID.at(colHitCount));
          Hit_True_GenType.push_back(ThisPType);

          int thisMarleyIndex=-1;
          int MainTrID=Hit_True_MainTrID.at(colHitCount);
          if(ThisPType==kMarl && MainTrID!=0){
            auto const it=trkIDToMarleyIndex.find(MainTrID);
            if(it==trkIDToMarleyIndex.end()){
              mf::LogDebug(fname) << "Track ID " << MainTrID << " is not in Marley index map" << std::endl;
            }
            else{
              thisMarleyIndex=it->second;
            }
          }
          Hit_True_MarleyIndex.push_back(thisMarleyIndex);

          if(Hit_True_MainTrID[colHitCount] == -1)
          {
            Hit_True_X     .push_back(-1);
            Hit_True_Y     .push_back(-1);
            Hit_True_Z     .push_back(-1);
            Hit_True_Energy.push_back(-1);
            Hit_True_nElec .push_back(-1);
          }
          else
          {
            for(unsigned int i = 0; i < ThisSimIDE.size(); i++)
            {
              if(ThisSimIDE.at(i)->trackID==Hit_True_MainTrID[colHitCount])
              {
                Hit_True_X     .push_back(ThisSimIDE.at(i)->x           );
                Hit_True_Y     .push_back(ThisSimIDE.at(i)->y           );
                Hit_True_Z     .push_back(ThisSimIDE.at(i)->z           );
                Hit_True_Energy.push_back(ThisSimIDE.at(i)->energy      );
                Hit_True_nElec .push_back(ThisSimIDE.at(i)->numElectrons);
                break;
              }
            }
          }

          if      (ThisPType == 0)  { ColHits_Oth      .push_back( ThisHit ); }
          else if (ThisPType == 1)  { ColHits_Marl     .push_back( ThisHit ); }
          else if (ThisPType == 2)  { ColHits_Ar39_Lar .push_back( ThisHit ); }
          else if (ThisPType == 3)  { ColHits_Ar42_Lar .push_back( ThisHit ); }
          else if (ThisPType == 4)  { ColHits_K42_Lar  .push_back( ThisHit ); }
          else if (ThisPType == 5)  { ColHits_Kr85_Lar .push_back( ThisHit ); }
          else if (ThisPType == 6)  { ColHits_Rn222_Lar.push_back( ThisHit ); }
          else if (ThisPType == 7)  { ColHits_Po218_Lar.push_back( ThisHit ); }
          else if (ThisPType == 8)  { ColHits_Pb214_Lar.push_back( ThisHit ); }
          else if (ThisPType == 9)  { ColHits_Bi214_Lar.push_back( ThisHit ); }
          else if (ThisPType == 10) { ColHits_Pb210_Lar.push_back( ThisHit ); }
          else if (ThisPType == 11) { ColHits_K40_CPA  .push_back( ThisHit ); }
          else if (ThisPType == 12) { ColHits_Ar42_CPA .push_back( ThisHit ); }
          else if (ThisPType == 13) { ColHits_U238_CPA .push_back( ThisHit ); }
          else if (ThisPType == 14) { ColHits_Po218_CPA.push_back( ThisHit ); }
          else if (ThisPType == 15) { ColHits_Pb214_CPA.push_back( ThisHit ); }
          else if (ThisPType == 16) { ColHits_Bi214_CPA.push_back( ThisHit ); }
          else if (ThisPType == 17) { ColHits_Pb210_CPA.push_back( ThisHit ); }
          else if (ThisPType == 18) { ColHits_Bi210_CPA.push_back( ThisHit ); }
          else if (ThisPType == 19) { ColHits_Co60_APA .push_back( ThisHit ); }
          else if (ThisPType == 20) { ColHits_U238_APA .push_back( ThisHit ); }
          else if (ThisPType == 21) { ColHits_Rn222_PDS.push_back( ThisHit ); }
          else if (ThisPType == 22) { ColHits_Neut     .push_back( ThisHit ); }
          else if (ThisPType == 23) { ColHits_CavGam   .push_back( ThisHit ); }

          colHitCount++;
        }
      }
      mf::LogInfo(fname) << "Total of:\n"
                         << " - Other    : " << ColHits_Oth.size() << " col plane hits\n"
                         << " - Marl     : " << TotGen_Marl        << " true parts\t| " << ColHits_Marl     .size() << " col plane hits\n"
                         << " - Ar39_Lar : " << TotGen_Ar39_Lar    << " true parts\t| " << ColHits_Ar39_Lar .size() << " col plane hits\n"
                         << " - Ar42_Lar : " << TotGen_Ar42_Lar    << " true parts\t| " << ColHits_Ar42_Lar .size() << " col plane hits\n"
                         << " - K42_Lar  : " << TotGen_K42_Lar     << " true parts\t| " << ColHits_K42_Lar  .size() << " col plane hits\n"
                         << " - Kr85_Lar : " << TotGen_Kr85_Lar    << " true parts\t| " << ColHits_Kr85_Lar .size() << " col plane hits\n"
                         << " - Rn222_Lar: " << TotGen_Rn222_Lar   << " true parts\t| " << ColHits_Rn222_Lar.size() << " col plane hits\n"
                         << " - Po218_Lar: " << TotGen_Po218_Lar   << " true parts\t| " << ColHits_Po218_Lar.size() << " col plane hits\n"
                         << " - Pb214_Lar: " << TotGen_Pb214_Lar   << " true parts\t| " << ColHits_Pb214_Lar.size() << " col plane hits\n"
                         << " - Bi214_Lar: " << TotGen_Bi214_Lar   << " true parts\t| " << ColHits_Bi214_Lar.size() << " col plane hits\n"
                         << " - Pb210_Lar: " << TotGen_Pb210_Lar   << " true parts\t| " << ColHits_Pb210_Lar.size() << " col plane hits\n"
                         << " - K40_CPA  : " << TotGen_K40_CPA     << " true parts\t| " << ColHits_K40_CPA  .size() << " col plane hits\n"
                         << " - Ar42_CPA : " << TotGen_Ar42_CPA    << " true parts\t| " << ColHits_Ar42_CPA .size() << " col plane hits\n"
                         << " - U238_CPA : " << TotGen_U238_CPA    << " true parts\t| " << ColHits_U238_CPA .size() << " col plane hits\n"
                         << " - Po218_CPA: " << TotGen_Po218_CPA   << " true parts\t| " << ColHits_Po218_CPA.size() << " col plane hits\n"
                         << " - Pb214_CPA: " << TotGen_Pb214_CPA   << " true parts\t| " << ColHits_Pb214_CPA.size() << " col plane hits\n"
                         << " - Bi214_CPA: " << TotGen_Bi214_CPA   << " true parts\t| " << ColHits_Bi214_CPA.size() << " col plane hits\n"
                         << " - Pb210_CPA: " << TotGen_Pb210_CPA   << " true parts\t| " << ColHits_Pb210_CPA.size() << " col plane hits\n"
                         << " - Bi210_CPA: " << TotGen_Bi210_CPA   << " true parts\t| " << ColHits_Bi210_CPA.size() << " col plane hits\n"
                         << " - Co60_APA : " << TotGen_Co60_APA    << " true parts\t| " << ColHits_Co60_APA .size() << " col plane hits\n"
                         << " - U238_APA : " << TotGen_U238_APA    << " true parts\t| " << ColHits_U238_APA .size() << " col plane hits\n"
                         << " - Rn222_PDS: " << TotGen_Rn222_PDS   << " true parts\t| " << ColHits_Rn222_PDS.size() << " col plane hits\n"
                         << " - Neut     : " << TotGen_Neut        << " true parts\t| " << ColHits_Neut     .size() << " col plane hits\n"
                         << " - CavGam   : " << TotGen_CavGam      << " true parts\t| " << ColHits_CavGam   .size() << " col plane hits\n";

    } else {
      mf::LogError(fname) << "Requested to save wire hits, but cannot load any wire hits";
      throw art::Exception(art::errors::NotFound) << "Requested to save wire hits, but cannot load any wire hits\n";
    }
  }

  if (fSavePDS) {

    auto OpHitHandle = evt.getHandle< std::vector< recob::OpHit > >(fOpHitModuleLabel);
    std::vector<art::Ptr<recob::OpHit> > ophitlist;
    std::map<PType, std::vector<art::Ptr<recob::OpHit> > > map_of_ophit;

    if (OpHitHandle) {
      art::fill_ptr_vector(ophitlist, OpHitHandle);

      mf::LogDebug(fname) << "There are " << ophitlist.size() << " optical hits in the event." << std::endl;

      for(size_t i = 0; i < ophitlist.size(); ++i)
      {
        std::vector<sim::TrackSDP> vec_tracksdp = pbt_serv->OpHitToTrackSDPs(ophitlist.at(i));
        PType gen = kUnknown;

        std::sort(vec_tracksdp.begin(), vec_tracksdp.end(),
                  [](const sim::TrackSDP& a, const sim::TrackSDP& b) -> bool { return a.energyFrac > b.energyFrac; });

        for (size_t iSDP=0; iSDP<vec_tracksdp.size(); ++iSDP) {
          PDS_OpHit_True_TrackIDAll.push_back(vec_tracksdp[iSDP].trackID);
          PDS_OpHit_True_GenTypeAll.push_back(WhichParType(vec_tracksdp[iSDP].trackID));
          PDS_OpHit_True_EnergyAll .push_back(vec_tracksdp[iSDP].energy);
          PDS_OpHit_True_IndexAll  .push_back((int)i);
        }

        if (vec_tracksdp.size()>0){
          int MainTrID = vec_tracksdp[0].trackID;
          PDS_OpHit_True_TrackID.push_back(vec_tracksdp[0].trackID);
          gen = WhichParType(vec_tracksdp[0].trackID);
          PDS_OpHit_True_GenType.push_back(gen);
          int thisMarleyIndex;
          if(gen==kMarl && MainTrID!=0){
            auto const it=trkIDToMarleyIndex.find(MainTrID);
            if(it==trkIDToMarleyIndex.end()){
              mf::LogDebug(fname) << "Track ID " << MainTrID << " is not in Marley index map" << std::endl;
            }
            else{
              thisMarleyIndex=it->second;
            }
          }
          PDS_OpHit_True_Index  .push_back(thisMarleyIndex);
          PDS_OpHit_True_Energy .push_back(vec_tracksdp[0].energy);
        } else {
          PDS_OpHit_True_TrackID.push_back(-1);
          PDS_OpHit_True_GenType.push_back(kUnknown);
          PDS_OpHit_True_Energy .push_back(-1);
        }

        map_of_ophit[gen].push_back(ophitlist.at(i));

        double xyz_optdet[3]={0,0,0};
        double xyz_world [3]={0,0,0};

        geo->OpDetGeoFromOpChannel(ophitlist[i]->OpChannel()).LocalToWorld(xyz_optdet,xyz_world);
        PDS_OpHit_OpChannel   .push_back(ophitlist[i]->OpChannel());
        PDS_OpHit_X           .push_back(xyz_world[0]);
        PDS_OpHit_Y           .push_back(xyz_world[1]);
        PDS_OpHit_Z           .push_back(xyz_world[2]);
        PDS_OpHit_PeakTimeAbs .push_back(ophitlist[i]->PeakTimeAbs());
        PDS_OpHit_PeakTime    .push_back(ophitlist[i]->PeakTime());
        PDS_OpHit_Frame       .push_back(ophitlist[i]->Frame());
        PDS_OpHit_Width       .push_back(ophitlist[i]->Width());
        PDS_OpHit_Area        .push_back(ophitlist[i]->Area());
        PDS_OpHit_Amplitude   .push_back(ophitlist[i]->Amplitude());
        PDS_OpHit_PE          .push_back(ophitlist[i]->PE());
        PDS_OpHit_FastToTotal .push_back(ophitlist[i]->FastToTotal());
      }
      mf::LogInfo(fname) << "Total of :\n"
                         << " - Other    : " << map_of_ophit[kUnknown].size() << " opt hits\n"
                         << " - Marl     : " << TotGen_Marl      << " true parts\t| " << map_of_ophit[kMarl]     .size() << " opt hits\n"
                         << " - Ar39_Lar : " << TotGen_Ar39_Lar  << " true parts\t| " << map_of_ophit[kAr39_Lar] .size() << " opt hits\n"
                         << " - Ar42_Lar : " << TotGen_Ar42_Lar  << " true parts\t| " << map_of_ophit[kAr42_Lar] .size() << " opt hits\n"
                         << " - K42_Lar  : " << TotGen_K42_Lar   << " true parts\t| " << map_of_ophit[kK42_Lar]  .size() << " opt hits\n"
                         << " - Kr85_Lar : " << TotGen_Kr85_Lar  << " true parts\t| " << map_of_ophit[kKr85_Lar] .size() << " opt hits\n"
                         << " - Rn222_Lar: " << TotGen_Rn222_Lar << " true parts\t| " << map_of_ophit[kRn222_Lar].size() << " opt hits\n"
                         << " - Po218_Lar: " << TotGen_Po218_Lar << " true parts\t| " << map_of_ophit[kPo218_Lar].size() << " opt hits\n"
                         << " - Pb214_Lar: " << TotGen_Pb214_Lar << " true parts\t| " << map_of_ophit[kPb214_Lar].size() << " opt hits\n"
                         << " - Bi214_Lar: " << TotGen_Bi214_Lar << " true parts\t| " << map_of_ophit[kBi214_Lar].size() << " opt hits\n"
                         << " - Pb210_Lar: " << TotGen_Pb210_Lar << " true parts\t| " << map_of_ophit[kPb210_Lar].size() << " opt hits\n"
                         << " - K40_CPA  : " << TotGen_K40_CPA   << " true parts\t| " << map_of_ophit[kK40_CPA]  .size() << " opt hits\n"
                         << " - Ar42_CPA : " << TotGen_Ar42_CPA  << " true parts\t| " << map_of_ophit[kAr42_CPA] .size() << " opt hits\n"
                         << " - U238_CPA : " << TotGen_U238_CPA  << " true parts\t| " << map_of_ophit[kU238_CPA] .size() << " opt hits\n"
                         << " - Po218_CPA: " << TotGen_Po218_CPA << " true parts\t| " << map_of_ophit[kPo218_CPA] .size() << " opt hits\n"
                         << " - Pb214_CPA: " << TotGen_Pb214_CPA << " true parts\t| " << map_of_ophit[kPb214_CPA] .size() << " opt hits\n"
                         << " - Bi214_CPA: " << TotGen_Bi214_CPA << " true parts\t| " << map_of_ophit[kBi214_CPA] .size() << " opt hits\n"
                         << " - Pb210_CPA: " << TotGen_Pb210_CPA << " true parts\t| " << map_of_ophit[kPb210_CPA] .size() << " opt hits\n"
                         << " - Bi210_CPA: " << TotGen_Bi210_CPA << " true parts\t| " << map_of_ophit[kBi210_CPA] .size() << " opt hits\n"
                         << " - Co60_APA : " << TotGen_Co60_APA  << " true parts\t| " << map_of_ophit[kCo60_APA] .size() << " opt hits\n"
                         << " - U238_APA : " << TotGen_U238_APA  << " true parts\t| " << map_of_ophit[kU238_APA] .size() << " opt hits\n"
                         << " - Rn222_PDS: " << TotGen_Rn222_PDS << " true parts\t| " << map_of_ophit[kRn222_PDS].size() << " opt hits\n"
                         << " - Neut     : " << TotGen_Neut      << " true parts\t| " << map_of_ophit[kNeut]     .size() << " opt hits\n"
                         << " - CavGam   : " << TotGen_CavGam    << " true parts\t| " << map_of_ophit[kCavGam]   .size() << " opt hits\n";

    }
    else {
      mf::LogError(fname) << "Requested to save optical hits, but cannot load any ophits";
      throw art::Exception(art::errors::NotFound) << "Requested to save optical hits, but cannot load any optical hits\n";
    }
  }

  if(fSaveIDEs) SaveIDEs(evt);


  fSNAnaTree->Fill();
}

void SNAna::FillTruth(const art::FindManyP<simb::MCParticle> Assn,
                      const art::Handle<std::vector<simb::MCTruth>>& Hand,
                      const PType type) {

  for(size_t i1=0; i1<Hand->size(); ++i1) {
    for ( size_t i2=0; i2 < Assn.at(i1).size(); ++i2 ) {
      const simb::MCParticle ThisPar = (*Assn.at(i1).at(i2));
      True_Bck_Mode      .push_back(type);
      True_Bck_PDG       .push_back(ThisPar.PdgCode   ());
      True_Bck_ID        .push_back(ThisPar.TrackId   ());
      True_Bck_Mother    .push_back(ThisPar.Mother    ());
      True_Bck_P         .push_back(ThisPar.P         ());
      True_Bck_VertX     .push_back(ThisPar.Vx        ());
      True_Bck_VertY     .push_back(ThisPar.Vy        ());
      True_Bck_VertZ     .push_back(ThisPar.Vz        ());
      True_Bck_Time      .push_back(ThisPar.T         ());
      True_Bck_Energy    .push_back(ThisPar.E         ());
      True_Bck_EndX      .push_back(ThisPar.EndX      ());
      True_Bck_EndY      .push_back(ThisPar.EndY      ());
      True_Bck_EndZ      .push_back(ThisPar.EndZ      ());
      True_Bck_EndT      .push_back(ThisPar.EndT      ());
      True_Bck_EndE      .push_back(ThisPar.EndE      ());
      True_Bck_EndProcess.push_back(ThisPar.EndProcess());
      True_Bck_Process   .push_back(ThisPar.Process   ());
    }
  }
}

void SNAna::FillMyMaps(std::map< int, simb::MCParticle> &MyMap,
                       art::FindManyP<simb::MCParticle> Assn,
                       art::Handle< std::vector<simb::MCTruth> > Hand,
                       std::map<int, int>* indexMap) {
  for ( size_t L1=0; L1 < Hand->size(); ++L1 ) {
    for ( size_t L2=0; L2 < Assn.at(L1).size(); ++L2 ) {
      const simb::MCParticle ThisPar = (*Assn.at(L1).at(L2));
      MyMap[ThisPar.TrackId()] = ThisPar;
      if(indexMap) indexMap->insert({ThisPar.TrackId(), L1});
    }
  }
  return;
}


PType SNAna::WhichParType( int TrID )
{
  PType ThisPType = kUnknown;
  auto const& it=trkIDToPType.find(TrID);
  if(it!=trkIDToPType.end()){
    ThisPType=it->second;
  }
  return ThisPType;

}


bool SNAna::InMyMap( int TrID, std::map< int, simb::MCParticle> ParMap )
{
  std::map<int, simb::MCParticle>::iterator ParIt;
  ParIt = ParMap.find( TrID );
  if ( ParIt != ParMap.end() ) {
    return true;
  } else
    return false;
}

//......................................................
void SNAna::SaveIDEs(art::Event const & evt)
{
  auto allParticles = evt.getValidHandle<std::vector<simb::MCParticle> >(fGEANTLabel);
  art::FindMany<simb::MCTruth> assn(allParticles,evt,fGEANTLabel);
  std::map<int, const simb::MCTruth*> idToTruth;
  for(size_t i=0; i<allParticles->size(); ++i){
    const simb::MCParticle& particle=allParticles->at(i);
    const std::vector<const simb::MCTruth*> truths=assn.at(i);
    if(truths.size()==1){
      idToTruth[particle.TrackId()]=truths[0];
    }
    else{
      mf::LogDebug("DAQSimAna") << "Particle " << particle.TrackId() << " has " << truths.size() << " truths";
      idToTruth[particle.TrackId()]=nullptr;
    }
  }

  // Get the SimChannels so we can see where the actual energy depositions were
  auto& simchs=*evt.getValidHandle<std::vector<sim::SimChannel>>("largeant");

  for(auto&& simch: simchs){
    // We only care about collection channels
    if(geo->SignalType(simch.Channel())!=geo::kCollection) continue;

    // The IDEs record energy depositions at every tick, but
    // mostly we have several depositions at contiguous times. So
    // we're going to save just one output IDE for each contiguous
    // block of hits on a channel. Each item in vector is a list
    // of (TDC, IDE*) for contiguous-in-time IDEs
    std::vector<std::vector<std::pair<int, const sim::IDE*> > > contigIDEs;
    int prevTDC=0;
    for (const auto& TDCinfo: simch.TDCIDEMap()) {
      // Do we need to start a new group of IDEs? Yes if this is
      // the first IDE in this channel. Yes if this IDE is not
      // contiguous with the previous one
      if(contigIDEs.empty() || TDCinfo.first-prevTDC>5){
        contigIDEs.push_back(std::vector<std::pair<int, const sim::IDE*> >());
      }
      std::vector<std::pair<int, const sim::IDE*> >& currentIDEs=contigIDEs.back();

      // Add all the current tick's IDEs to the list
      for (const sim::IDE& ide: TDCinfo.second) {
        currentIDEs.push_back(std::make_pair(TDCinfo.first, &ide));
      }
      prevTDC=TDCinfo.first;
    }

    for(auto const& contigs : contigIDEs){
      float energy=0;
      float electrons=0;
      int startTime=99999;
      int endTime=0;
      std::map<PType, float> ptypeToEnergy;
      for(auto const& timeide : contigs){
        const int tdc=timeide.first;
        startTime=std::min(tdc, startTime);
        endTime=std::max(tdc, endTime);
        const sim::IDE& ide=*timeide.second;
        const float thisEnergy=ide.energy;
        const PType thisPType=WhichParType(std::abs(ide.trackID));
        energy+=thisEnergy;
        electrons+=ide.numElectrons;
        ptypeToEnergy[thisPType]+=thisEnergy;
      }
      float bestEnergy=0;
      PType bestPType=kUnknown;
      for(auto const& it : ptypeToEnergy){
        if(it.second>bestEnergy){
          bestEnergy=it.second;
          bestPType=it.first;
        }
      }
      // Ignore anything past the end of the readout window
      if(startTime<4492){
        IDEChannel.push_back(simch.Channel());
        IDEStartTime.push_back(startTime);
        IDEEndTime.push_back(endTime);
        IDEEnergy.push_back(energy);
        IDEElectrons.push_back(electrons);
        IDEParticle.push_back(bestPType);
      }
    } // loop over our compressed IDEs
  } // loop over SimChannels
}


void SNAna::SaveNeighbourADC(int channel,
                             art::Handle< std::vector<raw::RawDigit> >rawDigitsVecHandle,
                             std::set<int>badChannels,
                             recob::Hit const& ThisHit) {

  Hit_AdjM5SADC.push_back(0);
  Hit_AdjM2SADC.push_back(0);
  Hit_AdjM1SADC.push_back(0);
  Hit_AdjP1SADC.push_back(0);
  Hit_AdjP2SADC.push_back(0);
  Hit_AdjP5SADC.push_back(0);
  Hit_AdjM5Chan.push_back(0);
  Hit_AdjM2Chan.push_back(0);
  Hit_AdjM1Chan.push_back(0);
  Hit_AdjP1Chan.push_back(0);
  Hit_AdjP2Chan.push_back(0);
  Hit_AdjP5Chan.push_back(0);
  int colHitCount = Hit_AdjM1Chan.size()-1;

  // A vector for us to uncompress the rawdigits into. Create it
  // outside the loop here so that we only have to allocate it once
  raw::RawDigit::ADCvector_t ADCs((*rawDigitsVecHandle)[0].Samples());

  std::vector<geo::WireID> wids = geo->ChannelToWire(channel);
  for(size_t i=0; i<rawDigitsVecHandle->size(); ++i)
  {
    int rawWireChannel=(*rawDigitsVecHandle)[i].Channel();
    const int chanDiff=rawWireChannel-channel;
    if(abs(chanDiff)!=1 && abs(chanDiff)!=2 && abs(chanDiff)!=5) continue;

    if(badChannels.find(rawWireChannel)!=badChannels.end()) continue;

    switch(chanDiff)
    {
    case -5:
      Hit_AdjM5SADC[colHitCount] = 0;
      Hit_AdjM5Chan[colHitCount] = rawWireChannel;
      raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                      (*rawDigitsVecHandle)[i].Compression());
      for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
        Hit_AdjM5SADC[colHitCount]+=ADCs[i];
      break;
    case -2:
      Hit_AdjM2SADC[colHitCount] = 0;
      Hit_AdjM2Chan[colHitCount] = rawWireChannel;
      raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                      (*rawDigitsVecHandle)[i].Compression());
      for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
        Hit_AdjM2SADC[colHitCount]+=ADCs[i];
      break;
    case -1:
      Hit_AdjM1SADC[colHitCount] = 0;
      Hit_AdjM1Chan[colHitCount] = rawWireChannel;
      raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                      (*rawDigitsVecHandle)[i].Compression());
      for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
        Hit_AdjM1SADC[colHitCount]+=ADCs[i];
      break;
    case  1:
      Hit_AdjP1SADC[colHitCount] = 0;
      Hit_AdjP1Chan[colHitCount] = rawWireChannel;
      raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                      (*rawDigitsVecHandle)[i].Compression());
      for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
        Hit_AdjP1SADC[colHitCount]+=ADCs[i];
      break;
    case  2:
      Hit_AdjP2SADC[colHitCount] = 0;
      Hit_AdjP2Chan[colHitCount] = rawWireChannel;
      raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                      (*rawDigitsVecHandle)[i].Compression());
      for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
        Hit_AdjP2SADC[colHitCount]+=ADCs[i];
      break;
    case  5:
      Hit_AdjP5SADC[colHitCount] = 0;
      Hit_AdjP5Chan[colHitCount] = rawWireChannel;
      raw::Uncompress((*rawDigitsVecHandle)[i].ADCs(), ADCs,
                      (*rawDigitsVecHandle)[i].Compression());
      for(int i=ThisHit.StartTick(); i<=ThisHit.EndTick();++i)
        Hit_AdjP5SADC[colHitCount]+=ADCs[i];
      break;
    default:
      break;
    }
  }
}


void SNAna::endJob()
{
  mf::LogDebug(fname) << firstCatch << " " << secondCatch << " " << thirdCatch << std::endl;
}

DEFINE_ART_MODULE(SNAna)
