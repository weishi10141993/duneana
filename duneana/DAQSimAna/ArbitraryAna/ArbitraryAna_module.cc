#include <functional>

//ROOT includes
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TMath.h"
#include "TGeoMatrix.h"

//LArSoft includes
#include "larcore/Geometry/Geometry.h"

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
//#include "art/Framework/Services/Optional/TFileDirectory.h"
//#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

class ArbitraryAna : public art::EDAnalyzer {

public:
  explicit ArbitraryAna(fhicl::ParameterSet const & p);

  ArbitraryAna(ArbitraryAna const &) = delete;
  ArbitraryAna(ArbitraryAna &&) = delete;
  ArbitraryAna & operator = (ArbitraryAna const &) = delete;
  ArbitraryAna & operator = (ArbitraryAna &&) = delete;

  void analyze(art::Event const & evt) override;
  void reconfigure(fhicl::ParameterSet const & p);
  void beginJob() override;
  void endJob() {
    std::cout << "first catch  " << firstCatch << std::endl;
    std::cout << "second catch " << secondCatch << std::endl;
    std::cout << "third catch  " << thirdCatch << std::endl;
  };

private:
  void ResetVariables();

  int firstCatch;
  int secondCatch;
  int thirdCatch;

  std::string fHitLabel;
  std::string fOpHitModuleLabel;

  std::string fGEANTLabel;
  std::string fGenLabel;
  
  std::map<int, bool> TrackIDGenerator;
  std::map<int, int>  TrackIDGeneratorIndex;

  bool fSaveTPC;
  bool fSavePDS;

  TTree* fArbitraryAnaTree;

  int Run;
  int SubRun;
  int Event;

  int NTotHit;
  int NColHit;
  int NIndHit;
  int NHitNoBT;

  std::vector<int>    Hit_View;
  std::vector<int>    Hit_Size;
  std::vector<int>    Hit_TPC;
  std::vector<int>    Hit_Chan;
  std::vector<double> Hit_X_start;
  std::vector<double> Hit_Y_start;
  std::vector<double> Hit_Z_start;
  std::vector<double> Hit_X_end;
  std::vector<double> Hit_Y_end;
  std::vector<double> Hit_Z_end;
  std::vector<double> Hit_Time;
  std::vector<double> Hit_RMS;
  std::vector<double> Hit_SADC;
  std::vector<double> Hit_Int;
  std::vector<double> Hit_Peak;
  
  std::vector<bool>   Hit_True_GenType;
  std::vector<int>    Hit_True_TrackID;
  std::vector<int>    Hit_True_GenIndex;
  std::vector<double> Hit_True_Energy;
  std::vector<double> Hit_True_nElec;
  std::vector<int>    Hit_True_nIDEs;

  std::vector<int>    PDS_OpHit_OpChannel;
  std::vector<double> PDS_OpHit_X;
  std::vector<double> PDS_OpHit_Y;
  std::vector<double> PDS_OpHit_Z;
  std::vector<double> PDS_OpHit_PeakTimeAbs;
  std::vector<double> PDS_OpHit_PeakTime;
  std::vector<int>    PDS_OpHit_Frame;
  std::vector<double> PDS_OpHit_Width;
  std::vector<double> PDS_OpHit_Area;
  std::vector<double> PDS_OpHit_Amplitude;
  std::vector<double> PDS_OpHit_PE;
  std::vector<double> PDS_OpHit_FastToTotal;

  std::vector<bool>   PDS_OpHit_True_GenType;
  std::vector<int>    PDS_OpHit_True_TrackID;
  std::vector<int>    PDS_OpHit_True_GenIndex;
  std::vector<double> PDS_OpHit_True_Energy;
  std::vector<double> PDS_OpHit_True_PE;

  std::vector<int>    True_PDG;
  std::vector<int>    True_ID;
  std::vector<double> True_E;
  std::vector<double> True_VertX;
  std::vector<double> True_VertY;
  std::vector<double> True_VertZ;
  std::vector<double> True_VertexT;
  std::vector<double> True_Px;
  std::vector<double> True_Py;
  std::vector<double> True_Pz;
  std::vector<double> True_Dirx;
  std::vector<double> True_Diry;
  std::vector<double> True_Dirz;
  std::vector<double> True_Time;

  std::vector<int>    True_Geant4_PDG;
  std::vector<int>    True_Geant4_ID;
  std::vector<int>    True_Geant4_MotherID;
  std::vector<double> True_Geant4_E;
  std::vector<double> True_Geant4_VertX;
  std::vector<double> True_Geant4_VertY;
  std::vector<double> True_Geant4_VertZ;
  std::vector<double> True_Geant4_Time;

  int TotGen;

  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::PhotonBackTrackerService> pbt_serv;

};

ArbitraryAna::ArbitraryAna(fhicl::ParameterSet const & p): EDAnalyzer(p) {
  this->reconfigure(p);
  firstCatch  = 0;
  secondCatch = 0;
  thirdCatch  = 0;
}


void ArbitraryAna::reconfigure(fhicl::ParameterSet const & p) {
  
  fHitLabel           = p.get<std::string>("HitLabel",         "hitfd");
  fOpHitModuleLabel   = p.get<std::string>("OpHitModuleLabel", "ophit");

  fGEANTLabel = p.get<std::string> ("GEANT4Label", "largeant");
  fGenLabel   = p.get<std::string> ("GenLabel",    "generator");
  fSaveTPC    = p.get<bool>("SaveTPC", 1);
  fSavePDS    = p.get<bool>("SavePDS", 1);
  
}


void ArbitraryAna::ResetVariables() {

  TrackIDGenerator.clear();
  TrackIDGeneratorIndex.clear();
  
  Hit_View.clear();
  Hit_Size.clear();
  Hit_TPC.clear();
  Hit_Chan.clear();
  Hit_X_start.clear();
  Hit_Y_start.clear();
  Hit_Z_start.clear();
  Hit_X_end.clear();
  Hit_Y_end.clear();
  Hit_Z_end.clear();
  Hit_Time.clear();
  Hit_RMS.clear();
  Hit_SADC.clear();
  Hit_Int.clear();
  Hit_Peak.clear();
  
  Hit_True_GenType.clear();
  Hit_True_TrackID.clear();
  Hit_True_GenIndex.clear();
  Hit_True_Energy.clear();
  Hit_True_nElec.clear();
  Hit_True_nIDEs.clear();

  PDS_OpHit_OpChannel.clear();
  PDS_OpHit_X.clear();
  PDS_OpHit_Y.clear();
  PDS_OpHit_Z.clear();
  PDS_OpHit_PeakTimeAbs.clear();
  PDS_OpHit_PeakTime.clear();
  PDS_OpHit_Frame.clear();
  PDS_OpHit_Width.clear();
  PDS_OpHit_Area.clear();
  PDS_OpHit_Amplitude.clear();
  PDS_OpHit_PE.clear();
  PDS_OpHit_FastToTotal.clear();

  PDS_OpHit_True_GenType.clear();
  PDS_OpHit_True_TrackID.clear();
  PDS_OpHit_True_GenIndex.clear();
  PDS_OpHit_True_Energy.clear();
  PDS_OpHit_True_PE.clear();

  True_PDG.clear();
  True_ID.clear();
  True_E.clear();
  True_VertX.clear();
  True_VertY.clear();
  True_VertZ.clear();
  True_VertexT.clear();
  True_Px.clear();
  True_Py.clear();
  True_Pz.clear();
  True_Dirx.clear();
  True_Diry.clear();
  True_Dirz.clear();
  True_Time.clear();

  True_Geant4_PDG.clear();
  True_Geant4_ID.clear();
  True_Geant4_MotherID.clear();
  True_Geant4_E.clear();
  True_Geant4_VertX.clear();
  True_Geant4_VertY.clear();
  True_Geant4_VertZ.clear();
  True_Geant4_Time.clear();

  Run = SubRun = Event = -1;
  TotGen = 0;

  NTotHit  = 0;
  NColHit  = 0;
  NIndHit  = 0;
  NHitNoBT = 0;

}


void ArbitraryAna::beginJob() {

  art::ServiceHandle<art::TFileService> tfs;

  fArbitraryAnaTree = tfs->make<TTree>("SNSimTree","SN simulation analysis tree");

  fArbitraryAnaTree->Branch("Run"     , &Run     , "Run/I"     );
  fArbitraryAnaTree->Branch("SubRun"  , &SubRun  , "SubRun/I"  );
  fArbitraryAnaTree->Branch("Event"   , &Event   , "Event/I"   );
  fArbitraryAnaTree->Branch("NTotHit" , &NTotHit , "NTotHit/I" );
  fArbitraryAnaTree->Branch("NColHit" , &NColHit , "NColHit/I" );
  fArbitraryAnaTree->Branch("NIndHit" , &NIndHit , "NIndHit/I" );
  fArbitraryAnaTree->Branch("NHitNoBT", &NHitNoBT, "NHitNoBT/I");

  if (fSaveTPC) {
    fArbitraryAnaTree->Branch("Hit_View"   , &Hit_View   );
    fArbitraryAnaTree->Branch("Hit_Size"   , &Hit_Size   );
    fArbitraryAnaTree->Branch("Hit_TPC"    , &Hit_TPC    );
    fArbitraryAnaTree->Branch("Hit_Chan"   , &Hit_Chan   );
    fArbitraryAnaTree->Branch("Hit_X_start", &Hit_X_start);
    fArbitraryAnaTree->Branch("Hit_Y_start", &Hit_Y_start);
    fArbitraryAnaTree->Branch("Hit_Z_start", &Hit_Z_start);
    fArbitraryAnaTree->Branch("Hit_X_end"  , &Hit_X_end  );
    fArbitraryAnaTree->Branch("Hit_Y_end"  , &Hit_Y_end  );
    fArbitraryAnaTree->Branch("Hit_Z_end"  , &Hit_Z_end  );
    fArbitraryAnaTree->Branch("Hit_Time"   , &Hit_Time   );
    fArbitraryAnaTree->Branch("Hit_RMS"    , &Hit_RMS    );
    fArbitraryAnaTree->Branch("Hit_SADC"   , &Hit_SADC   );
    fArbitraryAnaTree->Branch("Hit_Int"    , &Hit_Int    );
    fArbitraryAnaTree->Branch("Hit_Peak"   , &Hit_Peak   );

    fArbitraryAnaTree->Branch("Hit_True_GenType" , &Hit_True_GenType );
    fArbitraryAnaTree->Branch("Hit_True_TrackID" , &Hit_True_TrackID );
    fArbitraryAnaTree->Branch("Hit_True_GenIndex", &Hit_True_GenIndex);
    fArbitraryAnaTree->Branch("Hit_True_Energy"  , &Hit_True_Energy  );
    fArbitraryAnaTree->Branch("Hit_True_nElec"   , &Hit_True_nElec   );
    fArbitraryAnaTree->Branch("Hit_True_nIDEs"   , &Hit_True_nIDEs   );
  }

  if (fSavePDS) {
    fArbitraryAnaTree->Branch("PDS_OpHit_OpChannel"  , &PDS_OpHit_OpChannel  );
    fArbitraryAnaTree->Branch("PDS_OpHit_X"          , &PDS_OpHit_X          );
    fArbitraryAnaTree->Branch("PDS_OpHit_Y"          , &PDS_OpHit_Y          );
    fArbitraryAnaTree->Branch("PDS_OpHit_Z"          , &PDS_OpHit_Z          );
    fArbitraryAnaTree->Branch("PDS_OpHit_PeakTimeAbs", &PDS_OpHit_PeakTimeAbs);
    fArbitraryAnaTree->Branch("PDS_OpHit_PeakTime"   , &PDS_OpHit_PeakTime   );
    fArbitraryAnaTree->Branch("PDS_OpHit_Frame"      , &PDS_OpHit_Frame      );
    fArbitraryAnaTree->Branch("PDS_OpHit_Width"      , &PDS_OpHit_Width      );
    fArbitraryAnaTree->Branch("PDS_OpHit_Area"       , &PDS_OpHit_Area       );
    fArbitraryAnaTree->Branch("PDS_OpHit_Amplitude"  , &PDS_OpHit_Amplitude  );
    fArbitraryAnaTree->Branch("PDS_OpHit_PE"         , &PDS_OpHit_PE         );
    fArbitraryAnaTree->Branch("PDS_OpHit_FastToTotal", &PDS_OpHit_FastToTotal);

    fArbitraryAnaTree->Branch("PDS_OpHit_True_GenType" , &PDS_OpHit_True_GenType );
    fArbitraryAnaTree->Branch("PDS_OpHit_True_TrackID" , &PDS_OpHit_True_TrackID );
    fArbitraryAnaTree->Branch("PDS_OpHit_True_GenIndex", &PDS_OpHit_True_GenIndex);
    fArbitraryAnaTree->Branch("PDS_OpHit_True_Energy"  , &PDS_OpHit_True_Energy  );
    fArbitraryAnaTree->Branch("PDS_OpHit_True_PE"      , &PDS_OpHit_True_PE      );
  }
  
  fArbitraryAnaTree->Branch("True_PDG"    , &True_PDG    );
  fArbitraryAnaTree->Branch("True_ID"     , &True_ID     );
  fArbitraryAnaTree->Branch("True_E"      , &True_E      );
  fArbitraryAnaTree->Branch("True_VertX"  , &True_VertX  );
  fArbitraryAnaTree->Branch("True_VertY"  , &True_VertY  );
  fArbitraryAnaTree->Branch("True_VertZ"  , &True_VertZ  );
  fArbitraryAnaTree->Branch("True_VertexT", &True_VertexT);
  fArbitraryAnaTree->Branch("True_Px"     , &True_Px     );
  fArbitraryAnaTree->Branch("True_Py"     , &True_Py     );
  fArbitraryAnaTree->Branch("True_Pz"     , &True_Pz     );
  fArbitraryAnaTree->Branch("True_Dirx"   , &True_Dirx   );
  fArbitraryAnaTree->Branch("True_Diry"   , &True_Diry   );
  fArbitraryAnaTree->Branch("True_Dirz"   , &True_Dirz   );
  fArbitraryAnaTree->Branch("True_Time"   , &True_Time   );

  fArbitraryAnaTree->Branch("True_Geant4_PDG"     , &True_Geant4_PDG     );
  fArbitraryAnaTree->Branch("True_Geant4_ID"      , &True_Geant4_ID      );
  fArbitraryAnaTree->Branch("True_Geant4_MotherID", &True_Geant4_MotherID);
  fArbitraryAnaTree->Branch("True_Geant4_E"       , &True_Geant4_E       );
  fArbitraryAnaTree->Branch("True_Geant4_VertX"   , &True_Geant4_VertX   );
  fArbitraryAnaTree->Branch("True_Geant4_VertY"   , &True_Geant4_VertY   );
  fArbitraryAnaTree->Branch("True_Geant4_VertZ"   , &True_Geant4_VertZ   );
  fArbitraryAnaTree->Branch("True_Geant4_Time"    , &True_Geant4_Time    );
  
}


void ArbitraryAna::analyze(art::Event const & evt) {
  ResetVariables();

  Run    = evt.run();
  SubRun = evt.subRun();
  Event  = evt.event();

  //GET INFORMATION ABOUT THE DETECTOR'S GEOMETRY.
  auto const* geo  = lar::providerFrom<geo::Geometry>();
//  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  art::Handle<std::vector<simb::MCTruth>> GenTrue;
  if (evt.getByLabel(fGenLabel, GenTrue)) {
    art::FindManyP<simb::MCParticle> GenAssn(GenTrue,evt,fGEANTLabel);
    for (size_t igen=0; igen<GenTrue->size(); ++igen) {
      for (int ipart=0; ipart<(*GenTrue)[igen].NParticles(); ++ipart) {
        const simb::MCParticle& part = (*GenTrue)[igen].GetParticle(ipart);
        TrackIDGenerator     [part.TrackId()] = true;
        TrackIDGeneratorIndex[part.TrackId()] = igen;
          
        True_PDG    .push_back(part.PdgCode());
        True_ID     .push_back(part.TrackId());
        True_E      .push_back(part.E());
        True_VertX  .push_back(part.Vx());
        True_VertY  .push_back(part.Vy());
        True_VertZ  .push_back(part.Vz());
        True_VertexT.push_back(part.T());
        True_Px     .push_back(part.Px());
        True_Py     .push_back(part.Py());
        True_Pz     .push_back(part.Pz());
        True_Dirx   .push_back(part.Px() / part.P());
        True_Diry   .push_back(part.Py() / part.P());
        True_Dirz   .push_back(part.Pz() / part.P());
      }
      
      try {
        for (size_t ipart=0; ipart<GenAssn.at(igen).size(); ++ipart) {
          const simb::MCParticle& part = (*GenAssn.at(igen).at(ipart));
          TrackIDGenerator     [part.TrackId()] = true;
          TrackIDGeneratorIndex[part.TrackId()] = igen;
            
          True_Geant4_PDG     .push_back(part.PdgCode());
          True_Geant4_ID      .push_back(part.TrackId());
          True_Geant4_MotherID.push_back(part.Mother());
          True_Geant4_E       .push_back(part.E());
          True_Geant4_VertX   .push_back(part.Vx());
          True_Geant4_VertY   .push_back(part.Vy());
          True_Geant4_VertZ   .push_back(part.Vz());
          True_Geant4_Time    .push_back(part.T());
        }
      } catch (...) {
        //mf::LogWarning("ArbitraryAnaModule");
        std::cout <<"GEANT4 has not been run on the MCTruth " << (*GenTrue)[igen] << std::endl;
      }
    }
  } else {
    mf::LogWarning("ArbitraryAnaModule") << "MC truth cannot be loaded!";
  }
  
  
  if (fSaveTPC) {
    
    art::Handle<std::vector<recob::Hit>> reco_hits;
    
    if(evt.getByLabel(fHitLabel, reco_hits)) {
      
      std::vector<recob::Hit> ColHits_Gen;
      NTotHit = reco_hits->size();
      int LoopHits = NTotHit;

      for(int hit=0; hit<LoopHits; ++hit) {
        recob::Hit const& ThisHit = reco_hits->at(hit);
        if(ThisHit.View() == geo::kU || ThisHit.View() == geo::kV) {
          ++NIndHit;
        } else {
          ++NColHit;
        }
      }

      std::cout<<" Hello!! NColHit "<< NColHit<<"\n";

      auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

      for(int hit=0; hit<LoopHits; ++hit) {

        recob::Hit const& ThisHit = reco_hits->at(hit);

	//std::cout<<" Hello!! View "<< ThisHit.View()<<"\n";
	  
        if (ThisHit.View() == 2) {
          std::vector< sim::TrackIDE>  ThisHitIDE;
          std::vector<const sim::IDE*> ThisSimIDE;
          try {
            // HitToTrackIDEs opens a specific window around the hit. I want a
            // wider one, because the filtering can delay the hit. So this bit
            // is a copy of HitToTrackIDEs from the backtracker, with some
            // modification
            const double start = ThisHit.PeakTime()-20;
            const double end   = ThisHit.PeakTime()+ThisHit.RMS()+20;
            ThisHitIDE = bt_serv->ChannelToTrackIDEs(clockData, ThisHit.Channel(), start, end);
          } catch(...) {
            firstCatch++;
            try {
              ThisSimIDE = bt_serv->HitToSimIDEs_Ps(clockData, ThisHit);
            } catch(...) {
              secondCatch++;
              continue;
            }
            continue;
          }

          // Get the simIDEs.
	  /*
          try {
            ThisSimIDE = bt_serv->HitToSimIDEs_Ps(clockData, ThisHit);
          } catch(...) {
            thirdCatch++;
            continue;
          }
	  */

          Hit_View.push_back(ThisHit.View());
          Hit_Size.push_back(ThisHit.EndTick() - ThisHit.StartTick());
          Hit_TPC .push_back(ThisHit.WireID().TPC);
          int channel = ThisHit.Channel();
          Hit_Chan.push_back(channel);

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

	  std::cout<<" Hello!! "<<ThisHit.WireID().TPC<<" "<< channel <<" "<<ThisHit.SummedADC()<<"\n";

          if (ThisHitIDE.size()==0) NHitNoBT++;

          //WHICH PARTICLE CONTRIBUTED MOST TO THE HIT.
          double TopEFrac = -DBL_MAX;

          Hit_True_TrackID.push_back(0);
          for (size_t ideL=0; ideL < ThisHitIDE.size(); ++ideL) {

            if (ThisHitIDE[ideL].energyFrac > TopEFrac) {
              TopEFrac = ThisHitIDE[ideL].energyFrac;
              Hit_True_TrackID.back() = ThisHitIDE[ideL].trackID;
            }
          }

          Hit_True_GenType .push_back(TrackIDGenerator     [Hit_True_TrackID.back()]);
          Hit_True_GenIndex.push_back(TrackIDGeneratorIndex[Hit_True_TrackID.back()]);
        }
      }
      size_t nbacktracked=0;
      std::cout << fHitLabel << " found " << Hit_True_GenType.size() << " hits (";
      for (auto const& it: Hit_True_GenType) {
        if (it) nbacktracked++;
      }
      std::cout << nbacktracked << " were backtracked to signal)." << std::endl;
    } else {
      mf::LogWarning("ArbitraryAnaModule") << "Requested to save wire hits, but cannot load any wire hits";
    }
  }

  if (fSavePDS) {

    art::Handle< std::vector< recob::OpHit > > OpHitHandle;

    if (evt.getByLabel(fOpHitModuleLabel, OpHitHandle)) {

      for (size_t hit=0; hit<OpHitHandle->size(); ++hit) {
        const recob::OpHit& reco_hit = OpHitHandle->at(hit);
        std::vector<sim::TrackSDP> vec_tracksdp = pbt_serv->OpHitToTrackSDPs(reco_hit);

        std::sort(vec_tracksdp.begin(), vec_tracksdp.end(),
                  [](const sim::TrackSDP& a, const sim::TrackSDP& b) -> bool { return a.energyFrac > b.energyFrac; });

        if (vec_tracksdp.size()>0){
          int MainTrID = vec_tracksdp[0].trackID;
          PDS_OpHit_True_TrackID .push_back(MainTrID);
          PDS_OpHit_True_GenType .push_back(TrackIDGenerator     [MainTrID]);
          PDS_OpHit_True_GenIndex.push_back(TrackIDGeneratorIndex[MainTrID]);
        } else {
          PDS_OpHit_True_TrackID .push_back(-1);
          PDS_OpHit_True_GenType .push_back(false);
          PDS_OpHit_True_GenIndex.push_back(-1);
        }
        
        double xyz_optdet[3]={0,0,0};
        double xyz_world [3]={0,0,0};

        geo->OpDetGeoFromOpChannel(reco_hit.OpChannel()).LocalToWorld(xyz_optdet,xyz_world);
        PDS_OpHit_OpChannel   .push_back(reco_hit.OpChannel());
        PDS_OpHit_X           .push_back(xyz_world[0]);
        PDS_OpHit_Y           .push_back(xyz_world[1]);
        PDS_OpHit_Z           .push_back(xyz_world[2]);
        PDS_OpHit_PeakTimeAbs .push_back(reco_hit.PeakTimeAbs());
        PDS_OpHit_PeakTime    .push_back(reco_hit.PeakTime());
        PDS_OpHit_Frame       .push_back(reco_hit.Frame());
        PDS_OpHit_Width       .push_back(reco_hit.Width());
        PDS_OpHit_Area        .push_back(reco_hit.Area());
        PDS_OpHit_Amplitude   .push_back(reco_hit.Amplitude());
        PDS_OpHit_PE          .push_back(reco_hit.PE());
        PDS_OpHit_FastToTotal .push_back(reco_hit.FastToTotal());
      }
    }
    else {
      mf::LogWarning("ArbitraryAnaModule") << "Requested to save optical hits, but cannot load any ophits";
    }
  }

  fArbitraryAnaTree->Fill();
}

DEFINE_ART_MODULE(ArbitraryAna)
