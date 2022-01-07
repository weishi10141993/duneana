////////////////////////////////////////////////////////////////////////
//
// \file CAFMaker_module.cc
//
// Chris Marshall's version
// Largely based on historical FDSensOpt/CAFMaker_module.cc
//
///////////////////////////////////////////////////////////////////////

#ifndef CAFMaker_H
#define CAFMaker_H

// Generic C++ includes
#include <iostream>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRGlobal.h"

//#include "Utils/AppInit.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "dunereco/FDSensOpt/FDSensOptData/MVASelectPID.h"
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dunereco/CVN/func/InteractionType.h"
#include "dunereco/CVN/func/Result.h"
#include "dunereco/RegCNN/func/RegCNNResult.h"

// dunerw stuff
#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/exceptions.hh"
//#include "systematicstools/utility/md5.hh"

// root
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

// pdg
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"

// genie
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"

namespace caf {

  class CAFMaker : public art::EDAnalyzer {

    public:

      explicit CAFMaker(fhicl::ParameterSet const& pset);
      virtual ~CAFMaker();
      void beginJob() override;
      void endJob() override;
      void beginSubRun(const art::SubRun& sr) override;
      void endSubRun(const art::SubRun& sr) override;
      void analyze(art::Event const & evt) override;


    private:
      std::string fMVASelectLabel;
      std::string fMVASelectNueLabel;
      std::string fMVASelectNumuLabel;

      std::string fCVNLabel;
      std::string fRegCNNLabel;

      std::string fEnergyRecoNueLabel;
      std::string fEnergyRecoNumuLabel;
      std::string fMVAMethod;

      TFile* fOutFile;
      TTree* fTree;
      TTree* fMetaTree;

      double meta_pot;
      int meta_run, meta_subrun, meta_version;

      systtools::provider_list_t fSystProviders;

  }; // class CAFMaker


  //------------------------------------------------------------------------------
  CAFMaker::CAFMaker(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset), fOutFile(0)
  {
    fMVASelectLabel = pset.get<std::string>("MVASelectLabel");
    fMVASelectNueLabel = pset.get<std::string>("MVASelectNueLabel");
    fMVASelectNumuLabel = pset.get<std::string>("MVASelectNumuLabel");
    fCVNLabel = pset.get<std::string>("CVNLabel");
    fRegCNNLabel = pset.get<std::string>("RegCNNLabel");

    fEnergyRecoNueLabel = pset.get<std::string>("EnergyRecoNueLabel");
    fEnergyRecoNumuLabel = pset.get<std::string>("EnergyRecoNumuLabel");

    // Get DUNErw stuff from its fhicl, which should be included on the CAFMaker config file
    if( !pset.has_key("generated_systematic_provider_configuration") ) {
     std::cout << "[ERROR]: Could not find producer key: "
                  "\"generated_systematic_provider_configuration\". This should "
                  "contain a list of configured systematic providers generated by "
                  "GenerateSystProviderConfig." << std::endl;
     return;
    }

    fhicl::ParameterSet syst_provider_config = pset.get<fhicl::ParameterSet>("generated_systematic_provider_configuration");

    fSystProviders = systtools::ConfigureISystProvidersFromParameterHeaders(syst_provider_config);
  }

  //------------------------------------------------------------------------------
  caf::CAFMaker::~CAFMaker()
  {
  }

  //------------------------------------------------------------------------------
  void CAFMaker::beginJob()
  {
    fOutFile = new TFile("caf.root", "RECREATE");
    fTree = new TTree("cafTree", "cafTree");
    fMetaTree = new TTree("meta", "meta");

    // Create the branch. We will update the address before we write the tree
    caf::StandardRecord* rec = 0;
    fTree->Branch("rec", &rec);

    fMetaTree->Branch("pot", &meta_pot, "pot/D");
    fMetaTree->Branch("run", &meta_run, "run/I");
    fMetaTree->Branch("subrun", &meta_subrun, "subrun/I");
    fMetaTree->Branch("version", &meta_version, "version/I");

    caf::SRGlobal global;

    // make DUNErw variables
    for( auto &sp : fSystProviders ) {
      for(const systtools::SystParamHeader& head: sp->GetSystMetaData()){
        std::cout << "Adding reweight " << head.systParamId << " for " << head.prettyName << " with " << head.paramVariations.size() << " shifts" << std::endl;

        caf::SRSystParamHeader hdr;
        hdr.nshifts = head.paramVariations.size();
        hdr.name = head.prettyName;
        hdr.id = head.systParamId; // TODO is this necessary?

        global.wgts.params.push_back(hdr);
      }
    }

    fOutFile->cd();
    TTree* globalTree = new TTree("globalTree", "globalTree");
    caf::SRGlobal* pglobal = &global;
    TBranch* br = globalTree->Branch("global", "caf::SRGlobal", &pglobal);
    if(!br) abort();
    globalTree->Fill();
    globalTree->Write();

    meta_pot = 0.;
    meta_version = 1;
  }

  //------------------------------------------------------------------------------
  void CAFMaker::beginSubRun(const art::SubRun& sr)
  {
    auto pots = sr.getHandle< sumdata::POTSummary >("generator");
    if( pots ) meta_pot += pots->totpot;
  }

  //------------------------------------------------------------------------------
  void CAFMaker::analyze(art::Event const & evt)
  {
    caf::StandardRecord sr;
    caf::StandardRecord* psr = &sr;
    fTree->SetBranchAddress("rec", &psr);

    auto pidin = evt.getHandle<dunemva::MVASelectPID>(fMVASelectLabel);
    auto pidinnue = evt.getHandle<dunemva::MVASelectPID>(fMVASelectNueLabel);
    auto pidinnumu = evt.getHandle<dunemva::MVASelectPID>(fMVASelectNumuLabel);
    art::InputTag itag1(fCVNLabel, "cvnresult");
    auto cvnin = evt.getHandle<std::vector<cvn::Result>>(itag1);
    art::InputTag itag2(fRegCNNLabel, "regcnnresult");
    auto regcnnin = evt.getHandle<std::vector<cnn::RegCNNResult>>(itag2);
    auto ereconuein = evt.getHandle<dune::EnergyRecoOutput>(fEnergyRecoNueLabel);
    auto ereconumuin = evt.getHandle<dune::EnergyRecoOutput>(fEnergyRecoNumuLabel);

    sr.run = evt.id().run();
    sr.subrun = evt.id().subRun();
    sr.event = evt.id().event();
    meta_run = sr.run;
    meta_subrun = sr.subrun;

    if( !pidin.failedToGet() ) {
      sr.mvaresult = pidin->pid;

      //Fill MVA reco stuff
      sr.Ev_reco_nue     = ereconuein->fNuLorentzVector.E();
      sr.RecoLepEnNue    = ereconuein->fLepLorentzVector.E();
      sr.RecoHadEnNue    = ereconuein->fHadLorentzVector.E();
      sr.RecoMethodNue   = ereconuein->recoMethodUsed;
      sr.Ev_reco_numu    = ereconumuin->fNuLorentzVector.E();
      sr.RecoLepEnNumu   = ereconumuin->fLepLorentzVector.E();
      sr.RecoHadEnNumu   = ereconumuin->fHadLorentzVector.E();
      sr.RecoMethodNumu  = ereconumuin->recoMethodUsed;
      sr.LongestTrackContNumu = ereconumuin->longestTrackContained;
      sr.TrackMomMethodNumu   = ereconumuin->trackMomMethod;
    }

    if( !pidinnue.failedToGet() ) {
      sr.mvanue = pidinnue->pid;
    }

    if( !pidinnumu.failedToGet() ) {
      sr.mvanumu = pidinnumu->pid;
    }

    if( !cvnin.failedToGet() ) {
      //using i = cvn::Interaction;
      //if(cvnin->empty() || (*cvnin)[0].fOutput.size() <= i::kNutauOther){
      if(cvnin->empty()){
        sr.CVNResultIsAntineutrino = sr.CVNResultNue = sr.CVNResultNumu = sr.CVNResultNutau = sr.CVNResultNC = \
        sr.CVNResult0Protons = sr.CVNResult1Protons = sr.CVNResult2Protons = sr.CVNResultNProtons = \
        sr.CVNResult0Pions = sr.CVNResult1Pions = sr.CVNResult2Pions = sr.CVNResultNPions = \
        sr.CVNResult0Pizeros = sr.CVNResult1Pizeros = sr.CVNResult2Pizeros = sr.CVNResultNPizeros = \
        sr.CVNResult0Neutrons = sr.CVNResult1Neutrons = sr.CVNResult2Neutrons = sr.CVNResultNNeutrons = -3;
      }
      else{
        //const std::vector<float>& v = (*cvnin)[0].fOutput;
        //sr.CVNResultNue = v[i::kNueQE] + v[i::kNueRes] + v[i::kNueDIS] + v[i::kNueOther];
        //sr.CVNResultNumu = v[i::kNumuQE] + v[i::kNumuRes] + v[i::kNumuDIS] + v[i::kNumuOther];
        //sr.CVNResultNutau = v[i::kNutauQE] + v[i::kNutauRes] + v[i::kNutauDIS] + v[i::kNutauOther]

        sr.CVNResultIsAntineutrino = (*cvnin)[0].GetIsAntineutrinoProbability();

        sr.cvnnue = (*cvnin)[0].GetNueProbability();
        sr.cvnnumu = (*cvnin)[0].GetNumuProbability();
        sr.cvnnutau = (*cvnin)[0].GetNutauProbability();
        sr.cvnnc = (*cvnin)[0].GetNCProbability();

        sr.CVNResult0Protons = (*cvnin)[0].Get0protonsProbability();
        sr.CVNResult1Protons = (*cvnin)[0].Get1protonsProbability();
        sr.CVNResult2Protons = (*cvnin)[0].Get2protonsProbability();
        sr.CVNResultNProtons = (*cvnin)[0].GetNprotonsProbability();

        sr.CVNResult0Pions = (*cvnin)[0].Get0pionsProbability();
        sr.CVNResult1Pions = (*cvnin)[0].Get1pionsProbability();
        sr.CVNResult2Pions = (*cvnin)[0].Get2pionsProbability();
        sr.CVNResultNPions = (*cvnin)[0].GetNpionsProbability();

        sr.CVNResult0Pizeros = (*cvnin)[0].Get0pizerosProbability();
        sr.CVNResult1Pizeros = (*cvnin)[0].Get1pizerosProbability();
        sr.CVNResult2Pizeros = (*cvnin)[0].Get2pizerosProbability();
        sr.CVNResultNPizeros = (*cvnin)[0].GetNpizerosProbability();

        sr.CVNResult0Neutrons = (*cvnin)[0].Get0neutronsProbability();
        sr.CVNResult1Neutrons = (*cvnin)[0].Get1neutronsProbability();
        sr.CVNResult2Neutrons = (*cvnin)[0].Get2neutronsProbability();
        sr.CVNResultNNeutrons = (*cvnin)[0].GetNneutronsProbability();
      }
    }

    sr.RegCNNNueE = -1.;  // initializing
    if(!regcnnin.failedToGet()){
      if (!regcnnin->empty()){
        const std::vector<float>& v = (*regcnnin)[0].fOutput;
        sr.RegCNNNueE = v[0];
      }
    }

    std::vector< art::Ptr<simb::MCTruth> > truth;
    auto mct = evt.getHandle< std::vector<simb::MCTruth> >("generator");
    if ( mct )
      art::fill_ptr_vector(truth, mct);
    else
      mf::LogWarning("CAFMaker") << "No MCTruth.";

    std::vector< art::Ptr<simb::MCFlux> > flux;
    auto mcf = evt.getHandle< std::vector<simb::MCFlux> >("generator");
    if ( mcf )
      art::fill_ptr_vector(flux, mcf);
    else
      mf::LogWarning("CAFMaker") << "No MCFlux.";
/*
    std::vector< art::Ptr<simb::GTruth> > gtru;
    auto gt = evt.getHandle< std::vector<simb::GTruth> >("generator");
    if ( gt )
      art::fill_ptr_vector(gtru, gt);
    else
      mf::LogWarning("CAFMaker") << "No GTruth.";
*/

    for(size_t i=0; i<truth.size(); i++){

      if(i>1){
        mf::LogWarning("CAFMaker") << "Skipping MC truth index " << i;
        continue;
      }

      sr.isFD   = 1; // always FD
      sr.isFHC  = 999; // don't know how to get this?
      sr.isCC   = !(truth[i]->GetNeutrino().CCNC());  // ccnc is 0=CC 1=NC
      sr.nuPDG  = truth[i]->GetNeutrino().Nu().PdgCode();
      sr.nuPDGunosc = flux[i]->fntype;
      sr.mode   = truth[i]->GetNeutrino().Mode(); //0=QE/El, 1=RES, 2=DIS, 3=Coherent production; this is different than mode in ND
      sr.Ev     = truth[i]->GetNeutrino().Nu().E();
      sr.Q2     = truth[i]->GetNeutrino().QSqr();
      sr.W      = truth[i]->GetNeutrino().W();
      sr.X      = truth[i]->GetNeutrino().X();
      sr.Y      = truth[i]->GetNeutrino().Y();
      sr.NuMomX = truth[i]->GetNeutrino().Nu().Momentum().X();
      sr.NuMomY = truth[i]->GetNeutrino().Nu().Momentum().Y();
      sr.NuMomZ = truth[i]->GetNeutrino().Nu().Momentum().Z();

      sr.vtx_x  = truth[i]->GetNeutrino().Lepton().Vx();
      sr.vtx_y  = truth[i]->GetNeutrino().Lepton().Vy();
      sr.vtx_z  = truth[i]->GetNeutrino().Lepton().Vz();

      //Lepton stuff
      sr.LepPDG   = truth[i]->GetNeutrino().Lepton().PdgCode();
      sr.LepMomX  = truth[i]->GetNeutrino().Lepton().Momentum().X();
      sr.LepMomY  = truth[i]->GetNeutrino().Lepton().Momentum().Y();
      sr.LepMomZ  = truth[i]->GetNeutrino().Lepton().Momentum().Z();
      sr.LepE     = truth[i]->GetNeutrino().Lepton().Momentum().T();
      sr.LepNuAngle = truth[i]->GetNeutrino().Nu().Momentum().Vect().Angle(truth[i]->GetNeutrino().Lepton().Momentum().Vect());

      sr.nP     = 0;
      sr.nN     = 0;
      sr.nipip  = 0;
      sr.nipim  = 0;
      sr.nipi0  = 0;
      sr.nikp   = 0;
      sr.nikm   = 0;
      sr.nik0   = 0;
      sr.niem   = 0;
      sr.niother = 0;
      sr.nNucleus = 0;
      sr.nUNKNOWN = 0;

      sr.eP = 0.;
      sr.eN = 0.;
      sr.ePip = 0.;
      sr.ePim = 0.;
      sr.ePi0 = 0.;
      sr.eOther = 0.;

      for( int p = 0; p < truth[i]->NParticles(); p++ ) {
        if( truth[i]->GetParticle(p).StatusCode() == genie::kIStHadronInTheNucleus ) {

          int pdg = truth[i]->GetParticle(p).PdgCode();
          double ke = truth[i]->GetParticle(p).E() - truth[i]->GetParticle(p).Mass();
          if     ( pdg == genie::kPdgProton ) {
            sr.nP++;
            sr.eP += ke;
          } else if( pdg == genie::kPdgNeutron ) {
            sr.nN++;
            sr.eN += ke;
          } else if( pdg == genie::kPdgPiP ) {
            sr.nipip++;
            sr.ePip += ke;
          } else if( pdg == genie::kPdgPiM ) {
            sr.nipim++;
            sr.ePim += ke;
          } else if( pdg == genie::kPdgPi0 ) {
            sr.nipi0++;
            sr.ePi0 += ke;
          } else if( pdg == genie::kPdgKP ) {
            sr.nikp++;
            sr.eOther += ke;
          } else if( pdg == genie::kPdgKM ) {
            sr.nikm++;
            sr.eOther += ke;
          } else if( pdg == genie::kPdgK0 || pdg == genie::kPdgAntiK0 || pdg == genie::kPdgK0L || pdg == genie::kPdgK0S ) {
            sr.nik0++;
            sr.eOther += ke;
          } else if( pdg == genie::kPdgGamma ) {
            sr.niem++;
            sr.eOther += ke;
          } else if( genie::pdg::IsHadron(pdg) ) {
            sr.niother++; // charm mesons, strange and charm baryons, antibaryons, etc.
            sr.eOther += ke;
          } else if( genie::pdg::IsIon(pdg) ) {
            sr.nNucleus++;
          } else {
            sr.nUNKNOWN++;
          }

        }
      }

      // Reweighting variables

      // Consider
      //systtools::ScrubUnityEventResponses(er);

      sr.total_xsSyst_cv_wgt = 1;

      for(auto &sp : fSystProviders ) {
        std::unique_ptr<systtools::EventAndCVResponse> syst_resp = sp->GetEventVariationAndCVResponse(evt);
        if( !syst_resp ) {
          std::cout << "[ERROR]: Got nullptr systtools::EventResponse from provider "
                    << sp->GetFullyQualifiedName();
          abort();
        }

        // The iteration order here is the same as how we filled SRGlobal, so
        // no need to do any work to make sure they align.
        //
        // NB this will all go wrong if we ever support more than one MCTruth
        // per event.
        for(const std::vector<systtools::VarAndCVResponse>& resp: *syst_resp){
          for(const systtools::VarAndCVResponse& it: resp){
            // Need begin/end to convert double to float
            sr.xsSyst_wgt.emplace_back(it.responses.begin(), it.responses.end());
            sr.cvwgt.push_back(it.CV_response);
            sr.total_xsSyst_cv_wgt *= it.CV_response;
          }
        }
      }
    } // loop through MC truth i

    fTree->Fill();
  }

  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  void CAFMaker::endSubRun(const art::SubRun& sr){
  }

  void CAFMaker::endJob()
  {
    fMetaTree->Fill();

    fOutFile->cd();
    fTree->Write();
    fMetaTree->Write();
    fOutFile->Close();
  }

  DEFINE_ART_MODULE(CAFMaker)

} // namespace caf

#endif // CAFMaker_H
