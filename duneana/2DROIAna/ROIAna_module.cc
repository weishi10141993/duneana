#include "ROIAna_module.h"
//Constructor for cnn struct


roiana::ROIAna::ROIAna(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset}  ,
  fWireProducerLabel(pset.get< art::InputTag >("InputWireProducerLabel", "caldata")),
  fRawProducerLabel(pset.get< art::InputTag >("InputRawProducerLabel", "tpcrawdecoder:daq")),
  fSimChannelLabel(pset.get< art::InputTag >("SimChannelLabel", "elecDrift")),
  fSimulationProducerLabel(pset.get< art::InputTag >("SimulationProducerLabel", "largeant"))
{
  fLogLevel           = pset.get<int>("LogLevel", 10);
  fNChanPerApa        = pset.get<int>("ChannelPerApa", 2560);
  fNTicksPerWire      = pset.get<unsigned int>("TickesPerWire", 6000);
  //auto const* geo = lar::providerFrom<geo::Geometry>();
  geo = lar::providerFrom<geo::Geometry>();
  //geo = art::ServiceHandle<geo::Geometry>::get();
  fNPlanes = geo->Nplanes();

  fROI_Peak  = pset.get<float>("ROIPeak", 20);
  fROI_Range = pset.get<int>("ROIRange", 120);
  fROI_CH    = pset.get<int>("ROICH", 10);

  fTreeName          = pset.get<std::string>("TREENAME", "wireana");
  fECMin             = pset.get<float>("fECMin", -1e-8); //settting energy and charge accumulation to start at this value --> ensures true background, i.e 0 energy/0 charge falls into the underflow bin
  fHistEnergyMax     = pset.get<float>("fHistEnergyMax", 1000); 
  fHistChargeMax     = pset.get<float>("fHistChargeMax", 1000); 

}

void roiana::ROIAna::analyze(art::Event const & evt) {


  //get detector property
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

  //get event data
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  std::cout<<"########## EvtNo."<<event<<std::endl;
  std::cout<<"Let's TRY!"<<std::endl;

  MC = !evt.isRealData();
  /// FillTruthInfo to internal data objects
  /// i.e. trkid_to_label_map

  ////////////////////////////////////////////////////////////
  //Build Wire SimChanel List

  art::Handle<std::vector<sim::SimChannel>> simChannelListHandle;
  std::vector<art::Ptr<sim::SimChannel>> channellist;
  if (evt.getByLabel(fSimChannelLabel, simChannelListHandle)) 
    art::fill_ptr_vector(channellist, simChannelListHandle);
  SortWirePtrByChannel( channellist, true );

  ////////////////////////////////////////////////////////////
  //Build RawDigit List
  art::Handle<std::vector<raw::RawDigit>> rawListHandle, noiselessRawListHandle;
  std::vector<art::Ptr<raw::RawDigit>> rawList, noiselessRawList;

  if (evt.getByLabel(fRawProducerLabel, rawListHandle)) art::fill_ptr_vector(rawList, rawListHandle);
  SortWirePtrByChannel( rawList, true );

  // //Get channel-> wire,simchannel map
  if( fLogLevel >= 3 ) std::cout<<"Fill ch_w_sc"<<std::endl;
  for( auto w: rawList ) 
    ch_w_sc[ w->Channel() ].first = w;
  for( auto w: channellist) 
    ch_w_sc[ w->Channel() ].second= w;

  TruthFilter();

  //Creating ROI

}


std::map<int,std::vector<std::pair<int,int>>> roiana::ROIAna::ProcessROI(std::vector<art::Ptr<raw::RawDigit>>& rawList)
{

  std::map<int,std::vector<std::pair<int,int>>> ret;


  int num_channels = rawList.size();
  std::vector<std::vector<float>> ROI_array(num_channels);
  std::vector<std::vector<float>> old_array(num_channels);

  const int adc_size = rawList.front()->NADC();
  std::vector<int> num_store_ch(adc_size,0);
  //int num_store_ch[adc_size] = {};

  int lowest_ch = 0;

  for( auto rawdigit : rawList )
  {
    int channel =  rawdigit->Channel();
    std::vector<float> charges;
    for( auto adc : rawdigit->ADCs() )
    {
      charges.push_back( adc  );
    }
    std::vector<float> newcharge(adc_size,0.0);
    std::vector<float> oldcharge(adc_size,0.0);
    auto chargessort = charges;

    std::nth_element(chargessort.begin(), chargessort.begin()+chargessort.size()/2, chargessort.end());
    float median = chargessort[ chargessort.size()/2];

    //Time ROI
    std::vector<std::pair<int,int>> ranges;
    
    for (int bin = 0; bin < (int) charges.size(); ++bin)
    {
      float central_value = charges[bin]- median;
      oldcharge.at(bin) = central_value;
      // log->debug("RegionOfInterestFilter: carica nel bin {} = {}, median {}, Cvalue {}", bin, charges[bin], median, central_value );  

      if(central_value<-fROI_Peak or central_value>fROI_Peak)
      {

        int rmin = bin-fROI_Range;
        int rmax = bin+fROI_Range-1;
        bool make_newRange = true;
        for( auto it = ranges.begin(); it!= ranges.end(); ++it)
        {
          if (rmin > it->first and rmin < it->second) 
          {
            make_newRange = false;
            if (rmax > it->second) it->second = rmax;
            break;
          }
          if (rmax > it->first and rmax < it->second) 
          {
            make_newRange = false;
            if (rmin < it->first ) it->first= rmin;
            break;
          }
        }
        if (make_newRange) ranges.push_back({rmin,rmax});
        //log->debug("RegionOfInterestFilter: peak in the bin {} = {}, median {}, Cvalue {}, ispeak {}", bin, charges[bin], median, central_value, ispeak(charges[bin]) );
        //

        for(int delta = -fROI_Range; delta < fROI_Range; ++delta)
        {
          int newbin = bin+delta;
          if(newbin>-1 and newbin<(int)charges.size())
            newcharge.at(newbin) = charges[newbin]- median;
        }
      }
    }
    ROI_array[channel-lowest_ch] = newcharge;
    old_array[channel-lowest_ch] = oldcharge;
    ret[channel]=ranges;
  }

  // Channel ROI
  // Tejin Note: I changed ch_ind to begin from 1 to 0
  for (int ch_ind = 0; ch_ind < num_channels; ch_ind++) 
  {
    for(int bin=0; bin<adc_size; bin++)
    {
      // Lower channel ROI
      // peak in channel and zero in channel-1 (start of track)
      if(ispeak(ROI_array[ch_ind].at(bin)) and isZero(ROI_array[ch_ind-1].at(bin)))
      {
        // fill lower channel ROI
        for(int j=-fROI_CH; j<0; j++)
        {
          int update_channel = ch_ind+j;
          if(update_channel>-1)
          {
            ROI_array[update_channel].at(bin) = old_array[update_channel].at(bin);
          }
        }
      }

      // Upper channel ROI
      // zero in channel and peak in channel-1 (end of track)
      if(num_store_ch[bin]==0 and isZero(ROI_array[ch_ind].at(bin)) and ispeak(ROI_array[ch_ind-1].at(bin)))
      {
        // fill upper channel ROI
        ROI_array[ch_ind].at(bin) = old_array[ch_ind].at(bin);
        num_store_ch[bin] = fROI_CH;
      }
      else if (num_store_ch[bin]>0)
        ROI_array[ch_ind].at(bin) = old_array[ch_ind].at(bin);

      // iterate num_store_ch
      if(num_store_ch[bin]>1)
        num_store_ch[bin] -= 1; // >0 means store that many more channels in ROI
      else if(num_store_ch[bin]==1)
        num_store_ch[bin] = -1; // -1 means end of upper ROI
      else if(num_store_ch[bin]==-1)
        num_store_ch[bin] = 0; // 0 means ready to find next end of track
    }
  }
  return ret;
}

void roiana::ROIAna::beginJob() {

  gROOT->SetBatch(1);

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>(fTreeName.c_str() ,fTreeName.c_str() );

  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("MC", &MC);

  std::string name1 = Form("TrueEnergyChargeDeposited_%s",fTreeName.c_str() );
  std::string name2 = Form("TrueEnergyChargeDepositedInROI_%s",fTreeName.c_str() );
  name1 = Form("TrueEnergyDeposited_%s",fTreeName.c_str() );
  name2 = Form("TrueEnergyDepositedInROI_%s",fTreeName.c_str() );
  TrueEnergyDeposited = tfs->make<TH1F>( name1.c_str(), "Energy vs Charge; E (MeV)",100,0,fHistEnergyMax);
  TrueEnergyDepositedInROI =  tfs->make<TH1F>( name2.c_str(), "Energy vs Charge; E (MeV)",100,0,fHistEnergyMax);

  name1 = Form("TrueChargeDeposited_%s",fTreeName.c_str() );
  name2 = Form("TrueChargeDepositedInROI_%s",fTreeName.c_str() );
  TrueChargeDeposited = tfs->make<TH1F>( name1.c_str(), "Charge vs Charge; E (MeV)",100,0,fHistChargeMax);
  TrueChargeDepositedInROI =  tfs->make<TH1F>( name2.c_str(), "Charge vs Charge; E (MeV)",100,0,fHistChargeMax);

}

void roiana::ROIAna::endJob()
{
}

template<class T>
void roiana::ROIAna::SortWirePtrByChannel( std::vector<art::Ptr<T>> &vec, bool increasing )
{
  if( fLogLevel >= 10 ) 
  {
    std::cout<<"Entering SortWirePtrByChannel, sorting "<<vec.size()<<"channels."<<std::endl;
  }
  if (increasing)
  {
    std::sort(vec.begin(), vec.end(), [](art::Ptr<T> &a, art::Ptr<T> &b) { return a->Channel() < b->Channel(); });
  }
  else
  {
    std::sort(vec.begin(), vec.end(), [](art::Ptr<T> &a, art::Ptr<T> &b) { return a->Channel() > b->Channel(); });
  }
}

void roiana::ROIAna::TruthFilter()
{
  for( auto m: ch_w_sc)
  {
    //auto channel = m.first;
    //auto rawdigit = m.second.first;
    auto sim = m.second.second;

    //accumulate energy and charge for all channels
    for( auto &tdcide: sim->TDCIDEMap() )
    {
      std::vector<float> energies(partTypes.size(),fECMin );
      std::vector<float> charges(partTypes.size(),fECMin );
      std::vector<float> energiesNeut(partTypes.size(),fECMin );
      std::vector<float> chargesNeut(partTypes.size(),fECMin );
      std::vector<float> energiesRad(partTypes.size(),fECMin );
      std::vector<float> chargesRad(partTypes.size(),fECMin );

      for( auto &ide: tdcide.second )
      {
        
        if( fLogLevel >= 3 ) std::cout<<"ide.trackID: "<<ide.trackID<<std::endl;
        bool isSignal = trkid_to_label_map[ ide.trackID ] == "NuEScatter" || trkid_to_label_map[ ide.trackID ] == "marley";
        int pdg = PIS->TrackIdToParticle_P( ide.trackID )->PdgCode();
        float energy = ide.energy;
        float numElectrons = ide.numElectrons;
        energies[kAll]+=energy;
        charges[kAll]+=numElectrons;
        int partType = -1;
        if( abs(pdg) == 11 || abs(pdg) == 13 || abs(pdg) == 15 )
        {
          partType = kElectron;
        } else if( abs(pdg) == 2212)
        {
          partType = kProton;
        } else if(abs(pdg) == 2112)
        {
          partType = kNeutron;
        } else if(abs(pdg) == 22)
        {
          partType = kPhoton;
        } else
        {
          partType = kNuc;
        }
        //parsed particle, accumulate energy
        energies[partType]+=energy; charges[partType]+=numElectrons;
        if( isSignal )
        {
          energiesNeut[partType]+=energy; chargesNeut[partType]+=numElectrons;
        }
        else
        {
          energiesRad[partType]+=energy; chargesRad[partType]+=numElectrons;
        }
      }

      if( fLogLevel >= 3 ) std::cout<<"FillHistogram: begin"<<std::endl;
      //Some function to fill histogram
      TrueEnergyDeposited->Fill(energies[0]);
      TrueChargeDeposited->Fill(charges[0]);
      if( fLogLevel >= 3 ) std::cout<<"FillHistogram: end"<<std::endl;
    }
  }

}


DEFINE_ART_MODULE(roiana::ROIAna)
