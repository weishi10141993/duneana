////////////////////////////////////////////////////////////////////////
// Class:       HFChainValidation
// Module Type: analyzer
// File:        HFChainValidation_module.cc
//
// Output waveforms from each HF stage for validation  
////////////////////////////////////////////////////////////////////////


#include <fstream>
#include <algorithm> // std::transform
#include <numeric>   // std::accumulate
// Framework includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
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
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"



class HFChainValidation : public art::EDAnalyzer {

public:

  explicit HFChainValidation(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  HFChainValidation(HFChainValidation const &) = delete;
  HFChainValidation(HFChainValidation &&) = delete;
  HFChainValidation & operator = (HFChainValidation const &) = delete;
  HFChainValidation & operator = (HFChainValidation &&) = delete;

  void analyze(art::Event const & evt) override;

private:

  std::string m_HitLabel;
  std::string m_RawDigitLabel; 

  int          m_nChan; //number of channels to output (files can get pretty big) 
  bool         m_SaveRawDigit; // Save raw digit waveforms for nChannels? 

  //Options associated with pedsub/filtering
  short              m_frugalNContig;
  short              m_ped0; 
  std::vector<short> m_filterTaps;
  short              m_fir_shift;

  bool               m_SavePedSub;  // Save pedestal subtracted waveforms? 
  bool               m_SaveFiltered;  // Save FIR filtered waveforms? 
  bool               m_SaveHits; //Keep trigger primitives produced?


  std::vector<short> findPedestal      (const std::vector<short>& orig, short ncontig, short ped0 );
  std::vector<short> filter            (const std::vector<short>& orig);
  
};

//......................................................
HFChainValidation::HFChainValidation(fhicl::ParameterSet const & p)
  :  EDAnalyzer(p), 
     m_HitLabel           (p.get<std::string>       ("HitLabel")),
     m_RawDigitLabel      (p.get<std::string>       ("RawDigitLabel")),
     m_nChan              (p.get<int>               ("nChanelsToAnalyse",                             10)),
     m_SaveRawDigit       (p.get<bool>              ("SaveRawDigit",                                false)),
     m_frugalNContig      (p.get<short>             ("FrugalPedestalNContig",                          10)),
     m_ped0               (p.get<short>             ("PedInitialValue"      ,                        1000)),
     m_filterTaps         (p.get<std::vector<short>>("FilterCoeffs"         , {2,  9, 23, 31, 23,  9,  2})),
     m_fir_shift          (p.get<short>             ("FirShift",                                          6)),
     m_SavePedSub         (p.get<bool>              ("SavePedSub")),
     m_SaveFiltered       (p.get<bool>              ("SaveFiltered")),
     m_SaveHits           (p.get<bool>              ("SaveHits"))
{
  
}


std::vector<short> HFChainValidation::findPedestal(const std::vector<short>& raw_in, const short ncontig, short m_ped0 ){


  short median = m_ped0; //raw_in[0];                                                           
  std::vector<short> pedestal(raw_in.size(), 0);
  short accumulator=0;
  
  for(size_t i=0; i<raw_in.size(); ++i){
    short sample = raw_in[i]; // The current sample                                                                       
    
    if(sample>median) ++accumulator;
    if(sample<median) --accumulator;

    if(accumulator > ncontig){
      ++median;
      accumulator=0;
    }
    if(accumulator < -1*ncontig){
      --median;
      accumulator=0;
    }
    
    pedestal[i]=median;
  }
  
  return pedestal;
}

std::vector<short> HFChainValidation::filter(const std::vector<short>& pedsub) {
  

  const size_t ntaps = m_filterTaps.size();
  const short*  taps = m_filterTaps.data();
  
  std::vector<short> filtered(pedsub.size(), 0);
  for(size_t i=0; i<pedsub.size(); ++i){

    int temp = 0; //temporary value of filtered adc (to stop it from overflowing)              
                            
    //loop over taps                                                                                                                          
    for(size_t j=0; j<ntaps; ++j){
      const size_t index=i>j ? i-j : 0;
      temp+=pedsub[index]*taps[j];
    }
    //include the FIR bit shift                                                                                                                                                                                    
    filtered[i] = (temp >> m_fir_shift);
  }
  return filtered;
}


void HFChainValidation::analyze(art::Event const & evt){

  std::ofstream PDDump("PedSubWaveform_Dump.txt");
  std::ofstream FirDump("FilteredWaveform_Dump.txt");
  std::ofstream HitDump("TP_Dump.txt");
  
  auto const& digits_handle=evt.getValidHandle<std::vector<raw::RawDigit>>(m_RawDigitLabel);
  auto& digits_in =*digits_handle;
  
  std::vector< std::vector<short>> adc_samples;
  std::vector<unsigned int> channels;

  int counter = 0 ; 
  for(auto&& digit: digits_in){
    counter++; 
    channels.push_back(digit.Channel());
    adc_samples.push_back(digit.ADCs());
    if (counter == m_nChan) break; 
  }


  if (m_SaveRawDigit){ 
    
    std::ofstream RDDump ("RawDigitWaveform_Dump.txt");
    
    for (size_t j = 0; j < channels.size(); j++){
      RDDump << evt.event() << " "
	     << channels[j] << " ";
      for (auto const& adc : adc_samples[j]){
	RDDump << adc << " ";
      }
      RDDump << std::endl; 
    }
    RDDump.close(); 
  }
  
  //Output pedestal-subtracted/ filtered waveforms for checking 
  for (size_t ich=0; ich<adc_samples.size(); ++ich) { 

    const std::vector<short>& waveform = adc_samples[ich];
    std::vector<short> pedestal  = findPedestal(waveform, m_frugalNContig, m_ped0);
    std::vector<short> pedsub(waveform.size(), 0);
    for(size_t i=0; i<pedsub.size(); ++i)   
      pedsub[i]=waveform[i]-pedestal[i];
    std::vector<short> filtered  = filter(pedsub);
    
    
    if (m_SavePedSub){
      PDDump << evt.event() << " "
	     << channels[ich] << " ";
      for (auto const& adc: pedsub){
        PDDump << adc << " ";
      }
      PDDump << std::endl ;    
    }
    
    if (m_SaveFiltered){
      FirDump << evt.event() << " "
	      << channels[ich] << " ";
      for (auto const& adc: filtered){
	FirDump << adc << " ";
      }
      FirDump << std::endl ;    
    }
  }// loop over channels

  if (m_SaveHits) {
    auto reco_hits = evt.getValidHandle<std::vector<recob::Hit> >(m_HitLabel);
 
    for(size_t hit = 0; hit < reco_hits->size(); ++hit) {
      recob::Hit const& ThisHit = reco_hits->at(hit);   // current hit 

      //only ouput hits for the channels of interest 
      for (auto chan : channels){
	if (chan == ThisHit.Channel()){

	  HitDump << ThisHit.StartTick()  << " " << ThisHit.EndTick() - ThisHit.StartTick() 
		  << " " << ThisHit.PeakTime()  << " " << ThisHit.Channel()
		  << " " << ThisHit.SummedADC() << " "<< ThisHit.PeakAmplitude() << std::endl; 
	}
      }
    }
  }
  PDDump.close();
  FirDump.close();
  HitDump.close();
} // Analyze HFChainValidation.


//......................................................
DEFINE_ART_MODULE(HFChainValidation)
