////////////////////////////////////////////////////////////////////////
// Class:       SWTriggerPrimitiveFinder
// File:        SWTriggerPrimitiveFinder_tool.cc
////////////////////////////////////////////////////////////////////////

#include "duneana/DAQSimAna/DTPValidation/SWHitFinding/SWTriggerPrimitiveFinderPass1.h"

#include <algorithm> // for std::transform
#include <numeric> // for std::accumulate

//fcl-configurable parameters
SWTriggerPrimitiveFinderPass1::SWTriggerPrimitiveFinderPass1(fhicl::ParameterSet const & p)
  : m_threshold          (p.get<unsigned int>      ("Threshold"            )),
    m_frugalNContig      (p.get<short>             ("FrugalPedestalNContig",                          10)),
    m_ped0               (p.get<short>             ("PedInitialValue"      ,                        1000)),
    m_filterTaps         (p.get<std::vector<short>>("FilterCoeffs"         , {2,  4,  6,  7,  9,  11,  12,  13,
	                                                                      13,  12, 11,  9,  7,  6,  4,  2})),
    m_fir_shift          (p.get<short>             ("FIRShift"             ,                           6))
{

}



//Pedestal subtraction
std::vector<short> SWTriggerPrimitiveFinderPass1::findPedestal(const std::vector<short>& raw_in,
							       const int ncontig,
							       short m_ped0)
{

  short median = m_ped0; //raw_in[0];                                                                                                                 
  std::vector<short> pedestal(raw_in.size(), 0);
  int accumulator=0;

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


//Filtering
std::vector<short> SWTriggerPrimitiveFinderPass1::filter(const std::vector<short>& pedsub) {

 
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

//Hit Finding
void
SWTriggerPrimitiveFinderPass1::hitFinding(const std::vector<short>& waveform,
				    std::vector<AbsRunningSumTPFinderTool::Hit>& hits,
				    int channel) {


  bool is_hit  = false;
  bool was_hit = false;
  std::vector<int> hit_charge; 
  //initialize the hit 
  AbsRunningSumTPFinderTool::Hit hit(channel, 0, 0, 0, 0);
  
  for(size_t isample=0; isample<waveform.size()-1; ++isample){
    short adc         = waveform[isample];
    //ignore first ~0 ticks to let the pedestal stabilise    
    if (isample > 0) {
      is_hit = adc >  (short)m_threshold;
      if(is_hit && !was_hit) {
	hit_charge.push_back(adc); 
	hit.startTime         = isample;
	hit.SADC              = adc;
	hit.timeOverThreshold = 1;
      }
      if(is_hit && was_hit) {
	hit.SADC              += adc;
	hit.timeOverThreshold += 1;
	hit_charge.push_back(adc);
      }
      if(!is_hit && was_hit) {
	hit.peakCharge = *std::max_element(hit_charge.begin(), hit_charge.end()); 
	hits.push_back(hit);
	hit_charge.clear();
      }
    }
    was_hit = is_hit; 
  } //run over time samples for the waveform
}


std::vector<AbsRunningSumTPFinderTool::Hit>
SWTriggerPrimitiveFinderPass1::findHits(const std::vector<unsigned int>& channel_numbers, 
				  const std::vector<std::vector<short>>& adc_samples) {

  auto hits = std::vector<AbsRunningSumTPFinderTool::Hit>();
  std::cout << "findHits called with "      << adc_samples.size()
	    << " channels. First chan has " << adc_samples[0].size() << " samples" << std::endl;

  for(size_t ich=0; ich<adc_samples.size(); ++ich){
    const std::vector<short>& waveform = adc_samples[ich];
    std::vector<short> pedestal  = findPedestal(waveform, m_frugalNContig ,m_ped0);
    std::vector<short> pedsub(waveform.size(), 0);
    for(size_t i=0; i<pedsub.size(); ++i)
      pedsub[i]=waveform[i]-pedestal[i];
    std::vector<short> filtered  = filter(pedsub);
    hitFinding(filtered, hits, channel_numbers[ich]);
  }
  std::cout << "Returning " << hits.size() << " hits for a threshold of " << m_threshold << std::endl;
  return hits;
}

DEFINE_ART_CLASS_TOOL(SWTriggerPrimitiveFinderPass1)
