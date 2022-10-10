#ifndef SWTRIGGERPRIMITIVEFINDERPASS1_H
#define SWTRIGGERPRIMITIVEFINDERPASS1_H

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "duneana/DAQSimAna/AbsRunningSumHitFinder/AbsRunningSumTPFinderTool.h"

class SWTriggerPrimitiveFinderPass1 : public AbsRunningSumTPFinderTool {
 public:
  explicit SWTriggerPrimitiveFinderPass1(fhicl::ParameterSet const & p);
 
  
  virtual std::vector<AbsRunningSumTPFinderTool::Hit> findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& adc_samples);
    
  //protected:
  std::vector<short> findPedestal      (const std::vector<short>& orig, const int ncontig, short m_ped0);
  std::vector<short> filter            (const std::vector<short>& orig);


  void hitFinding(const std::vector<short>& waveform, std::vector<AbsRunningSumTPFinderTool::Hit>& hits, int channel);

  unsigned int       m_threshold;
  short              m_frugalNContig;
  short              m_ped0; 
  std::vector<short> m_filterTaps;
  short              m_fir_shift;
 
};

#endif
