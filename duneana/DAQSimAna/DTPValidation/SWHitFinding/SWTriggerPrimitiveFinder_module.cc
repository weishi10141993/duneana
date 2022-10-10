////////////////////////////////////////////////////////////////////////
// Class:       SWTriggerPrimitiveFinder
// File:        SWTriggerPrimitiveFinder_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::View_t, geo::SignalType ... 
#include "duneana/DAQSimAna/AbsRunningSumHitFinder/AbsRunningSumTPFinderTool.h" // use absRS Tool since we also want SADC & peak Charge 

#include <memory>

class SWTriggerPrimitiveFinder;


class SWTriggerPrimitiveFinder : public art::EDProducer {
public:
  explicit SWTriggerPrimitiveFinder(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SWTriggerPrimitiveFinder(SWTriggerPrimitiveFinder const &) = delete;
  SWTriggerPrimitiveFinder(SWTriggerPrimitiveFinder &&) = delete;
  SWTriggerPrimitiveFinder & operator = (SWTriggerPrimitiveFinder const &) = delete;
  SWTriggerPrimitiveFinder & operator = (SWTriggerPrimitiveFinder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:
  // The module name of the raw digits we're reading in
  std::string m_inputTag;

  // The actual Service that's doing the trigger primitive finding
  std::unique_ptr<AbsRunningSumTPFinderTool> m_finder;
};


SWTriggerPrimitiveFinder::SWTriggerPrimitiveFinder(fhicl::ParameterSet const & p)
  : EDProducer{p}, 
  m_inputTag(p.get<std::string>("InputTag", "convertedrds")),
  m_finder{art::make_tool<AbsRunningSumTPFinderTool>(p.get<fhicl::ParameterSet>("finder"))}
{
    produces<std::vector<recob::Hit>>();
    produces<art::Assns<raw::RawDigit, recob::Hit>>();
}

void SWTriggerPrimitiveFinder::produce(art::Event & e)
{


  std::vector< std::vector<short>> adc_samples;
  std::vector<unsigned int> channels;
  std::map<raw::ChannelID_t, const raw::RawDigit*> ChanToDigit;
  
  
  auto const& digits_handle=e.getValidHandle<std::vector<raw::RawDigit>>(m_inputTag);
  auto& digits_in =*digits_handle;
  
  for(auto&& digit: digits_in){
    
    ChanToDigit[digit.Channel()]=&digit;
    channels.push_back(digit.Channel());
    adc_samples.push_back(digit.ADCs());
    
  }
  

    // Pass the full list of collection channels to the hit finding algorithm
    std::vector<AbsRunningSumTPFinderTool::Hit> hits_col=m_finder->findHits(channels, adc_samples);
   
    // Loop over the returned trigger primitives and turn them into recob::Hits
    recob::HitCollectionCreator hcol(e, false /* doWireAssns */, true /* doRawDigitAssns */);
    for(auto const& hit : hits_col){

      const raw::RawDigit* digit=ChanToDigit[hit.channel];
      geo::WireID wireID = geo::WireID(); //FIX ME : empty constructor for WIREID since working with unspecified geo *for now* 

      recob::HitCreator lar_hit(*digit,                                //RAW DIGIT REFERENCE.
				wireID,                                       //WIRE ID.
				hit.startTime,                             //START TICK.
				hit.startTime+hit.timeOverThreshold,       //END TICK. 
				hit.timeOverThreshold,                     //RMS.
				hit.startTime + hit.timeOverThreshold*0.5, //PEAK_TIME.
				0,                                         //SIGMA_PEAK_TIME.
				hit.peakCharge,                            //PEAK_AMPLITUDE.
				0,                                         //SIGMA_PEAK_AMPLITUDE.
				hit.SADC,                                  //HIT_INTEGRAL.
				0,                                         //HIT_SIGMA_INTEGRAL.
				hit.SADC,                                  //SUMMED CHARGE. 
				0,                                         //MULTIPLICITY.
				0,                                         //LOCAL_INDEX.
				0,                                         //WIRE ID. (?)
				0                                          //DEGREES OF FREEDOM.
				);
      hcol.emplace_back(std::move(lar_hit), art::Ptr<raw::RawDigit>{digits_handle, 0});
    }
    hcol.put_into(e);
}

DEFINE_ART_MODULE(SWTriggerPrimitiveFinder)
