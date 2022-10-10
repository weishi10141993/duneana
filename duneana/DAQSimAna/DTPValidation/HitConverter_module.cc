////////////////////////////////////////////////////////////////////////
// Class:       HitConverter
// Plugin Type: producer
// File:        RawDigitConverter_module.cc
// Author:      Klaudia Wawrowska (K.Wawrowska@sussex.ac.uk)
//
// Convert firmware TP data from pandas dfs (in csv format) to recob::Hit data objects.
////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>
#include <fstream>
#include <chrono>
#include <functional>
#include <iostream>
#include <string.h>
#include <sstream>
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::View_t, geo::SignalType ... 
#include <memory>


//struct FWHit : recob<Hit>{ }




class HitConverter : public art::EDProducer {

public:

  explicit HitConverter(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  HitConverter(HitConverter const &) = delete;
  HitConverter(HitConverter &&) = delete;
  HitConverter & operator = (HitConverter const &) = delete;
  HitConverter & operator = (HitConverter &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  std::string m_RawDataFrame; 
  std::string m_TPDataFrame;
  std::vector< std::vector<std::string> > ReadData (std::string FPGAInputFileName);

};


HitConverter::HitConverter(fhicl::ParameterSet const & p)
  : EDProducer(p),
    m_RawDataFrame(p.get<std::string>("RawInputFile")),
    m_TPDataFrame(p.get<std::string>("TPInputFile"))
{
  produces< std::vector <recob::Hit> >(); 
}

 
//Read in the file for processing 
std::vector< std::vector<std::string> > HitConverter::ReadData(std::string FPGAInputFileName)
{
  std::ifstream FPGAInputFile;  
  FPGAInputFile.open(FPGAInputFileName);
  
  std::vector< std::vector<std::string> > content; 
  std::vector<std::string> row; 
  std::string line, word; 

  if (FPGAInputFile.is_open() ){
    while (getline(FPGAInputFile, line)){
      row.clear();
      std::stringstream str(line);
      while (getline(str, word, ',')) row.push_back(word);
      content.push_back(row);
    }
  }
  FPGAInputFile.close(); 
  return content; 
}

void HitConverter::produce(art::Event & e)
{

  // --- Little block of code to remove the hit time ambiguity in LArSoft 
  //     due to firmware hit finding being done in terms of timestamps --------
 
  //Raw firmware data 
  std::vector< std::vector<std::string> > raw_content = ReadData(m_RawDataFrame);

  int ticks_per_sample = 32;
  int fir_delay = 16; 
  size_t n_tsamples = raw_content.size() -1 ; 
  //min tstamp in raw data
  std::uint64_t min_ts = std::stoull( raw_content[1][0] );
  //max tstamp in raw data
  std::uint64_t max_ts = std::stoull( raw_content[n_tsamples][0] ); // or min_ts + (n_tsamples * ticks_per_sample.)

  // -------                                                      -------
  
  //Get the firmware hits 
  std::vector< std::vector<std::string> > content = ReadData(m_TPDataFrame);

  size_t n_hits = content.size(); // 1 hit/row

  std::cout << "Converting " << n_hits -1 << " firmware hits to recob::Hit data objects." << std::endl; 

  //initialize vector where we will store our firmware hits in recob::Hit format 
  std::unique_ptr<std::vector <recob::Hit> > firmware_hits( new std::vector<recob::Hit> );
  
  //run over hits (rows)  --> 1st row = column labels so skip 
  for (size_t i = 1; i < n_hits; i++){

    std::uint64_t timestamp = std::stoull(content[i][1]); 

    //only want hits within raw adc span for 1:1 SW-FW comparison 
    if (timestamp > max_ts || timestamp < min_ts) continue;  

    //column structure : 
    // index, ts, offline ch, crate no., slot no., fiber no., wire no., flags, median, accumulator, 
    // start time, end time, peak time, peak adc, hit continue, tp flags, sum adc

    std::uint64_t t_offset = (timestamp - min_ts) / ticks_per_sample; 

    raw::ChannelID_t channel = std::stof(content[i][2]);
    raw::TDCtick_t start_tick = (std::stoi(content[i][10]) - fir_delay) + (int)t_offset;
    raw::TDCtick_t end_tick =  (std::stoi(content[i][11]) - fir_delay) + (int)t_offset; 
    float peak_time =  (std::stoi(content[i][12]) - fir_delay) + (int)t_offset;
    float sigma_peak_time =  0;
    float rms  = 0;
    float peak_amplitude = std::stof(content[i][13]);
    float sigma_peak_amplitude = 0; 
    float summedADC = std::stof(content[i][16]);
    float hit_integral = std::stof(content[i][16]);
    float hit_sigma_integral = 0;
    short int multiplicity = 0; 
    short int local_index = 0; 
    float goodness_of_fit = 0; 
    int dof = 0; //degrees of freedom
    //FIX ME : assuming all hits are collection for now
    geo::View_t view = geo::kUnknown; 
    geo::SigType_t signal_type = geo::kCollection; 
    geo::WireID wireID = geo::WireID(); //(std::strtoul(content[i][5]); // FIX ME

    //Call the recob::Hit constructor 
    recob::Hit thisHit(channel,
		       start_tick,
		       end_tick,
		       peak_time,
		       sigma_peak_time,
		       rms,
		       peak_amplitude,
		       sigma_peak_amplitude,
		       summedADC,
		       hit_integral,
		       hit_sigma_integral,
		       multiplicity,
		       local_index,
		       goodness_of_fit,
		       dof,
		       view,
		       signal_type,
		       wireID); 

    firmware_hits->push_back(thisHit); 
  }
  e.put(std::move(firmware_hits));    
}


DEFINE_ART_MODULE(HitConverter) 
