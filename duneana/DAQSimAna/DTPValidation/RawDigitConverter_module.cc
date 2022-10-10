////////////////////////////////////////////////////////////////////////
// Class:       RawDigitConverter
// Plugin Type: producer (art v2_11_03)
// File:        RawDigitConverter_module.cc
// Author:      Klaudia Wawrowska (K.Wawrowska@sussex.ac.uk)
//
// Convert firmware wib data from pandas dfs (in csv format) to raw::RawDigit.
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
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "dunecore/DuneInterface/Service/SimChannelExtractService.h"
#include "dunecore/DuneInterface/Service/AdcDistortionService.h"
#include "dunecore/DuneInterface/Service/AdcCompressService.h"
#include "dunecore/DuneInterface/Service/AdcSuppressService.h"
#include "dunecore/DuneInterface/Service/ChannelNoiseService.h"
#include "dunecore/DuneInterface/Service/PedestalAdditionService.h"
#include <memory>

/*
struct FWRawDigit{ 

  std::vector< raw::RawDigit > RawDigits; 
  std::vector< std::uint_64t > timestamps; 

};
*/

class RawDigitConverter : public art::EDProducer {

public:

  explicit RawDigitConverter(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  RawDigitConverter(RawDigitConverter const &) = delete;
  RawDigitConverter(RawDigitConverter &&) = delete;
  RawDigitConverter & operator = (RawDigitConverter const &) = delete;
  RawDigitConverter & operator = (RawDigitConverter &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:
  std::string m_RawDataFrame;
  //art::ServiceHandle<AdcCompressService> m_pcmp;
  std::vector< std::vector<std::string> > ReadData (std::string FPGAInputFileName);
};


RawDigitConverter::RawDigitConverter(fhicl::ParameterSet const & p)
  : EDProducer(p),
    m_RawDataFrame(p.get<std::string>("RawInputFile"))
{
  produces< std::vector <raw::RawDigit> >(); 
}

 
//Read in the file for processing 
std::vector< std::vector<std::string> > RawDigitConverter::ReadData(std::string FPGAInputFileName)
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

void RawDigitConverter::produce(art::Event & e)
{


  std::vector< std::vector<std::string> > content = ReadData(m_RawDataFrame); 

  //need to accout for the headers/indices 
  size_t n_channels = content[0].size() -1; 
  size_t n_tsamples = content.size()-1; 

  size_t chan_col_index = 0; //channel numbers make up the header of the file
  size_t ts_row_index   = 0; // timestamp numbers = row indices 

  std::vector<std::vector<short>> adc_samples( n_channels, std::vector<short>(n_tsamples,0) );
  std::vector<unsigned int>       channel_numbers;
  std::vector<std::uint64_t>      timestamps; 

  std::cout << "Converting " << n_tsamples << " adc samples for " << n_channels << " channels." << std::endl; 

  //run over time samples/packets (rows)  
  for (size_t i = 0; i < n_tsamples+1; i++){

    //run over channels (columns)            l                                                                        
    for (size_t j = 0; j < n_channels+1; j++){

      // channel numbers stored as headers, timestamps stored as row indices 
      if ( (i == chan_col_index)  && (j >  ts_row_index)) channel_numbers.push_back(std::stoi(content[i][j]));
      if ( (i >  chan_col_index)  && (j == ts_row_index)) timestamps.push_back(std::stoull(content[i][j]));
      if ( (i >  chan_col_index)  && (j >  ts_row_index)) adc_samples[j-1][i-1] = (short)std::stoi(content[i][j]);
      //include a guard in case something overflows 
    }
  }

  //convert the data into raw::RawDigits
  std::unique_ptr<std::vector<raw::RawDigit>> rawDigits( new std::vector<raw::RawDigit>); 
  //pedestal 
  float pedval = 0.0;
  float pedrms = 0.0;  
  ULong64_t samples = n_tsamples; 
  raw::Compress_t compression = raw::kNone; 

  for (size_t ch = 0; ch < n_channels; ch++){
    raw::ChannelID_t channel = channel_numbers[ch];
    raw::RawDigit::ADCvector_t adclist = adc_samples[ch]; 

    raw::RawDigit RD(channel, 
		     samples, 
		     adclist, 
		     compression); 

    RD.SetPedestal(pedval, pedrms); 

    rawDigits->push_back( RD ); 
  }
  e.put(std::move(rawDigits));    
}


DEFINE_ART_MODULE(RawDigitConverter) 
