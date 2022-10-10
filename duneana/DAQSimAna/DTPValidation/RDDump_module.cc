////////////////////////////////////////////////////////////////////////
// Class:       RDDump
// Plugin Type: producer (art v2_10_03)
// File:        RDDump_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"

#include <memory>
#include <fstream>

class RDDump;


class RDDump : public art::EDAnalyzer {
public:
  explicit RDDump(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  RDDump(RDDump const &) = delete;
  RDDump(RDDump &&) = delete;
  RDDump & operator = (RDDump const &) = delete;
  RDDump & operator = (RDDump &&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
    // The module name of the raw digits we're reading in
    std::string m_inputTag;
    std::string m_outputFilename;
    std::ofstream m_outputFile;
};


RDDump::RDDump(fhicl::ParameterSet const & p)
    : EDAnalyzer(p),
      m_inputTag(p.get<std::string>("InputTag", "daq")), 
      m_outputFilename(p.get<std::string>("OutputFile")),
      m_outputFile(m_outputFilename)
{
}

void RDDump::analyze(art::Event const& e)
{
    auto const& digits_handle=e.getValidHandle<std::vector<raw::RawDigit>>(m_inputTag);
    auto& digits_in =*digits_handle;

    for(auto&& digit: digits_in){
        m_outputFile << e.event() << " "
                     << digit.Channel() << " ";
        for(auto const& adc: digit.ADCs()){
            m_outputFile << adc << " ";
        }
        m_outputFile << std::endl;
    }
}

DEFINE_ART_MODULE(RDDump)
