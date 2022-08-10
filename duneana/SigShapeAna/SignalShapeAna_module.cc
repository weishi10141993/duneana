// SignalShapeAna_module.cc
// written by T. Kosc
// kosc.thomas@gmail.com
//
// Module aimed at analyzing codbox vertical drift geometry data
// Take a root input file with larsoft oriented reco/dataprep objects
// which must include :
//   * recob::Tracks
//   * recob::Hit
//   * recob::Wire
// The module loops through reco track, select some which are "of interest" and then goes in detail
// into the signal at the Wire level.
// Each individual waveform of wires in which there is raw signal associated to the track is written
// The module also adds coherently all individual waveforms to produce a 'mean' signal.
/* ---
 
** assumes vertical drift geometry and drift along X direction.


--- */
//

#ifndef SignalShapeAna_module
#define SignalShapeAna_module

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
///#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
///#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
//#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"
#include "lardata/Utilities/AssociationUtil.h"
///#include "lardata/DetectorInfoServices/LArPropertiesService.h"
///#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
///#include "art/Framework/Principal/Handle.h"
///#include "art/Framework/Principal/View.h"
///#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
///#include "art/Utilities/make_tool.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Assns.h"
///#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "cetlib_except/exception.h"
///#include "nusimdata/SimulationBase/MCParticle.h"


#include "larcore/CoreUtils/ServiceUtil.h"

///#include "tools/IHitEfficiencyHistogramTool.h"
#include "TTree.h"
#include "TH1D.h"

// C++ Includes
#include <map>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <numeric>

// ------------------------------------------------------------------------------------------

 // track_of_interest <double> list of characteristics to recover from tracks
struct track_of_interest {
	size_t run;
	size_t subrun;
	size_t event;
	size_t plane;
	std::array<std::unordered_map<size_t, double>, 3 > wiremin;
	std::array<std::unordered_map<size_t, double>, 3 > wiremax;
	size_t tickmin;
	size_t tickmax;
	size_t Npoints;
	double theta;
	double phi;
	double length;
      }; 

// ------------------------------------------------------------------------------------------

bool Tolerance(double a, double b)
  {
  double diff = std::abs(a-b);
  return diff < 0.1;
  } 

// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------

namespace SignalShapeAna
{
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition
  class SignalShapeAna : public art::EDAnalyzer
  {
  public:
    // Standard constructor and destructor for an ART module.
    explicit SignalShapeAna(fhicl::ParameterSet const& pset);
    virtual ~SignalShapeAna();
    // This method is called once, at the start of the job. 
    void beginJob();
    void endJob();
    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    void beginRun(const art::Run& run);
    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);
    // The analysis routine, called once per event.
    void analyze (const art::Event& evt);

private:

    /* ---
 	Method that calculates the spherical coordinates of the track, assuming X as the vertical coordinate
 	Params:
		* recob::Track : larsoft track object
	Returns :
		* tuple<double>(theta, phi) 
    --- */ 
    std::tuple<double, double> GetThetaPhiAngles(recob::Track);
    std::tuple<double, double> GetThetaPhiAngles(art::Ptr<recob::Track> track){return GetThetaPhiAngles(*track);};
  
    /* ---
    Decides whether a track should be tagged of "quality"
    In practival, takes the length in X coordinate and compare it to a input parameter
 	Params:
		* recob::Track : larsoft track object
	Returns bool
    --- */ 
    bool TrackQualityCut(recob::Track track);


    /* ---
    Finds the matching index of the track (1st argument) in the list of track of interests (second argument)
    If no match is found, the function returns -1.
    In practival, takes the length in X coordinate and compare it to a input parameter
 	Params:
		* recob::Track : larsoft track object
		* track_of_interest hand-made struct that holds the main characteristics of the tracks needed in this analyssis
	Returns matching index of the track in the vector of track_of_interest. If no match is found, returns -1
    --- */ 
    int FindMatchingTrackOfInterest(art::Ptr<recob::Track>, std::vector<track_of_interest>);


    /* ---
    Add description of the function
    Updates the element track_of_interest depending on the input hit info
    argument : 
       * m_hit : art::Ptr<recob::Hit> larsoft object
       * track_of_interest object to be updated by the function. 
    no return
 * ---- */
    void UpdateTrackOfInterest(art::Ptr<recob::Hit>, track_of_interest&);

    /* ---
    Returns a triplet of (planeID, volTPC, wireID)
      Params :
         * str_wid is the char converted to a string returned by the method recob::Hit->WireID() 
      Return triplet
    Caveats : str_wid  also contains the cryostat ID, which is ignored here.
    ...
 * ---- */
   std::tuple<size_t, size_t, size_t> GetWireVolTPCViewandID(std::string str_wid);


    /* ---
    For a given track_of_interest object, digs into recob::Wire objects and store signal wires into histograms
    Histos X axis ranges in tick units from track_of_interest.tickmin to track_of_interest.tickmax with a bit of extension
    to make sure all waveforms are visible.
    Waveforms and centered waveforms are written.
      Params:
	* track_of_interest : strcuture contaning relevant information to describe a track to study
	* std::vector<recob::Wire> larsoft Wire objjects.
	* directory to write output objects inside root file 
      Returns an array of 3 vectors of TH1F which are the centered waveforms.
      1st element of array = vector of TH1 for induction 1, etc.
    ...
 * ---- */
   std::array<std::vector<TH1F*>, 3> StoreWireRegionsOfTracksOfInterest(track_of_interest &, std::vector<recob::Wire>, std::string);

    /* ---
    Write the content of input recob::Wire object associated to track_of_oiniterest into a histogram, and updates vector of TH1F (centered waveform)
    adding the current waveform.
      Params :
	* track_of_interest : strcuture contaning relevant information to describe a track to study
	* size_t : voltpc index
	* recob::Wire object (larsoft)
	* std::vector<TH1F*> : vector of cenetred waveform which will be updated by the method, adding the current waveform at play
	* art::TFileDirectory : directory to write histogram in current root file
      No return
 * ---- */
   void StoreWireRegionOfInterestAsHistogram(track_of_interest &, size_t, recob::Wire, std::vector<TH1F*>&, art::TFileDirectory);


    /* ---
    Add together all stored waveforms associated to track_of_interest
    Resulting histogram is stored
      Params:
	* track_of_interest : strcuture contaning relevant information to describe a track to study
	* std::array<std::vector<TH1F*>, 3> : array of vectors of centered waveforms :
		--> 1st array element is vector of centered waveforms for induction 1
		--> 2nd array element is vector of centered waveforms for induction 2
		--> 3rd array element is vector of centered waveforms for collection
	The method will check that the binning is the same for all the hists
	* std::string : directory in root file to write objects
	  if not, exception raised.
      No return
    ...
 * ---- */
  void PerformCoherentAdditionOfWaveforms(track_of_interest, std::array<std::vector<TH1F*>, 3>, std::string);


    // The variables that will go into the n-tuple.
    size_t fNumEvents; // local event counter, i.e 1st event analyzed = 1
    size_t fEvent;   // global event number
    size_t fRun;   // run number
    size_t fSubRun; // subrun number


    // parameter to read from fcl file   

    bool fVerbose; // boolean for the module to be a little talktative.
    std::string fPandoraTrackLabel; // the label of recob::Track objects
    std::string fTrackHitAssnsLabel; // the label of art::Assns<recob::Track,recob::Hit,void>
    std::string fWireLabel; // label associated to recob::Wire objects, i.e rawdigits with coherent noise removed
    double fTrackMinXExtension; // length in cm, projected along X ccoordinate, any track should have to be considered in the analysis.
    std::vector<double> fThetaRange; // [deg] theta_min and theta_max (X is the referent axis). Tracks are kept if theta_min < theta < theta_max
    size_t fPlane; // 0 = ind1, 1 = ind2, and 2 = col 

    // Some internal variables needed globally
    int mRiseWindowSize = 25; // for waveforms signals : time window desired in tick unit before the max of the waveform 
    int mDescendWindowSize = 50; // for waveforms signals : time window desired in tick unit after the max of the waveform
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;
    std::unordered_map<size_t, size_t> NtrackperEvent; // mapping of type (ievent, Ntracks) to keep track of event number in which track(s) have been saved, and how many have. 

 
    // Keep track of the angles of the tracks 
    // theta, phi angles of the traxks, spherical coordinate assuming X as referent axis (instead of z in canonical sph. coord.)
    TH1D* hTheta;
    TH1D* hPhi; 
    
    // Might need the geometry at some point
    const geo::GeometryCore*           fGeometry;       // pointer to Geometry service

  }; // end class SignalShapeAna def


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  //-------------******--------------*******-------------------------------
  //----------**--^^^--***---------*** ^^^-***-----------------------------
  //-----------^^^---^^^-----------^^^---^^^-------------------------------
  //---------^^^-------^^^-------^^^-------^^^-----------------------------
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  //---------///----------------------////---------------------------------
  //----------/////-------------//////-------------------------------------
  //---------------/////////////-------------------------------------------
  //-----------------------------------------------------------------------



  // class implementation
  //-----------------------------------------------------------------------
  // Constructor
  SignalShapeAna::SignalShapeAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {

   fGeometry = lar::providerFrom<geo::Geometry>();
   std::string strGDML = fGeometry->GDMLFile();
   // add an exception if someone tries to run this module with
   // a geometry different from vertical drift coldbox CRP1
   if (strGDML.find("vdcb1") == strGDML.npos)
   {
   throw cet::exception("SignalShapeAna") << "You're apparently trying to run this module with a geometry different from VD coldbox crp1. It is a risky attempt since this module will surely fail for any other geometry." << std::endl;
   }

    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------

  // Destructor
  SignalShapeAna::~SignalShapeAna()
  {}

  //-----------------------------------------------------------------------
  void SignalShapeAna::beginJob()
  {

    // store angle information
    art::TFileDirectory dir = tfs->mkdir("");
    hTheta = dir.make<TH1D>("thetaX", "thetaX", 91, -0.5, 180.5);
    hPhi = dir.make<TH1D>("phiX", "phiX", 91, -180.5, 180.5);

    // zero out the event counter
    fNumEvents = 0;
  } // end method beginjob()

  //-----------------------------------------------------------------------

  void SignalShapeAna::beginRun(const art::Run& /*run*/)
  {
  }
  //-----------------------------------------------------------------------
  void SignalShapeAna::reconfigure(fhicl::ParameterSet const& p)
  {
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    // Implement the tools for handling the responses
//    const std::vector<fhicl::ParameterSet>& hitHistogramToolVec = p.get<std::vector<fhicl::ParameterSet>>("HitEfficiencyHistogramToolList");


fhicl::Sequence<double> _ThetaRange  { fhicl::Name("ThetaRange" ), fhicl::Comment("") };


    fVerbose		 = p.get<bool>("Verbose");	
    fPandoraTrackLabel	 = p.get<std::string>("PandoraTrackLabel");
    fTrackHitAssnsLabel	 = p.get<std::string>("TrackHitAssnsLabel");
    fWireLabel 		 = p.get<std::string>("WireLabel");
    fTrackMinXExtension	 = p.get<double>("TrackMinXExtension");
    fThetaRange		 = p.get<std::vector<double> >("ThetaRange");
    fPlane		 = p.get<size_t>("Plane");

    if (fThetaRange.size() != 2) throw cet::exception("SignalShapeAna") << "fThetaRange range must be 2 (thetamin, thetamax). Current size = " << fThetaRange.size() << std::endl;
    if (fPlane != 0 && fPlane != 1 && fPlane != 2) throw cet::exception("SignalShapeAna") << "Plane must be 0, 1 or 2. Current value = " << fPlane << std::endl;
//    if (fPlane != 2) throw cet::exception("SignalShapeAna") << "module not ready yet to study induction type signals." << std::endl;

    return;
  }


  //-----------------------------------------------------------------------


  void SignalShapeAna::analyze(const art::Event& event)
  {
    // init a vector of track_of_interest objects
    std::vector<track_of_interest> fvec_Toi;

    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event();
    fRun    = event.run();
    fSubRun = event.subRun();
    fNumEvents++;



   //////////////////////////////////////////////////////////////////////////
   //////////////// 1st part of the analysis ////////////////////////////////
   //////////////////////////////////////////////////////////////////////////
   // Loop over recob::tracks objects and store basic features of interest useful
   // for the follwing

    art::InputTag tag_pandoraTrack{fPandoraTrackLabel};

   // Retrieve recob::Track and art::Assns objects
   auto const& pandoraTrack = *event.getValidHandle< std::vector<recob::Track> >(tag_pandoraTrack);

   // Get if data or sim
   // I think for the moment the module is fine, but I keep track of the call in case I need it at some point
///   bool isMC = !event.isRealData();

//   if (fVerbose) std::cout << "Ntracks to study = " << pandoraTrack.size() << std::endl;
  
   // Loop over tracks and select interesting ones
   for (size_t itrack=0; itrack<pandoraTrack.size(); itrack++){		

	// track must pass a quality cut
	if (!TrackQualityCut(pandoraTrack[itrack])) continue;
	

	// Calculate spherical angles using X as the vertical axis
        std::tuple<double, double> trackSphericalX = GetThetaPhiAngles(pandoraTrack[itrack]);	
        double thetaX = std::get<0>(trackSphericalX);
        double phiX = std::get<1>(trackSphericalX);
	hTheta->Fill(thetaX);
	hPhi->Fill(phiX);

	// track reco thetaX angle must lie in range thetamin --> thetamax
	if ( thetaX < fThetaRange[0] || thetaX > fThetaRange[1] ) continue;

      // a Track of interest is defined according to its number of points and the reco angle values
      // assume that 2 tracks will unlikely have these same characteristics in a given event
      // init easy characteristics
      // more complicated characteristics will be updated later
  
        track_of_interest ToI;
	ToI.run 	= fRun;
	ToI.subrun 	= fSubRun;
	ToI.event	= fEvent;
	ToI.Npoints	= pandoraTrack[itrack].NPoints();
	ToI.theta	= thetaX;
	ToI.phi		= phiX;
	ToI.length	= pandoraTrack[itrack].Length();
        fvec_Toi.push_back(ToI);

/*        std::vector<track_of_interest> ToI(3);
	for (size_t i = 0; i < 3; i++)
	  {
	  ToI[i].run	 	= fRun;
	  ToI[i].subrun 	= fSubRun;
	  ToI[i].event		= fEvent;
	  ToI[i].Npoints	= pandoraTrack[itrack].NPoints();
	  ToI[i].theta		= thetaX;
	  ToI[i].phi		= phiX;
	  ToI[i].length		= pandoraTrack[itrack].Length();
	  }
        fvec_Toi.push_back(ToI);
*/
	if (fVerbose) std::cout << "Adding track with NPoints = " << ToI.Npoints << " & thetaX = " << ToI.theta << "° & phiX " << ToI.phi << "° & length = " << ToI.length << " cm." << std::endl;

      } // end for over pandoraTrack
   

   if (fvec_Toi.size() < 1){ std::cout << "No track of interest found" << std::endl; return; }





 
   //////////////////////////////////////////////////////////////////////////
   //////////////// 2nd part of the analysis ////////////////////////////////
   //////////////////////////////////////////////////////////////////////////

   // Loop over art::Assns<recob::Track,recob::Hit> to update the track of interest
   // features info : wireID domain, voltpc, time region (in tick units) 

    // Retrieve art::Assns<recob::Track,recob::Hit>  object
    art::InputTag tag_TrackHitAssns{fTrackHitAssnsLabel};
    auto const& artTrackHit = *event.getValidHandle< art::Assns<recob::Track,recob::Hit,void> >(tag_TrackHitAssns);

	// now loop over tracks - hits associations to define the track of interest more precisely
	// including wire range and tick range 
	for (auto it = artTrackHit.begin(); it != artTrackHit.end(); ++it)
		{
		size_t i = std::distance(artTrackHit.begin(), it);
		auto m_track = std::get<0>(artTrackHit[i]);
		int i_Toi = FindMatchingTrackOfInterest(m_track, fvec_Toi);

// Here I should find a way to speed up the finding of the matching. 
// Typically atrTrackHit is sorted, so I should make the best out of that.

		if (i_Toi < 0) continue;
		if (i_Toi > static_cast<int>((fvec_Toi.size()-1)) ) throw cet::exception("SignalShapeAna") << "matching index of ToI vector greater than size of track_of_interest vector" << std::endl;;

		auto m_hit = std::get<1>(artTrackHit[i]);	
		// updates the ToI informations
		UpdateTrackOfInterest(m_hit, fvec_Toi[i_Toi]);

		} // end for over artTrackHit




	// check wiremin and wiremax members of track_of_interest must have same key values
	// WORK : add this in a function and perform a test, it'll be more explicit
	for (size_t k = 0; k < fvec_Toi.size(); k++)
	  {
	  // loop over planes 0, 1 and 2
	  for (size_t iplane = 0; iplane < 3; iplane++)
	    {
	    for (auto & it : fvec_Toi.at(k).wiremin.at(iplane))
	      { 
	      size_t vtpc = it.first; 
	      if (fvec_Toi.at(k).wiremax.at(iplane).find(vtpc) == fvec_Toi.at(k).wiremax.at(iplane).end()) throw cet::exception("SignalShapeAna") << "found a voltpc index in wiremin map but not in wiremax map" << std::endl;
	      } // end for over wiremin
	    for (auto & it : fvec_Toi.at(k).wiremax.at(iplane)) // switch roles of wiremin and wiremax
	      { 
	      size_t vtpc = it.first; 
	      if (fvec_Toi.at(k).wiremin.at(iplane).find(vtpc) == fvec_Toi.at(k).wiremin.at(iplane).end()) throw cet::exception("SignalShapeAna") << "found a voltpc index in wiremax map but not in wiremin map" << std::endl;
	      } // end for over wiremax
	    } // end loop over planes
	  } // end check loop


	// make some printing regarding the track of interest found int the current event
	if (fVerbose){
	  // check track of interest
	  for (size_t k = 0; k < fvec_Toi.size(); k++)
	    {
	    std::cout << "--- track of interest #" << k+1 << " --- \n";
	    std::cout << "\t--> run " << fvec_Toi.at(k).run << "; subrun " << fvec_Toi.at(k).subrun << "; event " << fvec_Toi.at(k).event << "\n"; 
	    std::cout << "\twire domain : \n";
	    for (size_t i_plane = 0; i_plane < 3; i_plane++)
	      {
	      std::cout << "\tview plane #" << i_plane << "\n";
	      for (auto& it_wmin: fvec_Toi.at(k).wiremin.at(i_plane)){
	        size_t vtpc = it_wmin.first;
	        std::cout << "\t\t\tvtpc " << vtpc;
	        std::cout << " w in [ " << it_wmin.second << " --> " << fvec_Toi.at(k).wiremax.at(i_plane)[vtpc] << " ]\n";
      	        }
	      } // end loop over view planes
	    std::cout << "\t tick domain : " << fvec_Toi.at(k).tickmin << " --> " << fvec_Toi.at(k).tickmax << "\n";
//	    std::cout << "\t reco info : thetaX = " << fvec_Toi.at(k).theta << "° ; phiX = " << fvec_Toi.at(k).phi << "° ; length = " << fvec_Toi.at(k).length << "\n";
	    std::cout << "-----------------------" << std::endl;
	  } // end for
	} // end if Verbose



   	//////////////////////////////////////////////////////////////////////////
	//////////////// 3rd part of the analysis ////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	// Retrieve recob::Wire objects in the {vtpc, wireID, ticks} regions
	// associated to tracks of interest
	// and store the wire signals in histograms
	
        art::InputTag tag_Wire{fWireLabel};
	auto const& Wires       = *event.getValidHandle< std::vector<recob::Wire> >(tag_Wire);


	// Startr loop over tracks of interest
	for (size_t kToi = 0; kToi < fvec_Toi.size(); kToi++)
	  {
	  // skip track of interest too close from time window boundaries
	  if (fvec_Toi[kToi].tickmin < 50 || fvec_Toi[kToi].tickmax > (Wires.at(0).NSignal() - 50) ) continue;

    	  // init the directory direction to store the wires
//    	  size_t mEvent = fvec_Toi.at(kToi).event; 
      	  std::string strdir = "event" + std::to_string(fEvent);
    	  if (NtrackperEvent.find(fEvent) != NtrackperEvent.end()) strdir += "_" + std::to_string(NtrackperEvent.at(fEvent));
	  NtrackperEvent[fEvent]++;

	  std::array<std::vector<TH1F*>, 3> mCenteredWaveforms = StoreWireRegionsOfTracksOfInterest(fvec_Toi[kToi], Wires, strdir);

	if (fVerbose) std::cout << "--- track of interest #" << kToi+1 << " --- \n";
	  PerformCoherentAdditionOfWaveforms(fvec_Toi[kToi], mCenteredWaveforms, strdir);
	if (fVerbose) std::cout << "-----------------------" << std::endl;
	  } // end for over tracks of interest fvec_Toi




    return;

  } // end method analyze

/*----------------------------------------------------------------------------------------------------*/

  void SignalShapeAna::endJob()
  {
    std::cout << "endjob running" << std::endl;
    return;
  }

/*-----------------------------------------------------------------------------------------------------*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------------------------------------ Non art methods ------------------------------------------------
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::tuple<double, double> SignalShapeAna::GetThetaPhiAngles(recob::Track track)
  {
  // get default theta and phi values of the track
  double thetaZ = track.Theta();
  double phiZ = track.Phi();
  // get start and end coordinates of track
  double startX = track.Start().X();
  double endX = track.End().X();
  // apply correction if track is reconstructed as travelling from bottom to top (ie low X to great X) 
  if ( (startX-endX)>0 ){ thetaZ = TMath::Pi() - thetaZ; phiZ = phiZ+TMath::Pi(); }

  // build track direction in cartesian system
  std::vector<double> trackDir = { TMath::Sin(thetaZ)*TMath::Cos(phiZ), TMath::Sin(thetaZ)*TMath::Sin(phiZ), TMath::Cos(thetaZ) };
  // calculate spherical angles taking X as the referent axis
  double invpi = 1. / TMath::Pi();
  double phiX = TMath::ATan2(trackDir[2], trackDir[1]) * 180. * invpi;
  double norm = TMath::Sqrt(trackDir[0]*trackDir[0]+trackDir[1]*trackDir[1]+trackDir[2]*trackDir[2]);
  double thetaX = TMath::ACos(trackDir[0] / norm ) * 180. * invpi;

  // checks
  //std::cout << "TrackDir = { " << trackDir[0] << " ; " << trackDir[1] << " ; " << trackDir[2] << " } " << std::endl; 
  //std::cout << "(thetaX ; phiX) = ( " << thetaX << " ; " << phiX << " )." << std::endl;

  return std::make_tuple(thetaX, phiX);

}	

// --------------------------------------------------------------------------------------------------------------------------

bool SignalShapeAna::TrackQualityCut(recob::Track track)
{
  double Lx = TMath::Abs(track.Start().X() - track.End().X());
  if (Lx < fTrackMinXExtension) return false;
  return true;
}

// --------------------------------------------------------------------------------------------------------------------------

int SignalShapeAna::FindMatchingTrackOfInterest(art::Ptr<recob::Track> Track, std::vector<track_of_interest>  ToIs){
  /*
 * 	find corresponding index of the track myTrack in the vector myToIs
 * 	Returns index position if found a match.
 * 	returns -1 if not 
 * 	*/
	int imatch = -1;
	for (size_t i=0; i<ToIs.size(); i++)
		{
		// calculate theta and phi associated to the current track
		std::tuple<double, double> trackSphericalX = GetThetaPhiAngles(Track);	
		double thetaX = std::get<0>(trackSphericalX);
		double phiX = std::get<1>(trackSphericalX);
		//std::cout << "Current track info : " << myTrack->NPoints() << "  " << thetaX << " " << phiX << std::endl;
		if ( Track->NPoints() == ToIs[i].Npoints && Tolerance(thetaX, ToIs[i].theta) && Tolerance(phiX, ToIs[i].phi) )	
			{
		 	imatch = i;
		 	break; 
		 	} // end if
		 } // end for
		
   return imatch;	

  }


// --------------------------------------------------------------------------------------------------------------------------

 std::tuple<size_t, size_t, size_t> SignalShapeAna::GetWireVolTPCViewandID(std::string str_wid)
 {
   // the str_wid must be formatted as C:0 T:2 P:1 W:95 (numbers are just for illustration)
   // otherwise the various erase functions below will crash, and we don't want that.
   if ( str_wid.find("C:")==str_wid.npos || str_wid.find("T:")==str_wid.npos || str_wid.find("P:")==str_wid.npos || str_wid.find("W:")==str_wid.npos)
	{
	throw cet::exception("SignalShapeAna") << "Method GetWireVolTPCViewandID(str_wid) probably called out of its application domain, please check." << std::endl;
	}

   // Get voltpc and turn it to an integer 
   char c_voltpc =  str_wid.at(str_wid.find("T:")+2);
   std::string str_voltpc; str_voltpc.push_back(c_voltpc);
   size_t i_voltpc = std::stoi(str_voltpc);

   // same operatn with plane
   char c_plane =  str_wid.at(str_wid.find("P:")+2);
   std::string str_plane; str_plane.push_back(c_plane);
   size_t i_plane = std::stoi(str_plane);

   // same operatn with wireID
   std::string str_wire =  str_wid.erase(0, str_wid.find("W:")+2);
   size_t i_wire = std::stoi(str_wire);
 
 return std::make_tuple(i_plane, i_voltpc, i_wire);
 }

// --------------------------------------------------------------------------------------------------------------------------



void SignalShapeAna::UpdateTrackOfInterest(art::Ptr<recob::Hit> m_hit, track_of_interest & _ToI)
  {
  /*
  Updates the element track_of_interest depending on the input hit info
  argument : 
  -- m_hit : art::Ptr<recob::Hit> larsoft object
  -- _ToI, : track_of_interest object to be updated
  */

	// store current hit info in variables
	std::tuple<size_t, size_t, size_t> hitIDInfo = GetWireVolTPCViewandID( (std::string) m_hit->WireID());
	size_t plane = std::get<0>(hitIDInfo);
///	if (plane != fPlane) return;
	size_t vtpc = std::get<1>(hitIDInfo);
	size_t wireID = std::get<2>(hitIDInfo);
	float peaktime = m_hit->PeakTime();

if (plane != 0 && plane !=1 && plane !=2) throw cet::exception("SignalShapeAna") << "Found a plane ID not in 0, 1 or 2." << std::endl;
 
//std::cout << "Current hit info (vtpc plane wireID tick ) : " << " " << vtpc << " " << plane << " " << wireID << " " << m_hit->PeakTime() << std::endl;
	
	// check whether vecROIs[i_ROI] has already been updated
	bool isFirst = false;
	if (_ToI.wiremin.at(plane).size() == 0) isFirst = true;

	// if _ToI has never been updated, simply init its wire/time infos
	if (isFirst) 
		{
		if (_ToI.wiremin.at(plane).size() > 0 || _ToI.wiremax.at(plane).size()) {
		  throw cet::exception("SignalShapeAna::UpdateTrackOfInterest") << "Vec member wiremin or wiremax of struct track_of_interest has size > 0 before first element added" << std::endl;
		  }
		_ToI.plane = plane; // obsolete
		_ToI.wiremin.at(plane)[vtpc] = wireID;
		_ToI.wiremax.at(plane)[vtpc] = wireID;
		_ToI.tickmin = static_cast<size_t>(peaktime+0.5); // approx float to closest integer
		_ToI.tickmax = static_cast<size_t>(peaktime+0.5);
		}
	else 
		{
		// update tick time
		if (_ToI.tickmin > peaktime)	_ToI.tickmin = static_cast<size_t>(peaktime+0.5);;
		if (_ToI.tickmax < peaktime)	_ToI.tickmax = static_cast<size_t>(peaktime+0.5);

		// update wireticks
		// check if any wire has already been added in the current vtpc
		auto it = _ToI.wiremin.at(plane).find(vtpc);
		bool isVtpcExist = (it != _ToI.wiremin.at(plane).end());

		// if current vtpc exists, update its extremas
		if (isVtpcExist)
		   {
		   if (_ToI.wiremin.at(plane)[vtpc] > wireID) _ToI.wiremin.at(plane)[vtpc] = wireID;
		   if (_ToI.wiremax.at(plane)[vtpc] < wireID) _ToI.wiremax.at(plane)[vtpc] = wireID;
		   }
		else // simply add the maps corresponding to key value vtpc
		   {
		   _ToI.wiremin.at(plane)[vtpc] = wireID;
		   _ToI.wiremax.at(plane)[vtpc] = wireID;
		   }

		} // end else (if isFirst)



  return;


  } // end method UpdateTrackofInterest

// --------------------------------------------------------------------------------------------------------------------------


// --------------------------------------------------------------------------------------------------------------------------


std::array<std::vector<TH1F*>, 3> SignalShapeAna::StoreWireRegionsOfTracksOfInterest(track_of_interest & mToI, std::vector<recob::Wire> mWires, std::string mainwritedir)
  {

  // update tick window to make it a bit larger
  mToI.tickmin -= 100.;
  if (mToI.tickmin < 0) mToI.tickmin = 0;
  mToI.tickmax += 100.;
  if (mToI.tickmax > (float) mWires.at(0).NSignal()) mToI.tickmax = mWires.at(0).NSignal() - 1; 

  // output vector declaration
  std::array<std::vector<TH1F*>, 3> arrayvecCenteredWaveforms; // vec of histograms each containing centered waveforms


  // loop over plane views
  for (int plane = 0; plane < 3; plane++)
    {
    std::string writedir = mainwritedir;
    if (plane == 0) writedir += "/ind1";
    else if (plane == 1) writedir += "/ind2";
    else if (plane == 2) writedir += "/col";
    else writedir += "/unknown";
    art::TFileDirectory dir = tfs->mkdir(writedir.c_str());

    // Start loop over wiremin objects
    // which is a bad way of saying that we actually loop over voltpc IDs
    for(auto& it_wmin : mToI.wiremin.at(plane))
      {
      size_t vtpc = it_wmin.first;
      // I check that the voltpc index also exists in wiremax map
      if (mToI.wiremax.at(plane).find(vtpc) == mToI.wiremax.at(plane).end()) throw cet::exception("SignalShapeAna") << "found a voltpc in wiremin map but not in wiremax map" << std::endl; 

      size_t eff_wmin = it_wmin.second;
      size_t eff_wmax = mToI.wiremax.at(plane)[vtpc];

      // eff_wmin and eff_wmax are the index of Wires in recob::Wire object
      // A correction regarding the plane ID and voltpc ID must be made
      // for instance, in coldbox v1,  induction 2 Wires ranges in [| 384 ; 1023 |]
      // Thus there are 640 Wires, 320 in each voltpc. 
      // [| 384 ; 703 |] for voltpc 0 and [| 704 ; 1023 |] for voltpc 2 (top electronics)
      // for the moment corrections applied are only valid for coldbox CRP1
      if (vtpc == 0 || vtpc == 1)
         {
         if (plane == 1) {eff_wmin += 384; eff_wmax += 384;}
         if (plane == 2) {eff_wmin += 1024; eff_wmax += 1024;}
         }
      if (vtpc == 2 || vtpc == 3)
         {
         if (plane == 0) {eff_wmin += 128; eff_wmax += 128;}
         if (plane == 1) {eff_wmin += 384 + 320; eff_wmax += 384 + 320;}
         if (plane == 2) {eff_wmin += 1024 + 288; eff_wmax += 1024 + 288;}
         }

      // check that effective wireIDs are not out of range, which for CRP1 coldbox is [0; 1599]
      if (eff_wmin < 0 || eff_wmin > 1599) throw::cet::exception("SignalShapeAna") << "effective wiremin index is out of range ! value = " << eff_wmin << std::endl;
      if (eff_wmax < 0 || eff_wmax > 1599) throw::cet::exception("SignalShapeAna") << "effective wiremax index is out of range ! value = " << eff_wmax << std::endl;

      // Start loop from effective wiremin to effective wiremax (both included)
      for (size_t iw = eff_wmin; iw <= eff_wmax; iw++)
        {
        StoreWireRegionOfInterestAsHistogram(mToI, vtpc, mWires.at(iw), arrayvecCenteredWaveforms.at(plane), dir);
        } // end loop over wires
      } // end for over map wiremin elements which are of type map({voltpc, min wire index})
  } // end loop over plane views


  return arrayvecCenteredWaveforms;
  }


// --------------------------------------------------------------------------------------------------------------------------


void SignalShapeAna::StoreWireRegionOfInterestAsHistogram(track_of_interest & mToI, size_t voltpc, recob::Wire mWire, std::vector<TH1F*>& _vecCenteredWaveforms, art::TFileDirectory dir)
  {

  size_t tmin = (size_t) mToI.tickmin;
  size_t tmax = (size_t) mToI.tickmax;

  
  // check that tmax didn't exceed the wire signal region
  // if tmax is greater than Nticks, set its value to NSignals() - 1, i.e last tick
  if (tmax >= mWire.NSignal()) tmax = mWire.NSignal() - 1;

  size_t Nbins = tmax - tmin + 1;


//    std::string str_plane = "";
//    if (mToI.plane == 0) str_plane = "ind1";
//    else if (mToI.plane == 1) str_plane = "ind2";
//    else if (mToI.plane == 2) str_plane = "col";
//    else str_plane = "unknown";
    std::string name = "vtpc" + std::to_string(voltpc) + "_" + std::to_string( static_cast<int>(mWire.Channel()) );
    TH1F * h = dir.make<TH1F>(name.c_str(), "Wire signal;ticks;ADC", Nbins, (double) tmin - 0.5, (double) tmax + 0.5);

    // Fill histo
    std::vector<float>  signal = mWire.Signal();
    for (int i=1; i<=h->GetNbinsX(); i++)
       {
       int ii = tmin + i - 1;
       h->SetBinContent(i, signal[ii]);
       }

    // add an histogram for which the waveform is centered around the max
    // caution with the number of bins to include
    int Nbins_centered = mDescendWindowSize + mRiseWindowSize + 1; // default case when the window is large enough to contain entirely the centered signal
    int refbin = h->GetMaximumBin(); // work to do

    name += "_centered";
    TH1F * hCentered = dir.make<TH1F>(name.c_str(), "Centered wire signal;tick", Nbins_centered, -0.5, Nbins_centered - 0.5);
    // keep track of this histogram
    _vecCenteredWaveforms.push_back(hCentered);

    // fill the histogram centered on the maximum
    for (int k = 1; k <= hCentered->GetNbinsX(); k++)
      {
      int effbin = refbin - mRiseWindowSize + k - 1;
      // test whether the effective bin does fall in the range [1; Nbins]
      if (effbin < 1 || effbin > h->GetNbinsX()) hCentered->SetBinContent(k, 0.);
      else hCentered->SetBinContent(k, h->GetBinContent(effbin));
      }


  return;
  } // endmethod StoreWireRegionOfInterestAsHistogram


// --------------------------------------------------------------------------------------------------------------------------



// note : the method assumes that all individual waveforms have the same binning convention
void SignalShapeAna::PerformCoherentAdditionOfWaveforms(track_of_interest mToI, std::array<std::vector<TH1F*>, 3> vecCenteredWaveforms, std::string mainwritedir)
  {


  // loop over plane views
  for (size_t plane = 0; plane < 3; plane++)
    {
    if (vecCenteredWaveforms.at(plane).size() == 0) continue;

    std::string writedir = mainwritedir;
    if (plane == 0) writedir += "/ind1";
    else if (plane == 1) writedir += "/ind2";
    else if (plane == 2) writedir += "/col";
    else writedir += "/unknown";
    art::TFileDirectory dir = tfs->mkdir(writedir.c_str());
/*
  //  art::ServiceHandle<art::TFileService> tfs;
    std::string strdir = "event" + std::to_string(mToI.event);
    if (plane == 0) strdir += "/ind1";
    else if (plane == 1) strdir += "/ind2";
    else if (plane == 2) strdir += "/col";
    art::TFileDirectory dir = tfs->mkdir(strdir.c_str());
*/
    // init coherent addition histogram
    size_t nbins = vecCenteredWaveforms.at(plane)[0]->GetNbinsX();
    double xlow = vecCenteredWaveforms.at(plane)[0]->GetBinLowEdge(1);
    double xup = vecCenteredWaveforms.at(plane)[0]->GetBinLowEdge(nbins) + vecCenteredWaveforms.at(plane)[0]->GetBinWidth(nbins);;
  
    // name of histogram
    // first get approx value of angles theta and phi that will go into the name
    int approxtheta;
    int approxphi;
    if (mToI.theta < 0) approxtheta = static_cast<int>(mToI.theta-0.5);
    else approxtheta = static_cast<int>(mToI.theta+0.5);
    if (mToI.phi < 0) approxphi = static_cast<int>(mToI.phi-0.5);
    else approxphi = static_cast<int>(mToI.phi+0.5);
    std::string name = "CoherentAddition_theta_" + std::to_string(approxtheta) + "_phi_" + std::to_string(approxphi);
    TH1D * hCoherent = dir.make<TH1D>(name.c_str(), "Signal coherently added", nbins, xlow, xup);


    if (fVerbose) std::cout << "will add coherently " << vecCenteredWaveforms.at(plane).size() << " waveforms on plane " << plane << "." << std::endl;



    std::vector<float> content(nbins, 0.);

    // start loop over centered waveforms
    int n = vecCenteredWaveforms.at(plane)[0]->GetNbinsX();
    int refbin = vecCenteredWaveforms.at(plane)[0]->GetMaximumBin(); // work to do
    for (size_t i = 0; i < vecCenteredWaveforms.at(plane).size(); i++)
       {
       // check that all waveforms have the same number of bins
       if (vecCenteredWaveforms.at(plane)[i]->GetNbinsX() != n) throw cet::exception("SignalShapeAna") << "Waveform histogram must have the same number of bins to be coherently added." << std::endl;
       // check that all waveforms are centered the same way
       if (vecCenteredWaveforms.at(plane)[i]->GetMaximumBin() != refbin) throw cet::exception("SignalShapeAna") << "Waveform histograms must be centered the same way." << std::endl;
       // retrieve content of waveform and add it to vector of size nbin
       for (int k = 0; k < vecCenteredWaveforms.at(plane)[i]->GetNbinsX(); k++) content[k] += vecCenteredWaveforms.at(plane)[i]->GetBinContent(k+1);
       } // end for over waveform histograms vector

    // Finally fill the histogram of coherent addition
    for (int k = 1; k <= hCoherent->GetNbinsX(); k++) hCoherent->SetBinContent(k, content[k-1]);

// check
//std::cout << "adding together " << vecCenteredWaveforms.at(plane).size() << " waveforms." << std::endl; 
//for (unsigned size_t k = 0; k <  content.size(); k++) std::cout << k << "\t" << content[k] << std::endl;

  } // end loop over plane views

  return;
  } // end method PerformCoherentAdditionOfWaveforms


// --------------------------------------------------------------------------------------------------------------------------



  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see HitEfficiencyAna.fcl for more information.
  DEFINE_ART_MODULE(SignalShapeAna)


} // namespace SignalShapeAna
#endif // HitEfficiencyAna_module



// TO DO LIST
/* 
	* make the module compatible with more greometries
	* there's a place with room for optimization, I don't remember where
	* in output root file, waveforms are not sorted by voltpc, make that happen
*/
