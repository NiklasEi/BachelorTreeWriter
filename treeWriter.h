#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <set>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3I.h"
#include "TPRegexp.h"

#include "SusyEvent.h"
#include "TreeObjects.h"

using namespace std;

class TreeWriter {
  public :
    TreeWriter( int nFiles, char** fileList, std::string const& );
    virtual ~TreeWriter();
    virtual void Loop( );

    // Set tigger and input Files
    void SetTriggerPaths( std::vector<const char*> const & tp ) { triggerNames = tp; }
    void SetJsonFile( TString const & filename );
    void SetPileUpWeightHisto( TH1F histo ) { pileupHisto = histo; }
    void fillGenParticles();
    
    
  private:
    bool loggingVerbosity;

    bool isSignalScan;
    bool passTrigger();
    bool isGoodLumi() const;

    float getPileUpWeight();
    void fillJets();

    // Additional information for producing the output
    TH1F pileupHisto;

    std::map<unsigned, std::set<unsigned> > goodLumiList;
    std::vector<const char*> triggerNames;

    TChain inputTree;
    susy::Event event;


    TFile outFile;
    TTree outputTree;

    std::map< std::string, TH2F > hist2D;
    std::map< std::string, TH1F > hist1D;

    
	unsigned int nTracksPV;
    unsigned int nVertex;

    float ht;
    float met;
    float weight;
    
    std::vector<tree::Jet> jets;
    std::vector<tree::Photon> photons;
    std::vector<tree::Photon> electrons;
    std::vector<tree::Particle> genPhotons;
    std::vector<tree::Particle> genElectrons;

};
