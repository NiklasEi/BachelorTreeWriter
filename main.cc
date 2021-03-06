#include "treeWriter.h"
#include "TreeObjects.h"

template <class HIST>
HIST getHisto( std::string const & filename, std::string const & histname ) {
  TFile file( filename.c_str() );
  if( file.IsZombie() )
    std::cerr << "ERROR: Could not open file " << filename << std::endl;

  if( file.GetListOfKeys()->Contains( histname.c_str() ) )
    return *((HIST*) file.Get( histname.c_str() ));
  else
    std::cerr << "ERROR: Could not extract " << histname
      << " histogram from " << filename << std::endl;
  return HIST();
}

int main( int argc, char** argv ) {

  if( argc < 3 ) {
    std::cout << "usage: ./execute outputFileName.root input1.root input2.root ..." << std::endl;
    return 1;
  }


  TreeWriter tw( argc-2, argv+2, argv[1] );


  const std::string lumiJsonName = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt";
  if( access( lumiJsonName.c_str(), F_OK ) != -1 ) {
    tw.SetJsonFile( lumiJsonName );
  } else {
    // try to get local cert file
    tw.SetJsonFile( "Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt" );
  }

  gSystem->Load("libHistPainter"); // to avoid waring and errors when reading th2 from file


  std::string pileupScenario = "S10";

  std::string puWeightFile = "pileUpReweighting/puWeights.root";

  tw.SetPileUpWeightHisto( getHisto<TH1F>( puWeightFile, "pileupWeight"+pileupScenario ) );

  std::vector<const char*> triggerNames;
  triggerNames.push_back( "HLT_Photon135_v" );
  tw.SetTriggerPaths( triggerNames );

  double start_time = time(NULL);
  tw.Loop();
  double end_time = time(NULL);

  std::cout << "Job needed " << 1.*(end_time - start_time)/3600 << " h real time." << std::endl;
  return 0;
}

