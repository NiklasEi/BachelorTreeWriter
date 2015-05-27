#include "treeWriter.h"


float effectiveAreaElectron( float eta ) {
   /** Returns the effective area for the isolation criteria for electrons.
    * See https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection
    * only for Delta R = 0.4 on 2012 Data
    */
    eta = fabs( eta );
    float ea;

   if( eta < 1.0 ) ea = 0.208;
   else if( eta < 1.479 ) ea = 0.209;
   else if( eta < 2.0 ) ea = 0.115;
   else if( eta < 2.2 ) ea = 0.143;
   else if( eta < 2.3 ) ea = 0.183;
   else if( eta < 2.4 ) ea = 0.194;
   else ea = 0.261;

   return ea;
}

float chargedHadronIso_corrected(const susy::Photon& gamma, float rho) {
   /** Correct isolation for photons,
    * see https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonID2012
    */
   float eta = fabs(gamma.caloPosition.Eta());
   float ea;

   if(eta < 1.0) ea = 0.012;
   else if(eta < 1.479) ea = 0.010;
   else if(eta < 2.0) ea = 0.014;
   else if(eta < 2.2) ea = 0.012;
   else if(eta < 2.3) ea = 0.016;
   else if(eta < 2.4) ea = 0.020;
   else ea = 0.012;

   float iso = gamma.chargedHadronIso;
   iso = std::max(iso - rho*ea, (float)0.);

   return iso;
}

float neutralHadronIso_corrected(const susy::Photon& gamma, float rho) {
   float eta = fabs(gamma.caloPosition.Eta());
   float ea;

   if(eta < 1.0) ea = 0.030;
   else if(eta < 1.479) ea = 0.057;
   else if(eta < 2.0) ea = 0.039;
   else if(eta < 2.2) ea = 0.015;
   else if(eta < 2.3) ea = 0.024;
   else if(eta < 2.4) ea = 0.039;
   else ea = 0.072;

   float iso = gamma.neutralHadronIso;
   iso = std::max(iso - rho*ea, (float)0.);

   return iso;
}

float photonIso_corrected(const susy::Photon& gamma, float rho) {
   float eta = fabs(gamma.caloPosition.Eta());
   float ea;

   if(eta < 1.0) ea = 0.148;
   else if(eta < 1.479) ea = 0.130;
   else if(eta < 2.0) ea = 0.112;
   else if(eta < 2.2) ea = 0.216;
   else if(eta < 2.3) ea = 0.262;
   else if(eta < 2.4) ea = 0.260;
   else ea = 0.266;

   float iso = gamma.photonIso;
   iso = std::max(iso - rho*ea, (float)0.);

   return iso;
}

bool isLooseJet( const susy::PFJet& jet ) {
   /**
    * \brief Apply loose cut on jets.
    *
    * See https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_7_TeV_data_a
    * for more information.
    */

   double energy = jet.momentum.E();
   return (jet.neutralHadronEnergy+jet.HFHadronEnergy) / energy < 0.99
      && jet.neutralEmEnergy / energy < 0.99
      && jet.nConstituents > 1
      && ( std::abs(jet.momentum.Eta()) >= 2.4
      || ( jet.chargedHadronEnergy / energy > 0
      && jet.chargedMultiplicity > 0
      && jet.chargedEmEnergy / energy < 0.99 ) );
}

bool goodVertex( const susy::Vertex& vtx ) {
   /** Definition of a good vertex. Returns true if the vertex is good.
   */
   return (!vtx.isFake() &&
      vtx.ndof >= 4 &&
      std::abs((vtx.position).z()) < 24.0 &&
      std::abs((vtx.position).Perp()) < 2.0 );
}

unsigned int numberOfGoodVertexInCollection( const std::vector<susy::Vertex>& vertexVector ) {
   // Counts the number of good vertices in the vertex Vector
   
   unsigned int number = 0;
   for( std::vector<susy::Vertex>::const_iterator vtx = vertexVector.begin();
      vtx != vertexVector.end(); ++vtx ) {
      if( goodVertex( *vtx ) )
         number++;
   }
   return number;
}

unsigned int nTrackPrimaryVertex( const std::vector<susy::Vertex>& vertexVector ) {
   // Tracks coming from the first good vertex 
   for( std::vector<susy::Vertex>::const_iterator vtx = vertexVector.begin();
      vtx != vertexVector.end(); ++vtx ) {
      if( goodVertex( *vtx ) )
         return vtx->tracksSize;
   }
   // if no valid vertex was found
   return 0;
}



template <typename VectorClass>
int indexOfnearestParticle( const tree::Particle& thisParticle, const std::vector<VectorClass>& particleVector,
      float deltaR_=.3, float ptRelMin_=-1e6, float ptRelMax_=1e6, TH2F* hist=NULL ) {
   /* Compare 'thisParticle' to each particle in the vector
    *
    * If a particle in the vector satisfies the cut values, its index is returned.
    * If no particle if found, -1 is returned.
    * If serveral particles satisfy the requirements, the nearest in deltaR is
    * returned.
    */
   int index = -1; //default value
   std::map<float,int> map_dr_i;
   float dR, relPt;
   //cout<<"test1 "<<particleVector.size()<<endl;
   for( typename std::vector<VectorClass>::const_iterator it = particleVector.begin();
         it != particleVector.end(); ++it ) {

      relPt = it->pt / thisParticle.pt;
      dR = thisParticle.DeltaR( *it );

      if( hist ) hist->Fill( dR, relPt );

      if( dR > deltaR_ ||  relPt > ptRelMax_ || relPt < ptRelMin_ ) continue;
      index = std::distance( particleVector.begin(), it );
      map_dr_i[dR] = index;
   }

   // search nearest if several particles found
   if( map_dr_i.size() > 1 ) {
      index = map_dr_i.begin()->second;
   }

   return index;
}

///////////////////////////////////////////////////////////////////////////////
// Here the class implementation begins ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

TreeWriter::TreeWriter( int nFiles, char** fileList, std::string const& outputName ) :
   loggingVerbosity(0),
   inputTree("susyTree"),
   event(),
   outFile( outputName.c_str(), "recreate" ),
   outputTree("myTree","Tree for single photon analysis")
{
   for( int i = 0; i<nFiles; ++i ){
      inputTree.Add( fileList[i] );
   }

   event.setInput( inputTree );

   outputTree.Branch("jets", &jets);
   outputTree.Branch("jetphotons", &jetphotons);
   outputTree.Branch("photons", &photons);
   outputTree.Branch("electrons", &electrons);
   outputTree.Branch("met", &met, "met/F");
   outputTree.Branch("ht", &ht, "ht/F");
   outputTree.Branch("weight", &weight, "weight/F");
   outputTree.Branch("nVertex", &nVertex, "nVertex/I");
   outputTree.Branch("nTracksPV", &nTracksPV, "nTracksPV/i");

   hist2D["matchGenPhoton"]   = TH2F("matchGenPhoton", ";#DeltaR;p_{T}^{gen} / p_{T}", 1000, 0, .5, 200, 0, 2 );
   hist2D["matchGenElectron"] = TH2F("matchGenElectron", ";#DeltaR;p_{T}^{gen} / p_{T}", 1000, 0, .5, 200, 0, 2 );
   hist1D["nGen"] = TH1F("", ";nGen;", 1, 0, 1 );
}

TreeWriter::~TreeWriter() {
   /** Deconstructor
    * Event has to be deleted before the deletion of the tree.
    */
   event.releaseTree(inputTree);
}

void TreeWriter::SetJsonFile(TString const& filename) {
   /** Read a Json file which contains good runNumbers and Lumi-sections.
    * The content will be stored in the class variable 'goodLumiList'.
    */

   if( goodLumiList.size() ) {
      std::cout << "WARNING: Good lumi sections already defined. Now overwritten by " << filename << std::endl;
      return;
   }

   std::ifstream inputFile( filename );
   if( !inputFile.is_open() )
      std::cerr << "ERROR: Cannot open JSON file " << filename << std::endl;

   std::string line;
   TString jsonText;
   while(inputFile.good()){
      getline(inputFile, line);
      jsonText += line;
   }
   inputFile.close();

   TPRegexp runBlockPat("\"([0-9]+)\":[ ]*\\[((?:\\[[0-9]+,[ ]*[0-9]+\\](?:,[ ]*|))+)\\]");
   TPRegexp lumiBlockPat("\\[([0-9]+),[ ]*([0-9]+)\\]");

   TArrayI positions(2);
   positions[1] = 0;
   while(runBlockPat.Match(jsonText, "g", positions[1], 10, &positions) == 3){
      TString runBlock(jsonText(positions[0], positions[1] - positions[0]));
      TString lumiPart(jsonText(positions[4], positions[5] - positions[4]));

      unsigned run(TString(jsonText(positions[2], positions[3] - positions[2])).Atoi());
      std::set<unsigned>& lumis(goodLumiList[run]);

      TArrayI lumiPos(2);
      lumiPos[1] = 0;
      while(lumiBlockPat.Match(lumiPart, "g", lumiPos[1], 10, &lumiPos) == 3){
         TString lumiBlock(lumiPart(lumiPos[0], lumiPos[1] - lumiPos[0]));
         int begin(TString(lumiPart(lumiPos[2], lumiPos[3] - lumiPos[2])).Atoi());
         int end(TString(lumiPart(lumiPos[4], lumiPos[5] - lumiPos[4])).Atoi());
         
         for(int lumi(begin); lumi <= end; ++lumi)
            lumis.insert(lumi);
      }
   }
   if( loggingVerbosity > 1 )
      std::cout << "JSON file for filtering included." << std::endl;
}



void TreeWriter::fillGenParticles() {
  /**fill in generated particles for matching
   * 
   */

   genPhotons.clear();
   genElectrons.clear();

   // genParticles
   tree::Particle thisGenParticle;
   for( susy::ParticleCollection::const_iterator it = event.genParticles.begin();
     it != event.genParticles.end(); ++it ) {

      // status 3: particles in matrix element
      // status 2: intermediate particles
      // status 1: final particles (but can decay in geant, etc)
      if( it->momentum.Pt() < 40 || it->status != 1) continue;

      thisGenParticle.pt = it->momentum.Pt();
      thisGenParticle.eta = it->momentum.Eta();
      thisGenParticle.phi = it->momentum.Phi();
      //thisGenParticle.bitFlag = 0;
      int pdgId = std::abs(it->pdgId);
      switch( pdgId ) {
         case 22: // photon
            genPhotons.push_back( thisGenParticle );
            break;
         case 11: // electron
            genElectrons.push_back( thisGenParticle );
            break;
      }
   }

   if( loggingVerbosity > 1 )
      std::cout << "Found " << genPhotons.size() << " generated photons and "
      << genElectrons.size() << " generated electrons" << std::endl;
}


bool TreeWriter::passTrigger() {
  /**
   * Checks if event passes the HLT trigger paths.
   *
   * If the a trigger path contains one of the triggers defined in triggerNames,
   * true will be returned. If no match is found, false is returned.
   */

   if( !triggerNames.size() ) {
      std::cout << "WARNING: No triggers selected" << std::endl;
      return true;
   }

   for( std::vector<const char*>::const_iterator it = triggerNames.begin();
     it != triggerNames.end(); ++it ) {
      for( susy::TriggerMap::const_iterator tm = event.hltMap.begin();
        tm != event.hltMap.end(); ++tm ) {
         if ( tm->first.Contains( *it ) && (int(tm->second.second))) {
            return true;
            if( loggingVerbosity > 1 )
               std::cout << "Pass trigger requirement." << std::endl;
         }
      }
   }
   if( loggingVerbosity > 1 )
      std::cout << "Fail trigger requirement." << std::endl;
   return false;
}

bool TreeWriter::isGoodLumi() const {
  /**
   * Check if current event is in json file added by 'IncludeAJson(TString )'
   */
   if( !goodLumiList.size() ) {
      std::cout << "WARNING: No json file for filtering found" << std::endl;
      return true;
   }
   bool goodLumi = false;
   unsigned run = event.runNumber;
   unsigned lumi = event.luminosityBlockNumber;
   std::map<unsigned, std::set<unsigned> >::const_iterator rItr(goodLumiList.find(run));
   if(rItr != goodLumiList.end()){
      std::set<unsigned>::const_iterator lItr(rItr->second.find(lumi));
      if(lItr != rItr->second.end()) goodLumi = true;
   }
   if( loggingVerbosity > 1 && goodLumi )
      std::cout << "Event is in a good lumi section." << std::endl;
   if( loggingVerbosity > 1 && !goodLumi )
      std::cout << "Event is not in a good lumi section." << std::endl;
   return goodLumi;
}

float TreeWriter::getPileUpWeight(){
  /**
   * If a pileup weight histogram has been added, the pile-up weight for the
   * current event is computed.
   */

   float thisWeight = 0;
   for( susy::PUSummaryInfoCollection::const_iterator iBX = event.pu.begin();
      iBX != event.pu.end(); ++iBX) {
      if (iBX->BX == 0) { // find bunch crossing for this event
         float trueNumInteractions = iBX->trueNumInteractions;
         thisWeight = pileupHisto.GetBinContent( pileupHisto.FindBin( trueNumInteractions ) );
         break;
      }
   }

   if( loggingVerbosity > 2 )
      std::cout << "Pile-up weight = " << thisWeight << std::endl;
   return thisWeight;
}


float TreeWriter::getHt() const {
	/* HT is sum jet + sum photon + sum photonJet + sum photonElectron
	 * For the jet sum, jets have to have a good Id or was used as match for a
	 * photon/photonJet/photonElectron.
	 * For the sum of photonObjects, the pt is only added in case the pt of the
	 * matched jet was not added.
	 */

	float returnedHt = 0;
	for(std::vector<tree::Jet>::const_iterator jet = jets.begin();
			jet != jets.end(); ++jet ) {

		if( jet->pt < 40 || std::abs(jet->eta) > 3. ) continue;
		//if( !jet->isStatus( tree::kJetId ) ) continue;
		//std::cout << " add jet to HT " << jet->pt << std::endl;

		returnedHt += jet->pt;
	}

	for( std::vector<tree::Photon>::const_iterator photon = photons.begin();
			photon != photons.end(); ++photon ) {
		if( photon->ptMJet == 0 ) {
			//std::cout << " add photon to HT " << photon->pt << std::endl;
			returnedHt += photon->pt;
		}
	}
	for( std::vector<tree::Photon>::const_iterator photon = jetphotons.begin();
			photon != jetphotons.end(); ++photon ) {
		if( photon->ptMJet == 0 ) {
			//std::cout << " add photonJet to HT " << photon->pt << std::endl;
			returnedHt += photon->pt;
		}
	}
	for( std::vector<tree::Photon>::const_iterator photon = electrons.begin();
			photon != electrons.end(); ++photon ) {
		if( photon->ptMJet == 0 ) {
			//std::cout << " add photonElectron to HT " << photon->pt << std::endl;
			returnedHt += photon->pt;
		}
	}
	return returnedHt;
}




void TreeWriter::fillJets() {
  /** Read the jets from susyEvent and save them to jet vector.
   * All jets for HT calculation, and photon-jet matching are saved.
   * This are not the final jets in the analysis.
   * All jets are corrected
   */
   jets.clear();
   tree::Jet jetToTree;

   std::vector<susy::PFJet> jetVector = event.pfJets.find("ak5chs")->second;
   for(std::vector<susy::PFJet>::const_iterator it = jetVector.begin();
     it != jetVector.end(); ++it) {

      TLorentzVector corrP4 = it->jecScaleFactors.at("L1FastL2L3") * it->momentum;

      if( std::abs(corrP4.Eta()) > 3 ) continue;
      if( corrP4.Pt() < 40 ) continue;
      if( !isLooseJet( *it ) ) continue;


      jetToTree.pt = corrP4.Pt();
      jetToTree.eta = corrP4.Eta();
      jetToTree.phi = corrP4.Phi();
      jets.push_back( jetToTree );

      if( loggingVerbosity > 2 )
         std::cout << " p_T, jet = " << jetToTree.pt << std::endl;
   }// for jet

   if( loggingVerbosity > 1 ){
      std::cout << "Found " << jets.size() << " uncleaned jets" << std::endl;
   }
}

void TreeWriter::Loop() {
  /**
   * \brief Loops over input chain and fills tree
   *
   * This is the major function of treeWriter, which initializes the output, loops
   * over all events and fills the tree. In the end, the tree is saved to the
   * output File
   */

// Declaration for objects saved in Tree
tree::Photon photonToTree;

   for (long jentry=0; jentry < inputTree.GetEntries(); ++jentry) {

      //if( jentry > 1000 ) break;



      event.getEntry(jentry);

      // For data, the weight is 1. Else take the pileup weight.
      weight = event.isRealData ? 1. : getPileUpWeight();
      hist1D["nGen"].Fill( 0 );
      if ( event.isRealData && !isGoodLumi() ) continue;
      if ( event.isRealData && !passTrigger() ) continue;
      
      
      //Jets einlesen
      fillJets();
      //Generierte Teilchen einlesen
      fillGenParticles();
      susy::MET pfMet = event.metMap["pfMet"];
      met = pfMet.met();
      
      
      // vertices
	  nVertex = numberOfGoodVertexInCollection( event.vertices );
	  if( !nVertex ){
         continue;
      }
	  
      
      nTracksPV = nTrackPrimaryVertex( event.vertices );
      if( loggingVerbosity > 2 ){
         std::cout << " nTracksPV = " << nTracksPV << std::endl;
      }
			

      photons.clear();
      electrons.clear();
      jetphotons.clear();
      //photons
      std::vector<susy::Photon> photonVector = event.photons["photons"];
      for(std::vector<susy::Photon>::iterator it = photonVector.begin();
        it != photonVector.end(); ++it ) {
		

         //dont use endcap
         if( std::abs( it->momentum.Eta() ) > susy::etaGapBegin ) continue;

         //define variables to save in the tree
         photonToTree.chargedIso = chargedHadronIso_corrected(*it, event.rho);
         photonToTree.neutralIso = neutralHadronIso_corrected(*it, event.rho);
         photonToTree.photonIso = photonIso_corrected(*it, event.rho);
         photonToTree.pt = it->momentum.Pt();
         photonToTree.eta = it->momentum.Eta();
         photonToTree.phi = it->momentum.Phi();
         photonToTree.r9 = it->r9;
         photonToTree.sigmaIetaIeta = it->sigmaIetaIeta;
         photonToTree.sigmaIphiIphi = it->sigmaIphiIphi;
         photonToTree.hadTowOverEm = it->hadTowOverEm;
         photonToTree.pixelseed = it->nPixelSeeds;
         //photonToTree.conversionSafeVeto = it->passelectronveto;
         //
         photonToTree.genPhoton = false;
         photonToTree.genElectron = false;
         

         int jetIndex = indexOfnearestParticle<tree::Jet>( photonToTree, jets, .2, .8, 3 );
         photonToTree.ptMJet = jetIndex>-1 ? jets.at(jetIndex).pt : 0.;
         photonToTree.etaMJet = jetIndex>-1 ? jets.at(jetIndex).eta : 0.;
         photonToTree.phiMJet = jetIndex>-1 ? jets.at(jetIndex).phi : 0.;
         photonToTree.matchedJetIndex = jetIndex;
         //cout<<"matchedJetIndex: "<<photonToTree.matchedJetIndex<<"       jetIndex:"<<jetIndex<<endl;
         //if(!photonToTree.ptMJet==0)
         //   cout<<"ptMJet1: "<<photonToTree.ptMJet<<endl;


         //photon definition barrel
         bool isPhotonOrElectron =
         ( std::abs(photonToTree.eta) <= susy::etaGapBegin
         && photonToTree.ptJet() > 145                                            // Min pt 
         && photonToTree.hadTowOverEm < 0.05
         && photonToTree.sigmaIetaIeta < 0.012
         && photonToTree.chargedIso < 2.6
         && photonToTree.neutralIso < 3.5+0.04*photonToTree.pt
         && photonToTree.photonIso < 1.3+0.005*photonToTree.pt
         )
         // and the endcap definition, which is not used now
         || ( std::abs( photonToTree.eta ) >= susy::etaGapEnd
         && std::abs( photonToTree.eta ) <= susy::etaMax
         && photonToTree.hadTowOverEm < 0.05
         && photonToTree.sigmaIetaIeta < 0.034
         && photonToTree.chargedIso < 2.3
         && photonToTree.neutralIso < 2.9+0.04*photonToTree.pt
         );
         
         //deferring between Photons and electrons via pixel-veto
         bool isPhoton = isPhotonOrElectron && !photonToTree.pixelseed;
         bool isPhotonElectron = isPhotonOrElectron && photonToTree.pixelseed;



         // photonJet definition
         bool isPhotonJet = !isPhotonOrElectron
         && photonToTree.ptJet() > 145
         && !photonToTree.pixelseed
         && photonToTree.hadTowOverEm < 0.05
         && photonToTree.sigmaIetaIeta < 0.012
         && photonToTree.chargedIso < 26 && photonToTree.chargedIso > 0.26
         && photonToTree.neutralIso < 35+0.4*photonToTree.pt && photonToTree.neutralIso > 0.35+0.004*photonToTree.pt
         && photonToTree.photonIso < 13+0.05*photonToTree.pt && photonToTree.photonIso > 0.13+0.0005*photonToTree.pt;


         if( loggingVerbosity > 2 ){ 
            if( isPhoton ) std::cout << " photon pT = " << photonToTree.pt << std::endl;
         }

         //matching generated particles to detected ones
         if( indexOfnearestParticle<tree::Particle>( photonToTree, genPhotons, .1, 0.9, 1.1, &hist2D["matchGenPhoton"] ) > -1 ){
            photonToTree.genPhoton = true;
         }
         if( indexOfnearestParticle<tree::Particle>( photonToTree, genElectrons, .1, -1e6, 1e6, &hist2D["matchGenElectron"] ) > -1 ){
            photonToTree.genElectron = true;
         }


	     //pushing candidates back
         if( isPhoton ){
            photons.push_back( photonToTree );
         }
         if( isPhotonElectron ){
            electrons.push_back( photonToTree );
         }
         if( isPhotonJet ){
            jetphotons.push_back( photonToTree );
         }
         

      }//finished pushing back photons and electrons (already matched if possible)
      if ( !photons.size() && !electrons.size() && !jetphotons.size()) continue;
      /*
      if(photons.size()==1){
		  if(!photons[0].ptMJet==0)
		     cout<<"ptMJet2: "<<photons[0].ptMJet<<"       #################################"<<endl;
		  
	  }
	  */
      if( loggingVerbosity > 1 ){
         std::cout << "Found " << photons.size() << " photons, "<<std::endl;
      }

      if( !event.passMetFilters() || !event.passMetFilter( susy::kEcalLaserCorr) ){
         continue;
      }
		  

      //calculate ht
      ht = getHt();
      //double oldht=0;
      //for(auto Jet : jets){
	  //   oldht+= Jet.pt;
      //}
      //cout<<"new ht: "<<ht<<"     old ht: "<<oldht<<endl;
   
      outputTree.Fill();
      
   }//finished loop over all tree-entries


   outFile.cd();
   outputTree.Write();
    
   for( auto m: hist2D){
	   m.second.SetName(m.first.c_str());
	   m.second.Write();
   }
   for( auto m: hist1D){
	   m.second.SetName(m.first.c_str());
	   m.second.Write();
   }

}


