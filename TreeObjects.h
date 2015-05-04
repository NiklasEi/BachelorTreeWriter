#include <TLorentzVector.h>
#include <TVector3.h>

#ifndef TreeObjects_h
#define TreeObjects_h

namespace tree {
// In this namespace classes for the trees are defined.

enum genParticles {
	kGenPhoton,
	kGenElectron
};

class Particle {
	public:
		// functions
		float DeltaR( const Particle &p2 ) const;
		float DeltaR( const TLorentzVector &vec2 ) const;
		float DeltaPhi( float phi2 ) const;
		void setStatus( int status );
		bool isStatus( int status ) const;

		// variables
		float pt, eta, phi;
		short bitFlag;
		//variables for matching
		bool genPhoton;
		bool genElectron;
};

class Photon : public Particle {
	public:
		float sigmaIphiIphi;
		float r9, sigmaIetaIeta, hadTowOverEm;
		float chargedIso, neutralIso, photonIso;
		int pixelseed;
};

class Jet : public Particle{
	public:
		float bCSV;
		float chargedHadronEnergy,
			neutralHadronEnergy,
			photonEnergy,
			electronEnergy,
			muonEnergy,
			HFHadronEnergy,
			HFEMEnergy,
			chargedEmEnergy,
			chargedMuEnergy,
			neutralEmEnergy;
};

bool EtGreater(const tree::Particle, const tree::Particle);

} // end namespace definition

#endif



