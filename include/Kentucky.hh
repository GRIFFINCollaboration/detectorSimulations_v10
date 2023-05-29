////////////////////////////////////////////////////////////////////////////////
// this class reads the (pre-formatted) ENDF files to generate angular distributions
// the files can either contain legendre coefficients (default)
// or a set of cosine(theta_cm), cross section pairs for a set of beam energies
////////////////////////////////////////////////////////////////////////////////
#ifndef __KENTUCKY_HH
#define __KENTUCKY_HH

#include <string>
#include <vector>

#include "G4String.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"

#include "TF1.h"

#define NEUTRON_MASS   (939.565379)
#define PROTON_MASS    (938.272046)
#define DEUTERON_MASS (1875.612859)
#define TRITON_MASS   (2808.921005)
#define HE3_MASS      (2808.391482)
#define HE4_MASS      (3727.379240)
#define C12_MASS     (11177.92873)

class Kentucky {
public:
	enum class EReaction { kPT, kDT, kDD, kUndefined };

	Kentucky(const G4double& minimumPhi, const G4double& maximumPhi, const G4double& minimumTheta, const G4double& maximumTheta);
	~Kentucky() { delete fCrossSectionFunction; };

	// setter
	void Energy(double val);
	void Reaction(G4String reaction);
	void VerbosityLevel(int val) { fVerboseLevel = val; }

	// function that calculates direction and energy and sets it on particle gun
	void DirectionAndEnergy(G4ParticleGun* particleGun) const;

	double operator() (double* x, double*) { return CrossSection(x[0]); }

	void Print() const;

private:
	double Angle(const double&, const double&) const;
	double SolidAngle(const double&, const double&) const;

	EReaction fReaction;

	bool Open(const char* fileName, bool tabulated = false);
	bool Open(std::string fileName, bool tabulated = false) {
		return Open(fileName.c_str(), tabulated);
	};

	double CrossSection(double) const;
	double CrossSection(double,std::vector<double>) const;

	int fVerboseLevel;

	std::string fFilename;

	TF1* fCrossSectionFunction;
	double fMax; // maximum value of cross section

	double fMinimumPhi;
	double fMaximumPhi;
	double fMinimumTheta;
	double fMaximumTheta;

	double fZeroDegreeEnergy;
	double fBeamEnergy;

	double fProjectileMass;
	double fTargetMass;
	double fRecoilMass;
	double fEjectileMass;
	std::vector<double> fEnergy;
	std::vector<std::vector<double> > fCosThCm;
	std::vector<std::vector<double> > fCoeff;
	std::vector<double> fCrossSection;
};

#endif
