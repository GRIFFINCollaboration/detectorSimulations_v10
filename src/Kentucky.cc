#include "Kentucky.hh"

#include <iostream>
#include <fstream>
#include <sstream>

#include <gsl/gsl_sf_legendre.h>

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "CLHEP/Units/PhysicalConstants.h"

#include "PhysicsUtilities.hh"

Kentucky::Kentucky(const G4double& minimumPhi, const G4double& maximumPhi, const G4double& minimumTheta, const G4double& maximumTheta)
	: fMinimumPhi(minimumPhi), fMaximumPhi(maximumPhi), fMinimumTheta(minimumTheta), fMaximumTheta(maximumTheta)
{
	fReaction = EReaction::kUndefined;
	fVerboseLevel = 0;

	fFilename = "";
	fCrossSectionFunction = new TF1("CrossSection", this, 0.*degree, 180.*degree, 0);
	fMax = 0.;

	fZeroDegreeEnergy = 0.;
	fBeamEnergy = 0.;

	fProjectileMass = 0.;
	fTargetMass = 0.;
	fRecoilMass = 0.;
	fEjectileMass = 0.;
	fEnergy.clear();
	fCosThCm.clear();
	fCoeff.clear();
	fCrossSection.clear();
}

bool Kentucky::Open(const char* fileName, bool tabulated)
{
	fFilename = fileName;
	std::ifstream input(fileName);

	if(!input.is_open()) {
		return false;
	}

	std::string line;
	std::istringstream stream;
	double tmpDouble;

	//the first line should have projectile, target, recoil, and ejectile masses (after any comments)
	while(input.good()) {
		getline(input,line);
		if(line[0] == '#') {
			continue;
		}
		stream.clear();
		stream.str(line);

		stream>>fProjectileMass>>fTargetMass>>fEjectileMass>>fRecoilMass;
		if(fVerboseLevel > 1) {
			G4cout<<__PRETTY_FUNCTION__<<"file '"<<fileName<<"': projectile mass = "<<fProjectileMass<<", target mass = "<<fTargetMass<<", recoil mass = "<<fRecoilMass<<", and ejectile mass = "<<fEjectileMass<<G4endl;
		}
		break;
	}

	//now the file itself follows
	while(input.good()) {
		getline(input,line);
		if(!tabulated) {
			if(line[0] == '#') {
				continue;
			}
			//line has beam energy, coeff (variable), cross section
			stream.clear();
			stream.str(line);

			fEnergy.push_back(0.);  
			fCoeff.push_back(std::vector<double>());

			stream>>fEnergy.back();
			fEnergy.back() *= CLHEP::eV; // ENDF is in eV
			while(stream.good()) {
				stream>>tmpDouble;
				fCoeff.back().push_back(tmpDouble);
			}
			//the last "coefficient" we read in is actually the cross section
			fCrossSection.push_back(fCoeff.back().back());
			fCoeff.back().pop_back();

			if(fVerboseLevel > 3) {
				G4cout<<__PRETTY_FUNCTION__<<fEnergy.back()<<": read legendre coeff. ";
				for(size_t i = 0; i < fCoeff.back().size(); ++i) {
					if(i != fCoeff.back().size()-1) {
						G4cout<<fCoeff.back()[i]<<", ";
					} else {
						G4cout<<"; cross-sec. "<<fCrossSection.back()<<G4endl;
					}
				}
			}
		} else {
			if(line.find("Incident Neutron Energy: ") == 0) {
				//start of new energy, energy given at the end of this line
				stream.clear();
				stream.str(line.substr(25));
				fEnergy.push_back(0.);  
				fCosThCm.push_back(std::vector<double>());
				fCoeff.push_back(std::vector<double>());
				stream>>fEnergy.back();
				fEnergy.back() *= CLHEP::eV; // ENDF is in eV
				if(fVerboseLevel > 3) {
					G4cout<<__PRETTY_FUNCTION__<<fEnergy.back()<<": looking for coeff.:"<<G4endl;
				}
			} else if(line.find("Cosine Value    Normalized Probability") == 0) {
				//start of a new block of angle/probability pairs
				//ended by a line starting with "ENDF/B-VII.1 Angular Distribution of Neutrons:"
				while(input.good()) {
					getline(input,line);
					if(line.find("ENDF/B-VII.1 Angular Distribution of Neutrons:") == 0) {
						break;
					}
					//line has angle, coeff
					stream.clear();
					stream.str(line);
					fCosThCm.back().push_back(0.);
					fCoeff.back().push_back(0.);
					stream>>fCosThCm.back().back()>>fCoeff.back().back();
					if(fVerboseLevel > 3) {
						G4cout<<__PRETTY_FUNCTION__<<fCosThCm.back().size()<<": theta = "<<acos(fCosThCm.back().back())/degree<<", cos(theta) = "<<fCosThCm.back().back()<<", coeff = "<<fCoeff.back().back()<<G4endl;
					}
				}
			}
		}//else of !tabulated
	}//while input is good

	if(fEnergy.size() < 2) {
		G4cout<<__PRETTY_FUNCTION__<<"failed to read at least two entries (got "<<fEnergy.size()<<" entries)!"<<G4endl;
		return false;
	}

	if(fVerboseLevel > 0) {
		G4cout<<__PRETTY_FUNCTION__<<"done reading "<<fEnergy.size()<<" entries ranging from "<<fEnergy[0]/CLHEP::keV<<" to "<<fEnergy.back()/CLHEP::keV<<" keV"<<G4endl;
	}

	return true;
}

//function with beam energy as parameter, returns cm-angle as function of lab-angle
double Kentucky::Angle(const double& x, const double& par) const
{
	if(fVerboseLevel > 3) {
		G4cout<<G4endl<<__PRETTY_FUNCTION__<<" parameters "<<x<<", "<<par<<G4endl;
		G4cout<<__PRETTY_FUNCTION__<<" "<<CLHEP::pi/rad<<" - acos(Physics::CosineThetaEjectileCm("<<par<<", "<<fProjectileMass<<", "<<fTargetMass<<", "<<fEjectileMass<<", "<<fRecoilMass<<", 0., "<<x/rad<<")) = "<<CLHEP::pi/rad<<" - acos("<<Physics::CosineThetaEjectileCm(par, fProjectileMass, fTargetMass, fEjectileMass, fRecoilMass, 0., x/rad)<<") = "<<CLHEP::pi/rad<<" - "<<acos(Physics::CosineThetaEjectileCm(par, fProjectileMass, fTargetMass, fEjectileMass, fRecoilMass, 0., x/rad))<<" = "<<CLHEP::pi/rad - acos(Physics::CosineThetaEjectileCm(par, fProjectileMass, fTargetMass, fEjectileMass, fRecoilMass, 0., x/rad))<<G4endl;
	}
	return CLHEP::pi/rad - acos(Physics::CosineThetaEjectileCm(par, fProjectileMass, fTargetMass, fEjectileMass, fRecoilMass, 0., x/rad));//excitation energy is always zero!
}

//function with beam energy as parameter, returns dOmega_cm/dOmega_lab as function of cm-angle
double Kentucky::SolidAngle(const double& x, const double& par) const
{
	return Physics::SolidAngleConversionEjectile(par, fProjectileMass, fTargetMass, fEjectileMass, fRecoilMass, 0., x/rad);
}

double Kentucky::CrossSection(double angle) const
{
	double conversion = 1.;
	//convert lab-frame angle to cm-angle
	if(fVerboseLevel > 3) {
		G4cout<<__PRETTY_FUNCTION__<<" beam energy = "<<fBeamEnergy/CLHEP::MeV<<" MeV, projectile mass = "<<fProjectileMass<<", target mass = "<<fTargetMass<<", recoil mass = "<<fRecoilMass<<", and ejectile mass = "<<fEjectileMass<<", lab-angle = "<<angle/deg<<"/"<<angle/rad;
	}
	//angle = Angle(angle/rad, fBeamEnergy/CLHEP::MeV)*rad;
	if(fVerboseLevel > 3) {
		G4cout<<" => cm-angle = "<<angle/deg<<G4endl;
	}
	//if we're in the lab-frame we also need to convert dsigma/dOmega
	conversion = 1.;//Physics::SolidAngleConversionEjectile(fBeamEnergy/CLHEP::MeV, fProjectileMass, fTargetMass, fEjectileMass, fRecoilMass, 0., angle/rad);

	//check energy 
	if(fEnergy.size() == 0) {
		if(fVerboseLevel > 1) {
			G4cout<<__PRETTY_FUNCTION__<<" energy range is zero: "<<fEnergy.size()<<G4endl;
		}
		return 0.;
	}

	//interpolate the coefficents, first find the bracketing pair of entries
	size_t index;
	for(index = 0; index < fEnergy.size()-1; ++index) {
		if(fBeamEnergy < fEnergy[index+1]) {
			break;
		}
	}

	if(fCosThCm.size() == 0) {
		//linear interpolation of all coefficients
		std::vector<double> coeff = fCoeff[index];
		if(fVerboseLevel > 2) {
			G4cout<<__PRETTY_FUNCTION__<<" "<<index<<": "<<coeff.size()<<", "<<fEnergy.size()<<", "<<fCoeff.size()<<G4endl;
		}
		for(size_t i = 0; i < coeff.size(); ++i) {
			coeff[i] += (fBeamEnergy - fEnergy[index]) * (fCoeff[index+1][i] - fCoeff[index][i])/(fEnergy[index+1] - fEnergy[index]);
		}
		//linear interpolation of cross section
		double crossSection = fCrossSection[index] + (fBeamEnergy - fEnergy[index]) * (fCrossSection[index+1] - fCrossSection[index])/(fEnergy[index+1] - fEnergy[index]);

		//this gives the cross section in barns/steradian
		if(fVerboseLevel > 2) {
			G4cout<<__PRETTY_FUNCTION__<<" "<<index<<": crossSection = "<<fCrossSection[index]<<" + "<<(fBeamEnergy - fEnergy[index])<<" * "<<(fCrossSection[index+1] - fCrossSection[index])<<"/"<<(fEnergy[index+1] - fEnergy[index])<<" = "<<fCrossSection[index] + (fBeamEnergy - fEnergy[index]) * (fCrossSection[index+1] - fCrossSection[index])/(fEnergy[index+1] - fEnergy[index])<<G4endl;
			G4cout<<__PRETTY_FUNCTION__<<" crossSection/(2.*CLHEP::pi)*CrossSection("<<angle/deg<<", coeff) = "<<crossSection<<"/("<<2.*CLHEP::pi<<")*"<<CrossSection(angle/rad,coeff)<<" = "<<crossSection/(2.*CLHEP::pi)*CrossSection(angle/rad,coeff)<<G4endl;
			G4cout<<__PRETTY_FUNCTION__<<": returning "<<crossSection/(2.*CLHEP::pi)*CrossSection(angle/rad,coeff)/conversion<<" = "<<crossSection<<"/(2.*CLHEP::pi)*"<<CrossSection(angle/rad,coeff)<<"/"<<conversion<<G4endl;
		}
		return crossSection/(2.*CLHEP::pi)*CrossSection(angle/rad,coeff)/conversion;//conversion is dOmega_lab/dOmega_cm
	} else {
		//variable different number of angle/probability pairs for each energy
		//so we first get the interpolated probability for this energy and the next energy and then do the interpolation between the energies
		size_t angInd1;
		for(angInd1 = 0; angInd1 < fCosThCm[index].size()-1; ++angInd1) {
			if(angle/rad < acos(fCosThCm[index][angInd1+1])) {
				break;
			}
		}
		double prob1 = fCoeff[index][angInd1] + (angle/rad - acos(fCosThCm[index][angInd1]))*(fCoeff[index][angInd1+1] - fCoeff[index][angInd1])/(acos(fCosThCm[index][angInd1+1]) - acos(fCosThCm[index][angInd1]));

		size_t angInd2;
		for(angInd2 = 0; angInd2 < fCosThCm[index+1].size()-1; ++angInd2) {
			if(angle/rad < acos(fCosThCm[index+1][angInd2+1])) {
				break;
			}
		}
		double prob2 = fCoeff[index+1][angInd2] + (angle/rad - acos(fCosThCm[index+1][angInd2]))*(fCoeff[index+1][angInd2+1] - fCoeff[index+1][angInd2])/(acos(fCosThCm[index+1][angInd2+1]) - acos(fCosThCm[index+1][angInd2]));

		if(fVerboseLevel > 3) {
			G4cout<<__PRETTY_FUNCTION__<<"index "<<index<<", angInd1 = "<<angInd1<<", angInd2 = "<<angInd2<<G4endl
				<<__PRETTY_FUNCTION__<<fEnergy[index]<<": "<<prob1<<" = "<<fCoeff[index][angInd1]<<" + "<<(angle/rad - acos(fCosThCm[index][angInd1]))*(fCoeff[index][angInd1+1] - fCoeff[index][angInd1])/(acos(fCosThCm[index][angInd1+1]) - acos(fCosThCm[index][angInd1]))<<" = "<<fCoeff[index][angInd1]<<" + "<<(angle/rad - acos(fCosThCm[index][angInd1]))<<"*"<<(fCoeff[index][angInd1+1] - fCoeff[index][angInd1])<<"/"<<(acos(fCosThCm[index][angInd1+1]) - acos(fCosThCm[index][angInd1]))<<G4endl
				<<__PRETTY_FUNCTION__<<fEnergy[index+1]<<": "<<prob2<<" = "<<fCoeff[index+1][angInd2]<<" + "<<(angle/rad - acos(fCosThCm[index+1][angInd2]))*(fCoeff[index+1][angInd2+1] - fCoeff[index+1][angInd2])/(acos(fCosThCm[index+1][angInd2+1]) - acos(fCosThCm[index+1][angInd2]))<<" = "<<fCoeff[index+1][angInd2]<<" + "<<(angle/rad - acos(fCosThCm[index+1][angInd2]))<<"*"<<(fCoeff[index+1][angInd2+1] - fCoeff[index+1][angInd2])<<"/"<<(acos(fCosThCm[index+1][angInd2+1]) - acos(fCosThCm[index+1][angInd2]))<<G4endl
				<<__PRETTY_FUNCTION__<<"=> "<<conversion*(prob1 + (fBeamEnergy - fEnergy[index]) * (prob2 - prob1)/(fEnergy[index+1] - fEnergy[index]))<<" = "<<conversion<<"*("<<prob1<<" + "<<(fBeamEnergy - fEnergy[index])<<" * "<<(prob2 - prob1)<<"/"<<(fEnergy[index+1] - fEnergy[index])<<")"<<G4endl;
			G4cout<<__PRETTY_FUNCTION__<<": returning "<<conversion*(prob1 + (fBeamEnergy - fEnergy[index]) * (prob2 - prob1)/(fEnergy[index+1] - fEnergy[index]))<<G4endl;
		}

		return conversion*(prob1 + (fBeamEnergy - fEnergy[index]) * (prob2 - prob1)/(fEnergy[index+1] - fEnergy[index]));
	}
}

//this function returns the cross section "per unit-cosine"
//to get an absolute value it has to be multiplied by the integrated cross section over two pi
double Kentucky::CrossSection(double angle, std::vector<double> coeff) const
{
	//reverse direction of theta (why?)
	angle = CLHEP::pi/rad - angle;

	double result = 0.5;

	for(size_t l = 0; l < coeff.size(); ++l) {
		result += (2.*(l+1.)+1.)/2.*coeff[l]*gsl_sf_legendre_Pl(l+1.,cos(angle/rad));
	}

	return result;
}

void Kentucky::Reaction(G4String reaction)
{
	if(fVerboseLevel > 0) {
		G4cout<<__PRETTY_FUNCTION__<<": '"<<reaction<<"'"<<G4endl;
	}
	reaction.toLower();
	if(reaction == "p,t") {
		fReaction = EReaction::kPT;
		if(!Open("PTCrossSections.dat")) {
			G4cout<<__PRETTY_FUNCTION__<<": failed to open 'PTCrossSections.dat'"<<G4endl;
		}
	} else if(reaction == "d,t") {
		fReaction = EReaction::kDT;
		if(!Open("DTCrossSections.dat")) {
			G4cout<<__PRETTY_FUNCTION__<<": failed to open 'DTCrossSections.dat'"<<G4endl;
		}
	} else if(reaction == "d,d") {
		fReaction = EReaction::kDD;
		if(!Open("DDCrossSections.dat")) {
			G4cout<<__PRETTY_FUNCTION__<<": failed to open 'DDCrossSections.dat'"<<G4endl;
		}
	} else {
		G4cerr<<"Unknown reaction '"<<reaction<<"', please use 'p,t', 'd,t', or 'd,d'. No whitespaces!"<<G4endl;
		return;
	}

	if(fVerboseLevel > 0) {
		G4cout<<__PRETTY_FUNCTION__<<": "<<fZeroDegreeEnergy<<", masses "<<fProjectileMass<<", "<<fTargetMass<<", "<<fEjectileMass<<", "<<fRecoilMass<<G4endl;
	}
	// if the energy was already set, set it again to calculate the beam energy
	if(fZeroDegreeEnergy > 0.) Energy(fZeroDegreeEnergy);
}

void Kentucky::Energy(double val)
{
	// this is the energy of the neutron at zero degree
	fZeroDegreeEnergy = val;
	if(fVerboseLevel > 0) {
		G4cout<<__PRETTY_FUNCTION__<<": "<<fZeroDegreeEnergy<<", masses "<<fProjectileMass<<", "<<fTargetMass<<", "<<fEjectileMass<<", "<<fRecoilMass<<G4endl;
	}
	// only proceed if the masses have been set (via Kentucky::Reaction)
	if(fProjectileMass == 0. || fTargetMass == 0. || fEjectileMass == 0. || fRecoilMass == 0.) return;
	// calculate energy of incoming beam
	fBeamEnergy = Physics::BeamEnergy(fProjectileMass, fTargetMass, fEjectileMass, fRecoilMass, fZeroDegreeEnergy);
	if(fVerboseLevel > 1) {
		G4cout<<__PRETTY_FUNCTION__<<" BeamEnergy("<<fProjectileMass<<", "<<fTargetMass<<", "<<fEjectileMass<<", "<<fRecoilMass<<", "<<fZeroDegreeEnergy<<") = "<<fBeamEnergy<<G4endl;
	}
	// check that the energy is in the range of the data file
	if(fBeamEnergy < fEnergy[0] || fBeamEnergy > fEnergy.back()) {
		if(fVerboseLevel > 0) {
			G4cout<<__PRETTY_FUNCTION__<<" energy "<<fBeamEnergy/CLHEP::keV<<" keV out of range "<<fEnergy[0]<<" - "<<fEnergy.back()<<G4endl;
		}
		if(fBeamEnergy < fEnergy[0]) {
			fBeamEnergy = fEnergy[0];
		}
		if(fBeamEnergy > fEnergy.back()) {
			fBeamEnergy = fEnergy.back();
		}
	}
	// calculate maximum cross section
	fMax = fCrossSectionFunction->GetMaximum(fMinimumTheta, fMaximumTheta);
	if(fVerboseLevel > 0) {
		G4cout<<__PRETTY_FUNCTION__<<" fMax = "<<fMax<<G4endl;
	}
	Print();
}

void Kentucky::DirectionAndEnergy(G4ParticleGun* particleGun) const
{
	/// this function calculates the direction and energy of the neutron
	// first calculate direction: phi is random from 0 - 2 pi, theta based on the cross section
	// quick and dirty method
	double phi = fMinimumPhi + (fMaximumPhi-fMinimumPhi)*G4UniformRand();
	double theta = fMinimumTheta + (fMaximumTheta-fMinimumTheta)*G4UniformRand();
	double y = fMax*G4UniformRand();
	int counter = 0;
	if(fVerboseLevel > 1) {
		G4cout<<__PRETTY_FUNCTION__<<": "<<counter<<" theta = "<<theta<<", y = "<<y<<", fMax = "<<fMax<<G4endl;
	}
	while(y > fCrossSectionFunction->Eval(theta)) {
		if(fVerboseLevel > 1) {
			G4cout<<__PRETTY_FUNCTION__<<": "<<counter<<" discarding theta "<<theta/deg<<" from cross section "<<fCrossSectionFunction->Eval(theta)<<G4endl;
		}
		theta = fMinimumTheta + (fMaximumTheta-fMinimumTheta)*G4UniformRand();
		y = fMax*G4UniformRand();
		++counter;
		if(fVerboseLevel > 1) {
			G4cout<<__PRETTY_FUNCTION__<<": "<<counter<<" theta = "<<theta<<", y = "<<y<<", fMax = "<<fMax<<G4endl;
		}
		if(counter > 1000) {
			G4cout<<__PRETTY_FUNCTION__<<": unable to find valid theta after "<<counter<<" iterations."<<G4endl;
			return;
		}
	}
	// create unity vector pointing along z-axis
	G4ThreeVector dir(0., 0., 1.);
	dir.rotateY(theta);
	dir.rotateZ(phi);
	particleGun->SetParticleMomentumDirection(dir);
	// now calculate the energy (excitation energy is zero for all reactions)
	if(fVerboseLevel > 0) {
		G4cout<<__PRETTY_FUNCTION__<<": calculating neutron energy from eex 0., theta "<<theta<<", beam energy "<<fBeamEnergy<<", and masses "<<fProjectileMass<<", "<<fTargetMass<<", "<<fEjectileMass<<", "<<fRecoilMass<<G4endl;
	}
	double energy = Physics::TRecoilLab(0., theta, fBeamEnergy, fProjectileMass, fTargetMass, fEjectileMass, fRecoilMass);
	if(fVerboseLevel > 0) {
		G4cout<<__PRETTY_FUNCTION__<<": energy "<<energy<<G4endl;
	}
	particleGun->SetParticleEnergy(energy);
}

void Kentucky::Print() const
{
	G4cout<<"Beam energy "<<fBeamEnergy<<" (zero degree energy "<<fZeroDegreeEnergy<<") and masses "<<fProjectileMass<<", "<<fTargetMass<<", "<<fEjectileMass<<", "<<fRecoilMass<<":"<<G4endl;
	for(double theta = 0.; theta < 180.*deg; theta += 1.*deg) {
		G4cout<<theta/deg<<"\t "<<fCrossSectionFunction->Eval(theta)<<G4endl;
	}
}
