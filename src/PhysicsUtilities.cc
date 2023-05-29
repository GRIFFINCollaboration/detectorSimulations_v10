#include "PhysicsUtilities.hh"

#include <cmath>
#include "globals.hh"
#include "G4SystemOfUnits.hh"

double Physics::bwe(int lam, int A, int outputflag)
{
  /* output */
	if(outputflag==1) {
		G4cout<<G4endl
			<<"bwe:"<<G4endl
			<<"*****"<<G4endl
			<<"calculates the Weisskopf single particle value for given massnumber"<<G4endl
			<<"and multipolarity for electric multipole transitions in e^2  fm^2lam"<<G4endl;
	}

	/* calculate Weisskopf single particle value */
	double bwe = -1;
	bwe = 1./(4.*(CLHEP::pi)) * pow((3./(double)(lam+3.)),2.) * pow(1.2,2.*lam) * pow(A,2.*lam/3.);

	/* output */
	if(outputflag==1) {
		G4cout<<G4endl
			<<"Weisskopf unit for A = "<<A<<" is: B_w(E"<<lam<<") = "<<bwe<<" e^2 fm^"<<2.*lam<<G4endl
			<<G4endl;
	}

	return bwe;
}

double Physics::SolidAngleConversionRecoil(double BetaCm, double RecoilBetaCm, double RecoilThetaCm)
{
	//calculate conversion of solid angle from lab to cm (dOmega_lab/dOmega_cm)
	//from heiko's phd thesis:
	//(dsigma/dOmega)_cm = (dsigma/dOmega)_lab * dOmega_lab/dOmega_cm
	//                   = (dsigma/dOmega)_lab * d(cos theta_lab)/d(cos theta_cm)
	//d(cos theta_lab)/d(cos theta_cm) = gamma*(1 + cos(theta_cm)*beta/beta_cm) / (gamma^2*(cos(theta_cm)+beta/beta_cm)^2 + (sin(theta_cm)^2))^3/2
	//beta_cm and theta_cm are those of the ejectile in the cm-frame and beta and gamma or those of the cm-frame in the lab-frame
	//but I want to calculate (dsigma/dOmega)_cm from (dsigma/dOmega)_lab of the recoil => switch ejectile with recoil
	//below everything w/o Recoil means the cm-frame
	RecoilThetaCm = CLHEP::pi - RecoilThetaCm;//different definitions of theta
	return SolidAngleConversionEjectile(BetaCm, RecoilBetaCm, RecoilThetaCm);
}

double Physics::SolidAngleConversionRecoil(double beamEnergy, double projectileMass, double targetMass, double ejectileMass, double recoilMass, double eEx, double recoilThetaCm)
{
	recoilThetaCm = CLHEP::pi - recoilThetaCm;//different definitions of theta
	return SolidAngleConversionEjectile(beamEnergy, projectileMass, targetMass, ejectileMass, recoilMass, eEx, recoilThetaCm);
}

double Physics::SolidAngleConversionEjectile(double BetaCm, double ejectileBetaCm, double ejectileThetaCm)
{
	//calculate conversion of solid angle from lab to cm (dOmega_lab/dOmega_cm)
	//from heiko's phd thesis:
	//(dsigma/dOmega)_cm = (dsigma/dOmega)_lab * dOmega_lab/dOmega_cm
	//                   = (dsigma/dOmega)_lab * d(cos theta_lab)/d(cos theta_cm)
	//d(cos theta_lab)/d(cos theta_cm) = gamma*(1 + cos(theta_cm)*beta/beta_cm) / (gamma^2*(cos(theta_cm)+beta/beta_cm)^2 + (sin(theta_cm)^2))^3/2
	//beta_cm and theta_cm are those of the ejectile in the cm-frame and beta and gamma or those of the cm-frame in the lab-frame
	//below everything w/o ejectile means the cm-frame
	double GammaCmSquared = 1/(1 - pow(BetaCm,2));

	double Denominator = pow(GammaCmSquared*pow(cos(ejectileThetaCm)+BetaCm/ejectileBetaCm,2)+pow(sin(ejectileThetaCm),2),3./2.);
	//double Denominator = pow(pow(cos(ejectileThetaCm)+BetaCm/ejectileBetaCm,2)+pow(sin(ejectileThetaCm),2),3./2.);

	if(Denominator != 0)
		return sqrt(GammaCmSquared) * (1 + cos(ejectileThetaCm)*BetaCm/ejectileBetaCm) / Denominator;
	else
		return -1.;
}

double Physics::SolidAngleConversionEjectile(double beamEnergy, double projectileMass, double targetMass, double ejectileMass, double recoilMass, double eEx, double ejectileThetaCm)
{
	double betaCm = BetaCm(beamEnergy, projectileMass, targetMass);
	double eCm = ECm(beamEnergy, projectileMass, targetMass);
	double ejectileBetaCm = BetaEjectileCm(eCm, eEx, ejectileMass, recoilMass);
	//calculate conversion of solid angle from lab to cm (dOmega_cm/dOmega_lab)
	//from heiko's phd thesis:
	//(dsigma/dOmega)_cm = (dsigma/dOmega)_lab * dOmega_lab/dOmega_cm
	//                   = (dsigma/dOmega)_lab * d(cos theta_lab)/d(cos theta_cm)
	//d(cos theta_lab)/d(cos theta_cm) = gamma*(1 + cos(theta_cm)*beta/beta_cm) / (gamma^2*(cos(theta_cm)+beta/beta_cm)^2 + (sin(theta_cm)^2))^3/2
	//beta_cm and theta_cm are those of the ejectile in the cm-frame and beta and gamma or those of the cm-frame in the lab-frame
	//below everything w/o ejectile means the cm-frame
	double gammaCmSquared = 1/(1 - pow(betaCm,2));

	double denominator = pow(gammaCmSquared*pow(cos(ejectileThetaCm)+betaCm/ejectileBetaCm,2)+pow(sin(ejectileThetaCm),2),3./2.);

	if(denominator != 0)
		return sqrt(gammaCmSquared) * (1 + cos(ejectileThetaCm)*betaCm/ejectileBetaCm) / denominator;
	else
		return -1.;
}

double Physics::SolidAngleConversionRecoil(double* x, double* par)
{
	//calculate conversion of solid angle from lab to cm (dOmega_cm/dOmega_lab)
	//parameters: 0 = beta of cm system, 1 = beta of recoil in cm system

	return SolidAngleConversionRecoil(par[0], par[1], x[0]);
}

double Physics::SolidAngleConversionEjectile(double* x, double* par)
{
	//calculate conversion of solid angle from lab to cm (dOmega_cm/dOmega_lab)
	//parameters: 0 = beta of cm system, 1 = beta of recoil in cm system

	return SolidAngleConversionEjectile(par[0], par[1], x[0]);
}


//WARNING: theta_lab is 180 degree for zero degree cm scattering (=> recoil)
double Physics::ThetaLab(double BetaCm, double RecoilBetaCm, double RecoilThetaCm)
{
	//calculate theta_lab
	double CosineThetaLab = 0.;

	//from heiko's phd thesis:
	//cos(theta_lab) = gamma*(cos(theta_cm)+beta/beta_cm)/sqrt(sin^2(theta_cm)+gamma^2*(cos(theta_cm)+beta/beta_cm)^2)
	//beta_cm and theta_cm are those of the ejectile in the cm-frame and beta and gamma or those of the cm-frame in the lab-frame
	double GammaCmSquared = 1/(1 - pow(BetaCm,2));

	double CosinePlusBetaRatio = cos(RecoilThetaCm)+BetaCm/RecoilBetaCm;

	double Denominator = sqrt(pow(sin(RecoilThetaCm),2)+GammaCmSquared*pow((CosinePlusBetaRatio),2));

	if(Denominator != 0)
		CosineThetaLab = sqrt(GammaCmSquared)*(CosinePlusBetaRatio) / Denominator;
	else
		return -1.;

	return CLHEP::pi - acos(CosineThetaLab);
}

double Physics::ThetaLab(double beamEnergy, double projectileMass, double targetMass, double ejectileMass, double recoilMass, double eEx, double recoilThetaCm)
{
	double betaCm = BetaCm(beamEnergy, projectileMass, targetMass);
	double eCm = ECm(beamEnergy, projectileMass, targetMass);
	double recoilBetaCm = BetaRecoilCm(eCm, eEx, ejectileMass, recoilMass);
	//calculate theta_lab
	double cosineThetaLab = 0.;

	//from heiko's phd thesis:
	//cos(theta_lab) = gamma*(cos(theta_cm)+beta/beta_cm)/sqrt(sin^2(theta_cm)+gamma^2*(cos(theta_cm)+beta/beta_cm)^2)
	//beta_cm and theta_cm are those of the ejectile in the cm-frame and beta and gamma or those of the cm-frame in the lab-frame
	double gammaCmSquared = 1/(1 - pow(betaCm,2));

	double cosinePlusBetaRatio = cos(recoilThetaCm)+betaCm/recoilBetaCm;

	double denominator = sqrt(pow(sin(recoilThetaCm),2)+gammaCmSquared*pow((cosinePlusBetaRatio),2));

	if(denominator != 0)
		cosineThetaLab = sqrt(gammaCmSquared)*(cosinePlusBetaRatio)/denominator;
	else
		return -1.;

	return CLHEP::pi - acos(cosineThetaLab);
}

//what's the difference to ThetaRecoilCm??? this one gives unreasonable results???
double Physics::ThetaCm(double BetaCm, double RecoilBetaCm, double RecoilThetaLab)
{
	//double Resolution = 1.0e-3;
	////calculate theta_cm by recursion
	//double BetaRatio = BetaCm/RecoilBetaCm;
	//double GammaCmSquared = 1/(1 - pow(BetaCm,2));
	//
	//double ThetaLabMinimum = 0.;
	//
	////check whether the given theta_lab can be reached with the given beta ratio
	//if(BetaRatio == 1)//maximum theta_lab is 90 degree
	//  {
	//    ThetaLabMinimum = CLHEP::pi/2.;
	//  }
	//else if(BetaRatio > 1)
	//  {
	//    ThetaLabMinimum = CLHEP::pi - atan2(RecoilBetaCm,sqrt(GammaCmSquared*(pow(BetaCm,2)-pow(RecoilBetaCm,2))));
	//
	//    G4cerr<<"Warning, theta_cm is not unambiguous, will only return one value!"<<G4endl;
	//  }
	//
	//if(RecoilThetaLab < ThetaLabMinimum)
	//  {
	//    G4cerr<<__PRETTY_FUNCTION__<<": this theta_lab can't be reached with the given beta's"<<G4endl;
	//
	//    return -1;
	//  }
	//
	//double ThetaCmLow = 0.;
	//double ThetaCmHigh = CLHEP::pi;
	//double ThetaCm = CLHEP::pi/2.;
	//
	//while(abs(ThetaLab(BetaCm, RecoilBetaCm, ThetaCm) - RecoilThetaLab) > Resolution)
	//  {
	//    //check which two value pairs bracket the desired theta_lab, ThetaCm and ThetaCmLow or ThetaCm and ThetaCmHigh
	//    if((ThetaLab(BetaCm, RecoilBetaCm, ThetaCm) - RecoilThetaLab)*(ThetaLab(BetaCm, RecoilBetaCm, ThetaCmLow) - RecoilThetaLab) < 0)
	//	{
	//	  ThetaCmHigh = ThetaCm;
	//	  ThetaCm = (ThetaCmLow + ThetaCmHigh)/2.;
	//	}
	//    else if((ThetaLab(BetaCm, RecoilBetaCm, ThetaCm) - RecoilThetaLab)*(ThetaLab(BetaCm, RecoilBetaCm, ThetaCmHigh) - RecoilThetaLab) < 0)
	//	{
	//	  ThetaCmLow = ThetaCm;
	//	  ThetaCm = (ThetaCmLow + ThetaCmHigh)/2.;
	//	}
	//    else
	//	{
	//	  G4cerr<<"bracketing method failed, return average of high and low bracket"<<G4endl;
	//
	//	  return (ThetaLab(BetaCm, RecoilBetaCm, ThetaCmLow) + ThetaLab(BetaCm, RecoilBetaCm, ThetaCmHigh))/2.;
	//	}
	//  }
	//
	//return ThetaLab(BetaCm, RecoilBetaCm, ThetaCm);

	//x4 = betacm/betarecoil
	//thetacm = acos((-x4*gammacm^2*tan(thetalab_recoil)^2+-sqrt(1+gammacm^2*tan(thetalab_recoil)^2*(1-x4^2)))/(1+gammacm^2*tan(thetalab_recoil)^2));
	double BetaRatio = BetaCm/RecoilBetaCm;
	double GammaCmSquared = 1/(1 - pow(BetaCm,2));
	double GammaTangensSquared = GammaCmSquared*pow(RecoilThetaLab,2);

	//check whether the denominator will be zero
	if(GammaTangensSquared == -1) {
		G4cerr<<__PRETTY_FUNCTION__<<": GammaTangensSquared is -1"<<G4endl;
		return -1.;
	}

	double CosineThetaCm;

	//WARNING: this is for recoils i.e. theta_lab = 180 degree => theta_cm = 0 degree!
	if(RecoilThetaLab > CLHEP::pi/2.) {
		CosineThetaCm = (-BetaRatio*GammaTangensSquared+sqrt(1+GammaTangensSquared*(1-pow(BetaRatio,2)))) / (1+GammaTangensSquared);
	} else {
		CosineThetaCm = (-BetaRatio*GammaTangensSquared-sqrt(1+GammaTangensSquared*(1-pow(BetaRatio,2)))) / (1+GammaTangensSquared);
	}

	return acos(CosineThetaCm);
}

double Physics::RelativisticGamma(double Beta)
{
	if(Beta < 0 || Beta > 1) {
		G4cerr<<__PRETTY_FUNCTION__<<": Beta is not in range 0-1: "<<Beta<<G4endl;
		return -1;
	}

	return 1./sqrt(1.-pow(Beta,2));
}

//total kinetic energy in the ingoing channel of the cm-system
double Physics::TiCm(double BeamEnergy, double ProjectileMass, double TargetMass)
{
	return sqrt(2.*(ProjectileMass+BeamEnergy)*TargetMass + pow(ProjectileMass, 2) + pow(TargetMass, 2)) - ProjectileMass - TargetMass;
}

//total kinetic energy in the outgoing channel of the cm-system
double Physics::TfCm(double ECm, double Eex, double EjectileMass, double RecoilMass)
{
	return ECm - Eex - (EjectileMass+RecoilMass);
}

//kinetic energy of the projectile in the cm-system
double Physics::TProjectileCm(double BeamEnergy, double ProjectileMass, double TargetMass)
{
	if(ECm(BeamEnergy, ProjectileMass, TargetMass) > 0.) {
		return TiCm(BeamEnergy, ProjectileMass, TargetMass)/2.*(TiCm(BeamEnergy, ProjectileMass, TargetMass) + 2.*TargetMass)/ECm(BeamEnergy, ProjectileMass, TargetMass);
	}
	return 0.;
}

//kinetic energy of the recoil in the cm-system
double Physics::TRecoilCm(double ECm, double Eex, double EjectileMass, double RecoilMass)
{
	double tfCm = TfCm(ECm,Eex,EjectileMass,RecoilMass);

	if(ECm == 0.) {
		G4cerr<<__PRETTY_FUNCTION__<<": ECm is 0"<<G4endl;
		return -1.;
	}

	//cout<<"tfCm = "<<tfCm<<", EjectileMass = "<<EjectileMass<<", Eex = "<<Eex<<", ECm = "<<ECm<<" => TRecoilCm = "<<tfCm*(tfCm+2*(EjectileMass+Eex))/(2.*ECm)<<G4endl;

	return tfCm*(tfCm+2.*(EjectileMass+Eex))/(2.*ECm);
}

//kinetic energy of the ejectile in the cm-system
double Physics::TEjectileCm(double ECm, double Eex, double EjectileMass, double RecoilMass)
{
	double tfCm = TfCm(ECm,Eex,EjectileMass,RecoilMass);

	if(ECm == 0.) {
		G4cerr<<__PRETTY_FUNCTION__<<": ECm is 0"<<G4endl;
		return -1.;
	}

	//cout<<"tfCm = "<<tfCm<<", EjectileMass = "<<EjectileMass<<", Eex = "<<Eex<<", ECm = "<<ECm<<" => TRecoilCm = "<<tfCm*(tfCm+2*(EjectileMass+Eex))/(2.*ECm)<<G4endl;

	return tfCm*(tfCm+2.*(RecoilMass+Eex))/(2.*ECm);
}

//beta of the recoil in the cm-system
double Physics::BetaRecoilCm(double ECm, double Eex, double EjectileMass, double RecoilMass)
{
	double tRecoilCm = TRecoilCm(ECm,Eex,EjectileMass,RecoilMass);

	//check that neither denominator gets zero or sqrt gets imaginary
	//this is fulfilled if trecoilcm and RecoilMass are larger than zero (which they have to be anyway)
	if(tRecoilCm < 0 || RecoilMass < 0) {
		G4cerr<<__PRETTY_FUNCTION__<<": either the recoil mass ("<<RecoilMass<<") or the recoil kinetic energy in the cm system ("<<tRecoilCm<<") are 0"<<G4endl;
		return -1;
	}

	return sqrt(pow(tRecoilCm,2) + 2.*tRecoilCm*RecoilMass)/(tRecoilCm + RecoilMass);
}

//beta of the ejectile in the cm-system
double Physics::BetaEjectileCm(double ECm, double Eex, double EjectileMass, double RecoilMass)
{
	double tEjectileCm = TEjectileCm(ECm,Eex,EjectileMass,RecoilMass);

	//check that neither denominator gets zero or sqrt gets imaginary
	//this is fulfilled if tEjectileCm and EjectileMass are larger than zero (which they have to be anyway)
	if(tEjectileCm < 0 || EjectileMass < 0) {
		G4cerr<<__PRETTY_FUNCTION__<<": either the ejectile mass ("<<EjectileMass<<") or the ejectile kinetic energy in the cm system ("<<tEjectileCm<<") are 0"<<G4endl;
		return -1;
	}

	return sqrt(pow(tEjectileCm,2) + 2.*tEjectileCm*EjectileMass)/(tEjectileCm + EjectileMass);
}

//theta of recoil in cm-system, Beta = beta of cm-system, Eex = excitation energy, WORKS ONLY IF THE PARTILCE CAN REACH GIVEN LAB ANGLE (180 DEGREE)!!!
double Physics::ThetaRecoilCm(double Beta, double ECm, double Eex, double ThetaRecoilLab, double EjectileMass, double RecoilMass)
{
	double BetaRatio = Beta/BetaRecoilCm(ECm,Eex,EjectileMass,RecoilMass);
	double GammaTangensSquared = pow(tan(ThetaRecoilLab),2)/(1.-pow(Beta,2));

	//check whether the denominator will be zero or sqrt gets imaginary
	if(GammaTangensSquared == -1 || (1.+GammaTangensSquared*(1.-pow(BetaRatio,2))) < 0) {
		//G4cerr<<__PRETTY_FUNCTION__<<": either the denominator gets zero or the sqrt will be imaginary: "<<GammaTangensSquared<<", "<<(1.+GammaTangensSquared*(1.-pow(Beta,2)))<<G4endl;
		return -1.;
	}

	//this seems to be the right combination, but why?
	if(ThetaRecoilLab < CLHEP::pi/2.)
		return  acos((BetaRatio*GammaTangensSquared - sqrt(1.+GammaTangensSquared*(1.-pow(BetaRatio,2))))/(1+GammaTangensSquared));

	return acos((BetaRatio*GammaTangensSquared + sqrt(1.+GammaTangensSquared*(1.-pow(BetaRatio,2))))/(1+GammaTangensSquared));
}

//cosine of theta of recoil in cm-system, Beta = beta of cm-system, Eex = excitation energy, WORKS ONLY IF THE PARTICLE CAN REACH GIVEN LAB ANGLE (180 DEGREE)!!!
double Physics::CosineThetaRecoilCm(double Beta, double ECm, double Eex, double ThetaRecoilLab, double EjectileMass, double RecoilMass)
{
	double BetaRatio = Beta/BetaRecoilCm(ECm,Eex,EjectileMass,RecoilMass);
	double GammaTangensSquared = pow(tan(ThetaRecoilLab),2)/(1.-pow(Beta,2));

	//check whether the denominator will be zero or sqrt gets imaginary
	if(GammaTangensSquared == -1) {
		G4cerr<<__PRETTY_FUNCTION__<<": the denominator will be zero: "<<GammaTangensSquared<<G4endl;
		return -2.;
	}

	if((1.+GammaTangensSquared*(1.-pow(BetaRatio,2))) < 0) {
		//G4cerr<<__PRETTY_FUNCTION__<<": the sqrt will be imaginary: "<<(1.+GammaTangensSquared*(1.-pow(BetaRatio,2)))<<" (1+"<<GammaTangensSquared<<"*(1.-pow("<<BetaRatio<<",2)), Beta = "<<Beta<<")"<<G4endl;
		return -2.;
	}

	//G4cout<<__PRETTY_FUNCTION__<<": the sqrt will be real: "<<(1.+GammaTangensSquared*(1.-pow(BetaRatio,2)))<<" (1+"<<GammaTangensSquared<<"*(1.-pow("<<BetaRatio<<",2)), Beta = "<<Beta<<")"<<G4endl;

	//this seems to be the right combination, but why?
	if(ThetaRecoilLab < CLHEP::pi/2.) {
		return (BetaRatio*GammaTangensSquared - sqrt(1.+GammaTangensSquared*(1.-pow(BetaRatio,2))))/(1+GammaTangensSquared);
	}

	return (BetaRatio*GammaTangensSquared + sqrt(1.+GammaTangensSquared*(1.-pow(BetaRatio,2))))/(1+GammaTangensSquared);
}

//cosine of theta of ejectile in cm-system, Beta = beta of cm-system, Eex = excitation energy, WORKS ONLY IF THE PARTILCE CAN REACH GIVEN LAB ANGLE (180 DEGREE)!!!
double Physics::CosineThetaEjectileCm(double Beta, double ECm, double Eex, double ThetaEjectileLab, double EjectileMass, double RecoilMass)
{
	double BetaRatio = Beta/BetaEjectileCm(ECm,Eex,EjectileMass,RecoilMass);
	double GammaTangensSquared = pow(tan(ThetaEjectileLab),2)/(1.-pow(Beta,2));

	//check whether the denominator will be zero or sqrt gets imaginary
	if(GammaTangensSquared == -1) {
		G4cerr<<__PRETTY_FUNCTION__<<": the denominator will be zero: "<<GammaTangensSquared<<G4endl;
		return -2.;
	}

	if((1.+GammaTangensSquared*(1.-pow(BetaRatio,2))) < 0) {
		//G4cerr<<__PRETTY_FUNCTION__<<": the sqrt will be imaginary: "<<(1.+GammaTangensSquared*(1.-pow(BetaRatio,2)))<<" (1+"<<GammaTangensSquared<<"*(1.-pow("<<BetaRatio<<",2)), Beta = "<<Beta<<")"<<G4endl;
		return -2.;
	}

	//this seems to be the right combination, but why?
	if(ThetaEjectileLab < CLHEP::pi/2.) {
		return (BetaRatio*GammaTangensSquared - sqrt(1.+GammaTangensSquared*(1.-pow(BetaRatio,2))))/(1+GammaTangensSquared);
	}

	return (BetaRatio*GammaTangensSquared + sqrt(1.+GammaTangensSquared*(1.-pow(BetaRatio,2))))/(1+GammaTangensSquared);
}

//cosine of theta of ejectile in cm-system, Beta = beta of cm-system, Eex = excitation energy, WORKS ONLY IF THE PARTILCE CAN REACH GIVEN LAB ANGLE (180 DEGREE)!!!
double Physics::CosineThetaEjectileCm(double beamEnergy, double projectileMass, double targetMass, double ejectileMass, double recoilMass, double eEx, double thetaEjectileLab)
{
	double betaCm = BetaCm(beamEnergy, projectileMass, targetMass);
	double eCm = ECm(beamEnergy, projectileMass, targetMass);
	double betaRatio = betaCm/BetaEjectileCm(eCm, eEx, ejectileMass, recoilMass);
	double gammaTangensSquared = pow(tan(thetaEjectileLab),2)/(1.-pow(betaCm,2));

	//check whether the denominator will be zero or sqrt gets imaginary
	if(gammaTangensSquared == -1) {
		G4cerr<<__PRETTY_FUNCTION__<<": the denominator will be zero: "<<gammaTangensSquared<<G4endl;
		return -2.;
	}

	if((1.+gammaTangensSquared*(1.-pow(betaRatio,2))) < 0) {
		//G4cerr<<__PRETTY_FUNCTION__<<": the sqrt will be imaginary: "<<(1.+gammaTangensSquared*(1.-pow(betaRatio,2)))<<" (1+"<<gammaTangensSquared<<"*(1.-pow("<<betaRatio<<",2)), Beta = "<<betaCm<<")"<<G4endl;
		return -2.;
	}

	//this seems to be the right combination, but why?
	if(thetaEjectileLab < CLHEP::pi/2.) {
		return (betaRatio*gammaTangensSquared - sqrt(1.+gammaTangensSquared*(1.-pow(betaRatio,2))))/(1+gammaTangensSquared);
	}

	return (betaRatio*gammaTangensSquared + sqrt(1.+gammaTangensSquared*(1.-pow(betaRatio,2))))/(1+gammaTangensSquared);
}

//beta of recoil in cm-system for given lab and cm-angle and beta of cm-system
double Physics::BetaRecoilCm(double Beta, double ThetaRecoilLab, double ThetaRecoilCm)
{
	//second solution: plus sign in numerator (does not give the right results)
	double Denominator = 1-pow(Beta,2)-pow(cos(ThetaRecoilCm),2)*(1-pow(Beta,2) + pow(tan(ThetaRecoilLab),2));

	//check whether the denominator will be zero
	if(Denominator == -1 || Denominator == 0 || Beta < 0 || Beta > 1) {
		G4cerr<<__PRETTY_FUNCTION__<<": Beta is not in range 0-1: "<<Beta<<" or denominator is wrong: "<<Denominator<<G4endl;
		return -1.;
	}

	return (Beta*cos(ThetaRecoilCm)*pow(tan(ThetaRecoilLab),2) - sin(ThetaRecoilCm)*tan(ThetaRecoilLab)*Beta*sqrt(1-pow(Beta,2))) / Denominator;
}

//excitation energy needed to reach certain lab and cm-angle combinations
double Physics::Eex(double ECm, double Beta, double EjectileMass, double RecoilMass, double ThetaRecoilLab, double ThetaRecoilCm)
{
	//alternative solutions: minus and plus signs before the two Sqrt terms (i.e. four combinations)
	//and of course the second solutions from BetaRecoilCm => total eight different results
	//this one looks best in gnuplot (somewhat matches orbits output) and this and the other BetaRecoilCm solution give the only reasonable Eex-values (i.e. neither around 1e6 nor around -1e7)

	if((1-pow(BetaRecoilCm(Beta, ThetaRecoilLab, ThetaRecoilCm),2)) < 0) {
		G4cerr<<__PRETTY_FUNCTION__<<": beta of recoil in cm is larger than 1: "<<BetaRecoilCm(Beta, ThetaRecoilLab, ThetaRecoilCm)<<G4endl;
		return -1.;
	}

	if((pow(ECm,2) + pow(RecoilMass,2) - 2*ECm*RecoilMass/sqrt(1-pow(BetaRecoilCm(Beta, ThetaRecoilLab, ThetaRecoilCm),2))) < 0) {
		G4cerr<<__PRETTY_FUNCTION__<<": something is wrong"<<G4endl;
		return -1.;
	}

	return (- EjectileMass + sqrt(pow(ECm,2) + pow(RecoilMass,2) - 2*ECm*RecoilMass/sqrt(1-pow(BetaRecoilCm(Beta, ThetaRecoilLab, ThetaRecoilCm),2))));
}

//theta of recoil in lab-system, Beta = beta of cm-system, Eex = excitation energy
double Physics::ThetaRecoilLab(double Beta, double ECm, double Eex, double ThetaRecoilCm, double EjectileMass, double RecoilMass)
{
	double BetaRatio = Beta/BetaRecoilCm(ECm,Eex,EjectileMass,RecoilMass);

	return atan((sqrt(1.-pow(Beta,2))*sin(ThetaRecoilCm))/(cos(ThetaRecoilCm)+BetaRatio));
}

double Physics::TRecoilLab(double Eex, double ThetaRecoilLab, double BeamEnergy, double ProjectileMass, double TargetMass, double EjectileMass, double RecoilMass)
{
	double betaCm = BetaCm(BeamEnergy, ProjectileMass, TargetMass);
	double eCm = ECm(BeamEnergy, ProjectileMass, TargetMass);
	//G4cout<<__PRETTY_FUNCTION__<<": betaCm "<<betaCm<<", eCm "<<eCm<<G4endl;
	return TRecoilLab(betaCm, eCm, Eex, ThetaRecoilLab, EjectileMass, RecoilMass);
}

//kinetic energy of recoil in the lab system
double Physics::TRecoilLab(double Beta, double ECm, double Eex, double ThetaRecoilLab, double EjectileMass, double RecoilMass)
{
	double gamma = RelativisticGamma(Beta);
	double tRecoilCm = TRecoilCm(ECm,Eex,EjectileMass,RecoilMass);
	double cosThetaRecoilCm = CosineThetaRecoilCm(Beta, ECm, Eex, ThetaRecoilLab, EjectileMass, RecoilMass);

	//G4cout<<__PRETTY_FUNCTION__<<": gamma "<<gamma<<", tRecoilCm "<<tRecoilCm<<", cosThetaRecoilCm "<<cosThetaRecoilCm<<G4endl;
	//this happens only if the computation of ThetaRecoilCm failed
	if(cosThetaRecoilCm == -2.) {
		return 0.;
	}

	if(tRecoilCm < 0. || RecoilMass < 0. || cosThetaRecoilCm < -1. || cosThetaRecoilCm > 1.) {
		G4cerr<<__PRETTY_FUNCTION__<<": either the recoil mass ("<<RecoilMass<<") or the recoil kinetic energy in the cm system ("<<tRecoilCm<<") are 0 or the computation of theta of recoil in cm failed: "<<cosThetaRecoilCm<<G4endl;
		return -1;
	}

	return (gamma-1)*RecoilMass + gamma*tRecoilCm - gamma*Beta*sqrt(tRecoilCm*(tRecoilCm+2.*RecoilMass))*cosThetaRecoilCm;
}

//kinetic energy of recoil in the lab system
double Physics::TRecoilLab2(double Beta, double ECm, double Eex, double ThetaRecoilCm, double EjectileMass, double RecoilMass)
{
	double gamma = RelativisticGamma(Beta);
	double tRecoilCm = TRecoilCm(ECm,Eex,EjectileMass,RecoilMass);
	double cosThetaRecoilCm = cos(ThetaRecoilCm);

	if(tRecoilCm < 0. || RecoilMass < 0. || cosThetaRecoilCm < -1. || cosThetaRecoilCm > 1.) {
		G4cerr<<__PRETTY_FUNCTION__<<": either the recoil mass ("<<RecoilMass<<") or the recoil kinetic energy in the cm system ("<<tRecoilCm<<") are 0 or the computation of theta of recoil in cm failed: "<<cosThetaRecoilCm<<G4endl;
		return -1;
	}

	return (gamma-1)*RecoilMass + gamma*tRecoilCm - gamma*Beta*sqrt(tRecoilCm*(tRecoilCm+2.*RecoilMass))*cosThetaRecoilCm;
}

//kinetic energy of recoil in the lab system - doesn't work???
double Physics::TEjectileLab(double Beta, double ECm, double Eex, double ThetaEjectileLab, double EjectileMass, double RecoilMass)
{
	double gamma = RelativisticGamma(Beta);
	double tEjectileCm = TEjectileCm(ECm,Eex,EjectileMass,RecoilMass);
	double cosThetaEjectileCm = CosineThetaEjectileCm(Beta, ECm, Eex, ThetaEjectileLab, EjectileMass, RecoilMass);

	//this happens only if the computation of ThetaEjectileCm failed
	if(cosThetaEjectileCm == -2.) {
		return 0.;
	}

	if(tEjectileCm < 0. || RecoilMass < 0. || cosThetaEjectileCm < -1. || cosThetaEjectileCm > 1.) {
		G4cerr<<__PRETTY_FUNCTION__<<": either the recoil mass ("<<RecoilMass<<") or the recoil kinetic energy in the cm system ("<<tEjectileCm<<") are 0 or the computation of theta of recoil in cm failed: "<<cosThetaEjectileCm<<G4endl;
		return -1;
	}

	return (gamma-1)*RecoilMass + gamma*tEjectileCm - gamma*Beta*sqrt(tEjectileCm*(tEjectileCm+2.*RecoilMass))*cosThetaEjectileCm;
}

//energy in the cm-system
double Physics::ECm(double BeamEnergy, double ProjectileMass, double TargetMass)
{
	return sqrt(2.*(ProjectileMass+BeamEnergy)*TargetMass + pow(ProjectileMass, 2) + pow(TargetMass, 2));
}

//beta of cm-system
double Physics::BetaCm(double BeamEnergy, double ProjectileMass, double TargetMass)
{
	return sqrt(2*ProjectileMass*BeamEnergy + pow(BeamEnergy,2))/(ProjectileMass + TargetMass + BeamEnergy);
}

double Physics::BeamEnergy(double ProjectileMass, double TargetMass, double LabAngle, double LabEnergy)
{
	//LabAngle has to be in rad, everything else (masses and energy) in the same unit as the result

	//from mathematica (fourth out of four solution, first two give imaginary results, third is always negative)
	//(-8*m1 - 4*m2 + t4Lab*pow(Sec(aLab),2) + Sec(aLab)*sqrt(16*pow(m2,2)*pow(cos(aLab),2) + t4Lab*(t4Lab*pow(Sec(aLab),2)*pow(cos(aLab) - sin(aLab),4) + 8*m2*(3 + sin(2*aLab)))) - 2*t4Lab*tan(aLab) + sqrt((pow(Sec(aLab),4)*(24*pow(m2,2)*t4Lab*pow(cos(aLab),2) + 32*pow(m1,2)*m2*pow(cos(aLab),4) + 2*m2*pow(t4Lab,2)*pow(cos(aLab) - sin(aLab),4) + 16*pow(m2,2)*t4Lab*pow(cos(aLab),3)*sin(aLab) - 8*pow(m1,2)*t4Lab*pow(cos(aLab),2)*(-1 + sin(2*aLab)) + 2*pow(m1,2)*cos(3*aLab)*sqrt(2*(4*pow(m2,2) + 12*m2*t4Lab + pow(t4Lab,2)) + (8*pow(m2,2) - 2*pow(t4Lab,2))*cos(2*aLab) + t4Lab*(t4Lab*pow(Sec(aLab),2) + 8*m2*sin(2*aLab) - 4*t4Lab*tan(aLab))) + 2*cos(aLab)*(3*pow(m1,2) + m2*t4Lab - m2*t4Lab*sin(2*aLab))*sqrt(2*(4*pow(m2,2) + 12*m2*t4Lab + pow(t4Lab,2)) + (8*pow(m2,2) - 2*pow(t4Lab,2))*cos(2*aLab) + t4Lab*(t4Lab*pow(Sec(aLab),2) + 8*m2*sin(2*aLab) - 4*t4Lab*tan(aLab)))))/m2))/8.

	double t2 = pow(TargetMass,2);
	double p2 = pow(ProjectileMass,2);
	double c2 = pow(cos(LabAngle),2);
	double e2 = pow(LabEnergy,2);
	double te = TargetMass*LabEnergy;

	return (-8.*ProjectileMass - 4.*TargetMass + LabEnergy/c2 - 2.*LabEnergy*tan(LabAngle) 
			+ sqrt(16.*t2*c2 + LabEnergy*(LabEnergy/c2*pow(cos(LabAngle) - sin(LabAngle),4) + 8.*TargetMass*(3. + sin(2.*LabAngle))))/cos(LabAngle) 
			+ sqrt((24.*t2*LabEnergy*c2 + 
					32.*p2*TargetMass*pow(cos(LabAngle),4) + 
					2.*TargetMass*e2*pow(cos(LabAngle) - sin(LabAngle),4) + 
					16.*t2*LabEnergy*pow(cos(LabAngle),3)*sin(LabAngle) - 
					8.*p2*LabEnergy*c2*(sin(2.*LabAngle) - 1) + 
					2.*p2*cos(3.*LabAngle)*sqrt(2.*(4.*t2 + 12.*te + e2) + 
						(8.*t2 - 2.*e2)*cos(2.*LabAngle) + 
						LabEnergy*(LabEnergy/c2 + 8.*TargetMass*sin(2.*LabAngle) - 4.*LabEnergy*tan(LabAngle))) + 
					2.*cos(LabAngle)*(3.*p2 + te - te*sin(2.*LabAngle))*sqrt(2.*(4.*t2 + 12.*te + e2) + (8.*t2 - 2.*e2)*cos(2.*LabAngle) + LabEnergy*(LabEnergy/c2 + 8.*TargetMass*sin(2.*LabAngle) - 4.*LabEnergy*tan(LabAngle))))/(pow(cos(LabAngle),4)*TargetMass)))/8.;
}

//-------------------- no excitation
//total kinetic energy in the outgoing channel of the cm-system
double Physics::TfCm(double ECm, double EjectileMass, double RecoilMass)
{
	return ECm - (EjectileMass+RecoilMass);
}

//theta of recoil in cm-system, Beta = beta of cm-system, WORKS ONLY IF THE PARTILCE CAN REACH GIVEN LAB ANGLE (180 DEGREE)!!!
double Physics::ThetaRecoilCm(double Beta, double ECm, double ThetaRecoilLab, double EjectileMass, double RecoilMass)
{
	double BetaRatio = Beta/BetaRecoilCm(ECm,0.,EjectileMass,RecoilMass);
	double GammaTangensSquared = pow(tan(ThetaRecoilLab),2)/(1.-pow(Beta,2));

	//check whether the denominator will be zero or sqrt gets imaginary
	if(GammaTangensSquared == -1 || (1.+GammaTangensSquared*(1.-pow(BetaRatio,2))) < 0) {
		//G4cerr<<__PRETTY_FUNCTION__<<": either the denominator gets zero or the sqrt will be imaginary: "<<GammaTangensSquared<<", "<<(1.+GammaTangensSquared*(1.-pow(Beta,2)))<<G4endl;
		return -1.;
	}

	//this seems to be the right combination, but why?
	if(ThetaRecoilLab < CLHEP::pi/2.)
		return  acos((BetaRatio*GammaTangensSquared - sqrt(1.+GammaTangensSquared*(1.-pow(BetaRatio,2))))/(1+GammaTangensSquared));

	return acos((BetaRatio*GammaTangensSquared + sqrt(1.+GammaTangensSquared*(1.-pow(BetaRatio,2))))/(1+GammaTangensSquared));
}

//cosine of theta of recoil in cm-system, Beta = beta of cm-system, WORKS ONLY IF THE PARTICLE CAN REACH GIVEN LAB ANGLE (180 DEGREE)!!!
double Physics::CosineThetaRecoilCm(double Beta, double ECm, double ThetaRecoilLab, double EjectileMass, double RecoilMass)
{
	double BetaRatio = Beta/BetaRecoilCm(ECm,0.,EjectileMass,RecoilMass);
	double GammaTangensSquared = pow(tan(ThetaRecoilLab),2)/(1.-pow(Beta,2));

	//check whether the denominator will be zero or sqrt gets imaginary
	if(GammaTangensSquared == -1) {
		G4cerr<<__PRETTY_FUNCTION__<<": the denominator will be zero: "<<GammaTangensSquared<<G4endl;
		return -2.;
	}

	if((1.+GammaTangensSquared*(1.-pow(BetaRatio,2))) < 0) {
		//G4cerr<<__PRETTY_FUNCTION__<<": the sqrt will be imaginary: "<<(1.+GammaTangensSquared*(1.-pow(BetaRatio,2)))<<" (1+"<<GammaTangensSquared<<"*(1.-pow("<<BetaRatio<<",2)), Beta = "<<Beta<<")"<<G4endl;
		return -2.;
	}

	//this seems to be the right combination, but why?
	if(ThetaRecoilLab < CLHEP::pi/2.) {
		return (BetaRatio*GammaTangensSquared - sqrt(1.+GammaTangensSquared*(1.-pow(BetaRatio,2))))/(1+GammaTangensSquared);
	}

	return (BetaRatio*GammaTangensSquared + sqrt(1.+GammaTangensSquared*(1.-pow(BetaRatio,2))))/(1+GammaTangensSquared);
}

//kinetic energy of the recoil in the cm-system
double Physics::TRecoilCm(double ECm, double EjectileMass, double RecoilMass)
{
	double tfCm = TfCm(ECm,EjectileMass,RecoilMass);

	if(ECm == 0.) {
		G4cerr<<__PRETTY_FUNCTION__<<": ECm is 0"<<G4endl;
		return -1.;
	}

	return tfCm*(tfCm+2*(EjectileMass))/(2.*ECm);
}

//kinetic energy of recoil in the lab system
double Physics::TRecoilLab(double Beta, double ECm, double ThetaRecoilLab, double EjectileMass, double RecoilMass)
{
	double gamma = RelativisticGamma(Beta);
	double tRecoilCm = TRecoilCm(ECm,EjectileMass,RecoilMass);
	double cosThetaRecoilCm = CosineThetaRecoilCm(Beta, ECm, ThetaRecoilLab, EjectileMass, RecoilMass);

	//this happens only if the computation of ThetaRecoilCm failed
	if(cosThetaRecoilCm == -2.) {
		return 0.;
	}

	if(tRecoilCm < 0. || RecoilMass < 0. || cosThetaRecoilCm < -1. || cosThetaRecoilCm > 1.) {
		G4cerr<<__PRETTY_FUNCTION__<<": either the recoil mass ("<<RecoilMass<<") or the recoil kinetic energy in the cm system ("<<tRecoilCm<<") are less than 0 or the computation of theta of recoil in cm failed: "<<cosThetaRecoilCm<<G4endl;
		return -1;
	}

	return (gamma-1)*RecoilMass + gamma*tRecoilCm - gamma*Beta*sqrt(tRecoilCm*(tRecoilCm+2.*RecoilMass))*cosThetaRecoilCm;
}

double Physics::BeamEnergy(double ProjectileMass, double TargetMass, double LabEnergy)
{
	//everything (masses and energy) in the same unit as the result

	//from mathematica (fourth out of four solution, first two give imaginary results, third is always negative)
	//(-8*m1 - 4*m2 + t4Lab*pow(Sec(aLab),2) + Sec(aLab)*sqrt(16*pow(m2,2)*pow(cos(aLab),2) + t4Lab*(t4Lab*pow(Sec(aLab),2)*pow(cos(aLab) - sin(aLab),4) + 8*m2*(3 + sin(2*aLab)))) - 2*t4Lab*tan(aLab) + sqrt((pow(Sec(aLab),4)*(24*pow(m2,2)*t4Lab*pow(cos(aLab),2) + 32*pow(m1,2)*m2*pow(cos(aLab),4) + 2*m2*pow(t4Lab,2)*pow(cos(aLab) - sin(aLab),4) + 16*pow(m2,2)*t4Lab*pow(cos(aLab),3)*sin(aLab) - 8*pow(m1,2)*t4Lab*pow(cos(aLab),2)*(-1 + sin(2*aLab)) + 2*pow(m1,2)*cos(3*aLab)*sqrt(2*(4*pow(m2,2) + 12*m2*t4Lab + pow(t4Lab,2)) + (8*pow(m2,2) - 2*pow(t4Lab,2))*cos(2*aLab) + t4Lab*(t4Lab*pow(Sec(aLab),2) + 8*m2*sin(2*aLab) - 4*t4Lab*tan(aLab))) + 2*cos(aLab)*(3*pow(m1,2) + m2*t4Lab - m2*t4Lab*sin(2*aLab))*sqrt(2*(4*pow(m2,2) + 12*m2*t4Lab + pow(t4Lab,2)) + (8*pow(m2,2) - 2*pow(t4Lab,2))*cos(2*aLab) + t4Lab*(t4Lab*pow(Sec(aLab),2) + 8*m2*sin(2*aLab) - 4*t4Lab*tan(aLab)))))/m2))/8.

	double t2 = pow(TargetMass,2);
	double p2 = pow(ProjectileMass,2);
	double e2 = pow(LabEnergy,2);
	double te = TargetMass*LabEnergy;

	return (-8.*ProjectileMass - 4.*TargetMass + LabEnergy
			+ sqrt(16.*t2 + LabEnergy*(LabEnergy + 8.*TargetMass*3.))
			+ sqrt((24.*t2*LabEnergy + 
					32.*p2*TargetMass + 
					2.*TargetMass*e2 - 
					8.*p2*LabEnergy + 
					2.*p2*sqrt(2.*(4.*t2 + 12.*te + e2) + 
						(8.*t2 - 2.*e2) + 
						pow(LabEnergy,2)) + 
					2.*(3.*p2 + te)*sqrt(2.*(4.*t2 + 12.*te + e2) + (8.*t2 - 2.*e2) + pow(LabEnergy,2)))/TargetMass))/8.;
}

//no excitation and theta_lab of recoil = 0 degree; two solutions: +- Sqrt, +Sqrt yields a much higher beam energty (too high)
double Physics::BeamEnergy(double ProjectileMass, double TargetMass, double EjectileMass, double RecoilMass, double LabEnergy)
{
	double nominator = (-pow(TargetMass,3) + TargetMass*pow(EjectileMass,2) + 3.*pow(TargetMass,2)*RecoilMass - pow(EjectileMass,2)*RecoilMass - 3.*TargetMass*pow(RecoilMass,2) + pow(RecoilMass,3) + 
			3.*pow(TargetMass,2)*LabEnergy - pow(EjectileMass,2)*LabEnergy - 4.*TargetMass*RecoilMass*LabEnergy + pow(RecoilMass,2)*LabEnergy - 2.*TargetMass*pow(LabEnergy,2) + 
			pow(ProjectileMass,2)*(-TargetMass + RecoilMass + LabEnergy) - 2.*ProjectileMass*(pow(TargetMass,2) + pow(RecoilMass,2) - 2.*TargetMass*(RecoilMass + LabEnergy)) - 
			sqrt(LabEnergy*(2.*RecoilMass + LabEnergy)*(pow(ProjectileMass,4) + pow((pow(TargetMass,2) - pow(EjectileMass,2) + pow(RecoilMass,2) - 2.*TargetMass*(RecoilMass + LabEnergy)),2) - 
					2.*pow(ProjectileMass,2)*(pow(TargetMass,2) + pow(EjectileMass,2) + pow(RecoilMass,2) - 2.*TargetMass*(RecoilMass + LabEnergy)))));
	double denominator = 2.*(pow(TargetMass,2) + pow(RecoilMass,2) - 2.*TargetMass*(RecoilMass + LabEnergy));

	if(denominator == 0.) {
		return 0.;
	}

	return nominator/denominator;
}

//-------------------- geometric functions (solid angle and such) --------------------
double Physics::SinusTheta(double X, double Y, double Z, double Factor)
{
	return Factor/sqrt(1 + pow(Z,2)/(pow(X,2)+pow(Y,2)));
}

double Physics::SinusThetaZInt(double X, double Y, double Z1, double Z2, double Factor)
{
	return - Factor * sqrt(pow(X,2) + pow(Y,2)) * log((Z1 + sqrt(pow(X,2) + pow(Y,2) + pow(Z1,2)))/(Z2 + sqrt(pow(X,2) + pow(Y,2) + pow(Z2,2))));
}

double Physics::SolidAngle(double x, double y, double z)
{
	//solid angle at point x,y,z
	return y/pow(pow(x,2)+pow(y,2)+pow(z,2),3./2.);
}

double Physics::IntegratedSolidAngle(double x, double y, double z, double d, double w)
{
	//integration of solid angle over a plane at position y from x-d/2 to x+d/2 and z-w/2 to z+w/2
	//division by d*w is not necessary (this would normalize the result to be independent from d and w)
	double x2 = pow(x,2);
	double y2 = pow(y,2);
	double z2 = pow(z,2);
	double d2 = pow(d,2);
	double w2 = pow(w,2);

	double xl = pow(x - d/2.,2);
	double xh = pow(x + d/2.,2);
	double zl = pow(z - w/2.,2);
	double zh = pow(z + w/2.,2);

	return (2*(atan((4.*d*w*y)/(-(d2 + w2 - 4.*(x2 + y2 + z2))*sqrt(y2 + xl + zl) + (-w2 + pow(d - 2.*x,2) + 4.*(y2 + z2))*sqrt(y2 + xh + zl)
						+ (-d2 + 4.*(x2 + y2) + pow(w - 2.*z,2))*sqrt(y2 + xl + zh) + 4.*sqrt((y2 + xl + zl)*(y2 + xh + zl)*(y2 + xl + zh))))
				+ atan((4.*d*w*y)/((-d2 + 4.*(x2 + y2) + pow(w + 2.*z,2))*sqrt(y2 + xh + zl) + (-w2 + pow(d + 2.*x,2) + 4.*(y2 + z2))*sqrt(y2 + xl + zh)
						- (d2 + w2 - 4.*(x2 + y2 + z2))*sqrt(y2 + xh + zh) + 4.*sqrt((y2 + xh + zl)*(y2 + xl + zh)*(y2 + xh + zh))))));
}

double Physics::GammaEfficiency(double* x, double* par)
{
	if(x[0] > 0)
		return par[0]*log(x[0]) + par[1]*log(x[0])/x[0] + par[2]*pow(log(x[0]),2)/x[0] + par[3]*pow(log(x[0]),4)/x[0] + par[4]*pow(log(x[0]),5)/x[0];

	return 1.;
}

double Physics::GammaEfficiencyError(double* x, double* par)
{
	//parameter 0 is the number of parameters that were fitted
	//parameters 1 to par[0]^2 are the entries of the variance-covariance matrix
	double result = 0.;

	if(par[0] > 5 || par[0] < 0) {
		G4cerr<<"Error, par[0] in GammaEfficiencyError should be between 0 and 5!"<<G4endl;

		return 0.;
	}

	if(x[0] <= 0.) {
		return 0.;
	}

	double partialDerivative[5] = {log(x[0]), log(x[0])/x[0], pow(log(x[0]),2)/x[0], pow(log(x[0]),4)/x[0], pow(log(x[0]),5)/x[0]};

	for(int i = 0; i < par[0]; i++) {
		for(int j = 0; j < par[0]; j++) {
			result += partialDerivative[i]*partialDerivative[j]*par[i*((int)par[0])+j+1];
		}
	}

	if(result >= 0.) {
		return sqrt(result);
	}

	//cout<<G4endl<<"x[0] = "<<x[0]<<G4endl;
	//
	//for(int i = 0; i < par[0]; i++)
	//  {
	//    for(int j = 0; j < par[0]; j++)
	//	{
	//	  cout<<setw(16)<<partialDerivative[i]*partialDerivative[j]*par[i*((int)par[0])+j+1]<<" ";
	//	}
	//    cout<<G4endl;
	//  }
	//
	//for(int i = 0; i < par[0]; i++)
	//  {
	//    for(int j = 0; j < par[0]; j++)
	//	{
	//	  cout<<partialDerivative[i]<<"*"<<partialDerivative[j]<<"*"<<par[i*((int)par[0])+j+1]<<" ";
	//	}
	//    cout<<G4endl;
	//  }
	//
	//cout<<"===================="<<G4endl;

	return 0.;
}

double Physics::GammaEfficiencyLowerBound(double* x, double* par)
{
	//parameters 0-4 are passed to GammaEfficiency and the remaining parameters to GammaEfficiencyError
	return GammaEfficiency(x,par)-GammaEfficiencyError(x,par+5);
}

double Physics::GammaEfficiencyUpperBound(double* x, double* par)
{
	//parameters 0-4 are passed to GammaEfficiency and the remaining parameters to GammaEfficiencyError
	return GammaEfficiency(x,par)+GammaEfficiencyError(x,par+5);
}

//ratio between two gamma efficiencies (e.g. cluster to core :)
double Physics::GammaEfficiencyRatio(double* x, double* par)
{
	//0-4:denominator parameters, 5-9: divisor parameters
	if(x[0] > 0) {
		if((par[5]*x[0] + par[6] + par[7]*log(x[0]) + par[8]*pow(log(x[0]),3) + par[9]*pow(log(x[0]),4)) != 0) {
			return (par[0]*x[0] + par[1] + par[2]*log(x[0]) + par[3]*pow(log(x[0]),3) + par[4]*pow(log(x[0]),4))/
				(par[5]*x[0] + par[6] + par[7]*log(x[0]) + par[8]*pow(log(x[0]),3) + par[9]*pow(log(x[0]),4));
		}
	}

	return 1.;
}

double Physics::ReleaseCurve(double* x, double* par)
{
	//parameter 0: number of components
	//parameter 1+n*5: rise time of nth component
	//parameter 2+n*5: fast release time of nth component
	//parameter 3+n*5: slow release time of nth component
	//parameter 4+n*5: ratio fast/slow of nth component
	//parameter 5+n*5: amplitude of nth component

	return ReleaseCurve(x[0], par);
}

double Physics::ReleaseCurve(double x, double* par)
{
	//parameter 0: number of components
	//parameter 1+n*5: rise time of nth component
	//parameter 2+n*5: fast release time of nth component
	//parameter 3+n*5: slow release time of nth component
	//parameter 4+n*5: ratio fast/slow of nth component
	//parameter 5+n*5: amplitude of nth component

	double result = 0;
	for(int i = 0; i < (int) par[0]; i++) {
		if(((par[4+i*5]*pow(par[2+i*5],2)/(par[2+i*5]+par[1+i*5]) + (1.-par[4+i*5])*pow(par[3+i*5],2)/(par[3+i*5]+par[1+i*5]))/log(2.)) != 0) {
			result += par[5+i*5]*(1. - exp(-log(2.)*x/abs(par[1+i*5]))) * (par[4+i*5]*exp(-log(2.)*x/abs(par[2+i*5])) + 
					(1.-par[4+i*5])*exp(-log(2.)*x/abs(par[3+i*5])))
				/((par[4+i*5]*pow(par[2+i*5],2)/(par[2+i*5]+par[1+i*5]) + 
							(1.-par[4+i*5])*pow(par[3+i*5],2)/(par[3+i*5]+par[1+i*5]))/log(2.));
		} else {
			G4cerr<<"Error, normalization is zero: "<<((par[4+i*5]*pow(par[2+i*5],2)/(par[2+i*5]+par[1+i*5]) + (1.-par[4+i*5])*pow(par[3+i*5],2)/(par[3+i*5]+par[1+i*5]))/log(2.))<<" = (("<<par[4+i*5]*pow(par[2+i*5],2)/(par[2+i*5]+par[1+i*5])<<" + "<<(1.-par[4+i*5])*pow(par[3+i*5],2)<<"/"<<(par[3+i*5]+par[1+i*5])<<")/log(2.))"<<G4endl;
			return 0.;
		}
	}

	return result;
}

double Physics::IsoldeReleaseCurve(double* x, double* par)
{
	//parameters 0-2: probability to get each 1-3th proton pulse
	//parameter 3: number of components
	//plus 5 parameters for each component => 5*(# of components)+4 parameters
	double result = ReleaseCurve(x[0],par+3) + 
		par[0]*ReleaseCurve(x[0]+1200.,par+3) + 
		par[1]*ReleaseCurve(x[0]+2400.,par+3) + 
		par[2]*ReleaseCurve(x[0]+3600.,par+3) +
		ReleaseCurve(x[0]+4800.,par+3);

	if(x[0] < 1200.) {
		return result;
	}

	for(int i = 0; i < 3; i++) {
		if(x[0] < 2400.+i*1200.) {
			return (1.-par[i])*result;
		}
	}

	return 0.;
}
