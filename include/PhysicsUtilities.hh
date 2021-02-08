#ifndef PHYSICS_UTILITIES_HH 
#define PHYSICS_UTILITIES_HH 

#include <stdio.h>
#include <iostream>
#include <iomanip>

namespace Physics 
{
  double bwe(int, int, int);
  
  double SolidAngleConversionRecoil(double, double, double);
  double SolidAngleConversionEjectile(double, double, double);

  double SolidAngleConversionRecoil(double beamEnergy, double projectileMass, double targetMass, double ejectileMass, double recoilMass, double eEx, double recoilThetaCm);
  double SolidAngleConversionEjectile(double beamEnergy, double projectileMass, double targetMass, double ejectileMass, double recoilMass, double eEx, double ejectileThetaCm);

  double SolidAngleConversionRecoil(double*, double*);
  double SolidAngleConversionEjectile(double*, double*);

  double ThetaLab(double, double, double);
  double ThetaLab(double beamEnergy, double projectileMass, double targetMass, double ejectileMass, double recoilMass, double eEx, double recoilThetaCm);
  
  double ThetaCm(double, double, double);
  
  double RelativisticGamma(double Beta);

  double TiCm(double BeamEnergy, double ProjectileMass, double TargetMass);

  double TfCm(double ECm, double Eex, double m3, double m4);
  double TfCm(double ECm, double m3, double m4);
  
  double TProjectileCm(double BeamEnergy, double ProjectileMass, double TargetMass);
  double TRecoilCm(double ECm, double Eex, double m3, double m4);
  double TRecoilCm(double ECm, double m3, double m4);
  
  double BetaRecoilCm(double ECm, double Eex, double m3, double m4);
  double BetaRecoilCm(double Beta, double ThetaRecoilLab, double ThetaRecoilCm);
  
  double ThetaRecoilCm(double Beta, double ECm, double Eex, double ThetaRecoilLab, double m3, double m4);
  double ThetaRecoilCm(double Beta, double ECm, double ThetaRecoilLab, double m3, double m4);
  
  double CosineThetaRecoilCm(double Beta, double ECm, double Eex, double ThetaRecoilLab, double m3, double m4);
  double CosineThetaRecoilCm(double Beta, double ECm, double ThetaRecoilLab, double m3, double m4);
    
  double Eex(double ECm, double Beta, double m3, double m4, double ThetaRecoilLab, double ThetaRecoilCm);
  
  double ThetaRecoilLab(double Beta, double ECm, double Eex, double ThetaRecoilCm, double EjectileMass, double RecoilMass);

  double TRecoilLab(double Eex, double ThetaRecoilLab, double BeamEnergy, double ProjectileMass, double TargetMass, double EjectileMass, double RecoilMass);
  double TRecoilLab(double Beta, double ECm, double Eex, double ThetaRecoilLab, double EjectileMass, double RecoilMass);
  double TRecoilLab2(double Beta, double ECm, double Eex, double ThetaRecoilCm, double EjectileMass, double RecoilMass);
  double TRecoilLab(double Beta, double ECm, double ThetaRecoilLab, double EjectileMass, double RecoilMass);
  
  double TEjectileCm(double ECm, double Eex, double m3, double m4);
  double BetaEjectileCm(double ECm, double Eex, double m3, double m4);
  double CosineThetaEjectileCm(double Beta, double ECm, double Eex, double ThetaEjectileLab, double m3, double m4);
  double CosineThetaEjectileCm(double beamEnergy, double projectileMass, double targetMass, double ejectileMass, double recoilMass, double eEx, double thetaEjectileLab);
  double TEjectileLab(double Beta, double ECm, double Eex, double ThetaEjectileLab, double EjectileMass, double RecoilMass);

  double ECm(double BeamEnergy, double ProjectileMass, double TargetMass);
  
  double BetaCm(double BeamEnergy, double ProjectileMass, double TargetMass);

  double BeamEnergy(double ProjectileMass, double TargetMass, double LabAngle, double LabEnergy);
  double BeamEnergy(double ProjectileMass, double TargetMass, double LabEnergy);
  double BeamEnergy(double ProjectileMass, double TargetMass, double EjectileMass, double RecoilMass, double LabEnergy);

  double SinusTheta(double X, double Y, double Z, double Factor = 1.);

  double SinusThetaZInt(double X, double Y, double Z1, double Z2, double Factor = 1.);

  double SolidAngle(double x, double y, double z);

  double IntegratedSolidAngle(double x, double y, double z, double d, double w);

  double GammaEfficiency(double* x, double* par);

  double GammaEfficiencyError(double* x, double* par);

  double GammaEfficiencyLowerBound(double* x, double* par);

  double GammaEfficiencyUpperBound(double* x, double* par);

  double GammaEfficiencyRatio(double* x, double* par);

  double ReleaseCurve(double x, double* par);

  double ReleaseCurve(double* x, double* par);

  double IsoldeReleaseCurve(double* x, double* par);
}

#endif
