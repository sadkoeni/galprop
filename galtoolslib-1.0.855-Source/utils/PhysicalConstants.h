#ifndef _utl_PhysicalConstants_h_
#define _utl_PhysicalConstants_h_

#include <Units.h>

#include <cmath>

namespace utl {

  constexpr double kPi = 3.14159265358979323846264338327950288419716939937510582097494;
  const double kExp = std::exp(1.); 
  const double kSqrtTwo = std::sqrt(2.);
  const double kLogTwo = std::log(2.);
  const double kSqrtThree = std::sqrt(3.);
  constexpr double kOneOnThree = 1./3.;
  constexpr double kOneOnTwo = 1./2.;
  constexpr double kFourPiOnThree = 4.*kPi/3.;
  constexpr double kPiOnTwo = kPi/2.;
  constexpr double kTwoPi = 2.*kPi;
  constexpr double kFourPi = 2.*kTwoPi;
  constexpr double kOneOnPi = 1./kPi;
  constexpr double kOneOnTwoPi = 1./kTwoPi;
  constexpr double kOneOnFourPi = 1./kFourPi;
  constexpr double kConvertDegreesToRadians = kPi/180.;
  constexpr double kConvertRadiansToDegrees = 1./kConvertDegreesToRadians;

  // All taken from PDG data tables (2002)

  // Physical constants 

  constexpr double kSpeedOfLight_SI = 299792458.0;
  constexpr double kSpeedOfLight = kSpeedOfLight_SI*m/s;
  constexpr double kPlanck_SI = 6.62606876e-34;
  constexpr double kPlanckReduced_SI = kPlanck_SI*kOneOnTwoPi;
  constexpr double kPlanck = kPlanck_SI*joule*s;
  constexpr double kPlanckReduced = kPlanckReduced_SI*joule*s;
  constexpr double kMuZero_SI = kFourPi*1.e-7;
  constexpr double kMuZero = kMuZero_SI*newton/(ampere*ampere); 
  constexpr double kBoltzmann_SI = 1.3806503e-23;
  constexpr double kBoltzmann = kBoltzmann_SI*joule/kelvin;

  constexpr double kPlanckTimesSpeedOfLight_SI = kPlanck_SI*kSpeedOfLight_SI;
  constexpr double kPlanckTimesSpeedOfLightSquared_SI = 
    kPlanckTimesSpeedOfLight_SI*kPlanckTimesSpeedOfLight_SI;
  
  // Particle and other masses

  constexpr double kMassConversion_SI = e_SI/(kSpeedOfLight_SI*kSpeedOfLight_SI);

  constexpr double kHydrogenMass_SI = 1.6735e-27;
  constexpr double kHydrogenMass = kHydrogenMass_SI*kg;

  constexpr double kElectronMass = 0.510998902*MeV; 
  constexpr double kElectronMass_SI = kElectronMass*kMassConversion_SI;
  constexpr double kMuonMass = 105.658357*MeV; 
  constexpr double kMuonMass_SI = kMuonMass*kMassConversion_SI;
  constexpr double kTauMass = 1776.99*MeV;
  constexpr double kTauMass_SI = kTauMass*kMassConversion_SI;

  constexpr double kProtonMass = 938.271998*MeV; 
  constexpr double kProtonMass_SI = kProtonMass*kMassConversion_SI;
  constexpr double kNeutronMass = 939.56533*MeV; 
  constexpr double kNeutronMass_SI = kNeutronMass*kMassConversion_SI;
  constexpr double kDeuteronMass = 1875.612762*MeV; 
  constexpr double kDeuteronMass_SI = kDeuteronMass*kMassConversion_SI;

  constexpr double kLambdaMass = 1115.683*MeV;
  constexpr double kLambdaMass_SI = kLambdaMass*kMassConversion_SI;
  constexpr double kSigmaZeroMass = 1192.642*MeV;
  constexpr double kSigmaZeroMass_SI = kSigmaZeroMass*kMassConversion_SI;
  constexpr double kSigmaPlusMass = 1189.37*MeV;
  constexpr double kSigmaPlusMass_SI = kSigmaPlusMass*kMassConversion_SI;
  constexpr double kSigmaMinusMass = 1197.449*MeV;
  constexpr double kSigmaMinusMass_SI = kSigmaMinusMass*kMassConversion_SI;
  constexpr double kXiZeroMass = 1314.83*MeV;
  constexpr double kXiZeroMass_SI = kXiZeroMass*kMassConversion_SI;
  constexpr double kXiMinusMass = 1321.31*MeV;
  constexpr double kXiMinusMass_SI = kXiMinusMass*kMassConversion_SI;
  constexpr double kOmegaMinusMass = 1672.45*MeV;
  constexpr double kOmegaMinusMass_SI = kOmegaMinusMass*kMassConversion_SI;

  constexpr double kPiZeroMass = 134.9766*MeV; 
  constexpr double kPiZeroMass_SI = kPiZeroMass*kMassConversion_SI;
  constexpr double kPiChargedMass = 139.57018*MeV; 
  constexpr double kPiChargedMass_SI = kPiChargedMass*kMassConversion_SI;
  constexpr double kKaonChargedMass = 493.677*MeV; 
  constexpr double kKaonChargedMass_SI = kKaonChargedMass*kMassConversion_SI;

  constexpr double kAtomicMassUnit_SI = 1.660538e-27;

  constexpr double kCarbonMass_SI = 12.0107*kAtomicMassUnit_SI;
  constexpr double kCarbonMass = kCarbonMass_SI*kg/g;
  constexpr double kOxygenMass_SI = 15.9994*kAtomicMassUnit_SI;
  constexpr double kOxygenMass = kOxygenMass_SI*kg/g;
  constexpr double kMagnesiumMass_SI = 24.3050*kAtomicMassUnit_SI;
  constexpr double kMagnesiumMass = kMagnesiumMass_SI*kg/g;
  constexpr double kSiliconMass_SI = 28.0855*kAtomicMassUnit_SI;
  constexpr double kSiliconMass = kSiliconMass_SI*kg/g;
  constexpr double kIronMass_SI = 55.845*kAtomicMassUnit_SI;
  constexpr double kIronMass = kIronMass_SI*kg/g;

  constexpr double kSilicateMass_SI = kMagnesiumMass_SI + kIronMass_SI + 
    kSiliconMass_SI + 4.*kOxygenMass_SI; // Silicate is MgFeSi0_4
  constexpr double kSilicateMass = kSilicateMass_SI*kg/g;

  constexpr double kGraphiteDensity = 2.24; // g cm^-3
  constexpr double kSilicateDensity = 3.50; // g cm^-3

  constexpr double kGraphiteDensity_SI = kGraphiteDensity*(g/kg)/(cm/m*cm/m*cm/m);//pow(cm/m, 3.); 
  constexpr double kSilicateDensity_SI = kSilicateDensity*(g/kg)/(cm/m*cm/m*cm/m);//pow(cm/m, 3.); 
  
  // Particle lifetimes

  constexpr double kMuonLifetime = 2.19703e-6*s;

  constexpr double kNeutronLifetime = 885.7*s;

  constexpr double kLambdaLifetime = 2.632e-10*s;
  constexpr double kSigmaZeroLifetime = 7.4e-20*s;
  constexpr double kSigmaPlusLifetime = 0.8018e-10*s;
  constexpr double kSigmaMinusLifetime = 1.479e-10*s;
  constexpr double kXiZeroLifetime = 2.9e-10*s;
  constexpr double kXiMinusLifetime = 1.639e-10*s;
  constexpr double kOmegaMinusLifetime = 0.821-10*s;

  constexpr double kPiZeroLifetime = 8.4e-17*s;
  constexpr double kPiChargedLifetime = 2.6033e-8*s;
  constexpr double kKaonChargedLifetime = 1.2384e-8*s;

  // Derived constants

  constexpr double kEpsilonZero_SI = 1.0/(kMuZero_SI*kSpeedOfLight_SI*kSpeedOfLight_SI);
  constexpr double kAlpha = (e_SI*e_SI)/
    (kFourPi*kEpsilonZero_SI*kPlanckReduced_SI*kSpeedOfLight_SI); 
  constexpr double kElectronRadius_SI = (e_SI*e_SI)/
    (kFourPi*kEpsilonZero_SI*kElectronMass_SI*
     kSpeedOfLight_SI*kSpeedOfLight_SI);
  constexpr double kThomsonCrossSection_SI = 
    8.0*kPi*kElectronRadius_SI*kElectronRadius_SI/3.0;

  // Distance conversions

  constexpr double kParsec = 3.0856775807e+16*m; 
  constexpr double kKiloParsec = kParsec*1.0e+3;
  constexpr double kMegaParsec = kKiloParsec*1.0e+3;

  constexpr double pc = kParsec;
  constexpr double kpc = kKiloParsec;
  constexpr double Mpc = kMegaParsec;

  // Some other conversions and constants

  constexpr double kYearToSec = 365.25*24.0*60.0*60.0;
  constexpr double kSolarLuminosity_SI = 3.846e26;
  constexpr double kSolarLuminosity = kSolarLuminosity_SI*watt;

  // Temperature, absolute magnitude at V, bolometric correction, magnitude, 
  // and mass for Sun 

  constexpr double kTemperatureSun = 5770;
  constexpr double kMVSun = 4.82;
  constexpr double kBCSun = -0.05;
  constexpr double kMBolometricSun = kMVSun + kBCSun;
  constexpr double kSurfaceGravitySun = 2.74e4; // cgs (cm s^-2)
  constexpr double kSolarMass_SI = 1.98844e30; 
  //constexpr double kSolarMass = kSolarMass_SI*kilogram/gram; // grams

  // Conversion from LSun/kpc^3 to eV cm^-3 sr^-1 kpc^-1
  // -> LSun (J s^-1) / (kpc/cm)^3 * kpc/cm / (4*Pi*c)
  constexpr double kWainscoatConversion = 
    (kSolarLuminosity_SI/e_SI)*(cm/kpc*cm/kpc*cm/kpc)*//pow(cm/kpc, 3.0)*
    (kpc/kSpeedOfLight_SI)*kOneOnFourPi;
    
  // Conversion from LSun/pc^3 to eV cm^-3 sr^-1 kpc^-1
  // -> LSun (J s^-1) / (pc/cm)^3 * kpc/cm / (4*Pi*c)
  constexpr double kMathisConversion = 
    (kSolarLuminosity_SI/e_SI)/(pc/cm*pc/cm*pc/cm)*//pow(pc/cm, 3.0)*
    (kpc/kSpeedOfLight_SI)*kOneOnFourPi;
    
  // Conversion from MJy sr^-1 kpc^-1 -> eV cm^-3 sr^-1 kpc^-1.
  // Take M == 10^6, Jy == 10^-26 W m^-2 Hz^-1 -> 10^-26/e_SI eV m^-2 Hz^-1 
  // 1/c converts cm^-2 s^-1 to cm^-3 if in kpc s^-1
  constexpr double kFreudenreichConversion = 
    1.e6*(1.e-26/e_SI)*(cm/kSpeedOfLight_SI)*(cm*cm);
  
  // Conversion of W sr^-1 H-atom^-1 H-atom cm^-3 -> eV cm^-3 sr^-1 kpc^-1
  // W -> J s^-1/e_SI -> eV s^-1
  // 1/c converts cm^-2 s^-1 -> cm^-2 kpc^-1 if c in kpc s^-1
  // So get W sr^-1 cm^-3 -> eV cm^-3 sr^-1 kpc^-1 
  // times factor 10^-32 associated with emissivity from table 

  constexpr double kSodroskiConversion = (1.e-32/e_SI)*(kpc/kSpeedOfLight_SI);

  constexpr double kFrequencyConstant = 
    2.*kPlanck_SI/(kSpeedOfLight_SI*kSpeedOfLight_SI);//pow(kSpeedOfLight_SI, 2.0);
  constexpr double kFrequencyExponentConstant = kPlanck_SI/kBoltzmann_SI;
  constexpr double kWavelengthConstant = 
    2.*kPlanck_SI*(kSpeedOfLight_SI*kSpeedOfLight_SI);//pow(kSpeedOfLight_SI, 2.0);
  constexpr double kWavelengthExponentConstant = 
    kPlanck_SI*kSpeedOfLight_SI/kBoltzmann_SI;
  constexpr double kPlanckIntegral = 2.*(kBoltzmann_SI*kPi*kBoltzmann_SI*kPi*kBoltzmann_SI*kPi*kBoltzmann_SI*kPi)//pow(kBoltzmann_SI*kPi, 4.0)
    /(15.*(kPlanck_SI*kSpeedOfLight_SI*kPlanck_SI*kSpeedOfLight_SI*kPlanck_SI));//pow(kPlanck_SI*kSpeedOfLight_SI, 2.0)*kPlanck_SI);
  constexpr double kStefanBoltzmann = kPlanckIntegral*kPi;

}				  

#endif 

