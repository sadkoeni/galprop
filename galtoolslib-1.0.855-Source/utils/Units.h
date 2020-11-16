#ifndef _utl_SystemOfUnits_h_
#define _utl_SystemOfUnits_h_

// Partition the global namespace in order to avoid ambiguity between 
// Geant4 units

namespace utl {

  /*
    The conversion factors defined in this file 
    convert your data into a set of base units, so that
    all dimensional quantities in the code are in a
    single system of units.  You can also 
    use the conversions defined here to, for example,
    display data with the unit of your choice.  For example:
    \code
    cout << "s = " << s/mm << " mm";
    \endcode
    The base units are : 
    - meter                   (meter)
    - nanosecond              (nanosecond)
    - electron Volt           (eV)
    - positron charge         (eplus)
    - degree Kelvin           (kelvin)
    - the amount of substance (mole)
    - luminous intensity      (candela)
    - radian                  (radian)
    - steradian               (steradian)
    
    Below is a non-exhaustive list of derived and pratical units
    (i.e. mostly the SI units).
    
    The SI numerical value of the positron charge is defined here,
    as it is needed for conversion factor : positron charge = e_SI (coulomb)
    
    This is a slightly modified version of the units definitions 
    written by the Geant4 collaboration
  
  */

  // 
  // Length [L]
  //
  constexpr double meter  = 1.0;
  constexpr double meter2 = meter*meter;
  constexpr double meter3 = meter*meter*meter;
  
  constexpr double millimeter  = 1.e-3*meter;                        
  constexpr double millimeter2 = millimeter*millimeter;
  constexpr double millimeter3 = millimeter*millimeter*millimeter;
  
  constexpr double centimeter  = 10.*millimeter;   
  constexpr double centimeter2 = centimeter*centimeter;
  constexpr double centimeter3 = centimeter*centimeter*centimeter;
  
  constexpr double kilometer = 1000.*meter;                   
  constexpr double kilometer2 = kilometer*kilometer;
  constexpr double kilometer3 = kilometer*kilometer*kilometer;
  
  constexpr double micrometer = 1.e-6*meter;             
  constexpr double micron     = 1.e-6*meter;
  constexpr double nanometer  = 1.e-9*meter;
  constexpr double angstrom   = 1.e-10*meter;
  constexpr double fermi      = 1.e-15*meter;
  
  constexpr double      barn = 1.e-28*meter2;
  constexpr double millibarn = 1.e-3 *barn;
  constexpr double microbarn = 1.e-6 *barn;
  constexpr double  nanobarn = 1.e-9 *barn;
  constexpr double  picobarn = 1.e-12*barn;
  
  // symbols
  constexpr double mm  = millimeter;                        
  constexpr double mm2 = millimeter2;
  constexpr double mm3 = millimeter3;
  
  constexpr double cm  = centimeter;   
  constexpr double cm2 = centimeter2;
  constexpr double cm3 = centimeter3;
  
  constexpr double m  = meter;                  
  constexpr double m2 = meter2;
  constexpr double m3 = meter3;
  
  constexpr double km  = kilometer;                   
  constexpr double km2 = kilometer2;
  constexpr double km3 = kilometer3;
  
  //
  // Angle
  //
  constexpr double radian      = 1.;                  
  constexpr double milliradian = 1.e-3*radian;
  constexpr double degree = (3.14159265358979323846/180.)*radian;
  
  constexpr double   steradian = 1.;
  
  // symbols
  constexpr double rad  = radian;	
  constexpr double mrad = milliradian;
  constexpr double sr   = steradian;
  constexpr double deg  = degree;
  
  //
  // Time [T]
  //
  constexpr double nanosecond  = 1.;
  constexpr double second      = 1.e+9 *nanosecond;
  constexpr double millisecond = 1.e-3 *second;
  constexpr double microsecond = 1.e-6 *second;
  constexpr double  picosecond = 1.e-12*second;
  constexpr double minute      = 60*second;
  constexpr double hour        = 60*minute;
  constexpr double day         = 24*hour;

  constexpr double hertz = 1./second;
  constexpr double kilohertz = 1.e+3*hertz;
  constexpr double megahertz = 1.e+6*hertz;
  
  // symbols
  constexpr double ns = nanosecond;			
  constexpr double  s = second;
  constexpr double ms = millisecond;

  //
  // Electric charge [Q]
  //
  constexpr double eplus = 1. ;		// positron charge
  constexpr double e_SI  = 1.602176462e-19;	// positron charge in coulomb
  constexpr double coulomb = eplus/e_SI;	// coulomb = 6.24150 e+18*eplus
  
  //
  // Energy [E]
  //
  constexpr double     electronvolt = 1.;
  constexpr double megaelectronvolt = 1.e+6*electronvolt;
  constexpr double kiloelectronvolt = 1.e+3*electronvolt;
  constexpr double gigaelectronvolt = 1.e+9*electronvolt;
  constexpr double teraelectronvolt = 1.e+12*electronvolt;
  constexpr double petaelectronvolt = 1.e+15*electronvolt;
  constexpr double exaelectronvolt  = 1.e+18*electronvolt;
  constexpr double zettaelectronvolt= 1.e+21*electronvolt;
  
  constexpr double joule = electronvolt/e_SI; // joule = 6.24150 e+12 * MeV
  constexpr double erg   = 1.0e-7*joule;
  
  // symbols
  constexpr double MeV = megaelectronvolt;
  constexpr double  eV = electronvolt;
  constexpr double keV = kiloelectronvolt;
  constexpr double GeV = gigaelectronvolt;
  constexpr double TeV = teraelectronvolt;
  constexpr double PeV = petaelectronvolt;
  constexpr double EeV = exaelectronvolt;
  constexpr double ZeV = zettaelectronvolt;
  
  //
  // Mass [E][T^2][L^-2]
  //
  constexpr double  kilogram = joule*second*second/(meter*meter);   
  constexpr double      gram = 1.e-3*kilogram;
  constexpr double milligram = 1.e-3*gram;
  
  // symbols
  constexpr double  kg = kilogram;
  constexpr double   g = gram;
  constexpr double  mg = milligram;
  
  //
  // Power [E][T^-1]
  //
  constexpr double watt = joule/second; // watt = 6.24150 e+3 * MeV/ns
  
  //
  // Force [E][L^-1]
  //
  constexpr double newton = joule/meter; // newton = 6.24150 e+9 * MeV/mm
  
  //
  // Pressure [E][L^-3]
  //
  constexpr double hep_pascal = newton/m2; // pascal = 6.24150 e+3 * MeV/mm3
  constexpr double bar = 100000*hep_pascal; // bar = 6.24150 e+8 * MeV/mm3
  constexpr double atmosphere = 101325*hep_pascal; // atm = 6.32420e+8*MeV/mm3
  
  //
  // Electric current [Q][T^-1]
  //
  constexpr double ampere = coulomb/second; // ampere = 6.24150e+9*eplus/ns
  constexpr double milliampere = 1.e-3*ampere;
  constexpr double microampere = 1.e-6*ampere;
  constexpr double  nanoampere = 1.e-9*ampere;
  
  //
  // Electric potential [E][Q^-1]
  //
  constexpr double megavolt = megaelectronvolt/eplus;
  constexpr double kilovolt = 1.e-3*megavolt;
  constexpr double     volt = 1.e-6*megavolt;
  
  //
  // Electric resistance [E][T][Q^-2]
  //
  constexpr double ohm = volt/ampere; // ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)
  
  //
  // Electric capacitance [Q^2][E^-1]
  //
  constexpr double farad = coulomb/volt; // farad = 6.24150e+24 * eplus/Megavolt
  constexpr double millifarad = 1.e-3*farad;
  constexpr double microfarad = 1.e-6*farad;
  constexpr double  nanofarad = 1.e-9*farad;
  constexpr double  picofarad = 1.e-12*farad;
  
  //
  // Magnetic Flux [T][E][Q^-1]
  //
  constexpr double weber = volt*second; // weber = 1000*megavolt*ns

  //
  // Magnetic Field [T][E][Q^-1][L^-2]
  //
  constexpr double tesla = volt*second/meter2; // tesla =0.001*megavolt*ns/mm2
  
  constexpr double gauss     = 1.e-4*tesla;
  constexpr double kilogauss = 1.e-1*tesla;
  
  //
  // Inductance [T^2][E][Q^-2]
  //
  constexpr double henry = weber/ampere; // henry = 1.60217e-7*MeV*(ns/eplus)**2
  
  //
  // Temperature
  //
  constexpr double kelvin = 1.;
  
  //
  // Amount of substance
  //
  constexpr double mole = 1.;
  
  //
  // Activity [T^-1]
  //
  constexpr double becquerel = 1./second ;
  constexpr double curie = 3.7e+10 * becquerel;
  
  //
  // Absorbed dose [L^2][T^-2]
  //
  constexpr double gray = joule/kilogram ;
  
  //
  // Luminous intensity [I]
  //
  constexpr double candela = 1.;
  
  //
  // Luminous flux [I]
  //
  constexpr double lumen = candela*steradian;
  
  //
  // Illuminance [I][L^-2]
  //
  constexpr double lux = lumen/meter2;
  
  //
  // Miscellaneous
  //
  constexpr double perCent     = 0.01 ;
  constexpr double perThousand = 0.001;
  constexpr double perMillion  = 0.000001;
  
} // end of utl namespace

#endif 
