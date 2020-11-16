
#include "spiralarms.h"

#include <vector>
#include <map>
#include <limits>
#include <Registry.h>

#define POW2(x) (x)*(x)

namespace GalacticStructure {

//Register the scale functions
static utl::Registry0<ScaleFunction>::Registrar<ExpFunction> registrarExpFunction("Exponential");
static utl::Registry0<ScaleFunction>::Registrar<GaussianFunction> registrarGaussFunction("Gaussian");
static utl::Registry0<ScaleFunction>::Registrar<Sech2Function> registrarSech2Function("Sech2");
static utl::Registry0<ScaleFunction>::Registrar<SechFunction> registrarSechFunction("Sech");
static utl::Registry0<ScaleFunction>::Registrar<Sech1_2Function> registrarSech1_2Function("Sech1_2");
static utl::Registry0<ScaleFunction>::Registrar<StepFunction> registrarStepFunction("Step");

#ifdef HAVE_OPENCL
std::string ScaleFunction::createFunction(const std::string &name, const std::string &returnCommand) const
{
   std::ostringstream codeStr;
   codeStr<<"inline float "<<name<<"( float x )\n"
      <<"{\n"
      <<"return "<<returnCommand<<";\n"
      <<"}\n";
   return codeStr.str();
}

std::string ExpFunction::getOpenCLFunction(std::string &name) const
{
   name = "ExpFunction_"+name;
   return createFunction(name, "exp( -x )");
}

std::string GaussianFunction::getOpenCLFunction(std::string &name) const
{
   name = "GaussianFunction_"+name;
   return createFunction(name, "exp( -x*x )");
}

std::string Sech2Function::getOpenCLFunction(std::string &name) const
{
   name = "Sech2Function_"+name;
   return createFunction(name, "4./(exp(2*x) + 2 + exp(-2*x))");
}

std::string SechFunction::getOpenCLFunction(std::string &name) const
{
   name = "SechFunction_"+name;
   return createFunction(name, "2./(exp(x) + exp(-x))");
}

std::string Sech1_2Function::getOpenCLFunction(std::string &name) const
{
   name = "Sech1_2Function_"+name;
   return createFunction(name, "sqrt((float)2.0)/sqrt(exp(x) + exp(-x))");
}

std::string StepFunction::getOpenCLFunction(std::string &name) const
{
   name = "StepFunction_"+name;
   return createFunction(name, "x <= 1 ? 1 : 0;");
}

#endif

ArmFunction::ArmFunction ( double a, double rMin, double phiMin, double phiExtent, double width,  std::unique_ptr<ScaleFunction> scale ) :
fA(a),
fRMin(rMin),
fPhiMin(phiMin),
fPhiExtent(phiExtent),
fWidth(width),
fscale(std::move(scale))
{
  armDistPointer = &GalacticStructure::ArmFunction::armDistLinear;
  fOneOverWidth = 1./fWidth;
  fCosAlpha = cos(atan(1./a));
  calculateRPdata();
  calculateXYdata();
}

double ArmFunction::operator() ( double radius, double phi ) const {
   const double dist = (this->*armDistPointer)( radius, phi );
   return (*fscale)( dist * fOneOverWidth );
}

void ArmFunction::calculateRPdata()
{
  //Use .5 degree steps, i.e. 720 steps per circle, and at least 2
  const size_t nsteps = fabs(fPhiExtent)*360/M_PI+2;
  fphiStep = fPhiExtent/(nsteps-1);
  
  nStepsPerRing = 2*M_PI/fabs(fphiStep);
  
  fRPArmData.resize(nsteps);
  
  for ( size_t i(0); i < nsteps; ++i ) {
    fRPArmData[i].second = fPhiMin + fphiStep*i;
    fRPArmData[i].first = fRMin * exp( fphiStep*i/fA );
  }
}

void ArmFunction::calculateXYdata()
{
  fXYArmData.resize(fRPArmData.size());

  for ( size_t i(0); i < fXYArmData.size(); ++i ) {
    fXYArmData[i].first = fRPArmData[i].first*cos(fRPArmData[i].second);
    fXYArmData[i].second = fRPArmData[i].first*sin(fRPArmData[i].second);
  }
}

void ArmFunction::CompareArmDist( double r, double phi, double &acc, double &app ) const
{
   acc = armDistBF(r, phi);
   app = armDistLinear(r, phi);
}

//Radial distance to arm
double ArmFunction::armDistRadial ( double rad, double phi ) const
{
   const double armPhi = angle(rad);

   if (armPhi < fPhiMin || armPhi > fPhiExtent) 
      return armDistLinear ( rad, phi );
   
   double angDist = phi - armPhi;

   //Make sure the distance is within -pi and pi
   while (angDist < -M_PI)
      angDist += 2*M_PI;
   while (angDist > M_PI)
      angDist -= 2*M_PI;

   return rad*fabs(angDist);
}


//Should work perfectly for infinite arms
//Also only works for pure logarithmic spiral arms.
double ArmFunction::armDistLinear ( double rad, double phi ) const
{

   //Find the closest distance, test phi values with 2 pi offset
   double distance(std::numeric_limits<double>::max());

   double mphi(phi);

   while ( mphi <= fPhiExtent + fPhiMin ) {

      if (mphi >= fPhiMin) {

         const double newDist = fabs(rad - radius(mphi));

         if (newDist < distance)
            distance = newDist;
      }

      mphi += 2*M_PI;

   }

   mphi = phi - M_PI;

   while ( mphi <= fPhiExtent + fPhiMin ) {

      if (mphi >= fPhiMin) {

         const double newDist = rad + radius(mphi);

         if (newDist < distance)
            distance = newDist;
      }

      mphi += 2*M_PI;

   }


   return distance * fCosAlpha;

}


double ArmFunction::armDist ( double rad, double phi ) const
{

  const double x = rad*cos(phi);
  const double y = rad*sin(phi);

  //Find the smallest distance to the spiral arm.  Try to be smart and use
  //the points that lie on the same angle as starting points.  Due the 
  //nature of spiral arms, look towards negative theta for larger R and 
  //positive theta for smaller R

  //Use points on the same angle as starting points and find the 
  //largest radius smaller than current and smallest larger than current
  //Always test end points
  int i1(0);
  size_t iphi(0);

  size_t ismall(0), ilarge(fRPArmData.size()-1);
  double rsmall(fRPArmData[ismall].first), rlarge(fRPArmData[ilarge].first);

  while (true) {

     i1 = int(floor( ( (phi + iphi*M_PI) - fPhiMin ) / fphiStep));
     iphi += 2;

     if (i1 < 0)
        continue;

     if (size_t(i1) > fRPArmData.size())
        break;

     if ( fRPArmData[i1].first < rad ) {
        ismall = i1;
        rsmall = fRPArmData[i1].first;
     } else {
        ilarge = i1;
        rlarge = fRPArmData[i1].first;
        break;
     }

  }


  //Use function to find nearest point closest to the two values, rmin and rlarge
  double sqrsmin = std::numeric_limits<double>::max();
  size_t imin(0);

  findBestPoint(imin, sqrsmin, ismall, rsmall, rad, x, y);
  findBestPoint(imin, sqrsmin, ilarge, rlarge, rad, x, y);


  // Find the nearest neighbour with smaller values and use linear interpolation between the points
  const double sqrsl = (imin == 0) ? std::numeric_limits< double >::max() : POW2(x-fXYArmData[imin-1].first) + POW2(y-fXYArmData[imin-1].second);
  const double sqrsu = (imin == fXYArmData.size()-1) ? std::numeric_limits< double >::max() : POW2(x-fXYArmData[imin+1].first) + POW2(y-fXYArmData[imin+1].second);
  
  const size_t il = (sqrsl < sqrsu) ? imin-1 : imin;
  
  // Scale the interpolation such that it goes from il to il+1 when xvar goes from 0 to 1.
  const double a = (fXYArmData[il+1].second - fXYArmData[il].second);
  const double y0 = fXYArmData[il].second;
  
  double xmin = ( x + a*y - a*y0 ) / ( 1 + POW2(a) );
  double ymin(0);
  
  if ( xmin <= 0 ) {
    
    xmin = fXYArmData[il].first;
    ymin = fXYArmData[il].second;
    
    
  } else if ( xmin >= 1 ) {
    
    xmin = fXYArmData[il+1].first;
    ymin = fXYArmData[il+1].second;
    
  } else {
    
    ymin = a*xmin + y0;
    xmin = fXYArmData[il].first + xmin*(fXYArmData[il+1].first - fXYArmData[il].first);
    
  }
  
  return sqrt( POW2(x-xmin) + POW2(y-ymin) );
}

void ArmFunction::findBestPoint(size_t &imin, double &sqrsmin, size_t ibegin, double rbegin, double rad, double x, double y) const
{

   const size_t pow2(3);
   size_t absStep = 1<<pow2;

   const int step = rbegin < rad ? absStep : -absStep;

   size_t i(ibegin),icmin(ibegin);  //Use underflow as a check for out of bounds

   //Start by finding the best index in steps of 8, then 4, then 2
   double cmin = std::numeric_limits<double>::max(); //A break criteria for the loop
   while (i < fXYArmData.size()) {

      const double sqrs = POW2(x-fXYArmData[i].first) + POW2(y-fXYArmData[i].second);
      
      if (sqrs < cmin) {
         cmin = sqrs;
         icmin = i;
      } else {
         break;
      }

      i += step;

   }

   //Now we test the points around the best one and half the search each time until the step is 2
   while (absStep > 2) {

      absStep = absStep >> 1; //Divide by 2;

      i = icmin + absStep;
      if ( i < fXYArmData.size() ) {

         const double sqrs = POW2(x-fXYArmData[i].first) + POW2(y-fXYArmData[i].second);

         if (sqrs < cmin) {
            cmin = sqrs;
            icmin = i;
            continue;
         } 

      }

      i = icmin - absStep;
      if ( i < fXYArmData.size() ) {

         const double sqrs = POW2(x-fXYArmData[i].first) + POW2(y-fXYArmData[i].second);

         if (sqrs < cmin) {
            cmin = sqrs;
            icmin = i;
            continue;
         } 

      }

   }

   //Update if needed
   if (cmin < sqrsmin) {
      sqrsmin = cmin;
      imin = icmin;
   }

}


double ArmFunction::armDistOld ( double rad, double phi ) const
{
  
  const double x = rad*cos(phi);
  const double y = rad*sin(phi);
  
  // Start by stepping in size of 10 and then search around that in steps of 2
  // Store the 2nd best value as well and explore
  double sqrtsmin = std::numeric_limits<double>::max();
  size_t imin(0);

  for (size_t i(0); i < fXYArmData.size(); i += 10) { 

    const double sqrs = POW2(x-fXYArmData[i].first) + POW2(y-fXYArmData[i].second);
    
    if (sqrs < sqrtsmin) {
      sqrtsmin = sqrs;
      imin = i;
    }
    
  }

  //Do the last bin to make sure the ends are sampled properly
  const double sqrs =  POW2(x-fXYArmData[fXYArmData.size()-1].first) + POW2(y-fXYArmData[fXYArmData.size()-1].second);
  if (sqrs < sqrtsmin) {
    sqrtsmin = sqrs;
    imin = fXYArmData.size()-1;
  }
  
  //Create the search limits
  const size_t ilow = imin < 10 ? 0 : imin-9;
  const size_t ihigh = imin + 10 >= fXYArmData.size() ? fXYArmData.size() : imin + 10;


  for (size_t i(ilow); i < ihigh; i += 2) {  //Step by 2 to speed things up

    const double sqrs = POW2(x-fXYArmData[i].first) + POW2(y-fXYArmData[i].second);
    
    if (sqrs < sqrtsmin) {
      sqrtsmin = sqrs;
      imin = i;
    }
    
  }
  
  
  // Find the nearest neighbour with smaller values and use linear interpolation between the points
  const double sqrsl = (imin == 0) ? std::numeric_limits< double >::max() : POW2(x-fXYArmData[imin-1].first) + POW2(y-fXYArmData[imin-1].second);
  const double sqrsu = (imin == fXYArmData.size()-1) ? std::numeric_limits< double >::max() : POW2(x-fXYArmData[imin+1].first) + POW2(y-fXYArmData[imin+1].second);
  
  const size_t il = (sqrsl < sqrsu) ? imin-1 : imin;
  
  // Scale the interpolation such that it goes from il to il+1 when xvar goes from 0 to 1.
  const double a = (fXYArmData[il+1].second - fXYArmData[il].second);
  const double y0 = fXYArmData[il].second;
  
  double xmin = ( x + a*y - a*y0 ) / ( 1 + POW2(a) );
  double ymin(0);
  
  if ( xmin <= 0 ) {
    
    xmin = fXYArmData[il].first;
    ymin = fXYArmData[il].second;
    
    
  } else if ( xmin >= 1 ) {
    
    xmin = fXYArmData[il+1].first;
    ymin = fXYArmData[il+1].second;
    
  } else {
    
    ymin = a*xmin + y0;
    xmin = fXYArmData[il].first + xmin*(fXYArmData[il+1].first - fXYArmData[il].first);
    
  }
  
  return sqrt( POW2(x-xmin) + POW2(y-ymin) );

}


double ArmFunction::armDistBF ( double radius, double phi ) const
{
  
  const double x = radius*cos(phi);
  const double y = radius*sin(phi);

  
  // Use brute force method to find minimum
  double sqrtsmin = std::numeric_limits<double>::max();
  size_t imin(0);

  for (size_t i(0); i < fXYArmData.size(); i += 2) {  //Step by 2 to speed things up
  
    const double sqrs = POW2(x-fXYArmData[i].first) + POW2(y-fXYArmData[i].second);
    
    if (sqrs < sqrtsmin) {
      sqrtsmin = sqrs;
      imin = i;
    }
    
  }

  
  // Find the nearest neighbour with smaller values and use linear interpolation between the points
  const double sqrsl = (imin == 0) ? std::numeric_limits< double >::max() : POW2(x-fXYArmData[imin-1].first) + POW2(y-fXYArmData[imin-1].second);
  const double sqrsu = (imin == fXYArmData.size()-1) ? std::numeric_limits< double >::max() : POW2(x-fXYArmData[imin+1].first) + POW2(y-fXYArmData[imin+1].second);

  const size_t il = (sqrsl < sqrsu) ? imin-1 : imin;

  // Scale the interpolation such that it goes from il to il+1 when xvar goes from 0 to 1.
  const double a = (fXYArmData[il+1].second - fXYArmData[il].second);
  const double y0 = fXYArmData[il].second;
  
  double xmin = ( x + a*y - a*y0 ) / ( 1 + POW2(a) );
  double ymin(0);
  
  if ( xmin <= 0 ) {

    xmin = fXYArmData[il].first;
    ymin = fXYArmData[il].second;
  
    
  } else if ( xmin >= 1 ) {
    
    xmin = fXYArmData[il+1].first;
    ymin = fXYArmData[il+1].second;

  } else {
    
    ymin = a*xmin + y0;
    xmin = fXYArmData[il].first + xmin*(fXYArmData[il+1].first - fXYArmData[il].first);
    
  }

    
  return sqrt( POW2(x-xmin) + POW2(y-ymin) );
}


double ArmFunction::a() const
{
  return fA;
}


double ArmFunction::angle ( const double radius ) const
{
  return fA*std::log(radius/fRMin) + fPhiMin;
}

double ArmFunction::radius ( double angle ) const
{
  return fRMin * exp( ( angle - fPhiMin ) / fA );
}

double ArmFunction::phiExtent() const
{
  return fPhiExtent;
}

double ArmFunction::phiMin() const
{
  return fPhiMin;
}

double ArmFunction::rMin() const
{
  return fRMin;
}

double ArmFunction::rMax() const
{
   return radius(fPhiExtent + fPhiMin);
}

double ArmFunction::width() const
{
  return fWidth;
}

void ArmFunction::setA ( double a )
{
  fA = a;
  fCosAlpha = cos(atan(1./a));
  calculateRPdata();
  calculateXYdata();
}

void ArmFunction::setPhiExtent ( double phiExtent )
{
  fPhiExtent = phiExtent;
  calculateRPdata();
  calculateXYdata();
}

void ArmFunction::setPhiExtentFromRMax ( double rMax )
{
   if (rMax <= fRMin) {
      std::cerr<<"Cannot have rMax <= rMin!!!"<<std::endl;
   } else {
      fPhiExtent = angle(rMax) - fPhiMin;
      calculateRPdata();
      calculateXYdata();
   }
}

void ArmFunction::setPhiMin ( double phiMin )
{
  fPhiMin = phiMin;
  calculateRPdata();
  calculateXYdata();
}

void ArmFunction::setRMin ( double rMin )
{
  fRMin = rMin;
  calculateRPdata();
  calculateXYdata();
}

void ArmFunction::setWidth ( double width )
{
  fWidth = width;
  fOneOverWidth = 1./fWidth;
}

const std::vector< std::pair< double, double > >& ArmFunction::RPArmData() const
{
  return fRPArmData;
}

void ArmFunction::setRPArmData ( const std::vector< std::pair< double, double > >& RPArmData )
{
   //Now use armdist, because that allows arbitrary arm shape
   armDistPointer = &GalacticStructure::ArmFunction::armDist;
  fRPArmData = RPArmData;
  
  // Find out how many steps per ring, assuming linear step size.
  nStepsPerRing = 2*M_PI/fabs(RPArmData[1].second - RPArmData[0].second);
  
  calculateXYdata();
}

#ifdef HAVE_OPENCL
std::string ArmFunction::getOpenCLFunction(std::string &name, size_t parIndex) const
{
   //fix the name
   name += "LogSpiral_"+name;

   //Get the scale function
   std::string scName = name;
   std::string scCode = fscale->getOpenCLFunction(scName);

   //Now the function
   std::ostringstream codeStr;
   codeStr<<scCode
      <<"inline float "<<name<<"_radius(const float t, __constant float* pars)\n"
      <<"{\n"
      <<"  return pars["<<parIndex+3<<"]*exp( (t - pars["<<parIndex+1<<"])*pars["<<parIndex<<"]);\n"
      <<"}\n"
      <<"float "<<name<<"(const float r, const float t, __constant float* pars)\n"
      <<"{\n"
      <<"  float dmin = MAXFLOAT;\n"
      <<"  float th, d;"
      <<"  for (th = t; th < pars["<<parIndex+1<<"]+pars["<<parIndex+2<<"]; th += 2*M_PI_F) \n"
      <<"  {\n"
      <<"    if (th < pars["<<parIndex+1<<"]) continue;\n"
      <<"    d = fabs(r-"<<name<<"_radius(th,pars))*pars["<<parIndex+5<<"];\n"
      <<"    if (d < dmin) dmin = d;\n"
      <<"  }\n"
      <<"  for (th = t-M_PI_F; th < pars["<<parIndex+1<<"]+pars["<<parIndex+2<<"]; th += 2*M_PI_F) \n"
      <<"  {\n"
      <<"    if (th < pars["<<parIndex+1<<"]) continue;\n"
      <<"    d = (r+"<<name<<"_radius(th,pars))*pars["<<parIndex+5<<"];\n"
      <<"    if (d < dmin) dmin = d;\n"
      <<"  }\n"
      <<"  return "<<scName<<"(dmin*pars["<<parIndex+4<<"]);\n"
      <<"}\n";

   return codeStr.str();
}

size_t ArmFunction::getOpenCLNPars() const
{
   return 6;
}

std::vector<cl_float> ArmFunction::getOpenCLPars() const
{
   std::vector<cl_float> pars(6);
   pars[0] = 1./fA;
   pars[1] = fPhiMin;
   pars[2] = fPhiExtent;
   pars[3] = fRMin;
   pars[4] = fOneOverWidth;
   pars[5] = fCosAlpha;
   return pars;
}
#endif

}
