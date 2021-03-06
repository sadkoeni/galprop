#include <vector>
#include <string>

void B_field_3D_model
(const std::string &name, const std::vector<double> &parameters,
 double x, double y, double z,
 int options,
 double &Breg,  double &Bregx, double &Bregy, double &Bregz,
 double &Bran,  double &Branx, double &Brany, double &Branz,
 int debug=0);

double Bperp(double x,double y,double z, double Bx,double By,double Bz, double x0,double y0,double z0);


typedef enum {B_TOTAL,B_RANDOM,B_REGULAR} B_field_component;

double B_field_3D_model_tot
(const std::string &name, const std::vector<double> &parameters,
 double x, double y, double z,
 int debug=0 );

double B_field_3D_model_abs
(const std::string &name, const std::vector<double> &parameters,
 double x, double y, double z,
 int debug=0, B_field_component bcomp=B_TOTAL );


double B_field_3D_model_tot
(const std::string &name, const std::vector<double> &parameters,
 double R, double z,
 int debug=0 );

double B_field_3D_model_abs
(const std::string &name, const std::vector<double> &parameters,
 double R, double z,
 int debug=0, B_field_component bcomp=B_TOTAL  );    
