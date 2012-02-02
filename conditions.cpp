#include <Fepic/CustomEigen>
#include <cmath>
#include <iostream>
using namespace Eigen;
using namespace std;

// space
typedef Matrix<double, Dynamic,1,0,3,1>              Vector;
typedef Matrix<double, Dynamic,Dynamic,RowMajor,3,3> Tensor;

const double pi  = 3.141592653589793;
const double pi2 = pi*pi;
inline double sqr(double v) {return v*v;}

double pho(Vector const& X, int tag);
double gama(Vector const& X, double t, int tag);
double niu(double t, int tag);
Vector force(Vector const& X, double t, int tag);
Vector u_boundary(Vector const& X, double t, int tag);
Vector traction(Vector const& X, double t, int tag);
double pressure_exact(Vector const& X, double t, int tag);
Vector grad_p_exact(Vector const& X, double t, int tag);
Tensor grad_u_exact(Vector const& X, double t, int tag);
Vector u_initial(Vector const& X, int tag);
double p_initial(Vector const& X, int tag);
Vector solid_normal(Vector const& X, double t, int tag);

#define PROBLEM_TYPE 8


#define CAVITY_2D_3D 1 // navier-stokes
#define KOVASZNAY_2D 2 // navier-stokes
#define STOKES_JD_2D 3 // see Jean-Donea
#define STATIC_BB    4 // static bubble
#define SLOSH_2D     5
#define ANYTHING_3D  6
#define OSCI_BB      7 // osci bubble
#define RAMP2D3D     8
#define RUSSA2D      9
#define ANGLE2D      10
#define RAMP3D       11

#if (PROBLEM_TYPE==CAVITY_2D_3D)
double pho(Vector const& X, int tag)
{
  return 1.0;
}
double cos_theta0()
{
  return 0.5;
}

double zeta(double u_norm, double angle)
{
  return 5.e-1;
}

double beta_diss()
{
  return 1.e-4;
}

double gama(Vector const& X, double t, int tag)
{
  return 1;
}
double niu(double t, int tag)
{
  return 1;
}
Vector force(Vector const& X, double t, int tag)
{
  Vector r(X.size());
  r.setZero();
  
  return r;
}
Vector u_boundary(Vector const& X, double t, int tag)
{
  Vector r(X.size());
  r.setZero();
  
  if (tag != 1)
    return r;
  r(0) = 1;
    return r;
}
Vector traction(Vector const& X, double t, int tag)
{
  Vector r(X.size());
  r.setZero();
  
  return r;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  return 1;
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  Vector r(X.size());
  r.setZero();
  
  return r;
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  Tensor r(X.size(), X.size());
  r.setZero();
  
  return r;
}
Vector u_initial(Vector const& X, int tag)
{
  Vector r(X.size());
  r.setZero();
  
  return r;
}
double p_initial(Vector const& X, int tag)
{
  return 1;
}
#endif

#if (PROBLEM_TYPE==KOVASZNAY_2D)
double pho(Vector const& X, int tag)
{
  return 1.0;
}
double cos_theta0()
{
  return 0.5;
}

double zeta(double u_norm, double angle)
{
  return 0*5.e-1;
}

double beta_diss()
{
  return 0*1.e-4;
}

double gama(Vector const& X, double t, int tag)
{
  return 1;
}
double niu(double t, int tag)
{
  return  0.025;
}
Vector force(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
Vector u_boundary(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(2));  
  
  const double a = 1./(2.*niu(t,tag)) - sqrt(1./(4.*sqr(niu(t,tag))) + 4.*pi2);
  v(0) = 1-exp(a*x)*cos(2.*pi*y);
  v(1) = a*exp(a*x)*sin(2.*pi*y)/(2.*pi);
  return v;
}
Vector traction(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector T(Vector::Zero(X.size()));

  const double a = 1./(2.*niu(t,tag)) - sqrt(1./(4.*sqr(niu(t,tag))) + 4.*pi2);
  T(0) = -2*a*niu(t,tag)*exp(a*x)*cos(2*pi*y) - (1-exp(2*a*x))/2.;
  T(1) = (a*a*niu(t,tag)*exp(a*x)*sin(2*pi*y))/(2*pi) + 2*pi*niu(t,tag)*exp(a*x)*sin(2*pi*y);
  return T;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  const double a = 1./(2.*niu(t,tag)) - sqrt(1./(4.*sqr(niu(t,tag))) + 4.*pi2);
  //return  0.5*(1.-exp(2.*a*x));
  return  (1.-exp(2.*a*x))/2.;
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(Vector::Zero(2));
  
  const double a = 1./(2.*niu(t,tag)) - sqrt(1./(4.*sqr(niu(t,tag))) + 4.*pi2);
  dxP(0) = -a*exp(2.*a*x);
  return dxP;
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(), X.size()));

  const double a = 1./(2.*niu(t,tag)) - sqrt(1./(4.*sqr(niu(t,tag))) + 4.*pi2);
  dxU(0,0) = -a*exp(a*x)*cos(2*pi*y);
  dxU(0,1) = 2*pi*exp(a*x)*sin(2*pi*y);
  dxU(1,0) = (pow(a,2)*exp(a*x)*sin(2*pi*y))/(2*pi);
  dxU(1,1) = a*exp(a*x)*cos(2*pi*y);
  return dxU;
}
Vector u_initial(Vector const& X, int tag)
{
  return Vector::Zero(X.size());
}
double p_initial(Vector const& X, int tag)
{
  return 0;
}
 
#endif

#if (PROBLEM_TYPE==STOKES_JD_2D)
double pho(Vector const& X, int tag)
{
  return 1.0;
}
double cos_theta0()
{
  return 0.5;
}

double zeta(double u_norm, double angle)
{
  return 0*5.e-1;
}

double beta_diss()
{
  return 0*1.e-4;
}

double gama(Vector const& X, double t, int tag)
{
  return 1;
}
double niu(double t, int tag)
{
  return 1;
}
Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector f(2);

  // 2d test - jean donea
  double x2 = x*x, x3=x*x*x, x4=x*x*x*x;
  double y2 = y*y, y3=y*y*y, y4=y*y*y*y;
  f(0) = (-48.*x2+48.*x-8.)*y3+(72.*x2-72.*x+12.)*y2+(-24.*x4+48.*x3-48.*x2+24.*x-4.)*y+12.*x4-24.*x3+12.*x2-2.*x+1.;
  f(1) = (24.*x-12.)*y4+(24.-48.*x)*y3+(48.*x3-72.*x2+48.*x-12.)*y2+(-48.*x3+72.*x2-24.*x)*y+8.*x3-12.*x2+4.*x;
  return f;
}
Vector u_boundary(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(2);

  // 2d test - jean donea
  double x2 = x*x, x3=x*x*x, x4=x*x*x*x;
  double y2 = y*y, y3=y*y*y, y4=y*y*y*y;
  v(0) = x2*(1.-x)*(1.-x)*(2.*y-6.*y2+4.*y3);
  v(1) = -y2*(1.-y)*(1.-y)*(2.*x-6.*x2+4.*x3);
  return v;
}
Vector traction(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector T(2);

  // in tag = 1

  // jean donea
  double x2 = x*x, x3=x*x*x, x4=x*x*x*x;
  double y2 = y*y, y3=y*y*y, y4=y*y*y*y;
  T(0) = pow(1.-x,2)*x2*(12.*y2-12.*y+2.)-(12.*x2-12.*x+2.)*pow(1.-y,2)*y2;
  T(1) = 4.*(4.*x3-6.*x2+2.*x)*(1.-y)*y2-4.*(4*x3-6.*x2+2.*x)*pow(1.-y,2)*y-(1.-x)*x;
  return T;

}
double pressure_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return x*(1-x) + 666;
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);  
  Vector r(2);
  
  r(0) = 1.-2.*x;
  r(1) = 0;
  
  
  return r;
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(2,2);

  //2d test - jean donea
  double x2 = x*x, x3=x*x*x, x4=x*x*x*x;
  double y2 = y*y, y3=y*y*y, y4=y*y*y*y;
  dxU(0,0) = 2.*(1.-x)*(1.-x)*x*(4.*y3-6.*y2+2.*y)-2.*(1.-x)*x2*(4.*y3-6.*y2+2.*y);
  dxU(0,1) = (1.-x)*(1.-x)*x2*(12.*y2-12.*y+2.);
  dxU(1,0) = -(12.*x2-12.*x+2.)*(1.-y)*(1.-y)*y2;
  dxU(1,1) = 2.*(4.*x3-6.*x2+2.*x)*(1.-y)*y2-2.*(4.*x3-6.*x2+2.*x)*(1.-y)*(1.-y)*y;
  return dxU;

}
Vector u_initial(Vector const& X, int tag)
{
  Vector r(X.size());
  r.setZero();
  
  return r;
}
double p_initial(Vector const& X, int tag)
{
  return 1;
}
#endif

#if (PROBLEM_TYPE==STATIC_BB)
double pho(Vector const& X, int tag)
{
  return 0.0;
}
double cos_theta0()
{
  return 0.5;
}

double zeta(double u_norm, double angle)
{
  return 0*5.e-1;
}

double beta_diss()
{
  return 0*1.e-4;
}

double gama(Vector const& X, double t, int tag)
{
  return 1;
}
double niu(double t, int tag)
{
  return 1;
}
Vector force(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
Vector u_boundary(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
Vector traction(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
double pressure_exact(Vector const& X, double t, int tag)
{
  return X.size()==2 ? 1 : 2;
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  return Tensor::Zero(X.size(),X.size());
}
Vector u_initial(Vector const& X, int tag)
{
  Vector r(X.size());
  r.setZero();
  
  return r;
}
double p_initial(Vector const& X, int tag)
{
  return 2;
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N = X;
  N(0) = 1;
  N(1) = 1;
  N(2) = 1;
  N.normalize();
  
  return N;
}
#endif

#if (PROBLEM_TYPE==SLOSH_2D)
double pho(Vector const& X, int tag)
{
  return 1.0;
}
double cos_theta0()
{
  return 0.;
}

double zeta(double u_norm, double angle)
{
  return 0*5.e-1;
}

double beta_diss()
{
  return 0*1.e-4;
}

double gama(Vector const& X, double t, int tag)
{
  return 0;
}
double niu(double t, int tag)
{
  return 0.01;
}
Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector f(Vector::Zero(X.size()));
  f(1)=-1;
  return f;
}
Vector u_boundary(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(2));

  return v;
}
Vector traction(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector T(Vector::Zero(X.size()));
  return T;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return 2;
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);  
  Vector r(Vector::Zero(X.size()));
  
  return r;
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(),X.size()));

  return dxU;
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(X.size());
  if (tag == 3 || tag == 20)
    N << 1,0;
  else
  if (tag == 4 || tag == 21)
    N << -1,0;
  else
  {
    cout << "TAGGG ERROR: " << tag << endl;
    throw;
  }
  return N;
}


Vector u_initial(Vector const& X, int tag)
{
  Vector r(X.size());
  r.setZero();
  
  return r;
}
double p_initial(Vector const& X, int tag)
{
  return 2;
}
#endif

#if (PROBLEM_TYPE==ANYTHING_3D)
double pho(Vector const& X, int tag)
{
  return 1.0;
}

double cos_theta0()
{
  return 0.5;
}

double zeta(double u_norm, double angle)
{
  return 0*5.e-1;
}

double beta_diss()
{
  return 0*1.e-4;
}

double gama(Vector const& X, double t, int tag)
{
  return 1;
}
double niu(double t, int tag)
{
  return 0.025;
}
Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double z = X(2);
  Vector f(3);
  f(0) = exp(-y*z-x*z-y-x)*(2*pi*y*z*exp(y*z+x*z+y+x)*cos(x*y*z)+(niu(t,tag)*exp(y+x)*pow(z,2)+niu(t,tag)*pow(y,2)*exp(y+x)+y)*exp(x*z)+exp(y+x)*z);
  f(1) = exp(-y*z-x*z-y-x)*(2*pi*x*z*exp(y*z+x*z+y+x)*cos(x*y*z)+(-niu(t,tag)*exp(y+x)*pow(z,2)-niu(t,tag)*pow(x,2)*exp(y+x)-x)*exp(y*z)+exp(y+x)*z);
  f(2) = exp(-y*z-x*z-y-x)*(2*pi*x*y*exp(y*z+x*z+y+x)*cos(x*y*z)+(-2*niu(t,tag)*exp(x*z)-1)*exp(y*z)+exp(x*z));

  return f;
}
Vector u_boundary(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double z = X(2);
  Vector v(3);

  v(0) = -exp(-y*z);
  v(1) = exp(-x*z);
  v(2) = exp(-x-y);

  //v(0) = y;
  //v(1) = 0;
  //v(2) = 0;

  return v;
}
Vector traction(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double z = X(2);
  Vector T(Vector::Zero(X.size()));

  T(0) = niu(t,tag)*z*exp(-y*z)-niu(t,tag)*z*exp(-x*z);
  T(1) = -2*pi*sin(x*y*z);
  T(2) = -niu(t,tag)*x*exp(-x*z)-niu(t,tag)*exp(-y-x);  
  
  return T;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double z = X(2);

  return 2.*pi*sin(x*y*z);
  //return 0;
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double z = X(2);
  Vector r(3);

  r(0) = 2.*pi*y*z*cos(x*y*z);
  r(1) = 2.*pi*x*z*cos(x*y*z);
  r(2) = 2.*pi*x*y*cos(x*y*z);
  
  return r;
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double z = X(2);
  Tensor dxU(3,3);

  dxU(0,0) = 0;
  dxU(0,1) = z*exp(-y*z);
  dxU(0,2) = y*exp(-y*z);
  dxU(1,0) = -z*exp(-x*z);
  dxU(1,1) = 0;
  dxU(1,2) = -x*exp(-x*z);
  dxU(2,0) = -exp(-y-x);
  dxU(2,1) = -exp(-y-x);
  dxU(2,2) = 0;

  return dxU;
}
Vector u_initial(Vector const& X, int tag)
{
  Vector r(X.size());
  r.setZero();
  
  return r;
}
double p_initial(Vector const& X, int tag)
{
  return 2;
}
#endif

#if (PROBLEM_TYPE==OSCI_BB)
double pho(Vector const& X, int tag)
{
  return 1.0;
}

double cos_theta0()
{
  return 0*0.5;
}

double zeta(double u_norm, double angle)
{
  return 0*5.e-1;
}

double beta_diss()
{
  return 0*1.e-4;
}

double gama(Vector const& X, double t, int tag)
{
  return 1;
}
double niu(double t, int tag)
{
  return 0.01;
}
Vector force(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
Vector u_boundary(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
Vector traction(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
double pressure_exact(Vector const& X, double t, int tag)
{
  return X.size()==2 ? 1 : 2;
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  return Tensor::Zero(X.size(),X.size());
}
Vector u_initial(Vector const& X, int tag)
{
  Vector r(X.size());
  r.setZero();
  
  return r;
}
double p_initial(Vector const& X, int tag)
{
  return 2;
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(X.size());
  if (tag == 3)
    N << 1,0;
  else
  if (tag == 4)
    N << 0,1;
  else
  if (tag == 5)
    N << -1,0;
  
  N = X;
  N.normalize();
  
  return N;
}


#endif

#if (PROBLEM_TYPE==RAMP2D3D)
double pho(Vector const& X, int tag)
{
  return 0.1;
}

double cos_theta0()
{
  return -sqrt(2)/2;
}

double zeta(double u_norm, double angle)
{
  return 0*1.e-2;
}

double beta_diss()
{
  return 0*1.e-2;
}

double gama(Vector const& X, double t, int tag)
{
  return 1.;
}
double niu(double t, int tag)
{
  return 1;
}
Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector f(Vector::Zero(X.size()));
  //f(1)=-1;
  //f(0)=1;
  //f /= 4*sqrt(2);
  
  
  return f;
}
Vector u_boundary(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}
Vector traction(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector T(Vector::Zero(X.size()));
  return T;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return 2;
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);  
  Vector r(Vector::Zero(X.size()));
  
  return r;
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(),X.size()));

  return dxU;
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  //if (true) //(X(0) <=1 )//|| X(0) >= 5)
  //{
    //if (tag == 3 || tag == 20 || tag == 21)
    //{
      //N << 1,1;
      //N /= sqrt(2.);
    //}
    //else
    //{
      //cout << "TAGGG ERROR: " << tag << endl;
      //throw;
    //}
  //}
  //else
  //{
    //if (tag == 3 || tag == 20 || tag == 21)
    //{
      ////N << -2.*X(0)/3.-7/3.,1;
      //N << 1-(X(0)-1),1; 
      //N.normalize();
    //}
  //}  
  
  //N << 0,1;
  N(1) = 1;
  return N;
}

Vector u_initial(Vector const& X, int tag)
{
  Vector r(X.size());
  r.setZero();
  
  return r;
}
double p_initial(Vector const& X, int tag)
{
  return 2;
}
#endif

#if (PROBLEM_TYPE==RUSSA2D)
double pho(Vector const& X, int tag)
{
  return 1.0;
}

double cos_theta0()
{
  return 0.5;
}

double zeta(double u_norm, double angle)
{
  return 1.e-1;
}

double beta_diss()
{
  return 1.e-5;
}

double gama(Vector const& X, double t, int tag)
{
  return 1;
}
double niu(double t, int tag)
{
  return 1.0;
}
Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector f(Vector::Zero(X.size()));
  f(1)=-1;
  f(0)=0;
  //f /= sqrt(2);
  return f;
}
Vector u_boundary(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(2));

  return v;
}
Vector traction(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector T(Vector::Zero(X.size()));
  return T;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return 2;
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);  
  Vector r(Vector::Zero(X.size()));
  
  return r;
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(),X.size()));

  return dxU;
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  double const alpha = 0.5;
  double const beta  = 0.7;
  double const w = 1;
  int const k = 1;
  
  double const b = acos(alpha/(w*beta))/w + 2.*k*pi;
  
  Vector N(X.size());
  
  
  if (tag == 3 || tag == 20 || tag == 21)
  {
    double const s = (X(0)<=b? -1 : +1);
    double const x = s<0 ? X(0) : 2.*b - X(0);
    
    N(0) =  s*(   -alpha + (beta-0.00*x*cos(w*x)) *w*cos(w*x) +  sin(w*x)*0.00*(x*w*sin(w*x) - cos(w*x))  );
    N(1) = 1;
    N.normalize();
  }
  else
  {
    cout << "TAGGG ERROR: " << tag << endl;
    throw;
  }

  
  //N << 0,1;
  return N;
}

Vector u_initial(Vector const& X, int tag)
{
  Vector r(X.size());
  r.setZero();
  
  return r;
}
double p_initial(Vector const& X, int tag)
{
  return 2;
}
#endif


#if (PROBLEM_TYPE==ANGLE2D)
double pho(Vector const& X, int tag)
{
  return 1.0;
}

double cos_theta0()
{
  return 0.5;
}

double zeta(double u_norm, double angle)
{
  return 0*1.e-1;
}

double beta_diss()
{
  return 0*1.e-1;
}

double gama(Vector const& X, double t, int tag)
{
  return 1;
}
double niu(double t, int tag)
{
  return 0.1;
}
Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector f(Vector::Zero(X.size()));
  f(1)=0;
  f(0)=0;
  //f /= sqrt(2);
  return f;
}
Vector u_boundary(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(2));

  return v;
}
Vector traction(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector T(Vector::Zero(X.size()));
  return T;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return 2;
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);  
  Vector r(Vector::Zero(X.size()));
  
  return r;
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(),X.size()));

  return dxU;
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(X.size());
  //if (X(0) <=1 )//|| X(0) >= 5)
  //{
  //  if (tag == 3 || tag == 20 || tag == 21)
  //  {
  //    N << 1,1;
  //    N /= sqrt(2.);
  //  }
  //  else
  //  {
  //    cout << "TAGGG ERROR: " << tag << endl;
  //    throw;
  //  }
  //}
  //else
  //{
  //  if (tag == 3 || tag == 20 || tag == 21)
  //  {
  //    //N << -2.*X(0)/3.-7/3.,1;
  //    N << 1-(X(0)-1),1; 
  //    N.normalize();
  //  }
  //}  
  
  N << 0,1;
  return N;
}

Vector u_initial(Vector const& X, int tag)
{
  Vector r(X.size());
  r.setZero();
  
  return r;
}
double p_initial(Vector const& X, int tag)
{
  return 2;
}
#endif









