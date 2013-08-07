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
double muu(int tag);
Vector force(Vector const& X, double t, int tag);
Vector u_exact(Vector const& X, double t, int tag);
Vector traction(Vector const& X, Vector const& normal, double t, int tag);
double pressure_exact(Vector const& X, double t, int tag);
Vector grad_p_exact(Vector const& X, double t, int tag);
Tensor grad_u_exact(Vector const& X, double t, int tag);
Vector u_initial(Vector const& X, int tag);
double p_initial(Vector const& X, int tag);
Vector solid_normal(Vector const& X, double t, int tag);
Tensor feature_proj(Vector const& X, double t, int tag);

inline double sign(double a) {a<0 ? -1 : 1;};

#define CAVITY_2D_3D   1 // navier-stokes
#define KOVASZNAY_2D   2 // navier-stokes
#define STOKES_JD_2D   3 // see Jean-Donea
#define STATIC_BB      4 // static bubble
#define SLOSH_2D       5
#define ANYTHING_3D    6
#define OSCI_BB        7 // osci bubble
#define RAMP2D3D       8
#define RUSSA2D        9
#define ANGLE2D        10
#define RAMP3D         11
#define COUETTE        12
#define PSEUDO_OSC2D   13
#define QUASI_STOKES   14
#define GERBEAU        15
#define ZETA3D         16
#define MICRO2D        17
#define TRACO          18
#define RUSSA_SIN2D    19

#define PROBLEM_TYPE 16


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
double muu(int tag)
{
  return 1;
}
Vector force(Vector const& X, double t, int tag)
{
  Vector r(X.size());
  r.setZero();

  return r;
}
Vector u_exact(Vector const& X, double t, int tag)
{
  Vector r(X.size());
  r.setZero();

  if (tag != 1)
    return r;
  r(0) = 1;
    return r;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
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
Vector v_exact(Vector const& , double , int ) //(X,t,tag)
{

}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  N(1) = 1;
  return N;
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
double muu(int tag)
{
  return  0.025;
}
Vector force(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(2));

  const double a = 1./(2.*muu(tag)) - sqrt(1./(4.*sqr(muu(tag))) + 4.*pi2);
  v(0) = 1-exp(a*x)*cos(2.*pi*y);
  v(1) = a*exp(a*x)*sin(2.*pi*y)/(2.*pi);
  return v;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector T(Vector::Zero(X.size()));

  const double a = 1./(2.*muu(tag)) - sqrt(1./(4.*sqr(muu(tag))) + 4.*pi2);
  T(0) = -2*a*muu(tag)*exp(a*x)*cos(2*pi*y) - (1-exp(2*a*x))/2.;
  T(1) = (a*a*muu(tag)*exp(a*x)*sin(2*pi*y))/(2*pi) + 2*pi*muu(tag)*exp(a*x)*sin(2*pi*y);
  return T;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  const double a = 1./(2.*muu(tag)) - sqrt(1./(4.*sqr(muu(tag))) + 4.*pi2);
  //return  0.5*(1.-exp(2.*a*x));
  return  (1.-exp(2.*a*x))/2.;
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(Vector::Zero(2));

  const double a = 1./(2.*muu(tag)) - sqrt(1./(4.*sqr(muu(tag))) + 4.*pi2);
  dxP(0) = -a*exp(2.*a*x);
  return dxP;
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(), X.size()));

  const double a = 1./(2.*muu(tag)) - sqrt(1./(4.*sqr(muu(tag))) + 4.*pi2);
  dxU(0,0) = -a*exp(a*x)*cos(2*pi*y);
  dxU(0,1) = 2*pi*exp(a*x)*sin(2*pi*y);
  dxU(1,0) = (pow(a,2)*exp(a*x)*sin(2*pi*y))/(2*pi);
  dxU(1,1) = a*exp(a*x)*cos(2*pi*y);
  return dxU;
}
Vector u_initial(Vector const& X, int tag)
{
  return 0*u_exact(X,0,tag);
}
double p_initial(Vector const& X, int tag)
{
  return 0;
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));

  return N;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  N(1) = 1;
  return N;
}


Vector v_exact(Vector const& X, double t, int ) //(X,t,tag)
{

  Vector v(X.size());
  v(0) = 0*sin(pi*t)*(X(0)+0.5)*(X(0)-1.)/2.;
  v(1) = 0*cos(pi*t)*(X(1)+0.5)*(X(1)-.5)/2.;
  return v;
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
double muu(int tag)
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
Vector u_exact(Vector const& X, double t, int tag)
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
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
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

  return x*(1-x);
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

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));

  return N;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  N(1) = 1;
  return N;
}


Vector v_exact(Vector const& X, double , int ) //(X,t,tag)
{
  return Vector::Zero(X.size());
}
#endif

#if (PROBLEM_TYPE==STATIC_BB)
double pho(Vector const& X, int tag)
{
  return 1.0;
}
double cos_theta0()
{
  return 0.0;
}

double zeta(double u_norm, double angle)
{
  return 0;
}

double beta_diss()
{
  return 0;
}


double gama(Vector const& X, double t, int tag)
{
  return 1;
}
double muu(int tag)
{
  return 0.001;
}
Vector force(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
Vector u_exact(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  return -gama(X,t,tag)*normal*(X.size()==2 ? 1 : 2);
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

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  N(1) = 1;
  return N;
}


Vector v_exact(Vector const& X, double , int ) //(X,t,tag)
{
  return Vector::Zero(X.size());
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
double muu(int tag)
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
Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(2));

  return v;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
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
Vector v_exact(Vector const& , double , int ) //(X,t,tag)
{

}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  N(1) = 1;
  return N;
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
double muu(int tag)
{
  return 0.025;
}
Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double z = X(2);
  Vector f(3);
  f(0) = exp(-y*z-x*z-y-x)*(2*pi*y*z*exp(y*z+x*z+y+x)*cos(x*y*z)+(muu(tag)*exp(y+x)*pow(z,2)+muu(tag)*pow(y,2)*exp(y+x)+y)*exp(x*z)+exp(y+x)*z);
  f(1) = exp(-y*z-x*z-y-x)*(2*pi*x*z*exp(y*z+x*z+y+x)*cos(x*y*z)+(-muu(tag)*exp(y+x)*pow(z,2)-muu(tag)*pow(x,2)*exp(y+x)-x)*exp(y*z)+exp(y+x)*z);
  f(2) = exp(-y*z-x*z-y-x)*(2*pi*x*y*exp(y*z+x*z+y+x)*cos(x*y*z)+(-2*muu(tag)*exp(x*z)-1)*exp(y*z)+exp(x*z));

  return f;
}
Vector u_exact(Vector const& X, double t, int tag)
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
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  double z = X(2);
  Vector T(Vector::Zero(X.size()));

  T(0) = muu(tag)*z*exp(-y*z)-muu(tag)*z*exp(-x*z);
  T(1) = -2*pi*sin(x*y*z);
  T(2) = -muu(tag)*x*exp(-x*z)-muu(tag)*exp(-y-x);

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
Vector v_exact(Vector const& X, double , int ) //(X,t,tag)
{
  return Vector::Zero(X.size());
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));

  return N;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  N(1) = 1;
  return N;
}


#endif

#if (PROBLEM_TYPE==OSCI_BB)
double pho(Vector const& X, int tag)
{
  return 1.0;
}

double cos_theta0()
{
  return 0;
}

double zeta(double u_norm, double angle)
{
  return 0;
}

double beta_diss()
{
  return 0;
}

double gama(Vector const& X, double t, int tag)
{
  return 1;
}
double muu(int tag)
{
  return 0.01;
}
Vector force(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
Vector u_exact(Vector const& X, double t, int tag)
{
  return Vector::Zero(X.size());
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
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
  Vector N(Vector::Zero(X.size()));

  if (tag==3)
    N(0) = -1;
  else
    N(1) = -1;

  return N;
}


Vector v_exact(Vector const& X, double , int ) //(X,t,tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}


#endif

#if (PROBLEM_TYPE==RAMP2D3D)
double pho(Vector const& X, int tag)
{
  return 0*2.e-7;
}

double cos_theta0()
{
  return sqrt(2.)/2.;
}

double zeta(double u_norm, double angle)
{
  return 0*gama(Vector::Zero(3),0,0)*cos_theta0()/4200.;
}

double beta_diss()
{
  return 1.e-5;
}

double gama(Vector const& X, double t, int tag)
{
  return 0.075;
}
double muu(int tag)
{
  return 1.e-5;
}
Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector f(Vector::Zero(X.size()));
  
  //f(1)=4;
  //f(0)=1;
  //f /= 4*sqrt(2);


  return f;
}
Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector T(Vector::Zero(X.size()));

  //T = -X;

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

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
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
  return X.size()-1;
}

Vector v_exact(Vector const& X, double , int ) //(X,t,tag)
{
  Vector v(X.size());
  v.setZero();
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
double muu(int tag)
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
Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(2));

  return v;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
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
Vector v_exact(Vector const& , double , int ) //(X,t,tag)
{

}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  N(1) = 1;
  return N;
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
double muu(int tag)
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
Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(2));

  return v;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
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

Vector v_exact(Vector const& , double , int ) //(X,t,tag)
{

}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  N(1) = 1;
  return N;
}


#endif

#if (PROBLEM_TYPE==COUETTE && true)

double const w_ = 1;
double const a__= 0.4;

double pho(Vector const& X, int tag)
{
  return 1.0;
}
double cos_theta0()
{
  return 0.0;
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
double muu(int tag)
{
  return  1.;
}
Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
  Tensor dxU(grad_u_exact(X,t,tag));

  //f = dxU*u_exact(X,t,tag);
  ////f(0) += -a__*y*w_*sin(w_*t);
  ////f(1) += +a__*x*w_*cos(w_*t);
  //f(1) += +a__*x;
  f(0) = t;


  
  return f;
}
Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector v(Vector::Zero(X.size()));

  //v(0) = a__*y*cos(w_*t);
  //v(1) = a__*x*sin(w_*t);
  //v(0) = a__*y*cos(w_*t);
  //v(1) = a__*x*cos(w_*t);
  //v(0) = a__*y*sin(w_*t);
  //v(1) = a__*x*sin(w_*t);
  //v(0) = 0;
  //v(1) = a__*x*t;
  v(0) = 1.-y*y;
  v(1) = 0;
  //v(0) = y*t;
  //v(1) = 0;


  
  
  return v;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  //double x = X(0);
  //double y = X(1);

  //return 0;
  return  -2*X(0) + X(0)*t;
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(2);
  dxP(0) = 0;
  dxP(1) = 0;

  return dxP;
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(), X.size()));

  dxU(0,0) = 0;
  dxU(0,1) = -2*y;//*cos(w_*t);
  dxU(1,0) = 0;//sin(w_*t);
  dxU(1,1) = 0;
  
  //dxU(0,0) = 0;
  //dxU(0,1) = 0;
  //dxU(1,0) = a__*t;
  //dxU(1,1) = 0;
  
  return dxU;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));

  //T(0) = -pressure_exact(X,t,tag);
  //T(1) = muu(tag)*(cos(w_*t) + sin(w_*t));

  Tensor dxU(grad_u_exact(X,t,tag));
  Tensor I(Tensor::Identity(2,2));

  T = (- pressure_exact(X,t,tag)*I +  muu(tag)*(dxU + dxU.transpose()))*normal;

  return T;
}
Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}
double p_initial(Vector const& X, int tag)
{
  return pressure_exact(X,0,tag);
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));

  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const tet = 10. * (pi/180.);
  double const x = X(0);
  double const y = X(1);
  Vector v(Vector::Zero(X.size()));
  //v(0) = +.1*( y );//*(X(0)+0.5)*(X(0)-0.5)/2.;
  //v(1) = +.1*( 0 );//*(X(1)+0.5)*(X(1)-0.5)/2.;
  //v(0) = +0.1*(  1 );//*(X(0)+0.5)*(X(0)-0.5)/2.;
  //v(1) = +0.1*(  0 );//*(X(1)+0.5)*(X(1)-0.5)/2.;
  //return v * (tag!=1 && tag!=2);
  v(0) = t*(1-x*x)*(1+y)/32.;
  v(1) = t*(1-y*y)*(x + t*(1-x*x)/32. + 1)/32.;
  return v;
  //return u_exact(X,t,tag);
}

// posição do contorno
Vector x_exact(Vector const& X, double t, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  N(1) = 1;
  return N;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  return Tensor::Zero(2,2);
}

#endif

#if (false)

double const w_ = 1;
double const a__= 0.4;

double pho(Vector const& X, int tag)
{
  return 1.0;
}
double cos_theta0()
{
  return 0.0;
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
double muu(int tag)
{
  return  1.;
}
Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));
  Tensor dxU(grad_u_exact(X,t,tag));

  f = dxU*u_exact(X,t,tag);
  f(0) += -a__*y*w_*sin(w_*t);
  f(1) += +a__*x*w_*cos(w_*t);
  //f(1) += +a__*x;
  
  return f;
}
Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector v(Vector::Zero(X.size()));

  v(0) = a__*y*cos(w_*t);
  v(1) = a__*x*sin(w_*t);
  //v(0) = a__*y*cos(w_*t);
  //v(1) = a__*x*cos(w_*t);
  //v(0) = a__*y*sin(w_*t);
  //v(1) = a__*x*sin(w_*t);
  //v(0) = 0;
  //v(1) = a__*x*t;
  //v(0) = 0;
  //v(1) = a__*x*t;
  //v(0) = 1.-y*y;
  //v(1) = 0;
  //v(0) = y*t;
  //v(1) = 0;
  
  
  return v;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  //double x = X(0);
  //double y = X(1);

  return 0;
  //return  -2*X(0);
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector dxP(2);
  dxP(0) = 0;
  dxP(1) = 0;

  return dxP;
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(), X.size()));

  //dxU(0,0) = 0;
  //dxU(0,1) = -2*y;//*cos(w_*t);
  //dxU(1,0) = 0;//sin(w_*t);
  //dxU(1,1) = 0;
  
  dxU(0,0) = 0;
  dxU(0,1) = a__*cos(w_*t);
  dxU(1,0) = a__*sin(w_*t);
  dxU(1,1) = 0;
  
  return dxU;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));

  //T(0) = -pressure_exact(X,t,tag);
  //T(1) = muu(tag)*(cos(w_*t) + sin(w_*t));

  Tensor dxU(grad_u_exact(X,t,tag));
  Tensor I(Tensor::Identity(2,2));

  T = (- pressure_exact(X,t,tag)*I +  muu(tag)*(dxU + dxU.transpose()))*normal;

  return T;
}
Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}
double p_initial(Vector const& X, int tag)
{
  return pressure_exact(X,0,tag);
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));

  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const tet = 10. * (pi/180.);
  double const x = X(0);
  double const y = X(1);
  Vector v(Vector::Zero(X.size()));
  v(0) = +.1*( -cos(t) );//*(X(0)+0.5)*(X(0)-0.5)/2.;
  v(1) = +.1*( +cos(t) );//*(X(1)+0.5)*(X(1)-0.5)/2.;
  //v(0) = +0.1*(  1 );//*(X(0)+0.5)*(X(0)-0.5)/2.;
  //v(1) = +0.1*(  0 );//*(X(1)+0.5)*(X(1)-0.5)/2.;
  //return v * (tag!=1 && tag!=2);
  //v(0) = t*(1-x*x)*(1+y)/32.;
  //v(1) = t*(1-y*y)*(x + t*(1-x*x)/32. + 1)/32.;
  return v;
  //return u_exact(X,t,tag);
}

// posição do contorno
Vector x_exact(Vector const& X, double t, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  N(1) = 1;
  return N;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  return Tensor::Zero(2,2);
}

#endif


#if (PROBLEM_TYPE==PSEUDO_OSC2D)

double const w_ = pi;
double const a__= 1;

double pho(Vector const& X, int tag)
{
  return 1.0;
}
double cos_theta0()
{
  return 0.0;
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
double muu(int tag)
{
  return  1.;
}

Vector force(Vector const& X, double t, int tag) // original
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));

  // linear case
  f(0) = a__* (  -w_*sin(w_*t)*x + a__*cos(w_*t)*cos(w_*t)*x );
  f(1) = a__* (  +w_*sin(w_*t)*y + a__*cos(w_*t)*cos(w_*t)*y );

  //// quadratic case
  //Tensor const dxU( grad_u_exact(X,t,tag) );
  //Vector const u(u_exact(X,t,tag));
  //f(0) =-a__*w_*sin(w_*t)* ( x*x/2. + y*x   ) ;
  //f(1) =+a__*w_*sin(w_*t)* ( y*y/2. + x*y   ) ;
  //f += dxU*u;

  return f;
}
Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector v(Vector::Zero(X.size()));

  // linear case
  v(0) = a__*x*cos(w_*t);
  v(1) =-a__*y*cos(w_*t);
  
  ////quadratic case
  //v(0) = a__*cos(w_*t)* ( x*x/2. + y*x   );
  //v(1) =-a__*cos(w_*t)* ( y*y/2. + x*y   );
  
  return v;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  //// linear case
  return 0;
  
  //// quadratic
  //return  a__*cos(w_*t)*(x - y);
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  //double x = X(0);
  //double y = X(1);
  Vector dxP(Vector::Zero(X.size()));
  
  // linear
  // nothing
  
  //// quadratic
  //dxP(0) = a__*cos(w_*t);
  //dxP(1) = a__*cos(w_*t);

  return dxP;
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(), X.size()));

  // linear case
  dxU(0,0) = a__*cos(w_*t);
  dxU(0,1) = 0;
  dxU(1,0) = 0;
  dxU(1,1) =-a__*cos(w_*t);
  
  //// quadratic case
  //dxU(0,0) = a__*cos(w_*t)*(x+y);
  //dxU(0,1) = a__*cos(w_*t)*x;
  //dxU(1,0) = a__*cos(w_*t)*(-y);
  //dxU(1,1) = a__*cos(w_*t)*(-x-y);
  
  return dxU;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));

  Tensor dxU(X.size(),X.size());
  Tensor I(Tensor::Identity(X.size(),X.size()));
  dxU = grad_u_exact(X,t,tag);
  double const p = pressure_exact(X,t,tag);

  T = (-p*I + dxU+dxU.transpose())*normal;

  return T;
}
Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}
double p_initial(Vector const& X, int tag)
{
  return 0;
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));

  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const w_ = pi;
  Vector v(Vector::Zero(X.size()));
  v(0) = -0.0*(cos(w_*t)) - 0.0*sin(w_*t);//*(X(0)+0.5)*(X(0)-0.5)/2.;
  v(1) = +0.0*(cos(w_*t)) + 0.0*sin(w_*t);//*(X(1)+0.5)*(X(1)-0.5)/2.;
  return v * (tag!=1 && tag!=2);
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  N(1) = 1;
  return N;
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
double muu(int tag)
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
Vector u_exact(Vector const& X, double t, int tag)
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
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
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

  return x*(1-x);
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

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));

  return N;
}

Vector v_exact(Vector const& X, double , int ) //(X,t,tag)
{
  return Vector::Zero(X.size());
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  N(1) = 1;
  return N;
}


#endif

#if (PROBLEM_TYPE==GERBEAU)
double pho(Vector const& X, int tag)
{
  return 0.81;
}

double cos_theta0()
{
  return 0;
}

double zeta(double u_norm, double angle)
{
  return 0;
}

double beta_diss()
{
  return 1.5;
}

double gama(Vector const& X, double t, int tag)
{
  return 5.5;
}

double muu(int tag)
{
  return 1.95;
}

Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector f(Vector::Zero(X.size()));
  
  return f;
}
Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(X.size()));
  double h = 13.6;
  double V = 0.2;

  v(1) = 0;
  v(0) = (y/h-1./2)*2*V;

  return v;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector T(Vector::Zero(X.size()));

  //T = -X;

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

  if (X(1)<0.5)
    N(1) = 1;
  else
    N(1) = -1;

  return N;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  if (X(1)<0.5)
    N(0) = -0.25;
  else
    N(0) = 0.25;
  return N;
}


Vector u_initial(Vector const& X, int /*tag*/)
{
  Vector r(Vector::Zero(X.size()));
  double h = 13.6;
  double L = 27.2;
  double V0 = 0.21;

  r(0) = (2*X(1)/h - 1)*V0;

  return r;
}
double p_initial(Vector const& X, int tag)
{
  return X.size()-1;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  Vector v(Vector::Zero(X.size()));
  
  if (tag != 4 && tag != 5)
    v(0) = sin(20*t);
  
  return v;
}


#endif

#if (PROBLEM_TYPE==ZETA3D)
double pho(Vector const& X, int tag)
{
  return 0;
}

double cos_theta0()
{
  return 0*sqrt(2.)/2.;
}

double zeta(double u_norm, double angle)
{
  return 0*1.e-4;
}

double beta_diss()
{
  return 0*1.e-4;
}

double gama(Vector const& X, double t, int tag)
{
  return 0*0.075;
}
double muu(int tag)
{
  return 1.e-0;
}
Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector f(Vector::Zero(X.size()));
  
  //f(1)=-4;
  f(0)=1;
  //f /= 4*sqrt(2);


  return f;
}
Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(X.size()));

  v(0) = 0;

  return v;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector T(Vector::Zero(X.size()));

  //T = -X;

  return T;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  //return 2*gama(X,t,tag);
  return x;
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
  
  double x = X(0);
  double y = X(1);
  double z = X(2);
  
  if (tag == 2) { // contact line
    if (x < 1.e-9)
      N(0) = 1;
    else if (y < 1.e-9)
      N(1) = 1;
    else if (z < 1.e-9)
      N(2) = 1;
    return N;
  }
  
  if (tag == 4)
    N(1) = 1;
  else if (tag == 6 )
    N(0) = 1;
  else
    N(2) = 1;
  
  return N;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}


Vector u_initial(Vector const& X, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return u_exact(X,0,tag);
}
double p_initial(Vector const& X, int tag)
{
  return pressure_exact(X,0,tag);
}

Vector v_exact(Vector const& X, double , int ) //(X,t,tag)
{
  Vector v(Vector::Zero(X.size()));
  return v;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  
  Vector A = solid_normal(X,t,tag);
  
  if (tag == 7)
  {
    f(0,0) = 1;
  }
  else if (tag == 8)
  {
    f(1,1) = 1;
  }
  else
  {
    f(2,2) = 1;
  }
  
  return f;
}
#endif

#if (PROBLEM_TYPE==MICRO2D)
double pho(Vector const& X, int tag)
{
  return 0.81;
}

double cos_theta0()
{
  return 0;
}

double zeta(double u_norm, double angle)
{
  return 0;
}

double beta_diss()
{
  return 1.5;
}

double gama(Vector const& X, double t, int tag)
{
  return 5.5;
}

double muu(int tag)
{
  return 1.95;
}

Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector f(Vector::Zero(X.size()));
  
  return f;
}
Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(X.size()));
  double h = 13.6;
  double V = 0.2;

  v(1) = 0;
  v(0) = (y/h-1./2)*2*V;

  return v;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector T(Vector::Zero(X.size()));

  //T = -X;

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

  if (X(1)<0.5)
    N(1) = 1;
  else
    N(1) = -1;

  return N;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  if (X(1)<0.5)
    N(0) = -0.25;
  else
    N(0) = 0.25;
  return N;
}


Vector u_initial(Vector const& X, int /*tag*/)
{
  Vector r(Vector::Zero(X.size()));
  double h = 13.6;
  double L = 27.2;
  double V0 = 0.21;

  //r(0) = (2*X(1)/h - 1)*V0;

  return r;
}
double p_initial(Vector const& X, int tag)
{
  return X.size()-1;
}

Vector v_exact(Vector const& X, double , int ) //(X,t,tag)
{
  Vector v(X.size());
  v.setZero();
  return v;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  return Tensor::Zero(X.size(), X.size());
}

#endif

#if (PROBLEM_TYPE==TRACO)

double const w_ = pi;
double const a__= 1;

double pho(Vector const& X, int tag)
{
  return 1.0;
}
double cos_theta0()
{
  return 0.0;
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
double muu(int tag)
{
  return  1.;
}

Vector force(Vector const& X, double t, int tag) // original
{
  double x = X(0);
  double y = X(1);

  Vector f(Vector::Zero(X.size()));

  // linear case
  //f(0) = 2.*a__*sin(x)*cos(y);
  //f(1) =-2.*a__*cos(x)*sin(y);

  //// quadratic case
  //f(0) = -2*x*y*y - 2.*x*x*x/3.;
  //f(1) = +2*y*x*x + 2.*y*y*y/3.;

  // se rho = 1
  f += grad_u_exact(X,t,tag)*u_exact(X,t,tag);
  f(0) += x+y;
  f(1) += x-y;
  
  //f += dxU*u;

  return f;
}
Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  Vector v(Vector::Zero(X.size()));

  //// linear case
  //v(0) = a__*sin(x)*cos(y);
  //v(1) =-a__*cos(x)*sin(y);
  
  //////quadratic case
  //v(0) =  x*x*x*y*y/3;
  //v(1) = -y*y*y*x*x/3;
  
  v(0) =  t*(x+y);
  v(1) =  t*(x-y);
  
  return v;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  //// linear case
  return 0;
  
  //// quadratic
  //return  a__*cos(w_*t)*(x - y);
}
Vector grad_p_exact(Vector const& X, double t, int tag)
{
  //double x = X(0);
  //double y = X(1);
  Vector dxP(Vector::Zero(X.size()));
  
  // linear
  // nothing
  
  //// quadratic
  //dxP(0) = a__*cos(w_*t);
  //dxP(1) = a__*cos(w_*t);

  return dxP;
}
Tensor grad_u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Tensor dxU(Tensor::Zero(X.size(), X.size()));

  // linear case
  //dxU(0,0) =  a__*cos(x)*cos(y);
  //dxU(0,1) = -a__*sin(x)*sin(y);
  //dxU(1,0) =  a__*sin(x)*sin(y);
  //dxU(1,1) = -a__*cos(x)*cos(y);;
  
  //dxU(0,0) = x*x*y*y;
  //dxU(0,1) = 2.*x*x*x*y/3;
  //dxU(1,0) = -2.*y*y*y*x/3;
  //dxU(1,1) = -y*y*x*x;
  
  dxU(0,0) = t*1;
  dxU(0,1) = t*1;
  dxU(1,0) = t*1;
  dxU(1,1) = -t*1;  
  
  return dxU;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  Vector T(Vector::Zero(X.size()));

  Tensor dxU(X.size(),X.size());
  Tensor I(Tensor::Identity(X.size(),X.size()));
  dxU = grad_u_exact(X,t,tag);
  double const p = pressure_exact(X,t,tag);

  T = (-p*I + dxU+dxU.transpose())*normal;

  return T;
}
Vector u_initial(Vector const& X, int tag)
{
  return u_exact(X,0,tag);
}
double p_initial(Vector const& X, int tag)
{
  return 0;
}

Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));

  return N;
}

Vector v_exact(Vector const& X, double t, int tag) //(X,t,tag)
{
  double const w_ = pi;
  Vector v(Vector::Zero(X.size()));
  //v = u_exact(X,t,tag);
  return v ;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  N(1) = 1;
  return N;
}


Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  return f;
}

#endif

#if (PROBLEM_TYPE==RUSSA_SIN2D)

const double a = 0.1;
const double w = 45.0;
const double R = 10.0;


double pho(Vector const& X, int tag)
{
  return 0*1.e+2;
}

double cos_theta0()
{
  return -0*sqrt(2)/2;
}

double zeta(double u_norm, double angle)
{
  return 1.e-4;
}

double beta_diss()
{
  return 1.e-3;
}

double gama(Vector const& X, double t, int tag)
{
  return 1e+1;
}
double muu(int tag)
{
  return 1.e+2;
}
Vector force(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector f(Vector::Zero(X.size()));
  
  //f(1)=max(-1.0,-t);
  f(1)=-1.e+1;
  //f(1)=-1.e+2;
  //f(0)=-1;
  //f /= 4*sqrt(2);


  return f;
}
Vector u_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector v(Vector::Zero(X.size()));

  return v;
}
Vector traction(Vector const& X, Vector const& normal, double t, int tag)
{
  double x = X(0);
  double y = X(1);
  Vector T(Vector::Zero(X.size()));

  //T = -X;

  return T;
}
double pressure_exact(Vector const& X, double t, int tag)
{
  double x = X(0);
  double y = X(1);

  return 0;
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

/*   Don't forget to normalize the vector  */
Vector solid_normal(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  
  double x = X(0);
  double y = X(1);
  
  double const r = sqrt(x*x + y*y);
  double const r2= r*r;
  
  const double tet = pi + acos(-x/r);
  
  N(0) = -x/r + a*w*sin(w*tet)*y/r2;
  N(1) = -y/r - a*w*sin(w*tet)*x/r2;

  //N(0) = 1;
  //N(1) = 1;

  N.normalize();

  return N;
}

Vector solid_veloc(Vector const& X, double t, int tag)
{
  Vector N(Vector::Zero(X.size()));
  return N;
}


Vector u_initial(Vector const& X, int tag)
{
  Vector r(Vector::Zero(X.size()));
  return r;
}
double p_initial(Vector const& X, int tag)
{
  return 0;
}

Vector v_exact(Vector const& X, double , int ) //(X,t,tag)
{
  Vector v(X.size());
  v.setZero();
  return v;
}

Tensor feature_proj(Vector const& X, double t, int tag)
{
  Tensor f(Tensor::Zero(X.size(), X.size()));
  return f;
}

#endif


