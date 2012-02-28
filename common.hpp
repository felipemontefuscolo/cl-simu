#define EIGEN_NO_AUTOMATIC_RESIZING
#include <Fepic/Mesh>
#include <Fepic/Quadrature>
#include <Fepic/DofHandler>
#include <Fepic/Shape>
#include <Eigen/LU>
#include "petscsnes.h"
#include <iostream>
#include <tr1/memory>
#include <tr1/array>
#include "mypetsc.hpp"


//#ifdef EIGEN_HAS_OPENMP
//#undef EIGEN_HAS_OPENMP
//#endif

using namespace std;
using namespace Eigen;
using namespace tr1;

// por célula
#define MAX_DOFS_U 33 // P2hp
#define MAX_DOFS_P 4  // P1
#define MAX_NNODES 10 // P2

enum VarNumber {
  // DH_UNKS
  VAR_U = 0,
  VAR_P = 1,
  VAR_L = 2,
  
  // DH_MESH
  VAR_M = 0
};

enum DofHandlerNumber {
  DH_UNKS = 0,
  DH_MESH = 1
};

enum PairSpace {

  P1P1      = 1,
  P1bP1_c   = 2,
  P2bPm1_c  = 3,
  P2P1      = 4,
  P1bP1     = 5,
  P2P0      = 6,
  P2bPm1    = 7
  
};

const double pi  = 3.141592653589793;
const double pi2 = pi*pi;

typedef Matrix<VectorXd, Dynamic,1>  VecOfVec;
typedef Matrix<MatrixXd, Dynamic,1>  VecOfMat;

// space
typedef Matrix<double, Dynamic,1,0,3,1>              Vector;
typedef Matrix<double, Dynamic,Dynamic,RowMajor,3,3> Tensor;


template<class Vec, class T>
bool is_in(T value, Vec const& v)
{
  for (int i = 0; i < (int)v.size(); i++)
  {
    if (value==v[i])
      return true;
  }
  return false;
}

/* stabilization type */
enum Behaviors {

  BH_bble_condens_PnPn = 0x01,
  BH_GLS               = 0x02,
  BH_Press_grad_elim   = 0x04,
  BH_bble_condens_CR   = 0x08
};

class AppCtx;
class Statistics;

double pho(Vector const& X, int tag);
double gama(Vector const& X, double t, int tag);
double cos_theta0();
double zeta(double u_norm, double angle);
double beta_diss();
double muu(int tag);
Vector force(Vector const& X, double t, int tag);
Vector u_exact(Vector const& X, double t, int tag);
Vector traction(Vector const& X, double t, int tag);
double pressure_exact(Vector const& X, double t, int tag);
Vector grad_p_exact(Vector const& X, double t, int tag);
Tensor grad_u_exact(Vector const& X, double t, int tag);
Vector u_initial(Vector const& X, int tag);
double p_initial(Vector const& X, int tag);
Vector solid_normal(Vector const& X, double t, int tag);
Vector v_exact(Vector const& X, double t, int tag);


inline double sqr(double v) {return v*v;}

void inline inverseAndDet(Tensor const& a, int dim, Tensor& inv, double& det)
{
  if (dim==2)
  {
    det = a(0,0)*a(1,1)-a(0,1)*a(1,0);

    inv(0,0) = a(1,1)/det;
    inv(0,1) = -a(0,1)/det;
    inv(1,0) = -a(1,0)/det;
    inv(1,1) = a(0,0)/det;
  }
  else if (dim==3)
  {
    det = a(0,0)*(a(1,1)*a(2,2)-a(1,2)*a(2,1))+a(0,1)*(a(1,2)*a(2,0)-a(1,0)*a(2,2))+a(0,2)*(a(1,0)*a(2,1)-a(1,1)*a(2,0));

    inv(0,0) = ( a(1,1)*a(2,2)-a(1,2)*a(2,1) )/det;
    inv(0,1) = ( a(0,2)*a(2,1)-a(0,1)*a(2,2) )/det;
    inv(0,2) = ( a(0,1)*a(1,2)-a(0,2)*a(1,1) )/det;
    inv(1,0) = ( a(1,2)*a(2,0)-a(1,0)*a(2,2) )/det;
    inv(1,1) = ( a(0,0)*a(2,2)-a(0,2)*a(2,0) )/det;
    inv(1,2) = ( a(0,2)*a(1,0)-a(0,0)*a(1,2) )/det;
    inv(2,0) = ( a(1,0)*a(2,1)-a(1,1)*a(2,0) )/det;
    inv(2,1) = ( a(0,1)*a(2,0)-a(0,0)*a(2,1) )/det;
    inv(2,2) = ( a(0,0)*a(1,1)-a(0,1)*a(1,0) )/det;

  }
  else
  {
    printf("invalid dim\n");
    throw;
  }


}

double inline determinant(Tensor const& a, int dim)
{
  if (dim==1)
    return a(0,0);
  else
  if (dim==2)
    return a(0,0)*a(1,1)-a(0,1)*a(1,0);
  else
  if (dim==3)
    return a(0,0)*(a(1,1)*a(2,2)-a(1,2)*a(2,1))+a(0,1)*(a(1,2)*a(2,0)-a(1,0)*a(2,2))+a(0,2)*(a(1,0)*a(2,1)-a(1,1)*a(2,0));
  else
  {
    printf("double determinant(Tensor const& a, int dim): invalid dim, get %d\n", dim);
    throw;
  }
}

template<class TensorType>
void invert(TensorType & a, int dim)
{
  if (dim==1)
  {
    a(0,0)=1./a(0,0);
  }
  else
  if (dim==2)
  {
    double const det = a(0,0)*a(1,1)-a(0,1)*a(1,0);

    double const inv00 = a(1,1)/det;
    double const inv01 = -a(0,1)/det;
    double const inv10 = -a(1,0)/det;
    double const inv11 = a(0,0)/det;

    a(0,0) = inv00;
    a(0,1) = inv01;
    a(1,0) = inv10;
    a(1,1) = inv11;
  }
  else if (dim==3)
  {
    double const det = a(0,0)*(a(1,1)*a(2,2)-a(1,2)*a(2,1))+a(0,1)*(a(1,2)*a(2,0)-a(1,0)*a(2,2))+a(0,2)*(a(1,0)*a(2,1)-a(1,1)*a(2,0));

    double const inv00 = ( a(1,1)*a(2,2)-a(1,2)*a(2,1) )/det;
    double const inv01 = ( a(0,2)*a(2,1)-a(0,1)*a(2,2) )/det;
    double const inv02 = ( a(0,1)*a(1,2)-a(0,2)*a(1,1) )/det;
    double const inv10 = ( a(1,2)*a(2,0)-a(1,0)*a(2,2) )/det;
    double const inv11 = ( a(0,0)*a(2,2)-a(0,2)*a(2,0) )/det;
    double const inv12 = ( a(0,2)*a(1,0)-a(0,0)*a(1,2) )/det;
    double const inv20 = ( a(1,0)*a(2,1)-a(1,1)*a(2,0) )/det;
    double const inv21 = ( a(0,1)*a(2,0)-a(0,0)*a(2,1) )/det;
    double const inv22 = ( a(0,0)*a(1,1)-a(0,1)*a(1,0) )/det;

    a(0,0) = inv00;
    a(0,1) = inv01;
    a(0,2) = inv02;
    a(1,0) = inv10;
    a(1,1) = inv11;
    a(1,2) = inv12;
    a(2,0) = inv20;
    a(2,1) = inv21;
    a(2,2) = inv22;

  }
  else
  {
    printf("invalid dim, try to run again dumb \n");
    throw;
  }


}


template<class AnyVector>
void cross(AnyVector & a, AnyVector const& b)
{
  double const r0 = a(1)*b(2) - a(2)*b(1);
  double const r1 = a(2)*b(0) - a(0)*b(2);
  double const r2 = a(0)*b(1) - a(1)*b(0);

  a(0) = r0;
  a(1) = r1;
  a(2) = r2;

}




class Statistics
{
public:
  Statistics()
  {
    this->reset();
  }
  
  void reset()
  {
    p_L2_norm=0;
    u_L2_norm=0;
    grad_u_L2_norm=0;
    grad_p_L2_norm=0;

    Np=0;
    Nu=0;
    Ngp=0;
    Ngu=0;    
  }
  
  void add_p_L2(double x)
  {
    p_L2_norm += x;
    //Np++;
  }
  void add_u_L2(double x)
  {
    u_L2_norm += x;
    //Nu ++;
  }
  void add_grad_u_L2(double x)
  {
    grad_u_L2_norm +=x;
    //Ngu++;
  }
  void add_grad_p_L2(double x)
  {
    grad_p_L2_norm += x;
    //Ngp++;
  }
  
  double mean_p_L2_norm()
  {
    return p_L2_norm;
  }
  
  double mean_u_L2_norm()
  {
    return u_L2_norm;
  }
  
  double mean_grad_u_L2_norm()
  {
    return grad_u_L2_norm;
  }
  
  double mean_grad_p_L2_norm()
  {
    return grad_p_L2_norm;
  }  
  
  double p_L2_norm;
  double u_L2_norm;
  double grad_u_L2_norm;
  double grad_p_L2_norm;
  
  
  int Np, Nu, Ngp, Ngu;
  
};



//
// User Class
//
class AppCtx
{
public:
  AppCtx(int argc, char **argv, bool &help_return, bool &erro);
  
  bool err_checks();
  void loadMesh();

  void loadDofs();
  void setUpDefaultOptions();
  bool getCommandLineOptions(int argc, char **/*argv*/);
  bool createFunctionsSpace();
  void createQuadrature();
  void meshAliasesUpdate();
  void dofsCreate();
  void dofsUpdate();
  PetscErrorCode allocPetscObjs();
  void matrixColoring();
  void printMatlabLoader();
  // must be called after loadDofs
  void evaluateQuadraturePts();
  void computeDirichletEntries();
  void onUpdateMesh();
  void setInitialConditions();
  PetscErrorCode checkSnesConvergence(SNES snes, PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason);
  PetscErrorCode solveTimeProblem();
  PetscErrorCode formJacobian(SNES /*snes*/,Vec x,Mat *Mat_Jac, Mat* /*prejac*/, MatStructure * /*flag*/);
  PetscErrorCode formFunction(SNES /*snes*/, Vec x, Vec f);
  // form the residue of the cell
  void formCellFunction(cell_iterator &cell,
                                  VectorXi &mapU_c,  VectorXi &/*mapP_c*/, // mappers
                                  MatrixXd &u_coefs_c_new,  VectorXd &p_coefs_c, // coefficients
                                  VectorXd &FUloc, VectorXd &FPloc); // output: local residue
  // form the residue of the facet
  void formFacetFunction(facet_iterator &facet,
                         VectorXi const&/*mapU_f*/,  VectorXi const&/*mapP_f*/, // mappers
                         MatrixXd &u_coefs_f,  VectorXd &/*p_coefs_f*/, // coefficients
                         VectorXd &FUloc); // output: local residue
  // form the residue of the contact line
  void formCornerFunction(corner_iterator &corner,
                          VectorXi const&/*mapU_r*/,  VectorXi const&/*mapP_r*/, // mappers
                          MatrixXd &u_coefs_r, // coefficients
                          VectorXd &FUloc);
                          
  
  // get points dofs even if he lie at an edge
  /// @param[out] dofs
  void getNodeDofs(Point const* point, DofHandlerNumber DH_, VarNumber VAR_, int * dofs) const
  {
    if (mesh->isVertex(point))
    {
      dof_handler[DH_].getVariable(VAR_).getVertexAssociatedDofs(dofs, point);
    }
    else
    {
      const int m = point->getPosition() - mesh->numVerticesPerCell();
      Cell const* cell = mesh->getCell(point->getIncidCell());
      if (dim==3)
      {
        const int edge_id = cell->getCornerId(m);
        dof_handler[DH_].getVariable(VAR_).getCornerAssociatedDofs(dofs, mesh->getCorner(edge_id));
      }
      else
      {
        const int edge_id = cell->getFacetId(m);
        dof_handler[DH_].getVariable(VAR_).getFacetAssociatedDofs(dofs, mesh->getFacet(edge_id));
      }
    }
  }

  double getMeshVolume();
  double getMaxVelocity();
  void printContactAngle(bool _print);
  
  
  void computeError(Vec const& Vec_x, Vec &Vec_up_1, double tt);
  void getVecNormals(Vec const* Vec_x_1, Vec & Vec_normal_);
  void smoothsMesh(Vec & Vec_normal_, Vec &Vec_x);
  void copyMesh2Vec(Vec &Vec_xmsh);
  void copyVec2Mesh(Vec const& Vec_xmsh);
  void swapMeshWithVec(Vec & Vec_xmsh);
  // @param[in] Vec_up_1 unknows vector with fluid velocity
  // @param[out] u_mesh
  PetscErrorCode calcMeshVelocity(Vec const& Vec_x_0, Vec const& Vec_x_1, Vec &Vec_v_mid);
  PetscErrorCode moveMesh(Vec const& Vec_x_0, Vec const& Vec_up_0, Vec const& Vec_up_1, double const vtheta, double tt, Vec & Vec_x_new);
  void freePetscObjs();
  
  
  // global settings
  int         dim; // = space dim
  int         mesh_cell_type;   // ECellType
  int         function_space;  // P1P1, P2P1, P2P0, etc ...
  int         behaviors;
  //double      Re;
  PetscBool   has_convec;
  PetscBool   unsteady;
  PetscBool   renumber_dofs;
  PetscBool   boundary_smoothing;
  PetscBool   print_to_matlab;
  PetscBool   force_dirichlet;
  PetscBool   full_diriclet;
  PetscBool   force_mesh_velocity;
  PetscBool   ale;
  PetscBool   plot_exact_sol;
  PetscBool   family_files;
  PetscBool   fprint_ca; // print contact angle
  int         converged_times;
  double      dt;
  double      steady_tol;
  double      utheta;
  double      vtheta;
  int         maxts;
  bool        force_pressure;
  bool        solve_the_sys;
  int         quadr_degree_cell;
  int         quadr_degree_facet;
  int         quadr_degree_corner;
  int         quadr_degree_err;
  bool        pres_pres_block;
  float       grow_factor;
  string      filename;
  string      filename_out;

  std::vector<int> dirichlet_tags;   // vetor com os tags que indicam contorno dirichlet
  std::vector<int> neumann_tags;     // vetor com os tags que indicam contorno neumann
  std::vector<int> interface_tags;
  std::vector<int> solid_tags;
  std::vector<int> triple_tags;

  DofHandler                   dof_handler[2];
  MeshIoMsh                    msh_reader;
  MeshIoVtk                    vtk_printer;
  shared_ptr<Mesh>             mesh;
  // velocity
  shared_ptr<ShapeFunction>    shape_phi_c; // cell
  shared_ptr<ShapeFunction>    shape_phi_f; // facet
  shared_ptr<ShapeFunction>    shape_phi_r; // corner
  // pressure
  shared_ptr<ShapeFunction>    shape_psi_c; // cell
  shared_ptr<ShapeFunction>    shape_psi_f; // facet
  shared_ptr<ShapeFunction>    shape_psi_r; // corner
  // mesh
  shared_ptr<ShapeFunction>    shape_qsi_c; // cell
  shared_ptr<ShapeFunction>    shape_qsi_f; // facet
  shared_ptr<ShapeFunction>    shape_qsi_r; // corner

  shared_ptr<ShapeFunction>    shape_bble;

  shared_ptr<Quadrature>       quadr_cell;
  shared_ptr<Quadrature>       quadr_facet;
  shared_ptr<Quadrature>       quadr_corner;
  shared_ptr<Quadrature>       quadr_err;   // to compute error
  int                          n_qpts_cell;
  int                          n_qpts_facet;
  int                          n_qpts_corner;
  int                          n_qpts_err;
  std::vector<int>             dir_entries;   // as linhas da matriz onde são forçadas as c.c.

  // velocity
  VecOfVec                     phi_c;         // shape function evaluated at quadrature points
  VecOfVec                     phi_f;         // shape function evaluated at quadrature points (facet)
  VecOfVec                     phi_r;         // shape function evaluated at quadrature points (corner)
  // pressure
  VecOfVec                     psi_c;         // shape function evaluated at quadrature points
  VecOfVec                     psi_f;         // shape function evaluated at quadrature points (facet)
  VecOfVec                     psi_r;         // shape function evaluated at quadrature points (corner)
  // mesh
  VecOfVec                     qsi_c;         // shape function evaluated at quadrature points
  VecOfVec                     qsi_f;         // shape function evaluated at quadrature points (facet)
  VecOfVec                     qsi_r;         // shape function evaluated at quadrature points (corner)
  VectorXd                     qsi_c_at_center; // shape function evaluated at quadrature points
  // velocity
  VecOfMat                     dLphi_c;       // matriz de gradiente no elemento unitário
  VecOfMat                     dLphi_f;       // matriz de gradiente no elemento unitário (facet)
  VecOfMat                     dLphi_r;       // matriz de gradiente no elemento unitário (corner)
  // pressure
  VecOfMat                     dLpsi_c;       // matriz de gradiente no elemento unitário
  VecOfMat                     dLpsi_f;       // matriz de gradiente no elemento unitário (facet)
  VecOfMat                     dLpsi_r;       // matriz de gradiente no elemento unitário (corner)
  // mesh
  VecOfMat                     dLqsi_c;       // matriz de gradiente no elemento unitário
  VecOfMat                     dLqsi_f;       // matriz de gradiente no elemento unitário (facet)
  VecOfMat                     dLqsi_r;       // matriz de gradiente no elemento unitário (corner)
  // normal
  VecOfMat                     dLphi_nf;       // matriz de gradiente nos nós de uma facet
  // velocity, cell
  VectorXd                     bble;          // bubble function evaluated at quadrature points
  VecOfVec                     dLbble;

  //
  //          to compute error in each cell
  // velocity
  VecOfVec                     phi_err;         // shape function evaluated at quadrature points
  VecOfMat                     dLphi_err;       // matriz de gradiente no elemento unitário
  // pressure
  VecOfVec                     psi_err;         // shape function evaluated at quadrature points
  VecOfMat                     dLpsi_err;       // matriz de gradiente no elemento unitário
  // mesh
  VecOfVec                     qsi_err;         // shape function evaluated at quadrature points
  VecOfMat                     dLqsi_err;       // matriz de gradiente no elemento unitário


  // dofs
  int           n_unknowns;
  int           n_dofs_u_mesh;
  int           n_dofs_u_per_cell;
  int           n_dofs_u_per_facet;
  int           n_dofs_u_per_corner;
  int           n_dofs_p_per_cell;
  int           n_dofs_p_per_facet;
  int           n_dofs_p_per_corner;

  // mesh alias
  int           n_nodes;
  int           n_cells;
  int           n_facets;
  int           n_corners;
  int           n_nodes_total;
  int           n_cells_total;
  int           n_facets_total;
  int           n_corners_total;
  int           nodes_per_cell;
  int           nodes_per_facet;
  int           nodes_per_corner;

  // solving ...
  double        current_time;
  int           time_step;
  int           print_step;
  double        beta1, beta2;

  Timer         timer;
  Statistics    Stats;

  // petsc vectors
  Vec                 Vec_res, Vec_up_0, Vec_up_1, Vec_v_mid, Vec_x_0, Vec_x_1, Vec_normal/**/;
  Mat                 Mat_Jac;
  SNES                snes;         /* nonlinear solver context */
  KSP    			        ksp;
  PC	   			        pc;

};






















