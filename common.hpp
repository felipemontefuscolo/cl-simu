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

double pho(Vector const& X, int tag);
double gama(Vector const& X, double t, int tag);
double cos_theta0();
double zeta(double u_norm, double angle);
double beta_diss();
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
  PetscErrorCode solveTimeProblem();
  PetscErrorCode formJacobian(SNES /*snes*/,Vec x,Mat *Jac, Mat* /*prejac*/, MatStructure * /*flag*/);
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
                          
  //// R size: ndofu x ndofu
  template<class AnyStaticMAtrix, class AnyStaticVector>
  void getRotationMatrix(AnyStaticMAtrix & R, AnyStaticVector/*VectorXi*/ const& nodes, int const n_nodes) const
  {
    //const int n_nodes = nodes.size();
    //const bool mesh_has_edge_nodes  = mesh->numNodesPerCell() > mesh->numVerticesPerCell();
    //bool const u_has_edge_assoc_dof = shape_phi_c->numDofsAssociatedToFacet() +
    //                                  shape_phi_c->numDofsAssociatedToCorner() > 0;
    Vector X(dim);
    Vector normal(dim);
    Vector tmp(dim);
    Point const* point;
    //Cell const* cell;
    //int m;
    int tag=0;
    int nodeid;
    double aux;

    R.setIdentity();

    // NODES
    for (int i = 0; i < n_nodes; ++i)
    {
      nodeid = nodes[i];
      point = mesh->getNode(nodeid);
      tag = point->getTag();
      //m = point->getPosition() - mesh->numVerticesPerCell();
      //cell = mesh->getCell(point->getIncidCell());

      if (!is_in(tag,solid_tags) && !is_in(tag,triple_tags))
        continue;


      point->getCoord(X.data());

      normal = solid_normal(X,current_time,tag);

      R.block(i*dim,i*dim,1,dim)  = normal.transpose();
      tmp.setZero();
      if (dim==2)
      {
        tmp(0) = -normal(1);
        tmp(1) =  normal(0);
        tmp.normalize();
        R.block(i*dim+1,i*dim,1,dim)  = tmp.transpose();
      }
      else
      if (dim==3)
      {
        aux = normal(0)*normal(0) + normal(1)*normal(1);
        if ( aux < 1.e-10)
        {
          aux = sqrt(normal(0)*normal(0) + normal(2)*normal(2));

          tmp(0) = normal(2)/aux;
          tmp(1) = 0;
          tmp(2) = -normal(0)/aux;

          R.block(i*dim+1,i*dim,1,dim)  = tmp.transpose();

          tmp(0) =  -normal(1)*normal(0)/aux;
          tmp(1) =  aux;
          tmp(2) =  -normal(1)*normal(2)/aux;
        }
        else
        {
          aux = sqrt(aux);

          tmp(0) = -normal(1)/aux;
          tmp(1) =  normal(0)/aux;
          tmp(2) = 0;

          R.block(i*dim+1,i*dim,1,dim)  = tmp.transpose();

          tmp(0) =  -normal(2)*normal(0)/aux;
          tmp(1) =  -normal(2)*normal(1)/aux;
          tmp(2) =  aux;
        }

        R.block(i*dim+2,i*dim,1,dim)  = tmp.transpose();
      }




    } // end nodes


  }

  // A <- R^T *A
  template<class AnyStaticMatrix, class OtherStaticMatrix, class TempStaticMatrix>
  void rotate_RtA(AnyStaticMatrix const& R, OtherStaticMatrix & M, TempStaticMatrix & tmp)
  {
    if (R.cols() != M.rows()) // então M é tratado como vetor
    {
      int const n = R.cols();
      tmp.resize(n,1);
      Map<VectorXd> v(M.data(),n);
      tmp = R.transpose()*v;
      v=tmp;
    }
    else
    {
      int const m = R.rows();
      int const n = M.cols();
      tmp.resize(m,n);
      tmp = R.transpose()*M;
      M = tmp;
    }
  }

  // A <- R *A
  template<class AnyStaticMatrix, class OtherStaticMatrix, class TempStaticMatrix>
  void rotate_RA(AnyStaticMatrix const& R, OtherStaticMatrix & M, TempStaticMatrix & tmp)
  {
    if (R.cols() != M.rows()) // então M é tratado como vetor
    {
      int const n = R.cols();
      tmp.resize(n,1);
      Map<VectorXd> v(M.data(),n);
      tmp = R*v;
      v=tmp;
    }
    else
    {
      int const m = R.rows();
      int const n = M.cols();
      tmp.resize(m,n);
      tmp = R*M;
      M = tmp;
    }
  }

  // A <- R *A *R^T
  template<class AnyStaticMatrix, class OtherStaticMatrix, class TempStaticMatrix>
  void rotate_RARt(AnyStaticMatrix const& R, OtherStaticMatrix & M, TempStaticMatrix & tmp)
  {
    int const m = R.rows();
    tmp.resize(m,m);
    tmp = R*M*R.transpose();
    M = tmp;
  }

  template<class AnyStaticMatrix, class OtherStaticMatrix, class TempStaticMatrix>
  void rotate_ARt(AnyStaticMatrix const& R, OtherStaticMatrix & M, TempStaticMatrix & tmp)
  {
    int const m = M.rows();
    int const n = R.cols();
    tmp.resize(m,n);
    tmp = M*R.transpose();
    M = tmp;
  }

  void computeError(Vec &qq, double tt);
  void updateNormals(Vec *x_mesh);
  void smoothsMesh(Vec &x_mesh);
  // @param[in] q unknows vector with fluid velocity
  // @param[out] u_mesh
  PetscErrorCode calcMeshVelocity(Vec const& q);
  PetscErrorCode moveMesh();
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
  double      dt;
  double      steady_tol;
  double      theta;
  int         maxts;
  bool        force_pressure;
  PetscBool   print_to_matlab;
  PetscBool   force_dirichlet;
  PetscBool   full_diriclet;
  PetscBool   ale;
  bool        solve_the_sys;
  PetscBool   plot_exact_sol;
  int         quadr_degree_cell;
  int         quadr_degree_facet;
  int         quadr_degree_corner;
  int         quadr_degree_err;
  bool        pres_pres_block;
  float       grow_factor;
  string      filename;
  string      filename_out;
  PetscBool  family_files;

  std::vector<int> dirichlet_tags;   // vetor com os tags que indicam contorno dirichlet
  std::vector<int> neumann_tags;     // vetor com os tags que indicam contorno neumann
  std::vector<int> interface_tags;
  std::vector<int> solid_tags;
  std::vector<int> triple_tags;

  DofHandler                   dof_handler_vars;
  DofHandler                   dof_handler_mesh;
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
  std::vector<int>             dir_vertices;  // os vertices dirichlet
  std::vector<int>             dir_corners;   // os vertices dirichlet
  std::vector<int>             dir_facets;
  std::vector<int>             dir_normal_vertices;  // os vertices dirichlet
  std::vector<int>             dir_normal_corners;   // os vertices dirichlet
  std::vector<int>             dir_normal_facets;
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

  // petsc vectors
  Vec                 q0, q, res, u_mesh, x_mesh, nml_mesh/**/;
  Mat                 Jac;
  SNES                snes;         /* nonlinear solver context */
  KSP    			        ksp;
  PC	   			        pc;

};


























