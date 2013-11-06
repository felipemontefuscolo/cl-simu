//static char help[] = "Navier-Stokes.\n\n";

//
// ALE
//
// pho*( (Ut + (U-Umsh) · nabla)U ) + grad p = muu* div grad U + force
//
// LIMITAÇÕES:
// * triangulo e tetrahedro apenas.
// * no maximo 1 grau de liberdade associado a uma aresta
// * condição de dirichlet na normal só pode valer 0 .. procure por UNORMAL_AQUI
// * elementos isoparamétricos
//
// OBS:
// * a normal no ponto de contato é igual ao limite da normal pela superfície livre
//
//
//


#include "common.hpp"
#include <iomanip>

PetscErrorCode FormJacobian(SNES snes,Vec Vec_up_1,Mat *Mat_Jac, Mat *prejac, MatStructure *flag, void *ptr);
PetscErrorCode FormFunction(SNES snes, Vec Vec_up_1, Vec Vec_fun, void *ptr);
PetscErrorCode CheckSnesConvergence(SNES snes, PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx);

PetscErrorCode FormJacobian_mesh(SNES snes,Vec Vec_up_1,Mat *Mat_Jac, Mat *prejac, MatStructure *flag, void *ptr);
PetscErrorCode FormFunction_mesh(SNES snes, Vec Vec_up_1, Vec Vec_fun, void *ptr);

class AppCtx;
class Statistics;

class GetDataVelocity : public DefaultGetDataVtk
{
public:
  GetDataVelocity(double *q_array_, AppCtx const& user_) : user(user_), q_array(q_array_){}
  //double get_data_r(int nodeid) const;
  void get_vec(int id, Real * vec_out) const;
  int vec_ncomps() const { return user.mesh->spaceDim(); }
  AppCtx const& user;
  double *q_array;
  virtual ~GetDataVelocity() {}
};

class GetDataPressure : public DefaultGetDataVtk
{
public:
  GetDataPressure(double *q_array_, AppCtx const& user_) : user(user_),q_array(q_array_){}
  double get_data_r(int nodeid) const;
  AppCtx const& user;
  double *q_array;
  virtual ~GetDataPressure() {}
};

class GetDataPressCellVersion : public DefaultGetDataVtk
{
public:
  GetDataPressCellVersion(double *q_array_, AppCtx const& user_) : user(user_),q_array(q_array_){}
  double get_data_r(int cellid) const;
  AppCtx const& user;
  double *q_array;
  virtual ~GetDataPressCellVersion() {}
};


class GetDataNormal : public DefaultGetDataVtk
{
public:
  GetDataNormal(double *q_array_, AppCtx const& user_) : user(user_), q_array(q_array_){}
  void get_vec(int id, Real * vec_out) const;
  int vec_ncomps() const { return user.mesh->spaceDim(); }
  AppCtx const& user;
  double *q_array;
  virtual ~GetDataNormal() {}
};


class GetDataMeshVel : public DefaultGetDataVtk
{
public:
  GetDataMeshVel(double *q_array_, AppCtx const& user_) : user(user_), q_array(q_array_){}
  void get_vec(int id, Real * vec_out) const;
  int vec_ncomps() const { return user.mesh->spaceDim(); }
  AppCtx const& user;
  double *q_array;
  virtual ~GetDataMeshVel() {}
};

class GetDataCellTag : public DefaultGetDataVtk
{
public:
  GetDataCellTag(AppCtx const& user_) : user(user_){}
  int get_data_i(int cellid) const;
  AppCtx const& user;
  virtual ~GetDataCellTag() {}
};




AppCtx::AppCtx(int argc, char **argv, bool &help_return, bool &erro)
{
  setUpDefaultOptions();
  if (getCommandLineOptions(argc, argv) )
    help_return = true;
  else
    help_return = false;


  // create some other things
  erro = createFunctionsSpace(); if (erro) return;
  createQuadrature();

}

bool AppCtx::err_checks()
{
  const bool mesh_has_edge_nodes  = mesh->numNodesPerCell() > mesh->numVerticesPerCell();
  bool const u_has_edge_assoc_dof = shape_phi_c->numDofsAssociatedToFacet() +
                                    shape_phi_c->numDofsAssociatedToCorner() > 0;
  if (!mesh_has_edge_nodes && u_has_edge_assoc_dof)
  {
    printf("ERROR: cant be superparametric element\n");
    return true;
  }
  return false;
}

void AppCtx::loadMesh()
{
  mesh.reset( Mesh::create(ECellType(mesh_cell_type),dim) );
  msh_reader.readFileMsh(filename.c_str(), mesh.get());
  vtk_printer.attachMesh(mesh.get());
  vtk_printer.isFamily(family_files);
  vtk_printer.setOutputFileName(filename_out.c_str());
  vtk_printer.setBinaryOutput( true);

  meshAliasesUpdate();

}

void AppCtx::loadDofs()
{
  dofsCreate();
  timer.restart();
  dofsUpdate();
  timer.elapsed("CuthillMcKeeRenumber");
  n_dofs_u_per_cell   = dof_handler[DH_UNKS].getVariable(VAR_U).numDofsPerCell();
  n_dofs_u_per_facet  = dof_handler[DH_UNKS].getVariable(VAR_U).numDofsPerFacet();
  n_dofs_u_per_corner = dof_handler[DH_UNKS].getVariable(VAR_U).numDofsPerCorner();
  n_dofs_p_per_cell   = dof_handler[DH_UNKS].getVariable(VAR_P).numDofsPerCell();
  n_dofs_p_per_facet  = dof_handler[DH_UNKS].getVariable(VAR_P).numDofsPerFacet();
  n_dofs_p_per_corner = dof_handler[DH_UNKS].getVariable(VAR_P).numDofsPerCorner();
  n_dofs_v_per_cell   = dof_handler[DH_UNKS].getVariable(VAR_M).numDofsPerCell();
  n_dofs_v_per_facet  = dof_handler[DH_UNKS].getVariable(VAR_M).numDofsPerFacet();
  n_dofs_v_per_corner = dof_handler[DH_UNKS].getVariable(VAR_M).numDofsPerCorner();

}

void AppCtx::setUpDefaultOptions()
{
/* global settings */
  dim                    = 2;
  mesh_cell_type         = TRIANGLE3;
  function_space         = 1; // P1P1
  behaviors              = BH_GLS;
  //Re                   = 0.0;
  dt                     = 0.1;
  unsteady               = PETSC_TRUE;
  boundary_smoothing     = PETSC_TRUE;
  steady_tol             = 1.e-6;
  utheta                 = 1;      // time step, theta method (momentum)
  vtheta                 = 1;      // time step, theta method (mesh velocity)
  maxts                  = 10;   // max num of time steps
  finaltime              = -1;
  force_pressure         = PETSC_FALSE;  // elim null space (auto)
  print_to_matlab        = PETSC_FALSE;  // imprime o jacobiano e a função no formato do matlab
  force_dirichlet        = PETSC_TRUE;   // impõe cc ? (para debug)
  full_diriclet          = PETSC_TRUE;
  force_mesh_velocity   = PETSC_FALSE;
  solve_the_sys          = true;   // for debug
  plot_exact_sol         = PETSC_FALSE;
  quadr_degree_cell      = (dim==2) ? 3 : 3;   // ordem de quadratura
  quadr_degree_facet     = 3;
  quadr_degree_corner    = 3;
  quadr_degree_err       = 8; // to compute error
  grow_factor            = 0.05;
  print_step             = 1;
  family_files           = PETSC_TRUE;
  has_convec             = PETSC_TRUE;
  renumber_dofs          = PETSC_FALSE;
  fprint_ca              = PETSC_FALSE;
  nonlinear_elasticity   = PETSC_FALSE;
  mesh_adapt             = PETSC_TRUE;

  filename = (dim==2 ? "malha/cavity2d-1o.msh" : "malha/cavity3d-1o.msh");
}

bool AppCtx::getCommandLineOptions(int argc, char **/*argv*/)
{
  PetscBool          flg_fin, flg_fout;
  char               finaux[PETSC_MAX_PATH_LEN];
  char               foutaux[PETSC_MAX_PATH_LEN];
  PetscBool          ask_help;

  if (argc == 1)
  {
    cout << "\nusage:\n";
    cout << "\n\t./main `cat args`\n\n";
    cout << "where args is a file with the command line parameters, or\n\n";
    cout << "\t./main -help\n\n";
    cout << "to show options.\n\n";
    return true;
  }

  /* opções do usuário */
  PetscOptionsBegin(PETSC_COMM_WORLD, "", "Options for the Navier-Stokes", "ahn?");
  PetscOptionsInt("-dim", "space dimension", "main.cpp", dim, &dim, PETSC_NULL);
  PetscOptionsInt("-mesh_type", "mesh type", "main.cpp", mesh_cell_type, &mesh_cell_type, PETSC_NULL);
  PetscOptionsInt("-function_space", "function_space", "main.cpp", 1, &function_space, PETSC_NULL);
  PetscOptionsInt("-maxts", "maximum number of time steps", "main.cpp", maxts, &maxts, PETSC_NULL);
  //PetscOptionsScalar("-Re", "Reynolds number", "main.cpp", Re, &Re, PETSC_NULL);
  PetscOptionsScalar("-dt", "time step", "main.cpp", dt, &dt, PETSC_NULL);
  PetscOptionsScalar("-utheta", "utheta value", "main.cpp", utheta, &utheta, PETSC_NULL);
  PetscOptionsScalar("-vtheta", "vtheta value", "main.cpp", vtheta, &vtheta, PETSC_NULL);
  PetscOptionsScalar("-sst", "steady state tolerance", "main.cpp", steady_tol, &steady_tol, PETSC_NULL);
  PetscOptionsBool("-print_to_matlab", "print jacobian to matlab", "main.cpp", print_to_matlab, &print_to_matlab, PETSC_NULL);
  PetscOptionsBool("-force_dirichlet", "force dirichlet bound cond", "main.cpp", force_dirichlet, &force_dirichlet, PETSC_NULL);
  PetscOptionsBool("-force_pressure", "force_pressure", "main.cpp", force_pressure, &force_pressure, PETSC_NULL);
  PetscOptionsBool("-plot_es", "plot exact solution", "main.cpp", plot_exact_sol, &plot_exact_sol, PETSC_NULL);
  PetscOptionsBool("-family_files", "plot family output", "main.cpp", family_files, &family_files, PETSC_NULL);
  PetscOptionsBool("-has_convec", "convective term", "main.cpp", has_convec, &has_convec, PETSC_NULL);
  PetscOptionsBool("-unsteady", "unsteady problem", "main.cpp", unsteady, &unsteady, PETSC_NULL);
  PetscOptionsBool("-boundary_smoothing", "boundary_smoothing", "main.cpp", boundary_smoothing, &boundary_smoothing, PETSC_NULL);
  PetscOptionsBool("-force_mesh_velocity", "force_mesh_velocity", "main.cpp", force_mesh_velocity, &force_mesh_velocity, PETSC_NULL);
  PetscOptionsBool("-renumber_dofs", "renumber dofs", "main.cpp", renumber_dofs, &renumber_dofs, PETSC_NULL);
  PetscOptionsBool("-fprint_ca", "print contact angle", "main.cpp", fprint_ca, &fprint_ca, PETSC_NULL);
  PetscOptionsBool("-nonlinear_elasticity", "put a non-linear term in the elasticity problem", "main.cpp", nonlinear_elasticity, &nonlinear_elasticity, PETSC_NULL);
  PetscOptionsBool("-mesh_adapt", "adapt the mesh during simulation", "main.cpp", mesh_adapt, &mesh_adapt, PETSC_NULL);
  PetscOptionsInt("-quadr_e", "quadrature degree (for calculating the error)", "main.cpp", quadr_degree_err, &quadr_degree_err, PETSC_NULL);
  PetscOptionsInt("-quadr_c", "quadrature degree", "main.cpp", quadr_degree_cell, &quadr_degree_cell, PETSC_NULL);
  PetscOptionsInt("-quadr_f", "quadrature degree (facet)", "main.cpp", quadr_degree_facet, &quadr_degree_facet, PETSC_NULL);
  PetscOptionsInt("-quadr_r", "quadrature degree (corner)", "main.cpp", quadr_degree_corner, &quadr_degree_corner, PETSC_NULL);
  PetscOptionsInt("-print_step", "print_step", "main.cpp", print_step, &print_step, PETSC_NULL);
  PetscOptionsScalar("-beta1", "par vel do fluido", "main.cpp", beta1, &beta1, PETSC_NULL);
  PetscOptionsScalar("-beta2", "par vel elastica", "main.cpp", beta2, &beta2, PETSC_NULL);
  PetscOptionsScalar("-finaltime", "the simulation ends at this time.", "main.cpp", finaltime, &finaltime, PETSC_NULL);
  PetscOptionsBool("-ale", "mesh movement", "main.cpp", ale, &ale, PETSC_NULL);
  PetscOptionsGetString(PETSC_NULL,"-fin",finaux,PETSC_MAX_PATH_LEN-1,&flg_fin);
  PetscOptionsGetString(PETSC_NULL,"-fout",foutaux,PETSC_MAX_PATH_LEN-1,&flg_fout);
  PetscOptionsHasName(PETSC_NULL,"-help",&ask_help);

  //is_bdf2 = PETSC_FALSE;
  is_bdf2 = PETSC_TRUE;
  if (is_bdf2 && utheta!=1)
  {
    cout << "ERROR:    BDF2 with utheta!=1" << endl;
    throw;
  }

  if (finaltime < 0)
    finaltime = maxts*dt;
  else
    maxts = 1 + static_cast<int> (round ( finaltime/dt ));

  // get boundary conditions tags (dirichlet)
  dirichlet_tags.resize(16);
  PetscBool flg_tags;
  int nmax = dirichlet_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-dir_tags", dirichlet_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    dirichlet_tags.resize(nmax);
  else
    dirichlet_tags.clear();

  neumann_tags.resize(16);
  nmax = neumann_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-neum_tags", neumann_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    neumann_tags.resize(nmax);
  else
    neumann_tags.clear();

  interface_tags.resize(16);
  nmax = interface_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-interf_tags", interface_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    interface_tags.resize(nmax);
  else
    interface_tags.clear();

  solid_tags.resize(16);
  nmax = solid_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-solid_tags", solid_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    solid_tags.resize(nmax);
  else
    solid_tags.clear();

  triple_tags.resize(16);
  nmax = triple_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-triple_tags", triple_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    triple_tags.resize(nmax);
  else
    triple_tags.clear();

  periodic_tags.resize(16);
  nmax = periodic_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-periodic_tags", periodic_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    periodic_tags.resize(nmax);
  else
    periodic_tags.clear();
  if (periodic_tags.size() % 2 != 0)
  {
    std::cout << "Invalid periodic tags\n";
    throw;
  }

  feature_tags.resize(16);
  nmax = feature_tags.size();
  PetscOptionsGetIntArray(PETSC_NULL, "-feature_tags", feature_tags.data(), &nmax, &flg_tags);
  if (flg_tags)
    feature_tags.resize(nmax);
  else
    feature_tags.clear();



  PetscOptionsEnd();

  switch (function_space)
  {
    case P1P1:
    {
      if(dim==2) mesh_cell_type = TRIANGLE3;
      else       mesh_cell_type = TETRAHEDRON4;
      break;
    }
    case P1bP1_c:
    {
      if(dim==2) mesh_cell_type = TRIANGLE3;
      else       mesh_cell_type = TETRAHEDRON4;
      break;
    }
    case P2bPm1_c:
    {
      if(dim==2) mesh_cell_type = TRIANGLE6;
      else       mesh_cell_type = TETRAHEDRON10;
      break;
    }
    case P2P1:
    {
      if(dim==2) mesh_cell_type = TRIANGLE6;
      else       mesh_cell_type = TETRAHEDRON10;
      break;
    }
    case P1bP1:
    {
      if(dim==2) mesh_cell_type = TRIANGLE3;
      else       mesh_cell_type = TETRAHEDRON4;
      break;
    }
    case P2P0:
    {
      if(dim==2) mesh_cell_type = TRIANGLE6;
      else       mesh_cell_type = TETRAHEDRON10;
      break;
    }
    case P2bPm1:
    {
      if(dim==2) mesh_cell_type = TRIANGLE6;
      else       mesh_cell_type = TETRAHEDRON10;
      break;
    }
    case P1P1unstable:
    {
      if(dim==2) mesh_cell_type = TRIANGLE3;
      else       mesh_cell_type = TETRAHEDRON4;
      break;
    }
    case P2bP1:
    {
      if(dim==2) mesh_cell_type = TRIANGLE6;
      else       mesh_cell_type = TETRAHEDRON10;
      break;
    }
    default:
    {
      cout << "invalid function space " << function_space << endl;
      throw;
    }
  };


  if(ask_help)
    return true;

  if (flg_fin)
    filename.assign(finaux);

  if (flg_fout)
    filename_out.assign(foutaux);

  //if (neumann_tags.size() + interface_tags.size() == 0 || force_pressure)
  if (force_pressure)
  {
    full_diriclet = PETSC_TRUE;
    force_pressure = PETSC_TRUE;
  }
  else
    force_pressure = PETSC_FALSE;

  if (!unsteady)
  {
    if (ale) {
      printf("ale with steady problem?\n");
      //throw;
    }
    //dt = 1.e50;

    if (utheta != 1)
    {
      cout << "WARNING!!!\n";
      cout << ">>> steady problem .. setting utheta to 1" << endl;
      utheta = 1;
    }

  }

  return false;
}

bool AppCtx::createFunctionsSpace()
{
  EShapeType velo_shape, pres_shape;
  bool is_simplex = ctype2cfamily(ECellType(mesh_cell_type)) == SIMPLEX;
  ECellType cell_type = ECellType(mesh_cell_type);


  switch (function_space)
  {
    case 1: // P1P1 (or Q1Q1) GLS stabilization
      {
        behaviors = BH_GLS;
        velo_shape = is_simplex ? P1 : Q1;
        pres_shape = is_simplex ? P1 : Q1;
      }
      break;
    case 2: // P1+P1 with bubble condensation
      {
        behaviors = BH_bble_condens_PnPn;
        velo_shape = P1;
        pres_shape = P1;
      }
      break;
    case 3: // P2+Pm1 with bubble condensation and pressure gradient elimination
      {
        behaviors = BH_bble_condens_CR;
        velo_shape = P2;
        pres_shape = P0;
      }
      break;
    case 4: // P2P1 (or Q2Q1)
      {
        behaviors = 0;
        velo_shape = is_simplex ? P2 : Q2;
        pres_shape = is_simplex ? P1 : Q1;
      }
      break;
    case 5: // P1+P1 traditional
      {
        behaviors = 0;
        velo_shape = P1ph ;
        pres_shape = P1;
      }
      break;
    case 6: // P2P0
      {
        behaviors = 0;
        velo_shape = P2;
        pres_shape = P0;
      }
      break;
    case 7: // P2+Pm1 full
      {
        behaviors = 0;
        velo_shape = P2ph;
        pres_shape = Pm1;
      }
      break;
    case 8: // P1P1 unstable
      {
        behaviors = 0;
        velo_shape = P1;
        pres_shape = P1;
      }
      break;
    case 9: // P2bP1 bubble condensation
      {
        behaviors = BH_bble_condens_PnPn;
        velo_shape = is_simplex ? P2 : Q2;
        pres_shape = is_simplex ? P1 : Q1;
      }
      break;
    default:
      {
        behaviors = BH_GLS;
        velo_shape = is_simplex ? P1 : Q1;
        pres_shape = is_simplex ? P1 : Q1;
      }

  }


  shape_phi_c.reset(ShapeFunction::create(cell_type, velo_shape));
  shape_psi_c.reset(ShapeFunction::create(cell_type, pres_shape));
  shape_qsi_c.reset(ShapeFunction::create(cell_type));
  shape_phi_f.reset(ShapeFunction::create(facetof(cell_type), facetof(velo_shape)));
  shape_psi_f.reset(ShapeFunction::create(facetof(cell_type), facetof(pres_shape)));
  shape_qsi_f.reset(ShapeFunction::create(facetof(cell_type)));
  shape_phi_r.reset(ShapeFunction::create(facetof(facetof(cell_type)), facetof(facetof(velo_shape))));
  shape_psi_r.reset(ShapeFunction::create(facetof(facetof(cell_type)), facetof(facetof(pres_shape))));
  shape_qsi_r.reset(ShapeFunction::create(facetof(facetof(cell_type))));

  shape_bble.reset(ShapeFunction::create(cell_type, BUBBLE));

  pres_pres_block = false;
  if (behaviors & (BH_bble_condens_PnPn | BH_GLS))
    pres_pres_block = true;


  return false;
}

void AppCtx::createQuadrature()
{
  ECellType ct = ECellType(mesh_cell_type);

  quadr_cell.reset( Quadrature::create(ECellType(ct)) );
  quadr_cell->setOrder(quadr_degree_cell);
  n_qpts_cell = quadr_cell->numPoints();

  quadr_facet.reset( Quadrature::create(facetof(ECellType(ct))) );
  quadr_facet->setOrder(quadr_degree_facet);
  n_qpts_facet = quadr_facet->numPoints();

  quadr_corner.reset( Quadrature::create(facetof(facetof(ECellType(ct)))) );
  quadr_corner->setOrder(quadr_degree_corner);
  n_qpts_corner = quadr_corner->numPoints();


  quadr_err.reset( Quadrature::create(ECellType(ct)) );
  quadr_err->setOrder(quadr_degree_err);
  n_qpts_err = quadr_err->numPoints();

}

void AppCtx::meshAliasesUpdate()
{
  n_nodes  = mesh->numNodes();
  n_cells  = mesh->numCells();
  n_facets = mesh->numFacets();
  n_corners= mesh->numCorners();
  n_nodes_total  = mesh->numNodesTotal();   // includes disableds
  n_cells_total  = mesh->numCellsTotal();   // includes disableds
  n_facets_total = mesh->numFacetsTotal();  // includes disableds
  n_corners_total= mesh->numCornersTotal(); // includes disableds
  nodes_per_cell  = mesh->numNodesPerCell();
  nodes_per_facet = mesh->numNodesPerFacet();
  nodes_per_corner= mesh->numNodesPerCorner();
}

void AppCtx::dofsCreate()
{
  // dof handler create
  dof_handler[DH_UNKS].setMesh(mesh.get());
  dof_handler[DH_UNKS].addVariable("velo",  shape_phi_c.get(), dim);
  dof_handler[DH_UNKS].addVariable("pres",  shape_psi_c.get(), 1);
  dof_handler[DH_UNKS].getVariable(VAR_P).setType(SPLITTED_BY_REGION_CELL,0,0);
  //Matrix<bool, Dynamic, Dynamic> blocks(2,2);
  //blocks.setOnes();
  //blocks(1,1)=pres_pres_block;
  //dof_handler[DH_UNKS].setVariablesRelationship(blocks.data());

  // mesh velocity
  dof_handler[DH_MESH].setMesh(mesh.get());
  dof_handler[DH_MESH].addVariable("mesh_veloc",  shape_qsi_c.get(), dim);
  //dof_handler[DH_MESH].setVariablesRelationship(blocks.data());
}

bool isPeriodic(Point const* p1, Point const* p2, int dim)
{
  double const TOL = 1.e-9;

  if (dim == 2)
  {
    if (fabs(p1->getCoord(0) - p2->getCoord(0)) < TOL)
      return true;
    if (fabs(p1->getCoord(1) - p2->getCoord(1)) < TOL)
      return true;
  }
  else // dim == 3
  {
    int c = 0;
    if (fabs(p1->getCoord(0) - p2->getCoord(0)) < TOL)
      c++;
    if (fabs(p1->getCoord(1) - p2->getCoord(1)) < TOL)
      c++;
    if (fabs(p1->getCoord(2) - p2->getCoord(2)) < TOL)
      c++;

    if (c==2)
      return true;
  }

  return false;
}

void AppCtx::dofsUpdate()
{
  dof_handler[DH_UNKS].SetUp();
  //if (renumber_dofs)
  //  dof_handler[DH_UNKS].CuthillMcKeeRenumber();

  std::vector<int> dofs1;
  std::vector<int> dofs2;

  int dofs_a[10];
  int dofs_b[10];

  // apply periodic boundary conditions here
  {
    if ((mesh_cell_type != TRIANGLE3) && (mesh_cell_type != TETRAHEDRON4) && !periodic_tags.empty() )
    {
      std::cout << "Periodic boundary conditions is not allowed with high order mesh\n";
      throw;
    }
    
      point_iterator point1 = mesh->pointBegin();
      point_iterator point1_end = mesh->pointEnd();
      point_iterator point2, point2_end;

      for (; point1 != point1_end; ++point1)
      {
        for (int i = 0; i < (int)periodic_tags.size(); i+=2)
        {

          int const tag1 = periodic_tags[i];
          int const tag2 = periodic_tags[i+1];

          if (tag1 != point1->getTag())
            continue;

          point2 = mesh->pointBegin();
          point2_end = mesh->pointEnd();
          for (; point2 != point2_end; ++point2)
          {
            if (tag2 != point2->getTag())
              continue;

            if (isPeriodic(&*point1, &*point2, dim))
            {
              for (int var = 0; var < 2; ++var)
              {
                int ndpv = dof_handler[DH_UNKS].getVariable(var).numDofsPerVertex();

                dof_handler[DH_UNKS].getVariable(var).getVertexDofs(dofs_a, &*point1);
                dof_handler[DH_UNKS].getVariable(var).getVertexDofs(dofs_b, &*point2);
                for (int k = 0; k < ndpv; ++k)
                {
                  dofs1.push_back(dofs_a[k]);
                  dofs2.push_back(dofs_b[k]);
                }
              }

            }
          }
        }
      }
  }
  dof_handler[DH_UNKS].linkDofs(dofs1.size(), dofs1.data(), dofs2.data());
  n_unknowns = dof_handler[DH_UNKS].numDofs();

  dof_handler[DH_MESH].SetUp();
  n_dofs_v_mesh = dof_handler[DH_MESH].numDofs();
}

PetscErrorCode AppCtx::allocPetscObjs()
{
  printf( "allocing petsc objs ... ");

  PetscErrorCode      ierr;
//  ierr = SNESCreate(PETSC_COMM_WORLD, &snes);                   CHKERRQ(ierr);

  //Vec Vec_res;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_res);                     CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_res, PETSC_DECIDE, n_unknowns);            CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_res);                                CHKERRQ(ierr);

  //Vec Vec_up_0;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_up_0);                      CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_up_0, PETSC_DECIDE, n_unknowns);             CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_up_0);                                 CHKERRQ(ierr);

  //Vec Vec_up_1;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_up_1);                       CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_up_1, PETSC_DECIDE, n_unknowns);              CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_up_1);                                  CHKERRQ(ierr);

  //Vec Vec_dup;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_dup);                       CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_dup, PETSC_DECIDE, n_unknowns);              CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_dup);                                  CHKERRQ(ierr);

  //Vec Vec_v_mid
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_v_mid);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_v_mid, PETSC_DECIDE, n_dofs_v_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_v_mid);                             CHKERRQ(ierr);

  //Vec Vec_v_1
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_v_1);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_v_1, PETSC_DECIDE, n_dofs_v_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_v_1);                             CHKERRQ(ierr);

  //Vec Vec_x_0;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_x_0);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_x_0, PETSC_DECIDE, n_dofs_v_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_x_0);                             CHKERRQ(ierr);

  //Vec Vec_x_1;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_x_1);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_x_1, PETSC_DECIDE, n_dofs_v_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_x_1);                             CHKERRQ(ierr);

  //Vec Vec_normal;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_normal);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_normal, PETSC_DECIDE, n_dofs_v_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_normal);                             CHKERRQ(ierr);

  //Vec Vec_res_m;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_res_m);                     CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_res_m, PETSC_DECIDE, n_dofs_v_mesh);         CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_res_m);                                CHKERRQ(ierr);


  std::vector<int> nnz;
  //if(false)
  {
    nnz.resize(n_unknowns,0);
    
    std::vector<SetVector<int> > table;
    dof_handler[DH_UNKS].getSparsityTable(table); // TODO: melhorar desempenho, função mt lenta


    //FEP_PRAGMA_OMP(parallel for)
    for (int i = 0; i < n_unknowns; ++i)
      nnz[i] = table[i].size();

    //// removendo a diagonal nula, pois o petsc da erro se for 0
    //if (!pres_pres_block)
    //{
    //  int const n_p_dofs_total = dof_handler[DH_UNKS].getVariable(VAR_P).totalSize();
    //  //FEP_PRAGMA_OMP(parallel for)
    //  for (int i = 0; i < n_p_dofs_total; ++i)
    //  {
    //    int const dof = dof_handler[DH_UNKS].getVariable(VAR_P).data()[i];
    //    if (dof >= 0)
    //      ++nnz[dof];
    //  }
    //}


  }

  // opções para nnz:
  // TETRA10 = 375

  //Mat Mat_Jac;
  ierr = MatCreate(PETSC_COMM_WORLD, &Mat_Jac);                                      CHKERRQ(ierr);
  ierr = MatSetSizes(Mat_Jac, PETSC_DECIDE, PETSC_DECIDE, n_unknowns, n_unknowns);   CHKERRQ(ierr);
  ierr = MatSetFromOptions(Mat_Jac);                                                 CHKERRQ(ierr);
  //ierr = MatSeqAIJSetPreallocation(Mat_Jac, 0, nnz.data());                          CHKERRQ(ierr);
  //max_nz = nnz.maxCoeff();
  ierr = MatSeqAIJSetPreallocation(Mat_Jac, 0, nnz.data());       CHKERRQ(ierr);
  //ierr = MatSetOption(Mat_Jac,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);                  CHKERRQ(ierr);
  ierr = MatSetOption(Mat_Jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE); CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD, &snes);                  CHKERRQ(ierr);
  ierr = SNESSetFunction(snes, Vec_res, FormFunction, this);      CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes, Mat_Jac, Mat_Jac, FormJacobian, this); CHKERRQ(ierr);
  //ierr = SNESSetJacobian(snes,Mat_Jac,Mat_Jac,SNESDefaultComputeJacobian,&user);  CHKERRQ(ierr);

  ierr = SNESSetConvergenceTest(snes,CheckSnesConvergence,this,PETSC_NULL); CHKERRQ(ierr);

  ierr = SNESGetKSP(snes,&ksp);                                                  CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);                                                      CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,Mat_Jac,Mat_Jac,SAME_NONZERO_PATTERN);                       CHKERRQ(ierr);
  //~ ierr = KSPSetType(ksp,KSPPREONLY);                                           CHKERRQ(ierr);
  //~ ierr = KSPSetType(ksp,KSPGMRES);                                               CHKERRQ(ierr);
  //~ ierr = PCSetType(pc,PCLU);                                                     CHKERRQ(ierr);
  //~ ierr = PCFactorSetMatOrderingType(pc, MATORDERINGNATURAL);                         CHKERRQ(ierr);
  //~ ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);  CHKERRQ(ierr);
  //~ ierr = SNESSetApplicationContext(snes,this);

//~ #ifdef PETSC_HAVE_MUMPS
  //~ PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
//~ #endif

  //ierr = SNESMonitorSet(snes, SNESMonitorDefault, 0, 0); CHKERRQ(ierr);
  //ierr = SNESMonitorSet(snes,Monitor,0,0);CHKERRQ(ierr);
  //ierr = SNESSetTolerances(snes,0,0,0,13,PETSC_DEFAULT);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
  //ierr = SNESLineSearchSet(snes, SNESLineSearchNo, &user); CHKERRQ(ierr);


  //  ------------------------------------------------------------------
  //  ----------------------- Mesh -------------------------------------
  //  ------------------------------------------------------------------
  nnz.clear();
  int n_mesh_dofs = dof_handler[DH_MESH].numDofs();
  {
    std::vector<SetVector<int> > table;
    dof_handler[DH_MESH].getSparsityTable(table); // TODO: melhorar desempenho, função mt lenta

    nnz.resize(n_mesh_dofs, 0);

    //FEP_PRAGMA_OMP(parallel for)
    for (int i = 0; i < n_mesh_dofs; ++i)
      nnz[i] = table[i].size();

  }


  //Mat Mat_Jac_m;
  ierr = MatCreate(PETSC_COMM_WORLD, &Mat_Jac_m);                                        CHKERRQ(ierr);
  ierr = MatSetType(Mat_Jac_m,MATSEQAIJ);                                                CHKERRQ(ierr);
  ierr = MatSetSizes(Mat_Jac_m, PETSC_DECIDE, PETSC_DECIDE, n_mesh_dofs, n_mesh_dofs);   CHKERRQ(ierr);
  //ierr = MatSetFromOptions(Mat_Jac_m);                                                 CHKERRQ(ierr);
  //ierr = MatSeqAIJSetPreallocation(Mat_Jac_m, 0, nnz.data());                          CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(Mat_Jac_m,  0, nnz.data());         CHKERRQ(ierr);
  //ierr = MatSetOption(Mat_Jac_m,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);                  CHKERRQ(ierr);
  ierr = MatSetOption(Mat_Jac_m,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);             CHKERRQ(ierr);
  ierr = MatSetOption(Mat_Jac_m,MAT_SYMMETRIC,PETSC_TRUE);                               CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD, &snes_m);                                          CHKERRQ(ierr);
  ierr = SNESSetFunction(snes_m, Vec_res_m, FormFunction_mesh, this);                      CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes_m, Mat_Jac_m, Mat_Jac_m, FormJacobian_mesh, this);         CHKERRQ(ierr);

  // descomente abaixo para sistemas lineares
  //ierr = SNESSetType(snes_m, SNESKSPONLY);                                               CHKERRQ(ierr);

  SNESGetLineSearch(snes_m,&linesearch);
  SNESLineSearchSetType(linesearch,SNESLINESEARCHBASIC);

  ierr = SNESGetKSP(snes_m,&ksp_m);                                                  CHKERRQ(ierr);
  ierr = KSPGetPC(ksp_m,&pc_m);                                                      CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp_m,Mat_Jac_m,Mat_Jac_m,SAME_NONZERO_PATTERN);            CHKERRQ(ierr);
  // ierr = KSPSetType(ksp_m,KSPPREONLY);                                            CHKERRQ(ierr);
  //ierr = KSPSetType(ksp_m,KSPGMRES);                                                 CHKERRQ(ierr);
  ierr = KSPSetType(ksp_m,KSPPREONLY);                                                 CHKERRQ(ierr);
  ierr = PCSetType(pc_m,PCLU);                                                   CHKERRQ(ierr);
  // ierr = PCFactorSetMatOrderingType(pc_m, MATORDERINGNATURAL);                         CHKERRQ(ierr);
  // ierr = KSPSetTolerances(ksp_m,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);  CHKERRQ(ierr);
  // ierr = SNESSetApplicationContext(snes_m,this);

  if(!nonlinear_elasticity)
  {
    ierr = SNESSetType(snes_m, SNESKSPONLY); CHKERRQ(ierr);
  }


//~ #ifdef PETSC_HAVE_MUMPS
  //~ PCFactorSetMatSolverPackage(pc_m,MATSOLVERMUMPS);
//~ #endif

  //ierr = SNESMonitorSet(snes_m, SNESMonitorDefault, 0, 0); CHKERRQ(ierr);
  //ierr = SNESMonitorSet(snes_m,Monitor,0,0);CHKERRQ(ierr);
  //ierr = SNESSetTolerances(snes_m,0,0,0,13,PETSC_DEFAULT);
  //ierr = SNESSetFromOptions(snes_m); CHKERRQ(ierr);
  //ierr = SNESLineSearchSet(snes_m, SNESLineSearchNo, &user); CHKERRQ(ierr);



  printf("done.\n");
  PetscFunctionReturn(0);
}

void AppCtx::matrixColoring()
{
  printf("matrix coloring ... ");

  VectorXi                     mapU_c(n_dofs_u_per_cell);
  VectorXi                     mapU_f(n_dofs_u_per_facet);
  VectorXi                     mapU_r(n_dofs_u_per_corner);
  VectorXi                     mapP_c(n_dofs_p_per_cell);
  VectorXi                     mapP_f(n_dofs_p_per_facet);
  VectorXi                     mapP_r(n_dofs_p_per_corner);
  // mesh velocity
  VectorXi                     mapM_c(dim*nodes_per_cell);
  VectorXi                     mapM_f(dim*nodes_per_facet);
  VectorXi                     mapM_r(dim*nodes_per_corner);

  MatrixXd Aloc = MatrixXd::Zero(n_dofs_u_per_cell, n_dofs_u_per_cell);
  MatrixXd Dloc = MatrixXd::Zero(n_dofs_p_per_cell, n_dofs_u_per_cell);
  MatrixXd Gloc = MatrixXd::Zero(n_dofs_u_per_cell, n_dofs_p_per_cell);
  MatrixXd Eloc = MatrixXd::Zero(n_dofs_p_per_cell, n_dofs_p_per_cell);

  cell_iterator cell = mesh->cellBegin();
  cell_iterator cell_end = mesh->cellEnd();
  for (; cell != cell_end; ++cell)
  {
    // mapeamento do local para o global:
    dof_handler[DH_UNKS].getVariable(VAR_U).getCellDofs(mapU_c.data(), &*cell);
    dof_handler[DH_UNKS].getVariable(VAR_P).getCellDofs(mapP_c.data(), &*cell);

    MatSetValues(Mat_Jac, mapU_c.size(), mapU_c.data(), mapU_c.size(), mapU_c.data(), Aloc.data(), ADD_VALUES);
    MatSetValues(Mat_Jac, mapU_c.size(), mapU_c.data(), mapP_c.size(), mapP_c.data(), Gloc.data(), ADD_VALUES);
    MatSetValues(Mat_Jac, mapP_c.size(), mapP_c.data(), mapU_c.size(), mapU_c.data(), Dloc.data(), ADD_VALUES);
    //if (pres_pres_block)
    MatSetValues(Mat_Jac, mapP_c.size(), mapP_c.data(), mapP_c.size(), mapP_c.data(), Eloc.data(), ADD_VALUES);

  }

  ////test
  //for (int i = 0; i < n_unknowns; ++i)
  //  for (int j = 0; j < n_unknowns; ++j)
  //    MatSetValue(Mat_Jac, i, j, 0.0, ADD_VALUES);

  //if (!pres_pres_block)
  //{
  //  int const n_p_dofs_total = dof_handler[DH_UNKS].getVariable(VAR_P).totalSize();
  //  for (int i = 0; i < n_p_dofs_total; ++i)
  //  {
  //    const double zero = 0.0;
  //    int const dof = dof_handler[DH_UNKS].getVariable(VAR_P).data()[i];
  //    if (dof>=0)
  //      MatSetValue(Mat_Jac, dof, dof, zero, ADD_VALUES);
  //  }
  //}

  Assembly(Mat_Jac);
  //MatSetOption(Mat_Jac,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
  //MatSetOption(Mat_Jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);
  printf("done \n");
}

void AppCtx::printMatlabLoader()
{
  FILE *fp = fopen("loadmat.m", "w");
  //fprintf(fp, "clear;\n"                   );
  fprintf(fp, "jacob;\n"                   );
  fprintf(fp, "clear zzz;\n"               );
  fprintf(fp, "B=Jac;\n"                   );
  fprintf(fp, "B(B!=0)=1;\n"               );
  fprintf(fp, "nU = %d;\n",dof_handler[DH_UNKS].getVariable(VAR_U).numPositiveDofs() );
  fprintf(fp, "nP = %d;\n",dof_handler[DH_UNKS].getVariable(VAR_P).numPositiveDofs() );
  fprintf(fp, "nT = nU + nP;\n"            );
  fprintf(fp, "K=Jac(1:nU,1:nU);\n"        );
  fprintf(fp, "G=Jac(1:nU,nU+1:nT);\n"     );
  fprintf(fp, "D=Jac(nU+1:nT,1:nU);\n"     );
  fprintf(fp, "E=Jac(nU+1:nT,nU+1:nT);\n"  );
  fprintf(fp, "\n"                        );
  fprintf(fp, "rhs;\n"                     );
  fprintf(fp, "f=res(1:nU);\n"             );
  fprintf(fp, "g=res(nU+1:nT);\n"          );
  fclose(fp);
}

// must be called after loadDofs
void AppCtx::evaluateQuadraturePts()
{
  // avaliando phi_c e as matrizes de transf nos pontos de quadratura.
  phi_c.resize(n_qpts_cell);
  psi_c.resize(n_qpts_cell);
  qsi_c.resize(n_qpts_cell);

  phi_f.resize(n_qpts_facet);
  psi_f.resize(n_qpts_facet);
  qsi_f.resize(n_qpts_facet);

  phi_r.resize(n_qpts_corner);
  psi_r.resize(n_qpts_corner);
  qsi_r.resize(n_qpts_corner);

  dLphi_c.resize(n_qpts_cell);
  dLpsi_c.resize(n_qpts_cell);
  dLqsi_c.resize(n_qpts_cell);

  dLphi_f.resize(n_qpts_facet);
  dLpsi_f.resize(n_qpts_facet);
  dLqsi_f.resize(n_qpts_facet);

  dLphi_r.resize(n_qpts_corner);
  dLpsi_r.resize(n_qpts_corner);
  dLqsi_r.resize(n_qpts_corner);

  bble.resize(n_qpts_cell);
  dLbble.resize(n_qpts_cell);


  for (int qp = 0; qp < n_qpts_cell; ++qp)
  {
    phi_c[qp].resize(n_dofs_u_per_cell/dim);
    psi_c[qp].resize(n_dofs_p_per_cell);
    qsi_c[qp].resize(nodes_per_cell);

    dLphi_c[qp].resize(n_dofs_u_per_cell/dim, dim);
    dLpsi_c[qp].resize(n_dofs_p_per_cell, dim);
    dLqsi_c[qp].resize(nodes_per_cell, dim);

    dLbble[qp].resize(dim);
    bble[qp] = shape_bble->eval(quadr_cell->point(qp), 0);

    for (int n = 0; n < n_dofs_u_per_cell/dim; ++n)
    {
      phi_c[qp][n] = shape_phi_c->eval(quadr_cell->point(qp), n);
      for (int d = 0; d < dim; ++d)
      {
        /* dLphi_c nao depende de qp no caso de funcoes lineares */
        dLphi_c[qp](n, d) = shape_phi_c->gradL(quadr_cell->point(qp), n, d);
      }
    }

    for (int n = 0; n < n_dofs_p_per_cell; ++n)
    {
      psi_c[qp][n] = shape_psi_c->eval(quadr_cell->point(qp), n);
      for (int d = 0; d < dim; ++d)
      {
        /* dLpsi_c nao depende de qp no caso de funcoes lineares */
        dLpsi_c[qp](n, d) = shape_psi_c->gradL(quadr_cell->point(qp), n, d);
      }
    }

    for (int n = 0; n < nodes_per_cell; ++n)
    {
      qsi_c[qp][n] = shape_qsi_c->eval(quadr_cell->point(qp), n);
      for (int d = 0; d < dim; ++d)
      {
        /* dLqsi_c nao depende de qp no caso de funcoes lineares */
        dLqsi_c[qp](n, d) = shape_qsi_c->gradL(quadr_cell->point(qp), n, d);
      }
    }

    for (int d = 0; d < dim; ++d)
    {
      dLbble[qp](d) = shape_bble->gradL(quadr_cell->point(qp), 0, d);
    }

  } // end qp

  // pts de quadratura no contorno
  for (int qp = 0; qp < n_qpts_facet; ++qp)
  {
    phi_f[qp].resize(n_dofs_u_per_facet/dim);
    psi_f[qp].resize(n_dofs_p_per_facet);
    qsi_f[qp].resize(n_dofs_v_per_facet/dim);

    dLphi_f[qp].resize(n_dofs_u_per_facet/dim, dim-1);
    dLpsi_f[qp].resize(n_dofs_p_per_facet, dim-1);
    dLqsi_f[qp].resize(n_dofs_v_per_facet/dim, dim-1);


    for (int n = 0; n < n_dofs_u_per_facet/dim; ++n)
    {
      phi_f[qp][n] = shape_phi_f->eval(quadr_facet->point(qp), n);
      for (int d = 0; d < dim-1; ++d)
      {
        /* dLphi_f nao depende de qp no caso de funcoes lineares */
        dLphi_f[qp](n, d) = shape_phi_f->gradL(quadr_facet->point(qp), n, d);
      }
    }

    for (int n = 0; n < n_dofs_p_per_facet; ++n)
    {
      psi_f[qp][n] = shape_psi_f->eval(quadr_facet->point(qp), n);
      for (int d = 0; d < dim-1; ++d)
      {
        /* dLpsi_f nao depende de qp no caso de funcoes lineares */
        dLpsi_f[qp](n, d) = shape_psi_f->gradL(quadr_facet->point(qp), n, d);
      }
    }

    for (int n = 0; n < nodes_per_facet; ++n)
    {
      qsi_f[qp][n] = shape_qsi_f->eval(quadr_facet->point(qp), n);
      for (int d = 0; d < dim-1; ++d)
      {
        /* dLqsi_f nao depende de qp no caso de funcoes lineares */
        dLqsi_f[qp](n, d) = shape_qsi_f->gradL(quadr_facet->point(qp), n, d);
      }
    }



  } // end qp



  // pts de quadratura no contorno do contorno
  for (int qp = 0; qp < n_qpts_corner; ++qp)
  {
    phi_r[qp].resize(n_dofs_u_per_corner/dim);
    psi_r[qp].resize(n_dofs_p_per_corner);
    qsi_r[qp].resize(nodes_per_corner);

    dLphi_r[qp].resize(n_dofs_u_per_corner/dim, 1 /* = dim-2*/);
    dLpsi_r[qp].resize(n_dofs_p_per_corner, 1 /* = dim-2*/);
    dLqsi_r[qp].resize(nodes_per_corner, 1 /* = dim-2*/);

    for (int n = 0; n < n_dofs_u_per_corner/dim; ++n)
    {
      phi_r[qp][n] = shape_phi_r->eval(quadr_corner->point(qp), n);
      for (int d = 0; d < dim-2; ++d)
      {
        /* dLphi_r nao depende de qp no caso de funcoes lineares */
        dLphi_r[qp](n, d) = shape_phi_r->gradL(quadr_corner->point(qp), n, d);
      }
    }

    for (int n = 0; n < n_dofs_p_per_corner; ++n)
    {
      psi_r[qp][n] = shape_psi_r->eval(quadr_corner->point(qp), n);
      for (int d = 0; d < dim-2; ++d)
      {
        /* dLpsi_r nao depende de qp no caso de funcoes lineares */
        dLpsi_r[qp](n, d) = shape_psi_r->gradL(quadr_corner->point(qp), n, d);
      }
    }

    for (int n = 0; n < nodes_per_corner; ++n)
    {
      qsi_r[qp][n] = shape_qsi_r->eval(quadr_corner->point(qp), n);
      for (int d = 0; d < dim-2; ++d)
      {
        /* dLqsi_r nao depende de qp no caso de funcoes lineares */
        dLqsi_r[qp](n, d) = shape_qsi_r->gradL(quadr_corner->point(qp), n, d);
      }
    }
  } // end qp

  //     Quadrature
  //     to compute the error
  //
  // velocity
  phi_err.resize(n_qpts_err);         // shape function evaluated at quadrature points
  dLphi_err.resize(n_qpts_err);       // matriz de gradiente no elemento unitário
  // pressure
  psi_err.resize(n_qpts_err);         // shape function evaluated at quadrature points
  dLpsi_err.resize(n_qpts_err);       // matriz de gradiente no elemento unitário
  // mesh
  qsi_err.resize(n_qpts_err);         // shape function evaluated at quadrature points
  dLqsi_err.resize(n_qpts_err);       // matriz de gradiente no elemento unitário

  for (int qp = 0; qp < n_qpts_err; ++qp)
  {
    phi_err[qp].resize(n_dofs_u_per_cell/dim);
    psi_err[qp].resize(n_dofs_p_per_cell);
    qsi_err[qp].resize(nodes_per_cell);

    dLphi_err[qp].resize(n_dofs_u_per_cell/dim, dim);
    dLpsi_err[qp].resize(n_dofs_p_per_cell, dim);
    dLqsi_err[qp].resize(nodes_per_cell, dim);

    for (int n = 0; n < n_dofs_u_per_cell/dim; ++n)
    {
      phi_err[qp][n] = shape_phi_c->eval(quadr_err->point(qp), n);
      for (int d = 0; d < dim; ++d)
      {
        /* dLphi_err nao depende de qp no caso de funcoes lineares */
        dLphi_err[qp](n, d) = shape_phi_c->gradL(quadr_err->point(qp), n, d);
      }
    }

    for (int n = 0; n < n_dofs_p_per_cell; ++n)
    {
      psi_err[qp][n] = shape_psi_c->eval(quadr_err->point(qp), n);
      for (int d = 0; d < dim; ++d)
      {
        /* dLpsi_err nao depende de qp no caso de funcoes lineares */
        dLpsi_err[qp](n, d) = shape_psi_c->gradL(quadr_err->point(qp), n, d);
      }
    }

    for (int n = 0; n < nodes_per_cell; ++n)
    {
      qsi_err[qp][n] = shape_qsi_c->eval(quadr_err->point(qp), n);
      for (int d = 0; d < dim; ++d)
      {
        /* dLqsi_err nao depende de qp no caso de funcoes lineares */
        dLqsi_err[qp](n, d) = shape_qsi_c->gradL(quadr_err->point(qp), n, d);
      }
    }


  } // end qp

  Vector Xc = Vector::Zero(dim);
  bool is_simplex = ctype2cfamily(ECellType(mesh_cell_type)) == SIMPLEX;
  if (is_simplex)
    for (int i = 0; i < dim; ++i)
      Xc(i) = 1./(dim+1);
  qsi_c_at_center.resize(nodes_per_cell);
  for (int i = 0; i < nodes_per_cell; ++i)
    qsi_c_at_center(i) = shape_qsi_c->eval(Xc.data(), i);

  // facets func derivatives on your nodes
  if (dim==2)
  {
    std::vector<double> parametric_pts;
    parametric_pts = genLineParametricPts(  ctypeDegree(ECellType(mesh_cell_type))  );

    dLphi_nf.resize(parametric_pts.size());

    for (int k = 0; k < (int)parametric_pts.size(); ++k)
    {
      dLphi_nf[k].resize(n_dofs_u_per_facet/dim, dim-1);
      for (int n = 0; n < n_dofs_u_per_facet/dim; ++n)
        for (int d = 0; d < dim-1; ++d)
          dLphi_nf[k](n, d) = shape_phi_f->gradL(&parametric_pts[k], n, d);
    }
  }
  else
  if(dim==3)
  {
    std::vector<Eigen::Vector2d> parametric_pts;
    parametric_pts = genTriParametricPts(  ctypeDegree(ECellType(mesh_cell_type))  );

    dLphi_nf.resize(parametric_pts.size());

    for (int k = 0; k < (int)parametric_pts.size(); ++k)
    {
      dLphi_nf[k].resize(n_dofs_u_per_facet/dim, dim-1);
      for (int n = 0; n < n_dofs_u_per_facet/dim; ++n)
        for (int d = 0; d < dim-1; ++d)
          dLphi_nf[k](n, d) = shape_phi_f->gradL(parametric_pts[k].data(), n, d);
    }
  }



}

void AppCtx::onUpdateMesh()
{
  allocPetscObjs();
  matrixColoring();
}

PetscErrorCode AppCtx::setInitialConditions()
{
  PetscErrorCode      ierr(0);

  Vector    Uf(dim);
  //double    pf;
  Vector    X(dim);
  Tensor    R(dim,dim);
  VectorXi  dofs(dim);
  //int       dof;
  int tag;

  current_time = 0;

  VecZeroEntries(Vec_res);
  VecZeroEntries(Vec_up_0);
  VecZeroEntries(Vec_up_1);
  VecZeroEntries(Vec_v_mid);
  VecZeroEntries(Vec_x_0);
  VecZeroEntries(Vec_x_1);
  VecZeroEntries(Vec_normal);


  copyMesh2Vec(Vec_x_0);
  copyMesh2Vec(Vec_x_1);

  // normals
  getVecNormals(&Vec_x_0, Vec_normal);

  // compute mesh sizes
  {
    mesh_sizes.reserve(n_nodes_total + 1*n_nodes_total/10);
    mesh_sizes.assign(n_nodes_total, 0.0);

    // stores how much edges were connecteds until now at each vertex
    std::vector<int> counter(n_nodes_total, 0);

    VectorXi edge_nodes(3);
    Vector Xa(dim), Xb(dim);
    Real h;

    const int n_edges_total = (dim==2) ? mesh->numFacetsTotal() : mesh->numCornersTotal();
    CellElement const* edge(NULL);

    //FEP_PRAGMA_OMP(for nowait)
    for (int a = 0; a < n_edges_total; ++a)
    {
      if (dim == 2)
        edge = mesh->getFacetPtr(a);
      else
        edge = mesh->getCornerPtr(a);

      if (edge==NULL || edge->isDisabled())
        continue;

      if (dim == 2)
        mesh->getFacetNodesId(a, edge_nodes.data());
      else
        mesh->getCornerNodesId(a, edge_nodes.data());

      mesh->getNodePtr(edge_nodes[0])->getCoord(Xa.data(),dim);
      mesh->getNodePtr(edge_nodes[1])->getCoord(Xb.data(),dim);
      h = (Xa-Xb).norm();

      mesh_sizes[edge_nodes[0]] = (mesh_sizes[edge_nodes[0]]*counter[edge_nodes[0]] + h)/(counter[edge_nodes[0]] + 1.0);
      mesh_sizes[edge_nodes[1]] = (mesh_sizes[edge_nodes[1]]*counter[edge_nodes[1]] + h)/(counter[edge_nodes[1]] + 1.0);

      ++counter[edge_nodes[0]];
      ++counter[edge_nodes[1]];

    } // end n_edges_total


  } // end compute mesh sizes



  // velocidade inicial e pressao inicial
  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    tag = point->getTag();

    point->getCoord(X.data(),dim);

    Uf = u_initial(X, tag);

    getNodeDofs(&*point,DH_UNKS,VAR_U,dofs.data());

    VecSetValues(Vec_up_0, dim, dofs.data(), Uf.data(), INSERT_VALUES);

    // press
    if (mesh->isVertex(&*point))
    {
      getNodeDofs(&*point,DH_UNKS,VAR_P,dofs.data());
      double p_in = p_initial(X, tag);
      VecSetValues(Vec_up_0, 1, dofs.data(), &p_in, INSERT_VALUES);
    }


  } // end point loop


  Assembly(Vec_up_0);
  VecCopy(Vec_up_0,Vec_up_1);

  // remember: Vec_normals follows the Vec_x_1

  //moveMesh(Vec_x_0, Vec_up_0, Vec_up_1, 1.0, tt, Vec_x_1);
  calcMeshVelocity(Vec_x_0, Vec_up_0, Vec_up_1, 1.0, Vec_v_mid, 0.0);
  // move the mesh
  VecWAXPY(Vec_x_1, dt, Vec_v_mid, Vec_x_0); // Vec_x_1 = Vec_v_mid*dt + Vec_x_0

  //point = mesh->pointBegin();
  //point_end = mesh->pointEnd();
  //for (; point != point_end; ++point)
  //{
  ////point->getCoord(X.data());
  ////cout << "HAAAA : " << X.transpose() << endl;
  //} // end point loop

  if (ale)
  {
    setUPInitialGuess();
    printf("Initial conditions:\n");
    for (int i = 0; i < 5; ++i)
    {
      printf("\tIterations %d\n", i);
      // * SOLVE THE SYSTEM *
      if (solve_the_sys)
      {
        //setUPInitialGuess();
        ierr = SNESSolve(snes,PETSC_NULL,Vec_up_1);  CHKERRQ(ierr);
      }
      // * SOLVE THE SYSTEM *

      // update
      if (ale)
      {
        //double tt = time_step==0? dt : current_time;
        //calcMeshVelocity(Vec_x_0, Vec_up_0, Vec_up_1, 1.0, Vec_v_mid, 0.0); // Euler (tem que ser esse no começo)
        calcMeshVelocity(Vec_x_0, Vec_up_0, Vec_up_1, 0.5, Vec_v_mid, 0.0); // Adams-Bashforth
        // move the mesh
        VecWAXPY(Vec_x_1, dt, Vec_v_mid, Vec_x_0); // Vec_x_1 = Vec_v_mid*dt + Vec_x_0

        //compute normal for the next time step, at n+1/2
        {
          Vec Vec_x_mid;
          int xsize;
          double *xarray;
          VecGetSize(Vec_x_0, &xsize);
          VecGetArray(Vec_res, &xarray);
          //prototipo no petsc-dev: VecCreateSeqWithArray(MPI_Comm comm,PetscInt bs,PetscInt n,const PetscScalar array[],Vec *V)
          //VecCreateSeqWithArray(MPI_COMM_SELF, xsize, xarray, &Vec_x_mid);
          VecCreateSeqWithArray(MPI_COMM_SELF, 1, xsize, xarray, &Vec_x_mid);
          Assembly(Vec_x_mid);

          VecCopy(Vec_x_0, Vec_x_mid);
          VecAXPY(Vec_x_mid,1.,Vec_x_1);
          VecScale(Vec_x_mid, 0.5);
          getVecNormals(&Vec_x_mid, Vec_normal);

          VecDestroy(&Vec_x_mid);
          VecRestoreArray(Vec_res, &xarray);
        }

        //// initial guess for the next time step; u(n+1) = 2*u(n) - u(n-1)
        //VecCopy(Vec_up_1, Vec_res);
        //VecScale(Vec_res,2.);
        //VecAXPY(Vec_res, -1., Vec_up_0);
        //VecCopy(Vec_res, Vec_up_1); // u(n+1) = 2*u(n) - u(n-1)

      }
    }
  }
  // normals
  getVecNormals(&Vec_x_1, Vec_normal);

  PetscFunctionReturn(0);
}

PetscErrorCode AppCtx::checkSnesConvergence(SNES snes, PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason)
{
  PetscErrorCode ierr;

  ierr = SNESConvergedDefault(snes,it,xnorm,pnorm,fnorm,reason,NULL); CHKERRQ(ierr);

  // se não convergiu, não terminou ou não é um método ale, retorna
  if (*reason<=0 || !ale)
  {
    return ierr;
  }
  else
  {
    // se não é a primeira vez que converge
    if (converged_times)
    {
      return ierr;
    }
    else
    {

      //Vec *Vec_up_k = &Vec_up_1;
      ////SNESGetSolution(snes,Vec_up_k);
      //
      //copyMesh2Vec(Vec_x_1);
      //calcMeshVelocity(*Vec_up_k, Vec_x_1, Vec_v_mid, current_time+dt);
      //moveMesh(Vec_x_1, Vec_vmsh_0, Vec_v_mid, 0.5);
      //
      //// mean mesh velocity
      //VecAXPY(Vec_v_mid,1,Vec_vmsh_0);
      //VecScale(Vec_v_mid, 0.5);
      //
      //*reason = SNES_CONVERGED_ITERATING;
      //++converged_times;

      return ierr;
    }
  }


/*
  converged
  SNES_CONVERGED_FNORM_ABS         =  2,  ||F|| < atol
  SNES_CONVERGED_FNORM_RELATIVE    =  3,  ||F|| < rtol*||F_initial||
  SNES_CONVERGED_PNORM_RELATIVE    =  4,  Newton computed step size small; || delta x || < stol
  SNES_CONVERGED_ITS               =  5,  maximum iterations reached
  SNES_CONVERGED_TR_DELTA          =  7,
   diverged
  SNES_DIVERGED_FUNCTION_DOMAIN    = -1,  the new x location passed the function is not in the domain of F
  SNES_DIVERGED_FUNCTION_COUNT     = -2,
  SNES_DIVERGED_LINEAR_SOLVE       = -3,  the linear solve failed
  SNES_DIVERGED_FNORM_NAN          = -4,
  SNES_DIVERGED_MAX_IT             = -5,
  SNES_DIVERGED_LINE_SEARCH        = -6,  the line search failed
  SNES_DIVERGED_LOCAL_MIN          = -8,  || J^T b || is small, implies converged to local minimum of F()
  SNES_CONVERGED_ITERATING         =  0
*/
}

PetscErrorCode AppCtx::setUPInitialGuess()
{
  // set u^n+1 b.c.

  VectorXi    u_dofs(dim);
  VectorXi    x_dofs(dim);
  int tag;

  Vector      X1(dim);
  Vector      U1(dim);

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    tag = point->getTag();

    getNodeDofs(&*point, DH_MESH, VAR_M, x_dofs.data());
    getNodeDofs(&*point, DH_UNKS, VAR_U, u_dofs.data());

    VecGetValues(Vec_x_1, dim, x_dofs.data(), X1.data());

    if (  is_in(tag, dirichlet_tags)  )
    {
      U1 = u_exact(X1, current_time+dt, tag);
      VecSetValues(Vec_up_1, dim, u_dofs.data(), U1.data(), INSERT_VALUES);
    }
    else
    if (is_in(tag, solid_tags) || is_in(tag, feature_tags) || is_in(tag, triple_tags))
    {
      U1.setZero();
      VecSetValues(Vec_up_1, dim, u_dofs.data(), U1.data(), INSERT_VALUES);
    }

  } // end for point

  Assembly(Vec_up_1);

  PetscFunctionReturn(0);
}

PetscErrorCode AppCtx::solveTimeProblem()
{
  PetscErrorCode      ierr(0);
  //cout << endl;
  //cout << "Convergence Reason: " << reason << endl;
  //cout << "Number of Newton iterations = " << its << endl;
  //cout << "number of linear iterations used by the nonlinear solver " << lits << endl;

  double initial_volume = getMeshVolume(), final_volume;
  Vector X(dim), Xe(dim);
  double x_error=0;

  if (print_to_matlab)
    printMatlabLoader();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Solve nonlinear system
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  printf("initial volume: %.15lf \n", initial_volume);
  current_time = 0;
  time_step = 0;
  double Qmax=0;
  double steady_error=1;
  setInitialConditions();
  int its;

  printf("Num of time iterations (maxts): %d\n",maxts);
  printf("starting time loop . . . \n");
  
  // BDF2
  if (is_bdf2)
  {
    current_time += dt;
    time_step += 1;
    VecCopy(Vec_up_1, Vec_dup);
    VecAXPY(Vec_dup,-1.0,Vec_up_0); // Vec_dup -= Vec_up_0
    VecScale(Vec_dup, 1./dt);

    copyVec2Mesh(Vec_x_1);

    VecScale(Vec_x_1, 2.0);
    VecAXPY(Vec_x_1,-1.0,Vec_x_0);
    copyMesh2Vec(Vec_x_0);

    calcMeshVelocity(Vec_x_0, Vec_up_0, Vec_up_1, 2.0, Vec_v_1, current_time);

    //VecWAXPY(Vec_x_1, dt, Vec_v_mid, Vec_x_0); // Vec_x_1 = Vec_v_mid*dt + Vec_x_0
    // Vec_x_1 = Vec_x_0 + (2/3 Vec_v_1  + 1/3 Vec_v_mid) * dt
    {
      VecCopy(Vec_v_1, Vec_x_1);
      VecScale(Vec_x_1, 2./3.);
      VecAXPY(Vec_x_1, 1./3.,Vec_v_mid);
      VecScale(Vec_x_1, dt);
      VecAXPY(Vec_x_1,1.,Vec_x_0);
    }
    //{
    //  VecCopy(Vec_v_1, Vec_x_1);
    //  VecScale(Vec_x_1, 2./2.);
    //  VecAXPY(Vec_x_1, 0./2.,Vec_v_mid);
    //  VecScale(Vec_x_1, dt);
    //  VecAXPY(Vec_x_1,1.,Vec_x_0);
    //}


    VecCopy(Vec_up_1, Vec_up_0);

  }

  for(;;)
  {


    cout << endl;
    cout << "current time: " << current_time << endl;
    cout << "time step: "    << time_step  << endl;
    cout << "steady error: " << steady_error << endl;

    //mesh->getNodePtr(0)->getCoord(X.data(), dim);
    //Xe(0) = 1.0*exp(sin(pi*current_time)/pi);
    //Xe(1) = 0;
    //x_error += (X-Xe).norm()/maxts;


    if (maxts == 0)
    {
      computeError(Vec_x_0, Vec_up_0,current_time);
      break;
    }
    
    bool const full_implicit = false;
    bool const try2 = false; // extrapolate geometry Vec_x_1 <- 2*Vec_x_1 - Vec_x_0
    
    // BDF2
    
    
    // * SOLVE THE SYSTEM *
    if (solve_the_sys)
    {
      setUPInitialGuess();
      
      for (int kk = 0 ; kk < 1+3*full_implicit; kk++)
      {
        printf("\tIterations %d\n", kk);
        
        ierr = SNESSolve(snes,PETSC_NULL,Vec_up_1);        CHKERRQ(ierr);
        
        // update
        if (full_implicit)
        {
          //double tt = time_step==0? dt : current_time;
          //calcMeshVelocity(Vec_x_0, Vec_up_0, Vec_up_1, 1.0, Vec_v_mid, 0.0); // Euler (tem que ser esse no começo)
          calcMeshVelocity(Vec_x_0, Vec_up_0, Vec_up_1, 0.5, Vec_v_mid, current_time); // Adams-Bashforth
          // move the mesh
          VecWAXPY(Vec_x_1, dt, Vec_v_mid, Vec_x_0); // Vec_x_1 = Vec_v_mid*dt + Vec_x_0

          //compute normal for the next time step, at n+1/2
          {
            Vec Vec_x_mid;
            int xsize;
            double *xarray;
            VecGetSize(Vec_x_0, &xsize);
            VecGetArray(Vec_res, &xarray);
            //prototipo no petsc-dev: VecCreateSeqWithArray(MPI_Comm comm,PetscInt bs,PetscInt n,const PetscScalar array[],Vec *V)
            //VecCreateSeqWithArray(MPI_COMM_SELF, xsize, xarray, &Vec_x_mid);
            VecCreateSeqWithArray(MPI_COMM_SELF, 1, xsize, xarray, &Vec_x_mid);
            Assembly(Vec_x_mid);

            VecCopy(Vec_x_0, Vec_x_mid);
            VecAXPY(Vec_x_mid,1.,Vec_x_1);
            VecScale(Vec_x_mid, 0.5);
            getVecNormals(&Vec_x_mid, Vec_normal);

            VecDestroy(&Vec_x_mid);
            VecRestoreArray(Vec_res, &xarray);
          }

          //// initial guess for the next time step; u(n+1) = 2*u(n) - u(n-1)
          //VecCopy(Vec_up_1, Vec_res);
          //VecScale(Vec_res,2.);
          //VecAXPY(Vec_res, -1., Vec_up_0);
          //VecCopy(Vec_res, Vec_up_1); // u(n+1) = 2*u(n) - u(n-1)

        }

        
        
      }
    }
    if (is_bdf2)
    {
      if (plot_exact_sol)
        computeError(Vec_x_0, Vec_up_0,current_time);
        VecCopy(Vec_up_1, Vec_dup);
        VecAXPY(Vec_dup,-1.0,Vec_up_0); // Vec_dup -= Vec_up_0
        VecScale(Vec_dup, 1./dt);
    }
    else
    {
      if (time_step == 0)
      {
        pressureTimeCorrection(Vec_up_0, Vec_up_1, 0., 1); // press(n) = press(n+1/2) - press(n-1/2)
        if (plot_exact_sol && maxts <= 1)
          computeError(Vec_x_0, Vec_up_0,current_time);
      }
      else
      {
        pressureTimeCorrection(Vec_up_0, Vec_up_1, .5, .5); // press(n) = press(n+1/2) - press(n-1/2)
        if (plot_exact_sol)
          computeError(Vec_x_0, Vec_up_0,current_time);
      }
    }
    
    bool must_print = false;
    if (is_bdf2)
    {
      if (((time_step-1)%print_step)==0 || time_step == (maxts-1))
        must_print = true;
    }
    else
      if ((time_step%print_step)==0 || time_step == (maxts-1))
        must_print = true;
        
    if (must_print)
    {
      if (family_files)
      {
        double  *q_array;
        double  *nml_array;
        double  *v_array;
        VecGetArray(Vec_up_0, &q_array);
        VecGetArray(Vec_normal, &nml_array);
        VecGetArray(Vec_v_mid, &v_array);
        vtk_printer.writeVtk();
        
        /* ---- nodes data ---- */
        vtk_printer.addNodeVectorVtk("u", GetDataVelocity(q_array, *this));
        //vtk_printer.addNodeVectorVtk("normal",  GetDataNormal(nml_array, *this));
        //vtk_printer.addNodeVectorVtk("v",  GetDataMeshVel(v_array, *this));
        vtk_printer.printPointTagVtk();
        
        if (!shape_psi_c->discontinuous())
          vtk_printer.addNodeScalarVtk("pressure", GetDataPressure(q_array, *this));
        else
          vtk_printer.addCellScalarVtk("pressure", GetDataPressCellVersion(q_array, *this));

        
        vtk_printer.addCellIntVtk("cell_tag", GetDataCellTag(*this));        
        

        //vtk_printer.printPointTagVtk("point_tag");
        VecRestoreArray(Vec_up_0, &q_array);
        VecRestoreArray(Vec_normal, &nml_array);
        VecRestoreArray(Vec_v_mid, &v_array);

        ierr = SNESGetIterationNumber(snes,&its);     CHKERRQ(ierr);
        cout << "num snes iterations: " << its << endl;
      }

    }

    printContactAngle(fprint_ca);
    //////if (plot_exact_sol) eh melhor depois da atualização
    //////  computeError(Vec_x_1, Vec_up_1,current_time+dt);
    current_time += dt;
    time_step += 1;

    // update
    if (ale)
    {
      //double tt = time_step==0? dt : current_time;

      // mesh adaption (it has topological changes) (2d only)
      // it destroys Vec_normal and Vec_v_mid
      //if(time_step == 1 || time_step == 3)
      if (mesh_adapt)
        meshAdapt();

      copyVec2Mesh(Vec_x_1);

      // MESH FLIPPING
      if (mesh_adapt)
      {
        if(time_step%1 == 0)
        {
          if (dim==2)
          {
            Facet *f;
            Cell* ca, *cb;
            // Delaunay
            for (int i = 0; i < n_facets_total; ++i)
            {
              f = mesh->getFacetPtr(i);
              if (f->isDisabled())
                continue;
              if (!MeshToolsTri::inCircle2d(f, &*mesh))
              {
                ca = mesh->getCellPtr(f->getIncidCell());
                cb = mesh->getCellPtr(ca->getIncidCell(f->getPosition()));

                if (cb)
                  if (ca->getTag() == cb->getTag())
                    MeshToolsTri::flipEdge(f, &*mesh, true);
              }
            }

          } // end dim==2
        } // end ts%1
      } // end mesh adapt

      if (full_implicit)
        VecCopy(Vec_x_1, Vec_x_0);
      else if (try2)
      {
        VecWAXPY(Vec_x_1, dt, Vec_v_mid, Vec_x_0);
        copyVec2Mesh(Vec_x_1);
        VecAXPBY(Vec_x_1,-1.0, 2.0,Vec_x_0); // Vec_x_1 <- 2*Vec_x_1 - Vec_x_0
        copyMesh2Vec(Vec_x_0);
      }
      else
      {
        VecScale(Vec_x_1, 2.0);
        VecAXPY(Vec_x_1,-1.0,Vec_x_0);
        copyMesh2Vec(Vec_x_0);
      }



      /* Compute Vec_x_1 and Vec_v_mid */
      /////////moveMesh(Vec_x_0, Vec_up_0, Vec_up_1, 1.0, current_time, Vec_x_1); // Euler
      ///////moveMesh(Vec_x_0, Vec_up_0, Vec_up_1, 1.5, current_time, Vec_x_1); // Adams-Bashforth
      /////////moveMesh(Vec_x_0, Vec_up_0, Vec_up_1, 0.5, current_time, Vec_x_1); // Alguma-coisa
      ///////calcMeshVelocity(Vec_x_0, Vec_up_0, Vec_up_1, 1.5, Vec_v_mid, current_time);
      if ( !full_implicit )
      {
        if (is_bdf2)
        {
          VecCopy(Vec_v_1, Vec_v_mid);
          calcMeshVelocity(Vec_x_0, Vec_up_0, Vec_up_1, 2.0, Vec_v_1, current_time);
          VecAXPBY(Vec_v_mid, .5, .5, Vec_v_1);
        }
        else
          calcMeshVelocity(Vec_x_0, Vec_up_0, Vec_up_1, 1.5, Vec_v_mid, current_time); // Adams-Bashforth
        //calcMeshVelocity(Vec_x_0, Vec_up_0, Vec_up_1, 1.0, Vec_v_mid, 0.0); // Euler
        // move the mesh
        if (!try2)
          VecWAXPY(Vec_x_1, dt, Vec_v_mid, Vec_x_0); // Vec_x_1 = Vec_v_mid*dt + Vec_x_0
      }


      //compute normal for the next time step, at n+1/2
      {
        Vec Vec_x_mid;
        int xsize;
        double *xarray;
        VecGetSize(Vec_x_0, &xsize);
        VecGetArray(Vec_res, &xarray);
        //prototipo no petsc-dev: VecCreateSeqWithArray(MPI_Comm comm,PetscInt bs,PetscInt n,const PetscScalar array[],Vec *V)
        //VecCreateSeqWithArray(MPI_COMM_SELF, xsize, xarray, &Vec_x_mid);
        VecCreateSeqWithArray(MPI_COMM_SELF, 1,xsize, xarray, &Vec_x_mid);
        Assembly(Vec_x_mid);

        VecCopy(Vec_x_0, Vec_x_mid);
        VecAXPY(Vec_x_mid,1.,Vec_x_1);
        VecScale(Vec_x_mid, 0.5);
        getVecNormals(&Vec_x_mid, Vec_normal);

        VecDestroy(&Vec_x_mid);
        VecRestoreArray(Vec_res, &xarray);
      }

      // initial guess for the next time step; u(n+1) = 2*u(n) - u(n-1)
      //VecCopy(Vec_up_1, Vec_res);
      //VecScale(Vec_res,2.);
      //VecAXPY(Vec_res, -1., Vec_up_0);
      VecCopy(Vec_up_1, Vec_up_0);
      //VecCopy(Vec_res, Vec_up_1); // u(n+1) = 2*u(n) - u(n-1)

    }
    else
    {
      VecCopy(Vec_up_1, Vec_up_0);
    }


    // compute steady error
    VecNorm(Vec_up_1, NORM_1, &Qmax);
    VecCopy(Vec_up_0, Vec_res);
    VecAXPY(Vec_res,-1.0,Vec_up_1);
    //steady_error = VecNorm(Vec_res, NORM_1)/(Qmax==0.?1.:Qmax);
    VecNorm(Vec_res, NORM_1, &steady_error);
    steady_error /= (Qmax==0.?1.:Qmax);



    if(time_step >= maxts) {
      cout << "\n==========================================\n";
      cout << "stop reason:\n";
      cout << "maximum number of iterations reached. \n";
      break;
    }
    if (steady_error <= steady_tol) {
      cout << "\n==========================================\n";
      cout << "stop reason:\n";
      cout << "steady state reached. \n";
      break;
    }
  }
  //
  // END TIME LOOP
  //

  // PRINT AGAIN
  if (false)
  {
    if (family_files)
    {
      double  *q_array;
      double  *nml_array;
      double  *v_array;
      VecGetArray(Vec_up_0, &q_array);
      VecGetArray(Vec_normal, &nml_array);
      VecGetArray(Vec_v_mid, &v_array);
      vtk_printer.writeVtk();
      
      /* ---- nodes data ---- */
      vtk_printer.addNodeVectorVtk("u", GetDataVelocity(q_array, *this));
      //vtk_printer.addNodeVectorVtk("normal",  GetDataNormal(nml_array, *this));
      //vtk_printer.addNodeVectorVtk("v",  GetDataMeshVel(v_array, *this));
      vtk_printer.printPointTagVtk();
      
      if (!shape_psi_c->discontinuous())
        vtk_printer.addNodeScalarVtk("pressure", GetDataPressure(q_array, *this));
      else
        vtk_printer.addCellScalarVtk("pressure", GetDataPressCellVersion(q_array, *this));

      
      vtk_printer.addCellIntVtk("cell_tag", GetDataCellTag(*this));        
      

      //vtk_printer.printPointTagVtk("point_tag");
      VecRestoreArray(Vec_up_0, &q_array);
      VecRestoreArray(Vec_normal, &nml_array);
      VecRestoreArray(Vec_v_mid, &v_array);

      ierr = SNESGetIterationNumber(snes,&its);     CHKERRQ(ierr);
      cout << "num snes iterations: " << its << endl;
    }

  }


  cout << endl;

  final_volume = getMeshVolume();
  printf("final volume: %.15lf \n", final_volume);
  printf("volume error 100*(f-i)/i: %.15lf per percent\n", 100*abs(final_volume-initial_volume)/initial_volume);
  printf("x error : %.15lf \n", x_error);


  SNESConvergedReason reason;
  ierr = SNESGetIterationNumber(snes,&its);     CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,&reason);  CHKERRQ(ierr);
  cout << "num snes iterations: " << its << endl;
  cout << "reason: " << reason << endl;
  switch (reason)
  {
  case( SNES_CONVERGED_FNORM_ABS       ): printf("SNES_CONVERGED_FNORM_ABS      /* ||F|| < atol */\n"); break;
  case( SNES_CONVERGED_FNORM_RELATIVE  ): printf("SNES_CONVERGED_FNORM_RELATIVE /* ||F|| < rtol*||F_initial|| */\n"); break;
  case( SNES_CONVERGED_SNORM_RELATIVE  ): printf("SNES_CONVERGED_PNORM_RELATIVE /* Newton computed step size small); break; || delta x || < stol */\n"); break;
  case( SNES_CONVERGED_ITS             ): printf("SNES_CONVERGED_ITS            /* maximum iterations reached */\n"); break;
  case( SNES_CONVERGED_TR_DELTA        ): printf("SNES_CONVERGED_TR_DELTA       \n"); break;
        /* diverged */
  case( SNES_DIVERGED_FUNCTION_DOMAIN  ): printf("SNES_DIVERGED_FUNCTION_DOMAIN /* the new x location passed the function is not in the domain of F */\n"); break;
  case( SNES_DIVERGED_FUNCTION_COUNT   ): printf("SNES_DIVERGED_FUNCTION_COUNT   \n"); break;
  case( SNES_DIVERGED_LINEAR_SOLVE     ): printf("SNES_DIVERGED_LINEAR_SOLVE    /* the linear solve failed */\n"); break;
  case( SNES_DIVERGED_FNORM_NAN        ): printf("SNES_DIVERGED_FNORM_NAN       \n"); break;
  case( SNES_DIVERGED_MAX_IT           ): printf("SNES_DIVERGED_MAX_IT          \n"); break;
  case( SNES_DIVERGED_LINE_SEARCH      ): printf("SNES_DIVERGED_LINE_SEARCH     /* the line search failed */ \n"); break;
  case( SNES_DIVERGED_LOCAL_MIN        ): printf("SNES_DIVERGED_LOCAL_MIN       /* || J^T b || is small, implies converged to local minimum of F() */\n"); break;
  case( SNES_CONVERGED_ITERATING       ): printf("SNES_CONVERGED_ITERATING      \n"); break;
  case( SNES_DIVERGED_INNER            ): printf("SNES_DIVERGED_INNER           \n"); break;
  }

  if (solve_the_sys)
    MatrixInfo(Mat_Jac);

  int lits;
  SNESGetLinearSolveIterations(snes,&lits);

  cout << "Greatest error reached during the simulation:" << endl;
  printf("%-21s %-21s %-21s %-21s %-21s %-21s %-21s %-21s %s\n", "# hmean", "u_L2_norm", "p_L2_norm", "grad_u_L2_norm", "grad_p_L2_norm", "u_L2_facet_norm", "u_inf_facet_norm", "u_inf_norm", "p_inf_norm" );
  printf("%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n\n",Stats.hmean, Stats.u_L2_norm, Stats.p_L2_norm, Stats.grad_u_L2_norm, Stats.grad_p_L2_norm, Stats.u_L2_facet_norm,  Stats.u_inf_facet_norm, Stats.u_inf_norm, Stats.p_inf_norm);

  ////if (unsteady)
  //{
  //  cout << "\nmean errors: \n";
  //  cout << "# hmean            u_L2_err         p_L2_err          grad_u_L2_err     grad_p_L2_err" << endl;
  //  printf( "%.15lf %.15lf %.15lf %.15lf %.15lf\n", Stats.mean_hmean() , Stats.mean_u_L2_norm(), Stats.mean_p_L2_norm(), Stats.mean_grad_u_L2_norm(), Stats.mean_grad_p_L2_norm());
  //}


  PetscFunctionReturn(0);
}

void AppCtx::computeError(Vec const& Vec_x, Vec &Vec_up, double tt)
{

  MatrixXd            u_coefs_c(n_dofs_u_per_cell/dim, dim);
  MatrixXd            u_coefs_c_trans(dim,n_dofs_u_per_cell/dim);
  MatrixXd            u_coefs_f(n_dofs_u_per_facet/dim, dim);
  MatrixXd            u_coefs_f_trans(dim,n_dofs_u_per_facet/dim);
  VectorXd            p_coefs_c(n_dofs_p_per_cell);
  MatrixXd            x_coefs_c(nodes_per_cell, dim);
  MatrixXd            x_coefs_c_trans(dim, nodes_per_cell);
  MatrixXd            x_coefs_f(nodes_per_facet, dim);
  MatrixXd            x_coefs_f_trans(dim, nodes_per_facet);
  Tensor              F_c(dim,dim), invF_c(dim,dim), invFT_c(dim,dim);
  Tensor              F_f(dim,dim-1), invF_f(dim-1,dim), fff(dim-1,dim-1);
  MatrixXd            dxphi_err(n_dofs_u_per_cell/dim, dim);
  MatrixXd            dxpsi_err(n_dofs_p_per_cell, dim);
  MatrixXd            dxqsi_err(nodes_per_cell, dim);
  Tensor              dxU(dim,dim); // grad u
  Vector              dxP(dim);     // grad p
  Vector              Xqp(dim);
  Vector              Uqp(dim);

  double              Pqp;
  VectorXi            cell_nodes(nodes_per_cell);
  double              J, JxW;
  double              weight;
  int                 tag;
  double              volume=0;


  double              p_L2_norm = 0.;
  double              u_L2_norm = 0.;
  double              grad_u_L2_norm = 0.;
  double              grad_p_L2_norm = 0.;
  double              p_inf_norm = 0.;
  double              u_inf_norm = 0.;
  double              hmean = 0.;
  double              u_L2_facet_norm = 0.;
  double              u_inf_facet_norm = 0.;


  VectorXi            mapU_c(n_dofs_u_per_cell);
  VectorXi            mapU_f(n_dofs_u_per_facet);
  VectorXi            mapU_r(n_dofs_u_per_corner);
  VectorXi            mapP_c(n_dofs_p_per_cell);
  VectorXi            mapP_r(n_dofs_p_per_corner);
  VectorXi            mapM_c(dim*nodes_per_cell);
  VectorXi            mapM_f(dim*nodes_per_facet);
  VectorXi            mapM_r(dim*nodes_per_corner);


  cell_iterator cell = mesh->cellBegin();
  cell_iterator cell_end = mesh->cellEnd();
  for (; cell != cell_end; ++cell)
  {
    tag = cell->getTag();

    // mapeamento do local para o global:
    //
    dof_handler[DH_UNKS].getVariable(VAR_U).getCellDofs(mapU_c.data(), &*cell);
    dof_handler[DH_UNKS].getVariable(VAR_P).getCellDofs(mapP_c.data(), &*cell);
    dof_handler[DH_MESH].getVariable(VAR_M).getCellDofs(mapM_c.data(), &*cell);

    /*  Pega os valores das variáveis nos graus de liberdade */
    VecGetValues(Vec_up, mapU_c.size(), mapU_c.data(), u_coefs_c.data());
    VecGetValues(Vec_up, mapP_c.size(), mapP_c.data(), p_coefs_c.data());
    VecGetValues(Vec_x,  mapM_c.size(), mapM_c.data(), x_coefs_c.data());

    u_coefs_c_trans = u_coefs_c.transpose();

    //mesh->getCellNodesId(&*cell, cell_nodes.data());
    //mesh->getNodesCoords(cell_nodes.begin(), cell_nodes.end(), x_coefs_c.data());
    x_coefs_c_trans = x_coefs_c.transpose();

    for (int qp = 0; qp < n_qpts_err; ++qp)
    {
      F_c    = x_coefs_c_trans * dLqsi_err[qp];
      J      = F_c.determinant();
      invF_c = F_c.inverse();
      invFT_c= invF_c.transpose();

      dxphi_err = dLphi_err[qp] * invF_c;
      dxpsi_err = dLpsi_err[qp] * invF_c;
      dxqsi_err = dLqsi_err[qp] * invF_c;

      dxP  = dxpsi_err.transpose() * p_coefs_c;
      dxU  = u_coefs_c_trans * dxphi_err;       // n+utheta

      Xqp  = x_coefs_c_trans * qsi_err[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
      Uqp  = u_coefs_c_trans * phi_err[qp];
      Pqp  = p_coefs_c.dot(psi_err[qp]);

      weight = quadr_err->weight(qp);

      JxW = J*weight;

      //  note: norm(u, H1)^2 = norm(u, L2)^2 + norm(gradLphi_c, L2)^2
      u_L2_norm        += (u_exact(Xqp, tt, tag) - Uqp).squaredNorm()*JxW;
      p_L2_norm        += sqr(pressure_exact(Xqp, tt, tag) - Pqp)*JxW;
      grad_u_L2_norm   += (grad_u_exact(Xqp, tt, tag) - dxU).squaredNorm()*JxW;
      grad_p_L2_norm   += (grad_p_exact(Xqp, tt, tag) - dxP).squaredNorm()*JxW;
      u_inf_norm       = max(u_inf_norm, (u_exact(Xqp, tt, tag) - Uqp).norm());
      p_inf_norm       = max(p_inf_norm, fabs(pressure_exact(Xqp, tt, tag) - Pqp));

      volume += JxW;
    } // fim quadratura

  } // end elementos

  facet_iterator facet = mesh->facetBegin();
  facet_iterator facet_end = mesh->facetEnd();
  for (; facet != facet_end; ++facet)
  {
    tag = facet->getTag();

    if (!(is_in(tag, dirichlet_tags) || is_in(tag, neumann_tags) || is_in(tag, interface_tags) ||
          is_in(tag, solid_tags) || is_in(tag, periodic_tags) || is_in(tag, feature_tags)))
      continue;

    // mapeamento do local para o global:
    //
    dof_handler[DH_UNKS].getVariable(VAR_U).getFacetDofs(mapU_f.data(), &*facet);
    dof_handler[DH_MESH].getVariable(VAR_M).getFacetDofs(mapM_f.data(), &*facet);

    /*  Pega os valores das variáveis nos graus de liberdade */
    VecGetValues(Vec_up, mapU_f.size(), mapU_f.data(), u_coefs_f.data());
    VecGetValues(Vec_x,  mapM_f.size(), mapM_f.data(), x_coefs_f.data());

    u_coefs_f_trans = u_coefs_f.transpose();
    x_coefs_f_trans = x_coefs_f.transpose();

    for (int qp = 0; qp < n_qpts_facet; ++qp)
    {
      F_f   = x_coefs_f_trans * dLqsi_f[qp];

      fff    = F_f.transpose()*F_f;
      J      = sqrt(fff.determinant());

      Xqp  = x_coefs_f_trans * qsi_f[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
      Uqp  = u_coefs_f_trans * phi_f[qp];

      weight = quadr_facet->weight(qp);

      JxW = J*weight;

      Vector U_exact = u_exact(Xqp, tt, tag);

      u_L2_facet_norm += (U_exact - Uqp).squaredNorm()*JxW;
      //u_L2_facet_norm += JxW;

      double const diff = (U_exact - Uqp).norm();

      if (diff > u_inf_facet_norm)
        u_inf_facet_norm = diff;


    } // fim quadratura

  }


  u_L2_norm      = sqrt(u_L2_norm     );
  p_L2_norm      = sqrt(p_L2_norm     );
  grad_u_L2_norm = sqrt(grad_u_L2_norm);
  grad_p_L2_norm = sqrt(grad_p_L2_norm);
  u_L2_facet_norm = sqrt(u_L2_facet_norm);

  // ASSUME QUE SÓ POSSA TER NO MÁXIMO 1 NÓ POR ARESTA


  VectorXi edge_nodes(3);
  Vector Xa(dim), Xb(dim);
  int n_edges=0;

  if (dim==2)
  //FEP_PRAGMA_OMP(parallel default(none) shared(hmean))
  {
    const int n_edges_total = mesh->numFacetsTotal();
    Facet const* edge(NULL);

    //FEP_PRAGMA_OMP(for nowait)
    for (int a = 0; a < n_edges_total; ++a)
    {
      edge = mesh->getFacetPtr(a);
      if (edge->isDisabled())
        continue;

      mesh->getFacetNodesId(&*edge, edge_nodes.data());

      mesh->getNodePtr(edge_nodes[0])->getCoord(Xa.data(),dim);
      mesh->getNodePtr(edge_nodes[1])->getCoord(Xb.data(),dim);
      hmean += (Xa-Xb).norm();
      ++n_edges;
    }
  }
  else
  if (dim==3)
  //FEP_PRAGMA_OMP(parallel default(none) shared(cout,hmean))
  {
    const int n_edges_total = mesh->numCornersTotal();
    Corner const* edge(NULL);

    //FEP_PRAGMA_OMP(for nowait)
    for (int a = 0; a < n_edges_total; ++a)
    {
      edge = mesh->getCornerPtr(a);
      if (edge->isDisabled())
        continue;

      mesh->getCornerNodesId(&*edge, edge_nodes.data());

      mesh->getNodePtr(edge_nodes[0])->getCoord(Xa.data(),dim);
      mesh->getNodePtr(edge_nodes[1])->getCoord(Xb.data(),dim);
      hmean += (Xa-Xb).norm();
      ++n_edges;
    }
  }
  hmean /= n_edges;

  //if (time_step==maxts)
  {
    cout << endl;
    //cout << "errors computed at last time step: "<< endl;
    //cout << "# hmean               u_L2_norm         p_L2_norm         grad_u_L2_norm    grad_p_L2_norm   u_L2_facet_norm    u_inf_facet_norm" << endl;
    printf("%-21s %-21s %-21s %-21s %-21s %-21s %s\n", "# hmean", "u_L2_norm", "p_L2_norm", "grad_u_L2_norm", "grad_p_L2_norm", "u_L2_facet_norm", "u_inf_facet_norm" );
    printf("%.15e %.15e %.15e %.15e %.15e %.15e %.15e\n\n",hmean, u_L2_norm, p_L2_norm, grad_u_L2_norm, grad_p_L2_norm, u_L2_facet_norm,  u_inf_facet_norm);
  }

  Stats.add_p_L2_norm        (p_L2_norm       );
  Stats.add_u_L2_norm        (u_L2_norm       );
  Stats.add_grad_u_L2_norm   (grad_u_L2_norm  );
  Stats.add_grad_p_L2_norm   (grad_p_L2_norm  );
  Stats.add_p_inf_norm       (p_inf_norm      );
  Stats.add_u_inf_norm       (u_inf_norm      );
  Stats.add_hmean            (hmean           );
  Stats.add_u_L2_facet_norm  (u_L2_facet_norm );
  Stats.add_u_inf_facet_norm (u_inf_facet_norm);


}


void AppCtx::pressureTimeCorrection(Vec &Vec_up_0, Vec &Vec_up_1, double a, double b) // p(n+1) = a*p(n+.5) + b* p(n)
{
  Vector    Uf(dim);
  //double    pf;
  Vector    X(dim);
  Tensor    R(dim,dim);
  int       dof;
  //int       dof;
  int tag;
  double P0, P1, P2;

  if (behaviors & BH_Press_grad_elim)
  {
    cout << "FIX ME: NOT SUPPORTED YET!!!!\n";
    throw;
  }

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    // press
    if (mesh->isVertex(&*point))
    {
      getNodeDofs(&*point,DH_UNKS,VAR_P,&dof);

      VecGetValues(Vec_up_1, 1, &dof, &P1);
      VecGetValues(Vec_up_0, 1, &dof, &P0);

      P2 = a*P1 + b*P0;

      VecSetValues(Vec_up_0, 1, &dof, &P2, INSERT_VALUES);
    }

  } // end point loop
  Assembly(Vec_up_0);
}

double AppCtx::getMaxVelocity()
{
  Vector     Uf(dim); // Ue := elastic velocity
  VectorXi   vtx_dofs_fluid(dim); // indices de onde pegar a velocidade
  double     Umax=0;


  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    if (!mesh->isVertex(&*point))
        continue;
    dof_handler[DH_UNKS].getVariable(VAR_U).getVertexDofs(vtx_dofs_fluid.data(), &*point);
    VecGetValues(Vec_up_1, dim, vtx_dofs_fluid.data(), Uf.data());
    Umax = max(Umax,Uf.norm());
  }
  return Umax;
}

double AppCtx::getMeshVolume()
{
  MatrixXd            x_coefs_c(nodes_per_cell, dim);
  MatrixXd            x_coefs_c_trans(dim, nodes_per_cell);
  Tensor              F_c(dim,dim), invF_c(dim,dim), invFT_c(dim,dim);
  VectorXi            cell_nodes(nodes_per_cell);
  double              Jx;
  //int                 tag;
  double              volume=0;

  cell_iterator cell = mesh->cellBegin();
  cell_iterator cell_end = mesh->cellEnd();
  for (; cell != cell_end; ++cell)
  {
    //tag = cell->getTag();

    mesh->getCellNodesId(&*cell, cell_nodes.data());
    mesh->getNodesCoords(cell_nodes.begin(), cell_nodes.end(), x_coefs_c.data());
    x_coefs_c_trans = x_coefs_c.transpose();

    for (int qp = 0; qp < n_qpts_err; ++qp)
    {
      F_c    = x_coefs_c_trans * dLqsi_err[qp];
      Jx     = F_c.determinant();
      volume += Jx * quadr_err->weight(qp);

    } // fim quadratura

  } // end elementos
  return volume;
}

// TODO Only for 2D !!!!!!!!!!!!!!!!!!
void AppCtx::printContactAngle(bool _print)
{
  if (!_print)
    return;

  //at file:
  // current_time theta_max theta_min

  int tag;
  double theta_max=0;
  double theta_min=pi;
  double tet;
  Vector X(dim);
  Vector normal_solid(dim);
  Vector normal_surf(dim);
  VectorXi vtx_dofs_mesh(dim);

  Point* plc = NULL;

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    tag = point->getTag();

    if (!is_in(tag,triple_tags))
      continue;

    point->getCoord(X.data(), dim);
    normal_solid = solid_normal(X, current_time, tag);
    dof_handler[DH_MESH].getVariable(VAR_M).getVertexDofs(vtx_dofs_mesh.data(), &*point);
    VecGetValues(Vec_normal, dim, vtx_dofs_mesh.data(), normal_surf.data());

    tet = asin(sqrt(1. - pow(normal_solid.dot(normal_surf),2)  ));

    if (tet < theta_min)
      theta_min = tet;
    if (tet > theta_max)
      theta_max = tet;

    if (X(0)>0)
      plc = &*point;
  }

  theta_min = theta_min*180./pi;
  theta_max = theta_max*180./pi;

  cout << "theta min: " << theta_min << "\ntheta max: " << theta_max << endl;


  // also print some energies
  double field_energy   = 0;
  double kinetic_energy = 0;
  double euler_dissip   = 0; // ver Gerbeau
  double viscous_power  = 0;
  double surface_energy = 0;
  double solid_power    = 0;
  double cl_power       = 0;
  double volume         = 0;

  // Field energy, Kinetic energy, Viscous power
  // The variables are at time step n by default
  {
    MatrixXd            u_coefs_c(n_dofs_u_per_cell/dim, dim);
    MatrixXd            u_coefs_c_trans(dim,n_dofs_u_per_cell/dim);
    MatrixXd            u_coefs_c_new(n_dofs_u_per_cell/dim, dim);
    MatrixXd            u_coefs_c_t_new(dim,n_dofs_u_per_cell/dim);
    VectorXd            p_coefs_c(n_dofs_p_per_cell);
    MatrixXd            x_coefs_c(nodes_per_cell, dim);
    MatrixXd            x_coefs_c_trans(dim, nodes_per_cell);
    Tensor              F_c(dim,dim), invF_c(dim,dim), invFT_c(dim,dim);
    MatrixXd            dxphi_c(n_dofs_u_per_cell/dim, dim);
    MatrixXd            dxpsi_c(n_dofs_p_per_cell, dim);
    MatrixXd            dxqsi_c(nodes_per_cell, dim);
    Tensor              dxU(dim,dim); // grad u
    Vector              dxP(dim);     // grad p
    Vector              Xqp(dim);
    Vector              Uqp(dim);
    Vector              Uqp_new(dim);

//    double              Pqp;
    VectorXi            cell_nodes(nodes_per_cell);
    double              Jx, JxW;
    double              weight;
    int                 tag;

    VectorXi            mapU_c(n_dofs_u_per_cell);
    VectorXi            mapU_r(n_dofs_u_per_corner);
    VectorXi            mapP_c(n_dofs_p_per_cell);
    VectorXi            mapP_r(n_dofs_p_per_corner);
    VectorXi            mapM_c(dim*nodes_per_cell);
    VectorXi            mapM_r(dim*nodes_per_corner);


    cell_iterator cell = mesh->cellBegin();
    cell_iterator cell_end = mesh->cellEnd();
    for (; cell != cell_end; ++cell)
    {
      tag = cell->getTag();

      // mapeamento do local para o global:
      //
      dof_handler[DH_UNKS].getVariable(VAR_U).getCellDofs(mapU_c.data(), &*cell);
      dof_handler[DH_UNKS].getVariable(VAR_P).getCellDofs(mapP_c.data(), &*cell);
      dof_handler[DH_MESH].getVariable(VAR_M).getCellDofs(mapM_c.data(), &*cell);

      /*  old coefs */
      VecGetValues(Vec_up_0, mapU_c.size(), mapU_c.data(), u_coefs_c.data());
      VecGetValues(Vec_up_0, mapP_c.size(), mapP_c.data(), p_coefs_c.data());
      VecGetValues(Vec_x_0,  mapM_c.size(), mapM_c.data(), x_coefs_c.data());

      // new coefs
      VecGetValues(Vec_up_1, mapU_c.size(), mapU_c.data(), u_coefs_c_new.data());

      u_coefs_c_trans = u_coefs_c.transpose();
      u_coefs_c_t_new = u_coefs_c_new.transpose();

      //mesh->getCellNodesId(&*cell, cell_nodes.data());
      //mesh->getNodesCoords(cell_nodes.begin(), cell_nodes.end(), x_coefs_c.data());
      x_coefs_c_trans = x_coefs_c.transpose();

      for (int qp = 0; qp < n_qpts_cell; ++qp)
      {
        F_c    = x_coefs_c_trans * dLqsi_c[qp];
        Jx     = F_c.determinant();
        invF_c = F_c.inverse();
        invFT_c= invF_c.transpose();

        dxphi_c = dLphi_c[qp] * invF_c;
        dxpsi_c = dLpsi_c[qp] * invF_c;
        dxqsi_c = dLqsi_c[qp] * invF_c;

        dxP  = dxpsi_c.transpose() * p_coefs_c;
        dxU  = u_coefs_c_trans * dxphi_c;       // n+utheta

        Xqp  = x_coefs_c_trans * qsi_c[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        Uqp  = u_coefs_c_trans * phi_c[qp];
        Uqp_new = u_coefs_c_t_new * phi_c[qp];
        //Pqp  = p_coefs_c.dot(psi_c[qp]);

        weight = quadr_cell->weight(qp);
        JxW = Jx*weight;

        field_energy   += force(Xqp,current_time,tag).dot(Xqp) * JxW;
        kinetic_energy += 0.5* pho(Xqp,tag) * Uqp.squaredNorm() * JxW;
        euler_dissip   += 0.5* pho(Xqp,tag) * (Uqp-Uqp_new).squaredNorm() * JxW / dt;
        viscous_power  += 0.5* muu(tag)* (dxU + dxU.transpose()).squaredNorm() * JxW;
        volume         += JxW;

      } // fim quadratura

    } // end elementos

  }


  // Surface energy
  {
    MatrixXd            u_coefs_f(n_dofs_u_per_facet/dim, dim);
    MatrixXd            u_coefs_f_trans(dim,n_dofs_u_per_facet/dim);
    VectorXd            p_coefs_f(n_dofs_p_per_facet);
    MatrixXd            x_coefs_f(nodes_per_facet, dim);
    MatrixXd            x_coefs_f_trans(dim, nodes_per_facet);
    Tensor             F_f(dim,dim-1);       //
    Tensor             invF_f(dim-1,dim);    //
    Tensor             fff_f(dim-1,dim-1);   //  fff = first fundamental form
    MatrixXd            dxphi_f(n_dofs_u_per_facet/dim, dim);
    MatrixXd            dxqsi_f(nodes_per_facet, dim);
    Vector             Xqp(dim);
    Vector             Uqp(dim);

    VectorXi            facet_nodes(nodes_per_facet);
    double              Jx, JxW;
    double              weight;
    int                 tag;
    bool                is_neumann;
    bool                is_surface;
    bool                is_solid;

    VectorXi            mapU_f(n_dofs_u_per_facet);
    VectorXi            mapP_f(n_dofs_p_per_facet);
    VectorXi            mapM_f(dim*nodes_per_facet);


    facet_iterator facet = mesh->facetBegin();
    facet_iterator facet_end = mesh->facetEnd();
    for (; facet != facet_end; ++facet)
    {
      tag = facet->getTag();

      is_neumann = (neumann_tags.end() != std::find(neumann_tags.begin(), neumann_tags.end(), tag));
      is_surface = (interface_tags.end() != std::find(interface_tags.begin(), interface_tags.end(), tag));
      is_solid = (solid_tags.end() != std::find(solid_tags.begin(), solid_tags.end(), tag));

      if ((!is_neumann) && (!is_surface) && (!is_solid))
        //PetscFunctionReturn(0);
        continue;

      // mapeamento do local para o global:
      //
      dof_handler[DH_UNKS].getVariable(VAR_U).getFacetDofs(mapU_f.data(), &*facet);
      dof_handler[DH_UNKS].getVariable(VAR_P).getFacetDofs(mapP_f.data(), &*facet);
      dof_handler[DH_MESH].getVariable(VAR_M).getFacetDofs(mapM_f.data(), &*facet);

      /*  Pega os valores das variáveis nos graus de liberdade */
      VecGetValues(Vec_up_0, mapU_f.size(), mapU_f.data(), u_coefs_f.data());
      VecGetValues(Vec_up_0, mapP_f.size(), mapP_f.data(), p_coefs_f.data());
      VecGetValues(Vec_x_0,  mapM_f.size(), mapM_f.data(), x_coefs_f.data());

      u_coefs_f_trans = u_coefs_f.transpose();

      //mesh->getFacetNodesId(&*facet, facet_nodes.data());
      //mesh->getNodesCoords(facet_nodes.begin(), facet_nodes.end(), x_coefs_f.data());
      x_coefs_f_trans = x_coefs_f.transpose();

      for (int qp = 0; qp < n_qpts_facet; ++qp)
      {
        F_f    = x_coefs_f_trans * dLqsi_f[qp];
        fff_f.resize(dim-1,dim-1);
        fff_f  = F_f.transpose()*F_f;
        Jx     = sqrt(fff_f.determinant());
        invF_f = fff_f.inverse()*F_f.transpose();

        dxphi_f = dLphi_f[qp] * invF_f;
        dxqsi_f = dLqsi_f[qp] * invF_f;

        Xqp  = x_coefs_f_trans * qsi_f[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        Uqp  = u_coefs_f_trans * phi_f[qp];

        weight = quadr_facet->weight(qp);
        JxW = Jx*weight;


        //if (is_neumann)
        //{
        //  //Vector no(Xqp);
        //  //no.normalize();
        //  //traction_ = utheta*(traction(Xqp,current_time+dt,tag)) + (1.-utheta)*traction(Xqp,current_time,tag);
        //  traction_ = traction(Xqp, normal, current_time + dt/2,tag);
        //  //traction_ = (traction(Xqp,current_time,tag) +4.*traction(Xqp,current_time+dt/2.,tag) + traction(Xqp,current_time+dt,tag))/6.;
        //
        //  for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
        //  {
        //    for (int c = 0; c < dim; ++c)
        //    {
        //      FUloc(i*dim + c) -= JxW_mid * traction_(c) * phi_f[qp][i] ; // força
        //    }
        //  }
        //}
        //
        if (is_surface)
        {
          surface_energy += JxW * gama(Xqp,current_time,tag);
        }

        if (is_solid)
        {
          // obs gama_s = - gama*cos(theta)
          //surface_energy -= JxW * gama(Xqp,current_time,tag) * cos_theta0(); // comente isso para ser igual do gerbeau
          solid_power +=  JxW * beta_diss() * Uqp.dot(Uqp - solid_veloc(Xqp,current_time,tag));
        }


      } // fim quadratura

    } // end elementos

  }


  // Cl energy (2D only)
  if (dim == 2)
  {
    Vector              line_normal(dim);
    Vector              Xqp(dim);
    Vector              aux(dim);
    Vector              normal(dim);
    Vector              U(dim);
    VectorXi            mapU_f(dim);
    int tag;
    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    for (; point != point_end; ++point)
    {
      tag = point->getTag();
      if (!is_in(tag, triple_tags))
        continue;

      // compute line normal
      {
        Point * sol_point = NULL;
        Point * sol_point_2;
        int iVs[FEPIC_MAX_ICELLS];
        int *iVs_end, *iVs_it;
        Vector aux(dim);

        iVs_end = mesh->connectedVtcs(&*point, iVs);

        // se esse nó está na linha, então existe um vértice vizinho que está no sólido
        for (iVs_it = iVs; iVs_it != iVs_end ; ++iVs_it)
        {
          sol_point_2 = mesh->getNodePtr(*iVs_it);
          if ( is_in(sol_point_2->getTag(), solid_tags) )
            if( !sol_point || (sol_point_2->getTag() > sol_point->getTag()) )
              sol_point = sol_point_2;
        }
        if (!sol_point)
        {
          //FEP_PRAGMA_OMP(critical)
          {
            printf("ERRO: ponto na linha tríplice não tem um vértice vizinho no sólido");
            throw;
          }
        }
        point->getCoord(Xqp.data(),dim);

        normal = solid_normal(Xqp, current_time, tag);

        sol_point->getCoord(aux.data(),dim);


        // choose a line_normal candidate
        line_normal(0) = -normal(1);
        line_normal(1) =  normal(0);

        // check orientation
        if (line_normal.dot(Xqp-aux) < 0)
          line_normal *= -1;
      }

      //
      dof_handler[DH_UNKS].getVariable(VAR_U).getVertexDofs(mapU_f.data(), &*point);

      /*  Pega os valores das variáveis nos graus de liberdade */
      VecGetValues(Vec_up_0, mapU_f.size(), mapU_f.data(), U.data());

      cl_power += gama(Xqp,current_time,tag)*cos_theta0()*line_normal.dot(U);
    }

  }

  ofstream File("ContactHistory", ios::app);

  plc = mesh->getNodePtr(0);
  if (plc!=NULL && mesh->isVertex(plc))
  {
    int dofs[3];
    dof_handler[DH_MESH].getVariable(VAR_M).getVertexDofs(dofs, &*plc);
    VecGetValues(Vec_x_0,  dim, dofs, X.data());
  }
  File.precision(12);
  File << current_time << " "
       << theta_min << " "
       << theta_max << " "
       << sqrt(X(0)*X(0) + X(1)*X(1)) <<" "
       << viscous_power << " "
       << kinetic_energy << " "
       << euler_dissip << " "
       << field_energy << " "
       << surface_energy << " "
       << solid_power << " "
       << cl_power << " "
       << endl;

  File.close();
}

void AppCtx::freePetscObjs()
{
  Destroy(Mat_Jac);
  Destroy(Mat_Jac_m);
  Destroy(Vec_res);
  Destroy(Vec_res_m);
  Destroy(Vec_up_0);
  Destroy(Vec_up_1);
  Destroy(Vec_dup);
  Destroy(Vec_v_mid);
  Destroy(Vec_v_1);
  Destroy(Vec_x_0);
  Destroy(Vec_x_1);
  Destroy(Vec_normal);
  //Destroy(ksp);
  //Destroy(snes);
  SNESDestroy(&snes);
  SNESDestroy(&snes_m);
}

void GetDataVelocity::get_vec(int id, Real * vec_out) const
{
  Point const* point = user.mesh->getNodePtr(id);
  std::vector<int> dofs(user.dim);

  user.getNodeDofs(&*point, DH_UNKS, VAR_U, dofs.data());

  for (int i = 0; i < user.dim; ++i)
    vec_out[i] = q_array[*(dofs.data()+i)];
}

double GetDataPressure::get_data_r(int nodeid) const
{
  Point const*const point = user.mesh->getNodePtr(nodeid);
  if (!user.mesh->isVertex(point))
  {
    int dofs[3];
    // position of the node at edge
    const int m = point->getPosition() - user.mesh->numVerticesPerCell();
    Cell const*const cell = user.mesh->getCellPtr(point->getIncidCell());
    if (user.dim==3)
    {
      const int edge_id = cell->getCornerId(m);
      user.dof_handler[DH_UNKS].getVariable(VAR_P).getCornerDofs(dofs, user.mesh->getCornerPtr(edge_id));
    }
    else
    {
      const int edge_id = cell->getFacetId(m);
      user.dof_handler[DH_UNKS].getVariable(VAR_P).getFacetDofs(dofs, user.mesh->getFacetPtr(edge_id));
    }
    return (q_array[dofs[0]] + q_array[dofs[1]])/2.;
  }
  int dof;
  user.dof_handler[DH_UNKS].getVariable(VAR_P).getVertexDofs(&dof, point);
  return q_array[dof];
}

double GetDataPressCellVersion::get_data_r(int cellid) const
{
  // assume que só há 1 grau de liberdade na célula
  int dof[user.dof_handler[DH_UNKS].getVariable(VAR_P).numDofsPerCell()];
  user.dof_handler[DH_UNKS].getVariable(VAR_P).getCellDofs(dof, user.mesh->getCellPtr(cellid));
  return q_array[dof[0]];
}

void GetDataNormal::get_vec(int id, Real * vec_out) const
{
  Point const* point = user.mesh->getNodePtr(id);
  vector<int> dofs(user.dim);

  user.getNodeDofs(&*point, DH_MESH, VAR_M, dofs.data());

  for (int i = 0; i < user.dim; ++i)
    vec_out[i] = q_array[*(dofs.data()+i)];
}

void GetDataMeshVel::get_vec(int id, Real * vec_out) const
{
  Point const* point = user.mesh->getNodePtr(id);
  vector<int> dofs(user.dim);

  user.getNodeDofs(&*point, DH_MESH, VAR_M, dofs.data());

  for (int i = 0; i < user.dim; ++i)
    vec_out[i] = q_array[*(dofs.data()+i)];
}

int GetDataCellTag::get_data_i(int cellid) const
{
  // assume que só há 1 grau de liberdade na célula
  return user.mesh->getCellPtr(cellid)->getTag();
}

/* petsc functions*/
extern PetscErrorCode FormJacobian(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode Monitor(SNES,PetscInt,PetscReal,void *);


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{

  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

  bool help_return;
  bool erro;
  AppCtx user(argc, argv, help_return, erro);

  if (help_return)
    return 0;

  if (erro)
    return 1;

#ifdef FEP_HAS_OPENMP
  {
    int nthreads = 0;
    #pragma omp parallel
    {
      #pragma omp critical
      nthreads = omp_get_num_threads();
    }
    printf("OpenMP version:\n");
    printf("num threads: %d\n", nthreads);
  }
#else
  {
    printf("Serial version.\n");
  }
#endif


  user.loadMesh();
  user.loadDofs();
  user.evaluateQuadraturePts();

  erro = user.err_checks(); if (erro) return 1;

  // print info
  cout << "mesh: " << user.filename << endl;
  user.mesh->printInfo();
  cout << "\n# velocity unknows: " << user.dof_handler[DH_UNKS].getVariable(VAR_U).numPositiveDofs();
  cout << "\n# preassure unknows: " << user.dof_handler[DH_UNKS].getVariable(VAR_P).numPositiveDofs() << endl;
  user.mesh->printStatistics();
  user.mesh->timer.printTimes();

  cout << "feature_tags =  ";
  for (int i = 0; i < (int)user.feature_tags.size(); ++i)
  {
    cout << user.feature_tags[i] << " ";
  }
  cout << endl;

  printf("on update mesh\n");
  user.onUpdateMesh();
  user.solveTimeProblem();

  cout << "\n";
  user.timer.printTimes();

  user.freePetscObjs();
  PetscFinalize();

  cout << "\a" << endl;
  return 0.;
}


#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
PetscErrorCode FormJacobian(SNES snes,Vec Vec_up_1,Mat *Mat_Jac, Mat *prejac, MatStructure *flag, void *ptr)
{
  AppCtx *user    = static_cast<AppCtx*>(ptr);
  user->formJacobian(snes,Vec_up_1,Mat_Jac,prejac,flag);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
PetscErrorCode FormFunction(SNES snes, Vec Vec_up_1, Vec Vec_fun, void *ptr)
{
  AppCtx *user    = static_cast<AppCtx*>(ptr);
  user->formFunction(snes,Vec_up_1,Vec_fun);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CheckSnesConvergence"
PetscErrorCode CheckSnesConvergence(SNES snes, PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
{
  AppCtx *user    = static_cast<AppCtx*>(ctx);
  user->checkSnesConvergence(snes, it, xnorm, pnorm, fnorm, reason);
  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "FormJacobian_mesh"
PetscErrorCode FormJacobian_mesh(SNES snes,Vec Vec_up_1,Mat *Mat_Jac, Mat *prejac, MatStructure *flag, void *ptr)
{
  AppCtx *user    = static_cast<AppCtx*>(ptr);
  user->formJacobian_mesh(snes,Vec_up_1,Mat_Jac,prejac,flag);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction_mesh"
PetscErrorCode FormFunction_mesh(SNES snes, Vec Vec_up_1, Vec Vec_fun, void *ptr)
{
  AppCtx *user    = static_cast<AppCtx*>(ptr);
  user->formFunction_mesh(snes,Vec_up_1,Vec_fun);
  PetscFunctionReturn(0);
}






