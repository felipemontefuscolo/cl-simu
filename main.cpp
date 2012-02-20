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


PetscErrorCode FormJacobian(SNES snes,Vec Vec_up_1,Mat *Mat_Jac, Mat *prejac, MatStructure *flag, void *ptr);
PetscErrorCode FormFunction(SNES snes, Vec Vec_up_1, Vec Vec_fun, void *ptr);
PetscErrorCode CheckSnesConvergence(SNES snes, PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx);

class AppCtx;
class Statistics;

template<int Coord>
class GetDataVelocity : public DefaultGetDataVtk
{
public:
  GetDataVelocity(double *q_array_, AppCtx const& user_) : user(user_), q_array(q_array_){}
  double get_data_r(int nodeid) const;
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


template<int Coord>
class GetDataNormal : public DefaultGetDataVtk
{
public:
  GetDataNormal(double *q_array_, AppCtx const& user_) : user(user_), q_array(q_array_){}
  double get_data_r(int nodeid) const;
  AppCtx const& user;
  double *q_array;
  virtual ~GetDataNormal() {}
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

  meshAliasesUpdate();

}

void AppCtx::loadDofs()
{
  dofsCreate();
  timer.restart();
  dofsUpdate();
  timer.elapsed("CuthillMcKeeRenumber");
  n_dofs_u_per_cell   = dof_handler_vars.getVariable(0).numDofsPerCell();
  n_dofs_u_per_facet  = dof_handler_vars.getVariable(0).numDofsPerFacet();
  n_dofs_u_per_corner = dof_handler_vars.getVariable(0).numDofsPerCorner();
  n_dofs_p_per_cell   = dof_handler_vars.getVariable(1).numDofsPerCell();
  n_dofs_p_per_facet  = dof_handler_vars.getVariable(1).numDofsPerFacet();
  n_dofs_p_per_corner = dof_handler_vars.getVariable(1).numDofsPerCorner();
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
  maxts                  = 1500;   // max num of time steps
  force_pressure         = false;  // elim null space (auto)
  print_to_matlab        = PETSC_FALSE;  // imprime o jacobiano e a função no formato do matlab
  force_dirichlet        = PETSC_TRUE;   // impõe cc ? (para debug)
  full_diriclet          = PETSC_TRUE;
  force_bound_mesh_vel   = PETSC_FALSE;
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
  PetscOptionsBool("-plot_es", "plot exact solution", "main.cpp", plot_exact_sol, &plot_exact_sol, PETSC_NULL);
  PetscOptionsBool("-family_files", "plot family output", "main.cpp", family_files, &family_files, PETSC_NULL);
  PetscOptionsBool("-has_convec", "convective term", "main.cpp", has_convec, &has_convec, PETSC_NULL);
  PetscOptionsBool("-unsteady", "unsteady problem", "main.cpp", unsteady, &unsteady, PETSC_NULL);
  PetscOptionsBool("-boundary_smoothing", "boundary_smoothing", "main.cpp", boundary_smoothing, &boundary_smoothing, PETSC_NULL);
  PetscOptionsBool("-force_bound_mesh_vel", "force_bound_mesh_vel", "main.cpp", force_bound_mesh_vel, &force_bound_mesh_vel, PETSC_NULL);
  PetscOptionsBool("-renumber_dofs", "renumber dofs", "main.cpp", renumber_dofs, &renumber_dofs, PETSC_NULL);
  PetscOptionsBool("-fprint_ca", "print contact angle", "main.cpp", fprint_ca, &fprint_ca, PETSC_NULL);
  PetscOptionsInt("-quadr_e", "quadrature degree (for calculating the error)", "main.cpp", quadr_degree_err, &quadr_degree_err, PETSC_NULL);
  PetscOptionsInt("-quadr_c", "quadrature degree", "main.cpp", quadr_degree_cell, &quadr_degree_cell, PETSC_NULL);
  PetscOptionsInt("-quadr_f", "quadrature degree (facet)", "main.cpp", quadr_degree_facet, &quadr_degree_facet, PETSC_NULL);
  PetscOptionsInt("-quadr_r", "quadrature degree (corner)", "main.cpp", quadr_degree_corner, &quadr_degree_corner, PETSC_NULL);
  PetscOptionsInt("-print_step", "print_step", "main.cpp", print_step, &print_step, PETSC_NULL);
  PetscOptionsScalar("-beta1", "par vel do fluido", "main.cpp", beta1, &beta1, PETSC_NULL);
  PetscOptionsScalar("-beta2", "par vel elastica", "main.cpp", beta2, &beta2, PETSC_NULL);
  PetscOptionsBool("-ale", "mesh movement", "main.cpp", ale, &ale, PETSC_NULL);
  PetscOptionsGetString(PETSC_NULL,"-fin",finaux,PETSC_MAX_PATH_LEN-1,&flg_fin);
  PetscOptionsGetString(PETSC_NULL,"-fout",foutaux,PETSC_MAX_PATH_LEN-1,&flg_fout);
  PetscOptionsHasName(PETSC_NULL,"-help",&ask_help);

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


  PetscOptionsEnd();

  if(ask_help)
    return true;

  if (flg_fin)
    filename.assign(finaux);

  if (flg_fout)
    filename_out.assign(foutaux);

  if (neumann_tags.size() + interface_tags.size() == 0 )
  {
    full_diriclet = PETSC_TRUE;
    force_pressure = true;
  }
  else
    force_pressure = false;

  if (!unsteady)
  {
    if (ale) {
      printf("ale with steady problem?\n");
      throw;
    }
    //dt = 1.e50;
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
  dof_handler_vars.setMesh(mesh.get());
  dof_handler_vars.addVariable("velo",  shape_phi_c.get(), dim);
  dof_handler_vars.addVariable("pres",  shape_psi_c.get(), 1);
  Matrix<bool, Dynamic, Dynamic> blocks(2,2);
  blocks.setOnes();
  blocks(1,1)=pres_pres_block;
  dof_handler_vars.setVariablesRelationship(blocks.data());

  // mesh velocity
  dof_handler_mesh.setMesh(mesh.get());
  dof_handler_mesh.addVariable("mesh_veloc",  shape_qsi_c.get(), dim);
  dof_handler_mesh.setVariablesRelationship(blocks.data());
}

void AppCtx::dofsUpdate()
{
  dof_handler_vars.SetUp();
  if (renumber_dofs)
    dof_handler_vars.CuthillMcKeeRenumber();
  n_unknowns = dof_handler_vars.numDofs();
  dof_handler_mesh.SetUp();
  n_dofs_u_mesh = dof_handler_mesh.numDofs();
}

PetscErrorCode AppCtx::allocPetscObjs()
{

  PetscErrorCode      ierr;
  ierr = SNESCreate(PETSC_COMM_WORLD, &snes);                   CHKERRQ(ierr);

  //Vec Vec_up_1;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_up_1);                       CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_up_1, PETSC_DECIDE, n_unknowns);              CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_up_1);                                  CHKERRQ(ierr);

  //Vec Vec_up_0;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_up_0);                      CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_up_0, PETSC_DECIDE, n_unknowns);             CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_up_0);                                 CHKERRQ(ierr);

  //Vec Vec_res;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_res);                     CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_res, PETSC_DECIDE, n_unknowns);            CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_res);                                CHKERRQ(ierr);

  //Vec Vec_vmsh_0
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_vmsh_0);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_vmsh_0, PETSC_DECIDE, n_dofs_u_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_vmsh_0);                             CHKERRQ(ierr);

  //Vec Vec_vmsh_med
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_vmsh_med);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_vmsh_med, PETSC_DECIDE, n_dofs_u_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_vmsh_med);                             CHKERRQ(ierr);

  //Vec Vec_xmsh_0;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_xmsh_0);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_xmsh_0, PETSC_DECIDE, n_dofs_u_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_xmsh_0);                             CHKERRQ(ierr);

  //Vec Vec_tmp;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_tmp);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_tmp, PETSC_DECIDE, n_dofs_u_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_tmp);                             CHKERRQ(ierr);

  //Vec Vec_nmsh;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_nmsh);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_nmsh, PETSC_DECIDE, n_dofs_u_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_nmsh);                             CHKERRQ(ierr);

  VectorXi nnz;
  {
    std::vector<std::set<int> > table;
    dof_handler_vars.getSparsityTable(table);

    nnz.resize(n_unknowns);

    //#pragma omp parallel for
    for (int i = 0; i < n_unknowns; ++i)
      nnz[i] = table[i].size();

    // removendo a diagonal nula
    if (!pres_pres_block)
    {
      int const n_p_dofs_total = dof_handler_vars.getVariable(1).totalSize();
      //#pragma omp parallel for
      for (int i = 0; i < n_p_dofs_total; ++i)
      {
        int const dof = dof_handler_vars.getVariable(1).data()[i];
        if (dof >= 0)
          ++nnz[dof];
      }
    }
  }

  //Mat Mat_Jac;
  ierr = MatCreate(PETSC_COMM_WORLD, &Mat_Jac);                                      CHKERRQ(ierr);
  ierr = MatSetSizes(Mat_Jac, PETSC_DECIDE, PETSC_DECIDE, n_unknowns, n_unknowns);   CHKERRQ(ierr);
  ierr = MatSetFromOptions(Mat_Jac);                                                 CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(Mat_Jac, 0, nnz.data());                          CHKERRQ(ierr);
  ierr = MatSetOption(Mat_Jac,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);                  CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD, &snes);                  CHKERRQ(ierr);
  ierr = SNESSetFunction(snes, Vec_res, FormFunction, this);      CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes, Mat_Jac, Mat_Jac, FormJacobian, this); CHKERRQ(ierr);
  //ierr = SNESSetJacobian(snes,Mat_Jac,Mat_Jac,SNESDefaultComputeJacobian,&user);  CHKERRQ(ierr);

  ierr = SNESSetConvergenceTest(snes,CheckSnesConvergence,this,PETSC_NULL); CHKERRQ(ierr);

  ierr = SNESGetKSP(snes,&ksp);                                                  CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);                                                      CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,Mat_Jac,Mat_Jac,SAME_NONZERO_PATTERN);                       CHKERRQ(ierr);
  //ierr = KSPSetType(ksp,KSPPREONLY);                                           CHKERRQ(ierr);
  //ierr = KSPSetType(ksp,KSPGMRES);                                               CHKERRQ(ierr);
  //ierr = PCSetType(pc,PCLU);                                                     CHKERRQ(ierr);
  //ierr = PCFactorSetMatOrderingType(pc, MATORDERINGNATURAL);                         CHKERRQ(ierr);
  //ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);  CHKERRQ(ierr);
  //ierr = SNESSetApplicationContext(snes,this);                 

  //ierr = SNESMonitorSet(snes, SNESMonitorDefault, 0, 0); CHKERRQ(ierr);
  //ierr = SNESMonitorSet(snes,Monitor,0,0);CHKERRQ(ierr);
  //ierr = SNESSetTolerances(snes,0,0,0,13,PETSC_DEFAULT);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
  //ierr = SNESLineSearchSet(snes, SNESLineSearchNo, &user); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

void AppCtx::matrixColoring()
{

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
    dof_handler_vars.getVariable(0).getCellDofs(mapU_c.data(), &*cell);
    dof_handler_vars.getVariable(1).getCellDofs(mapP_c.data(), &*cell);


    MatSetValues(Mat_Jac, mapU_c.size(), mapU_c.data(), mapU_c.size(), mapU_c.data(), Aloc.data(), ADD_VALUES);
    MatSetValues(Mat_Jac, mapU_c.size(), mapU_c.data(), mapP_c.size(), mapP_c.data(), Gloc.data(), ADD_VALUES);
    MatSetValues(Mat_Jac, mapP_c.size(), mapP_c.data(), mapU_c.size(), mapU_c.data(), Dloc.data(), ADD_VALUES);
    if (pres_pres_block)
      MatSetValues(Mat_Jac, mapP_c.size(), mapP_c.data(), mapP_c.size(), mapP_c.data(), Eloc.data(), ADD_VALUES);

  }

  ////test
  //for (int i = 0; i < n_unknowns; ++i)
  //  for (int j = 0; j < n_unknowns; ++j)
  //    MatSetValue(Mat_Jac, i, j, 0.0, ADD_VALUES);

  if (!pres_pres_block)
  {
    int const n_p_dofs_total = dof_handler_vars.getVariable(1).totalSize();
    for (int i = 0; i < n_p_dofs_total; ++i)
    {
      const double zero = 0.0;
      int const dof = dof_handler_vars.getVariable(1).data()[i];
      if (dof>=0)
        MatSetValue(Mat_Jac, dof, dof, zero, ADD_VALUES);
    }
  }



  Assembly(Mat_Jac);
  //MatSetOption(Mat_Jac,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
  //MatSetOption(Mat_Jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);

}

void AppCtx::printMatlabLoader()
{
  FILE *fp = fopen("loadmat.m", "w");
  fprintf(fp, "clear;\n"                   );
  fprintf(fp, "jacob;\n"                   );
  fprintf(fp, "clear zzz;\n"               );
  fprintf(fp, "B=Mat_Jac;\n"                   );
  fprintf(fp, "B(B!=0)=1;\n"               );
  fprintf(fp, "nU = %d;\n",dof_handler_vars.getVariable(0).numDofs() );
  fprintf(fp, "nP = %d;\n",dof_handler_vars.getVariable(1).numDofs() );
  fprintf(fp, "nT = nU + nP;\n"            );
  fprintf(fp, "K=Mat_Jac(1:nU,1:nU);\n"        );
  fprintf(fp, "G=Mat_Jac(1:nU,nU+1:nT);\n"     );
  fprintf(fp, "D=Mat_Jac(nU+1:nT,1:nU);\n"     );
  fprintf(fp, "E=Mat_Jac(nU+1:nT,nU+1:nT);\n"  );
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
    qsi_f[qp].resize(nodes_per_facet);

    dLphi_f[qp].resize(n_dofs_u_per_facet/dim, dim-1);
    dLpsi_f[qp].resize(n_dofs_p_per_facet, dim-1);
    dLqsi_f[qp].resize(nodes_per_facet, dim-1);

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
    parametric_pts = genLineParametricPts(  orderForCtype(ECellType(mesh_cell_type))  );

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
    parametric_pts = genTriParametricPts(  orderForCtype(ECellType(mesh_cell_type))  );

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

void AppCtx::computeDirichletEntries()
{
  // armazenando as entradas da jacobiana que são dirichlet
  int tag;
  int count_dir_entries=0;
  int count_dir_vertices=0;
  int count_dir_corners=0;
  int count_dir_facets=0;
  int count_dir_normal_vertices=0;
  int count_dir_normal_corners=0;
  int count_dir_normal_facets=0;

  bool const has_vertices_dofs = shape_phi_c->numDofsAssociatedToVertice() > 0;
  bool const has_corner_dofs   = shape_phi_c->numDofsAssociatedToCorner()  > 0;
  bool const has_facet_dofs    = shape_phi_c->numDofsAssociatedToFacet()   > 0;

  //#pragma omp parallel private(tag) shared(cout,count_dir_vertices,count_dir_corners,count_dir_facets,count_dir_normal_vertices,count_dir_normal_corners,count_dir_normal_facets) default(none)
  {
    int count_dir_vertices_local=0;
    int count_dir_corners_local=0;
    int count_dir_facets_local=0;
    int count_dir_normal_vertices_local=0;
    int count_dir_normal_corners_local=0;
    int count_dir_normal_facets_local=0;

    Point const* point;
    Facet const* facet;
    Corner const* corner;

    // counting dirichlet vertices
    if (has_vertices_dofs)
    {
      //#pragma omp for nowait
      for (int i = 0; i < n_nodes_total; i++)
      {
        point = mesh->getNode(i);
        if (point->disabled() || (!mesh->isVertex(point)))
          continue;

        tag = point->getTag();

        // dir vertices
        if ( is_in(tag, dirichlet_tags))
          ++count_dir_vertices_local;

        // normals
        if ( is_in(tag, solid_tags) || is_in(tag, triple_tags))
          ++count_dir_normal_vertices_local;

      }
    }

    // counting dirichlet edges (to compute higher order nodes later)
    if (has_corner_dofs)
    {
      //#pragma omp for nowait
      for (int i = 0; i < n_corners_total; i++)
      {
        corner = mesh->getCorner(i);
        if (corner->disabled())
          continue;

        tag = corner->getTag();

        // dir corners
        if ( is_in(tag, dirichlet_tags) )
          ++count_dir_corners_local;

        // normal corners
        if ( is_in(tag, solid_tags) || is_in(tag, triple_tags))
          ++count_dir_normal_corners_local;
      }
    }

    // counting dirichlet faces (to compute higher order nodes later)
    if (has_facet_dofs)
    {
      //#pragma omp for nowait
      ////#pragma omp single
      for (int i = 0; i < n_facets_total; i++)
      {
        facet = mesh->getFacet(i);
        if (facet->disabled())
          continue;

        tag = facet->getTag();

        // dir facets
        if ( is_in(tag, dirichlet_tags) )
          ++count_dir_facets_local;

        // normal corners
        if ( is_in(tag, solid_tags) || is_in(tag, triple_tags))
          ++count_dir_normal_facets_local;
      }
    }

    //#pragma omp critical
    {
      count_dir_vertices        += count_dir_vertices_local;
      count_dir_corners         += count_dir_corners_local;
      count_dir_facets          += count_dir_facets_local;
      count_dir_normal_vertices += count_dir_normal_vertices_local;
      count_dir_normal_corners  += count_dir_normal_corners_local;
      count_dir_normal_facets   += count_dir_normal_facets_local;
    }

  } // end parallel


  count_dir_entries = dim*(count_dir_vertices + count_dir_corners + count_dir_facets) +
                      count_dir_normal_vertices + count_dir_normal_corners + count_dir_normal_facets;
  if (force_pressure)
    ++count_dir_entries;

  // pure dirichlets
  if (count_dir_vertices > (int)dir_vertices.size())
    dir_vertices.reserve(dir_vertices.size()*(1.+grow_factor));
  dir_vertices.resize(count_dir_vertices);

  if (count_dir_corners > (int)dir_corners.size())
    dir_corners.reserve(dir_corners.size()*(1.+grow_factor));
  dir_corners.resize(count_dir_corners);

  if (count_dir_facets > (int)dir_facets.size())
    dir_facets.reserve(dir_facets.size()*(1.+grow_factor));
  dir_facets.resize(count_dir_facets);


  // normal dirichlet
  if (count_dir_normal_vertices > (int)dir_normal_vertices.size())
    dir_normal_vertices.reserve(dir_normal_vertices.size()*(1.+grow_factor));
  dir_normal_vertices.resize(count_dir_normal_vertices);

  if (count_dir_normal_corners > (int)dir_normal_corners.size())
    dir_normal_corners.reserve(dir_normal_corners.size()*(1.+grow_factor));
  dir_normal_corners.resize(count_dir_normal_corners);

  if (count_dir_normal_facets > (int)dir_normal_facets.size())
    dir_normal_facets.reserve(dir_normal_facets.size()*(1.+grow_factor));
  dir_normal_facets.resize(count_dir_normal_facets);


  // all entries
  if (count_dir_entries > (int)dir_entries.size())
    dir_entries.reserve(dir_entries.size()*(1.+grow_factor));
  dir_entries.resize(count_dir_entries);


  cout << "Dirichlet entries " << endl;
  cout << "dir_vertices.size()        : "<< count_dir_vertices        << endl;
  cout << "dir_corners.size()         : "<< count_dir_corners         << endl;
  cout << "dir_facets.size()          : "<< count_dir_facets          << endl;
  cout << "dir_normal_vertices.size() : "<< count_dir_normal_vertices << endl;
  cout << "dir_normal_corners.size()  : "<< count_dir_normal_corners  << endl;
  cout << "dir_normal_facets.size()   : "<< count_dir_normal_facets   << endl;
  cout << "total (dir_entries)        : "<< dir_entries.size() <<endl;


  // é possível paralelizar melhor isso daqui ....
  // store
  ////#pragma omp parallel private(tag) shared(cout,count_dir_vertices,count_dir_corners,count_dir_facets,  count_dir_normal_vertices,count_dir_normal_corners,count_dir_normal_facets) default(none)
  {
    Point const* point;
    Facet const* facet;
    Corner const* corner;
    int kk, mm;

    //#pragma omp sections nowait
    {
      //#pragma omp section
      // counting dirichlet vertices
      if (has_vertices_dofs)
      {
        kk=0;
        mm=0;
        for (int i = 0; i < n_nodes_total; i++)
        {
          point = mesh->getNode(i);
          if (point->disabled() || (!mesh->isVertex(point)))
            continue;
          tag = point->getTag();

          if ( is_in(tag, dirichlet_tags))
            dir_vertices.at(kk++) = i;

          if ( is_in(tag, solid_tags) || is_in(tag, triple_tags))
            dir_normal_vertices.at(mm++) = i;
        }
      }

      //#pragma omp section
      // counting dirichlet edges (to compute higher order nodes later)
      if (has_corner_dofs)
      {
        kk=0;
        mm=0;
        for (int i = 0; i < n_corners_total; i++)
        {
          corner = mesh->getCorner(i);
          if (corner->disabled())
            continue;

          tag = corner->getTag();

          if ( is_in(tag, dirichlet_tags) )
            dir_corners.at(kk++) = i;

          if ( is_in(tag, solid_tags) || is_in(tag, triple_tags))
            dir_normal_corners.at(mm++) = i;
        }
      }

      //#pragma omp section
      // counting dirichlet faces (to compute higher order nodes later)
      if (has_facet_dofs)
      {
        kk=0;
        mm=0;
        for (int i = 0; i < n_facets_total; i++)
        {
          facet = mesh->getFacet(i);
          if (facet->disabled())
            continue;

          tag = facet->getTag();

          if ( is_in(tag, dirichlet_tags) )
            dir_facets.at(kk++) = i;

          if ( is_in(tag, solid_tags) || is_in(tag, triple_tags))
            dir_normal_facets.at(mm++) = i;

        }
      }
    }

  } // end parallel


  int xadj[7];
  int const xadj_size = 7;

  xadj[0] = 0;
  xadj[1] = count_dir_vertices;
  xadj[2] = count_dir_corners         + xadj[1];
  xadj[3] = count_dir_facets          + xadj[2];
  xadj[4] = count_dir_normal_vertices + xadj[3];
  xadj[5] = count_dir_normal_corners  + xadj[4];
  xadj[6] = count_dir_normal_facets   + xadj[5];


  //#pragma omp parallel shared(cout,xadj) default(none)
  {
    Point const* point;
    Facet const* facet;
    Corner const* corner;

    // store dir entries
    std::vector<int> vertice_dofs(dim);
    std::vector<int> corner_dofs(dim);
    std::vector<int> facet_dofs(dim);

    // pure dirichlet
    //#pragma omp for nowait
    for (int i = xadj[0]; i < xadj[1]; i++)
    {
      point = mesh->getNode(dir_vertices[i]);
      dof_handler_vars.getVariable(0).getVertexAssociatedDofs(vertice_dofs.data(), point);
      for (int c = 0; c < dim; ++c)
      {
        dir_entries.at(i*dim + c) = vertice_dofs[c];
      }
    }

    //#pragma omp for nowait
    for (int i = xadj[1]; i < xadj[2]; i++)
    {
      corner = mesh->getCorner(dir_corners[i-xadj[1]]);
      dof_handler_vars.getVariable(0).getCornerAssociatedDofs(corner_dofs.data(), corner);
      for (int c = 0; c < dim; ++c)
        dir_entries.at(i*dim + c) = corner_dofs[c];
    }

    //#pragma omp for nowait
    for (int i = xadj[2]; i < xadj[3]; i++)
    {
      facet = mesh->getFacet(dir_facets[i-xadj[2]]);
      dof_handler_vars.getVariable(0).getFacetAssociatedDofs(facet_dofs.data(), facet);
      for (int c = 0; c < dim; ++c)
        dir_entries.at(i*dim + c) = facet_dofs[c];
    }

    const int shift = xadj[3]*(dim-1);
    // normal dirichlet
    //#pragma omp for nowait
    for (int i = xadj[3]; i < xadj[4]; i++)
    {
      point = mesh->getNode(dir_normal_vertices[i-xadj[3]]);
      dof_handler_vars.getVariable(0).getVertexAssociatedDofs(vertice_dofs.data(), point);
      dir_entries.at(shift + i) = vertice_dofs[0];
    }

    //#pragma omp for nowait
    for (int i = xadj[4]; i < xadj[5]; i++)
    {
      corner = mesh->getCorner(dir_normal_corners[i-xadj[4]]);
      dof_handler_vars.getVariable(0).getCornerAssociatedDofs(corner_dofs.data(), corner);
      dir_entries.at(shift + i) = corner_dofs[0];
    }

    //#pragma omp for nowait
    for (int i = xadj[5]; i < xadj[6]; i++)
    {
      facet = mesh->getFacet(dir_normal_facets[i-xadj[5]]);
      dof_handler_vars.getVariable(0).getFacetAssociatedDofs(facet_dofs.data(), facet);
      dir_entries.at(shift + i) = facet_dofs[0];
    }



  }

  if (force_pressure)
  {
    int p_dof=-1;
    // find first valid dof
    for (int i = 0; p_dof<0 && i<dof_handler_vars.getVariable(1).totalSize(); ++i)
    {
      p_dof = dof_handler_vars.getVariable(1).data()[i];
    }
    dir_entries.at(dim*xadj[xadj_size-1]) = p_dof;
  }



} // end computeDirichletEntries()

void AppCtx::onUpdateMesh()
{
  allocPetscObjs();
  matrixColoring();
  computeDirichletEntries();
}

void AppCtx::setInitialConditions()
{
  Vector    Uf(dim);
  Vector    X(dim);
  Tensor    R(dim,dim);
  VectorXi  dofs(dim);
  int nodeid;
  int tag;
  
  VecSet(Vec_up_0,0.);
  VecSet(Vec_up_1,0.);
  VecSet(Vec_vmsh_0,0.);
  VecSet(Vec_vmsh_med,0.);
  //VecSet(Vec_tmp,0.);
  
  //VecCopy(Vec_up_1,Vec_up_0); // Vec_up_0 <- Vec_up_1
  copyMesh2Vec(Vec_xmsh_0);
  
  updateNormals(&Vec_xmsh_0);
  
  // velocidade inicial
  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    nodeid = mesh->getPointId(&*point);
    tag = point->getTag();

    point->getCoord(X.data());

    Uf = u_initial(X, tag);

    getRotationMatrix(R,&nodeid,1);

    if (mesh->isVertex(&*point))
      dof_handler_vars.getVariable(0).getVertexDofs(dofs.data(), &*point);
    else
    {
      // pega os graus de liberdade
      const int m = point->getPosition() - mesh->numVerticesPerCell();
      Cell const* icell = mesh->getCell(point->getIncidCell());
      if (dim==3)
      {
        Corner *edge = mesh->getCorner(icell->getCornerId(m));
        dof_handler_vars.getVariable(0).getCornerAssociatedDofs(dofs.data(), &*edge);
      }
      else // dim=2
      {
        Facet *edge = mesh->getFacet(icell->getFacetId(m));
        dof_handler_vars.getVariable(0).getFacetAssociatedDofs(dofs.data(), &*edge);
      }

    }
    
    VecSetValues(Vec_up_0, dim, dofs.data(), Uf.data(), INSERT_VALUES);

  } // end point loop 
  Assembly(Vec_up_0);
  VecCopy(Vec_up_0,Vec_up_1);
  
  calcMeshVelocity(Vec_up_0, Vec_xmsh_0, Vec_vmsh_0, 0);
}

PetscErrorCode AppCtx::checkSnesConvergence(SNES snes, PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason)
{
  PetscErrorCode ierr;
  
  ierr = SNESDefaultConverged(snes,it,xnorm,pnorm,fnorm,reason,NULL); CHKERRQ(ierr);
  
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
      //copyMesh2Vec(Vec_tmp);
      //calcMeshVelocity(*Vec_up_k, Vec_tmp, Vec_vmsh_med, current_time+dt);
      //moveMesh(Vec_xmsh_0, Vec_vmsh_0, Vec_vmsh_med, 0.5);
      //
      //// mean mesh velocity
      //VecAXPY(Vec_vmsh_med,1,Vec_vmsh_0);
      //VecScale(Vec_vmsh_med, 0.5);
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

PetscErrorCode AppCtx::solveTimeProblem()
{
  PetscErrorCode      ierr(0);
  //cout << endl;
  //cout << "Convergence Reason: " << reason << endl;
  //cout << "Number of Newton iterations = " << its << endl;
  //cout << "number of linear iterations used by the nonlinear solver " << lits << endl;

  if (print_to_matlab)
    printMatlabLoader();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Solve nonlinear system
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setInitialConditions();
  int its;
    
  printf("initial volume: %.15lf \n", getMeshVolume());
  current_time = 0;
  time_step = 0;
  double Qmax=0;
  double steady_error=1;
  for(;;)
  {
    if ((time_step%print_step) == 0)
    {
      if (family_files)
      {
        double  *q_array;
        double  *nml_array;
        VecGetArray(Vec_up_1, &q_array);
        VecGetArray(Vec_nmsh, &nml_array);
        vtk_printer.writeVtk();
        vtk_printer.addNodeScalarVtk("ux",  GetDataVelocity<0>(q_array, *this));
        vtk_printer.addNodeScalarVtk("uy",  GetDataVelocity<1>(q_array, *this));
        if (dim==3)
          vtk_printer.addNodeScalarVtk("uz",   GetDataVelocity<2>(q_array, *this));
          
        vtk_printer.addNodeScalarVtk("nx",  GetDataNormal<0>(nml_array, *this));
        vtk_printer.addNodeScalarVtk("ny",  GetDataNormal<1>(nml_array, *this));
        if (dim==3)
          vtk_printer.addNodeScalarVtk("nz",   GetDataNormal<2>(nml_array, *this));
        
        if (shape_psi_c->discontinuous())
          vtk_printer.addCellScalarVtk("pressure", GetDataPressCellVersion(q_array, *this));
        else
          vtk_printer.addNodeScalarVtk("pressure", GetDataPressure(q_array, *this));


        //vtk_printer.printPointTagVtk("point_tag");
        //VecRestoreArray(Vec_up_1, &q_array);
        //VecRestoreArray(Vec_nmsh, &nml_array);
        
        ierr = SNESGetIterationNumber(snes,&its);     CHKERRQ(ierr);
        cout << "num snes iterations: " << its << endl;
      }


      cout << endl;
      cout << "current time: " << current_time << endl;
      cout << "time step: "    << time_step  << endl;
      cout << "steady error: " << steady_error << endl;
      
    }

    if (ale)
    {
      //calcMeshVelocity(Vec_up_0, Vec_xmsh_0, Vec_vmsh_0, current_time);
      //VecCopy(Vec_vmsh_0, Vec_vmsh_med);
      //moveMesh(Vec_xmsh_0, Vec_vmsh_0, Vec_vmsh_med, 0.0);
      //////calcMeshVelocity(Vec_up_0, Vec_xmsh_0, Vec_vmsh_med, current_time+dt);
      //calcMeshVelocity(Vec_up_0, Vec_xmsh_0, Vec_vmsh_0, current_time);
      //VecCopy(Vec_vmsh_0, Vec_vmsh_med);
      
      //VecCopy(Vec_vmsh_med, Vec_vmsh_0);
      //VecScale(Vec_vmsh_0, 0.5);
      
      copyMesh2Vec(Vec_xmsh_0);
      calcMeshVelocity(Vec_up_0, Vec_xmsh_0, Vec_vmsh_med, current_time);
      moveMesh(Vec_xmsh_0, Vec_vmsh_0, Vec_vmsh_med, 1);
      //swapMeshWithVec(Vec_xmsh_0);
      //calcMeshVelocity(Vec_up_0, Vec_xmsh_0, Vec_vmsh_med, current_time+dt);
    }
    
    // * SOLVE THE SYSTEM *
    if (solve_the_sys)
      ierr = SNESSolve(snes,PETSC_NULL,Vec_up_1);        CHKERRQ(ierr);
    // * SOLVE THE SYSTEM *
    
    if (plot_exact_sol)
      computeError(Vec_up_1,current_time+dt);
    


    //if (plot_exact_sol)
      //computeError(Vec_up_1,current_time+dt);

    current_time += dt;
    time_step += 1;

    // compute steady error
    Qmax = VecNorm(Vec_up_1, NORM_1);
    VecAXPY(Vec_up_0,-1.0,Vec_up_1);
    steady_error = VecNorm(Vec_up_0, NORM_1)/(Qmax==0.?1.:Qmax);
    
    //VecScale(Vec_up_0,0);
    // update
    VecCopy(Vec_up_1,Vec_up_0); // Vec_up_0 <- Vec_up_1
    
    if (ale)
    {
      copyMesh2Vec(Vec_xmsh_0);
      VecCopy(Vec_vmsh_med,Vec_vmsh_0); // Vec_up_0 <- Vec_up_1
    }
    
    printContactAngle(fprint_ca);
    
    
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
  cout << endl;

  if (family_files)
  {
    double  *q_array;
    double  *nml_array;
    VecGetArray(Vec_up_1, &q_array);
    VecGetArray(Vec_nmsh, &nml_array);
    vtk_printer.writeVtk();
    vtk_printer.addNodeScalarVtk("ux",  GetDataVelocity<0>(q_array, *this));
    vtk_printer.addNodeScalarVtk("uy",  GetDataVelocity<1>(q_array, *this));
    if (dim==3)
      vtk_printer.addNodeScalarVtk("uz",   GetDataVelocity<2>(q_array, *this));
      
    vtk_printer.addNodeScalarVtk("nx",  GetDataNormal<0>(nml_array, *this));
    vtk_printer.addNodeScalarVtk("ny",  GetDataNormal<1>(nml_array, *this));
    if (dim==3)
      vtk_printer.addNodeScalarVtk("nz",   GetDataNormal<2>(nml_array, *this));
    
    if (shape_psi_c->discontinuous())
      vtk_printer.addCellScalarVtk("pressure", GetDataPressCellVersion(q_array, *this));
    else
      vtk_printer.addNodeScalarVtk("pressure", GetDataPressure(q_array, *this));


    vtk_printer.printPointTagVtk("point_tag");
    VecRestoreArray(Vec_up_1, &q_array);
    VecRestoreArray(Vec_nmsh, &nml_array);
  }

  printf("final volume: %.15lf \n", getMeshVolume());


  SNESConvergedReason reason;
  ierr = SNESGetIterationNumber(snes,&its);     CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,&reason);  CHKERRQ(ierr);
  cout << "num snes iterations: " << its << endl;
  cout << "reason: " << reason << endl;
  switch (reason)
  {
  case( SNES_CONVERGED_FNORM_ABS       ): printf("SNES_CONVERGED_FNORM_ABS      /* ||F|| < atol */\n"); break;
  case( SNES_CONVERGED_FNORM_RELATIVE  ): printf("SNES_CONVERGED_FNORM_RELATIVE /* ||F|| < rtol*||F_initial|| */\n"); break;
  case( SNES_CONVERGED_PNORM_RELATIVE  ): printf("SNES_CONVERGED_PNORM_RELATIVE /* Newton computed step size small); break; || delta x || < stol */\n"); break;
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
  //case( SNES_DIVERGED_INNER            ): printf("SNES_DIVERGED_INNER           \n"); break;
  }

  if (solve_the_sys)
    MatrixInfo(Mat_Jac);

  int lits;
  SNESGetLinearSolveIterations(snes,&lits);

  if (plot_exact_sol)
    computeError(Vec_up_1,current_time);
  
  if (unsteady)
  {
    cout << "\nmean errors: \n";
    cout << "# hmean            u_L2_norm         p_L2_norm         grad_u_L2_norm    grad_p_L2_norm" << endl;
    printf( "<hmean>            %.15lf %.15lf %.15lf %.15lf\n", Stats.mean_u_L2_norm(), Stats.mean_p_L2_norm(), Stats.mean_grad_u_L2_norm(), Stats.mean_grad_p_L2_norm());
  }

  PetscFunctionReturn(0);
}


void AppCtx::computeError(Vec &qq, double tt)
{

  MatrixXd            u_coefs_c(n_dofs_u_per_cell/dim, dim);
  MatrixXd            u_coefs_c_trans(dim,n_dofs_u_per_cell/dim);
  VectorXd            p_coefs_c(n_dofs_p_per_cell);
  MatrixXd            x_coefs_c(nodes_per_cell, dim);
  MatrixXd            x_coefs_c_trans(dim, nodes_per_cell);
  Tensor              F_c(dim,dim), invF_c(dim,dim), invFT_c(dim,dim);
  MatrixXd            dxphi_err(n_dofs_u_per_cell/dim, dim);
  MatrixXd            dxpsi_err(n_dofs_p_per_cell, dim);
  MatrixXd            dxqsi_err(nodes_per_cell, dim);
  Tensor              dxU(dim,dim); // grad u
  Vector              dxP(dim);     // grad p
  Vector              Xqp(dim);
  Vector              Uqp(dim);

  double              Pqp;
  VectorXi            cell_nodes(nodes_per_cell);
  double              Jx;
  double              weight;
  int                 tag;

  double              p_L2_norm = 0.;
  double              u_L2_norm = 0.;
  double              grad_u_L2_norm = 0.;
  double              grad_p_L2_norm = 0.;

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
    dof_handler_vars.getVariable(0).getCellDofs(mapU_c.data(), &*cell);
    dof_handler_vars.getVariable(1).getCellDofs(mapP_c.data(), &*cell);

    /*  Pega os valores das variáveis nos graus de liberdade */
    VecGetValues(qq, mapU_c.size(), mapU_c.data(), u_coefs_c.data());
    VecGetValues(qq, mapP_c.size(), mapP_c.data(), p_coefs_c.data());

    u_coefs_c_trans = u_coefs_c.transpose();

    mesh->getCellNodesId(&*cell, cell_nodes.data());
    mesh->getNodesCoords(cell_nodes.begin(), cell_nodes.end(), x_coefs_c.data());
    x_coefs_c_trans = x_coefs_c.transpose();

    for (int qp = 0; qp < n_qpts_err; ++qp)
    {
      F_c    = x_coefs_c_trans * dLqsi_err[qp];
      Jx     = F_c.determinant();
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

      //  note: norm(u, H1)^2 = norm(u, L2)^2 + norm(gradLphi_c, L2)^2
      u_L2_norm        += (u_exact(Xqp, tt, tag) - Uqp).squaredNorm()*weight*Jx;
      p_L2_norm        += sqr(pressure_exact(Xqp, tt, tag) - Pqp)*weight*Jx;
      grad_u_L2_norm   += (grad_u_exact(Xqp, tt, tag) - dxU).squaredNorm()*weight*Jx;
      grad_p_L2_norm   += (grad_p_exact(Xqp, tt, tag) - dxP).squaredNorm()*weight*Jx;

    } // fim quadratura

  } // end elementos

  u_L2_norm      = sqrt(u_L2_norm     );
  p_L2_norm      = sqrt(p_L2_norm     );
  grad_u_L2_norm = sqrt(grad_u_L2_norm);
  grad_p_L2_norm = sqrt(grad_p_L2_norm);

  // compute hmean
  double hmean=0;
  // ASSUME QUE SÓ POSSA TER NO MÁXIMO 1 NÓ POR ARESTA


  VectorXi edge_nodes(3);
  Vector Xa(dim), Xb(dim);
  int n_edges=0;

  if (dim==2)
  //#pragma omp parallel default(none) shared(hmean)
  {
    const int n_edges_total = mesh->numFacetsTotal();
    Facet const* edge(NULL);

    //#pragma omp for nowait
    for (int a = 0; a < n_edges_total; ++a)
    {
      edge = mesh->getFacet(a);
      if (edge->disabled())
        continue;

      mesh->getFacetNodesId(&*edge, edge_nodes.data());

      mesh->getNode(edge_nodes[0])->getCoord(Xa.data());
      mesh->getNode(edge_nodes[1])->getCoord(Xb.data());
      hmean += (Xa-Xb).norm();
      ++n_edges;
    }
  }

  if (dim==3)
  //#pragma omp parallel default(none) shared(cout,hmean)
  {
    const int n_edges_total = mesh->numCornersTotal();
    Corner const* edge(NULL);

    //#pragma omp for nowait
    for (int a = 0; a < n_edges_total; ++a)
    {
      edge = mesh->getCorner(a);
      if (edge->disabled())
        continue;

      mesh->getCornerNodesId(&*edge, edge_nodes.data());

      mesh->getNode(edge_nodes[0])->getCoord(Xa.data());
      mesh->getNode(edge_nodes[1])->getCoord(Xb.data());
      hmean += (Xa-Xb).norm();
      ++n_edges;
    }
  }
  hmean /= n_edges;

  if (time_step==maxts)
  {
    cout << endl;
    cout << "errors computed at last time step: "<< endl;
    cout << "# hmean            u_L2_norm         p_L2_norm         grad_u_L2_norm    grad_p_L2_norm" << endl;
    printf("%.15lf %.15lf %.15lf %.15lf %.15lf\n",hmean, u_L2_norm, p_L2_norm, grad_u_L2_norm, grad_p_L2_norm);
  }
  
  Stats.add_u_L2(u_L2_norm/maxts);
  Stats.add_p_L2(p_L2_norm/maxts);
  Stats.add_grad_u_L2(grad_u_L2_norm/maxts);
  Stats.add_grad_p_L2(grad_p_L2_norm/maxts);
  
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
    dof_handler_vars.getVariable(0).getVertexDofs(vtx_dofs_fluid.data(), &*point);
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
  
  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    tag = point->getTag();
    
    if (!is_in(tag,triple_tags))
      continue;
    
    point->getCoord(X.data());
    normal_solid = solid_normal(X, current_time, tag);
    dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_mesh.data(), &*point);
    VecGetValues(Vec_nmsh, dim, vtx_dofs_mesh.data(), normal_surf.data());   
    
    tet = asin(sqrt(1. - pow(normal_solid.dot(normal_surf),2)  ));
    
    if (tet < theta_min)
      theta_min = tet;
    if (tet > theta_max)
      theta_max = tet;
  }
  
  theta_min = theta_min*180./pi;
  theta_max = theta_max*180./pi;
  
  cout << "theta min: " << theta_min << "\ntheta max: " << theta_max << endl;
  
  ofstream File("ContactHistory", ios::app);
  
  File.precision(12);
  File << current_time << " " << theta_min << " " << theta_max << endl;
  
  File.close();
}

void AppCtx::freePetscObjs()
{
  Destroy(Mat_Jac);
  Destroy(Vec_up_1);
  Destroy(Vec_up_0);
  Destroy(Vec_res);
  Destroy(Vec_vmsh_0);
  Destroy(Vec_vmsh_med);
  Destroy(Vec_xmsh_0);
  Destroy(Vec_tmp);
  Destroy(Vec_nmsh);
  //Destroy(ksp);
  //Destroy(snes);
  SNESDestroy(&snes);
}


template<int Coord>
double GetDataVelocity<Coord>::get_data_r(int nodeid) const
{
  Point const* point = user.mesh->getNode(nodeid);
  Tensor R(user.dim,user.dim);
  Vector Ur(user.dim);
  Vector tmp(user.dim);
  user.getRotationMatrix(R,&nodeid,1);
  int dofs[user.dim];
  if (!user.mesh->isVertex(point))
  {
    const int m = point->getPosition() - user.mesh->numVerticesPerCell();
    Cell const* cell = user.mesh->getCell(point->getIncidCell());
    if (user.dim==3)
    {
      const int edge_id = cell->getCornerId(m);
      user.dof_handler_vars.getVariable(0).getCornerAssociatedDofs(dofs, user.mesh->getCorner(edge_id));
    }
    else
    {
      const int edge_id = cell->getFacetId(m);
      user.dof_handler_vars.getVariable(0).getFacetAssociatedDofs(dofs, user.mesh->getFacet(edge_id));
    }

    for (int i = 0; i < user.dim; ++i)
      Ur(i) = q_array[dofs[i]];
    //return q_array[dofs[2*user.dim + Coord]];

  }
  else
  {
    //VectorXi dofs(user.dim);
    user.dof_handler_vars.getVariable(0).getVertexDofs(dofs, point);
    for (int i = 0; i < user.dim; ++i)
      Ur(i) = q_array[dofs[i]];
  }

  tmp = R.transpose()*Ur;
  return tmp(Coord);
}

double GetDataPressure::get_data_r(int nodeid) const
{
  Point const*const point = user.mesh->getNode(nodeid);
  if (!user.mesh->isVertex(point))
  {
    int dofs[3];
    // position of the node at edge
    const int m = point->getPosition() - user.mesh->numVerticesPerCell();
    Cell const*const cell = user.mesh->getCell(point->getIncidCell());
    if (user.dim==3)
    {
      const int edge_id = cell->getCornerId(m);
      user.dof_handler_vars.getVariable(1).getCornerDofs(dofs, user.mesh->getCorner(edge_id));
    }
    else
    {
      const int edge_id = cell->getFacetId(m);
      user.dof_handler_vars.getVariable(1).getFacetDofs(dofs, user.mesh->getFacet(edge_id));
    }
    return (q_array[dofs[0]] + q_array[dofs[1]])/2.;
  }
  int dof;
  user.dof_handler_vars.getVariable(1).getVertexDofs(&dof, point);
  return q_array[dof];
  return 0;
}

double GetDataPressCellVersion::get_data_r(int cellid) const
{
  // assume que só há 1 grau de liberdade na célula
  int dof[user.dof_handler_vars.getVariable(1).numDofsPerCell()];
  user.dof_handler_vars.getVariable(1).getCellDofs(dof, user.mesh->getCell(cellid));
  return q_array[dof[0]];
}

template<int Coord>
double GetDataNormal<Coord>::get_data_r(int nodeid) const
{
  Point const* point = user.mesh->getNode(nodeid);
  Vector Normal(user.dim);
  int dofs[user.dim];
  if (!user.mesh->isVertex(point))
  {
    const int m = point->getPosition() - user.mesh->numVerticesPerCell();
    Cell const* cell = user.mesh->getCell(point->getIncidCell());
    if (user.dim==3)
    {
      const int edge_id = cell->getCornerId(m);
      user.dof_handler_vars.getVariable(0).getCornerAssociatedDofs(dofs, user.mesh->getCorner(edge_id));
    }
    else
    {
      const int edge_id = cell->getFacetId(m);
      user.dof_handler_vars.getVariable(0).getFacetAssociatedDofs(dofs, user.mesh->getFacet(edge_id));
    }

    for (int i = 0; i < user.dim; ++i)
      Normal(i) = q_array[dofs[i]];
    //return q_array[dofs[2*user.dim + Coord]];

  }
  else
  {
    //VectorXi dofs(user.dim);
    user.dof_handler_mesh.getVariable(0).getVertexDofs(dofs, point);
    for (int i = 0; i < user.dim; ++i)
      Normal(i) = q_array[dofs[i]];
  }

  return Normal(Coord);
}


/* petsc functions*/
extern PetscErrorCode FormJacobian(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode Monitor(SNES,PetscInt,PetscReal,void *);


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{

  PetscInitialize(&argc,&argv);

  bool help_return;
  bool erro;
  AppCtx user(argc, argv, help_return, erro);

  if (help_return)
    return 0;

  if (erro)
    return 1;


  user.loadMesh();
  user.loadDofs();
  user.evaluateQuadraturePts();

  erro = user.err_checks(); if (erro) return 1;

  // print info
  cout << "mesh: " << user.filename << endl;
  user.mesh->printInfo();
  cout << "\n# velocity unknows: " << user.dof_handler_vars.getVariable(0).numDofs();
  cout << "\n# preassure unknows: " << user.dof_handler_vars.getVariable(1).numDofs() << endl;
  user.mesh->printStatistics();
  user.mesh->timer.printTimes();

  user.onUpdateMesh();
  user.solveTimeProblem();

  cout << "\n";
  user.timer.printTimes();

  user.freePetscObjs();
  PetscFinalize();

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











/*
  // ***********
  // form the residue of the contact line (2D)
  void formCornerFunction(Point * point,
                            Vector &Ur,
                            Vector &FUloc)
  {
    bool is_triple;
    int tag;
    //Point * point;
    Point * sol_point; // solid point
    int iVs[FEPIC_MAX_ICELLS];
    int *iVs_end;
    int *iVs_it;
    Vector line_normal(dim);
    Vector aux(dim);
    Vector normal(dim); // solid normal
    Vector X(dim);
    VectorXi vtx_dofs_fluid(dim);
    int nodeid;
    //Vector FUloc(dim);
    Vector U(dim);
    //Vector Ur(dim); // rotated
    Tensor R(dim,dim);
    Vector tmp(dim);

    tag = point->getTag();

    is_triple = is_in(tag, triple_tags);

    if (!is_triple)
      return;


    iVs_end = mesh->connectedVtcs(point, iVs);

    // se esse nó está na linha, então existe um vértice vizinho que está no sólido
    for (iVs_it = iVs; iVs_it != iVs_end ; ++iVs_it)
    {
      sol_point = mesh->getNode(*iVs_it);
      if ( is_in(sol_point->getTag(), solid_tags) )
        break;
    }
    if (iVs_it == iVs_end)
    {
      //#pragma omp critical
      {
        printf("ERRO: ponto na linha tríplice não tem um vértice vizinho no sólido");
        throw;
      }
    }

    normal = solid_normal(X, current_time, tag);

    point->getCoord(X.data());

    sol_point->getCoord(aux.data());


    // choose a line_normal candidate
    line_normal(0) = -normal(1);
    line_normal(1) =  normal(0);

    // check orientation
    if (line_normal.dot(X-aux) < 0)
      line_normal *= -1;

    dof_handler_vars.getVariable(0).getVertexDofs(vtx_dofs_fluid.data(), &*point);

    nodeid = mesh->getPointId(point);
    getRotationMatrix(R, &nodeid, 1);

    U = R.transpose()*Ur;

    FUloc(0) = (-gama(X, current_time, tag)*cos_theta0() + zeta(0,0)*line_normal.dot(U))*line_normal(0);
    FUloc(1) = (-gama(X, current_time, tag)*cos_theta0() + zeta(0,0)*line_normal.dot(U))*line_normal(1);

    nodeid = mesh->getPointId(&*point);

    rotate_RA(R,FUloc,tmp);
  }
*/

