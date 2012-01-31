//static char help[] = "Navier-Stokes.\n\n";

/*
 * ALE
 *
 *  pho*( (Ut + (U-Umsh) · nabla)U ) + grad p = niu* div grad U + force
 *
 *  LIMITAÇÕES:
 *  * triangulo e tetrahedro apenas.
 *  * no maximo 1 grau de liberdade associado a uma aresta
 *  * condição de dirichlet na normal só pode valer 0 .. procure por UNORMAL_AQUI
 *  * elementos isoparamétricos e subparamétricos
 */

#include "common.hpp"


PetscErrorCode FormJacobian(SNES snes,Vec x,Mat *Jac, Mat *prejac, MatStructure *flag, void *ptr);
PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ptr);

class AppCtx;

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
  //Re                     = 0.0;
  dt                     = 0.1;
  unsteady               = PETSC_TRUE;
  boundary_smoothing     = PETSC_TRUE;
  steady_tol             = 1.e-6;
  theta                  = 1;      // time step, theta method
  maxts                  = 1500;   // max num of time steps
  force_pressure         = false;  // elim null space (auto)
  print_to_matlab        = PETSC_FALSE;  // imprime o jacobiano e a função no formato do matlab
  force_dirichlet        = PETSC_TRUE;   // impõe cc ? (para debug)
  full_diriclet          = PETSC_TRUE;
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
  renumber_dofs          = PETSC_TRUE;

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
  PetscOptionsScalar("-theta", "theta value", "main.cpp", theta, &theta, PETSC_NULL);
  PetscOptionsScalar("-sst", "steady state tolerance", "main.cpp", steady_tol, &steady_tol, PETSC_NULL);
  PetscOptionsBool("-print_to_matlab", "print jacobian to matlab", "main.cpp", print_to_matlab, &print_to_matlab, PETSC_NULL);
  PetscOptionsBool("-force_dirichlet", "force dirichlet bound cond", "main.cpp", force_dirichlet, &force_dirichlet, PETSC_NULL);
  PetscOptionsBool("-plot_es", "plot exact solution", "main.cpp", plot_exact_sol, &plot_exact_sol, PETSC_NULL);
  PetscOptionsBool("-family_files", "plot family output", "main.cpp", family_files, &family_files, PETSC_NULL);
  PetscOptionsBool("-has_convec", "convective term", "main.cpp", has_convec, &has_convec, PETSC_NULL);
  PetscOptionsBool("-unsteady", "unsteady problem", "main.cpp", unsteady, &unsteady, PETSC_NULL);
  PetscOptionsBool("-boundary_smoothing", "boundary_smoothing", "main.cpp", boundary_smoothing, &boundary_smoothing, PETSC_NULL);
  PetscOptionsBool("-renumber_dofs", "renumber dofs", "main.cpp", renumber_dofs, &renumber_dofs, PETSC_NULL);
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
  dof_handler_mesh.addVariable("mesh_veloc",  shape_phi_c.get(), dim);
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

  //Vec q;
  ierr = VecCreate(PETSC_COMM_WORLD, &q);                       CHKERRQ(ierr);
  ierr = VecSetSizes(q, PETSC_DECIDE, n_unknowns);              CHKERRQ(ierr);
  ierr = VecSetFromOptions(q);                                  CHKERRQ(ierr);

  //Vec q0;
  ierr = VecCreate(PETSC_COMM_WORLD, &q0);                      CHKERRQ(ierr);
  ierr = VecSetSizes(q0, PETSC_DECIDE, n_unknowns);             CHKERRQ(ierr);
  ierr = VecSetFromOptions(q0);                                 CHKERRQ(ierr);

  //Vec res;
  ierr = VecCreate(PETSC_COMM_WORLD, &res);                     CHKERRQ(ierr);
  ierr = VecSetSizes(res, PETSC_DECIDE, n_unknowns);            CHKERRQ(ierr);
  ierr = VecSetFromOptions(res);                                CHKERRQ(ierr);

  //Vec u_mesh;
  ierr = VecCreate(PETSC_COMM_WORLD, &u_mesh);                  CHKERRQ(ierr);
  ierr = VecSetSizes(u_mesh, PETSC_DECIDE, n_dofs_u_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(u_mesh);                             CHKERRQ(ierr);

  //Vec x_mesh;
  ierr = VecCreate(PETSC_COMM_WORLD, &x_mesh);                  CHKERRQ(ierr);
  ierr = VecSetSizes(x_mesh, PETSC_DECIDE, n_dofs_u_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(x_mesh);                             CHKERRQ(ierr);

  //Vec nml_mesh;
  ierr = VecCreate(PETSC_COMM_WORLD, &nml_mesh);                  CHKERRQ(ierr);
  ierr = VecSetSizes(nml_mesh, PETSC_DECIDE, n_dofs_u_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(nml_mesh);                             CHKERRQ(ierr);

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

  //Mat Jac;
  ierr = MatCreate(PETSC_COMM_WORLD, &Jac);                                      CHKERRQ(ierr);
  ierr = MatSetSizes(Jac, PETSC_DECIDE, PETSC_DECIDE, n_unknowns, n_unknowns);   CHKERRQ(ierr);
  ierr = MatSetFromOptions(Jac);                                                 CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(Jac, 0, nnz.data());                          CHKERRQ(ierr);
  ierr = MatSetOption(Jac,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);                  CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD, &snes);                  CHKERRQ(ierr);
  ierr = SNESSetFunction(snes, res, FormFunction, this);      CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes, Jac, Jac, FormJacobian, this); CHKERRQ(ierr);
  //ierr = SNESSetJacobian(snes,Jac,Jac,SNESDefaultComputeJacobian,&user);  CHKERRQ(ierr);

  ierr = SNESGetKSP(snes,&ksp);                                                  CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);                                                      CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,Jac,Jac,SAME_NONZERO_PATTERN);                       CHKERRQ(ierr);
  //ierr = KSPSetType(ksp,KSPPREONLY);                                           CHKERRQ(ierr);
  //ierr = KSPSetType(ksp,KSPGMRES);                                               CHKERRQ(ierr);
  //ierr = PCSetType(pc,PCLU);                                                     CHKERRQ(ierr);
  //ierr = PCFactorSetMatOrderingType(pc, MATORDERINGNATURAL);                         CHKERRQ(ierr);
  //ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);  CHKERRQ(ierr);


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


    MatSetValues(Jac, mapU_c.size(), mapU_c.data(), mapU_c.size(), mapU_c.data(), Aloc.data(), ADD_VALUES);
    MatSetValues(Jac, mapU_c.size(), mapU_c.data(), mapP_c.size(), mapP_c.data(), Gloc.data(), ADD_VALUES);
    MatSetValues(Jac, mapP_c.size(), mapP_c.data(), mapU_c.size(), mapU_c.data(), Dloc.data(), ADD_VALUES);
    if (pres_pres_block)
      MatSetValues(Jac, mapP_c.size(), mapP_c.data(), mapP_c.size(), mapP_c.data(), Eloc.data(), ADD_VALUES);

  }

  ////test
  //for (int i = 0; i < n_unknowns; ++i)
  //  for (int j = 0; j < n_unknowns; ++j)
  //    MatSetValue(Jac, i, j, 0.0, ADD_VALUES);

  if (!pres_pres_block)
  {
    int const n_p_dofs_total = dof_handler_vars.getVariable(1).totalSize();
    for (int i = 0; i < n_p_dofs_total; ++i)
    {
      const double zero = 0.0;
      int const dof = dof_handler_vars.getVariable(1).data()[i];
      if (dof>=0)
        MatSetValue(Jac, dof, dof, zero, ADD_VALUES);
    }
  }



  Assembly(Jac);
  //MatSetOption(Jac,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
  //MatSetOption(Jac,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);

}

void AppCtx::printMatlabLoader()
{
  FILE *fp = fopen("loadmat.m", "w");
  fprintf(fp, "clear;\n"                   );
  fprintf(fp, "jacob;\n"                   );
  fprintf(fp, "clear zzz;\n"               );
  fprintf(fp, "B=Jac;\n"                   );
  fprintf(fp, "B(B!=0)=1;\n"               );
  fprintf(fp, "nU = %d;\n",dof_handler_vars.getVariable(0).numDofs() );
  fprintf(fp, "nP = %d;\n",dof_handler_vars.getVariable(1).numDofs() );
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

    dLphi_r[qp].resize(n_dofs_u_per_corner/dim, dim);
    dLpsi_r[qp].resize(n_dofs_p_per_corner, dim);
    dLqsi_r[qp].resize(nodes_per_corner, dim);

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
    std::vector<Eigen::Vector2d> tri_parametric_pts;
    tri_parametric_pts = genTriParametricPts(  orderForCtype(ECellType(mesh_cell_type))  );

    dLphi_nf.resize(tri_parametric_pts.size());

    for (int k = 0; k < (int)tri_parametric_pts.size(); ++k)
    {
      dLphi_nf[k].resize(n_dofs_u_per_facet/dim, dim-1);
      for (int n = 0; n < n_dofs_u_per_facet/dim; ++n)
        for (int d = 0; d < dim-1; ++d)
          dLphi_nf[k](n, d) = shape_phi_f->gradL(tri_parametric_pts[k].data(), n, d);
    }
  }
  else
  if(dim==3)
  {
    std::vector<Eigen::Vector3d> tet_parametric_pts;
    tet_parametric_pts = genTetParametricPts(  orderForCtype(ECellType(mesh_cell_type))  );

    dLphi_nf.resize(tet_parametric_pts.size());

    for (int k = 0; k < (int)tet_parametric_pts.size(); ++k)
    {
      dLphi_nf[k].resize(n_dofs_u_per_facet/dim, dim-1);
      for (int n = 0; n < n_dofs_u_per_facet/dim; ++n)
        for (int d = 0; d < dim-1; ++d)
          dLphi_nf[k](n, d) = shape_phi_f->gradL(tet_parametric_pts[k].data(), n, d);
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
  VecSet(q0,0.);
  VecSet(q,0.);
  VecSet(u_mesh,0.);
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
  current_time = 0;
  time_step = 0;
  double Umax=0;
  double steady_error=1;
  do
  {
    VecCopy(q,q0); // q0 <- q
    if (ale) {
      calcMeshVelocity(q);
    }


    if ((time_step%print_step) == 0)
    {
      if (family_files)
      {
        double  *q_array;
        double  *nml_array;
        VecGetArray(q, &q_array);
        VecGetArray(nml_mesh, &nml_array);
        vtk_printer.writeVtk();
        vtk_printer.addNodeScalarVtk("u_x",  GetDataVelocity<0>(q_array, *this));
        vtk_printer.addNodeScalarVtk("u_y",  GetDataVelocity<1>(q_array, *this));
        if (dim==3)
          vtk_printer.addNodeScalarVtk("u_z",   GetDataVelocity<2>(q_array, *this));
        if (shape_psi_c->discontinuous())
          vtk_printer.addCellScalarVtk("pressure", GetDataPressCellVersion(q_array, *this));
        else
          vtk_printer.addNodeScalarVtk("pressure", GetDataPressure(q_array, *this));

        vtk_printer.addNodeScalarVtk("n_x",  GetDataNormal<0>(nml_array, *this));
        vtk_printer.addNodeScalarVtk("n_y",  GetDataNormal<1>(nml_array, *this));
        if (dim==3)
          vtk_printer.addNodeScalarVtk("n_z",   GetDataNormal<2>(nml_array, *this));

        //vtk_printer.printPointTagVtk("point_tag");
        VecRestoreArray(q, &q_array);
        VecRestoreArray(nml_mesh, &nml_array);
      }


      cout << endl;
      cout << "current time: " << current_time << endl;
      cout << "time step: "    << time_step  << endl;
      cout << "steady error: " << steady_error << endl << endl;
    }

    if (solve_the_sys)
      ierr = SNESSolve(snes,PETSC_NULL,q);        CHKERRQ(ierr);

    if (ale)
      moveMesh();

    //formFunction(snes, q, res);


    current_time += dt;
    time_step += 1;

    Umax = VecNorm(q, NORM_1);
    VecAXPY(q0,-1.0,q);
    steady_error = VecNorm(q0, NORM_1)/(Umax==0.?1.:Umax);
  }
  while (time_step < maxts && steady_error > steady_tol);


  SNESConvergedReason reason;
  int its;
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
    MatrixInfo(Jac);

  int lits;
  SNESGetLinearSolveIterations(snes,&lits);

  if (plot_exact_sol)
    computeError(q,-1);

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
  MatrixXd            dxphi_err;
  MatrixXd            dxpsi_err;
  MatrixXd            dxqsi_err;
  Tensor              dxU(dim,dim); // grad u
  Vector              dxP(dim);     // grad p
  Vector              Xqp(dim);
  Vector              Uqp(dim);

  double              Pqp;
  VectorXi          cell_nodes(nodes_per_cell);
  double              Jx;
  double              weight;
  int                 tag;

  double              p_L2_norm = 0.;
  double              u_L2_norm = 0.;
  double              grad_u_L2_norm = 0.;
  double              grad_p_L2_norm = 0.;

  VectorXi               mapU_c(n_dofs_u_per_cell);
  VectorXi               mapU_r(n_dofs_u_per_corner);
  VectorXi               mapP_c(n_dofs_p_per_cell);
  VectorXi               mapP_r(n_dofs_p_per_corner);
  VectorXi               mapM_c(dim*nodes_per_cell);
  VectorXi               mapM_r(dim*nodes_per_corner);


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
      dxU  = u_coefs_c_trans * dxphi_err;       // n+theta

      Xqp  = x_coefs_c_trans * qsi_err[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
      Uqp  = u_coefs_c_trans * phi_err[qp];
      Pqp  = p_coefs_c.dot(psi_err[qp]);

      weight = quadr_err->weight(qp);

      //  note: norm(u, H1)^2 = norm(u, L2)^2 + norm(gradLphi_c, L2)^2
      u_L2_norm        += (u_boundary(Xqp, tt, tag) - Uqp).squaredNorm()*weight*Jx;
      p_L2_norm        += sqr(pressure_exact(Xqp, tt, tag) - Pqp)*weight*Jx;
      grad_u_L2_norm   += (grad_u_exact(Xqp, tt, tag) - dxU).squaredNorm()*weight*Jx;
      grad_p_L2_norm   += (grad_p_exact(Xqp, tt, tag) - dxP).squaredNorm()*weight*Jx;

    } // fim quadratura

  } // end elementos

  u_L2_norm      = sqrt(u_L2_norm     );
  p_L2_norm      = sqrt(p_L2_norm     );
  grad_u_L2_norm = sqrt(grad_u_L2_norm);
  grad_p_L2_norm = sqrt(grad_p_L2_norm);

  // compute hmin
  double hmin=99999999;
  const int high_order = mesh->numNodesPerCell() > mesh->numVerticesPerCell();
  // ASSUME QUE SÓ POSSA TER NO MÁXIMO 1 NÓ POR ARESTA

  if (dim==2)
  //#pragma omp parallel default(none) shared(hmin)
  {
    const int n_facets_total = mesh->numFacetsTotal();
    Facet const* facet(NULL);
    double hmin_local=99999999, dist(0);
    VectorXi facet_nodes(nodes_per_facet);
    Vector Xa(dim), Xb(dim);

    //#pragma omp for nowait
    for (int a = 0; a < n_facets_total; ++a)
    {
      facet = mesh->getFacet(a);
      if (facet->disabled())
        continue;

      mesh->getFacetNodesId(&*facet, facet_nodes.data());

      mesh->getNode(facet_nodes[0])->getCoord(Xa.data());

      for (int i = 1; i < nodes_per_facet-high_order; ++i)
      {
        mesh->getNode(facet_nodes[i])->getCoord(Xb.data());
        dist = (Xa-Xb).norm();
        if (hmin_local >  dist)
          hmin_local = dist;
      }
    }
    //#pragma omp critical
    {
      if (hmin_local < hmin)
        hmin = hmin_local;
    }
  }

  if (dim==3)
  //#pragma omp parallel default(none) shared(cout,hmin)
  {
    const int n_corners_total = mesh->numCornersTotal();
    Corner const* corner(NULL);
    double hmin_local=99999999, dist(0);
    VectorXi corner_nodes(nodes_per_corner);
    Vector Xa(dim), Xb(dim);

    //#pragma omp for nowait
    for (int a = 0; a < n_corners_total; ++a)
    {
      corner = mesh->getCorner(a);
      if (corner->disabled())
        continue;

      mesh->getCornerNodesId(&*corner, corner_nodes.data());

      mesh->getNode(corner_nodes[0])->getCoord(Xa.data());


      for (int i = 1; i < nodes_per_corner-high_order; ++i)
      {
        mesh->getNode(corner_nodes[i])->getCoord(Xb.data());
        dist = (Xa-Xb).norm();
        if (hmin_local >  dist)
          hmin_local = dist;

      }
    }
    //#pragma omp critical
    {
      if (hmin_local < hmin)
        hmin = hmin_local;
    }
  }

  cout << endl;
  cout << "# hmin            u_L2_norm         p_L2_norm         grad_u_L2_norm    grad_p_L2_norm" << endl;
  printf("%.15lf %.15lf %.15lf %.15lf %.15lf\n",hmin, u_L2_norm, p_L2_norm, grad_u_L2_norm, grad_p_L2_norm);

}

void AppCtx::getNormalsFromMesh(Vec *x_mesh)
{

  VectorXi          facet_nodes(nodes_per_facet);
  MatrixXd             x_coefs(nodes_per_facet, dim);                // coordenadas nodais da célula
  MatrixXd            x_coefs_trans(dim, nodes_per_facet);
  Vector              X;
  Vector              normal(dim);
  Tensor              F(dim,dim-1);
  VectorXi               map(n_dofs_u_per_facet);
  //bool                is_surface, is_solid;
  int                 tag;
  bool                virtual_mesh;

  if (x_mesh==NULL)
    virtual_mesh = false;
  else
    virtual_mesh = true;


  // LOOP NAS FACES DO CONTORNO
  facet_iterator facet = mesh->facetBegin();
  facet_iterator facet_end = mesh->facetEnd();
  for (; facet != facet_end; ++facet)
  {
    tag = facet->getTag();


    if (!mesh->inBoundary(&*facet))
      continue;

    //is_surface = is_in(tag, interface_tags);
    //is_solid   = is_in(tag, solid_tags);
    //
    //if ( !is_surface && !is_solid)
    //  continue;

    dof_handler_mesh.getVariable(0).getFacetDofs(map.data(), &*facet);

    mesh->getFacetNodesId(&*facet, facet_nodes.data());
    if (virtual_mesh)
      VecGetValues(*x_mesh, map.size(), map.data(), x_coefs.data());
    else
      mesh->getNodesCoords(facet_nodes.begin(), facet_nodes.end(), x_coefs.data());
    x_coefs_trans = x_coefs.transpose();

    // find the normal
    for (int k = 0; k < nodes_per_facet; ++k)
    {
      tag = mesh->getNode(facet_nodes(k))->getTag();

      if (is_in(tag, solid_tags))
        continue;

      F   = x_coefs_trans * dLphi_nf[k];

      if (dim==2)
      {
        normal(0) = +F(1,0);
        normal(1) = -F(0,0);
      }
      else
      {
        normal = F.col(0);
        X = F.col(1);
        // a = a x b
        cross(normal, X);
      }

      VecSetValues(nml_mesh, dim, map.data()+k*dim, normal.data(), ADD_VALUES);


    } // nodes

  }
  Assembly(nml_mesh);

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    tag = point->getTag();

    if (!mesh->inBoundary(&*point))
      continue;

    if (!mesh->isVertex(&*point))
    {
      const int m = point->getPosition() - mesh->numVerticesPerCell();
      Cell const* cell = mesh->getCell(point->getIncidCell());
      if (dim==3)
      {
        const int edge_id = cell->getCornerId(m);
        Corner const* corner = mesh->getCorner(edge_id);
        dof_handler_mesh.getVariable(0).getCornerAssociatedDofs(map.data(), corner);
      }
      else
      {
        const int edge_id = cell->getFacetId(m);
        Facet const* facet = mesh->getFacet(edge_id);
        dof_handler_mesh.getVariable(0).getFacetAssociatedDofs(map.data(), facet);
      }
    }
    else
      dof_handler_mesh.getVariable(0).getVertexDofs(map.data(), &*point);

    if (!is_in(tag, solid_tags) && !is_in(tag, triple_tags))
    {
      VecGetValues(nml_mesh, dim, map.data(),normal.data());
      normal.normalize();
      VecSetValues(nml_mesh, dim, map.data(),normal.data(), INSERT_VALUES);
      Assembly(nml_mesh);
    }
    else
    {
      if (virtual_mesh)
        VecGetValues(*x_mesh, dim, map.data(), X.data());
      else
        point->getCoord(X.data());
      normal = -solid_normal(X,current_time,tag);

      VecSetValues(nml_mesh, dim, map.data(), normal.data(), INSERT_VALUES);

    }

  }



}

// @param[in] q unknows vector with fluid velocity
// @param[out] u_mesh
PetscErrorCode AppCtx::calcMeshVelocity(Vec const& q)
{
  int        nodeid;
  //double     *x_mesh_array;
  bool const u_has_edge_assoc_dof = shape_phi_c->numDofsAssociatedToFacet()+shape_phi_c->numDofsAssociatedToCorner() > 0;

  Vector     Xm(dim); // X mean
  Vector     Xi(dim);
  Vector     dX(dim);
  Vector     normal(dim);
  Vector     tmp(dim), tmp2(dim);
  Vector     Uf(dim), Ue(dim), Umsh(dim); // Ue := elastic velocity
  int        tag;
  int        iVs[128], *iVs_end;
  VectorXi   vtx_dofs_umesh(dim);  // indices de onde pegar a velocidade
  VectorXi   vtx_dofs_fluid(dim); // indices de onde pegar a velocidade
  VectorXi edge_dofs_umesh(3*dim);
  VectorXi edge_dofs_fluid((2+u_has_edge_assoc_dof)*dim);
  VectorXi   edge_nodes(3);
  Tensor     R(dim,dim);
  bool       in_boundary;
  int        id;

  //VecGetArray(x_mesh, &x_mesh_array);

  //
  // Initial guess for laplacian smoothing.
  // Note that LS is done at timestep n+1.
  //
  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    if (!mesh->isVertex(&*point))
        continue;
    dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), &*point);
    dof_handler_vars.getVariable(0).getVertexDofs(vtx_dofs_fluid.data(), &*point);
    VecGetValues(q, dim, vtx_dofs_fluid.data(), Uf.data());
    id = mesh->getPointId(&*point);
    getRotationMatrix(R,&id,1);
    point->getCoord(Xi.data());
    Xi += R.transpose()*Uf*dt;
    VecSetValues(x_mesh, dim, vtx_dofs_umesh.data(), Xi.data(), INSERT_VALUES);
  }
  if (dim==2 && u_has_edge_assoc_dof)
  {
    facet_iterator edge = mesh->facetBegin();
    facet_iterator edge_end = mesh->facetEnd();
    for (; edge != edge_end; ++edge)
    {
      mesh->getFacetNodesId(&*edge, edge_nodes.data());
      dof_handler_mesh.getVariable(0).getFacetAssociatedDofs(edge_dofs_umesh.data(), &*edge);
      dof_handler_vars.getVariable(0).getFacetAssociatedDofs(edge_dofs_fluid.data(), &*edge);
      VecGetValues(q, dim, edge_dofs_fluid.data(), Uf.data());
      mesh->getNode(edge_nodes(2))->getCoord(Xi.data());
      getRotationMatrix(R,&edge_nodes(2),1);
      Xi += R.transpose()*Uf*dt;
      VecSetValues(x_mesh, dim, edge_dofs_umesh.data(), Xi.data(), INSERT_VALUES);
    }
  }
  else if (dim==3 && u_has_edge_assoc_dof)
  {
    corner_iterator edge = mesh->cornerBegin();
    corner_iterator edge_end = mesh->cornerEnd();
    for (; edge != edge_end; ++edge)
    {
      mesh->getCornerNodesId(&*edge, edge_nodes.data());
      dof_handler_mesh.getVariable(0).getCornerAssociatedDofs(edge_dofs_umesh.data(), &*edge);
      dof_handler_vars.getVariable(0).getCornerAssociatedDofs(edge_dofs_fluid.data(), &*edge);
      VecGetValues(q, dim, edge_dofs_fluid.data(), Uf.data());
      mesh->getNode(edge_nodes(2))->getCoord(Xi.data());
      getRotationMatrix(R,&edge_nodes(2),1);
      Xi += R.transpose()*Uf*dt;
      VecSetValues(x_mesh, dim, edge_dofs_umesh.data(), Xi.data(), INSERT_VALUES);
    }
  }

  Assembly(x_mesh);

  /* suavização laplaciana */
  point = mesh->pointBegin();
  point_end = mesh->pointEnd();
  for (int smooth_it = 0; smooth_it < 1; ++smooth_it)
  {
    getNormalsFromMesh(&x_mesh);

    for (; point != point_end; ++point)
    {

      in_boundary = mesh->inBoundary(&*point);

      if (!boundary_smoothing && in_boundary)
        continue;

      tag = point->getTag();

      if (mesh->isVertex(&*point))
      {
        //if (  is_in(tag,interface_tags) || is_in(tag,triple_tags) || is_in(tag,solid_tags) ||
        //    is_in(tag,dirichlet_tags) || is_in(tag,neumann_tags)  )
        ////if (is_in(tag,triple_tags))
        //  continue;

        Xm = Vector::Zero(dim);
        iVs_end = mesh->connectedVtcs(&*point, iVs);

        if (!in_boundary)
        {
          for (int *it = iVs; it != iVs_end ; ++it)
          {
            dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), mesh->getNode(*it));
            // debug
            VecGetValues(x_mesh, dim, vtx_dofs_umesh.data(), tmp.data());
            Xm += tmp;
          }
          Xm /= (iVs_end-iVs);
          dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), &*point);
          VecSetValues(x_mesh, dim, vtx_dofs_umesh.data(), Xm.data(), INSERT_VALUES);
        }
        else
        {
          int N=0;
          for (int *it = iVs; it != iVs_end ; ++it)
          {
            if (!mesh->inBoundary((mesh->getNode(*it))))
              continue;
            dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), mesh->getNode(*it));
            // debug
            VecGetValues(x_mesh, dim, vtx_dofs_umesh.data(), tmp.data());
            ++N;
            Xm += tmp;
          }
          dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), &*point);
          VecGetValues(x_mesh, dim, vtx_dofs_umesh.data(), Xi.data());

          if (dim==3)
            Xm = (N*Xi + 2*Xm)/(3*N);
            //Xm = Xm/N;
          else
            Xm = (N*Xi + Xm)/(2*N);

          dX = Xm - Xi;
          VecGetValues(nml_mesh, dim, vtx_dofs_umesh.data(), normal.data());
          dX -= normal.dot(dX)*normal;
          Xi += 0.1*dX;
          
          VecSetValues(x_mesh, dim, vtx_dofs_umesh.data(), Xi.data(), INSERT_VALUES);
        }

      }
      //
      // mid node
      //
      else
      {
        if (!in_boundary)
          continue;

        const int m = point->getPosition() - mesh->numVerticesPerCell();
        Cell const* icell = mesh->getCell(point->getIncidCell());
        if (dim==3)
        {
          Corner *edge = mesh->getCorner(icell->getCornerId(m));
          dof_handler_mesh.getVariable(0).getCornerDofs(edge_dofs_umesh.data(), &*edge);
          dof_handler_vars.getVariable(0).getCornerDofs(edge_dofs_fluid.data(), &*edge);
          mesh->getCornerNodesId(&*edge, edge_nodes.data());
        }
        else // dim=2
        {
          Facet *edge = mesh->getFacet(icell->getFacetId(m));
          dof_handler_mesh.getVariable(0).getFacetDofs(edge_dofs_umesh.data(), &*edge);
          dof_handler_vars.getVariable(0).getFacetDofs(edge_dofs_fluid.data(), &*edge);
          mesh->getFacetNodesId(&*edge, edge_nodes.data());
        }

        VecGetValues(x_mesh, dim, edge_dofs_umesh.data(), Xm.data());    // Umsh0
        VecGetValues(x_mesh, dim, edge_dofs_umesh.data()+dim, tmp.data());    // Umsh0
        VecGetValues(x_mesh, dim, edge_dofs_umesh.data()+2*dim, Xi.data());    // Umsh0

        //Xm = (Xm+tmp+2*Xi)/4.;
        Xm = (Xm+tmp)/2.;

        dX = Xm - Xi;
        VecGetValues(nml_mesh, dim, edge_dofs_umesh.data()+2*dim, normal.data());
        dX -= normal.dot(dX)*normal;
        Xi += dX;
        VecSetValues(x_mesh, dim, edge_dofs_umesh.data()+2*dim, Xi.data(), INSERT_VALUES);
        Assembly(x_mesh);

      }

    } // end point
  } // end smooth
  Assembly(x_mesh);
  


  /* calculando a velocidade da malha*/
  point = mesh->pointBegin();
  point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    nodeid = mesh->getPointId(&*point);
    in_boundary = mesh->inBoundary(&*point);
    tag = point->getTag();

    //
    //  vertex
    //
    if (mesh->isVertex(&*point))
    {
      dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), &*point);
      dof_handler_vars.getVariable(0).getVertexDofs(vtx_dofs_fluid.data(), &*point);

      // pega a velocidade do fluido
      VecGetValues(q, dim, vtx_dofs_fluid.data(), Uf.data());
      // pega a coordenada suavizada
      VecGetValues(x_mesh, dim, vtx_dofs_umesh.data(), Xm.data());
      // pega a coordenada atual
      point->getCoord(Xi.data());


      if ((!boundary_smoothing && in_boundary) || is_in(tag,triple_tags))
      {
        // se esta na sup., a vel. da malha é igual a do fluido
        Umsh = Uf;
        VecSetValues(u_mesh, dim, vtx_dofs_umesh.data(), Umsh.data(),INSERT_VALUES);
        continue;
      }
      Ue = (Xm - Xi)/dt;


      getRotationMatrix(R,&nodeid,1);

      //mesh.getNode(i)->setCoord(Xi + dt*(bet1*U+bet2*Ue));

      if (in_boundary)
        Umsh = R*Ue;
      else
        Umsh = beta1*Uf + beta2*Ue;

      VecSetValues(u_mesh, dim, vtx_dofs_umesh.data(), Umsh.data(), INSERT_VALUES);
    }
    else
    //
    // mid node
    //
    {
      const int m = point->getPosition() - mesh->numVerticesPerCell();
      Cell const* icell = mesh->getCell(point->getIncidCell());
      if (dim==3)
      {
        Corner *edge = mesh->getCorner(icell->getCornerId(m));
        dof_handler_mesh.getVariable(0).getCornerDofs(edge_dofs_umesh.data(), &*edge);
        dof_handler_vars.getVariable(0).getCornerDofs(edge_dofs_fluid.data(), &*edge);
        mesh->getCornerNodesId(&*edge, edge_nodes.data());
      }
      else // dim=2
      {
        Facet *edge = mesh->getFacet(icell->getFacetId(m));
        dof_handler_mesh.getVariable(0).getFacetDofs(edge_dofs_umesh.data(), &*edge);
        dof_handler_vars.getVariable(0).getFacetDofs(edge_dofs_fluid.data(), &*edge);
        mesh->getFacetNodesId(&*edge, edge_nodes.data());
      }

      VecGetValues(q, dim, edge_dofs_fluid.data()+2*dim, Uf.data());
      // umesh = ufluid
      if ((!boundary_smoothing && in_boundary) || is_in(tag,triple_tags))
      {
        // se esta na sup., a vel. da malha é igual a do fluido
        Umsh = Uf;
        VecSetValues(u_mesh, dim, edge_dofs_fluid.data()+2*dim, Umsh.data(),INSERT_VALUES);
        continue;
      }

      // se está na interface
      if (in_boundary )
      //if (is_in(tag,triple_tags) )
      {
        point->getCoord(Xi.data());
        VecGetValues(x_mesh, dim, edge_dofs_umesh.data()+2*dim, Xm.data());
        Ue = (Xm - Xi)/dt;
        getRotationMatrix(R,&nodeid,1);
        Umsh = R*Ue;
        VecSetValues(u_mesh, dim, edge_dofs_umesh.data()+2*dim, Umsh.data(),INSERT_VALUES);
      }
      else // se está no interior do domínio, arrasta o nó para o centro da aresta
      {

        VecGetValues(x_mesh, dim, edge_dofs_umesh.data(), Xm.data());    // Umsh0
        VecGetValues(x_mesh, dim, edge_dofs_umesh.data()+dim, Xi.data());    // Umsh0

        Xm = (Xm+Xi)/2.;

        point->getCoord(Xi.data());

        Umsh = (Xm-Xi)/dt;

        getRotationMatrix(R,edge_nodes.data()+2,1);                 // mid node
        Umsh = R*Umsh;

        VecSetValues(u_mesh, dim, edge_dofs_umesh.data()+2*dim, Umsh.data(), INSERT_VALUES);
      }
    }


  } // end point loop


  Assembly(u_mesh);
  //VecRestoreArray(x_mesh, &x_mesh_array);

  PetscFunctionReturn(0);
}



//
// MOVE MESH
//
PetscErrorCode AppCtx::moveMesh()
{
  Vector      Xi(dim), tmp(dim);
  Vector      U(dim), Umsh(dim);
  VectorXi    vtx_dofs_umesh(3); // indices de onde pegar a velocidade
  //int       vtx_dofs_fluid[3]; // indices de onde pegar a velocidade
  //int       tag;
  VectorXi    edge_dofs(dim);
  VectorXi    edge_nodes(3);
  Tensor      R(dim,dim);
  int         nodeid;
  //int         edge_id;
  //Cell const* cell;

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    if (mesh->isVertex(&*point))
    {
      dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), &*point);

      VecGetValues(u_mesh, dim, vtx_dofs_umesh.data(), Umsh.data());

      nodeid = mesh->getPointId(&*point);
      getRotationMatrix(R, &nodeid, 1);

      point->getCoord(Xi.data());
      tmp = Xi + dt*R.transpose()*Umsh;
      point->setCoord(tmp.data() );
    }
  }

  // higher order nodes
  bool const mesh_has_edge_node = mesh->numNodesPerCell()-mesh->numVerticesPerCell() > 0;
  if (mesh_has_edge_node)
  {
    //Point * point;

    if (dim==3) // MESMA COISA, trocando corner <-> corner
    {
      corner_iterator edge = mesh->cornerBegin();
      corner_iterator edge_end = mesh->cornerEnd();
      for (; edge != edge_end; ++edge)
      {
        mesh->getCornerNodesId(&*edge, edge_nodes.data());

        // mid node
        point = point_iterator(mesh.get(),mesh->getNode(edge_nodes[2]));
        //dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh, point);

        dof_handler_mesh.getVariable(0).getCornerAssociatedDofs(edge_dofs.data(), &*edge);

        VecGetValues(u_mesh, dim, edge_dofs.data(), Umsh.data());

        nodeid = mesh->getPointId(&*point);
        getRotationMatrix(R, &nodeid, 1);

        point->getCoord(Xi.data());
        tmp = Xi + dt*R.transpose()*Umsh;
        //tmp = Xi + dt*R*Umsh;
        point->setCoord(tmp.data() );

      } // end loop
    } // end dim==3
    else
    if (dim==2)
    {
      facet_iterator edge = mesh->facetBegin();
      facet_iterator edge_end = mesh->facetEnd();
      for (; edge != edge_end; ++edge)
      {
        mesh->getFacetNodesId(&*edge, edge_nodes.data());

        // mid node
        point = point_iterator(mesh.get(),mesh->getNode(edge_nodes[2]));
        //dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh, point);

        dof_handler_mesh.getVariable(0).getFacetAssociatedDofs(edge_dofs.data(), &*edge);

        VecGetValues(u_mesh, dim, edge_dofs.data(), Umsh.data());

        nodeid = mesh->getPointId(&*point);
        getRotationMatrix(R, &nodeid, 1);

        point->getCoord(Xi.data());
        tmp = Xi + dt*R.transpose()*Umsh;
        //tmp = Xi + dt*R*Umsh;
        point->setCoord(tmp.data() );

      } // end loop
    } // end dim==2


  }

  getNormalsFromMesh(NULL);

  PetscFunctionReturn(0);
}


void AppCtx::freePetscObjs()
{
  Destroy(Jac);
  Destroy(q);
  Destroy(res);
  Destroy(q0);
  Destroy(u_mesh);
  Destroy(x_mesh);
  Destroy(nml_mesh);
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
PetscErrorCode FormJacobian(SNES snes,Vec x,Mat *Jac, Mat *prejac, MatStructure *flag, void *ptr)
{
  AppCtx *user    = static_cast<AppCtx*>(ptr);
  user->formJacobian(snes,x,Jac,prejac,flag);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ptr)
{
  AppCtx *user    = static_cast<AppCtx*>(ptr);
  user->formFunction(snes,x,f);
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

