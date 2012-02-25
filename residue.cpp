#include "common.hpp"

// ******************************************************************************
//                            FORM FUNCTION
// ******************************************************************************
PetscErrorCode AppCtx::formFunction(SNES /*snes*/, Vec Vec_up_k, Vec Vec_fun)
{
  
  int iter;
  
  SNESGetIterationNumber(snes,&iter);
  
  if (!iter)
  {
    converged_times=0;
  }

  //PetscErrorCode      ierr;

  // LOOP NAS CÉLULAS
  VecZeroEntries(Vec_fun);
  //#pragma omp parallel default(none) shared(Vec_up_k,Vec_fun,cout)
  {
    VectorXd          FUloc(n_dofs_u_per_cell);
    VectorXd          FPloc(n_dofs_p_per_cell);

    MatrixXd          u_coefs_c_new(n_dofs_u_per_cell/dim, dim);
    VectorXd          p_coefs_c(n_dofs_p_per_cell);

    VectorXi          mapU_c(n_dofs_u_per_cell);
    VectorXi          mapP_c(n_dofs_p_per_cell);

    //const int tid = omp_get_thread_num();
    //const int nthreads = omp_get_num_threads();
    //const int n_cell_colors = mesh->numCellColors();

    //cell_iterator cell;
    //cell_iterator cell_end;

    //for (int color = 0; color < n_cell_colors; ++color)
    //{
      //cell = mesh->cellBegin(EColor(color),tid,nthreads);
      //cell_end = mesh->cellEnd(EColor(color),tid,nthreads);

      cell_iterator cell = mesh->cellBegin();
      cell_iterator cell_end = mesh->cellEnd();
      for (; cell != cell_end; ++cell)
      {

        // mapeamento do local para o global:
        //
        dof_handler[DH_UNKS].getVariable(VAR_U).getCellDofs(mapU_c.data(), &*cell);
        dof_handler[DH_UNKS].getVariable(VAR_P).getCellDofs(mapP_c.data(), &*cell);

        /*  Pega os valores das variáveis nos graus de liberdade */
        VecGetValues(Vec_up_k , mapU_c.size(), mapU_c.data(), u_coefs_c_new.data());
        VecGetValues(Vec_up_k , mapP_c.size(), mapP_c.data(), p_coefs_c.data());


        formCellFunction(cell, mapU_c, mapP_c, u_coefs_c_new, p_coefs_c, FUloc, FPloc);


        VecSetValues(Vec_fun, mapU_c.size(), mapU_c.data(), FUloc.data(), ADD_VALUES);
        VecSetValues(Vec_fun, mapP_c.size(), mapP_c.data(), FPloc.data(), ADD_VALUES);

      }

      //#pragma omp barrier
    //}


  }

  // LOOP NAS FACES DO CONTORNO
  //#pragma omp parallel default(none) shared(Vec_up_k,Vec_fun,cout)
  {
    VectorXd          FUloc(n_dofs_u_per_facet);
    //VectorXd          FPloc(n_dofs_p_per_facet); // don't need it

    MatrixXd             u_coefs_f(n_dofs_u_per_facet/dim, dim);
    VectorXd             p_coefs_f(n_dofs_p_per_facet);

    VectorXi               mapU_f(n_dofs_u_per_facet);
    VectorXi               mapP_f(n_dofs_p_per_facet);

    bool is_neumann;
    bool is_surface;
    bool is_solid;
    int tag;

    //const int tid = omp_get_thread_num();
    //const int nthreads = omp_get_num_threads();
    //const int n_facet_colors = mesh->numFacetColors();

    //facet_iterator facet;
    //facet_iterator facet_end;

    //for (int color = 0; color < n_facet_colors; ++color)
    //{
      //facet = mesh->facetBegin(EColor(color),tid,nthreads);
     //facet_end = mesh->facetEnd(EColor(color),tid,nthreads);

      // LOOP NAS FACES DO CONTORNO
      facet_iterator facet = mesh->facetBegin();
      facet_iterator facet_end = mesh->facetEnd();
      if (neumann_tags.size() != 0 || interface_tags.size() != 0 || solid_tags.size() != 0)
      for (; facet != facet_end; ++facet)
      {
        tag = facet->getTag();
        is_neumann = is_in(tag, neumann_tags);
        is_surface = is_in(tag, interface_tags);
        is_solid   = is_in(tag, solid_tags);

        if ((!is_neumann) && (!is_surface) && (!is_solid))
        //PetscFunctionReturn(0);
          continue;

        // mapeamento do local para o global:
        //
        dof_handler[DH_UNKS].getVariable(VAR_U).getFacetDofs(mapU_f.data(), &*facet);
        dof_handler[DH_UNKS].getVariable(VAR_P).getFacetDofs(mapP_f.data(), &*facet);

        VecGetValues(Vec_up_k , mapU_f.size(), mapU_f.data(), u_coefs_f.data());

        formFacetFunction(facet, mapU_f, mapP_f, u_coefs_f, p_coefs_f,FUloc);

        VecSetValues(Vec_fun, mapU_f.size(), mapU_f.data(), FUloc.data(), ADD_VALUES);

      }

      //#pragma omp barrier
    //} // end color

  } // end parallel


  // LINHA DE CONTATO
  //#pragma omp parallel shared(Vec_up_k,Vec_fun,cout) default(none)
  {
    int                 tag;
    bool                is_triple;

    MatrixXd             u_coefs_r(n_dofs_u_per_corner/dim, dim);
    VectorXd          FUloc(n_dofs_u_per_corner);

    VectorXi               mapU_r(n_dofs_u_per_corner);
    VectorXi               mapP_r(n_dofs_p_per_corner);

    //const int tid = omp_get_thread_num();
    //const int nthreads = omp_get_num_threads();
    //const int n_corner_colors = mesh->numCornerColors();


    // LOOP NAS ARESTAS DA LINHA TRIPLICE
    corner_iterator corner = mesh->cornerBegin();
    corner_iterator corner_end = mesh->cornerEnd();
    if (triple_tags.size() != 0)
    for (; corner != corner_end; ++corner)
    {
      tag = corner->getTag();
      is_triple = is_in(tag,triple_tags);
      if (!is_triple)
        continue;

      // mapeamento do local para o global:
      //
      dof_handler[DH_UNKS].getVariable(VAR_U).getCornerDofs(mapU_r.data(), &*corner);
      dof_handler[DH_UNKS].getVariable(VAR_P).getCornerDofs(mapP_r.data(), &*corner);

      VecGetValues(Vec_up_k , mapU_r.size(), mapU_r.data(), u_coefs_r.data());

      formCornerFunction(corner,mapU_r,mapP_r,u_coefs_r,FUloc);

      VecSetValues(Vec_fun, mapU_r.size(), mapU_r.data(), FUloc.data(), ADD_VALUES);
      //cout << FUloc.transpose() << endl;

    }



  }


  Assembly(Vec_fun);

  //#pragma omp parallel default(none) shared(Vec_up_k,Vec_fun,cout)
  {
    int       tag;
    Vector    X(dim);
    Vector    Ue(dim);
    Vector    U(dim);
    Vector    dU(dim);
    VectorXi  nodes_dof_fluid(dim);
    VectorXi  nodes_dof_mesh(dim);

    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    for (; point != point_end; ++point)
    {
      tag = point->getTag();
      
      if (!is_in(tag,dirichlet_tags))
        continue;
      
      getNodeDofs(&*point, DH_UNKS, VAR_U, nodes_dof_fluid.data());
      getNodeDofs(&*point, DH_MESH, VAR_M, nodes_dof_mesh.data());
      
      VecGetValues(Vec_up_k, dim, nodes_dof_fluid.data(), U.data());
      VecGetValues(Vec_x_1,  dim, nodes_dof_mesh.data(), X.data());
      
      Ue = u_exact(X,current_time+dt,tag);
      
      dU = U - Ue;
      
      VecSetValues(Vec_fun, dim, nodes_dof_fluid.data(), dU.data(), INSERT_VALUES);
    }

  } // end parallel

  if (force_pressure)
  { 
    int idx;
    double p, pe = 0;
    
    idx = dir_entries.back();
    VecGetValues(Vec_up_k, 1, &idx, &p);
    VecSetValue(Vec_fun, idx, p - pe, INSERT_VALUES);
    Assembly(Vec_fun);
  }

  if(print_to_matlab)
  {
    static bool ja_foi=false;
    if (!ja_foi) View(Vec_fun, "rhs.m","res");
    ja_foi = true;
  }

  Assembly(Vec_fun);

  PetscFunctionReturn(0);

} // END formFunction




// ***********
// form the residue of the cell
void AppCtx::formCellFunction(cell_iterator &cell,
                              VectorXi &mapU_c,  VectorXi &mapP_c, // mappers
                              MatrixXd &u_coefs_c_new,  VectorXd &p_coefs_c_new, // coefficients
                              VectorXd &FUloc, VectorXd &FPloc) // output: local residue
{

  /* local data */
  int                 tag;
  //MatrixXd            u_coefs_c(n_dofs_u_per_cell/dim, dim);
  MatrixXd            u_coefs_c_mid_trans(dim, n_dofs_u_per_cell/dim); // n+utheta
  MatrixXd            u_coefs_c_old(n_dofs_u_per_cell/dim, dim);        // n
  MatrixXd            u_coefs_c_old_trans(dim,n_dofs_u_per_cell/dim);   // n
  //MatrixXd           u_coefs_c_new(n_dofs_u_per_cell/dim, dim);       // n+1
  MatrixXd            u_coefs_c_new_trans(dim,n_dofs_u_per_cell/dim);   // n+1

  MatrixXd            v_coefs_c_trans(dim, nodes_per_cell);      // mesh velocity; n+utheta
  MatrixXd            v_coefs_c_mid(nodes_per_cell, dim);       // mesh velocity; n
  MatrixXd            v_coefs_c_mid_trans(dim,nodes_per_cell);  // mesh velocity; n

  //VectorXd           p_coefs_c_new(n_dofs_p_per_cell);  // n+1
  VectorXd            p_coefs_c_old(n_dofs_p_per_cell);    // n
  VectorXd            p_coefs_c_mid(n_dofs_p_per_cell);   // n+utheta
  
  //MatrixXd            x_coefs_c_mid(nodes_per_cell, dim);       // n+utheta
  MatrixXd            x_coefs_c_mid_trans(dim, nodes_per_cell); // n+utheta
  MatrixXd            x_coefs_c_new(nodes_per_cell, dim);       // n+1
  MatrixXd            x_coefs_c_new_trans(dim, nodes_per_cell); // n+1
  MatrixXd            x_coefs_c_old(nodes_per_cell, dim);       // n
  MatrixXd            x_coefs_c_old_trans(dim, nodes_per_cell); // n
  
  Tensor              F_c_mid(dim,dim);       // n+utheta
  Tensor              invF_c_mid(dim,dim);    // n+utheta
  Tensor              invFT_c_mid(dim,dim);   // n+utheta
  Tensor              F_c_new(dim,dim);       // n+1
  Tensor              invF_c_new(dim,dim);    // n+1
  Tensor              invFT_c_new(dim,dim);   // n+1
  Tensor              F_c_old(dim,dim);       // n
  Tensor              invF_c_old(dim,dim);    // n
  Tensor              invFT_c_old(dim,dim);   // n
  
  /* All variables are in (n+utheta) by default */
  
  MatrixXd            dxphi_c(n_dofs_u_per_cell/dim, dim);
  MatrixXd            dxpsi_c(n_dofs_p_per_cell, dim);
  MatrixXd            dxqsi_c(nodes_per_cell, dim);
  Vector              dxbble(dim);
  Tensor              dxU(dim,dim);   // grad u
  Tensor              dxUb(dim,dim);  // grad u bble
  Vector              dxP(dim);   // grad p
  Vector              Xqp(dim);
  Vector              Xc(dim);  // cell center; to compute CR element
  Vector              Uqp(dim);
  Vector              Uqp_new(dim);
  Vector              Uqp_old(dim);
  Vector              Ubqp(dim); // bble
  Vector              Vqp(dim);
  Vector              Uconv_qp(dim);
  Vector              dUdt(dim);
  double              Pqp=0;
  double              bble_integ=0;

  VectorXi            cell_nodes(nodes_per_cell);
  double              Jx_mid, Jx_new, Jx_old;
  double              weight=0;
  double              visc=-1; // viscosity
  double              cell_volume=0;
  double              hk2=0;
  double              tauk=0;
  double              delk=0;
  double              delta_cd=0;
  double              rho;

  //VectorXd          FUloc(n_dofs_u_per_cell); // subvetor da função f (parte de U)
  //VectorXd          FPloc(n_dofs_p_per_cell);     // subvetor da função f (parte de P)
  Tensor              iBbb(dim, dim);                               // BC, i : inverse ..it is not the inverse to CR element
  MatrixXd            Bnb(n_dofs_u_per_cell, dim);
  MatrixXd            Gbp(dim, n_dofs_p_per_cell);
  MatrixXd            Gnx(n_dofs_u_per_cell, dim);                  // CR ;; suffix x means p gradient
  Vector              FUb(dim);
  Vector              FPx(dim); // pressure gradient

  Vector              force_at_theta(dim);
  Vector              Res(dim);                                     // residue
  Tensor const        I(Tensor::Identity(dim,dim));
  Vector              vec(dim);     // temp
  Tensor              Ten(dim,dim); // temp  

  VectorXi            mapM_c(dim*nodes_per_cell); // mesh velocity



  // ----- computations ------

  tag = cell->getTag();

  dof_handler[DH_MESH].getVariable(VAR_M).getCellDofs(mapM_c.data(), &*cell);

  VecGetValues(Vec_v_mid,   mapM_c.size(), mapM_c.data(), v_coefs_c_mid.data());
  VecGetValues(Vec_x_0,     mapM_c.size(), mapM_c.data(), x_coefs_c_old.data());
  VecGetValues(Vec_x_1,     mapM_c.size(), mapM_c.data(), x_coefs_c_new.data());
  VecGetValues(Vec_up_0,    mapU_c.size(), mapU_c.data(), u_coefs_c_old.data());
  VecGetValues(Vec_up_0,    mapP_c.size(), mapP_c.data(), p_coefs_c_old.data());
  
  // get nodal coordinates of the old and new cell
  mesh->getCellNodesId(&*cell, cell_nodes.data());
  //mesh->getNodesCoords(cell_nodes.begin(), cell_nodes.end(), x_coefs_c.data());
  //x_coefs_c_trans = x_coefs_c.transpose();
  
  v_coefs_c_mid_trans = v_coefs_c_mid.transpose();
  x_coefs_c_old_trans = x_coefs_c_old.transpose();
  x_coefs_c_new_trans = x_coefs_c_new.transpose();
  u_coefs_c_old_trans = u_coefs_c_old.transpose();
  u_coefs_c_new_trans = u_coefs_c_new.transpose();

  u_coefs_c_mid_trans = utheta*u_coefs_c_new_trans + (1.-utheta)*u_coefs_c_old_trans;
  x_coefs_c_mid_trans = utheta*x_coefs_c_new_trans + (1.-utheta)*x_coefs_c_old_trans;
  p_coefs_c_mid       = utheta*p_coefs_c_new       + (1.-utheta)*p_coefs_c_old;

  visc = muu(current_time, tag);
  rho  = pho(Xqp,tag);
  FUloc.setZero();
  FPloc.setZero();
  
  /*
  if (behaviors & BH_bble_condens_PnPn)
  {
    iBbb.setZero();
    Bnb.setZero();
    Gbp.setZero();
    FUb.setZero();
  }

  if(behaviors & BH_GLS)
  {
    cell_volume = 0;
    for (int qp = 0; qp < n_qpts_cell; ++qp) {
      F_c_mid = x_coefs_c.transpose() * dLqsi_c[qp];
      Jx_mid = determinant(F_c_mid,dim);
      cell_volume += Jx_mid * quadr_cell->weight(qp);
    }

    hk2 = cell_volume / pi; // element size
    double const uconv = (u_coefs_c_old - v_coefs_c_mid).lpNorm<Infinity>();
    
    tauk = 4.*visc/hk2 + 2.*rho*uconv/sqrt(hk2);
    tauk = 1./tauk;
    if (dim==3)
      tauk *= 0.1;
    
    delk = 4.*visc + 2.*rho*uconv*sqrt(hk2);

  }
  if (behaviors & BH_bble_condens_CR)
  {
    bble_integ = 0;
    Gnx.setZero();
    iBbb.setZero();
    Bnb.setZero();
    FUb.setZero();
    FPx.setZero();

    cell_volume = 0;
    Xc.setZero();
    for (int qp = 0; qp < n_qpts_cell; ++qp) {
      F_c_mid = x_coefs_c_trans * dLqsi_c[qp];
      Jx_mid = determinant(F_c_mid,dim);
      Xqp  = x_coefs_c_trans * qsi_c[qp];
      cell_volume += Jx_mid * quadr_cell->weight(qp);
      Xc += Jx_mid * quadr_cell->weight(qp) * Xqp;
    }
    Xc /= cell_volume;
  }
  */

  // Quadrature
  for (int qp = 0; qp < n_qpts_cell; ++qp)
  {
    //F_c_new = x_coefs_c_new_trans * dLqsi_c[qp];
    //inverseAndDet(F_c_new,dim,invF_c_new,Jx_new);
    //invFT_c_new= invF_c_new.transpose();
    //
    //F_c_old = x_coefs_c_old_trans * dLqsi_c[qp];
    //inverseAndDet(F_c_old,dim,invF_c_old,Jx_old);
    //invFT_c_old= invF_c_old.transpose();
    
    F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];
    inverseAndDet(F_c_mid,dim,invF_c_mid,Jx_mid);
    invFT_c_mid= invF_c_mid.transpose();

    dxphi_c = dLphi_c[qp] * invF_c_mid;
    dxpsi_c = dLpsi_c[qp] * invF_c_mid;
    dxqsi_c = dLqsi_c[qp] * invF_c_mid;

    dxP  = dxpsi_c.transpose() * p_coefs_c_mid;
    dxU  = u_coefs_c_mid_trans * dxphi_c;       // n+utheta

    Xqp  = x_coefs_c_mid_trans * qsi_c[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
    Uqp  = u_coefs_c_mid_trans * phi_c[qp]; //n+utheta
    Uqp_new  = u_coefs_c_new_trans * phi_c[qp]; //n+utheta
    Uqp_old  = u_coefs_c_old_trans * phi_c[qp]; //n+utheta
    Pqp  = p_coefs_c_mid.dot(psi_c[qp]);
    Vqp  = v_coefs_c_mid_trans * qsi_c[qp];
    Uconv_qp = Uqp - Vqp;
    //Uconv_qp = Uqp_old;
    dUdt    = (Uqp_new-Uqp_old)/dt;

    force_at_theta = utheta*force(Xqp,current_time+dt,tag) + (1.-utheta)*force(Xqp,current_time,tag);

    weight = quadr_cell->weight(qp);

    if (Jx_mid < 1.e-10)
    {
      //#pragma omp critical
      //if (tid==0)
      {
        std::cout << "erro: jacobiana da integral não invertível: ";
        std::cout << "Jx_mid = " << Jx_mid << endl;
        cout << "trans matrix:\n" << F_c_mid << endl;
        cout << "x coefs mid:" << endl;
        cout << x_coefs_c_mid_trans.transpose() << endl;
        cout << "-----" << endl;
        cout << "cell id: " << mesh->getCellId(&*cell) << endl;
        cout << "cell nodes:\n" << cell_nodes.transpose() << endl;
        cout << "mapM :\n" << mapM_c.transpose() << endl;
        throw;
      }
    }

    for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
    {
      for (int c = 0; c < dim; ++c)
      {
        FUloc(i*dim + c) += Jx_mid*weight*
                ( rho*(dUdt(c) + has_convec*Uconv_qp.dot(dxU.row(c)))*phi_c[qp][i] + // aceleração
                  visc*dxphi_c.row(i).dot(dxU.row(c) + dxU.col(c).transpose()) - //rigidez
                  Pqp*dxphi_c(i,c) - // pressão
                  force_at_theta(c)*phi_c[qp][i]   ); // força

      }
    }
    for (int i = 0; i < n_dofs_p_per_cell; ++i)
      FPloc(i) -= Jx_mid*weight* dxU.trace()*psi_c[qp][i];

    /* ----------------
      * stabilization
      * ---------------- */
    /*
    if (behaviors & (BH_bble_condens_PnPn | BH_bble_condens_CR))
    {
      dxbble = invFT_c_mid * dLbble[qp];
    
      for (int c = 0; c < dim; c++)
      {
        for (int j = 0; j < n_dofs_u_per_cell/dim; j++)
        {
          for (int d = 0; d < dim; d++)
          {
            delta_cd = c==d;
            Bnb(j*dim + d, c) += Jx_mid*weight*
                                 ( has_convec*phi_c[qp][j]*utheta *pho(Xqp,tag)*(  delta_cd*Uconv_qp.dot(dxbble)  )   // convective
                                 + delta_cd*pho(Xqp,tag)*phi_c[qp][j]*bble[qp]/dt     // time derivative
                                 + utheta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) );    // rigidez

          }
        }
        if (behaviors & BH_bble_condens_PnPn)
          for (int j = 0; j < n_dofs_p_per_cell; ++j)
            Gbp(c, j) -= Jx_mid*weight*psi_c[qp][j]*dxbble(c);
      }

      for (int c = 0; c < dim; c++)
      {
        for (int d = 0; d < dim; d++)
        {
          delta_cd = c==d;
          iBbb(c, d) += Jx_mid*weight*
                        ( has_convec*bble[qp]*utheta *rho*( delta_cd*Uconv_qp.dot(dxbble) )   // convective
                        + delta_cd*pho(Xqp,tag)*bble[qp]*bble[qp]/dt     // time derivative
                        + utheta*visc*(delta_cd* dxbble.dot(dxbble)  + dxbble(d)*dxbble(c)) ); // rigidez
        }

        FUb(c) += Jx_mid*weight*
                  ( bble[qp]*rho*(dUdt(c) + has_convec*Uconv_qp.dot(dxU.row(c))) + // time derivative + convective
                    visc*dxbble.dot(dxU.row(c) + dxU.col(c).transpose()) - //rigidez
                    Pqp*dxbble(c) -                     // pressão
                    force_at_theta(c)*bble[qp]   ); // força
      }
    }
    else
    if(behaviors & BH_GLS)
    {
      Res = rho*( dUdt +  has_convec*dxU*Uconv_qp) + dxP - force_at_theta;
      
      for (int i = 0; i < n_dofs_u_per_cell/dim; i++)
      {
        vec = Jx_mid*weight*(  has_convec*tauk*  rho* Uconv_qp.dot(dxphi_c.row(i)) * Res    +   delk*dxU.trace()*dxphi_c.row(i).transpose() );
        
        for (int c = 0; c < dim; c++)
          FUloc(i*dim + c) += vec(c);
          
      }
      for (int i = 0; i < n_dofs_p_per_cell; ++i)
        FPloc(i) -= Jx_mid*weight *tauk* dxpsi_c.row(i).dot(Res);
    }

    if (behaviors & BH_bble_condens_CR)
    {
      bble_integ += Jx_mid*weight*bble[qp];

      for (int c = 0; c < dim; ++c)
      {
        for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
          for (int j = 0; j < dim; ++j) // pressure gradient
            Gnx(i*dim + c,j) -= Jx_mid*weight* (Xqp(j) - Xc(j))*dxphi_c(i,c);

        FPx(c) -= Jx_mid*weight* dxU_new.trace()*(Xqp(c) - Xc(c));
      }
    }
    */
  } // fim quadratura

  
  /* stabilization */
  /*
  if (behaviors & BH_bble_condens_PnPn)
  {
    //iBbb = iBbb.inverse().eval();
    invert(iBbb,dim);
    
    FUloc = FUloc - Bnb*iBbb*FUb;
    FPloc = FPloc - Gbp.transpose()*iBbb*FUb;
  }

  if(behaviors & BH_bble_condens_CR)
  {
    double const a = 1./(bble_integ*bble_integ);
    double const b = 1./bble_integ;
    
    FUloc +=  a*Gnx*iBbb*FPx - b*Bnb*FPx - b*Gnx*FUb;
  }
  */


  //PetscFunctionReturn(0);
  return;
} // end formCellFunction

// ***********
// form the residue of the facet
void AppCtx::formFacetFunction(facet_iterator &facet,
                       VectorXi const&/*mapU_f*/,  VectorXi const&/*mapP_f*/, // mappers
                       MatrixXd &u_coefs_f,  VectorXd &/*p_coefs_f*/, // coefficients
                       VectorXd &FUloc) // output: local residue
{
  int                 tag;
  bool                is_neumann;
  bool                is_surface;
  bool                is_solid;
  //MatrixXd           u_coefs_f(n_dofs_u_per_facet/dim, dim);
  //VectorXd           p_coefs_f(n_dofs_p_per_facet);
  MatrixXd            u_coefs_f_trans(dim, n_dofs_u_per_facet/dim);
  MatrixXd             x_coefs_f(nodes_per_facet, dim);
  MatrixXd            x_coefs_f_trans(dim, nodes_per_facet);
  Tensor              F_f(dim,dim-1);
  Tensor              invF_f(dim-1,dim);
  MatrixXd            dxphi_f(n_dofs_u_per_facet/dim, dim);
  Tensor              dxU_f(dim,dim);   // grad u
  Vector              Xqp(dim);
  Vector              Uqp(dim);
  //VectorXd          FUloc(n_dofs_u_per_facet);
  MatrixXd            Aloc_f(n_dofs_u_per_facet, n_dofs_u_per_facet);
  VectorXi            facet_nodes(nodes_per_facet);
  Vector              normal(dim);
  double              Jx_mid=0;
  double              weight=0;

  VectorXi               mapM_f(dim*nodes_per_facet);

  MatrixXd         R(n_dofs_u_per_facet,n_dofs_u_per_facet);
  MatrixXd         tmp;
  Vector           traction_(dim);

  // ----- computations ------
  FUloc.setZero();

  tag = facet->getTag();

  is_neumann = (neumann_tags.end() != std::find(neumann_tags.begin(), neumann_tags.end(), tag));
  is_surface = (interface_tags.end() != std::find(interface_tags.begin(), interface_tags.end(), tag));
  is_solid = (solid_tags.end() != std::find(solid_tags.begin(), solid_tags.end(), tag));

  if ((!is_neumann) && (!is_surface) && (!is_solid))
    //PetscFunctionReturn(0);
    return;

  mesh->getFacetNodesId(&*facet, facet_nodes.data());
  mesh->getNodesCoords(facet_nodes.begin(), facet_nodes.end(), x_coefs_f.data());
  x_coefs_f_trans = x_coefs_f.transpose();


  
  u_coefs_f_trans = u_coefs_f.transpose();


  for (int qp = 0; qp < n_qpts_facet; ++qp)
  {

    F_f   = x_coefs_f_trans * dLqsi_f[qp];

    tmp.resize(dim-1,dim-1);
    tmp = F_f.transpose()*F_f;
    Jx_mid = sqrt(tmp.determinant());
    tmp = tmp.inverse();
    invF_f = tmp*F_f.transpose();


    weight  = quadr_facet->weight(qp);
    Xqp     = x_coefs_f_trans * qsi_f[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
    dxphi_f = dLphi_f[qp] * invF_f;
    dxU_f   = u_coefs_f_trans * dxphi_f; // n+utheta
    Uqp  = u_coefs_f_trans * phi_f[qp];

    if (is_neumann)
    //traction_ = utheta*(traction(Xqp,current_time+dt,tag)) + (1.-utheta)*traction(Xqp,current_time,tag);
    traction_ = traction(Xqp,current_time + dt,tag);
    //traction_ = (traction(Xqp,current_time,tag) +4.*traction(Xqp,current_time+dt/2.,tag) + traction(Xqp,current_time+dt,tag))/6.;
    
      for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
      {
        for (int c = 0; c < dim; ++c)
        {
          FUloc(i*dim + c) -= Jx_mid*weight * traction_(c) * phi_f[qp][i] ; // força
        }
      }

    if (is_surface)
      for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
      {
        for (int c = 0; c < dim; ++c)
        {
          FUloc(i*dim + c) += Jx_mid*weight *gama(Xqp,current_time,tag)*
                                          (dxphi_f(i,c) + 0*(unsteady*dt) *dxU_f.row(c).dot(dxphi_f.row(i)));
        }
      }

    if (is_solid)
      for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
      {
        for (int c = 0; c < dim; ++c)
        {
          FUloc(i*dim + c) += Jx_mid*weight *beta_diss()*Uqp(c)*phi_f[qp][i];
        }
      }

  } // end quadratura



  //PetscFunctionReturn(0);
  //return;
} // end formFacetFunction


// ***********
// form the residue of the contact line (2D)
void AppCtx::formCornerFunction(corner_iterator &corner,
                        VectorXi const&/*mapU_r*/,  VectorXi const&/*mapP_r*/, // mappers
                        MatrixXd &u_coefs_r, // coefficients
                        VectorXd &FUloc)
{

  bool                gen_error = false;
  int                 tag;
  bool                is_triple;
  //MatrixXd             u_coefs_r(n_dofs_u_per_corner/dim, dim);
  MatrixXd            u_coefs_r_trans(dim, n_dofs_u_per_corner/dim);
  MatrixXd            x_coefs_r(nodes_per_corner, dim);
  MatrixXd            x_coefs_r_trans(dim, nodes_per_corner);
  Tensor              F_r(dim,dim-2);
  Tensor              invF_r(dim-2,dim);
  MatrixXd            dxphi_r(n_dofs_u_per_corner/dim, dim);
  Tensor              dxU_r(dim,dim);   // grad u
  Vector              Xqp(dim);
  Vector              Uqp(dim);
  //VectorXd          FUloc(n_dofs_u_per_corner);
  //MatrixXd         Aloc_r(n_dofs_u_per_corner, n_dofs_u_per_corner);
  VectorXi            corner_nodes(nodes_per_corner);
  Vector              normal(dim);
  Vector              line_normal(dim);
  Vector              solid_point(dim); // ponto na superfície do sólido .. ele é único
  Vector              point_a(dim); // ponto na linha triplice
  Vector              point_b(dim); // ponto na linha triplice
  Vector              ifacet_normal(dim); // ponto na linha triplice
  double              line_normal_sign; // +1 or -1
  double              Jx_mid=0;
  double              weight=0;
  int                 iCs[FEPIC_MAX_ICELLS];
  int                 eiCs[FEPIC_MAX_ICELLS];
  int                 *iCs_end;
  int                 *iCs_it;
  Cell                *fluid_cell;

  //VectorXi               mapU_r(n_dofs_u_per_corner);
  //VectorXi               mapP_r(n_dofs_p_per_corner);
  VectorXi               mapM_r(dim*nodes_per_corner);

  MatrixXd         R(n_dofs_u_per_corner,n_dofs_u_per_corner);
  MatrixXd         tmp;


  tag = corner->getTag();

  is_triple = (triple_tags.end() != std::find(triple_tags.begin(), triple_tags.end(), tag));

  FUloc.setZero();

  if (!is_triple)
    return;

  mesh->getCornerNodesId(&*corner, corner_nodes.data());

  if (dim==3)
  {
    // encontrando o ponto da superfície sólido. com ele, é calculado uma normal
    // que corrigi o sinal de line_normal
    iCs_end = mesh->edgeStar(&*corner, iCs, eiCs);
    if (iCs_end == iCs)
    {
      printf("ERROR!: no icell found\n");
      throw;
    }
    gen_error = true;
    for (iCs_it = iCs; iCs_it != iCs_end ; ++iCs_it)
    {
      fluid_cell = mesh->getCell(*iCs_it);

      for (int kk = 0; kk < mesh->numVerticesPerCell(); ++kk)
      {
        Point const* pp      = mesh->getNode(fluid_cell->getNodeId(kk) );
        const int    tag_aux = pp->getTag();

        if (is_in(tag_aux, solid_tags) && !is_in(tag_aux, triple_tags))
        {
          gen_error = false;
          pp->getCoord(solid_point.data());
          break;
        }
      }
      if (gen_error==false)
        break;

    }
    if (gen_error)
    {
      printf("ERROR!: solid point not found\n");
      throw;
    }


    mesh->getNode(corner_nodes(0))->getCoord(point_a.data());
    mesh->getNode(corner_nodes(1))->getCoord(point_b.data());

    // se (a-c) cross (b-c) dot solid_normal > 0, então line_normal_sign = 1, se não, =-1
    Xqp = (point_a+point_b)/2.;
    point_a -= solid_point;
    point_b -= solid_point;
    normal  = solid_normal(Xqp, current_time, tag);
    line_normal_sign = point_a(0)*point_b(1)*normal(2) + point_a(1)*point_b(2)*normal(0) + point_a(2)*point_b(0)*normal(1)
                      -point_a(0)*point_b(2)*normal(1) - point_a(1)*point_b(0)*normal(2) - point_a(2)*point_b(1)*normal(0);
    if (line_normal_sign>0)
      line_normal_sign = 1;
    else
      line_normal_sign = -1;

  }
  else // dim==2
  {
    Point * point = mesh->getNode(corner_nodes[0]);
    Point * sol_point;
    int iVs[FEPIC_MAX_ICELLS];
    int *iVs_end, *iVs_it;
    Vector aux(dim);

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
    point->getCoord(Xqp.data());

    normal = solid_normal(Xqp, current_time, tag);

    sol_point->getCoord(aux.data());


    // choose a line_normal candidate
    line_normal(0) = -normal(1);
    line_normal(1) =  normal(0);

    // check orientation
    if (line_normal.dot(Xqp-aux) < 0)
      line_normal *= -1;
  }


  //mesh->getCornerNodesId(&*corner, corner_nodes.data());
  mesh->getNodesCoords(corner_nodes.begin(), corner_nodes.end(), x_coefs_r.data());
  x_coefs_r_trans = x_coefs_r.transpose();

  
  u_coefs_r_trans = u_coefs_r.transpose();


  for (int qp = 0; qp < n_qpts_corner; ++qp)
  {
    if (dim==3)
    {
      F_r   = x_coefs_r_trans * dLqsi_r[qp];
      Jx_mid = F_r.norm();
    }
    else
    {
      Jx_mid = 1;
    }
    //invF_r = F_r.transpose()/(Jx_mid*Jx_mid);

    weight  = quadr_corner->weight(qp);
    Xqp = x_coefs_r_trans * qsi_r[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
    //dxphi_r = dLphi_r[qp] * invF_r;
    //dxU_r   = u_coefs_r_trans * dxphi_r; // n+utheta
    Uqp  = u_coefs_r_trans * phi_r[qp];

    if (dim==3)
    {
      normal  = solid_normal(Xqp, current_time, tag);
      line_normal(0)= F_r(0,0);line_normal(1)= F_r(1,0);line_normal(2)= F_r(2,0);
      line_normal *= line_normal_sign;
      // a = a cross b
      cross(line_normal, normal);
      line_normal.normalize();
    }

    for (int i = 0; i < n_dofs_u_per_corner/dim; ++i)
    {
      for (int j = 0; j < dim; ++j)
      {
        FUloc(i*dim + j) += Jx_mid*weight*(-gama(Xqp, current_time, tag)*cos_theta0() + zeta(0,0)*line_normal.dot(Uqp))*line_normal(j)*phi_r[qp][i];
      }
    }


  } // end quadratura

  ////#pragma omp critical
  //{
  //cout << "line_normal " << line_normal.size() << endl;
  //cout << "---------------------------------" << endl;
  //cout << line_normal << endl;
  //cout << "=================================" << endl;
  //cout << n_dofs_u_per_corner/dim << endl;
  //}


}




