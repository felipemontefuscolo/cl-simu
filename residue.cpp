#include "common.hpp"

template <typename Derived>
void getProjectorMatrix(MatrixBase<Derived> & P, int n_nodes, int const* nodes, Vec const& Vec_x_, double t, AppCtx const& app)
{
  int const dim = app.dim;
  Mesh const* mesh = &*app.mesh;
  //DofHandler const* dof_handler = &*app.dof_handler;
  std::vector<int> const& dirichlet_tags  = app.dirichlet_tags;
  //std::vector<int> const& neumann_tags    = app.neumann_tags  ;
  //std::vector<int> const& interface_tags  = app.interface_tags;
  std::vector<int> const& solid_tags      = app.solid_tags    ;
  std::vector<int> const& triple_tags     = app.triple_tags   ;
  //std::vector<int> const& periodic_tags   = app.periodic_tags ;
  std::vector<int> const& feature_tags    = app.feature_tags  ;
  //Vec const& Vec_normal = app.Vec_normal;
  
  P.setIdentity();

  Tensor I(dim,dim);
  Tensor Z(dim,dim);
  Vector X(dim);
  Vector normal(dim);
  int    dofs[dim];
  int    tag;
  Point const* point;

  I.setIdentity();
  Z.setZero();

  // NODES
  for (int i = 0; i < n_nodes; ++i)
  {
    point = mesh->getNodePtr(nodes[i]);
    tag = point->getTag();
    //m = point->getPosition() - mesh->numVerticesPerCell();
    //cell = mesh->getCellPtr(point->getIncidCell());

    if (is_in(tag,feature_tags))
    {
      app.getNodeDofs(&*point, DH_MESH, VAR_M, dofs);
      VecGetValues(Vec_x_, dim, dofs, X.data());
      P.block(i*dim,i*dim,dim,dim)  = feature_proj(X,t,tag);
      continue;
    }
    else
    if (is_in(tag,solid_tags) || is_in(tag, triple_tags))
    {
      app.getNodeDofs(&*point, DH_MESH, VAR_M, dofs);
      VecGetValues(Vec_x_, dim, dofs, X.data());

      normal = solid_normal(X,t,tag);

      P.block(i*dim,i*dim,dim,dim)  = I - normal*normal.transpose();      
      //P.block(i*dim,i*dim,dim,dim) = Z;

    }
    else
    if (is_in(tag,dirichlet_tags))
    {
      P.block(i*dim,i*dim,dim,dim) = Z;
    }


  } // end nodes
}



// ******************************************************************************
//                            FORM FUNCTION
// ******************************************************************************
PetscErrorCode AppCtx::formFunction(SNES /*snes*/, Vec Vec_up_k, Vec Vec_fun)
{
  //PetscErrorCode      ierr;

  int null_space_press_dof=-1;
  
  if (force_pressure)
  {
    if (behaviors & BH_Press_grad_elim)
    {
      cell_iterator cell = mesh->cellBegin();
      dof_handler[DH_UNKS].getVariable(VAR_P).getCellDofs(&null_space_press_dof, &*cell);
    }
    else
    {
      point_iterator point = mesh->pointBegin();
      while (!mesh->isVertex(&*point))
        ++point;
      dof_handler[DH_UNKS].getVariable(VAR_P).getVertexDofs(&null_space_press_dof, &*point);
    }
    // fix the initial guess
    VecSetValue(Vec_up_k, null_space_press_dof, 0.0, INSERT_VALUES);
    Assembly(Vec_up_k);
  }


  Mat *JJ = &Mat_Jac;

  int iter;

  SNESGetIterationNumber(snes,&iter);

  if (!iter)
  {
    converged_times=0;
  }

  //PetscErrorCode      ierr;

  // LOOP NAS CÉLULAS
  VecZeroEntries(Vec_fun);
  MatZeroEntries(*JJ);
  
#ifdef FEP_HAS_OPENMP
  FEP_PRAGMA_OMP(parallel default(none) shared(Vec_up_k,Vec_fun,cout))
#endif
  {
    VectorXd          FUloc(n_dofs_u_per_cell);
    VectorXd          FPloc(n_dofs_p_per_cell);

    /* local data */
    int                 tag;
    MatrixXd            u_coefs_c_mid_trans(dim, n_dofs_u_per_cell/dim);  // n+utheta
    MatrixXd            u_coefs_c_old(n_dofs_u_per_cell/dim, dim);        // n
    MatrixXd            u_coefs_c_old_trans(dim,n_dofs_u_per_cell/dim);   // n
    MatrixXd            u_coefs_c_new(n_dofs_u_per_cell/dim, dim);        // n+1
    MatrixXd            u_coefs_c_new_trans(dim,n_dofs_u_per_cell/dim);   // n+1

    MatrixXd            v_coefs_c_trans(dim, nodes_per_cell);      // mesh velocity; n+utheta
    MatrixXd            v_coefs_c_mid(nodes_per_cell, dim);        // mesh velocity; n
    MatrixXd            v_coefs_c_mid_trans(dim,nodes_per_cell);   // mesh velocity; n

    VectorXd            p_coefs_c_new(n_dofs_p_per_cell);  // n+1
    //VectorXd            p_coefs_c_old(n_dofs_p_per_cell);  // n

    //MatrixXd            x_coefs_c_mid(nodes_per_cell, dim);       // n+utheta
    MatrixXd            x_coefs_c_mid_trans(dim, nodes_per_cell); // n+utheta
    MatrixXd            x_coefs_c_new(nodes_per_cell, dim);       // n+1
    MatrixXd            x_coefs_c_new_trans(dim, nodes_per_cell); // n+1
    MatrixXd            x_coefs_c_old(nodes_per_cell, dim);       // n
    MatrixXd            x_coefs_c_old_trans(dim, nodes_per_cell); // n

    Tensor              F_c_mid(dim,dim);       // n+utheta
    Tensor              invF_c_mid(dim,dim);    // n+utheta
    Tensor              invFT_c_mid(dim,dim);   // n+utheta
    //Tensor              F_c_new(dim,dim);       // n+1
    //Tensor              invF_c_new(dim,dim);    // n+1
    //Tensor              invFT_c_new(dim,dim);   // n+1
    //Tensor              F_c_old(dim,dim);       // n
    //Tensor              invF_c_old(dim,dim);    // n
    //Tensor              invFT_c_old(dim,dim);   // n

    /* All variables are in (n+utheta) by default */

    MatrixXd            dxphi_c(n_dofs_u_per_cell/dim, dim);
    MatrixXd            dxpsi_c(n_dofs_p_per_cell, dim);
    MatrixXd            dxqsi_c(nodes_per_cell, dim);
    Vector              dxbble(dim);
    Tensor              dxU(dim,dim);   // grad u
    Tensor              dxUb(dim,dim);  // grad u bble
    Vector              dxP_new(dim);   // grad p
    Vector              Xqp(dim);
    Vector              Xc(dim);  // cell center; to compute CR element
    Vector              Uqp(dim);
    Vector              Ubqp(dim); // bble
    Vector              Uqp_old(dim);  // n
    Vector              Uqp_new(dim);  // n+1
    Vector              Vqp(dim);
    Vector              Uconv_qp(dim);
    Vector              dUdt(dim);
    double              Pqp_new;
    double              bble_integ=0;
    //VectorXd          FUloc(n_dofs_u_per_cell); // subvetor da função f (parte de U)
    //VectorXd          FPloc(n_dofs_p_per_cell);     // subvetor da função f (parte de P)
    VectorXi            cell_nodes(nodes_per_cell);
    double              J_mid;
    //double              J_new, J_old;
    double              JxW_mid;
    double              weight;
    double              visc=-1; // viscosity
    double              cell_volume;
    double              hk2;
    double              tauk=0;
    double              delk=0;
    double              delta_cd;
    double              rho;

    MatrixXd            Aloc(n_dofs_u_per_cell, n_dofs_u_per_cell);
    MatrixXd            Gloc(n_dofs_u_per_cell, n_dofs_p_per_cell);
    MatrixXd            Dloc(n_dofs_p_per_cell, n_dofs_u_per_cell);
    MatrixXd            Eloc(n_dofs_p_per_cell, n_dofs_p_per_cell);   // GSL, BC
    MatrixXd            Cloc(n_dofs_u_per_cell, n_dofs_p_per_cell);   // GSL
    Tensor              iBbb(dim, dim);                               // BC, i : inverse ..it is not the inverse to CR element
    MatrixXd            Bbn(dim, n_dofs_u_per_cell);                  // BC
    MatrixXd            Bnb(n_dofs_u_per_cell, dim);                  // BC
    MatrixXd            Dpb(n_dofs_p_per_cell, dim);                  // BC
    MatrixXd            Gbp(dim, n_dofs_p_per_cell);                  // BC
    MatrixXd            Gnx(n_dofs_u_per_cell, dim);                  // CR ;; suffix x means p gradient
    Vector              FUb(dim);                                     // BC
    Vector              FPx(dim); // pressure gradient

    Vector              force_at_mid(dim);
    Vector              Res(dim);                                     // residue
    Tensor              dResdu(dim,dim);                              // residue derivative
    Tensor const        I(Tensor::Identity(dim,dim));
    Vector              vec(dim);     // temp
    Tensor              Ten(dim,dim); // temp

    VectorXi            mapU_c(n_dofs_u_per_cell);
    VectorXi            mapU_r(n_dofs_u_per_corner);
    VectorXi            mapP_c(n_dofs_p_per_cell);
    VectorXi            mapP_r(n_dofs_p_per_corner);
    // mesh velocity
    VectorXi            mapM_c(dim*nodes_per_cell);
    VectorXi            mapM_f(dim*nodes_per_facet);
    VectorXi            mapM_r(dim*nodes_per_corner);

    MatrixXd            Prj(n_dofs_u_per_cell,n_dofs_u_per_cell); // projector matrix
    //VectorXi            cell_nodes(nodes_per_cell);


    const int tid = omp_get_thread_num();
    const int nthreads = omp_get_num_threads();

    cell_iterator cell = mesh->cellBegin(tid,nthreads);
    cell_iterator cell_end = mesh->cellEnd(tid,nthreads);

    //cell_iterator cell = mesh->cellBegin();
    //cell_iterator cell_end = mesh->cellEnd();
    for (; cell != cell_end; ++cell)
    {

      tag = cell->getTag();

      // mapeamento do local para o global:
      //
      dof_handler[DH_MESH].getVariable(VAR_M).getCellDofs(mapM_c.data(), &*cell);
      dof_handler[DH_UNKS].getVariable(VAR_U).getCellDofs(mapU_c.data(), &*cell);
      dof_handler[DH_UNKS].getVariable(VAR_P).getCellDofs(mapP_c.data(), &*cell);

      VecGetValues(Vec_v_mid,   mapM_c.size(), mapM_c.data(), v_coefs_c_mid.data());
      VecGetValues(Vec_x_0,     mapM_c.size(), mapM_c.data(), x_coefs_c_old.data());
      VecGetValues(Vec_x_1,     mapM_c.size(), mapM_c.data(), x_coefs_c_new.data());
      VecGetValues(Vec_up_0,    mapU_c.size(), mapU_c.data(), u_coefs_c_old.data());
      VecGetValues(Vec_up_k,    mapU_c.size(), mapU_c.data(), u_coefs_c_new.data());
      VecGetValues(Vec_up_k,    mapP_c.size(), mapP_c.data(), p_coefs_c_new.data());

      // get nodal coordinates of the old and new cell
      mesh->getCellNodesId(&*cell, cell_nodes.data());
      //mesh->getNodesCoords(cell_nodes.begin(), cell_nodes.end(), x_coefs_c.data());
      //x_coefs_c_trans = x_coefs_c_mid_trans;

      v_coefs_c_mid_trans = v_coefs_c_mid.transpose();
      x_coefs_c_old_trans = x_coefs_c_old.transpose();
      x_coefs_c_new_trans = x_coefs_c_new.transpose();
      u_coefs_c_old_trans = u_coefs_c_old.transpose();
      u_coefs_c_new_trans = u_coefs_c_new.transpose();

      u_coefs_c_mid_trans = utheta*u_coefs_c_new_trans + (1.-utheta)*u_coefs_c_old_trans;
      x_coefs_c_mid_trans = utheta*x_coefs_c_new_trans + (1.-utheta)*x_coefs_c_old_trans;
      //p_coefs_c_mid       = utheta*p_coefs_c_new       + (1.-utheta)*p_coefs_c_old;

      visc = muu(tag);
      rho  = pho(Xqp,tag);
      Aloc.setZero();
      Gloc.setZero();
      Dloc.setZero();      
      FUloc.setZero();
      FPloc.setZero();


      if (behaviors & BH_bble_condens_PnPn)
      {
        iBbb.setZero();
        Bnb.setZero();
        Gbp.setZero();
        FUb.setZero();
        Bbn.setZero();
        Dpb.setZero();
      }

      if(behaviors & BH_GLS)
      {
        cell_volume = 0;
        for (int qp = 0; qp < n_qpts_cell; ++qp) {
          F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];
          J_mid = determinant(F_c_mid,dim);
          cell_volume += J_mid * quadr_cell->weight(qp);
        }

        hk2 = cell_volume / pi; // element size
        double const uconv = (u_coefs_c_old - v_coefs_c_mid).lpNorm<Infinity>();

        tauk = 4.*visc/hk2 + 2.*rho*uconv/sqrt(hk2);
        tauk = 1./tauk;
        if (dim==3)
          tauk *= 0.1;

        delk = 4.*visc + 2.*rho*uconv*sqrt(hk2);
        //delk = 0;

        Eloc.setZero();
        Cloc.setZero();
      }
      if (behaviors & BH_bble_condens_CR)
      {
        bble_integ = 0;
        Gnx.setZero();
        iBbb.setZero();
        Bnb.setZero();
        FUb.setZero();
        FPx.setZero();
        Bbn.setZero();

        cell_volume = 0;
        Xc.setZero();
        for (int qp = 0; qp < n_qpts_cell; ++qp) {
          F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];
          J_mid = determinant(F_c_mid,dim);
          Xqp  = x_coefs_c_mid_trans * qsi_c[qp];
          cell_volume += J_mid * quadr_cell->weight(qp);
          Xc += J_mid * quadr_cell->weight(qp) * Xqp;
        }
        Xc /= cell_volume;
      }


      // Quadrature
      for (int qp = 0; qp < n_qpts_cell; ++qp)
      {
        //F_c_new = x_coefs_c_new_trans * dLqsi_c[qp];
        //inverseAndDet(F_c_new,dim,invF_c_new,J_new);
        //invFT_c_new= invF_c_new.transpose();

        //F_c_old = x_coefs_c_old_trans * dLqsi_c[qp];
        //inverseAndDet(F_c_old,dim,invF_c_old,J_old);
        //invFT_c_old= invF_c_old.transpose();

        F_c_mid = x_coefs_c_mid_trans * dLqsi_c[qp];
        inverseAndDet(F_c_mid,dim,invF_c_mid,J_mid);
        invFT_c_mid= invF_c_mid.transpose();

        dxphi_c = dLphi_c[qp] * invF_c_mid;
        dxpsi_c = dLpsi_c[qp] * invF_c_mid;
        dxqsi_c = dLqsi_c[qp] * invF_c_mid;

        dxP_new  = dxpsi_c.transpose() * p_coefs_c_new;
        dxU  = u_coefs_c_mid_trans * dxphi_c;       // n+utheta

        Xqp      = x_coefs_c_mid_trans * qsi_c[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        Uqp      = u_coefs_c_mid_trans * phi_c[qp]; //n+utheta
        Uqp_new  = u_coefs_c_new_trans * phi_c[qp]; //n+utheta
        Uqp_old  = u_coefs_c_old_trans * phi_c[qp]; //n+utheta
        Pqp_new  = p_coefs_c_new.dot(psi_c[qp]);
        Vqp      = v_coefs_c_mid_trans * qsi_c[qp];
        Uconv_qp = Uqp - Vqp;
        //Uconv_qp = Uqp_old;
        dUdt     = (Uqp_new-Uqp_old)/dt;

        force_at_mid = force(Xqp,current_time+utheta*dt,tag);

        weight = quadr_cell->weight(qp);
        JxW_mid = J_mid*weight;

        //~ if (mesh->getCellId(&*cell) == 0)
        //~ {
          //~ printf("cHEcKKKKKKKKKKK!!\n");
          //~ cout << "x coefs mid:" << endl;
          //~ cout << x_coefs_c_mid_trans.transpose() << endl;
        //~ }
        if (J_mid < 1.e-10)
        {
          //FEP_PRAGMA_OMP(critical)
          //if (tid==0)
          {
            printf("in formCellFunction:\n");
            std::cout << "erro: jacobiana da integral não invertível: ";
            std::cout << "J_mid = " << J_mid << endl;
            cout << "trans matrix:\n" << F_c_mid << endl;
            cout << "x coefs mid:" << endl;
            cout << x_coefs_c_mid_trans.transpose() << endl;
            cout << "-----" << endl;
            cout << "cell id: " << mesh->getCellId(&*cell) << endl;
            cout << "cell Contig id: " << mesh->getCellContigId( mesh->getCellId(&*cell) ) << endl;
            cout << "cell nodes:\n" << cell_nodes.transpose() << endl;
            cout << "mapM :\n" << mapM_c.transpose() << endl;
            throw;
          }
        }

        for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
        {
          for (int c = 0; c < dim; ++c)
          {
            FUloc(i*dim + c) += JxW_mid*
                    ( rho*(unsteady*dUdt(c) + has_convec*Uconv_qp.dot(dxU.row(c)))*phi_c[qp][i] + // aceleração
                      visc*dxphi_c.row(i).dot(dxU.row(c) + dxU.col(c).transpose()) - //rigidez
                      Pqp_new*dxphi_c(i,c) - // pressão
                      force_at_mid(c)*phi_c[qp][i]   ); // força

            for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
            {
              for (int d = 0; d < dim; ++d)
              {
                delta_cd = c==d;
                Aloc(i*dim + c, j*dim + d) += JxW_mid*
                                              ( has_convec*phi_c[qp][i]*utheta *rho*( delta_cd*Uconv_qp.dot(dxphi_c.row(j))  +  dxU(c,d)*phi_c[qp][j] )   // advecção
                                              + unsteady* delta_cd*rho*phi_c[qp][i]*phi_c[qp][j]/dt     // time derivative
                                              + utheta*visc*( delta_cd * dxphi_c.row(i).dot(dxphi_c.row(j)) + dxphi_c(i,d)*dxphi_c(j,c))   ); // rigidez

              }
            }
            for (int j = 0; j < n_dofs_p_per_cell; ++j)
              Gloc(i*dim + c,j) -= JxW_mid * psi_c[qp][j]* dxphi_c(i,c);

          }

          //FUloc.segment(i*dim,dim) += JxW_mid*
          //                        (   phi_c[qp][i]*rho* (  dUdt + has_convec* dxU*Uconv_qp  ) // aceleração
          //                          + visc* ( dxU + dxU.transpose() )*dxphi_c.row(i).transpose()       //rigidez
          //                          - Pqp_new*dxphi_c.row(i).transpose() // pressão
          //                          - phi_c[qp][i]* force_at_mid   ); // força

        }
        for (int i = 0; i < n_dofs_p_per_cell; ++i)
          FPloc(i) -= JxW_mid* dxU.trace()*psi_c[qp][i];

        // ----------------
        //
        //  STABILIZATION
        //
        //  ----------------
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
                
                Bbn(c, j*dim + d) += JxW_mid*
                                     ( has_convec*bble[qp]*utheta *rho*( delta_cd*Uconv_qp.dot(dxphi_c.row(j))  +  dxU(c,d)*phi_c[qp][j] )   // convective
                                     + unsteady*delta_cd*rho*bble[qp]*phi_c[qp][j]/dt     // time derivative
                                     + utheta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) );    // rigidez                
                
                Bnb(j*dim + d, c) += JxW_mid*
                                     ( has_convec*phi_c[qp][j]*utheta *rho*(  delta_cd*Uconv_qp.dot(dxbble)  )   // convective
                                     + delta_cd*rho*phi_c[qp][j]*bble[qp]/dt  * unsteady    // time derivative
                                     + utheta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) );    // rigidez

              }
            }
            if (behaviors & BH_bble_condens_PnPn)
              for (int j = 0; j < n_dofs_p_per_cell; ++j)
                Gbp(c, j) -= JxW_mid*psi_c[qp][j]*dxbble(c);
          }

          for (int c = 0; c < dim; c++)
          {
            for (int d = 0; d < dim; d++)
            {
              delta_cd = c==d;
              iBbb(c, d) += JxW_mid*
                            ( has_convec*bble[qp]*utheta *rho*( delta_cd*Uconv_qp.dot(dxbble) )   // convective
                            + delta_cd*rho*bble[qp]*bble[qp]/dt  * unsteady   // time derivative
                            + utheta*visc*(delta_cd* dxbble.dot(dxbble)  + dxbble(d)*dxbble(c)) ); // rigidez
            }

            FUb(c) += JxW_mid*
                      ( bble[qp]*rho*(dUdt(c)*unsteady + has_convec*Uconv_qp.dot(dxU.row(c))) + // time derivative + convective
                        visc*dxbble.dot(dxU.row(c) + dxU.col(c).transpose()) - //rigidez
                        Pqp_new*dxbble(c) -                     // pressão
                        force_at_mid(c)*bble[qp]   ); // força
          }
        }
        else
        if(behaviors & BH_GLS)
        {
          Res = rho*( dUdt * unsteady +  has_convec*dxU*Uconv_qp) + dxP_new - force_at_mid;

          for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
          {
            dResdu = unsteady*(rho*phi_c[qp][j]/dt)*I + has_convec*rho*utheta*(  phi_c[qp][j]*dxU + Uconv_qp.dot(dxphi_c.row(j))*I );

            for (int i = 0; i < n_dofs_p_per_cell; ++i)
            {
              vec = dxpsi_c.row(i).transpose();
              vec = dResdu.transpose()*vec;
              vec = -JxW_mid*tauk*  vec;
              for (int d = 0; d < dim; d++)
                Dloc(i, j*dim + d) += vec(d);

              // atençao nos indices
              vec = JxW_mid*tauk*  has_convec*rho*Uconv_qp.dot(dxphi_c.row(j))* dxpsi_c.row(i).transpose();
              for (int d = 0; d < dim; d++)
                Cloc(j*dim + d,i) += vec(d);
            }

            for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
            {
              // supg term
              Ten = JxW_mid*tauk* has_convec*(  utheta*rho*phi_c[qp][j]*Res*dxphi_c.row(i) + rho*Uconv_qp.dot(dxphi_c.row(i))*dResdu  );
              // divergence term
              Ten+= JxW_mid*delk*utheta*dxphi_c.row(i).transpose()*dxphi_c.row(j);

              for (int c = 0; c < dim; ++c)
                for (int d = 0; d < dim; ++d)
                  Aloc(i*dim + c, j*dim + d) += Ten(c,d);
            }
          }

          for (int i = 0; i < n_dofs_p_per_cell; ++i)
            for (int j = 0; j < n_dofs_p_per_cell; ++j)
              Eloc(i,j) -= tauk*JxW_mid * dxphi_c.row(i).dot(dxphi_c.row(j));


          for (int i = 0; i < n_dofs_u_per_cell/dim; i++)
          {
            vec = JxW_mid*(  has_convec*tauk*  rho* Uconv_qp.dot(dxphi_c.row(i)) * Res    +   delk*dxU.trace()*dxphi_c.row(i).transpose() );

            for (int c = 0; c < dim; c++)
              FUloc(i*dim + c) += vec(c);

          }
          for (int i = 0; i < n_dofs_p_per_cell; ++i)
            FPloc(i) -= JxW_mid *tauk* dxpsi_c.row(i).dot(Res);
        }

        if (behaviors & BH_bble_condens_CR)
        {
          bble_integ += JxW_mid*bble[qp];

          for (int c = 0; c < dim; ++c)
          {
            for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
              for (int j = 0; j < dim; ++j) // pressure gradient
                Gnx(i*dim + c,j) -= JxW_mid* (Xqp(j) - Xc(j))*dxphi_c(i,c);

            FPx(c) -= JxW_mid* dxU.trace()*(Xqp(c) - Xc(c));
          }
        }


      } // fim quadratura

      Dloc += utheta*Gloc.transpose();

      //
      // stabilization
      //
      if (behaviors & BH_bble_condens_PnPn)
      {
        //iBbb = iBbb.inverse().eval();
        invert(iBbb,dim);

        FUloc = FUloc - Bnb*iBbb*FUb;
        FPloc = FPloc - utheta*Gbp.transpose()*iBbb*FUb;
        
        Dpb  = utheta*Gbp.transpose();

        // correções com os coeficientes da bolha

        Ubqp = -utheta*iBbb*FUb; // U bolha no tempo n+utheta        
        
        for (int qp = 0; qp < n_qpts_cell; ++qp)
        {
          F_c_mid    = x_coefs_c_mid_trans * dLqsi_c[qp];
          inverseAndDet(F_c_mid, dim, invF_c_mid,J_mid);
          invFT_c_mid= invF_c_mid.transpose();

          Uqp    = u_coefs_c_mid_trans * phi_c[qp];     //n+utheta
          dxbble = invFT_c_mid * dLbble[qp];
          dxUb   = Ubqp*dxbble.transpose();

          weight = quadr_cell->weight(qp);
          JxW_mid = J_mid*weight;

          for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
          {
            for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
            {
              Ten = has_convec*JxW_mid*rho*utheta*phi_c[qp][i]* phi_c[qp][j] * dxUb;  // advecção

              for (int c = 0; c < dim; ++c)
                for (int d = 0; d < dim; ++d)
                  Aloc(i*dim + c, j*dim + d) += Ten(c,d);
            }
            Ten = has_convec*JxW_mid*  rho*utheta* bble[qp] * phi_c[qp][j] *dxUb; // advecção
            for (int c = 0; c < dim; ++c)
              for (int d = 0; d < dim; ++d)
                Bbn(c, j*dim + d) += Ten(c,d);
          }
        } // fim quadratura 2 vez

        Aloc -= Bnb*iBbb*Bbn;
        Gloc -= Bnb*iBbb*Gbp;
        Dloc -= Dpb*iBbb*Bbn;
        Eloc = -Dpb*iBbb*Gbp;        
        
      }
      if(behaviors & BH_GLS)
      {
        Gloc += Cloc;
      }
      if(behaviors & BH_bble_condens_CR)
      {
        Ubqp.setZero();
        //for (int i = 0; i < Gnx.cols(); ++i)
        //  for (int j = 0; j < Gnx.rows(); ++j)
        //    Ubqp(i) += Gnx(j,i)*u_coefs_c_new(j);
        //Ubqp /= -bble_integ;
        //Ubqp *= utheta;

        for (int c = 0; c < dim; ++c)
          for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
            for (int j = 0; j < dim; ++j) // pressure gradient
              //Gnx(i*dim + c,j) -= JxW_mid* (Xqp(j) - Xc(j))*dxphi_c(i,c);
            Ubqp(j) += Gnx(i*dim + c,j) * u_coefs_c_mid_trans(c,i);

        Ubqp /= -bble_integ;
        Ubqp *= utheta;


        //Ubqp = -Gnx.transpose()*u_coefs_c_new; // U bolha no tempo n+utheta

        for (int qp = 0; qp < n_qpts_cell; ++qp)
        {
          F_c_mid    = x_coefs_c_mid_trans * dLqsi_c[qp];
          inverseAndDet(F_c_mid, dim, invF_c_mid,J_mid);
          invFT_c_mid= invF_c_mid.transpose();

          Uqp    = u_coefs_c_mid_trans * phi_c[qp];     //n+utheta
          dxbble = invFT_c_mid * dLbble[qp];
          dxUb   = Ubqp*dxbble.transpose();

          weight = quadr_cell->weight(qp);
          JxW_mid = J_mid*weight;

          for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
          {
            for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
            {
              Ten = has_convec*JxW_mid*rho*utheta*phi_c[qp][i]* phi_c[qp][j] * dxUb;  // advecção

              for (int c = 0; c < dim; ++c)
                for (int d = 0; d < dim; ++d)
                  Aloc(i*dim + c, j*dim + d) += Ten(c,d);
            }
            Ten = has_convec*JxW_mid*  rho*utheta* bble[qp] * phi_c[qp][j] *dxUb; // advecção
            for (int c = 0; c < dim; ++c)
              for (int d = 0; d < dim; ++d)
                Bbn(c, j*dim + d) += Ten(c,d);
          }
        } // fim quadratura 2 vez

        double const a = 1./(bble_integ*bble_integ);
        double const b = 1./bble_integ;
        Aloc += utheta*a*Gnx*iBbb*Gnx.transpose() - utheta*b*Bnb*Gnx.transpose() - b*Gnx*Bbn;        
        
        FUloc +=  a*Gnx*iBbb*FPx - b*Bnb*FPx - b*Gnx*FUb;
      }



      // Projection - to force non-penetrarion bc
      mesh->getCellNodesId(&*cell, cell_nodes.data());
      getProjectorMatrix(Prj, nodes_per_cell, cell_nodes.data(), Vec_x_1, current_time+dt, *this);

      FUloc = Prj*FUloc;
      Aloc = Prj*Aloc*Prj;
      Gloc = Prj*Gloc;      
      Dloc = Dloc*Prj;
      
      if (force_pressure)
      {
        for (int i = 0; i < mapP_c.size(); ++i)
        {
          if (mapP_c(i) == null_space_press_dof)
          {
            Gloc.col(i).setZero();
            Dloc.row(i).setZero();
            FPloc(i) = 0;
            break;
          }
        }
      }
      
#ifdef FEP_HAS_OPENMP
      FEP_PRAGMA_OMP(critical)
#endif
      {
        VecSetValues(Vec_fun, mapU_c.size(), mapU_c.data(), FUloc.data(), ADD_VALUES);
        VecSetValues(Vec_fun, mapP_c.size(), mapP_c.data(), FPloc.data(), ADD_VALUES);
        MatSetValues(*JJ, mapU_c.size(), mapU_c.data(), mapU_c.size(), mapU_c.data(), Aloc.data(),  ADD_VALUES);
        MatSetValues(*JJ, mapU_c.size(), mapU_c.data(), mapP_c.size(), mapP_c.data(), Gloc.data(),  ADD_VALUES);
        MatSetValues(*JJ, mapP_c.size(), mapP_c.data(), mapU_c.size(), mapU_c.data(), Dloc.data(),  ADD_VALUES);
        if (pres_pres_block)
        {
          MatSetValues(*JJ, mapP_c.size(), mapP_c.data(), mapP_c.size(), mapP_c.data(), Eloc.data(),  ADD_VALUES);
        }        
      }
    }


  }

  // LOOP NAS FACES DO CONTORNO
  //~ FEP_PRAGMA_OMP(parallel default(none) shared(Vec_up_k,Vec_fun,cout))
  {
    int                 tag;
    bool                is_neumann;
    bool                is_surface;
    bool                is_solid;
    //MatrixXd           u_coefs_c_new(n_dofs_u_per_facet/dim, dim);
    //VectorXd           p_coefs_f(n_dofs_p_per_facet);
    MatrixXd            u_coefs_f_mid_trans(dim, n_dofs_u_per_facet/dim);  // n+utheta
    MatrixXd            u_coefs_f_old(n_dofs_u_per_facet/dim, dim);        // n
    MatrixXd            u_coefs_f_new(n_dofs_u_per_facet/dim, dim);        // n+1
    MatrixXd            u_coefs_f_old_trans(dim,n_dofs_u_per_facet/dim);   // n
    MatrixXd            u_coefs_f_new_trans(dim,n_dofs_u_per_facet/dim);   // n+1

    MatrixXd            x_coefs_f_mid_trans(dim, n_dofs_v_per_facet/dim); // n+utheta
    MatrixXd            x_coefs_f_new(n_dofs_v_per_facet/dim, dim);       // n+1
    MatrixXd            x_coefs_f_new_trans(dim, n_dofs_v_per_facet/dim); // n+1
    MatrixXd            x_coefs_f_old(n_dofs_v_per_facet/dim, dim);       // n
    MatrixXd            x_coefs_f_old_trans(dim, n_dofs_v_per_facet/dim); // n

    MatrixXd            noi_coefs_f_new(n_dofs_v_per_facet/dim, dim);  // normal interpolada em n+1
    MatrixXd            noi_coefs_f_new_trans(dim, n_dofs_v_per_facet/dim);  // normal interpolada em n+1

    Tensor              F_f_mid(dim,dim-1);       // n+utheta
    Tensor              invF_f_mid(dim-1,dim);    // n+utheta
    Tensor              fff_f_mid(dim-1,dim-1);   // n+utheta; fff = first fundamental form
    //Tensor              invFT_c_mid(dim,dim);   // n+utheta

    MatrixXd            Aloc_f(n_dofs_u_per_facet, n_dofs_u_per_facet);
    VectorXd            FUloc(n_dofs_u_per_facet);
                        
    MatrixXd            tmp(n_dofs_u_per_facet,n_dofs_u_per_facet);
                        
    VectorXi            mapU_f(n_dofs_u_per_facet);
    VectorXi            mapP_f(n_dofs_p_per_facet);
    VectorXi            mapM_f(dim*nodes_per_facet);
                        
    MatrixXd            Prj(n_dofs_u_per_facet,n_dofs_u_per_facet);

    MatrixXd            dxphi_f(n_dofs_u_per_facet/dim, dim);
    Tensor              dxU_f(dim,dim);   // grad u
    Vector              Xqp(dim);
    Vector              Xqp2(dim);
    Vector              Xqp_new(dim);
    Vector              Xqp_old(dim);
    Vector              Uqp(dim);
    Vector              Uqp_new(dim);
    Vector              Uqp_old(dim);
    //VectorXd          FUloc(n_dofs_u_per_facet);
    VectorXi            facet_nodes(nodes_per_facet);
    Vector              normal(dim);
    Vector              noi(dim); // normal interpolada
    Vector              some_vec(dim);
    double              J_mid=0,JxW_mid;
    double              weight=0;
    //double              visc;
    //double              rho;
    Vector              Uqp_solid(dim);

    Vector              traction_(dim);

    //~ const int tid = omp_get_thread_num();
    //~ const int nthreads = omp_get_num_threads();
//~
    //~ facet_iterator facet = mesh->facetBegin(tid,nthreads);
    //~ facet_iterator facet_end = mesh->facetEnd(tid,nthreads);

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
      dof_handler[DH_MESH].getVariable(VAR_M).getFacetDofs(mapM_f.data(), &*facet);

      VecGetValues(Vec_normal,  mapM_f.size(), mapM_f.data(), noi_coefs_f_new.data());
      VecGetValues(Vec_x_0,     mapM_f.size(), mapM_f.data(), x_coefs_f_old.data());
      VecGetValues(Vec_x_1,     mapM_f.size(), mapM_f.data(), x_coefs_f_new.data());
      VecGetValues(Vec_up_0,    mapU_f.size(), mapU_f.data(), u_coefs_f_old.data());
      VecGetValues(Vec_up_k ,   mapU_f.size(), mapU_f.data(), u_coefs_f_new.data());
      
      // get nodal coordinates of the old and new cell
      mesh->getFacetNodesId(&*facet, facet_nodes.data());

      x_coefs_f_old_trans = x_coefs_f_old.transpose();
      x_coefs_f_new_trans = x_coefs_f_new.transpose();
      u_coefs_f_old_trans = u_coefs_f_old.transpose();
      u_coefs_f_new_trans = u_coefs_f_new.transpose();
      noi_coefs_f_new_trans = noi_coefs_f_new.transpose();

      u_coefs_f_mid_trans = utheta*u_coefs_f_new_trans + (1.-utheta)*u_coefs_f_old_trans;
      x_coefs_f_mid_trans = utheta*x_coefs_f_new_trans + (1.-utheta)*x_coefs_f_old_trans;

      FUloc.setZero();
      Aloc_f.setZero();

      //visc = muu(tag);
      //rho  = pho(Xqp,tag);

      //noi_coefs_f_new_trans = x_coefs_f_mid_trans;

      for (int qp = 0; qp < n_qpts_facet; ++qp)
      {

        F_f_mid   = x_coefs_f_mid_trans * dLqsi_f[qp];

        if (dim==2)
        {
          normal(0) = +F_f_mid(1,0);
          normal(1) = -F_f_mid(0,0);
          normal.normalize();
        }
        else
        {
          normal = cross(F_f_mid.col(0), F_f_mid.col(1));
          normal.normalize();
        }

        fff_f_mid.resize(dim-1,dim-1);
        fff_f_mid  = F_f_mid.transpose()*F_f_mid;
        J_mid     = sqrt(fff_f_mid.determinant());
        invF_f_mid = fff_f_mid.inverse()*F_f_mid.transpose();

        weight  = quadr_facet->weight(qp);
        JxW_mid = J_mid*weight;
        Xqp     = x_coefs_f_mid_trans * qsi_f[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        dxphi_f = dLphi_f[qp] * invF_f_mid;
        dxU_f   = u_coefs_f_mid_trans * dxphi_f; // n+utheta
        Uqp     = u_coefs_f_mid_trans * phi_f[qp];
        noi     = noi_coefs_f_new_trans * qsi_f[qp];

        if (is_neumann)
        {
          //Vector no(Xqp);
          //no.normalize();
          //traction_ = utheta*(traction(Xqp,current_time+dt,tag)) + (1.-utheta)*traction(Xqp,current_time,tag);
          traction_ = traction(Xqp, normal, current_time + dt*utheta,tag);
          //traction_ = (traction(Xqp,current_time,tag) +4.*traction(Xqp,current_time+dt/2.,tag) + traction(Xqp,current_time+dt,tag))/6.;

          for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
          {
            for (int c = 0; c < dim; ++c)
            {
              FUloc(i*dim + c) -= JxW_mid * traction_(c) * phi_f[qp][i] ; // força
            }
          }
        }

        if (is_surface)
        {
          //Vector no(Xqp);
          //no.normalize();
          for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
          {
            for (int c = 0; c < dim; ++c)
            {
              //FUloc(i*dim + c) += JxW_mid *gama(Xqp,current_time,tag)*(dxphi_f(i,c) + 0*(unsteady*dt) *dxU_f.row(c).dot(dxphi_f.row(i))); // correto
              FUloc(i*dim + c) += JxW_mid *gama(Xqp,current_time,tag)*dxphi_f(i,c);
              //for (int d = 0; d < dim; ++d)
              //  FUloc(i*dim + c) += JxW_mid * gama(Xqp,current_time,tag)* ( (c==d?1:0) - noi(c)*noi(d) )* dxphi_f(i,d) ;
              //FUloc(i*dim + c) += JxW_mid * gama(Xqp,current_time,tag)* ( unsteady*dt *dxU_f.row(c).dot(dxphi_f.row(i)));
            }
          }
          
          if (false) // semi-implicit term
          {
            for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
              for (int j = 0; j < n_dofs_u_per_facet/dim; ++j)
                for (int c = 0; c < dim; ++c)
                  Aloc_f(i*dim + c, j*dim + c) += utheta*JxW_mid* (unsteady*dt) *gama(Xqp,current_time,tag)*dxphi_f.row(i).dot(dxphi_f.row(j));
          }          
          
        }

        if (is_solid)
        {
          Uqp_solid = solid_veloc(Xqp, current_time+utheta*dt, tag);

          for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
          {
            for (int c = 0; c < dim; ++c)
            {
              FUloc(i*dim + c) += JxW_mid *beta_diss()*(Uqp(c)-Uqp_solid(c))*phi_f[qp][i];
              //FUloc(i*dim + c) += x_coefs_f_old_trans.norm()*beta_diss()*Uqp(c)*phi_f[qp][i];
              
              for (int j = 0; j < n_dofs_u_per_facet/dim; ++j)
                  Aloc_f(i*dim + c, j*dim + c) += utheta*JxW_mid *beta_diss()*phi_f[qp][j]*phi_f[qp][i];                   
              
            }
          }
        }

      } // end quadratura
     

      // Projection - to force non-penetrarion bc
      mesh->getFacetNodesId(&*facet, facet_nodes.data());
      getProjectorMatrix(Prj, nodes_per_facet, facet_nodes.data(), Vec_x_1, current_time+dt, *this);

      FUloc = Prj*FUloc;
      Aloc_f = Prj*Aloc_f*Prj;

      //~ FEP_PRAGMA_OMP(critical)
      {
        VecSetValues(Vec_fun, mapU_f.size(), mapU_f.data(), FUloc.data(), ADD_VALUES);
        MatSetValues(*JJ, mapU_f.size(), mapU_f.data(), mapU_f.size(), mapU_f.data(), Aloc_f.data(),  ADD_VALUES);
      }

    }


  } // end parallel


  // LINHA DE CONTATO
  //FEP_PRAGMA_OMP(parallel shared(Vec_up_k,Vec_fun,cout) default(none))
  {
    // will be useful I hope
    //Real const eps = std::numeric_limits<Real>::epsilon();
    //Real const eps_root = pow(eps,1./3.);    
    //double            h;
    //volatile double   hh;
    
    int              tag;
    bool             is_triple;

    VectorXd         FUloc(n_dofs_u_per_corner);
    MatrixXd         Aloc_r(n_dofs_u_per_corner, n_dofs_u_per_corner);

    VectorXi         mapU_r(n_dofs_u_per_corner);
    VectorXi         mapP_r(n_dofs_p_per_corner);

    MatrixXd         Prj(n_dofs_u_per_corner,n_dofs_u_per_corner);
    VectorXi         corner_nodes(nodes_per_corner);

    bool                gen_error = false;
    //MatrixXd             u_coefs_r_mid(n_dofs_u_per_corner/dim, dim);
    MatrixXd            u_coefs_r_mid_trans(dim, n_dofs_u_per_corner/dim);  // n+utheta
    MatrixXd            u_coefs_r_old(n_dofs_u_per_corner/dim, dim);        // n
    MatrixXd            u_coefs_r_new(n_dofs_u_per_corner/dim, dim);        // n+1
    MatrixXd            x_coefs_r_mid_trans(dim, nodes_per_corner);
    MatrixXd            x_coefs_r_new(nodes_per_corner, dim);
    MatrixXd            x_coefs_r_old(nodes_per_corner, dim);
    Tensor              F_r_mid(dim,dim-2);
    Tensor              invF_r_mid(dim-2,dim);
    MatrixXd            dxphi_r(n_dofs_u_per_corner/dim, dim);
    Tensor              dxU_r(dim,dim);   // grad u
    Vector              Xqp(dim);
    Vector              Uqp(dim);
    //VectorXd          FUloc(n_dofs_u_per_corner);
    //MatrixXd         Aloc_r(n_dofs_u_per_corner, n_dofs_u_per_corner);
    Vector              normal(dim);
    Vector              line_normal(dim);
    Vector              solid_point(dim); // ponto na superfície do sólido .. ele é único
    Vector              point_a(dim); // ponto na linha triplice
    Vector              point_b(dim); // ponto na linha triplice
    Vector              ifacet_normal(dim); // ponto na linha triplice
    double              line_normal_sign; // +1 or -1
    double              J_mid=0, JxW_mid;
    double              weight=0;
    double              gama_mid;
    //double              visc;
    //double              rho;
    int                 iCs[FEPIC_MAX_ICELLS];
    int                 eiCs[FEPIC_MAX_ICELLS];
    int                 *iCs_end;
    int                 *iCs_it;
    Cell                *fluid_cell;

    //VectorXi            mapU_r(n_dofs_u_per_corner);
    //VectorXi            mapP_r(n_dofs_p_per_corner);
    VectorXi            mapM_r(dim*nodes_per_corner);

    //const int tid = omp_get_thread_num();
    //const int nthreads = omp_get_num_threads();
    //const int n_corner_colors = mesh->numCornerColors();


    // LOOP NAS ARESTAS DA LINHA TRIPLICE
    CellElement * corner;

    if (triple_tags.size() != 0)
    for (int _r = 0; _r < n_corners_total; ++_r)
    {
      if (dim==2)
      {
        corner = mesh->getNodePtr(_r);
        if (!mesh->isVertex(corner))
          continue;
      }
      else
        corner = mesh->getCornerPtr(_r);
      if (corner->isDisabled())
        continue;

      tag = corner->getTag();
      is_triple = is_in(tag,triple_tags);
      if (!is_triple)
        continue;

      FUloc.setZero();
      Aloc_r.setZero();

      mesh->getCornerNodesId(&*corner, corner_nodes.data());
     // mesh->getNodesCoords(corner_nodes.begin(), corner_nodes.end(), x_coefs_r_mid.data());

      dof_handler[DH_UNKS].getVariable(VAR_U).getCornerDofs(mapU_r.data(), &*corner);
      dof_handler[DH_UNKS].getVariable(VAR_P).getCornerDofs(mapP_r.data(), &*corner);
      dof_handler[DH_MESH].getVariable(VAR_M).getCornerDofs(mapM_r.data(), &*corner);

      VecGetValues(Vec_x_0,     mapM_r.size(), mapM_r.data(), x_coefs_r_old.data());
      VecGetValues(Vec_x_1,     mapM_r.size(), mapM_r.data(), x_coefs_r_new.data());
      VecGetValues(Vec_up_0,    mapU_r.size(), mapU_r.data(), u_coefs_r_old.data());
      VecGetValues(Vec_up_k,    mapU_r.size(), mapU_r.data(), u_coefs_r_new.data());

      u_coefs_r_mid_trans = utheta*u_coefs_r_new.transpose() + (1.-utheta)*u_coefs_r_old.transpose();
      x_coefs_r_mid_trans = utheta*x_coefs_r_new.transpose() + (1.-utheta)*x_coefs_r_old.transpose();

      //visc = muu(tag);
      //rho  = pho(Xqp,tag);

      if (dim==3)
      {
        // encontrando o ponto da superfície sólido. com ele, é calculado uma normal
        // que corrigi o sinal de line_normal
        iCs_end = mesh->edgeStar(corner, iCs, eiCs);
        if (iCs_end == iCs)
        {
          printf("ERROR!: no icell found\n");
          throw;
        }
        gen_error = true;
        for (iCs_it = iCs; iCs_it != iCs_end ; ++iCs_it)
        {
          fluid_cell = mesh->getCellPtr(*iCs_it);

          for (int kk = 0; kk < mesh->numVerticesPerCell(); ++kk)
          {
            Point const* pp      = mesh->getNodePtr(fluid_cell->getNodeId(kk) );
            const int    tag_aux = pp->getTag();

            if ((is_in(tag_aux, solid_tags) || is_in(tag_aux,feature_tags)) && !is_in(tag_aux, triple_tags))
            {
              gen_error = false;
              pp->getCoord(solid_point.data(), dim);
              break;
            }
          }
          if (gen_error==false)
            break;

        }
        if (gen_error)
        {
          printf("ERROR!: solid point not found\n");
          cout << "corner id: " << (mesh->getCellPtr(corner->getIncidCell())->getCornerId(corner->getPosition())) << endl;
          cout << "first icell : " << (*iCs) << endl;
          mesh->getNodePtr(corner_nodes(0))->getCoord(point_a.data(),dim);
          mesh->getNodePtr(corner_nodes(1))->getCoord(point_b.data(),dim);
          cout << "point a = " << point_a[0] << " " << point_a[1] << " " << point_a[2] << "\n";
          cout << "point b = " << point_b[0] << " " << point_b[1] << " " << point_b[2] << "\n";
          cout << "number of cells = " << (iCs_end-iCs) << endl;
          throw;
        }


        mesh->getNodePtr(corner_nodes(0))->getCoord(point_a.data(),dim);
        mesh->getNodePtr(corner_nodes(1))->getCoord(point_b.data(),dim);

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
        Point * point = mesh->getNodePtr(corner_nodes[0]);
        Point * sol_point = NULL;
        Point * sol_point_2;
        int iVs[FEPIC_MAX_ICELLS];
        int *iVs_end, *iVs_it;
        Vector aux(dim);

        iVs_end = mesh->connectedVtcs(point, iVs);

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



      for (int qp = 0; qp < n_qpts_corner; ++qp)
      {
        if (dim==3)
        {
          F_r_mid   = x_coefs_r_mid_trans * dLqsi_r[qp];
          J_mid = F_r_mid.norm();
          weight  = quadr_corner->weight(qp);
          JxW_mid = J_mid*weight;
        }
        else
        {
          J_mid = 1;
          weight  = 1;
          JxW_mid = 1;
        }
        //invF_r_mid = F_r_mid.transpose()/(J_mid*J_mid);


        Xqp = x_coefs_r_mid_trans * qsi_r[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        //dxphi_r = dLphi_r[qp] * invF_r_mid;
        //dxU_r   = u_coefs_r_mid_trans * dxphi_r; // n+utheta
        Uqp  = u_coefs_r_mid_trans * phi_r[qp];

        gama_mid = gama(Xqp, current_time+dt*utheta, tag);

        if (dim==3)
        {
          normal  = solid_normal(Xqp, current_time, tag);
          line_normal(0)= F_r_mid(0,0);
          line_normal(1)= F_r_mid(1,0);
          line_normal(2)= F_r_mid(2,0);
          line_normal *= line_normal_sign;
          line_normal = cross(line_normal, normal);
          line_normal.normalize();
        }

        double const Uqp_norm = abs(line_normal.dot(Uqp));
        for (int i = 0; i < n_dofs_u_per_corner/dim; ++i)
        {
          for (int c = 0; c < dim; ++c)
          {
            FUloc(i*dim + c) += JxW_mid*(-gama_mid*cos_theta0() + zeta(Uqp_norm,0)*line_normal.dot(Uqp))*line_normal(c)*phi_r[qp][i];
            //FUloc(i*dim + c) += JxW_mid*(-gama_mid*cos_theta0() )*line_normal(c)*phi_r[qp][i];
            //FUloc(i*dim + c) += JxW_mid*zeta(Uqp_norm,0)*Uqp(c)*phi_r[qp][i];
          
            for (int j = 0; j < n_dofs_u_per_corner/dim; ++j)
              for (int d = 0; d < dim; ++d)
                Aloc_r(i*dim + c, j*dim + d) += JxW_mid* zeta(Uqp_norm,0)*utheta*line_normal(d) *phi_r[qp][j]*line_normal(c)*phi_r[qp][i];
                //if (c==d)
                //  Aloc_r(i*dim + c, j*dim + d) += JxW_mid* zeta(Uqp_norm,0) *utheta*phi_r[qp][j]*phi_r[qp][i];
          
          }
          
          //FUloc.segment(i*dim, dim) += JxW_mid* phi_r[qp][i] * (-gama_mid*cos_theta0() + zeta(0,0)*line_normal.dot(Uqp))*line_normal;
        }


      } // end quadratura


     

      // Projection - to force non-penetrarion bc
      mesh->getCornerNodesId(&*corner, corner_nodes.data());
      getProjectorMatrix(Prj, nodes_per_corner, corner_nodes.data(), Vec_x_1, current_time+dt, *this);

      FUloc = Prj*FUloc;
      Aloc_r = Prj*Aloc_r*Prj;
      
      VecSetValues(Vec_fun, mapU_r.size(), mapU_r.data(), FUloc.data(), ADD_VALUES);
      MatSetValues(*JJ, mapU_r.size(), mapU_r.data(), mapU_r.size(), mapU_r.data(), Aloc_r.data(),  ADD_VALUES);
      //cout << FUloc.transpose() << endl;

    }



  }

  // boundary conditions on global Jacobian
    // solid & triple tags .. force normal
  if (force_dirichlet)
  {
    int      nodeid;
    int      u_dofs[dim];
    Vector   normal(dim);
    Tensor   A(dim,dim);
    Tensor   I(Tensor::Identity(dim,dim));
    int      tag;

    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    for ( ; point != point_end; ++point)
    {
      tag = point->getTag();
      if (!(is_in(tag,feature_tags)   ||
            is_in(tag,solid_tags)     ||
            is_in(tag,interface_tags) ||
            is_in(tag,triple_tags)    ||
            is_in(tag,dirichlet_tags) ||
            is_in(tag,neumann_tags)   ||
            is_in(tag,periodic_tags)     )  )
        continue;
      //dof_handler[DH_UNKS].getVariable(VAR_U).getVertexAssociatedDofs(u_dofs, &*point);
      getNodeDofs(&*point, DH_UNKS, VAR_U, u_dofs);
      
      nodeid = mesh->getPointId(&*point);
      getProjectorMatrix(A, 1, &nodeid, Vec_x_1, current_time+dt, *this);
      A = I - A;
      MatSetValues(*JJ, dim, u_dofs, dim, u_dofs, A.data(), ADD_VALUES);
    }
  }

  if (force_pressure)
  {
    double const p =1.0;
    MatSetValues(*JJ, 1, &null_space_press_dof, 1, &null_space_press_dof, &p, ADD_VALUES);
  }

  if(print_to_matlab)
  {
    static bool ja_foi=false;
    if (!ja_foi)
    { 
      View(Vec_fun, "rhs.m","res");
      View(*JJ,"jacob.m","Jac");
    }
    ja_foi = true;
    
  }

  Assembly(Vec_fun);
  Assembly(*JJ);

  PetscFunctionReturn(0);

} // END formFunction


PetscErrorCode AppCtx::formJacobian(SNES /*snes*/,Vec /*Vec_up_k*/,Mat* /*Mat_Jac*/, Mat* /*prejac*/, MatStructure * /*flag*/)
{
  PetscFunctionReturn(0);
}

// ====================================================================================================

// to apply boundary conditions on linear elasticity problem.
template <typename Derived>
void getProjectorBC(MatrixBase<Derived> & P, int n_nodes, int const* nodes, Vec const& Vec_x_, double t, AppCtx const& app)
{
  int const dim = app.dim;
  Mesh const* mesh = &*app.mesh;
  //DofHandler const* dof_handler = &*app.dof_handler;
  std::vector<int> const& dirichlet_tags  = app.dirichlet_tags;
  std::vector<int> const& neumann_tags    = app.neumann_tags  ;
  std::vector<int> const& interface_tags  = app.interface_tags;
  std::vector<int> const& solid_tags      = app.solid_tags    ;
  std::vector<int> const& triple_tags     = app.triple_tags   ;
  std::vector<int> const& periodic_tags   = app.periodic_tags ;
  std::vector<int> const& feature_tags    = app.feature_tags  ;
  Vec const& Vec_normal = app.Vec_normal;
  
  P.setIdentity();

  Tensor I(dim,dim);
  Tensor Z(dim,dim);
  Vector X(dim);
  Vector normal(dim);
  int    dofs[dim];
  int    tag;
  Point const* point;

  I.setIdentity();
  Z.setZero();

  // NODES
  for (int i = 0; i < n_nodes; ++i)
  {
    point = mesh->getNodePtr(nodes[i]);
    tag = point->getTag();
    //m = point->getPosition() - mesh->numVerticesPerCell();
    //cell = mesh->getCellPtr(point->getIncidCell());

    if (is_in(tag,feature_tags))
    {
      app.getNodeDofs(&*point, DH_MESH, VAR_M, dofs);
      VecGetValues(Vec_x_, dim, dofs, X.data());
      P.block(i*dim,i*dim,dim,dim)  = feature_proj(X,t,tag);
    }
    else
    if (is_in(tag,solid_tags) )
    {
      if (app.boundary_smoothing)
      {
        app.getNodeDofs(&*point, DH_MESH, VAR_M, dofs);
        VecGetValues(Vec_x_, dim, dofs, X.data());
        normal = -solid_normal(X,t,tag);
        P.block(i*dim,i*dim,dim,dim)  = I - normal*normal.transpose();      
      }
      else
      {
        P.block(i*dim,i*dim,dim,dim) = Z;  
      }
    }
    else
    if (is_in(tag,interface_tags))
    {
      if (app.boundary_smoothing)
      {
        app.getNodeDofs(&*point, DH_MESH, VAR_M, dofs);
        VecGetValues(Vec_normal, dim, dofs, X.data());
        P.block(i*dim,i*dim,dim,dim) = I - X*X.transpose();
      }
      else
      {
        P.block(i*dim,i*dim,dim,dim) = Z;
      }
    }
    else
    if (is_in(tag,triple_tags) || is_in(tag,dirichlet_tags) || is_in(tag,neumann_tags) || is_in(tag,periodic_tags))
    {
      P.block(i*dim,i*dim,dim,dim) = Z;
    }


  } // end nodes
}



// function to compute mesh velocity
PetscErrorCode AppCtx::formFunction_mesh(SNES /*snes_m*/, Vec Vec_v, Vec Vec_fun)
{
  Mat *JJ = &Mat_Jac_m;

  //SNESGetJacobian(snes_m, JJ, NULL, NULL, NULL);

  // NOTE: solve elasticity problem in the mesh at time step n
  // NOTE: The mesh used is the Vec_x_0
  // WARNING: this function assumes that the boundary conditions was already applied

  // LOOP NAS CÉLULAS
  VecZeroEntries(Vec_fun);
  MatZeroEntries(*JJ);


#ifdef FEP_HAS_OPENMP
  FEP_PRAGMA_OMP(parallel default(none) shared(Vec_v, Vec_fun, cout))
#endif
  {
    bool const non_linear = nonlinear_elasticity;
    
    Tensor            dxV(dim,dim);   // grad u
    Tensor            F_c(dim,dim);
    Tensor            invF_c(dim,dim);
    Tensor            invFT_c(dim,dim);
    Vector            Vqp(dim);
    MatrixXd          v_coefs_c_trans(dim, nodes_per_cell);      // mesh velocity;
    MatrixXd          v_coefs_c(nodes_per_cell, dim);
    MatrixXd          x_coefs_c_trans(dim, nodes_per_cell);
    MatrixXd          x_coefs_c(nodes_per_cell, dim);
    MatrixXd          dxqsi_c(nodes_per_cell, dim);
    double            J, weight, JxW;

    VectorXd          Floc(n_dofs_v_per_cell);
    MatrixXd          Aloc(n_dofs_v_per_cell, n_dofs_v_per_cell);

    VectorXi          mapV_c(n_dofs_u_per_cell);

    MatrixXd          Prj(n_dofs_v_per_cell, n_dofs_v_per_cell);
    VectorXi          cell_nodes(nodes_per_cell);

    double            sigma_ck;
    double            dsigma_ckjd;  // dphi, d_compd

    const int tid = omp_get_thread_num();
    const int nthreads = omp_get_num_threads();

    cell_iterator cell = mesh->cellBegin(tid,nthreads);
    cell_iterator cell_end = mesh->cellEnd(tid,nthreads);

    //cell_iterator cell = mesh->cellBegin();
    //cell_iterator cell_end = mesh->cellEnd();
    for (; cell != cell_end; ++cell)
    {

      // mapeamento do local para o global:
      //
      dof_handler[DH_MESH].getVariable(VAR_M).getCellDofs(mapV_c.data(), &*cell);

      /*  Pega os valores das variáveis nos graus de liberdade */
      VecGetValues(Vec_v ,  mapV_c.size(), mapV_c.data(), v_coefs_c.data());
      VecGetValues(Vec_x_0, mapV_c.size(), mapV_c.data(), x_coefs_c.data());
      
      v_coefs_c_trans = v_coefs_c.transpose();
      x_coefs_c_trans = x_coefs_c.transpose();
      
      Floc.setZero();
      Aloc.setZero();

      // Quadrature
      for (int qp = 0; qp < n_qpts_cell; ++qp)
      {
        F_c = x_coefs_c_trans * dLqsi_c[qp];
        inverseAndDet(F_c,dim,invF_c,J);
        invFT_c= invF_c.transpose();

        dxqsi_c = dLqsi_c[qp] * invF_c;

        dxV  = v_coefs_c_trans * dxqsi_c;       // n+utheta
        Vqp  = v_coefs_c_trans * qsi_c[qp];
        //Xqp      = x_coefs_c_trans * qsi_c[qp]; // coordenada espacial (x,y,z) do ponto de quadratura

        weight = quadr_cell->weight(qp);
        JxW = J*weight;


        for (int i = 0; i < n_dofs_v_per_cell/dim; ++i)
        {
          for (int c = 0; c < dim; ++c)
          {
            for (int k = 0; k < dim; ++k)
            {
              sigma_ck = dxV(c,k) + dxV(k,c);

              if (non_linear)
              {
                for (int l = 0; l < dim; ++l)
                {
                  sigma_ck += dxV(l,c)*dxV(l,k);

                  if (c==k)
                  {
                    sigma_ck -= dxV(l,l);
                    for (int m = 0; m < dim; ++m)
                      sigma_ck -=  dxV(l,m)*dxV(l,m); 
                  }

                }
              }

              Floc(i*dim + c) += sigma_ck*dxqsi_c(i,k); // (JxW/JxW) is to compiler not complain about unused variables

              for (int j = 0; j < n_dofs_v_per_cell/dim; ++j)
              {
                for (int d = 0; d < dim; ++d)
                {
                  dsigma_ckjd = 0;
                  
                  if (c==d)
                    dsigma_ckjd = dxqsi_c(j,k);

                  if (k==d)
                    dsigma_ckjd += dxqsi_c(j,c);

                  if (non_linear)
                  {
                    for (int l = 0; l < dim; ++l)
                    {
                      if (l==d)
                        dsigma_ckjd += dxqsi_c(j,c)*dxV(l,k) + dxV(l,c)*dxqsi_c(j,k);

                      if (c==k)
                      {
                        if (l==d)
                        {
                          dsigma_ckjd -= dxqsi_c(j,l);
                          for (int m = 0; m < dim; ++m)
                            dsigma_ckjd -= 2.*dxqsi_c(j,m)*dxV(l,m);
                        }
                      }
                    }
                  }
                  
                  Aloc(i*dim + c, j*dim + d) += dsigma_ckjd*dxqsi_c(i,k);

                } // end d

              } // end j

            } // end k

          }// end c
        } // endi


      } // fim quadratura


      // Projection - to force non-penetrarion bc
      mesh->getCellNodesId(&*cell, cell_nodes.data());
      getProjectorBC(Prj, nodes_per_cell, cell_nodes.data(), Vec_x_0, current_time, *this /*AppCtx*/);

      Floc = Prj*Floc;
      Aloc = Prj*Aloc*Prj;

#ifdef FEP_HAS_OPENMP
      FEP_PRAGMA_OMP(critical)
#endif
      {
        VecSetValues(Vec_fun, mapV_c.size(), mapV_c.data(), Floc.data(), ADD_VALUES);
        MatSetValues(*JJ, mapV_c.size(), mapV_c.data(), mapV_c.size(), mapV_c.data(), Aloc.data(),  ADD_VALUES);
      }
    } // end cell loop


  } // end parallel

  
  // boundary conditions on global Jacobian
    // solid & triple tags .. force normal
  if (force_dirichlet)
  {
    int      nodeid;
    int      v_dofs[dim];
    Vector   normal(dim);
    Tensor   A(dim,dim);
    Tensor   I(Tensor::Identity(dim,dim));
    int      tag;

    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    for ( ; point != point_end; ++point)
    {
      tag = point->getTag();
      if (!(is_in(tag,feature_tags)   ||
            is_in(tag,solid_tags)     ||
            is_in(tag,interface_tags) ||
            is_in(tag,triple_tags)    ||
            is_in(tag,dirichlet_tags) ||
            is_in(tag,neumann_tags)   ||
            is_in(tag,periodic_tags)     )  )
        continue;
      //dof_handler[DH_UNKS].getVariable(VAR_U).getVertexAssociatedDofs(v_dofs, &*point);
      getNodeDofs(&*point, DH_MESH, VAR_M, v_dofs);
      
      nodeid = mesh->getPointId(&*point);
      getProjectorBC(A, 1, &nodeid, Vec_x_0, current_time, *this);
      A = I - A;
      MatSetValues(*JJ, dim, v_dofs, dim, v_dofs, A.data(), ADD_VALUES);
    }
  }
  
  Assembly(*JJ);
  Assembly(Vec_fun);

  PetscFunctionReturn(0);
}


PetscErrorCode AppCtx::formJacobian_mesh(SNES /*snes*/,Vec /*Vec_up_k*/,Mat* /**Mat_Jac*/, Mat* /*prejac*/, MatStructure * /*flag*/)
{
  // jacobian matrix is done in the formFunction_mesh
  PetscFunctionReturn(0);
}














































