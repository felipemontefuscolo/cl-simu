#include "common.hpp"

// ******************************************************************************
//                            FORM JACOBIAN
// ******************************************************************************
PetscErrorCode AppCtx::formJacobian(SNES /*snes*/,Vec Vec_up_k,Mat *Mat_Jac, Mat* /*prejac*/, MatStructure * /*flag*/)
{

  Mat                 *JJ = Mat_Jac;
  PetscErrorCode      ierr;

  int iter;
  SNESGetIterationNumber(snes, &iter);



  MatZeroEntries(*JJ);


  //  LOOP NOS ELEMENTOS
  //
#if (FEP_HAS_OPENMP)
  #pragma omp parallel default(none) shared(JJ,Vec_up_k,cout)
#endif
  {

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
    //const int n_cell_colors = mesh->numCellColors();
    //
    //cell_iterator cell;
    //cell_iterator cell_end;

    //for (int color = 0; color < n_cell_colors; ++color)
    //{
      //cell = mesh->cellBegin(EColor(color),tid,nthreads);
      //cell_end = mesh->cellEnd(EColor(color),tid,nthreads);

      cell_iterator cell = mesh->cellBegin(tid,nthreads);
      cell_iterator cell_end = mesh->cellEnd(tid,nthreads);

      //cell_iterator cell = mesh->cellBegin();
      //cell_iterator cell_end = mesh->cellEnd();
      for (; cell != cell_end; ++cell)
      {
        tag = cell->getTag();

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
        //mesh->getCellNodesId(&*cell, cell_nodes.data());
        //mesh->getNodesCoords(cell_nodes.begin(), cell_nodes.end(), x_coefs_c.data());
        //x_coefs_c_trans = x_coefs_c.transpose();

        v_coefs_c_mid_trans = v_coefs_c_mid.transpose();
        x_coefs_c_old_trans = x_coefs_c_old.transpose();
        x_coefs_c_new_trans = x_coefs_c_new.transpose();
        u_coefs_c_old_trans = u_coefs_c_old.transpose();
        u_coefs_c_new_trans = u_coefs_c_new.transpose();

        u_coefs_c_mid_trans = utheta*u_coefs_c_new_trans + (1.-utheta)*u_coefs_c_old_trans;
        x_coefs_c_mid_trans = utheta*x_coefs_c_new_trans + (1.-utheta)*x_coefs_c_old_trans;

        //  SetUp
        visc = muu(tag);
        rho  = pho(Xqp,tag);
        Aloc.setZero();
        Gloc.setZero();
        Dloc.setZero();

        if (behaviors & BH_bble_condens_PnPn)
        {
          iBbb.setZero();
          Bnb.setZero();
          Bbn.setZero();
          Dpb.setZero();
          Gbp.setZero();
          FUb.setZero();
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
          iBbb.setZero(); // it is not the inverse to CR element
          Gnx.setZero();
          Bnb.setZero();
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

        for (int qp = 0; qp < n_qpts_cell; ++qp)
        {
          //F_c_new = x_coefs_c_new_trans * dLqsi_c[qp];
          //inverseAndDet(F_c_new,dim,invF_c_new,J_new);
          //invFT_c_new= invF_c_new.transpose();

          //F_c_old = x_coefs_c_old_trans * dLqsi_c[qp];
          //inverseAndDet(F_c_old,dim,invF_c_old,J_old);
          //invFT_c_old= invF_c_old.transpose();

          F_c_mid    = x_coefs_c_mid_trans * dLqsi_c[qp];
          inverseAndDet(F_c_mid, dim, invF_c_mid,J_mid);
          invFT_c_mid= invF_c_mid.transpose();

          dxphi_c = dLphi_c[qp] * invF_c_mid;
          dxpsi_c = dLpsi_c[qp] * invF_c_mid;
          dxqsi_c = dLqsi_c[qp] * invF_c_mid;

          dxP_new  = dxpsi_c.transpose() * p_coefs_c_new;
          dxU  = u_coefs_c_mid_trans * dxphi_c;

          Xqp      = x_coefs_c_mid_trans * qsi_c[qp];     // coordenada espacial (x,y,z) do ponto de quadratura
          Uqp      = u_coefs_c_mid_trans * phi_c[qp];     //n+utheta
          Uqp_old  = u_coefs_c_old_trans * phi_c[qp];
          Uqp_new  = u_coefs_c_new_trans * phi_c[qp]; //n+1
          Vqp      = v_coefs_c_mid_trans * qsi_c[qp];
          Uconv_qp = Uqp - Vqp;
          dUdt     = (Uqp_new-Uqp_old)/dt;
          Pqp_new  = p_coefs_c_new.dot(psi_c[qp]);

          force_at_mid = force(Xqp,current_time+utheta*dt,tag);

          weight = quadr_cell->weight(qp);
          JxW_mid = J_mid * weight;
          if (J_mid < 1.e-10)
          {
            ////#pragma omp critical
            {
              std::cout << "erro: jacobiana da integral não invertível: ";
              std::cout << "J_mid = " << J_mid << endl;
            }
            throw;
          }

          for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
          {
            for (int c = 0; c < dim; ++c)
            {
              for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
              {
                for (int d = 0; d < dim; ++d)
                {
                  delta_cd = c==d;
                  Aloc(i*dim + c, j*dim + d) += JxW_mid*
                                                ( has_convec*phi_c[qp][i]*utheta *rho*( delta_cd*Uconv_qp.dot(dxphi_c.row(j))  +  dxU(c,d)*phi_c[qp][j] )   // advecção
                                                + delta_cd*rho*phi_c[qp][i]*phi_c[qp][j]/dt     // time derivative
                                                + utheta*visc*( delta_cd * dxphi_c.row(i).dot(dxphi_c.row(j)) + dxphi_c(i,d)*dxphi_c(j,c))   ); // rigidez

                }
              }
              for (int j = 0; j < n_dofs_p_per_cell; ++j)
                Gloc(i*dim + c,j) -= JxW_mid * psi_c[qp][j]* dxphi_c(i,c);
            }

          }

          //+
          //+ ------------------------------
          //+        stabilization
          //+ ------------------------------
          //+
          if (behaviors & (BH_bble_condens_PnPn | BH_bble_condens_CR))
          {
            dxbble = invFT_c_mid * dLbble[qp];

            for (int c = 0; c < dim; ++c)
            {
              for (int j = 0; j < n_dofs_u_per_cell/dim; j++)
              {
                for (int d = 0; d < dim; d++)
                {
                  delta_cd = c==d;
                  Bbn(c, j*dim + d) += JxW_mid*
                                       ( has_convec*bble[qp]*utheta *rho*( delta_cd*Uconv_qp.dot(dxphi_c.row(j))  +  dxU(c,d)*phi_c[qp][j] )   // convective
                                       + delta_cd*rho*bble[qp]*phi_c[qp][j]/dt     // time derivative
                                       + utheta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) );    // rigidez

                  Bnb(j*dim + d, c) += JxW_mid*
                                       ( has_convec*phi_c[qp][j]*utheta *rho*(  delta_cd*Uconv_qp.dot(dxbble)  )   // convective
                                       + delta_cd*rho*phi_c[qp][j]*bble[qp]/dt     // time derivative
                                       + utheta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) );    // rigidez

                }
              }
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
                              + delta_cd*rho*bble[qp]*bble[qp]/dt     // time derivative
                              + utheta*visc*(delta_cd* dxbble.dot(dxbble)  + dxbble(d)*dxbble(c)) ); // rigidez
              }

              FUb(c) += JxW_mid*
                        ( bble[qp]*rho*(dUdt(c) + has_convec*Uconv_qp.dot(dxU.row(c))) + // time derivative + convective
                          visc*dxbble.dot(dxU.row(c) + dxU.col(c).transpose()) - //rigidez
                          Pqp_new*dxbble(c) -                     // pressão
                          force_at_mid(c)*bble[qp]   ); // força

            }

          }
          else
          if(behaviors & BH_GLS)
          {
            Res = rho*( dUdt +  has_convec*dxU*Uconv_qp) + dxP_new - force_at_mid;

            for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
            {
              dResdu = (rho*phi_c[qp][j]/dt)*I + has_convec*rho*utheta*(  phi_c[qp][j]*dxU + Uconv_qp.dot(dxphi_c.row(j))*I );

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
          }

          if (behaviors & BH_bble_condens_CR)
          {
            bble_integ += JxW_mid*bble[qp];

            for (int c = 0; c < dim; ++c)
              for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
                for (int j = 0; j < dim; ++j) // pressure gradient
                  Gnx(i*dim + c,j) -= JxW_mid* (Xqp(j) - Xc(j))*dxphi_c(i,c);
          }


        } // fim quadratura

        Dloc += utheta*Gloc.transpose();

        /* estabilização */

        if (behaviors & BH_bble_condens_PnPn)
        {
          invert(iBbb,dim);

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
        }


        // Projection - to force non-penetrarion bc
        mesh->getCellNodesId(&*cell, cell_nodes.data());
        getProjectorMatrix(Prj, nodes_per_cell, cell_nodes.data(), Vec_x_1, current_time+dt);

        Aloc = Prj*Aloc;
        Gloc = Prj*Gloc;
#if (FEP_HAS_OPENMP)
        #pragma omp critical
#endif
        {
          MatSetValues(*JJ, mapU_c.size(), mapU_c.data(), mapU_c.size(), mapU_c.data(), Aloc.data(),  ADD_VALUES);
          MatSetValues(*JJ, mapU_c.size(), mapU_c.data(), mapP_c.size(), mapP_c.data(), Gloc.data(),  ADD_VALUES);
          MatSetValues(*JJ, mapP_c.size(), mapP_c.data(), mapU_c.size(), mapU_c.data(), Dloc.data(),  ADD_VALUES);
          if (pres_pres_block)
          {
            MatSetValues(*JJ, mapP_c.size(), mapP_c.data(), mapP_c.size(), mapP_c.data(), Eloc.data(),  ADD_VALUES);
          }
        }
      }

      //#pragma omp barrier
    //} // endl color



  } // end parallel

  // LOOP NAS FACETS
  //#pragma omp parallel default(none)  shared(JJ,Vec_up_k,cout)
  {

    int                tag;
    bool               is_surface;
    bool               is_solid;
    bool               is_neumann;
    double             J_mid;
    double             weight;
    double             JxW_mid;
    MatrixXd           x_coefs_f_old(nodes_per_facet, dim);                // coordenadas nodais
    MatrixXd           x_coefs_f_new(nodes_per_facet, dim);                // coordenadas nodais
    MatrixXd           x_coefs_f_new_trans(dim, nodes_per_facet);
    MatrixXd           x_coefs_f_old_trans(dim, nodes_per_facet);
    MatrixXd           x_coefs_f_mid_trans(dim, nodes_per_facet);
    Vector             normal(dim);
    Tensor             F_f_mid(dim,dim-1);
    Tensor             invF_f_mid(dim-1,dim);
    MatrixXd           dxphi_f(n_dofs_u_per_facet/dim, dim);
    Vector             Xqp(dim);
    MatrixXd           Aloc_f(n_dofs_u_per_facet, n_dofs_u_per_facet);

    MatrixXd           tmp(n_dofs_u_per_facet,n_dofs_u_per_facet);

    VectorXi           mapU_f(n_dofs_u_per_facet);
    VectorXi           mapP_f(n_dofs_p_per_facet);
    VectorXi           mapM_f(dim*nodes_per_facet);

    MatrixXd           Prj(n_dofs_u_per_facet,n_dofs_u_per_facet);
    VectorXi           facet_nodes(nodes_per_facet);


    //~ const int tid = omp_get_thread_num();
    //~ const int nthreads = omp_get_num_threads();
//~ 
    //~ facet_iterator facet = mesh->facetBegin(tid,nthreads);
    //~ facet_iterator facet_end = mesh->facetEnd(tid,nthreads);

    // LOOP NAS FACES DO CONTORNO
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

      dof_handler[DH_UNKS].getVariable(VAR_U).getFacetDofs(mapU_f.data(), &*facet);
      dof_handler[DH_UNKS].getVariable(VAR_P).getFacetDofs(mapP_f.data(), &*facet);
      dof_handler[DH_MESH].getVariable(VAR_M).getFacetDofs(mapM_f.data(), &*facet);

      VecGetValues(Vec_x_0,     mapM_f.size(), mapM_f.data(), x_coefs_f_old.data());
      VecGetValues(Vec_x_1,     mapM_f.size(), mapM_f.data(), x_coefs_f_new.data());      
      x_coefs_f_old_trans = x_coefs_f_old.transpose();
      x_coefs_f_new_trans = x_coefs_f_new.transpose();
      x_coefs_f_mid_trans = utheta*x_coefs_f_new_trans + (1.-utheta)*x_coefs_f_old_trans;

      ///*  Pega os valores das variáveis nos graus de liberdade */
      //VecGetValues(Vec_up_k, mapUb.size(), mapUb.data(), halfU_trans.data());

      Aloc_f.setZero();

      //dxUb  = halfU_trans * dxphib; // n+1


      for (int qp = 0; qp < n_qpts_facet; ++qp)
      {
        F_f_mid   = x_coefs_f_mid_trans * dLqsi_f[qp];

        tmp.resize(dim-1,dim-1);
        tmp = F_f_mid.transpose()*F_f_mid;
        J_mid = sqrt(tmp.determinant());
        invF_f_mid = tmp.inverse()*F_f_mid.transpose();


        weight = quadr_facet->weight(qp);
        JxW_mid = J_mid * weight;
        Xqp  = x_coefs_f_mid_trans * qsi_f[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        dxphi_f = dLphi_f[qp] * invF_f_mid;

        if (false)
        if (is_surface) // semi-implicit term
        {
          for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
            for (int j = 0; j < n_dofs_u_per_facet/dim; ++j)
              for (int c = 0; c < dim; ++c)
                Aloc_f(i*dim + c, j*dim + c) += utheta*JxW_mid* (unsteady*dt) *gama(Xqp,current_time,tag)*dxphi_f.row(i).dot(dxphi_f.row(j));
        }

        if (is_solid) // dissipation
        {
          for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
            for (int j = 0; j < n_dofs_u_per_facet/dim; ++j)
              for (int c = 0; c < dim; ++c)
                Aloc_f(i*dim + c, j*dim + c) += utheta*JxW_mid *beta_diss()*phi_f[qp][j]*phi_f[qp][i];
                //Aloc_f(i*dim + c, j*dim + c) += x_coefs_f_old_trans.norm()*utheta*beta_diss()*phi_f[qp][j]*phi_f[qp][i];
        }



      } // end quadratura

      // Projection - to force non-penetrarion bc
      mesh->getFacetNodesId(&*facet, facet_nodes.data());
      getProjectorMatrix(Prj, nodes_per_facet, facet_nodes.data(), Vec_x_1, current_time+dt);

      Aloc_f = Prj*Aloc_f;
      //~ #pragma omp critical
      {
        MatSetValues(*JJ, mapU_f.size(), mapU_f.data(), mapU_f.size(), mapU_f.data(), Aloc_f.data(),  ADD_VALUES);
      }


    }

  } // end parallel

  //// LOOP CORNERS
  //#pragma omp parallel shared(Vec_up_k,JJ,cout) default(none)
  //if (false)
  {
    Real const eps = std::numeric_limits<Real>::epsilon();
    Real const eps_root = pow(eps,1./3.);

    int                 tag;
    bool                is_triple;

    MatrixXd          u_coefs_k(n_dofs_u_per_corner/dim, dim);
    MatrixXd          u_coefs_kp1(n_dofs_u_per_corner/dim, dim);
    MatrixXd          u_coefs_km1(n_dofs_u_per_corner/dim, dim);
    VectorXd          FUloc_km1(n_dofs_u_per_corner);
    VectorXd          FUloc_kp1(n_dofs_u_per_corner);
    MatrixXd          Aloc_r(n_dofs_u_per_corner, n_dofs_u_per_corner);
    VectorXi          mapU_r(n_dofs_u_per_corner);
    VectorXi          mapP_r(n_dofs_p_per_corner);
    double            h;
    volatile double   hh;

    MatrixXd          Prj(n_dofs_u_per_corner,n_dofs_u_per_corner);
    VectorXi          corner_nodes(nodes_per_corner);

    //const int tid = omp_get_thread_num();
    //const int nthreads = omp_get_num_threads();

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

      VecGetValues(Vec_up_k , mapU_r.size(), mapU_r.data(), u_coefs_k.data());

      //VecSetValues(f, mapU_r.size(), mapU_r.data(), FUloc.data(), ADD_VALUES);


      h = max(u_coefs_k.norm(),1.)*eps_root;

      for (int i = 0; i < n_dofs_u_per_corner; ++i)
      {
        u_coefs_kp1 = u_coefs_k;
        u_coefs_km1 = u_coefs_k;

        u_coefs_kp1(i) += h;
        u_coefs_km1(i) -= h;
        hh = u_coefs_kp1(i) - u_coefs_km1(i);

        formCornerFunction(corner,mapU_r,mapP_r,u_coefs_kp1,FUloc_kp1);
        formCornerFunction(corner,mapU_r,mapP_r,u_coefs_km1,FUloc_km1);

        for (int l=0; l<n_dofs_u_per_corner; ++l)
        {
          Aloc_r(l,i) = (FUloc_kp1(l)-FUloc_km1(l))/hh;
        }

      }

      // Projection - to force non-penetrarion bc
      mesh->getCornerNodesId(&*corner, corner_nodes.data());
      getProjectorMatrix(Prj, nodes_per_corner, corner_nodes.data(), Vec_x_1, current_time+dt);

      Aloc_r = Prj*Aloc_r;
      MatSetValues(*JJ, mapU_r.size(), mapU_r.data(), mapU_r.size(), mapU_r.data(), Aloc_r.data(),  ADD_VALUES);

    }
  } // end parallel

  // solid & triple tags .. force normal
  if (force_dirichlet)
  {
    int      nodeid;
    int      u_dofs[dim];
    Vector   normal(dim);
    Tensor   A(dim,dim);
    Tensor   I(Tensor::Identity(dim,dim));

    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    for ( ; point != point_end; ++point)
    {
      if (!is_in(point->getTag(),solid_tags) && !is_in(point->getTag(),triple_tags))
        continue;
      //dof_handler[DH_UNKS].getVariable(VAR_U).getVertexAssociatedDofs(u_dofs, &*point);
      getNodeDofs(&*point, DH_UNKS, VAR_U, u_dofs);
      
      nodeid = mesh->getPointId(&*point);
      getProjectorMatrix(A, 1, &nodeid, Vec_x_1, current_time+dt);
      A = I - A;
      MatSetValues(*JJ, dim, u_dofs, dim, u_dofs, A.data(), ADD_VALUES);
      //if (nodeid==23)
      //{
      //  cout << "OLHA AQUI FDPPPPPPPPPPP\n";
      //  for (int j = 0; j < dim; ++j)
      //    cout << u_dofs[j] << " ";
      //  cout << endl;
      //}
    }
  }

  Assembly(*JJ);

  // impondo as condições de contorno
  if (force_dirichlet)
  {
    ierr = MatZeroRows(*JJ, dir_entries.size(), dir_entries.data(), 1., NULL, NULL);  CHKERRQ(ierr);
  }
  Assembly(*JJ);

  //*flag = SAME_NONZERO_PATTERN;

  if (print_to_matlab)
  {
    static bool ja_foi=false;
    if (!ja_foi) {View(*JJ,"jacob.m","Jac"); cout << "imprimiu matlab" << dir_entries.size() <<  endl;}
    ja_foi = true;
  }
  PetscFunctionReturn(0);


  //for (int i = 0; i < n_unknowns; ++i)
  //{
  //  double const a = 10000000;
  //  MatSetValues(*JJ, 1, &i, 1, &i, &a,  INSERT_VALUES);
  //}
  //Assembly(*JJ);


} // end form jacobian



