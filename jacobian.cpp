#include "common.hpp"

// ******************************************************************************
//                            FORM JACOBIAN
// ******************************************************************************
PetscErrorCode AppCtx::formJacobian(SNES /*snes*/,Vec x,Mat *Jac, Mat* /*prejac*/, MatStructure * /*flag*/)
{

  Mat                 *JJ = Jac;
  PetscErrorCode      ierr;

  int iter;
  SNESGetIterationNumber(snes, &iter);


  MatZeroEntries(*JJ);


  //  LOOP NOS ELEMENTOS
  //
  //#pragma omp parallel default(none) shared(JJ,x,cout)
  {

    //MatrixXd            u_coefs_c(n_dofs_u_per_cell/dim, dim);
    MatrixXd            u_coefs_c_trans(dim,n_dofs_u_per_cell/dim);       // n+1 também
    MatrixXd            u_coefs_c_old(n_dofs_u_per_cell/dim, dim);        // n
    MatrixXd            u_coefs_c_old_trans(dim,n_dofs_u_per_cell/dim);   //
    MatrixXd            v_coefs_c_old(nodes_per_cell, dim);        // convective velocity; n
    MatrixXd            v_coefs_c_old_trans(dim,nodes_per_cell);   //
    MatrixXd            u_coefs_c_new(n_dofs_u_per_cell/dim, dim);        // n+1
    MatrixXd            u_coefs_c_new_trans(dim,n_dofs_u_per_cell/dim);   //
    VectorXd            p_coefs_c(n_dofs_p_per_cell);                  // valores de P na célula
    MatrixXd            x_coefs_c(nodes_per_cell, dim);                // coordenadas nodais da célula
    MatrixXd            x_coefs_c_trans(dim, nodes_per_cell);          // coordenadas nodais da célula
    Tensor              F_c(dim,dim);
    Tensor              invF_c(dim,dim);
    Tensor              invFT_c(dim,dim);
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
    Vector              Ubqp(dim); // bble
    //double              Pqp;
    Vector              Uqp_old(dim);
    Vector              Uconv_old(dim);
    VectorXi            cell_nodes(nodes_per_cell);
    double              Jx;
    double              weight;
    double              visc;
    double              cell_volume; // or area, 2d case
    double              hk2;
    double              tau=0;
    double              bble_integ=0;
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

    double              delta_cd;
    int                 tag;

    VectorXi            mapU_c(n_dofs_u_per_cell);
    VectorXi            mapU_r(n_dofs_u_per_corner);
    VectorXi            mapP_c(n_dofs_p_per_cell);
    VectorXi            mapP_r(n_dofs_p_per_corner);
    // mesh velocity    
    VectorXi            mapM_c(dim*nodes_per_cell);
    VectorXi            mapM_f(dim*nodes_per_facet);
    VectorXi            mapM_r(dim*nodes_per_corner);

    MatrixXd            R(n_dofs_u_per_cell,n_dofs_u_per_cell);
    MatrixXd            tmp(n_dofs_u_per_cell,n_dofs_u_per_cell);

    //const int tid = omp_get_thread_num();
    //const int nthreads = omp_get_num_threads();
    //const int n_cell_colors = mesh->numCellColors();
    //
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
        tag = cell->getTag();

        mesh->getCellNodesId(&*cell, cell_nodes.data());
        mesh->getNodesCoords(cell_nodes.begin(), cell_nodes.end(), x_coefs_c.data());
        x_coefs_c_trans = x_coefs_c.transpose();

        getRotationMatrix(R,cell_nodes,cell_nodes.size());

        // mapeamento do local para o global:
        //
        dof_handler_vars.getVariable(0).getCellDofs(mapU_c.data(), &*cell);
        dof_handler_vars.getVariable(1).getCellDofs(mapP_c.data(), &*cell);
        dof_handler_mesh.getVariable(0).getCellDofs(mapM_c.data(), &*cell);

        /*  Pega os valores das variáveis nos graus de liberdade */
        VecGetValues(x , mapU_c.size(), mapU_c.data(), u_coefs_c_new.data());
        VecGetValues(q0, mapU_c.size(), mapU_c.data(), u_coefs_c_old.data());
        VecGetValues(u_mesh, mapM_c.size(), mapM_c.data(), v_coefs_c_old.data());
        VecGetValues(x , mapP_c.size(), mapP_c.data(), p_coefs_c.data());

        // transformando para coordenada verdadeira
        rotate_RtA(R,u_coefs_c_new,tmp);
        rotate_RtA(R,u_coefs_c_old,tmp);
        rotate_RtA(R,v_coefs_c_old,tmp);

        v_coefs_c_old_trans = v_coefs_c_old.transpose();
        u_coefs_c_old_trans = u_coefs_c_old.transpose();
        u_coefs_c_trans = u_coefs_c_new.transpose();


        visc = niu(current_time, tag);

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
            F_c = x_coefs_c_trans * dLqsi_c[qp];
            Jx = determinant(F_c,dim);
            cell_volume += Jx * quadr_cell->weight(qp);
          }

          hk2 = cell_volume / pi; // element size
          tau = hk2/(4.*visc);

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
            F_c = x_coefs_c_trans * dLqsi_c[qp];
            Jx = determinant(F_c,dim);
            Xqp  = x_coefs_c_trans * qsi_c[qp];
            cell_volume += Jx * quadr_cell->weight(qp);
            Xc += Jx * quadr_cell->weight(qp) * Xqp;
          }
          Xc /= cell_volume;
        }
        //cout << "--------------------------------------------\n";
        for (int qp = 0; qp < n_qpts_cell; ++qp)
        {
          F_c    = x_coefs_c_trans * dLqsi_c[qp];
          inverseAndDet(F_c, dim, invF_c,Jx);
          invFT_c= invF_c.transpose();

          dxphi_c = dLphi_c[qp] * invF_c;
          dxpsi_c = dLpsi_c[qp] * invF_c;
          dxqsi_c = dLqsi_c[qp] * invF_c;

          dxP  = dxpsi_c.transpose() * p_coefs_c;
          dxU  = u_coefs_c_trans * dxphi_c;

          Xqp  = x_coefs_c_trans * qsi_c[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
          Uqp  = u_coefs_c_trans * phi_c[qp];
          //Pqp  = p_coefs_c.dot(psi_c[qp]);
          Uqp_old = u_coefs_c_old_trans * phi_c[qp];
          Uconv_old = Uqp_old - v_coefs_c_old_trans * qsi_c[qp];
          //Uconv_old = Uqp_old;

          weight = quadr_cell->weight(qp);
          if (Jx < 1.e-10)
          {
            ////#pragma omp critical
            {
              std::cout << "erro: jacobiana da integral não invertível: ";
              std::cout << "Jx = " << Jx << endl;
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
                  Aloc(i*dim + c, j*dim + d) += Jx*weight*
                                                ( delta_cd*phi_c[qp][i]* pho(Xqp,tag)*(phi_c[qp][j]/dt + has_convec*theta * Uconv_old.dot(dxphi_c.row(j)))   // advecção
                                                  +theta*visc*( delta_cd * dxphi_c.row(i).dot(dxphi_c.row(j)) + dxphi_c(i,d)*dxphi_c(j,c))   ); // rigidez

                }
              }
              for (int j = 0; j < n_dofs_p_per_cell; ++j)
                Gloc(i*dim + c,j) -= Jx*weight*psi_c[qp][j]*dxphi_c(i,c);
            }

          }


          /* ------------------------------
            *        stabilization
            * ------------------------------ */
          if (behaviors & (BH_bble_condens_PnPn | BH_bble_condens_CR))
          {
            dxbble = invFT_c * dLbble[qp];

            for (int c = 0; c < dim; ++c)
            {
              for (int j = 0; j < n_dofs_u_per_cell/dim; j++)
              {
                for (int d = 0; d < dim; d++)
                {
                  delta_cd = c==d;
                  Bbn(c, j*dim + d) += Jx*weight*
                                        (  delta_cd *bble[qp] *pho(Xqp,tag)*(phi_c[qp][j]/dt + has_convec*theta*Uconv_old.dot(dxphi_c.row(j)))  // advecção
                                        + theta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) );    // rigidez

                  Bnb(j*dim + d, c) += Jx*weight*
                                        (  delta_cd * phi_c[qp][j] *pho(Xqp,tag)*(bble[qp]/dt + has_convec*theta*Uconv_old.dot(dxbble))  // advecção
                                        + theta*visc*(delta_cd * dxphi_c.row(j).dot(dxbble) + dxphi_c(j,c)*dxbble(d)) );    // rigidez
                }
              }
              for (int j = 0; j < n_dofs_p_per_cell; ++j)
                Gbp(c, j) -= Jx*weight*psi_c[qp][j]*dxbble(c);
            }


            for (int c = 0; c < dim; c++)
            {
              for (int d = 0; d < dim; d++)
              {
                delta_cd = c==d;
                iBbb(c, d) += Jx*weight*
                            ( bble[qp]* delta_cd *pho(Xqp,tag)*(bble[qp]/dt+ has_convec*theta*Uconv_old.dot(dxbble) ) // advecção
                              +theta*visc*(delta_cd* dxbble.dot(dxbble)  + dxbble(d)*dxbble(c)) ); // rigidez
              }

            }

          }
          else
          if(behaviors & BH_GLS)
          {
            for (int i = 0; i < n_dofs_p_per_cell; i++)
              for (int j = 0; j < n_dofs_u_per_cell/dim; j++)
                for (int d = 0; d < dim; d++)
                {
                  Dloc(i, j*dim + d) -= Jx*weight* tau*dxpsi_c(i,d)*pho(Xqp,tag)*( phi_c[qp][j]/dt + has_convec*theta*Uconv_old.dot(dxphi_c.row(j)) );

                  Cloc(j*dim + d,i) += Jx*weight *tau* pho(Xqp,tag)*has_convec* Uconv_old.dot(dxphi_c.row(j)) * dxpsi_c(i,d);
                }

            for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
              for (int c = 0; c < dim; c++)
                for (int j = 0; j < n_dofs_u_per_cell/dim; ++j)
                  for (int d = 0; d < dim; ++d)
                  {
                    delta_cd = c==d;
                    Aloc(i*dim + c, j*dim + d) += Jx*weight*tau*delta_cd*has_convec*pho(Xqp,tag)*Uconv_old.dot(dxphi_c.row(i))*pho(Xqp,tag)*(phi_c[qp][j]/dt + theta*has_convec*Uconv_old.dot(dxphi_c.row(j)));
                  }

            for (int i = 0; i < n_dofs_p_per_cell; ++i)
              for (int j = 0; j < n_dofs_p_per_cell; ++j)
                Eloc(i,j) -= tau*Jx*weight * dxphi_c.row(i).dot(dxphi_c.row(j));
          }

          if (behaviors & BH_bble_condens_CR)
          {
            bble_integ += Jx*weight*bble[qp];

            for (int c = 0; c < dim; ++c)
              for (int i = 0; i < n_dofs_u_per_cell/dim; ++i)
                for (int j = 0; j < dim; ++j) // pressure gradient
                  Gnx(i*dim + c,j) -= Jx*weight* (Xqp(j) - Xc(j))*dxphi_c(i,c);
          }


        } // fim quadratura

        Dloc += Gloc.transpose();

        /* estabilização */
        if (behaviors & BH_bble_condens_PnPn)
        {
          invert(iBbb,dim);

          Dpb  = Gbp.transpose();

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
          Aloc += (  (Gnx*iBbb/bble_integ - bble_integ*Bnb)*Gnx.transpose() - Gnx*Bbn )/bble_integ;
        }

        rotate_RARt(R, Aloc, tmp);
        rotate_RA(R, Gloc, tmp);
        rotate_ARt(R, Dloc,tmp);
        MatSetValues(*JJ, mapU_c.size(), mapU_c.data(), mapU_c.size(), mapU_c.data(), Aloc.data(),  ADD_VALUES);
        MatSetValues(*JJ, mapU_c.size(), mapU_c.data(), mapP_c.size(), mapP_c.data(), Gloc.data(),  ADD_VALUES);
        MatSetValues(*JJ, mapP_c.size(), mapP_c.data(), mapU_c.size(), mapU_c.data(), Dloc.data(),  ADD_VALUES);
        if (pres_pres_block)
        {
          MatSetValues(*JJ, mapP_c.size(), mapP_c.data(), mapP_c.size(), mapP_c.data(), Eloc.data(),  ADD_VALUES);
        }

      }

      //#pragma omp barrier
    //} // endl color



  } // end parallel

  // LOOP NAS FACETS
  if (interface_tags.size() != 0)
  //#pragma omp parallel default(none)  shared(JJ,x,cout)
  {

    int                tag;
    PetscBool          is_surface;
    VectorXi           facet_nodes(nodes_per_facet);
    double             Jx;
    double             weight;
    MatrixXd           x_coefs_f(nodes_per_facet, dim);                // coordenadas nodais da célula
    MatrixXd           x_coefs_f_trans(dim, nodes_per_facet);
    Vector             normal(dim);
    Tensor             F_f(dim,dim-1);
    Tensor             invF_f(dim-1,dim);
    MatrixXd           dxphi_f(n_dofs_u_per_facet/dim, dim);
    Vector             Xqp(dim);
    MatrixXd           Aloc_f(n_dofs_u_per_facet, n_dofs_u_per_facet);

    VectorXi           mapU_f(n_dofs_u_per_facet);
    VectorXi           mapP_f(n_dofs_p_per_facet);

    MatrixXd           R(n_dofs_u_per_facet,n_dofs_u_per_facet);
    MatrixXd           tmp;

    // LOOP NAS FACES DO CONTORNO
    facet_iterator facet = mesh->facetBegin();
    facet_iterator facet_end = mesh->facetEnd();
    for (; facet != facet_end; ++facet)
    {
      tag = facet->getTag();

      //is_neumann = (neumann_tags.end() != std::find(neumann_tags.begin(), neumann_tags.end(), tag));
      is_surface = (PetscBool)is_in(tag, interface_tags);

      if ( (!is_surface) ) continue;

      dof_handler_vars.getVariable(0).getFacetDofs(mapU_f.data(), &*facet);
      dof_handler_vars.getVariable(1).getFacetDofs(mapP_f.data(), &*facet);

      mesh->getFacetNodesId(&*facet, facet_nodes.data());
      mesh->getNodesCoords(facet_nodes.begin(), facet_nodes.end(), x_coefs_f.data());
      x_coefs_f_trans = x_coefs_f.transpose();

      getRotationMatrix(R,facet_nodes,facet_nodes.size());

      ///*  Pega os valores das variáveis nos graus de liberdade */
      //VecGetValues(x, mapUb.size(), mapUb.data(), halfU_trans.data());

      Aloc_f.setZero();

      //dxUb  = halfU_trans * dxphib; // n+1


      for (int qp = 0; qp < n_qpts_facet; ++qp)
      {
        F_f   = x_coefs_f_trans * dLqsi_f[qp];

        tmp.resize(dim-1,dim-1);
        tmp = F_f.transpose()*F_f;
        Jx = sqrt(tmp.determinant());
        tmp = tmp.inverse();
        invF_f = tmp*F_f.transpose();


        weight = quadr_facet->weight(qp);
        Xqp  = x_coefs_f_trans * qsi_f[qp]; // coordenada espacial (x,y,z) do ponto de quadratura
        dxphi_f = dLphi_f[qp] * invF_f;


        for (int i = 0; i < n_dofs_u_per_facet/dim; ++i)
        {
          for (int j = 0; j < n_dofs_u_per_facet/dim; ++j)
          {
            for (int c = 0; c < dim; ++c)
            {
              Aloc_f(i*dim + c, j*dim + c) += Jx*weight* (unsteady*dt) *gama(Xqp,current_time,tag)*dxphi_f.row(i).dot(dxphi_f.row(j));
            }
          }
        }
      } // end quadratura

      rotate_RARt(R, Aloc_f, tmp);
      MatSetValues(*JJ, mapU_f.size(), mapU_f.data(), mapU_f.size(), mapU_f.data(), Aloc_f.data(),  ADD_VALUES);


    }

  } // end parallel


  //// LOOP CORNERS
  //#pragma omp parallel shared(x,JJ,cout) default(none)
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
    volatile    double hh;

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
      dof_handler_vars.getVariable(0).getCornerDofs(mapU_r.data(), &*corner);
      dof_handler_vars.getVariable(1).getCornerDofs(mapP_r.data(), &*corner);

      VecGetValues(x , mapU_r.size(), mapU_r.data(), u_coefs_k.data());

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
      MatSetValues(*JJ, mapU_r.size(), mapU_r.data(), mapU_r.size(), mapU_r.data(), Aloc_r.data(),  ADD_VALUES);

    }




  }



  Assembly(*JJ);

  // impondo as condições de contorno
  if (force_dirichlet) {
    ierr = MatZeroRows(*JJ, dir_entries.size(), dir_entries.data(), 1., NULL, NULL);  CHKERRQ(ierr);
  }
  Assembly(*JJ);


  //*flag = SAME_NONZERO_PATTERN;

  if (print_to_matlab) {
    static bool ja_foi=false;
    if (!ja_foi) {View(*JJ,"jacob.m","Jac"); cout << "imprimiu matlab" << dir_entries.size() <<  endl;}
    ja_foi = true;
  }
  PetscFunctionReturn(0);


} // end form jacobian



