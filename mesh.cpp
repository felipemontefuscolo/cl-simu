#include "common.hpp"






void AppCtx::getVecNormals(Vec const* Vec_x_1, Vec & Vec_normal_)
{

  VectorXi          facet_nodes(nodes_per_facet);
  MatrixXd          x_coefs(nodes_per_facet, dim);                // coordenadas nodais da célula
  MatrixXd          x_coefs_trans(dim, nodes_per_facet);
  Vector            X(dim);
  Vector            normal(dim);
  Tensor            F(dim,dim-1);
  VectorXi          map(n_dofs_u_per_facet);
  //bool               is_surface, is_solid;
  int               tag;
  //int               tag_other;
  bool              virtual_mesh;

  if (Vec_x_1==NULL)
    virtual_mesh = false;
  else
    virtual_mesh = true;

  VecSet(Vec_normal_,0);

  // LOOP NAS FACES DO CONTORNO
  facet_iterator facet = mesh->facetBegin();
  facet_iterator facet_end = mesh->facetEnd();
  for (; facet != facet_end; ++facet)
  {
    tag = facet->getTag();


    if (is_in(tag, solid_tags) || !mesh->inBoundary(&*facet))
      continue;


    //is_surface = is_in(tag, interface_tags);
    //is_solid   = is_in(tag, solid_tags);
    //
    //if ( !is_surface && !is_solid)
    //  continue;

    dof_handler[DH_MESH].getVariable(VAR_M).getFacetDofs(map.data(), &*facet);

    mesh->getFacetNodesId(&*facet, facet_nodes.data());
    if (virtual_mesh)
      VecGetValues(*Vec_x_1, map.size(), map.data(), x_coefs.data());
    else
      mesh->getNodesCoords(facet_nodes.begin(), facet_nodes.end(), x_coefs.data());
    x_coefs_trans = x_coefs.transpose();


    if (false) // método alternativo que serve para todos os métodos (normal consistente)
    {
      // find the normal
      for (int k = 0; k < nodes_per_facet; ++k)
      {
        //tag_other = mesh->getNode(facet_nodes(k))->getTag();

        F   = x_coefs_trans * dLphi_nf[k];

        if (dim==2)
        {
          normal(0) = +F(1,0);
          normal(1) = -F(0,0);
        }
        else
        {
          normal = cross(F.col(0),F.col(1));
        }
        //normal.normalize();
        VecSetValues(Vec_normal_, dim, map.data()+k*dim, normal.data(), ADD_VALUES);

      } // nodes
            
    }
    else if (mesh_cell_type == TETRAHEDRON4) // facet = triangle3
    {
      F   = x_coefs_trans * dLphi_nf[0];
      Vector a(F.col(0)), b(F.col(1)-F.col(0)), c(F.col(1));
      
      normal = cross(F.col(0), F.col(1));
      
      // 0
      X = normal;
      X /= a.dot(a) * c.dot(c);
      VecSetValues(Vec_normal_, dim, map.data()+0*dim, X.data(), ADD_VALUES);
      
      // 1
      X = normal;
      X /= a.dot(a) * b.dot(b);
      VecSetValues(Vec_normal_, dim, map.data()+1*dim, X.data(), ADD_VALUES);
      
      // 2
      X = normal;
      X /= b.dot(b) * c.dot(c);
      VecSetValues(Vec_normal_, dim, map.data()+2*dim, X.data(), ADD_VALUES);

    }
    else if(mesh_cell_type == TETRAHEDRON10) // facet = triangle6
    {
      Vector x0(x_coefs_trans.col(0));
      Vector x1(x_coefs_trans.col(1));
      Vector x2(x_coefs_trans.col(2));
      Vector x3(x_coefs_trans.col(3));
      Vector x4(x_coefs_trans.col(4));
      Vector x5(x_coefs_trans.col(5));
      
      //edges
      double e03 = (x0-x3).squaredNorm();
      double e05 = (x0-x5).squaredNorm();
      double e13 = (x1-x3).squaredNorm();
      double e14 = (x1-x4).squaredNorm();
      double e24 = (x2-x4).squaredNorm();
      double e25 = (x2-x5).squaredNorm();
      double e34 = (x3-x4).squaredNorm();
      double e35 = (x3-x5).squaredNorm();
      double e45 = (x4-x5).squaredNorm();
      
      // dividi o triangulo em 4 mini-triangulos
      
      // node 0
      cross(normal,x3-x0,x5-x0);
      normal /= e03*e05;
      VecSetValues(Vec_normal_, dim, map.data()+0*dim, normal.data(), ADD_VALUES);
      
      // node 1
      cross(normal,x4-x1,x3-x1);
      normal /= e14*e13;
      VecSetValues(Vec_normal_, dim, map.data()+1*dim, normal.data(), ADD_VALUES);
      
      // node 2
      cross(normal,x5-x2,x4-x2);
      normal /= e25*e24;
      VecSetValues(Vec_normal_, dim, map.data()+2*dim, normal.data(), ADD_VALUES);
      
      // node 3
      normal = cross(x1-x3,x4-x3)/(e13*e34) + cross(x4-x3,x5-x3)/(e34*e35) + cross(x5-x3,x0-x3)/(e35*e03);
      VecSetValues(Vec_normal_, dim, map.data()+3*dim, normal.data(), ADD_VALUES);
      
      // node 4
      normal = cross(x2-x4,x5-x4)/(e24*e45) + cross(x5-x4,x3-x4)/(e45*e34) + cross(x3-x4,x1-x4)/(e34*e14);
      VecSetValues(Vec_normal_, dim, map.data()+4*dim, normal.data(), ADD_VALUES);
      
      // node 5
      normal = cross(x0-x5,x3-x5)/(e05*e35) + cross(x3-x5,x4-x5)/(e35*e45) + cross(x4-x5,x2-x5)/(e45*e25);
      VecSetValues(Vec_normal_, dim, map.data()+5*dim, normal.data(), ADD_VALUES);
    }
    else if(mesh_cell_type == TRIANGLE3)
    {
      normal(0) = x_coefs_trans(1,1)-x_coefs_trans(1,0);
      normal(1) = x_coefs_trans(0,0)-x_coefs_trans(0,1);
      
      normal /= (x_coefs_trans.col(0)-x_coefs_trans.col(1)).squaredNorm();
      VecSetValues(Vec_normal_, dim, map.data()+0*dim, normal.data(), ADD_VALUES);
      VecSetValues(Vec_normal_, dim, map.data()+1*dim, normal.data(), ADD_VALUES);
    }
    else if(mesh_cell_type == TRIANGLE6) // dividi a aresta em duas partes
    {
      normal(0) = x_coefs_trans(1,2)-x_coefs_trans(1,0);
      normal(1) = x_coefs_trans(0,0)-x_coefs_trans(0,2);
      normal /= (x_coefs_trans.col(0)-x_coefs_trans.col(2)).squaredNorm();
      VecSetValues(Vec_normal_, dim, map.data()+0*dim, normal.data(), ADD_VALUES);
      VecSetValues(Vec_normal_, dim, map.data()+2*dim, normal.data(), ADD_VALUES);
      
      normal(0) = x_coefs_trans(1,1)-x_coefs_trans(1,2);
      normal(1) = x_coefs_trans(0,2)-x_coefs_trans(0,1);
      normal /= (x_coefs_trans.col(1)-x_coefs_trans.col(2)).squaredNorm();
      VecSetValues(Vec_normal_, dim, map.data()+1*dim, normal.data(), ADD_VALUES);
      VecSetValues(Vec_normal_, dim, map.data()+2*dim, normal.data(), ADD_VALUES);
    }
    

  }
  Assembly(Vec_normal_);

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    tag = point->getTag();

    if (!mesh->inBoundary(&*point))
      continue;

    getNodeDofs(&*point, DH_MESH, VAR_M, map.data());

    if (!is_in(tag, solid_tags))// && !is_in(tag, triple_tags))
    {
      VecGetValues(Vec_normal_, dim, map.data(),normal.data());
      normal.normalize();
      VecSetValues(Vec_normal_, dim, map.data(),normal.data(), INSERT_VALUES);
      Assembly(Vec_normal_);
    }
    else
    {
      if (virtual_mesh)
        VecGetValues(*Vec_x_1, dim, map.data(), X.data());
      else
        point->getCoord(X.data(),dim);
      normal = -solid_normal(X,current_time,tag);

      VecSetValues(Vec_normal_, dim, map.data(), normal.data(), INSERT_VALUES);
    }

  }


}

/// @param[out] Vec_normal_
/// @param[out] Vec_x_
void AppCtx::smoothsMesh(Vec & Vec_normal_, Vec &Vec_x_)
{
  //int        nodeid;
  //double     *Vec_x__array;
  bool const u_has_edge_assoc_dof = shape_phi_c->numDofsAssociatedToFacet()+shape_phi_c->numDofsAssociatedToCorner() > 0;

  Vector      Xm(dim); // X mean
  Vector      Xi(dim);
  Vector      dX(dim);
  Vector      normal(dim);
  Vector      tmp(dim), tmp2(dim);
  Vector      Uf(dim), Ue(dim), Umsh(dim); // Ue := elastic velocity
  int         tag;
  int         iVs[128], *iVs_end;
  int         iCs[128], viCs[128];
  VectorXi    vtx_dofs_umesh(dim);  // indices de onde pegar a velocidade
  VectorXi    vtx_dofs_fluid(dim); // indices de onde pegar a velocidade
  VectorXi    edge_dofs_umesh(3*dim);
  VectorXi    edge_dofs_fluid((2+u_has_edge_assoc_dof)*dim);
  VectorXi    edge_nodes(3);
  double      error;
  double      old_quality, new_quality;
  bool        in_boundary;
  //int        id;


  /* suavização laplaciana */
  for (int smooth_it = 0; smooth_it < 3; ++smooth_it)
  {
    error = 0;
    getVecNormals(&Vec_x_, Vec_normal_);

    // VERTICES
    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    for (; point != point_end; ++point)
    {

      in_boundary = mesh->inBoundary(&*point);

      // pula o caso em que o ponto está no contoro mas não tem boundary smoothing
      if (!boundary_smoothing && in_boundary)
        continue;

      tag = point->getTag();

      if (mesh->isVertex(&*point))
      {
        //if (  is_in(tag,interface_tags) || is_in(tag,triple_tags) || is_in(tag,solid_tags) ||
        //    is_in(tag,dirichlet_tags) || is_in(tag,neumann_tags)  )
        if (is_in(tag,triple_tags))
          continue;

        Xm = Vector::Zero(dim);
        //iVs_end = mesh->connectedVtcs(&*point, iVs);
        iVs_end = mesh->connectedVtcs(&*point, iVs, iCs, viCs);

        if (!in_boundary)
        {
          for (int *it = iVs; it != iVs_end ; ++it)
          {
            dof_handler[DH_MESH].getVariable(VAR_M).getVertexDofs(vtx_dofs_umesh.data(), mesh->getNode(*it));
            VecGetValues(Vec_x_, dim, vtx_dofs_umesh.data(), tmp.data());
            Xm += tmp;
          }
          Xm /= (iVs_end-iVs);

          dof_handler[DH_MESH].getVariable(VAR_M).getVertexDofs(vtx_dofs_umesh.data(), &*point);
          //// compute error
          VecGetValues(Vec_x_, dim, vtx_dofs_umesh.data(), tmp.data()); // old coord
          error += (tmp-Xm).norm();
          old_quality = getCellPatchQuality(Vec_x_, iCs);

          VecSetValues(Vec_x_, dim, vtx_dofs_umesh.data(), Xm.data(), INSERT_VALUES);

          new_quality = getCellPatchQuality(Vec_x_, iCs);

          // se a qualidade piorou, volta no que estava antes
          if (new_quality < old_quality)
          {
            VecSetValues(Vec_x_, dim, vtx_dofs_umesh.data(), tmp.data(), INSERT_VALUES);
          }

        }
        else
        {
          int N=0;
          for (int *it = iVs; it != iVs_end ; ++it)
          {
            Point const* viz_pt = mesh->getNode(*it);
            if (viz_pt->getTag()!=tag && !is_in( viz_pt->getTag(), triple_tags ))
              continue;
            dof_handler[DH_MESH].getVariable(VAR_M).getVertexDofs(vtx_dofs_umesh.data(), mesh->getNode(*it));
            // debug
            VecGetValues(Vec_x_, dim, vtx_dofs_umesh.data(), tmp.data());
            ++N;
            Xm += tmp;

          }

          dof_handler[DH_MESH].getVariable(VAR_M).getVertexDofs(vtx_dofs_umesh.data(), &*point);
          VecGetValues(Vec_x_, dim, vtx_dofs_umesh.data(), Xi.data());

          if (dim==3)
            Xm = (N*Xi + 2*Xm)/(3*N);
            //Xm = Xm/N;
          else
            Xm = (N*Xi + Xm)/(2*N);
            //Xm = Xm/N;
          //

          dX = Xm - Xi;
          VecGetValues(Vec_normal_, dim, vtx_dofs_umesh.data(), normal.data());
          dX -= normal.dot(dX)*normal;
          Xi += dX;

          // compute error
          VecGetValues(Vec_x_, dim, vtx_dofs_umesh.data(), tmp.data());
          error += (tmp-Xi).norm();

          old_quality = getCellPatchQuality(Vec_x_, iCs);
          VecSetValues(Vec_x_, dim, vtx_dofs_umesh.data(), Xi.data(), INSERT_VALUES);
          new_quality = getCellPatchQuality(Vec_x_, iCs);
          // se a qualidade piorou, volta no que estava antes
          if (new_quality < old_quality)
          {
            VecSetValues(Vec_x_, dim, vtx_dofs_umesh.data(), tmp.data(), INSERT_VALUES);
          }
        }

      }

    } // end point


    // MID NODES
    point = mesh->pointBegin();
    point_end = mesh->pointEnd();
    if (u_has_edge_assoc_dof)
    for (; point != point_end; ++point)
    {

      in_boundary = mesh->inBoundary(&*point);

      // pula o caso em que o ponto está no contoro mas não tem boundary smoothing
      if (!boundary_smoothing && in_boundary)
        continue;

      tag = point->getTag();

      if (!mesh->isVertex(&*point))
      {

        const int m = point->getPosition() - mesh->numVerticesPerCell();
        Cell const* icell = mesh->getCell(point->getIncidCell());
        if (dim==3)
        {
          Corner *edge = mesh->getCorner(icell->getCornerId(m));
          dof_handler[DH_MESH].getVariable(VAR_M).getCornerDofs(edge_dofs_umesh.data(), &*edge);
          dof_handler[DH_UNKS].getVariable(VAR_U).getCornerDofs(edge_dofs_fluid.data(), &*edge);
          mesh->getCornerNodesId(&*edge, edge_nodes.data());
        }
        else // dim=2
        {
          Facet *edge = mesh->getFacet(icell->getFacetId(m));
          dof_handler[DH_MESH].getVariable(VAR_M).getFacetDofs(edge_dofs_umesh.data(), &*edge);
          dof_handler[DH_UNKS].getVariable(VAR_U).getFacetDofs(edge_dofs_fluid.data(), &*edge);
          mesh->getFacetNodesId(&*edge, edge_nodes.data());
        }

        VecGetValues(Vec_x_, dim, edge_dofs_umesh.data(), Xm.data());
        VecGetValues(Vec_x_, dim, edge_dofs_umesh.data()+dim, tmp.data());
        VecGetValues(Vec_x_, dim, edge_dofs_umesh.data()+2*dim, Xi.data());    // mid

        if (in_boundary)
        {
          Xm = (Xm+tmp+2*Xi)/4.;

          //Xm = (Xm+tmp)/2.;

          dX = Xm - Xi;
          VecGetValues(Vec_normal_, dim, edge_dofs_umesh.data()+2*dim, normal.data());
          dX -= normal.dot(dX)*normal;
          Xi += dX;
        }
        else
        {
          Xi = (Xm+tmp)/2.;
        }

        // compute error
        VecGetValues(Vec_x_, dim, edge_dofs_umesh.data()+2*dim, tmp.data());
        error += (tmp-Xi).norm();

        VecSetValues(Vec_x_, dim, edge_dofs_umesh.data()+2*dim, Xi.data(), INSERT_VALUES);
        Assembly(Vec_x_);

      }

    } // end point


    error /= mesh->numNodes();
    //if (error < 1.e-10)
    //{
    //  cout << "STOP, " << smooth_it << " iterations\n";
    //  break;
    //}


  } // end smooth
  Assembly(Vec_x_);

  getVecNormals(&Vec_x_, Vec_normal_);
}

// copy mesh of the data structure to a Petsc Vector
void AppCtx::copyMesh2Vec(Vec & Vec_xmsh)
{
  //bool const u_has_edge_assoc_dof = shape_phi_c->numDofsAssociatedToFacet()+shape_phi_c->numDofsAssociatedToCorner() > 0;

  Vector     Xi(dim);
  VectorXi   node_dofs_mesh(dim);  // indices de onde pegar a velocidade

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    getNodeDofs(&*point, DH_MESH, VAR_M, node_dofs_mesh.data());
    point->getCoord(Xi.data(),dim);
    VecSetValues(Vec_xmsh, dim, node_dofs_mesh.data(), Xi.data(), INSERT_VALUES);
  }
  Assembly(Vec_xmsh);
}

// copy mesh of the data structure to a Petsc Vector
void AppCtx::copyVec2Mesh(Vec const& Vec_xmsh)
{
  //bool const u_has_edge_assoc_dof = shape_phi_c->numDofsAssociatedToFacet()+shape_phi_c->numDofsAssociatedToCorner() > 0;

  Vector     Xi(dim);
  VectorXi   node_dofs_mesh(dim);  // indices de onde pegar a velocidade

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    getNodeDofs(&*point, DH_MESH, VAR_M, node_dofs_mesh.data());
    VecGetValues(Vec_xmsh, dim, node_dofs_mesh.data(), Xi.data());
    point->setCoord(Xi.data(),dim);
  }
}


// copy mesh of the data structure to a Petsc Vector
void AppCtx::swapMeshWithVec(Vec & Vec_xmsh)
{
  //double     *Vec_xmsh_array;
  //bool const u_has_edge_assoc_dof = shape_phi_c->numDofsAssociatedToFacet()+shape_phi_c->numDofsAssociatedToCorner() > 0;

  Vector     tmp(dim); // X mean
  Vector     Xi(dim); // X mean
  VectorXi   node_dofs_umesh(dim);  // indices de onde pegar a velocidade

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    getNodeDofs(&*point, DH_MESH, VAR_M, node_dofs_umesh.data());

    point->getCoord(tmp.data(),dim);
    VecGetValues(Vec_xmsh, dim, node_dofs_umesh.data(), Xi.data());

    point->setCoord(Xi.data(),dim);
    VecSetValues(Vec_xmsh, dim, node_dofs_umesh.data(), tmp.data(), INSERT_VALUES);
  }

  Assembly(Vec_xmsh);
}


PetscErrorCode AppCtx::calcMeshVelocity(Vec const& Vec_x_0, Vec const& Vec_x_1, Vec &Vec_v_mid)
{
  PetscErrorCode ierr;
  ierr = VecCopy(Vec_x_1, Vec_v_mid);      CHKERRQ(ierr);
  ierr = VecAXPY(Vec_v_mid,-1.,Vec_x_0);   CHKERRQ(ierr);
  ierr = VecScale(Vec_v_mid, 1./dt);       CHKERRQ(ierr);
  Assembly(Vec_v_mid);

  PetscFunctionReturn(0);
}

/// @brief move the implicit mesh.
///
/// x1 = Vec_x + dt*(vtheta*Vec_up_1 + (1-vtheta)*Vec_up_0)
///
/// @param[in] Vec_x_1 initial mesh
/// @param[in] Vec_up_0
/// @param[in] Vec_up_1
/// @param[in] Vec_x_new initial guess
/// @param[out] Vec_x_new
/// @warning CHANGE THE VECTOR VEC_NORMAL
PetscErrorCode AppCtx::moveMesh(Vec const& Vec_x_0, Vec const& Vec_up_0, Vec const& Vec_up_1, double const vtheta, double tt, Vec & Vec_x_new)
{
  Vector      Xnew(dim);
  Vector      X0(dim);
  Vector      tmp(dim);
  Vector      U(dim);
  Vector      U0(dim);
  Vector      U1(dim);
  VectorXi    node_dofs_mesh(dim);
  VectorXi    node_dofs_fluid(dim);
  int         tag;
  VectorXi    edge_dofs(dim);
  VectorXi    edge_nodes(3);
  //int         edge_id;
  //Cell const* cell;

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    tag = point->getTag();

    getNodeDofs(&*point, DH_MESH, VAR_M, node_dofs_mesh.data());
    getNodeDofs(&*point, DH_UNKS, VAR_U, node_dofs_fluid.data());

    VecGetValues(Vec_x_0,   dim, node_dofs_mesh.data(), X0.data());

    if (force_mesh_velocity)
    {
      Vector k1 = dt * v_exact(X0,tt,tag);
      Vector k2 = dt * v_exact(X0+0.5*k1,tt+0.5*dt,tag);
      Vector k3 = dt * v_exact(X0+0.5*k2,tt+0.5*dt,tag);
      Vector k4 = dt * v_exact(X0+k3,tt+dt,tag);
      tmp =  (k1 + 2.*(k2+k3) + k4)/6.; // velocity
      //if (!mesh->isVertex(&*point)) // APAGAR TEMP ERASE-ME ... este codigo é só para não entortar a malha
      //{ 
      //  int vtcs[3];
      //  const int m = point->getPosition() - mesh->numVerticesPerCell();
      //  Cell const* cell = mesh->getCell(point->getIncidCell());
      //  if (dim==3)
      //    cell->getCornerVerticesId(m, vtcs);
      //  else
      //    cell->getFacetVerticesId(m, vtcs);
      //  if (mesh->inBoundary(mesh->getNode(vtcs[0])) || mesh->inBoundary(mesh->getNode(vtcs[1])))
      //    tmp = tmp / 2.;
      //  if (mesh->inBoundary(mesh->getNode(vtcs[0])) && mesh->inBoundary(mesh->getNode(vtcs[1])))
      //    tmp = tmp * 0.;
      //}
      tmp = X0 + tmp;
      VecSetValues(Vec_x_new, dim, node_dofs_mesh.data(), tmp.data(), INSERT_VALUES);
    }
    else
    {
      VecGetValues(Vec_up_0,  dim, node_dofs_fluid.data(), U0.data());
      VecGetValues(Vec_up_1,  dim, node_dofs_fluid.data(), U1.data());

      tmp = X0 + dt*(  vtheta*U1 + (1.-vtheta)*U0 );
      VecSetValues(Vec_x_new, dim, node_dofs_mesh.data(), tmp.data(), INSERT_VALUES);
    }
  }

  Assembly(Vec_x_new);

  if (!force_mesh_velocity)
    smoothsMesh(Vec_normal, Vec_x_new);

  PetscFunctionReturn(0);
}

/// @param Vec_x the mesh
/// @param cell_id cell id
/// @return quality number in range ]-inf,1] ... 1 is the best
double AppCtx::getCellQuality(Vec const& Vec_x_, int cell_id) const
{
  MatrixXd  x_coefs(nodes_per_cell, dim);
  VectorXi  mapM_c(dim*nodes_per_cell); // mesh velocity

  dof_handler[DH_MESH].getVariable(VAR_M).getCellDofs(mapM_c.data(), mesh->getCell(cell_id));

  VecGetValues(Vec_x_, mapM_c.size(), mapM_c.data(), x_coefs.data());

  if (dim==2)
  {
    double l_sqr; // l0^2 + l1^2 + l2^2

    l_sqr  = sqr(x_coefs(1,0)-x_coefs(0,0)) + sqr(x_coefs(1,1)-x_coefs(0,1));
    l_sqr += sqr(x_coefs(2,0)-x_coefs(1,0)) + sqr(x_coefs(2,1)-x_coefs(1,1));
    l_sqr += sqr(x_coefs(0,0)-x_coefs(2,0)) + sqr(x_coefs(0,1)-x_coefs(2,1));

    double const f = 6.92820323027551; // 4*sqrt(3)
    double const area = ((x_coefs(0,0)-x_coefs(2,0))*(x_coefs(1,1)-x_coefs(2,1))-(x_coefs(1,0)-x_coefs(2,0))*(x_coefs(0,1)-x_coefs(2,1)))/2;

    return f*area/l_sqr;
  }
  else // if dim==3
  {
    double l_rms_3; // pow( (l0^2 + l1^2 + l2^2 + l3^2 + l4^2 + l5^2 + l6^2)/6 , 1.5)

    l_rms_3  = sqr(x_coefs(0,0)-x_coefs(1,0))   +   sqr(x_coefs(0,1)-x_coefs(1,1))   +   sqr(x_coefs(0,2)-x_coefs(1,2));
    l_rms_3 += sqr(x_coefs(0,0)-x_coefs(2,0))   +   sqr(x_coefs(0,1)-x_coefs(2,1))   +   sqr(x_coefs(0,2)-x_coefs(2,2));
    l_rms_3 += sqr(x_coefs(0,0)-x_coefs(3,0))   +   sqr(x_coefs(0,1)-x_coefs(3,1))   +   sqr(x_coefs(0,2)-x_coefs(3,2));
    l_rms_3 += sqr(x_coefs(1,0)-x_coefs(2,0))   +   sqr(x_coefs(1,1)-x_coefs(2,1))   +   sqr(x_coefs(1,2)-x_coefs(2,2));
    l_rms_3 += sqr(x_coefs(1,0)-x_coefs(3,0))   +   sqr(x_coefs(1,1)-x_coefs(3,1))   +   sqr(x_coefs(1,2)-x_coefs(3,2));
    l_rms_3 += sqr(x_coefs(2,0)-x_coefs(3,0))   +   sqr(x_coefs(2,1)-x_coefs(3,1))   +   sqr(x_coefs(2,2)-x_coefs(3,2));
    l_rms_3 /= 6.;

    l_rms_3 = pow(l_rms_3,1.5);

    double const f = 8.48528137423857; // 6*sqrt(2)

    double x21 = x_coefs(1,0) - x_coefs(0,0);
    double x32 = x_coefs(2,0) - x_coefs(1,0);
    double x43 = x_coefs(3,0) - x_coefs(2,0);
    double y12 = x_coefs(0,1) - x_coefs(1,1);
    double y23 = x_coefs(1,1) - x_coefs(2,1);
    double y34 = x_coefs(2,1) - x_coefs(3,1);
    double z12 = x_coefs(0,2) - x_coefs(1,2);
    double z23 = x_coefs(1,2) - x_coefs(2,2);
    double z34 = x_coefs(2,2) - x_coefs(3,2);
    double volume = x21*(y23*z34 - y34*z23 ) + x32*(y34*z12 - y12*z34 ) + x43*(y12*z23 - y23*z12 );

    return f*volume/l_rms_3;
  }


}

/// @param Vec_x_ the mesh
///
double AppCtx::getCellPatchQuality(Vec const& Vec_x_, int const* cells) const
{
  double quality=99999999, aux;

  for (int i = 0; ; ++i)
  {
    if (cells[i] < 0)
      break;

    aux = getCellQuality(Vec_x_, cells[i]);
    if (aux < quality)
      quality = aux;
  }

  return quality;

}







































