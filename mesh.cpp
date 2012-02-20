#include "common.hpp"






void AppCtx::updateNormals(Vec *Vec_xmsh_0)
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

  if (Vec_xmsh_0==NULL)
    virtual_mesh = false;
  else
    virtual_mesh = true;

  VecSet(Vec_nmsh,0);

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

    dof_handler_mesh.getVariable(0).getFacetDofs(map.data(), &*facet);

    mesh->getFacetNodesId(&*facet, facet_nodes.data());
    if (virtual_mesh)
      VecGetValues(*Vec_xmsh_0, map.size(), map.data(), x_coefs.data());
    else
      mesh->getNodesCoords(facet_nodes.begin(), facet_nodes.end(), x_coefs.data());
    x_coefs_trans = x_coefs.transpose();

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
        normal = F.col(0);
        X = F.col(1);
        // a = a x b
        cross(normal, X);
      }

      VecSetValues(Vec_nmsh, dim, map.data()+k*dim, normal.data(), ADD_VALUES);


    } // nodes

  }
  Assembly(Vec_nmsh);

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

    if (!is_in(tag, solid_tags))// && !is_in(tag, triple_tags))
    {
      VecGetValues(Vec_nmsh, dim, map.data(),normal.data());
      normal.normalize();
      VecSetValues(Vec_nmsh, dim, map.data(),normal.data(), INSERT_VALUES);
      Assembly(Vec_nmsh);
    }
    else
    {
      if (virtual_mesh)
        VecGetValues(*Vec_xmsh_0, dim, map.data(), X.data());
      else
        point->getCoord(X.data());
      normal = -solid_normal(X,current_time,tag);

      VecSetValues(Vec_nmsh, dim, map.data(), normal.data(), INSERT_VALUES);
    }

  }


}

void AppCtx::smoothsMesh(Vec &Vec_xmsh)
{
  //int        nodeid;
  //double     *Vec_xmsh_array;
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
  VectorXi   edge_dofs_umesh(3*dim);
  VectorXi   edge_dofs_fluid((2+u_has_edge_assoc_dof)*dim);
  VectorXi   edge_nodes(3);
  Tensor     R(dim,dim);
  double     error;
  bool       in_boundary;
  //int        id;
  
  
  /* suavização laplaciana */
  for (int smooth_it = 0; smooth_it < 2; ++smooth_it)
  {
    error = 0;
    updateNormals(&Vec_xmsh);

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
        iVs_end = mesh->connectedVtcs(&*point, iVs);

        if (!in_boundary)
        {
          for (int *it = iVs; it != iVs_end ; ++it)
          {
            dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), mesh->getNode(*it));
            // debug
            VecGetValues(Vec_xmsh, dim, vtx_dofs_umesh.data(), tmp.data());
            Xm += tmp;
          }
          Xm /= (iVs_end-iVs);
          dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), &*point);
          
          // compute error
          VecGetValues(Vec_xmsh, dim, vtx_dofs_umesh.data(), tmp.data());
          error += (tmp-Xm).norm();
          
          VecSetValues(Vec_xmsh, dim, vtx_dofs_umesh.data(), Xm.data(), INSERT_VALUES);
        }
        else
        {
          int N=0;
          for (int *it = iVs; it != iVs_end ; ++it)
          {
            Point const* viz_pt = mesh->getNode(*it);
            if (viz_pt->getTag()!=tag && !is_in( viz_pt->getTag(), triple_tags ))
              continue;
            dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), mesh->getNode(*it));
            // debug
            VecGetValues(Vec_xmsh, dim, vtx_dofs_umesh.data(), tmp.data());
            ++N;
            Xm += tmp;
            
          }
          
          dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), &*point);
          VecGetValues(Vec_xmsh, dim, vtx_dofs_umesh.data(), Xi.data());
          
          if (dim==3)
            Xm = (N*Xi + 2*Xm)/(3*N);
            //Xm = Xm/N;
          else
            Xm = (N*Xi + Xm)/(2*N);
            //Xm = Xm/N;
          //

          dX = Xm - Xi;
          VecGetValues(Vec_nmsh, dim, vtx_dofs_umesh.data(), normal.data());
          dX -= normal.dot(dX)*normal;
          Xi += dX;

          // compute error
          VecGetValues(Vec_xmsh, dim, vtx_dofs_umesh.data(), tmp.data());
          error += (tmp-Xi).norm();          
          
          VecSetValues(Vec_xmsh, dim, vtx_dofs_umesh.data(), Xi.data(), INSERT_VALUES);
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

        VecGetValues(Vec_xmsh, dim, edge_dofs_umesh.data(), Xm.data());    // Umsh0
        VecGetValues(Vec_xmsh, dim, edge_dofs_umesh.data()+dim, tmp.data());    // Umsh0
        VecGetValues(Vec_xmsh, dim, edge_dofs_umesh.data()+2*dim, Xi.data());    // Umsh0

        if (in_boundary)
        {
          Xm = (Xm+tmp+2*Xi)/4.;
          
          //Xm = (Xm+tmp)/2.;

          dX = Xm - Xi;
          VecGetValues(Vec_nmsh, dim, edge_dofs_umesh.data()+2*dim, normal.data());
          dX -= normal.dot(dX)*normal;
          Xi += dX;
        }
        else
        {
          Xi = (Xm+tmp)/2.;
        }
        
        // compute error
        VecGetValues(Vec_xmsh, dim, edge_dofs_umesh.data()+2*dim, tmp.data());
        error += (tmp-Xi).norm();             
        
        VecSetValues(Vec_xmsh, dim, edge_dofs_umesh.data()+2*dim, Xi.data(), INSERT_VALUES);
        Assembly(Vec_xmsh);

      }

    } // end point

    
    error /= mesh->numNodes();
    //if (error < 1.e-10)
    //{
    //  cout << "STOP, " << smooth_it << " iterations\n";
    //  break;
    //}
    
    
  } // end smooth
  Assembly(Vec_xmsh);
  
}

// copy mesh of the data structure to a Petsc Vector
void AppCtx::copyMesh2Vec(Vec &Vec_xmsh)
{
  //double     *Vec_xmsh_array;
  bool const u_has_edge_assoc_dof = shape_phi_c->numDofsAssociatedToFacet()+shape_phi_c->numDofsAssociatedToCorner() > 0;

  Vector     Xm(dim); // X mean
  Vector     Xi(dim);
  //Vector     dX(dim);
  //Vector     normal(dim);
  //Vector     tmp(dim), tmp2(dim);
  //Vector     Uf(dim), Ue(dim), Umsh(dim); // Ue := elastic velocity
  //int        iVs[128], *iVs_end;
  VectorXi   vtx_dofs_umesh(dim);  // indices de onde pegar a velocidade
  //VectorXi   vtx_dofs_fluid(dim); // indices de onde pegar a velocidade
  VectorXi   edge_dofs_umesh(3*dim);
  VectorXi   edge_dofs_fluid((2+u_has_edge_assoc_dof)*dim);
  VectorXi   edge_nodes(3);
  //Tensor     R(dim,dim);
  //int        id;

  //VecGetArray(Vec_xmsh, &Vec_xmsh_array);

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    if (!mesh->isVertex(&*point))
        continue;
    dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), &*point);
    //dof_handler_vars.getVariable(0).getVertexDofs(vtx_dofs_fluid.data(), &*point);
    //VecGetValues(Vec_up_1, dim, vtx_dofs_fluid.data(), Uf.data());
    //id = mesh->getPointId(&*point);
    //getRotationMatrix(R,&id,1);
    point->getCoord(Xi.data());
    //Xi = Xi + dt*R.transpose()*Uf;
    VecSetValues(Vec_xmsh, dim, vtx_dofs_umesh.data(), Xi.data(), INSERT_VALUES);
  }
  if (dim==2 && u_has_edge_assoc_dof)
  {
    facet_iterator edge = mesh->facetBegin();
    facet_iterator edge_end = mesh->facetEnd();
    for (; edge != edge_end; ++edge)
    {
      mesh->getFacetNodesId(&*edge, edge_nodes.data());
      dof_handler_mesh.getVariable(0).getFacetAssociatedDofs(edge_dofs_umesh.data(), &*edge);
     // dof_handler_vars.getVariable(0).getFacetAssociatedDofs(edge_dofs_fluid.data(), &*edge);
      //VecGetValues(Vec_up_1, dim, edge_dofs_fluid.data(), Uf.data());
      mesh->getNode(edge_nodes(2))->getCoord(Xi.data());
      //getRotationMatrix(R,&edge_nodes(2),1);
      //Xi = Xi + dt*R.transpose()*Uf;
      VecSetValues(Vec_xmsh, dim, edge_dofs_umesh.data(), Xi.data(), INSERT_VALUES);
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
      //dof_handler_vars.getVariable(0).getCornerAssociatedDofs(edge_dofs_fluid.data(), &*edge);
      //VecGetValues(Vec_up_1, dim, edge_dofs_fluid.data(), Uf.data());
      mesh->getNode(edge_nodes(2))->getCoord(Xi.data());
      //getRotationMatrix(R,&edge_nodes(2),1);
      //Xi = Xi + dt*R.transpose()*Uf;
      VecSetValues(Vec_xmsh, dim, edge_dofs_umesh.data(), Xi.data(), INSERT_VALUES);
    }
  }

  Assembly(Vec_xmsh);  
}


// copy mesh of the data structure to a Petsc Vector
void AppCtx::copyVec2Mesh(Vec const& Vec_xmsh)
{
  bool const u_has_edge_assoc_dof = shape_phi_c->numDofsAssociatedToFacet()+shape_phi_c->numDofsAssociatedToCorner() > 0;

  Vector     Xm(dim); // X mean
  Vector     Xi(dim);
  VectorXi   vtx_dofs_umesh(dim);  // indices de onde pegar a velocidade
  VectorXi   edge_dofs_umesh(3*dim);
  VectorXi   edge_dofs_fluid((2+u_has_edge_assoc_dof)*dim);
  VectorXi   edge_nodes(3);


  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    if (!mesh->isVertex(&*point))
        continue;
    dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), &*point);

    VecGetValues(Vec_xmsh, dim, vtx_dofs_umesh.data(), Xi.data());
    point->setCoord(Xi.data());
  }
  if (dim==2 && u_has_edge_assoc_dof)
  {
    facet_iterator edge = mesh->facetBegin();
    facet_iterator edge_end = mesh->facetEnd();
    for (; edge != edge_end; ++edge)
    {
      mesh->getFacetNodesId(&*edge, edge_nodes.data());
      dof_handler_mesh.getVariable(0).getFacetAssociatedDofs(edge_dofs_umesh.data(), &*edge);

      VecGetValues(Vec_xmsh, dim, edge_dofs_umesh.data(), Xi.data());
      mesh->getNode(edge_nodes(2))->setCoord(Xi.data());
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

      VecGetValues(Vec_xmsh, dim, edge_dofs_umesh.data(), Xi.data());
      mesh->getNode(edge_nodes(2))->setCoord(Xi.data());
    }
  }

}


// copy mesh of the data structure to a Petsc Vector
void AppCtx::swapMeshWithVec(Vec & Vec_xmsh)
{
  //double     *Vec_xmsh_array;
  bool const u_has_edge_assoc_dof = shape_phi_c->numDofsAssociatedToFacet()+shape_phi_c->numDofsAssociatedToCorner() > 0;

  Vector     tmp(dim); // X mean
  Vector     Xi(dim);
  VectorXi   vtx_dofs_umesh(dim);  // indices de onde pegar a velocidade
  VectorXi   edge_dofs_umesh(3*dim);
  VectorXi   edge_dofs_fluid((2+u_has_edge_assoc_dof)*dim);
  VectorXi   edge_nodes(3);

  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
  for (; point != point_end; ++point)
  {
    if (!mesh->isVertex(&*point))
        continue;
    dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), &*point);
    
    point->getCoord(tmp.data());
    VecGetValues(Vec_xmsh, dim, vtx_dofs_umesh.data(), Xi.data());
    
    point->setCoord(Xi.data());
    VecSetValues(Vec_xmsh, dim, vtx_dofs_umesh.data(), tmp.data(), INSERT_VALUES);
  }
  if (dim==2 && u_has_edge_assoc_dof)
  {
    facet_iterator edge = mesh->facetBegin();
    facet_iterator edge_end = mesh->facetEnd();
    for (; edge != edge_end; ++edge)
    {
      mesh->getFacetNodesId(&*edge, edge_nodes.data());
      dof_handler_mesh.getVariable(0).getFacetAssociatedDofs(edge_dofs_umesh.data(), &*edge);
      
      mesh->getNode(edge_nodes(2))->getCoord(tmp.data());
      VecGetValues(Vec_xmsh, dim, edge_dofs_umesh.data(), Xi.data());
      
      mesh->getNode(edge_nodes(2))->setCoord(Xi.data());
      VecSetValues(Vec_xmsh, dim, edge_dofs_umesh.data(), tmp.data(), INSERT_VALUES);
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
      
      mesh->getNode(edge_nodes(2))->getCoord(tmp.data());
      VecGetValues(Vec_xmsh, dim, edge_dofs_umesh.data(), Xi.data());
      
      mesh->getNode(edge_nodes(2))->setCoord(Xi.data());
      VecSetValues(Vec_xmsh, dim, edge_dofs_umesh.data(), tmp.data(), INSERT_VALUES);
    }
  }

  Assembly(Vec_xmsh);  
}



/// @param[in] Vec_up vector with fluid velocity and pressure
/// @param[in] Vec_xmsh the mesh
/// @param[out] Vec_vmsh
PetscErrorCode AppCtx::calcMeshVelocity(Vec const& Vec_up, Vec const& Vec_xmsh, Vec & Vec_vmsh, double const current_time)
{
  int        nodeid;
  //double     *Vec_xmsh_array;
  bool const u_has_edge_assoc_dof = shape_phi_c->numDofsAssociatedToFacet()+shape_phi_c->numDofsAssociatedToCorner() > 0;

  Vector     Xm(dim); // X mean
  Vector     Xi(dim);
  Vector     dX(dim);
  Vector     normal(dim);
  Vector     tmp(dim), tmp2(dim);
  Vector     Uf(dim), Ue(dim), Umsh(dim); // Ue := elastic velocity
  int        tag;
  //int        iVs[128], *iVs_end;
  VectorXi   vtx_dofs_umesh(dim);  // indices de onde pegar a velocidade
  VectorXi   vtx_dofs_fluid(dim); // indices de onde pegar a velocidade
  VectorXi   edge_dofs_umesh(3*dim);
  VectorXi   edge_dofs_fluid((2+u_has_edge_assoc_dof)*dim);
  VectorXi   edge_nodes(3);
  Tensor     R(dim,dim);
  bool       in_boundary;
  //int        id;

  // compute elastic velocity
  VecCopy(Vec_xmsh, Vec_vmsh);
  smoothsMesh(Vec_vmsh);
  VecAXPY(Vec_vmsh,-1.,Vec_xmsh);
  VecScale(Vec_vmsh, 1./dt);

  /* calculando a velocidade da malha*/
  point_iterator point = mesh->pointBegin();
  point_iterator point_end = mesh->pointEnd();
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
      VecGetValues(Vec_up, dim, vtx_dofs_fluid.data(), Uf.data());
      // pega a velocidade elástica
      VecGetValues(Vec_vmsh, dim, vtx_dofs_umesh.data(), Ue.data());
      //// pega a coordenada atual
      //point->getCoord(Xi.data());

      //if (in_boundary && force_bound_mesh_vel)
      if (force_bound_mesh_vel)
      {
        VecGetValues(Vec_xmsh, dim, vtx_dofs_umesh.data(), Xi.data());
        
        
        point->getCoord(Xm.data());
        double dtt = this->dt;
        double k1=current_time;
        double k2=current_time+dtt/4.;
        double k3=current_time+dtt/2.;
        double k4=current_time+3.*dtt/4.;
        double k5=current_time+dtt;
        Umsh  =  7.*v_exact(Xm,k1,tag);
        Umsh += 32.*v_exact(Xm,k2,tag);
        Umsh += 12.*v_exact(Xm,k3,tag);
        Umsh += 32.*v_exact(Xm,k4,tag);
        Umsh +=  7.*v_exact(Xm,k5,tag);
        Umsh /= 90.;

        //////////////////////Umsh = v_exact(Xi,current_time,tag);
        VecSetValues(Vec_vmsh, dim, vtx_dofs_umesh.data(), Umsh.data(),INSERT_VALUES);
        continue;
      }
      else
      if ((!boundary_smoothing && in_boundary) || is_in(tag,triple_tags))
      {
        // se esta na sup., a vel. da malha é igual a do fluido
        Umsh = Uf;
        VecSetValues(Vec_vmsh, dim, vtx_dofs_umesh.data(), Umsh.data(),INSERT_VALUES);
        continue;
      }

      getRotationMatrix(R,&nodeid,1);

      if (in_boundary)
      {
        VecGetValues(Vec_nmsh, dim, vtx_dofs_umesh.data(), normal.data());
        
        Ue = Ue - normal.dot(Ue)*normal;
        Uf = R.transpose()*Uf;
        Uf = normal.dot(Uf)*normal;
        
        Umsh = R*(Uf + Ue);
        //Umsh = Uf;
      }
      else
        // se não está no contorno, não precisa rotacionar
        Umsh = beta1*Uf + beta2*Ue;

      VecSetValues(Vec_vmsh, dim, vtx_dofs_umesh.data(), Umsh.data(), INSERT_VALUES);
    }
    else
    //
    // mid node
    //
    {
      // pega os graus de liberdade
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

      // pega a velocidade do fluido
      VecGetValues(Vec_up,   dim, edge_dofs_fluid.data()+2*dim, Uf.data());
      // pega a velocidade elástica
      VecGetValues(Vec_vmsh, dim, edge_dofs_umesh.data()+2*dim, Ue.data());
      // pega a coordenada atual

      //if (in_boundary && force_bound_mesh_vel)
      if (force_bound_mesh_vel)
      {
        VecGetValues(Vec_xmsh, dim, edge_dofs_umesh.data()+2*dim, Xi.data());
        
        Umsh = v_exact(Xi,current_time,tag);
        VecSetValues(Vec_vmsh, dim, edge_dofs_umesh.data()+2*dim, Umsh.data(),INSERT_VALUES);
        continue;
      }
      else      
      if ((!boundary_smoothing && in_boundary) || is_in(tag,triple_tags))
      {
        // se esta na sup., a vel. da malha é igual a do fluido
        Umsh = Uf;
        VecSetValues(Vec_vmsh, dim, edge_dofs_fluid.data()+2*dim, Umsh.data(),INSERT_VALUES);
        continue;
      }

      // se está na interface
      if (in_boundary )
      {
        getRotationMatrix(R,&nodeid,1);
        
        VecGetValues(Vec_nmsh, dim, edge_dofs_umesh.data()+2*dim, normal.data());
        Ue = Ue - normal.dot(Ue)*normal;
        Uf = R.transpose()*Uf;
        Uf = normal.dot(Uf)*normal;
        Umsh = R*(Uf + Ue);

        VecSetValues(Vec_vmsh, dim, edge_dofs_umesh.data()+2*dim, Umsh.data(),INSERT_VALUES);
        
      }
      else // se está no interior do domínio, arrasta o nó para o centro da aresta
      {
        Umsh = R*Ue;

        VecSetValues(Vec_vmsh, dim, edge_dofs_umesh.data()+2*dim, Umsh.data(), INSERT_VALUES);
      }
    }


  } // end point loop


  Assembly(Vec_vmsh);
  //VecRestoreArray(Vec_xmsh, &Vec_xmsh_array);

  PetscFunctionReturn(0);
}

/// @brief move the implicit mesh.
///
/// x1 = Vec_xmsh_0 + dt*(vtheta*Vec_vmsh_med + (1-vtheta)*Vec_vmsh_0)
///
/// @param[in] Vec_xmsh_0 initial mesh
/// @param[in] Vec_vmsh_0
/// @param[in] Vec_vmsh_med
/// @param[out] implicit mesh
///
PetscErrorCode AppCtx::moveMesh(Vec const& Vec_xmsh_0, Vec const& Vec_vmsh_0, Vec const& Vec_vmsh_1, double const vtheta)
{
  Vector      Xi(dim);
  Vector      X0(dim);
  Vector      tmp(dim);
  Vector      U(dim);
  Vector      Umsh0(dim);
  Vector      Umsh1(dim);
  VectorXi    vtx_dofs_umesh(3); // indices de onde pegar a velocidade
  //int       vtx_dofs_fluid[3]; // indices de onde pegar a velocidade
  int         tag;
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
    tag = point->getTag();
    if (mesh->isVertex(&*point))
    {
      dof_handler_mesh.getVariable(0).getVertexDofs(vtx_dofs_umesh.data(), &*point);

      VecGetValues(Vec_xmsh_0, dim, vtx_dofs_umesh.data(), X0.data());
      VecGetValues(Vec_vmsh_0, dim, vtx_dofs_umesh.data(), Umsh0.data());
      VecGetValues(Vec_vmsh_1, dim, vtx_dofs_umesh.data(), Umsh1.data());

      nodeid = mesh->getPointId(&*point);
      getRotationMatrix(R, &nodeid, 1);

      //point->getCoord(Xi.data());
      ///////////tmp = X0 + dt*R.transpose()*(  vtheta*Umsh1 + (1.-vtheta)*Umsh0 );
      double dtt = this->dt;
      double k1=current_time;
      double k2=current_time+dtt/4.;
      double k3=current_time+dtt/2.;
      double k4=current_time+3.*dtt/4.;
      double k5=current_time+dtt;
      tmp  =  7.*v_exact(X0,k1,tag);
      tmp += 32.*v_exact(X0,k2,tag);
      tmp += 12.*v_exact(X0,k3,tag);
      tmp += 32.*v_exact(X0,k4,tag);
      tmp +=  7.*v_exact(X0,k5,tag);
      tmp *= dtt/90.;
      tmp += X0;
      point->setCoord(tmp.data() );
    }
  }

  // higher order nodes
  bool const mesh_has_edge_node = mesh->numNodesPerCell()-mesh->numVerticesPerCell() > 0;
  if (mesh_has_edge_node)
  {
    //Point * point;

    if (dim==3) // MESMA COISA, trocando corner <-> facet
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

        VecGetValues(Vec_xmsh_0, dim, edge_dofs.data(), X0.data());
        VecGetValues(Vec_vmsh_0, dim, edge_dofs.data(), Umsh0.data());
        VecGetValues(Vec_vmsh_1, dim, edge_dofs.data(), Umsh1.data());

        nodeid = mesh->getPointId(&*point);
        getRotationMatrix(R, &nodeid, 1);

        //point->getCoord(Xi.data());
        tmp = X0 + dt*R.transpose()*(  vtheta*Umsh1 + (1.-vtheta)*Umsh0 );
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

        VecGetValues(Vec_xmsh_0, dim, edge_dofs.data(), X0.data());
        VecGetValues(Vec_vmsh_0, dim, edge_dofs.data(), Umsh0.data());
        VecGetValues(Vec_vmsh_1, dim, edge_dofs.data(), Umsh1.data());

        nodeid = mesh->getPointId(&*point);
        getRotationMatrix(R, &nodeid, 1);

        //point->getCoord(Xi.data());
        tmp = X0 + dt*R.transpose()*(  vtheta*Umsh1 + (1.-vtheta)*Umsh0 );
        point->setCoord(tmp.data() );

      } // end loop
    } // end dim==2


  }

  updateNormals(NULL);

  PetscFunctionReturn(0);
}


















































