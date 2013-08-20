#include "common.hpp"
#include <tr1/tuple>

using std::tr1::tuple;
using std::tr1::make_tuple;
using std::tr1::get;

extern PetscErrorCode FormJacobian(SNES snes,Vec Vec_up_1,Mat *Mat_Jac, Mat *prejac, MatStructure *flag, void *ptr);
extern PetscErrorCode FormFunction(SNES snes, Vec Vec_up_1, Vec Vec_fun, void *ptr);
extern PetscErrorCode CheckSnesConvergence(SNES snes, PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx);

extern PetscErrorCode FormJacobian_mesh(SNES snes,Vec Vec_up_1,Mat *Mat_Jac, Mat *prejac, MatStructure *flag, void *ptr);
extern PetscErrorCode FormFunction_mesh(SNES snes, Vec Vec_up_1, Vec Vec_fun, void *ptr);

void checkConsistencyTri(Mesh *mesh)
{
  cell_iterator cell = mesh->cellBegin();
  cell_iterator cell_end = mesh->cellEnd();
  
  int const nfpc = mesh->numFacetsPerCell();
  //  int const nnpc = mesh->numNodesPerCell();
  int const nnpf = mesh->numNodesPerFacet();
  
  int f_nds[nnpf];
  
  //  Point *p;
  Facet *f;
  Cell  *c;
  
  for (; cell != cell_end; ++cell)
  {
    int myid = mesh->getCellId(&*cell);
    for (int i = 0; i < nfpc; ++i)
    {
      if (cell->getIncidCell(i) >= 0)
      {
        // verifica o vizinho
        c = mesh->getCellPtr(cell->getIncidCell(i));
        int pos = cell->getIncidCellPos(i);
        if(myid != c->getIncidCell(pos))
        {
          cout  << "myid=" << myid<<"; c->getIncidCell(pos)="<<c->getIncidCell(pos)<<"; i="<<i<<"; pos="<<pos ;
          throw;
        };
        
        // verifica a face
        f = mesh->getFacetPtr(cell->getFacetId(i));
        int icf = f->getIncidCell();
        if(!(icf==myid || icf==cell->getIncidCell(i)))
        {
          cout <<"myid="<<myid<<"; icf="<<icf<<"; i="<<i<<"; cell->getIncidCell(i)="<<cell->getIncidCell(i)<<"\n";
          throw;
        }
        if (icf==myid)
        {
          if(f->getPosition() != i)
          {
            cout << "myid=" << myid<<"; f->getPosition()="<<f->getPosition()<<"; i="<<i<<"\n";
            throw;
          }
        }
        else
        {
          if(f->getPosition() != pos)
          {
            cout << "myid=" << myid<<"; f->getPosition()="<<f->getPosition()<<"; pos="<<pos<<"; i="<<i<<"\n";
            throw;
          }
        }
      }
      else // bordo
      {
        // verifica a face
        f = mesh->getFacetPtr(cell->getFacetId(i));
        int icf = f->getIncidCell();
        // só pode ser o myid, pq do outro lado não tem ngm
        if(icf!=myid)
        {
          cout << "icf = " << icf << ", myid = " << myid << endl;
          throw;
        }
        
        // verifica se os nós da face estão no contorno
        mesh->getFacetNodesId(f, f_nds);
        for (int j = 0; j < nnpf; ++j)
        {
          if(!mesh->inBoundary(mesh->getNodePtr(f_nds[j])))
          {
            cout << "node="<<f_nds[j]<<"; icell="<<mesh->getNodePtr(f_nds[j])->getIncidCell()
                << "; pos="<< mesh->getNodePtr(f_nds[j])->getPosition();
            throw;
          }
        }
        
      }
    }
    
  }
  
  std::vector<int> ics;
  std::vector<int> cc_ids;
  std::vector<int>::iterator it;
  for (point_iterator point = mesh->pointBegin(); point != mesh->pointEnd(); ++point)
  {
    point->getAllIncidences(ics);
    
    cc_ids.clear();

    int myid = point.index();
    
    for (int i = 0; i < (int)ics.size()/2; ++i)
    {
      int ic = ics.at(2*i);
      int pos = ics.at(2*i+1);
      Cell *c = mesh->getCellPtr(ic);
      if(c==NULL)
      {
        cout << "c = NULL" << endl;
        throw;
      }
      if(myid != c->getNodeId(pos))
      {
        cout << "ic = " << ic << "; pos = " << pos;
        throw;
      }
      cc_ids.push_back(c->getConnectedComponentId());
    }
    
    std::sort(cc_ids.begin(), cc_ids.end());
    it = unique (cc_ids.begin(), cc_ids.end());
    
    // checks if all incident cells have distinct connected component id
    if( std::distance(cc_ids.begin(), it) != (int)cc_ids.size())
    {
      cout << "MERDA DE std::distance\n";
      throw;
    }
    
  }
  
  
}

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
  int               tag_0, tag_1, tag_2;
  bool              contrib_0, contrib_1, contrib_2;
  //int               tag_other;
  bool              virtual_mesh;
  bool              is_surface;
  bool              is_solid;
  bool              is_cl;
  int               sign_;
  int               pts_id[15];


  if (Vec_x_1==NULL)
    virtual_mesh = false;
  else
    virtual_mesh = true;

  Assembly(Vec_normal_);
  VecSet(Vec_normal_,0);
  Assembly(Vec_normal_);

  // LOOP NAS FACES DO CONTORNO
  facet_iterator facet = mesh->facetBegin();
  facet_iterator facet_end = mesh->facetEnd();
  for (; facet != facet_end; ++facet)
  {
    tag = facet->getTag();

    is_surface = is_in(tag, interface_tags);
    is_solid   = is_in(tag, solid_tags);
    is_cl      = is_in(tag, triple_tags);

    if ( !(is_surface || is_solid || is_cl || mesh->inBoundary(&*facet)) )
      continue;

    contrib_0 = true;
    contrib_1 = true;
    contrib_2 = true;
    // the solid normal doesn't contribute to the triple normal
    if (is_solid)
    {
      mesh->getFacetNodesId(&*facet, pts_id);

      tag_0 = mesh->getNodePtr(pts_id[0])->getTag();
      contrib_0 = !is_in(tag_0, triple_tags);

      tag_1 = mesh->getNodePtr(pts_id[1])->getTag();
      contrib_1 = !is_in(tag_1, triple_tags);

      if (dim==3)
      {
        tag_2 = mesh->getNodePtr(pts_id[2])->getTag();
        contrib_2 = !is_in(tag_2, triple_tags);
      }
    }

    dof_handler[DH_MESH].getVariable(VAR_M).getFacetDofs(map.data(), &*facet);

    mesh->getFacetNodesId(&*facet, facet_nodes.data());
    if (virtual_mesh)
      VecGetValues(*Vec_x_1, map.size(), map.data(), x_coefs.data());
    else
      mesh->getNodesCoords(facet_nodes.begin(), facet_nodes.end(), x_coefs.data());
    x_coefs_trans = x_coefs.transpose();

    // fix orientation in the case where the gas phase isn't passive
    {
      sign_ = 1;
      cell_handler cc = mesh->getCell(facet->getIncidCell());
      cell_handler oc = mesh->getCell(cc->getIncidCell(facet->getPosition()));
      if ( oc.isValid()  )
        if ( oc->getTag() > cc->getTag() )
          sign_ = -1;
    }

    if (false) // método alternativo que serve para todos os métodos (normal consistente)
    {
      // find the normal
      for (int k = 0; k < nodes_per_facet; ++k)
      {
        //tag_other = mesh->getNodePtr(facet_nodes(k))->getTag();

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
      normal *= sign_;

      // 0
      if (contrib_0)
      {
        X = normal;
        X /= a.dot(a) * c.dot(c);
        VecSetValues(Vec_normal_, dim, map.data()+0*dim, X.data(), ADD_VALUES);
      }

      // 1
      if (contrib_1)
      {
        X = normal;
        X /= a.dot(a) * b.dot(b);
        VecSetValues(Vec_normal_, dim, map.data()+1*dim, X.data(), ADD_VALUES);
      }

      // 2
      if (contrib_2)
      {
        X = normal;
        X /= b.dot(b) * c.dot(c);
        VecSetValues(Vec_normal_, dim, map.data()+2*dim, X.data(), ADD_VALUES);
      }

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
      if (contrib_0)
      {
        cross(normal,x3-x0,x5-x0);
        normal /= e03*e05;
        normal *= sign_;
        VecSetValues(Vec_normal_, dim, map.data()+0*dim, normal.data(), ADD_VALUES);
      }

      // node 1
      if (contrib_1)
      {
        cross(normal,x4-x1,x3-x1);
        normal /= e14*e13;
        normal *= sign_;
        VecSetValues(Vec_normal_, dim, map.data()+1*dim, normal.data(), ADD_VALUES);
      }

      // node 2
      if (contrib_2)
      {
        cross(normal,x5-x2,x4-x2);
        normal /= e25*e24;
        normal *= sign_;
        VecSetValues(Vec_normal_, dim, map.data()+2*dim, normal.data(), ADD_VALUES);
      }

      // node 3
      if (contrib_0 && contrib_1)
      {
        normal = cross(x1-x3,x4-x3)/(e13*e34) + cross(x4-x3,x5-x3)/(e34*e35) + cross(x5-x3,x0-x3)/(e35*e03);
        normal *= sign_;
        VecSetValues(Vec_normal_, dim, map.data()+3*dim, normal.data(), ADD_VALUES);
      }

      // node 4
      if (contrib_1 && contrib_2)
      {
        normal = cross(x2-x4,x5-x4)/(e24*e45) + cross(x5-x4,x3-x4)/(e45*e34) + cross(x3-x4,x1-x4)/(e34*e14);
        normal *= sign_;
        VecSetValues(Vec_normal_, dim, map.data()+4*dim, normal.data(), ADD_VALUES);
      }

      // node 5
      if (contrib_2 && contrib_0)
      {
        normal = cross(x0-x5,x3-x5)/(e05*e35) + cross(x3-x5,x4-x5)/(e35*e45) + cross(x4-x5,x2-x5)/(e45*e25);
        normal *= sign_;
        VecSetValues(Vec_normal_, dim, map.data()+5*dim, normal.data(), ADD_VALUES);
      }
    }
    else if(mesh_cell_type == TRIANGLE3)
    {
      normal(0) = x_coefs_trans(1,1)-x_coefs_trans(1,0);
      normal(1) = x_coefs_trans(0,0)-x_coefs_trans(0,1);

      normal /= (x_coefs_trans.col(0)-x_coefs_trans.col(1)).squaredNorm();
      normal *= sign_;
      if (contrib_0)
        VecSetValues(Vec_normal_, dim, map.data()+0*dim, normal.data(), ADD_VALUES);
      if (contrib_1)
        VecSetValues(Vec_normal_, dim, map.data()+1*dim, normal.data(), ADD_VALUES);
      //if (is_surface)
      //  cout << "normal = " << normal(0) << " " << normal(1) << endl;
    }
    else if(mesh_cell_type == TRIANGLE6) // dividi a aresta em duas partes
    {
      normal(0) = x_coefs_trans(1,2)-x_coefs_trans(1,0);
      normal(1) = x_coefs_trans(0,0)-x_coefs_trans(0,2);
      normal /= (x_coefs_trans.col(0)-x_coefs_trans.col(2)).squaredNorm();
      normal *= sign_;
      if (contrib_0)
        VecSetValues(Vec_normal_, dim, map.data()+0*dim, normal.data(), ADD_VALUES);
      VecSetValues(Vec_normal_, dim, map.data()+2*dim, normal.data(), ADD_VALUES);

      normal(0) = x_coefs_trans(1,1)-x_coefs_trans(1,2);
      normal(1) = x_coefs_trans(0,2)-x_coefs_trans(0,1);
      normal /= (x_coefs_trans.col(1)-x_coefs_trans.col(2)).squaredNorm();
      normal *= sign_;
      if (contrib_1)
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

    is_surface = is_in(tag, interface_tags);
    is_solid   = is_in(tag, solid_tags);
    is_cl      = is_in(tag, triple_tags);

    if ( !(is_surface || is_solid || is_cl || mesh->inBoundary(&*point)) )
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
      Assembly(Vec_normal_);
    }

  }
  Assembly(Vec_normal_);

}

/// @param[out] Vec_normal_
/// @param[out] Vec_x_
void AppCtx::smoothsMesh(Vec & Vec_normal_, Vec &Vec_x_)
{
  //int        nodeid;
  //double     *Vec_x__array;
  bool const u_has_edge_assoc_dof = shape_phi_c->numDofsAssociatedToFacet()+shape_phi_c->numDofsAssociatedToCorner() > 0;

  Vector      Xm(dim); // X mean
  Vector      X0(dim);
  Vector      Xi(dim);
  Vector      dX(dim);
  Vector      normal(dim);
  Vector      tmp(dim), tmp2(dim);
  Vector      Uf(dim), Ue(dim), Umsh(dim); // Ue := elastic velocity
  int         tag;
  int         viz_tag;
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
  bool        is_surface;
  bool        is_solid;

  int const n_smooths = 10;

  /* suavização laplaciana */
  for (int smooth_it = 0; smooth_it < n_smooths; ++smooth_it)
  {
    error = 0;
    getVecNormals(&Vec_x_, Vec_normal_);

    // VERTICES
    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    for (; point != point_end; ++point)
    {
      tag = point->getTag();

      is_surface = is_in(tag, interface_tags);
      is_solid   = is_in(tag, solid_tags);

      in_boundary =  is_surface || is_solid;

      // pula o caso em que o ponto está no contoro mas não tem boundary smoothing
      if (!boundary_smoothing && in_boundary)
        continue;

      if (is_in(tag,triple_tags) || is_in(tag,feature_tags))
          continue;

      if (is_in(tag,dirichlet_tags) || is_in(tag,neumann_tags) || is_in(tag,periodic_tags))
          continue;

      if (mesh->isVertex(&*point))
      {
        //if (  is_in(tag,interface_tags) || is_in(tag,triple_tags) || is_in(tag,solid_tags) ||
        //    is_in(tag,dirichlet_tags) || is_in(tag,neumann_tags)  )

        dof_handler[DH_MESH].getVariable(VAR_M).getVertexDofs(vtx_dofs_umesh.data(), &*point);
        VecGetValues(Vec_x_, dim, vtx_dofs_umesh.data(), X0.data()); // old coord

        Xm = Vector::Zero(dim);
        //iVs_end = mesh->connectedVtcs(&*point, iVs);
        iVs_end = mesh->connectedVtcs(&*point, iVs, iCs, viCs);

        if (!in_boundary)
        {
          int N=0;
          Point const* viz_pt;
          for (int *it = iVs; it != iVs_end ; ++it)
          {
            viz_pt = mesh->getNodePtr(*it);
            viz_tag = viz_pt->getTag();
            dof_handler[DH_MESH].getVariable(VAR_M).getVertexDofs(vtx_dofs_umesh.data(), viz_pt);
            VecGetValues(Vec_x_, dim, vtx_dofs_umesh.data(), tmp.data());
            ++N;
            Xm += tmp;
            //if (isFixedPoint(viz_tag))
            //{
            //  Xm += 5*X0;
            //  N += 5;
            //}
          }
          Xm /= N;

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
            Point const* viz_pt = mesh->getNodePtr(*it);
            //if (viz_pt->getTag()!=tag && !is_in( viz_pt->getTag(), triple_tags ))
            viz_tag = viz_pt->getTag();
            if (viz_tag!=tag && (is_in(viz_tag,interface_tags) || is_in(viz_tag,solid_tags)) )
              continue;
            if (viz_tag!=tag && !is_in(viz_tag,triple_tags) && !isFixedPoint(viz_tag) && !is_in(viz_tag,feature_tags))
              continue;
            dof_handler[DH_MESH].getVariable(VAR_M).getVertexDofs(vtx_dofs_umesh.data(), mesh->getNodePtr(*it));
            // debug
            VecGetValues(Vec_x_, dim, vtx_dofs_umesh.data(), tmp.data());
            ++N;
            Xm += tmp;
            //if (isFixedPoint(viz_tag))
            //{
            //  Xm += 3*X0;
            //  N += 3;
            //}
          }

          dof_handler[DH_MESH].getVariable(VAR_M).getVertexDofs(vtx_dofs_umesh.data(), &*point);
          VecGetValues(Vec_x_, dim, vtx_dofs_umesh.data(), Xi.data());

          if (dim==3)
            Xm = (N*Xi + 2*Xm)/(3*N);
            //Xm = Xm/N;
          else
            //Xm = (N*Xi + Xm)/(2*N);
            Xm = Xm/N;
          //

          dX = Xm - Xi;
          VecGetValues(Vec_normal_, dim, vtx_dofs_umesh.data(), normal.data());
          dX -= normal.dot(dX)*normal;
          Xi += dX;

          // compute error
          VecGetValues(Vec_x_, dim, vtx_dofs_umesh.data(), tmp.data());
          error += (tmp-Xi).norm();

          //old_quality = getCellPatchQuality(Vec_x_, iCs);
          VecSetValues(Vec_x_, dim, vtx_dofs_umesh.data(), Xi.data(), INSERT_VALUES);
          //new_quality = getCellPatchQuality(Vec_x_, iCs);
          //// se a qualidade piorou, volta no que estava antes
          //if (new_quality < old_quality)
          //{
          //  VecSetValues(Vec_x_, dim, vtx_dofs_umesh.data(), tmp.data(), INSERT_VALUES);
          //}
        }

      }

    } // end point


    // MID NODES
    point = mesh->pointBegin();
    point_end = mesh->pointEnd();
    if (u_has_edge_assoc_dof)
    for (; point != point_end; ++point)
    {
      tag = point->getTag();

      is_surface = is_in(tag, interface_tags);
      is_solid   = is_in(tag, solid_tags);

      in_boundary = mesh->inBoundary(&*point) || is_surface || is_solid;

      // pula o caso em que o ponto está no contoro mas não tem boundary smoothing
      if (!boundary_smoothing && in_boundary)
        continue;

      if (!mesh->isVertex(&*point))
      {

        const int m = point->getPosition() - mesh->numVerticesPerCell();
        Cell const* icell = mesh->getCellPtr(point->getIncidCell());
        if (dim==3)
        {
          Corner *edge = mesh->getCornerPtr(icell->getCornerId(m));
          dof_handler[DH_MESH].getVariable(VAR_M).getCornerDofs(edge_dofs_umesh.data(), &*edge);
          dof_handler[DH_UNKS].getVariable(VAR_U).getCornerDofs(edge_dofs_fluid.data(), &*edge);
          mesh->getCornerNodesId(&*edge, edge_nodes.data());
        }
        else // dim=2
        {
          Facet *edge = mesh->getFacetPtr(icell->getFacetId(m));
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


PetscErrorCode AppCtx::meshAdapt()
{
  //PetscFunctionReturn(0);
  PetscErrorCode      ierr;

  // only for 2d ... TODO for 3d too
  if (dim != 2)
    PetscFunctionReturn(0);
  // only for linear elements
  if (mesh->numNodesPerCell() > mesh->numVerticesPerCell())
    PetscFunctionReturn(0);

  const Real TOL = 0.6;

  typedef tuple<int, int, int> EdgeVtcs; // get<0> = mid node, get<1> = top node, get<2> = bot node


  std::list<EdgeVtcs> adde_vtcs;

  CellElement *edge;
  Real h;
  Real expected_h;
  VectorXi edge_nodes(3); // 3 nodes at most
  Vector Xa(dim), Xb(dim);
  int tag_a;
  int tag_b;
  int tag_e;

  bool mesh_was_changed = false;

  //printf("ENTRANDO NO LOOP is_splitting!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  // first make collapses, then make splittings.
  // NOTE: mark added points
  for (int is_splitting = 0; is_splitting < 2 ; ++is_splitting)
  {
    int const n_edges_total = (dim==2) ? mesh->numFacetsTotal() : mesh->numCornersTotal();

    //printf("ENTRANDO NO LOOP eid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    for (int eid = 0; eid < n_edges_total; ++eid)
    {

      if (dim == 2)
        edge = mesh->getFacetPtr(eid);
      else
        edge = mesh->getCornerPtr(eid);

      if (edge==NULL || edge->isDisabled())
        continue;

      if (dim == 2)
        mesh->getFacetNodesId(eid, edge_nodes.data());
      else
        mesh->getCornerNodesId(eid, edge_nodes.data());

      tag_a = mesh->getNodePtr(edge_nodes[0])->getTag();
      tag_b = mesh->getNodePtr(edge_nodes[1])->getTag();
      tag_e = edge->getTag();

      if (!( tag_a==tag_b && tag_b==tag_e ) && (is_splitting==0))
        continue;

      Point const* pt_a = mesh->getNodePtr(edge_nodes[0]);
      Point const* pt_b = mesh->getNodePtr(edge_nodes[1]);

      if (pt_a->isMarked() || pt_b->isMarked())
        continue;

      pt_a->getCoord(Xa.data(),dim);
      pt_b->getCoord(Xb.data(),dim);
      
      h = (Xa-Xb).norm();

      expected_h = .5*(mesh_sizes[edge_nodes[0]] + mesh_sizes[edge_nodes[1]]);

      //printf("mesh size  %lf!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", (h - expected_h)/expected_h);

      if (is_splitting )//&& (time_step%2==0))
      {
        if ((h - expected_h)/expected_h > TOL)
        {
          mesh_was_changed = true;
          int pt_id = MeshToolsTri::insertVertexOnEdge(edge->getIncidCell(), edge->getPosition(), 0.5, &*mesh);
          //printf("INSERTED %d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", pt_id);
          adde_vtcs.push_back(make_tuple(pt_id, edge_nodes[0], edge_nodes[1]));
          mesh->getNodePtr(pt_id)->setMarkedTo(true);
          if (pt_id < (int)mesh_sizes.size())
            mesh_sizes[pt_id] = expected_h;
          else
          {
            if (pt_id > (int)mesh_sizes.size())
            {
              printf("ERROR: Something with mesh_sizes is wrong!!\n");
              throw;
            }
            mesh_sizes.push_back(expected_h);
          }
        }
      }
      else if(!is_splitting )//&& !(time_step%2==0)) // is collapsing
      {
        if ((h - expected_h)/expected_h < -TOL)
        {
          mesh_was_changed = true;
          //printf("COLLAPSED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
          //printf("COLLAPSING  tags = %d  %d %d #################################################\n", tag_a, tag_b, tag_e);
          int pt_id = MeshToolsTri::collapseEdge2d(edge->getIncidCell(), edge->getPosition(), 0.0, &*mesh);
          //printf("COLLAPSED %d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", pt_id);
          //mesh->getNodePtr(pt_id)->setMarkedTo(true);
        }
      }

    } // end n_edges_total
  }

  if (!mesh_was_changed)
    PetscFunctionReturn(0);

  // atualiza tudo!
  meshAliasesUpdate();

  // free petsc matrices first and residues to save memory.
  // destroy only those who not depends on DofHandler
  Destroy(Mat_Jac);
  Destroy(Mat_Jac_m);
  Destroy(Vec_res);
  Destroy(Vec_res_m);
  Destroy(Vec_normal);
  Destroy(Vec_v_mid);
  SNESReset(snes);
  SNESReset(snes_m);
  KSPReset(ksp);
  KSPReset(ksp_m);
  PCReset(pc);
  PCReset(pc_m);
  SNESLineSearchReset(linesearch);

  DofHandler  dof_handler_tmp[2];
  
  // tranfers variables values from old to new mesh
  {
    // First fix the u-p unknows
    dof_handler_tmp[DH_MESH].copy(dof_handler[DH_MESH]);
    dof_handler_tmp[DH_UNKS].copy(dof_handler[DH_UNKS]);

    dofsUpdate();


    Vec *petsc_vecs[] = {&Vec_up_0, &Vec_up_1, &Vec_x_0, &Vec_x_1};
    int DH_t[]       = {DH_UNKS , DH_UNKS , DH_MESH, DH_MESH};

    std::vector<Real> temp;

    // NOTE: the mesh must not be changed in this loop
    for (int v = 0; v < static_cast<int>( sizeof(DH_t)/sizeof(int) ); ++v)
    {
      int vsize;

      VecGetSize(*petsc_vecs[v], &vsize);

      temp.assign(vsize, 0.0);

      double *array;

      // copy
      VecGetArray(*petsc_vecs[v], &array);
      for (int i = 0; i < vsize; ++i)
        temp[i] = array[i];
      VecRestoreArray(*petsc_vecs[v], &array);

      Destroy(*petsc_vecs[v]);

      ierr = VecCreate(PETSC_COMM_WORLD, petsc_vecs[v]);                              CHKERRQ(ierr);
      ierr = VecSetSizes(*petsc_vecs[v], PETSC_DECIDE, dof_handler[DH_t[v]].numDofs()); CHKERRQ(ierr);
      ierr = VecSetFromOptions(*petsc_vecs[v]);                                         CHKERRQ(ierr);

      int dofs_0[64]; // old
      int dofs_1[64]; // new

      // copy data from old mesh to new mesh
      VecGetArray(*petsc_vecs[v], &array);
      for (point_iterator point = mesh->pointBegin(), point_end = mesh->pointEnd(); point != point_end; ++point)
      {
        if (!mesh->isVertex(&*point))
          continue;

        if (!point->isMarked())
        {
          for (int k = 0; k < dof_handler_tmp[DH_t[v]].numVars(); ++k)
          {
            //printf("v=%d,  point = %d, k=%d @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n  ", v, mesh->getPointId(&*point), k);
            //dof_handler_tmp[DH_t[v]].getVariable(k).getVertexAssociatedDofs(dofs_0, &*point);
            dof_handler_tmp[DH_t[v]].getVariable(k).getVertexDofs(dofs_0, mesh->getPointId(&*point));
            dof_handler    [DH_t[v]].getVariable(k).getVertexAssociatedDofs(dofs_1, &*point);
            for (int j = 0; j < dof_handler_tmp[DH_t[v]].getVariable(k).numDofsPerVertex(); ++j)
              array[dofs_1[j]] = temp.at(dofs_0[j]);

          }
        }
      }
      // interpolate values at new points
      std::list<EdgeVtcs>::iterator it     = adde_vtcs.begin();
      std::list<EdgeVtcs>::iterator it_end = adde_vtcs.end();
      for (; it != it_end; ++it)
      {
        int const pt_id        =  get<0>(*it);
        int const a_id         =  get<1>(*it);
        int const b_id         =  get<2>(*it);
        Point      * point = mesh->getNodePtr(pt_id);
        Point const* pt_a  = mesh->getNodePtr(a_id);
        Point const* pt_b  = mesh->getNodePtr(b_id);
        Point const* link[] = {pt_a, pt_b};
        int const Nlinks = sizeof(link)/sizeof(Point*);
        double const weight = 1./Nlinks;
        //printf("a_id = b_id %d %d!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", a_id, b_id);
        for (int k = 0; k < dof_handler_tmp[DH_t[v]].numVars(); ++k)
        {
          dof_handler[DH_t[v]].getVariable(k).getVertexAssociatedDofs(dofs_1, &*point);

          for (int j = 0; j < Nlinks; ++j)
          {
            dof_handler_tmp[DH_t[v]].getVariable(k).getVertexAssociatedDofs(dofs_0, link[j]);
            //printf("v = %d  x = %lf  y = %lf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", v, temp[dofs_0[0]], temp[dofs_0[1]]);
            for (int c = 0; c < dof_handler_tmp[DH_t[v]].getVariable(k).numDofsPerVertex(); ++c)
              array[dofs_1[c]] += weight*temp[dofs_0[c]];
          }
        }
        
      }
      VecRestoreArray(*petsc_vecs[v], &array);
      Assembly(*petsc_vecs[v]);
    } // end loop in vectors

    std::list<EdgeVtcs>::iterator it     = adde_vtcs.begin();
    std::list<EdgeVtcs>::iterator it_end = adde_vtcs.end();
    for (; it != it_end; ++it)
    {
      int const pt_id =  get<0>(*it);
      mesh->getNodePtr(pt_id)->setMarkedTo(false);
    }

    //for (point_iterator point = mesh->pointBegin(), point_end = mesh->pointEnd(); point != point_end; ++point)
    //{
    //  if (point->isMarked())
    //  {
    //    printf("%d MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA\n", mesh->getPointId(&*point));
    //  }
    //}


  } // end tranf


  //Vec Vec_res;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_res);                     CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_res, PETSC_DECIDE, n_unknowns);            CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_res);                                CHKERRQ(ierr);

  //Vec Vec_v_mid
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_v_mid);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_v_mid, PETSC_DECIDE, n_dofs_v_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_v_mid);                             CHKERRQ(ierr);

  //Vec Vec_normal;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_normal);                  CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_normal, PETSC_DECIDE, n_dofs_v_mesh);      CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_normal);                             CHKERRQ(ierr);

  //Vec Vec_res_m;
  ierr = VecCreate(PETSC_COMM_WORLD, &Vec_res_m);                     CHKERRQ(ierr);
  ierr = VecSetSizes(Vec_res_m, PETSC_DECIDE, n_dofs_v_mesh);         CHKERRQ(ierr);
  ierr = VecSetFromOptions(Vec_res_m);

  //Mat Mat_Jac;
  ierr = MatCreate(PETSC_COMM_WORLD, &Mat_Jac);                                      CHKERRQ(ierr);
  ierr = MatSetSizes(Mat_Jac, PETSC_DECIDE, PETSC_DECIDE, n_unknowns, n_unknowns);   CHKERRQ(ierr);
  ierr = MatSetFromOptions(Mat_Jac);                                                 CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(Mat_Jac,  max_nz, NULL);                          CHKERRQ(ierr);
  ierr = MatSetOption(Mat_Jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);           CHKERRQ(ierr);

  //Mat Mat_Jac_m;
  int n_mesh_dofs = dof_handler[DH_MESH].numDofs();
  ierr = MatCreate(PETSC_COMM_WORLD, &Mat_Jac_m);                                        CHKERRQ(ierr);
  ierr = MatSetType(Mat_Jac_m,MATSEQAIJ);                                                CHKERRQ(ierr);
  ierr = MatSetSizes(Mat_Jac_m, PETSC_DECIDE, PETSC_DECIDE, n_mesh_dofs, n_mesh_dofs);   CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(Mat_Jac_m,  max_nz_m, NULL);                          CHKERRQ(ierr);
  ierr = MatSetOption(Mat_Jac_m,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);             CHKERRQ(ierr);
  ierr = MatSetOption(Mat_Jac_m,MAT_SYMMETRIC,PETSC_TRUE);                               CHKERRQ(ierr);

  ierr = SNESSetFunction(snes, Vec_res, FormFunction, this);                             CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes, Mat_Jac, Mat_Jac, FormJacobian, this);                    CHKERRQ(ierr);
  ierr = SNESSetConvergenceTest(snes,CheckSnesConvergence,this,PETSC_NULL);              CHKERRQ(ierr);
  ierr = SNESGetKSP(snes,&ksp);                                                          CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);                                                              CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,Mat_Jac,Mat_Jac,SAME_NONZERO_PATTERN);                      CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);                                                       CHKERRQ(ierr);

  ierr = SNESSetFunction(snes_m, Vec_res_m, FormFunction_mesh, this);                    CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes_m, Mat_Jac_m, Mat_Jac_m, FormJacobian_mesh, this);         CHKERRQ(ierr);
  ierr = SNESGetLineSearch(snes_m,&linesearch);                                      CHKERRQ(ierr);
  ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHBASIC);                          CHKERRQ(ierr);
  ierr = SNESGetKSP(snes_m,&ksp_m);                                                  CHKERRQ(ierr);
  ierr = KSPGetPC(ksp_m,&pc_m);                                                      CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp_m,Mat_Jac_m,Mat_Jac_m,SAME_NONZERO_PATTERN);            CHKERRQ(ierr);
  //ierr = KSPSetType(ksp_m,KSPCG);                                                    CHKERRQ(ierr);
  //ierr = PCSetType(pc_m,PCILU);                                                      CHKERRQ(ierr);

  if(!nonlinear_elasticity)
  {
    ierr = SNESSetType(snes_m, SNESKSPONLY); CHKERRQ(ierr);
  }


  // check
  if (false)
  {
    int const n_edges_total = (dim==2) ? mesh->numFacetsTotal() : mesh->numCornersTotal();
    
    //printf("ENTRANDO NO LOOP eid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    for (int eid = 0; eid < n_edges_total; ++eid)
    {

      if (dim == 2)
        edge = mesh->getFacetPtr(eid);
      else
        edge = mesh->getCornerPtr(eid);

      if (edge==NULL || edge->isDisabled())
        continue;

      int tag = edge->getTag();

      if (tag != 2 && tag != 3 && tag != 7)
      {
        printf(" MERDA  MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA\n");
        printf("%d\n", tag);
        throw;
      }
      
      if (tag != 2 && tag != 3 && mesh->inBoundary((Facet*)edge))
      {
        printf(" MERDA  MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA MERDA\n");
        printf("%d\n", tag);
        throw;
      }
      
    }    
    
  }

  checkConsistencyTri(&*mesh);

  PetscFunctionReturn(0);
}

//PetscErrorCode AppCtx::calcMeshVelocity(Vec const& Vec_x_0, Vec const& Vec_x_1, Vec &Vec_v_mid)
//{
//  PetscErrorCode ierr;
//  ierr = VecCopy(Vec_x_1, Vec_v_mid);      CHKERRQ(ierr);
//  ierr = VecAXPY(Vec_v_mid,-1.,Vec_x_0);   CHKERRQ(ierr);
//  ierr = VecScale(Vec_v_mid, 1./dt);       CHKERRQ(ierr);
//  Assembly(Vec_v_mid);
//
//  // DEBUG
//  //Vector      U0(dim);
//  //Vector      X0(dim);
//  //Vector      X1(dim);
//  //VectorXi    node_dofs_fluid(dim);
//  //getNodeDofs(&*mesh->getNodePtr(120), DH_UNKS, VAR_U, node_dofs_fluid.data());
//  //VecGetValues(Vec_v_mid,  dim, node_dofs_fluid.data(), U0.data());
//  //VecGetValues(Vec_x_0,  dim, node_dofs_fluid.data(), X0.data());
//  //VecGetValues(Vec_x_1,  dim, node_dofs_fluid.data(), X1.data());
//  //cout << "VELOCITY VELOCITY VELOCITY VELOCITY VELOCITY VELOCITY = " << U0(0) << " " << U0(1) << endl;
//  //cout << "VELOCITY VELOCITY VELOCITY VELOCITY VELOCITY VELOCITY = " << X0(0) << " " << X0(1) << endl;
//  //cout << "VELOCITY VELOCITY VELOCITY VELOCITY VELOCITY VELOCITY = " << X1(0) << " " << X1(1) << endl;
//  //cout.flush();
//
//
//  PetscFunctionReturn(0);
//}
//

// by elasticity
// @brief A mean between Vec_up_0 and Vec_up_1 is used as boundary conditions to an elasticity problem. This mean
// is computed by = vtheta*Vec_up_1 + (1-Vec_up_1)*Vec_up_0

PetscErrorCode AppCtx::calcMeshVelocity(Vec const& Vec_x_0, Vec const& Vec_up_0, Vec const& Vec_up_1, double vtheta, Vec &Vec_v_mid, double tt)
{
  PetscErrorCode ierr;

  // The normal used must be at time step n.

  // set boundary conditions on the initial guess
  {
    VectorXi    node_dofs_mesh(dim);
    VectorXi    node_dofs_fluid(dim);
    int tag;

    Vector      U0(dim);
    Vector      X0(dim);
    Vector      U1(dim);
    Vector      tmp(dim);
    Vector      k1(dim),k2(dim),k3(dim),k4(dim); //  RK4

    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    for (; point != point_end; ++point)
    {
      tag = point->getTag();

      getNodeDofs(&*point, DH_MESH, VAR_M, node_dofs_mesh.data());
      getNodeDofs(&*point, DH_UNKS, VAR_U, node_dofs_fluid.data());

      VecGetValues(Vec_x_0, dim, node_dofs_mesh.data(), X0.data());

      if (force_mesh_velocity)
      {
        k1 = v_exact(X0,tt,tag);
        k2 = v_exact(X0+0.5*k1*dt,tt+0.5*dt,tag);
        k3 = v_exact(X0+0.5*k2*dt,tt+0.5*dt,tag);
        k4 = v_exact(X0+k3*dt,tt+dt,tag);
        tmp =  (k1 + 2.*(k2+k3) + k4)/6.; // velocity
        //tmp = v_exact(X0,tt,tag);
        
        //if (!mesh->isVertex(&*point)) // APAGAR TEMP ERASE-ME ... este codigo é só para não entortar a malha
        //{
        //  int vtcs[3];
        //  const int m = point->getPosition() - mesh->numVerticesPerCell();
        //  Cell const* cell = mesh->getCellPtr(point->getIncidCell());
        //  if (dim==3)
        //    cell->getCornerVerticesId(m, vtcs);
        //  else
        //    cell->getFacetVerticesId(m, vtcs);
        //  if (mesh->inBoundary(mesh->getNodePtr(vtcs[0])) || mesh->inBoundary(mesh->getNodePtr(vtcs[1])))
        //    tmp = tmp / 2.;
        //  if (mesh->inBoundary(mesh->getNodePtr(vtcs[0])) && mesh->inBoundary(mesh->getNodePtr(vtcs[1])))
        //    tmp = tmp * 0.;
        //}
        VecSetValues(Vec_v_mid, dim, node_dofs_mesh.data(), tmp.data(), INSERT_VALUES);
      }
      else
      {
        if (  false &&  (  is_in(tag, neumann_tags) || is_in(tag, dirichlet_tags) || is_in(tag, periodic_tags)   )   )
        {
          tmp.setZero();
          VecSetValues(Vec_v_mid, dim, node_dofs_mesh.data(), tmp.data(), INSERT_VALUES);
        }
        else
        //if (is_in(tag, interface_tags) || is_in(tag, triple_tags) || is_in(tag, solid_tags) || is_in(tag,feature_tags))
        {
          VecGetValues(Vec_up_0,  dim, node_dofs_fluid.data(), U0.data());
          VecGetValues(Vec_up_1,  dim, node_dofs_fluid.data(), U1.data());

          tmp = vtheta*U1 + (1.-vtheta)*U0;
          VecSetValues(Vec_v_mid, dim, node_dofs_mesh.data(), tmp.data(), INSERT_VALUES);
        }

      } // if force_mesh_velocity


    } // end for point

  } // b.c.

  if (!force_mesh_velocity)
  {
    ierr = SNESSolve(snes_m, PETSC_NULL, Vec_v_mid);  CHKERRQ(ierr);
  }


  // don't let edges to be curved
  if (mesh->numNodesPerCell() > mesh->numVerticesPerCell())
  {
    int tag;
    int dofs[50];
    Vector Vm(dim);
    Vector V0(dim);
    Vector V1(dim);

    point_iterator point = mesh->pointBegin();
    point_iterator point_end = mesh->pointEnd();
    for (; point != point_end; ++point)
    {
      tag = point->getTag();

      if (mesh->isVertex(&*point))
        continue;

      if (is_in(tag, solid_tags) || is_in(tag, triple_tags) || is_in(tag, interface_tags) || is_in(tag, feature_tags))
        continue;


      const int m = point->getPosition() - mesh->numVerticesPerCell();
      Cell const* cell = mesh->getCellPtr(point->getIncidCell());
      if (dim==3)
      {
        const int edge_id = cell->getCornerId(m);
        dof_handler[DH_MESH].getVariable(VAR_M).getCornerDofs(dofs, mesh->getCornerPtr(edge_id));
      }
      else
      {
        const int edge_id = cell->getFacetId(m);
        dof_handler[DH_MESH].getVariable(VAR_M).getFacetDofs(dofs, mesh->getFacetPtr(edge_id));
      }
      VecGetValues(Vec_v_mid, dim, dofs + 0*dim, V0.data());
      VecGetValues(Vec_v_mid, dim, dofs + 1*dim, V1.data());
      Vm = .5*(V0+V1);
      VecSetValues(Vec_v_mid, dim, dofs + 2*dim, Vm.data(), INSERT_VALUES);

    }
    Assembly(Vec_v_mid);
  }


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
      //  Cell const* cell = mesh->getCellPtr(point->getIncidCell());
      //  if (dim==3)
      //    cell->getCornerVerticesId(m, vtcs);
      //  else
      //    cell->getFacetVerticesId(m, vtcs);
      //  if (mesh->inBoundary(mesh->getNodePtr(vtcs[0])) || mesh->inBoundary(mesh->getNodePtr(vtcs[1])))
      //    tmp = tmp / 2.;
      //  if (mesh->inBoundary(mesh->getNodePtr(vtcs[0])) && mesh->inBoundary(mesh->getNodePtr(vtcs[1])))
      //    tmp = tmp * 0.;
      //}
      tmp = X0 + tmp;
      VecSetValues(Vec_x_new, dim, node_dofs_mesh.data(), tmp.data(), INSERT_VALUES);
    }
    else
    {

      if (is_in(tag, neumann_tags) || is_in(tag, dirichlet_tags) || is_in(tag, periodic_tags))
      {
        VecSetValues(Vec_x_new, dim, node_dofs_mesh.data(), X0.data(), INSERT_VALUES);
      }
      else
      //if (is_in(tag, interface_tags) || is_in(tag, triple_tags) || is_in(tag, solid_tags) || is_in(tag,feature_tags))
      if (true)
      {
        VecGetValues(Vec_up_0,  dim, node_dofs_fluid.data(), U0.data());
        VecGetValues(Vec_up_1,  dim, node_dofs_fluid.data(), U1.data());

        tmp = X0 + dt*(  vtheta*U1 + (1.-vtheta)*U0 );
        VecSetValues(Vec_x_new, dim, node_dofs_mesh.data(), tmp.data(), INSERT_VALUES);
      }
      else // volume
      {
        VecGetValues(Vec_v_mid,   dim, node_dofs_mesh.data(), U0.data()); // get old vmesh

        //double const a = 0.001351644*X0(0)*(54.4 - X0(0));
        tmp = X0 + dt*U0 ;
        VecSetValues(Vec_x_new, dim, node_dofs_mesh.data(), tmp.data(), INSERT_VALUES);
      }
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

  dof_handler[DH_MESH].getVariable(VAR_M).getCellDofs(mapM_c.data(), mesh->getCellPtr(cell_id));

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







































