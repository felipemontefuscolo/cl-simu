// This file is part of FEPiC++, a toolbox for finite element codes.
//
// FEPiC++ is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// FEPiC++ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// FEPiC++. If not, see <http://www.gnu.org/licenses/>.

#ifndef FEPIC_MYPETSC_HPP
#define FEPIC_MYPETSC_HPP

#include "petscsnes.h"
#include "petscksp.h"
#include <iostream>
#include <fstream>

template<class Vec>
void printVector(Vec const& v)
{
	for (int i=0; i<v.size(); ++i)
		std::cout << v[i]<< " ";
	std::cout << std::endl;
}

void inline MatrixInfo(Mat& K) {
  MatInfo info;
  double  mal, nz_a, nz_u;
  MatGetInfo(K,MAT_LOCAL,&info);
  mal  = info.nz_allocated;
  nz_a = info.nz_used;
  nz_u = info.nz_unneeded;
  std::cout << mal << " allocated\t\t" << nz_a << " used\t" << nz_u << " unneeded" << std::endl;
}

void inline View(Mat &K, const char* file_name, const char* obj_name) {

  PetscViewer viewer;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name, &viewer);
  PetscObjectSetName((PetscObject)K,obj_name);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE);
  MatView(K, viewer);
  PetscViewerDestroy(&viewer);
}

void inline View(Vec &v, const char* file_name, const char* obj_name) {

  PetscViewer viewer;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name, &viewer);
  PetscObjectSetName((PetscObject)v,obj_name);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(v, viewer);
  PetscViewerDestroy(&viewer);
}

// OBSOLETO D+++
void inline ViewXPM(Mat &K, const char* name = "MATRIX") {
  std::ofstream file(name);
  int width, heigth;

  MatGetSize(K, &heigth, &width);

  file << "/* XPM */" << std::endl;
  file << "static char * blarg_xpm[] = {" << std::endl;
  file << "\"" <<  width << " " << heigth << " 3 1\",";
  file << "\"a c #000000\"," << std::endl;
  file << "\"b c #FFFFFF\"," << std::endl;
  file << "\"c c #FF0000\"," << std::endl;

  for (int i=0;i<heigth-1;i++) {
    file << "\"";
    double val;
    for (int j=0; j<width; j++) {
      MatGetValue(K,i,j,&val);
      if (val==1) file << "c";
      else if (val != 0) file << "b";
      else file << "a";
    }
    file << "\"," << std::endl;
  }

  file << "\"";
  for (int j=0; j<width; j++) {
    double val;
    MatGetValue(K,heigth-1,j,&val);
    if (val==1) file << "c";
    else if (val != 0) file << "b";
    else file << "a";
  }
  file << "\"" << std::endl;
  file << "}";

  file.close();

}

//inline
//void MyMatRestoreArray(Mat mat,PetscScalar *v[])
//{
//  (*mat->ops->restorearray)(mat,v);
//}






// REFERÊNCIA:
// PETSC_DIR/src/mat/impls/aij/seq
// a: jac_array
void inline MyMatSumValues(PetscScalar* a, PetscInt *ai, PetscInt *aj, PetscInt m,const PetscInt im[],PetscInt n,const PetscInt in[],PetscScalar v[])
{
  PetscInt     *rp,k,low,high,t,row,nrow,i,col,l;
  PetscScalar    *ap,*aa = a;



  for (k=0; k<m; k++) { /* loop over rows */
    row  = im[k];
    if (row < 0) {v += n; continue;} /* SETERRQ1(PETSC_ERR_ARG_OUTOFRANGE,"Negative row: %D",row); */
    //if (row >= A->rmap->n) SETERRQ2(PETSC_ERR_ARG_OUTOFRANGE,"Row too large: row %D max %D",row,A->rmap->n-1);
    rp   = aj + ai[row]; ap = aa + ai[row];
    //nrow = ailen[row];
    nrow = ai[row+1]-ai[row];
    for (l=0; l<n; l++) { /* loop over columns */
      if (in[l] < 0) {v++; continue;} /* SETERRQ1(PETSC_ERR_ARG_OUTOFRANGE,"Negative column: %D",in[l]); */
      //if (in[l] >= A->cmap->n) SETERRQ2(PETSC_ERR_ARG_OUTOFRANGE,"Column too large: col %D max %D",in[l],A->cmap->n-1);
      col = in[l] ;
      high = nrow; low = 0; /* assume unsorted */
      while (high-low > 5) {
        t = (low+high)/2;
        if (rp[t] > col) high = t;
        else             low  = t;
      }
      for (i=low; i<high; i++) {
        if (rp[i] > col) break;
        if (rp[i] == col) {
          ap[i] += *v++;
          goto finished;
        }
      }
      *v++ = 0.0;
      finished:;
    }
  }
}

void inline MyVecGetValues(double const* array,int mapsize, int const* mapdata, double * data)
{
  for (int i = 0; i < mapsize; ++i)
    data[i] = array[mapdata[i]];
}

void inline MyVecSetValues(double * array,int mapsize, int const* mapdata, double const* data)
{
  for (int i = 0; i < mapsize; ++i)
    array[mapdata[i]] = data[i];
}

void inline VecSumValues(double *f_array, int mapsize, int const* mapdata, double const* data)
{
  for (int i = 0; i < mapsize; ++i)
    f_array[mapdata[i]] += data[i];
}


void inline Assembly(Mat &K) {
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
};

void inline Assembly(Vec &v) {
  VecAssemblyBegin(v);
  VecAssemblyEnd(v);
}

void inline Destroy(Mat &K) {
  MatDestroy(&K);
}

void inline Destroy(Vec &v) {
  VecDestroy(&v);
}

void inline Destroy(PC &pc) {
  PCDestroy(&pc);
}

void inline Destroy(KSP &ksp) {
  KSPDestroy(&ksp);
}

void inline Destroy(SNES &snes) {
  SNESDestroy(&snes);
}



#endif

