/*----------------------------  vanka_matrix_vector_base.h
 * ---------------------------*/
/*      $Id:$                 */
#ifndef __vanka_matrix_vector_base_H
#define __vanka_matrix_vector_base_H
/*----------------------------   vanka_matrix_vector_base.h
 * ---------------------------*/

#include "columndiagstencil.h"
#include "fmatrixblock.h"
#include "gascoignehash.h"
#include "gascoignemesh2d.h"
#include "gascoignemesh3d.h"
#include "sparseblockmatrix.h"
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <map>
#include <set>

#include "iluinterface.h"

namespace Gascoigne {
class Vanka_Matrix_Vector_base {
private:
public:
  virtual void build_patches(const GascoigneMesh *M) {
    cout << "not written Vanka_Matrix_Vector_base" << endl;
    abort();
  }

  virtual void ReInit(const MatrixInterface *A_Matrix) {
    cout << "not written Vanka_Matrix_Vector_base" << endl;
    abort();
  }

  virtual Eigen::FullPivLU<Eigen::MatrixXd> &GetLU(int n) {
    cout << "not written Vanka_Matrix_Vector_base" << endl;
    abort();
  }

  virtual int patchindices_size() {
    cout << "not written Vanka_Matrix_Vector_base" << endl;
    abort();
  }

  virtual const std::vector<int> &Getpatchindices(int p) {
    cout << "not written Vanka_Matrix_Vector_base" << endl;
    abort();
  }

  virtual const HASHMAP<int, int> &GetINP(int p) {
    cout << "not written Vanka_Matrix_Vector_base" << endl;
    abort();
  }

  virtual const std::vector<int> &Get_Vector_Columns(int p) {
    cout << "not written Vanka_Matrix_Vector_base" << endl;
    abort();
  }

  virtual std::vector<int>::iterator Get_Vector_Columns_begin(int n) {
    cout << "not written Vanka_Matrix_Vector_base" << endl;
    abort();
  }

  virtual std::vector<int>::iterator Get_Vector_Columns_end(int n) {
    cout << "not written Vanka_Matrix_Vector_base" << endl;
    abort();
  }
  virtual void clear() {
    cout << "not written Vanka_Matrix_Vector_base" << endl;
    abort();
  }
};
} // namespace Gascoigne

/*----------------------------   vanka_matrix_vector_base.h
 * ---------------------------*/
/* end of #ifndef __vanka_matrix_vector_base_H */
#endif
/*----------------------------   vanka_matrix_vector_base.h
 * ---------------------------*/
