/*----------------------------   blockvanka.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __blockvanka_H
#define __blockvanka_H
/*----------------------------   blockvanka.h     ---------------------------*/

#include "iluinterface.h"

namespace Gascoigne {

template <int N>
class PatchVanka : public IluInterface {
public:
  BlockVanka() {
  }

  int          n() const                                                          = 0;
  std::string  GetName() const return "Block Vanka Smoother"                      = 0;
  virtual void ReInit(const SparseStructureInterface* A)                          = 0;
  virtual void ConstructStructure(const IntVector& perm, const MatrixInterface& A)= 0;

  virtual void modify(int c, double s) {
    std::cerr << "\"IluInterface::modify\" not written!" << std::endl;
    abort();
  }
  virtual void zero() {
    std::cerr << "\"IluInterface::zero\" not written!" << std::endl;
    abort();
  }
  virtual void compute_ilu() {
    std::cerr << "\"IluInterface::compute_ilu\" not written!" << std::endl;
    abort();
  }
  virtual void copy_entries(const MatrixInterface* A) {
    std::cerr << "\"IluInterface::copy_entries\" not written!" << std::endl;
    abort();
  }
  virtual void solve(GlobalVector& x) const {
    std::cerr << "\"IluInterface::solve\" not written!" << std::endl;
    abort();
  }
  virtual void solve_transpose(GlobalVector& x) const {
    std::cerr << "\"IluInterface::solve_transpose\" not written!" << std::endl;
    abort();
  }
  virtual std::ostream& Write(std::ostream& s) const {
    std::cerr << "\"IluInterface::Write\" not written!" << std::endl;
    abort();
  }
};

};  // namespace Gascoigne
}

/*----------------------------   blockvanka.h     ---------------------------*/
/* end of #ifndef __blockvanka_H */
#endif
/*----------------------------   blockvanka.h     ---------------------------*/
