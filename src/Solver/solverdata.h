#ifndef  __SolverData_h
#define  __SolverData_h

#include  "gascoigne.h"
#include  "paramfile.h"
#include  <string>

/*----------------------------------------------*/

namespace Gascoigne
{
class SolverData
{
  protected:

    int  exact_lu, enlarge, iter_pre, iter_post, iter_exact;
    int  bicgstab_pstep;
    DoubleVector     ilum;
    DoubleVector     vector_direction;
    IntVector        stream_direction, pfilter;
    double              omega;
    std::string             ilusort, linear_smooth, bicgstab_residual;

    // gibt an wieviele Iterationen mindestens gemacht werden, und
    // ab wann im bicgstab ein Abbruchskriterium greift.
    // iter/miniter wird immer iteriert.
    double              bicgstab_miniter;

  public:

    SolverData();
    void Init(const ParamFile* pf, int ncomp);
    ~SolverData();

    const IntVector& GetPfilter()const { return pfilter;}
    int    GetExactLu()             const { return exact_lu;}
    int    GetEnlarge()             const { return enlarge;}
    double GetIluModify(int c)      const { return ilum[c];}
    const DoubleVector& GetIluModify() const { return ilum;}
    double GetOmega()           const { return omega;}
    int    GetIterPre ()        const { return iter_pre;}
    int    GetIterPost ()       const { return iter_post;}
    int    GetIterExact ()      const { return iter_exact;}
    const std::string& GetIluSort() const { return ilusort;}
    const IntVector& GetStreamDirection() const
      { return stream_direction; }
    const DoubleVector& GetVectorDirection() const
      { return vector_direction; }
    
    int    GetBiCGStabPStep()    const { return bicgstab_pstep;}
    std::string GetBiCGStabResidual()const { return bicgstab_residual;}
    std::string GetLinearSmooth()    const { return linear_smooth;}
    double GetBiCGStabMinIter()  const { return bicgstab_miniter;}
};
}

#endif
