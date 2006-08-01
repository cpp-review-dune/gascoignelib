#include "piq2.h"

using namespace std;

namespace Gascoigne
{
/**********************************************************/

PiQ2::PiQ2() : _DI(NULL), _MP(NULL)
{
  _q2weight.resize(3,DoubleVector(5));
  _q2weight[0][0]=1.;
  _q2weight[0][1]=3./8.;
  _q2weight[0][2]=0;
  _q2weight[0][3]=-1./8.;
  _q2weight[0][4]=0;

  _q2weight[1][0]=0.;
  _q2weight[1][1]=3./4.;
  _q2weight[1][2]=1.;
  _q2weight[1][3]=3./4.;
  _q2weight[1][4]=0;

  _q2weight[2][0]=0.;
  _q2weight[2][1]=-1./8;
  _q2weight[2][2]=0;
  _q2weight[2][3]=3./8.;
  _q2weight[2][4]=1.;
}

/**********************************************************/

void PiQ2::Init(const SolverInterface *SI)
{
  _DI = SI->GetDiscretization();
  assert(_DI);
  _MP = dynamic_cast<const GascoigneMesh *>(SI->GetMesh());
  assert(_MP);
  assert(_MP->HasQ4Patch());
}

/**********************************************************/

void PiQ2::vmult(GlobalVector &pu, const GlobalVector &u) const
{
  int npatch=25,MAXZ=1,MAXIZ=1;
  if(_MP->dimension()==3)
  {
    npatch = 125;
    MAXZ   = 5;
    MAXIZ  = 3;
  }

  pu.zero();

  for(int p=0; p<_MP->nq4patches(); p++)
  {
    const IntVector &q4patch = *_MP->IndicesOfQ4Patch(p);
    assert(q4patch.size()==npatch);
    for(int jx=0; jx<5; jx++)
    {
      for(int jy=0; jy<5; jy++)
      {
        for(int jz=0; jz<MAXZ; jz++)
        {
          int n = 25*jz+jy*5+jx;
          int node = q4patch[n];
          if((jx%2==0)&&(jy&2==0)&&(jz%2==0))
          {
            pu.zero_node(node);
          }
          assert(node<pu.n());
          assert(node<u.n());
          pu.equ_node(node,1.,node,u);
          for(int ix=0; ix<3; ix++)
          {
            for(int iy=0; iy<3; iy++)
            {
              for(int iz=0; iz<MAXIZ; iz++)
              {
                int nodebig = q4patch[ix*2+iy*2*5+iz*2*25];
                double weight = _q2weight[ix][jx]*_q2weight[iy][jy];
                if(_MP->dimension()==3)
                {
                  weight *= _q2weight[iz][jz];
                }
                pu.add_node(node,-weight,nodebig,u);
              }
            }
          }
        }
      }
    }
  }
}

/**********************************************************/
}

