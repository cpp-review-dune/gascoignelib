#include  "dirichletdatabycolor.h"

using namespace std;

/*-----------------------------------------*/

namespace Gascoigne
{
  DirichletDataByColor::DirichletDataByColor(nvector<int> comps, const set<int>& cl, nvector<double> s)
    : DirichletData(), __cols(cl), __comps(comps), __scales(s)
  {
    assert(__comps.size()==__scales.size());
    assert(__comps.size()>0);
    assert(__cols.size()>0);
  }

  DirichletDataByColor::DirichletDataByColor(int comps, std::set<int>& cl, double s) : __cols(cl)
  {
    __comps .push_back(comps);
    __scales.push_back(s);
  }

  /*-----------------------------------------*/

  DirichletDataByColor::DirichletDataByColor(const vector<string>& args)
  {
    bool ok=true;
    
    int n = args.size(); 
    if (n<5) ok=false;
    int ncol = atoi(args[0].c_str());
    if (n<4+ncol) ok = false;
    for (int i=0;i<ncol;++i)
      __cols.insert(atoi(args[i+1].c_str()));
    int ncomp = atoi(args[ncol+1].c_str());
    if (n!=2*ncomp+ncol+2) ok = false;
    for (int i=0;i<ncol;++i)
      {
	__comps.push_back( atoi(args[2*i+ncol+2].c_str()));
	__scales.push_back(atof(args[2*i+ncol+3].c_str()));
      }
    
    if (!ok)
      {
	cerr << "DirichletDataByColor::DirichletDataByColor: Usage" << endl
	     << "\t ncols_[cols]_ncomps_[comp_weight] " << endl;
	abort();
      }
  }

  /*-----------------------------------------*/

  void DirichletDataByColor::operator()
    (DoubleVector& b, const Vertex2d& v, int color) const
  {
    b.zero();

    if(__cols.find(color)!=__cols.end()) 
      for (int i=0;i<__comps.size();++i)
	b[__comps[i]] = __scales[i];
  }

  /*-----------------------------------------*/

  void DirichletDataByColor::operator()
    (DoubleVector& b, const Vertex3d& v, int color) const
  {
    b.zero();
    if(__cols.find(color)!=__cols.end()) 
      for (int i=0;i<__comps.size();++i)
	b[__comps[i]] = __scales[i];
  }

  /*-----------------------------------------*/

  set<int> DirichletDataByColor::preferred_colors()const 
  {
    return __cols;
  }
}
