#ifndef __cgmatrix_h
#define __cgmatrix_h


template <class MT,class VT,class MEM, class INFO>
inline int CgMatrix(const MT& M, VT& x, const VT& b, MEM& mem, INFO& info)
{
  
  VT& g  = mem[0];
  VT& h  = mem[1];
  VT& d  = mem[2];
  VT& Ad = mem[3];

  int    it=0,reached = 0;
  double gh,alpha,beta;
  double     res;

  M.vmulteq(g,x);
  g.sequ(1.,-1.,b);
  res = g.norm();

  if(res==0.)      return 0;

  d = g;
  //M.precondition(d,g);

  gh  =  g*d;
  gh *=  -1.;

  for(it=0;!reached;it++)
  {
    M.vmulteq(Ad,d);

    alpha = d*Ad;
    alpha = gh/alpha;

    g.add(alpha,Ad);
    x.add(alpha,d );
    res = g.norm();

    reached = info.check_residual(it,res);
    if (reached) continue;
    
    h = g;
    //M.precondition(h,g);
    
    beta = gh;
    gh   = g*h;
    beta = gh/beta;
    
    d.sequ(beta,-1.,h);
  }
  if (reached<0) return 1;
  return 0;
}

#endif
