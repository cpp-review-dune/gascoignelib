#include  "extrapolator.h"

/*-----------------------------------------*/

Extrapolator::Extrapolator() {}
Extrapolator::~Extrapolator() {}

/*-----------------------------------------*/

void Extrapolator::Print()
{
  int n = order.size();
  std::cout << "extra\t";

  std::cout.precision(12);
  for (int i=0; i<n; i++)
    {
      std::cout << valextra[i];
      std::cout.precision(4);
      std::cout << " [" << order[i] << "]  ";
      std::cout.precision(12);
    }
  std::cout << std::endl;
}

/*-----------------------------------------*/

void Extrapolator::NewValues(const dvector& J)
{
  int n = J.size();
  if (!vals.size())
    {
      valextra.resize(n);
      order.resize(n);
      valextra = 0.;
      order = -1.;
    }
  vals.push_back(J);
  int nv = vals.size();
  if(nv>=3)
    {
      for(int i=0;i<n;i++)
        {
          double j0 = vals[nv-3][i];
          double j1 = vals[nv-2][i];
          double j2 = vals[nv-1][i];

          double d0 = j0 - j1;
          double d1 = j1 - j2;
          double alpha = d0/d1;
          double c = alpha*d0/(alpha-1.);
          double j = j0-c;
          
          order[i] = alpha;
          valextra[i] = j;
        }
    }
}
