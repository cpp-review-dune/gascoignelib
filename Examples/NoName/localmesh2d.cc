#include  "localmesh2d.h"

using namespace std;

/*-----------------------------------------*/

LocalMesh2d::LocalMesh2d() : GascoigneMesh2d() 
{
  _count=0;
  _q.resize(10);
  for(int k=0;k<_q.size();k++)
    {
      _q[k].reservesize(10);
      _q[k].zero();
      _q[k][k]=0.25;
//       for(int n=0;n<_q[k].size();n++)
// 	{
// 	  _q[k][n] = 0.25*pow(1./(n+1.),k);
// 	}
    }
}

/*-----------------------------------------*/

void LocalMesh2d::SetCoordinates(const LevelMesh2d* LM)
{
//   cerr << "$$$$$$$$$$$$$$$$$$$$ " << _count << endl;
//   cerr << "$$$$$$$$$$$$$$$$$$$$ " << _q[_count] << endl;

  assert(_count<_q.size());

  double rect_min = 0.75;
  double rect_max = 3.;
  for(int i=0;i<nnodes();i++)
    {
      const Vertex2d& v = LM->vertex2d(i);
      double x = v.x();
      double y = v.y();

      double rect = max(abs(x),abs(y));
      if( (rect<=rect_max) && (rect>=rect_min) )
	{
	  double p = (rect_max-rect)/(rect_max-rect_min);
	  double q = 1.-p;

	  double beta = atan2(y,x);
	  double alpha = beta;
	  if(y<0) alpha = 2.*M_PI+beta;

// 	  cerr << "xy alpha " << x << " " << y << "  :::  " << alpha << " *** " << sin(alpha) << endl;


	  double xnew=0.;
	  double ynew=0.;

	  {
	    //  Choix I
	    double r = rect;
	    for(int n=0;n<_q[_count].size();n++)
	      {
		r += _q[_count][n] * (cos((n+1)*alpha)-1.);
	      }
	    xnew = q * y  +  p * r*sin(alpha);
	    ynew = q * x  +  p * r*cos(alpha);
	  }

	  nx[i].y() = xnew;
	  nx[i].x() = ynew;
	}
    }
  _count++;
}
