#ifndef __leveljumper_h

#include <math.h>
#include <map>
#include "quad.h"
#include "hex.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class LevelJump
{
  int mini, maxi;

public:

  int  Min() const { return mini;}
  int  Max() const { return maxi;}
  int& Min()       { return mini;}
  int& Max()       { return maxi;}
  const int& Maxi() const  { return maxi;}

  void update(int i)           { maxi = std::max(maxi,i); mini=std::min(mini,i); }
  void reinit()   { maxi=0; mini=1000;}

  LevelJump() { reinit();}
  
};

/*---------------------------------------------------*/

class LevelJumper
{
  std::map<int,LevelJump>  Phi;

public:

  LevelJumper() : Phi() {} 

  void update(const Quad& Q)
    {
      for (int i=0; i<4; i++)
	{
	  Phi[ Q[i] ].update(Q.level());
	}
    }
  void update(const Hex& Q)
    {
      for (int i=0; i<8; i++)
	{
	  Phi[ Q[i] ].update(Q.level());
	}
    }
  bool check() const
    {
      for(std::map<int,LevelJump>::const_iterator p=Phi.begin();
	  p!=Phi.end();p++)
	{
	  if(p->second.Max()>p->second.Min()+1) return 1;
	}
      return 0;
    }
  LevelJump& operator()(int i)  { return Phi[i];}
  const LevelJump& operator()(int i) const { return Phi.find(i)->second;}

  bool VertexOK(const Quad& Q)
    {
      for(int ii=0;ii<4;ii++)
	{
	  if( Phi[ Q[ii] ].Maxi() > Q.level()+1 )
	    {
// 	      std::cerr << "§§§\t" << Q[ii] << " " << Phi[ Q[ii] ].Maxi() << " " << Q.level() << std::endl;
	      return 0;
	    }
	}
      return 1;
    }
  bool VertexOK(const Hex& Q)
    {
      for(int ii=0;ii<8;ii++)
	{
	  if( Phi[ Q[ii] ].Maxi() > Q.level()+1 )
	    {
	      return 0;
	    }
	}
      return 1;
    }


};
}

/*---------------------------------------------------*/

#endif
