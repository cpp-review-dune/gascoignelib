#ifndef __dgedge_h__
#define __dgedge_h__


namespace Gascoigne
{

  //////////////////////////////////////////////////
  // basic class to store the information of one edge for integration in dg discretizations
  //
  // master and slave are indices of the adjacent cells. It always holds master > slave
  class DGEdge
  {
  protected:
    int         _master,_slave;
    signed char _localmaster,_localslave;
  public:


  DGEdge() : _master(-1), _slave(-1), _localmaster(-1), _localslave(-1) {}
    ~DGEdge(){}
  
  DGEdge(int m, int s, signed char lm, signed char ls) : _master(m), _slave(s), _localmaster(lm), _localslave(ls)
    {
      sort_edge();
    }
  
    //////
    int master() const { return _master; }
    int slave()  const { return _slave; }
    int& master()      { return _master; }
    int& slave()       { return _slave; }

    signed char  localmaster() const { return _localmaster; }
    signed char  localslave()  const { return _localslave; }
    signed char& localmaster()      { return _localmaster; }
    signed char& localslave()       { return _localslave; }

    void sort_edge()
    {
      assert(_master!=_slave);
      if (_master<_slave)
	{
	  std::swap(_master      , _slave);
	  std::swap(_localmaster , _localslave);
	}
    }
    bool operator=(const DGEdge& E) const
      {
	if (! ((_master==E.master())&&(_slave==E.slave())) ) return false;
	assert(_localmaster==E.localmaster());
	assert(_localslave==E.localslave());
	return true;
      }
    
    
  };
}

  
#endif
