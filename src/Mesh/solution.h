#ifndef __solution_h
#define __solution_h

/***************************************************************/

class Solution
{
 protected:

 public:

  virtual ~Solution(){}

  // fier interpolarion
  virtual void equ(int, double, int, double, int) {}
  virtual void equ(int, double, int, double, int, double, int, double, int) {}
  virtual void equ(int, double, int, int, int, int, int, int, int, int) {}


  // fuer visualization
  virtual int    comp()     const {}
  virtual int    n()        const {}
  virtual double u(int i,int c) const {}
  virtual double z(int i,int c) const {}

  // fuer backup
  virtual double& u(int i,int c) {}
  virtual double& z(int i,int c) {}
};

/***************************************************************/

#endif
