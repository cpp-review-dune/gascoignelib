#ifndef  __FunctionalManager_h
#define  __FunctionalManager_h

#include  "functional.h"
#include  "dirichletdata.h"
#include  "righthandsidedata.h"
#include  "equation.h"
#include  "nvector.h"

/*-----------------------------------------*/


class FunctionalManager
{
protected:

  typedef nvector<double>     Vector;
  typedef nvector<int>        ivector;

  std::vector<Functional*>         FF;
  std::vector<std::string>             names, gnames;

  virtual Functional* ConstructFunctional
    (const std::string& name, const std::string& params, const Equation& EQ);

  void Construct(int i, const std::vector<std::string>& functional, const Equation& EQ);

public:

  FunctionalManager();
  ~FunctionalManager();

  void Print(std::ostream& os) const;

  void ConstructSet(const std::string& paramfile, const Equation& EQ);

  const Functional*   GetFunctional  (const std::string& name) const;

  const std::vector<std::string>& GetGridFunctionalNames() const {return gnames;}
  const std::vector<std::string>& GetFunctionalNames    () const {return names;}
};


#endif
