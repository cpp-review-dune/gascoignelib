#ifndef __nlinfo_h
#define __nlinfo_h

#include "cginfo.h"

/*************************************************************/

class NLStatisticData : public StatisticData
{
  int  _newmatrix;
  int  _totalmatrix;

 public:

  int& totalmatrix()       { return _totalmatrix; }
  int  totalmatrix() const { return _totalmatrix; }
  int& newmatrix()       { return _newmatrix; }
  int  newmatrix() const { return _newmatrix; }

  void reset();

  friend std::ostream& operator<<(std::ostream &s, const NLStatisticData& A);
};

/*************************************************************/

class NLControlData : public ControlData
{
  int _relax, _newmatrix, _laststepbad, _matrixmustbebuild;

 public:

  NLControlData();
  int  relax()     const { return _relax;}
  int& relax()           { return _relax;}
  int  newmatrix() const { return _newmatrix;}
  int& newmatrix()       { return _newmatrix;}
  int  laststepbad() const { return _laststepbad;}
  int& laststepbad()       { return _laststepbad;}

  int  matrixmustbebuild() const { return _matrixmustbebuild;}
  int& matrixmustbebuild()       { return _matrixmustbebuild;}

  void reset();

  friend std::ostream& operator<<(std::ostream &s, const NLControlData& A);
};

/*************************************************************/

class NLUserData : public UserData
{
  double _rho, _linrho, _maxresincrease;
  int    _maxrelax;

 public:

  NLUserData();

  double  rho()      const { return _rho;}
  double& rho()            { return _rho;}
  double  linrho()      const { return _linrho;}
  double& linrho()            { return _linrho;}
  int     maxrelax() const { return _maxrelax;}
  int&    maxrelax()       { return _maxrelax;}

  double  maxresincrease() const { return _maxresincrease;}
  double& maxresincrease()       { return _maxresincrease;}

  friend std::ostream& operator<<(std::ostream &s, const NLUserData& A);
};

/*************************************************************/

class NLInfo
{
  CGInfo&        Linfo;

  NLStatisticData SD;
  NLControlData   CD;
  NLUserData      UD;

  void compute_reduction_rate();
  void matrix_control();

 public:

  NLInfo(CGInfo&, double f = 1.e-8, double r=1.e-15, int p = 10, 
	 int m = 5000, const std::string& txt = "NLIt: ");

  CGInfo& GetLinearInfo() {return Linfo;}
  const CGInfo& GetLinearInfo() const {return Linfo;}

  friend std::ostream& operator<<(std::ostream &s, const NLInfo& A);

  const NLStatisticData& statistics() const { return SD; }
        NLStatisticData& statistics()       { return SD; }
  const NLControlData&   control()    const { return CD; }
        NLControlData&   control()          { return CD; }
  const NLUserData&      user()       const { return UD; }
        NLUserData&      user()             { return UD; }

 bool   check         (double res, double cor = 0.);
 bool   check         (int iter, double res, double cor = 0.);
 std::string check_damping (int, double);
 void reset();
 void new_matrix();
};

/*************************************************************/

#endif
