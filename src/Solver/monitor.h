#ifndef __monitor_h
#define __monitor_h

#include  <vector>
#include  <iostream>

#include  "cginfo.h"
#include  "nlinfo.h"
#include  "paramfile.h"
#include  <string>


/*******************************************************************/

class Monitor
{
  int    prec;
  int    m_length[10], m_counter;
  int    ps, _newmatrix;

 protected:

  typedef std::vector<double> dvector;
  typedef std::vector<int>    ivector;

  std::vector<std::string>  header;
  std::string          format;

  std::string directory;

  std::string protokoll, aos, bos;
  std::string texfile, numfile;

  void  error_io(const std::string&) const;
  void  matrix_info(int);
  void  print_message ();

  void  PrintAscii (std::ostream& os, const std::string&) const;
  void  PrintAscii (std::ostream& os, int    i) const {os <<i;} 
  void  PrintAscii (std::ostream& os, double d) const {os.precision(prec);os <<d;} 
  void  PrintAscii (std::ostream& os, const ivector&) const;
  void  PrintAscii (std::ostream& os, const dvector&) const;
  void  PrintHeader(std::ostream& os) const;

 public:

  int     control;
  char    message[400];

  Monitor();

  void SetAos(const std::string& s) {aos=s;}
  void SetBos(const std::string& s) {bos=s;}

  std::string&  aim_of_simulation()   { return aos;}
  int   set_print_step   (int i) { return ps = i;}

  void set_directory(const std::string& dir);

  void  failed_step   ();
  void  pre_monitor   (char*);
  void  post_monitor  ();
  void  init          (const Gascoigne::ParamFile* pf, int);
  void  mesh          (int,int);
  void  pre_nonlinear (int);
  void  post_nonlinear(const dvector&, double, int,int,int);
  void  nonlinear_step(const CGInfo&, const NLInfo&);
  int&  new_matrix()  { return _newmatrix;};

  
  void  Print          (const std::string& s, std::string se="\n") const;

  void precision(int n) {prec=n;}

  void  PrintResults   (const std::string& s="\n") const;
  void  PrintResults   (double) const;
  void  PrintResults   (int) const;
  void  PrintResults   (const ivector& iv) const;
  void  PrintResults   (const dvector& dv) const;

  void PrintInfoSummary(const NLInfo& nlinfo) const;
  void PrintInfoSummary(const CGInfo& nlinfo) const;
};

/*******************************************************************/

#endif
