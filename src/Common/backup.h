#ifndef __backup_h
#define __backup_h

#include  "gascoigne.h"
#include  "compvector.h"
#include  <string>

/********************************************************************/

class WriteBackUp
{
 public:

  WriteBackUp(const Gascoigne::GlobalVector&, const std::string&);
};

/********************************************************************/

class WriteBackUpBinary
{
 public:

  WriteBackUpBinary(const Gascoigne::GlobalVector&, const std::string&);
};

/********************************************************************/

class ReadBackUp
{
 public:

  ReadBackUp(const std::string&, int&, int&);
  ReadBackUp(Gascoigne::GlobalVector&, const std::string&);
};

/********************************************************************/

class ReadBackUpResize
{
 public:

  ReadBackUpResize(Gascoigne::GlobalVector&, const std::string&);
};

/********************************************************************/

class ReadBackUpBinary
{
 public:

  ReadBackUpBinary(Gascoigne::GlobalVector&, const std::string&);
};

/********************************************************************/

#endif
