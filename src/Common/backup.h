#ifndef __backup_h
#define __backup_h

#include  "gascoigne.h"
#include  "compvector.h"
#include  <string>

using namespace Gascoigne;

/********************************************************************/

class WriteBackUp
{
 public:

  WriteBackUp(const GlobalVector&, const std::string&);
};

/********************************************************************/

class WriteBackUpBinary
{
 public:

  WriteBackUpBinary(const GlobalVector&, const std::string&);
};

/********************************************************************/

class ReadBackUp
{
 public:

  ReadBackUp(const std::string&, int&, int&);
  ReadBackUp(GlobalVector&, const std::string&);
};

/********************************************************************/

class ReadBackUpResize
{
 public:

  ReadBackUpResize(GlobalVector&, const std::string&);
};

/********************************************************************/

class ReadBackUpBinary
{
 public:

  ReadBackUpBinary(GlobalVector&, const std::string&);
};

/********************************************************************/

#endif
