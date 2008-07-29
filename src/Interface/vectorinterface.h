#ifndef  __VectorInterface_h
#define  __VectorInterface_h

#include  <string>
#include  <iostream>
#include  <cassert>


namespace Gascoigne
{

  /////////////////////////////////////////////
  ////
  ////@brief
  ////  ... comments VectorInterface
  ////
  ////
  /////////////////////////////////////////////

  class VectorInterface : public std::string
  {
    private:
      std::string _type;

    protected:

    public:
      VectorInterface(const std::string& name) : std::string(name), _type("node")  { }
      VectorInterface(const std::string& name, const std::string& type) : std::string(name) {
        assert(type=="node" || type=="cell" || type=="parameter");
        GetType()=type;
      }
      VectorInterface(const VectorInterface& v) {
        SetName(v.GetName());
        SetType(v.GetType());
      }
      virtual ~VectorInterface() { }

      void SetName(const std::string& name) { GetName()=name; }
      std::string& GetName() { return *this; }
      const std::string& GetName() const { return *this; }

      void SetType(const std::string& type) {
        assert(type=="node" || type=="cell" || type=="parameter");
        GetType()=type;
      }
      std::string& GetType() { return _type; }
      const std::string& GetType() const { return _type; }

      friend std::ostream& operator<<(std::ostream& os, const VectorInterface& g) {
        os << "Name: '" << g.GetName() << "' ";
        os << "Type: '" << g.GetType() << "'" << std::endl;
        return os;
      }

  };
}

#endif
