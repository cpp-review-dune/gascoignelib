#ifndef  __StencilInterface_h
#define  __StencilInterface_h


namespace Gascoigne
{

  /////////////////////////////////////////////
  ////
  ////@brief
  ////  ... comments StencilInterface

  ////
  ////
  /////////////////////////////////////////////

  class StencilInterface
  {
    private:

    protected:

    public:
      //
      ////  Con(De)structor 
      //
      StencilInterface() {}
      virtual ~StencilInterface() {}
      virtual int n() const =0;
  };
}

#endif
