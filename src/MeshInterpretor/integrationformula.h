#ifndef __integrationformula_h
#define __integrationformula_h

#include  "integrationformulainterface.h"

/*------------------------------------------------------------*/

template<int DIM>
class IntegrationFormula : public IntegrationFormulaInterface
{
protected:

  typedef Vertex<DIM>   VERTEX;

  std::vector<VERTEX> ic;

public:

  IntegrationFormula();
  IntegrationFormula(int n);
  ~IntegrationFormula()    {}

  IntegrationFormula(const IntegrationFormula& IF) : 
    IntegrationFormulaInterface(IF), ic(IF.c()) {}

  void init(int n);

  void xi(VERTEX& v, int k)  const { v = ic[k];}
  const std::vector<VERTEX>& c()  const { return ic;}
};

/*------------------------------------------------------------*/

typedef IntegrationFormula<1>  IntegrationFormula1d;
typedef IntegrationFormula<2>  IntegrationFormula2d;
typedef IntegrationFormula<3>  IntegrationFormula3d;

/*------------------------------------------------------------*/

class LineMidPoint : public IntegrationFormula1d{
public:  LineMidPoint();};

/*------------------------------------------------------------*/

class LineTrapez : public IntegrationFormula1d{
public:  LineTrapez();};

/*------------------------------------------------------------*/

class LineSimpson : public IntegrationFormula1d{
public:  LineSimpson();};

/*------------------------------------------------------------*/

class LineGauss1 : public IntegrationFormula1d{
public:  LineGauss1();};

/*------------------------------------------------------------*/

class LineGauss2 : public IntegrationFormula1d{
public:  LineGauss2();};

/*------------------------------------------------------------*/

class LineGauss3 : public IntegrationFormula1d{
public:  LineGauss3();};

/*------------------------------------------------------------*/

class LineGauss4 : public IntegrationFormula1d{
public:  LineGauss4();};

/*------------------------------------------------------------*/

class LineGauss5 : public IntegrationFormula1d{
public:  LineGauss5();};

/*------------------------------------------------------------*/

class LineGauss6 : public IntegrationFormula1d{
public:  LineGauss6();};

/*------------------------------------------------------------*/

class LineGauss7 : public IntegrationFormula1d{
public:  LineGauss7();};

/*------------------------------------------------------------*/

class LineGauss8 : public IntegrationFormula1d{
public:  LineGauss8();};

/*------------------------------------------------------------*/

class LineGauss9 : public IntegrationFormula1d{
public:  LineGauss9();};

/*------------------------------------------------------------*/

class LineGauss10 : public IntegrationFormula1d{
public:  LineGauss10();};

/*------------------------------------------------------------*/

template<int N, class INT>
class PatchFormula2d : public IntegrationFormula2d{
public: PatchFormula2d();};

/*------------------------------------------------------------*/

template<int N, class LineFormula>
class TensorFormula2d : public IntegrationFormula2d
{
public: TensorFormula2d();};

/*------------------------------------------------------------*/

template<int N, class Line>
class TensorFormula3d : public IntegrationFormula3d{
public: TensorFormula3d();};

/*------------------------------------------------------------*/

template<int N, class INT>
class PatchFormula3d : public IntegrationFormula3d{
public: PatchFormula3d();};

/*------------------------------------------------------------*/

typedef TensorFormula2d<1,LineMidPoint> QuadMidPoint;
typedef TensorFormula2d<2,LineTrapez>   QuadTrapez;
typedef TensorFormula2d<3,LineSimpson>  QuadSimpson;

typedef TensorFormula2d<1,LineGauss1>   QuadGauss1; 
typedef TensorFormula2d<2,LineGauss2>   QuadGauss4; 
typedef TensorFormula2d<3,LineGauss3>   QuadGauss9;
typedef TensorFormula2d<4,LineGauss4>   QuadGauss16; 
typedef TensorFormula2d<5,LineGauss5>   QuadGauss25; 
typedef TensorFormula2d<6,LineGauss6>   QuadGauss36; 
typedef TensorFormula2d<7,LineGauss7>   QuadGauss49; 
typedef TensorFormula2d<8,LineGauss8>   QuadGauss64; 
typedef TensorFormula2d<9,LineGauss9>   QuadGauss81; 
typedef TensorFormula2d<10,LineGauss10> QuadGauss100; 

/*------------------------------------------------------------*/

typedef TensorFormula3d<2,LineTrapez> HexTrapez;
typedef TensorFormula3d<2,LineGauss2> HexGauss8; 
typedef TensorFormula3d<3,LineGauss3> HexGauss27; 
typedef TensorFormula3d<4,LineGauss4> HexGauss64; 

/*------------------------------------------------------------*/

#endif
