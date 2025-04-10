#ifndef MULTIGRID_H_
#define MULTIGRID_H_

#include <cmath>
#include <iostream>
#include "Array.h"
#include "InputOutput.h"
#include "Nutrients.h"
using namespace std;

class Multigrid;

class Multigrid
{
public:
    Multigrid(int levelsAga, int levelsAgarMG, int nx, int ny, int nzAgar);
    ~Multigrid();
    void resetZero();
    void relaxationAgarAtLevel(int l);
    void relaxationAgarElementwise(int l, int i, int j, int k, int nx, int ny, int nz);
    void residueRestrictionAga(int l);
    void restrictionAgaElementwise(int l, int i, int j, int k, int nx, int ny, int nz, int nxp, int nyp, int nzp, double h2, int d);
    void interpolationAga(int l);
    void interpolationAgaElementwise(int l, int i, int j, int k, int nx, int ny, int nz, int nxp, int nyp, int nzp);
    void updateSolAga(int l);
    void updateResAga(int l);
    void NestSolAga(int l);
    void NestSolAgaElementwise(int l, int i, int j, int k, int nx, int ny, int nz, int nxp, int nyp, int nzp, double h2);
    void updateSolAgaElementwise(int l, int i, int j, int k, int nx, int ny, int nz, int nxp, int nyp, int nzp, double h2);
    void updateResAgaElementwise(int l, int i, int j, int k, int nx, int ny, int nz, int nxp, int nyp, int nzp, double h2);
    void rewriteSolBdryAga(int l);
    void setup();
    void Vcycle(int mu1, int mu2);
    void solve(Array3D<LocalAga>** CurrentAgar, Array2D<LocalAga>** CurrentWall, Array3D<LocalAga>** PreviousAgar, Array2D<LocalAga>** PreviousWall, int mu1, int mu2, double Dt);
    void reset();
    void resetSolution();
    double evalUpdateErr();
    double evalUpdateRes();
    AgaArray3D** Aga;
    AgaArray3D** AgaRes;
    void copyToMG(Array3D<LocalAga>** CurrentAgar, Array2D<LocalAga>** CurrentWall);
    void copyFromMG(Array3D<LocalAga>** CurrentAgar, Array2D<LocalAga>** CurrentWall, Array3D<LocalAga>** PreviousAgar, Array2D<LocalAga>** PreviousWall);
    void setInsideColony(Array3D<double>* insideColony){insideColonyDen = insideColony;};
    int getLevels() {return agarLevels;};
private:
    // private variables.
    AgaArray3D** AgaErr;
    Array3D<double>* insideColonyDen;
    double h;
    int agarLevels;
    int agarMGLevels;
    int BoxX;
    int BoxY;
    int BoxZAgar;
    double dt;
    
    // private member functions.
    int ipow(int base, int exp);
    double insideColonyAtLevel(int l, int i, int j, int k){return insideColonyDen[0](ipow(2,l)*i, ipow(2,l)*j, ipow(2,l)*k);};
};


int UpdateEnvArrayMG_OxygenAcetate(Array3D<LocalEnv>* CurrentColony, Array3D<LocalEnv>* PreviousColony, Array3D<LocalAga>** CurrentAgar, Array3D<LocalAga>** PreviousAgar, Array2D<LocalAga>** CurrentWall, Array2D<LocalAga>** PreviousWall, Array3D<double>& DensityShiftP, Array3D<double>& Density1ShiftP, Array3D<double>& Density2ShiftP, Array2D<double>& WallDensityShiftP, Array2D<double>& WallDensity1ShiftP, Array2D<double>& WallDensity2ShiftP, int minX, int maxX, int minY, int maxY, int maxH, Array2D<double>& Height, Array3D<double>& insideColonyDen, Multigrid& mg);

#endif
