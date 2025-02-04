#include <float.h>
#include "Multigrid.h"
#include <omp.h>

Multigrid::Multigrid(int levelsAga, int levelsAgarMG, int nx, int ny, int nzAgar)
{
    agarLevels = levelsAga;
    agarMGLevels = levelsAgarMG;
    BoxX = nx;
    BoxY = ny;
    BoxZAgar = nzAgar;
    
    Aga = new AgaArray3D*[levelsAgarMG];
    AgaRes = new AgaArray3D*[levelsAgarMG];
    AgaErr = new AgaArray3D*[levelsAgarMG];
    nx = BoxX;
    ny = BoxY;
    for (int i=0; i<levelsAga; i++)
    {
        Aga[i] = new AgaArray3D(nx+1,ny,nzAgar+1);
        AgaRes[i] = new AgaArray3D(nx+1,ny,nzAgar+1);
        AgaErr[i] = new AgaArray3D(nx+1,ny,nzAgar+1);
    }
    for (int i=levelsAga; i<levelsAgarMG; i++)
    {
        nx = nx/2;
        nzAgar = nzAgar/2;
        Aga[i] = new AgaArray3D(nx+1,ny,nzAgar+1);
        AgaRes[i] = new AgaArray3D(nx+1,ny,nzAgar+1);
        AgaErr[i] = new AgaArray3D(nx+1,ny,nzAgar+1);
    }
    
};

Multigrid::~Multigrid()
{
    for (int i=0; i<agarMGLevels; i++)
    {
        delete Aga[i];
        delete AgaRes[i];
        delete AgaErr[i];
    }
    delete []Aga;
    delete []AgaRes;
    delete []AgaErr;
    
}

void Multigrid::resetZero()
{
    for (int l=0;l<agarMGLevels;l++)
    {
        int d = ipow(2,l);
        int d0 = ipow(2,agarLevels);
        int nx = min(BoxX, BoxX*d0/d/2);
        int nz = min(BoxZAgar, BoxZAgar*d0/d/2);
        
        int i, k;
#pragma omp parallel for default(shared) private(i,k)
        for (i=0; i<=nx; i++)
        {
            for (k=0; k<=nz; k++)
            {
                AgaErr[l][0](i,0,k).CarbonAgar = 0;
                AgaErr[l][0](i,0,k).OxygenAgar = 0;
                AgaErr[l][0](i,0,k).AcetateAgar = 0;
            }
        }
        if (l>0)
        {
#pragma omp parallel for default(shared) private(i,k)
            for (i=0; i<=nx; i++)
            {
                for (k=0; k<=nz; k++)
                {
                    AgaRes[l][0](i,0,k).CarbonAgar = 0;
                    AgaRes[l][0](i,0,k).OxygenAgar = 0;
                    AgaRes[l][0](i,0,k).AcetateAgar = 0;
                }
            }
        }
    }
}

void Multigrid::resetSolution()
{
    for (int l=0;l<agarMGLevels;l++)
    {
        int d = ipow(2,l);
        int d0 = ipow(2,agarLevels);
        int nx = min(BoxX, BoxX*d0/d/2);
        int nz = min(BoxZAgar, BoxZAgar*d0/d/2);
        
        int i, k;
#pragma omp parallel for default(shared) private(i,k)
        for (i=0; i<=nx; i++)
        {
            for (k=0; k<=nz; k++)
            {
                AgaErr[l][0](i,0,k).CarbonAgar = 0;
                AgaErr[l][0](i,0,k).OxygenAgar = 0;
                AgaErr[l][0](i,0,k).AcetateAgar = 0;
            }
        }
        if (l>0)
        {
#pragma omp parallel for default(shared) private(i,k)
            for (i=0; i<=nx; i++)
            {
                for (k=0; k<=nz; k++)
                {
                    AgaRes[l][0](i,0,k).CarbonAgar = 0;
                    AgaRes[l][0](i,0,k).OxygenAgar = 0;
                    AgaRes[l][0](i,0,k).AcetateAgar = 0;
                }
            }
        }
    }
    
}

int Multigrid::ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp % 2)
            result *= base;
        exp /= 2;
        base *= base;
    }
    return result;
}

void Multigrid::relaxationAgarAtLevel(int l)
{
    int d = ipow(2,l);
    int d0 = ipow(2,agarLevels);
    int nx = min(BoxX, BoxX*d0/d/2);
    int ny = BoxY;
    int nz = min(BoxZAgar, BoxZAgar*d0/d/2);
    
    int i, k;
#pragma omp parallel for default(shared) private(i,k)
    for (i=0; i<nx+1; i+=2)
        for (k=0; k<nz+1; k++)
        {
            relaxationAgarElementwise(l, i, 0, k, nx, ny, nz);
        }
    
#pragma omp parallel for default(shared) private(i,k)
    for (i=1; i<nx+1; i+=2)
        for (k=0; k<nz+1; k++)
        {
            relaxationAgarElementwise(l, i, 0, k, nx, ny, nz);
        }
    
}

void Multigrid::relaxationAgarElementwise(int l, int i, int j, int k, int nx, int ny, int nz)
{
    double h2 = BoxLength*BoxLength*ipow(4,l);
    int d = ipow(2,l);
    if (k==0)
    {
        int i0 = min(max((i-nx/2)*d+BoxX/2,0),BoxX-1);
        if (insideColonyDen[0](i0,0,0)>0)
        {
            AgaErr[l][0](i,0,0).CarbonAgar = AgaRes[l][0](i,0,0).CarbonAgar;
            AgaErr[l][0](i,0,0).OxygenAgar = AgaRes[l][0](i,0,0).OxygenAgar;
            AgaErr[l][0](i,0,0).AcetateAgar = AgaRes[l][0](i,0,0).AcetateAgar;
        }
        else if ((i==0 || i==nx))
        {
            if (l>=agarLevels-1)
            {
                if (i==0)
                {
                    AgaErr[l][0](i,0,0).CarbonAgar = (AgaRes[l][0](i,0,0).CarbonAgar*h2/dt/DiffAgarCarbon+(AgaErr[l][0](i+1,0,0).CarbonAgar+AgaErr[l][0](i,0,1).CarbonAgar))/(h2/dt/DiffAgarCarbon+2);
                    AgaErr[l][0](i,0,0).OxygenAgar = AgaRes[l][0](i,0,0).OxygenAgar;
                    AgaErr[l][0](i,0,0).AcetateAgar = (AgaRes[l][0](i,0,0).AcetateAgar*h2/dt/DiffAgarAce+(AgaErr[l][0](i+1,0,0).AcetateAgar+AgaErr[l][0](i,0,1).AcetateAgar))/(h2/dt/DiffAgarAce+2);
                }
                else
                {
                    AgaErr[l][0](i,0,0).CarbonAgar = (AgaRes[l][0](i,0,0).CarbonAgar*h2/dt/DiffAgarCarbon+(AgaErr[l][0](i-1,0,0).CarbonAgar+AgaErr[l][0](i,0,1).CarbonAgar))/(h2/dt/DiffAgarCarbon+2);
                    AgaErr[l][0](i,0,0).OxygenAgar = AgaRes[l][0](i,0,0).OxygenAgar;
                    AgaErr[l][0](i,0,0).AcetateAgar = (AgaRes[l][0](i,0,0).AcetateAgar*h2/dt/DiffAgarAce+(AgaErr[l][0](i-1,0,0).AcetateAgar+AgaErr[l][0](i,0,1).AcetateAgar))/(h2/dt/DiffAgarAce+2);
                }
            }
        }
        else
        {
            AgaErr[l][0](i,0,0).CarbonAgar = (AgaRes[l][0](i,0,0).CarbonAgar*h2/dt/DiffAgarCarbon + AgaErr[l][0](i,0,1).CarbonAgar + 0.5*AgaErr[l][0](i-1,0,0).CarbonAgar + 0.5*AgaErr[l][0](i+1,0,0).CarbonAgar)/(h2/dt/DiffAgarCarbon+2);
            AgaErr[l][0](i,0,0).OxygenAgar = AgaRes[l][0](i,0,0).OxygenAgar;
            AgaErr[l][0](i,0,0).AcetateAgar = (AgaRes[l][0](i,0,0).AcetateAgar*h2/dt/DiffAgarAce + AgaErr[l][0](i,0,1).AcetateAgar + 0.5*AgaErr[l][0](i-1,0,0).AcetateAgar + 0.5*AgaErr[l][0](i+1,0,0).AcetateAgar)/(h2/dt/DiffAgarAce+2);
        }
    }
    else if (k==nz)
    {
        if (l>=agarLevels-1)
        {
            if (i==0 || i==nx)
            {
                if (i==0)
                {
                    AgaErr[l][0](i,0,k).CarbonAgar = (AgaRes[l][0](i,0,k).CarbonAgar*h2/dt/DiffAgarCarbon+(AgaErr[l][0](i+1,0,k).CarbonAgar+AgaErr[l][0](i,0,k-1).CarbonAgar))/(h2/dt/DiffAgarCarbon+2);
                    AgaErr[l][0](i,0,k).OxygenAgar = AgaRes[l][0](i,0,k).OxygenAgar;
                    AgaErr[l][0](i,0,k).AcetateAgar = (AgaRes[l][0](i,0,k).AcetateAgar*h2/dt/DiffAgarAce+(AgaErr[l][0](i+1,0,k).AcetateAgar+AgaErr[l][0](i,0,k-1).AcetateAgar))/(h2/dt/DiffAgarAce+2);
                }
                else
                {
                    AgaErr[l][0](i,0,k).CarbonAgar = (AgaRes[l][0](i,0,k).CarbonAgar*h2/dt/DiffAgarCarbon+(AgaErr[l][0](i-1,0,k).CarbonAgar+AgaErr[l][0](i,0,k-1).CarbonAgar))/(h2/dt/DiffAgarCarbon+2);
                    AgaErr[l][0](i,0,k).OxygenAgar = AgaRes[l][0](i,0,k).OxygenAgar;
                    AgaErr[l][0](i,0,k).AcetateAgar = (AgaRes[l][0](i,0,k).AcetateAgar*h2/dt/DiffAgarAce+(AgaErr[l][0](i-1,0,k).AcetateAgar+AgaErr[l][0](i,0,k-1).AcetateAgar))/(h2/dt/DiffAgarAce+2);
                }
            }
            else
            {
                AgaErr[l][0](i,0,k).CarbonAgar = (AgaRes[l][0](i,0,k).CarbonAgar*h2/dt/DiffAgarCarbon + AgaErr[l][0](i,0,k-1).CarbonAgar + 0.5*AgaErr[l][0](i+1,0,k).CarbonAgar + 0.5*AgaErr[l][0](i-1,0,k).CarbonAgar)/(h2/dt/DiffAgarCarbon+2);
                AgaErr[l][0](i,0,k).OxygenAgar = (AgaRes[l][0](i,0,k).OxygenAgar*h2/dt/DiffAgarO2 + AgaErr[l][0](i,0,k-1).OxygenAgar + 0.5*AgaErr[l][0](i+1,0,k).OxygenAgar + 0.5*AgaErr[l][0](i-1,0,k).OxygenAgar)/(h2/dt/DiffAgarO2+2);
                AgaErr[l][0](i,0,k).AcetateAgar = (AgaRes[l][0](i,0,k).AcetateAgar*h2/dt/DiffAgarAce + AgaErr[l][0](i,0,k-1).AcetateAgar + 0.5*AgaErr[l][0](i+1,0,k).AcetateAgar + 0.5*AgaErr[l][0](i-1,0,k).AcetateAgar)/(h2/dt/DiffAgarAce+2);
            }
        }
    }
    else if (i==0 || i==nx)
    {
        if (l>=agarLevels-1)
        {
            if (i==0)
            {
	        AgaErr[l][0](i,0,k).CarbonAgar = (AgaRes[l][0](i,0,k).CarbonAgar*h2/dt/DiffAgarCarbon + 0.5*AgaErr[l][0](i,0,k-1).CarbonAgar + 0.5*AgaErr[l][0](i,0,k+1).CarbonAgar + AgaErr[l][0](i+1,0,k).CarbonAgar)/(h2/dt/DiffAgarCarbon+2);
                AgaErr[l][0](i,0,k).OxygenAgar = AgaRes[l][0](i,0,k).OxygenAgar;
                AgaErr[l][0](i,0,k).AcetateAgar = (AgaRes[l][0](i,0,k).AcetateAgar*h2/dt/DiffAgarAce + 0.5*AgaErr[l][0](i,0,k-1).AcetateAgar + 0.5*AgaErr[l][0](i,0,k+1).AcetateAgar + AgaErr[l][0](i+1,0,k).AcetateAgar)/(h2/dt/DiffAgarAce+2);
            }
            else
            {
	        AgaErr[l][0](i,0,k).CarbonAgar = (AgaRes[l][0](i,0,k).CarbonAgar*h2/dt/DiffAgarCarbon + 0.5*AgaErr[l][0](i,0,k-1).CarbonAgar + 0.5*AgaErr[l][0](i,0,k+1).CarbonAgar + AgaErr[l][0](i-1,0,k).CarbonAgar)/(h2/dt/DiffAgarCarbon+2);
                AgaErr[l][0](i,0,k).OxygenAgar = AgaRes[l][0](i,0,k).OxygenAgar;
                AgaErr[l][0](i,0,k).AcetateAgar = (AgaRes[l][0](i,0,k).AcetateAgar*h2/dt/DiffAgarAce + 0.5*AgaErr[l][0](i,0,k-1).AcetateAgar + 0.5*AgaErr[l][0](i,0,k+1).AcetateAgar + AgaErr[l][0](i-1,0,k).AcetateAgar)/(h2/dt/DiffAgarAce+2);
            }
        }
    }
    else
    {
        AgaErr[l][0](i,0,k).CarbonAgar = (AgaRes[l][0](i,0,k).CarbonAgar*h2/dt/DiffAgarCarbon + 0.5*AgaErr[l][0](i,0,k-1).CarbonAgar + 0.5*AgaErr[l][0](i,0,k+1).CarbonAgar + 0.5*AgaErr[l][0](i-1,0,k).CarbonAgar + 0.5*AgaErr[l][0](i+1,0,k).CarbonAgar)/(h2/dt/DiffAgarCarbon+2);
        AgaErr[l][0](i,0,k).OxygenAgar = (AgaRes[l][0](i,0,k).OxygenAgar*h2/dt/DiffAgarO2 + 0.5*AgaErr[l][0](i,0,k-1).OxygenAgar + 0.5*AgaErr[l][0](i,0,k+1).OxygenAgar + 0.5*AgaErr[l][0](i-1,0,k).OxygenAgar + 0.5*AgaErr[l][0](i+1,0,k).OxygenAgar)/(h2/dt/DiffAgarO2+2);
        AgaErr[l][0](i,0,k).AcetateAgar = (AgaRes[l][0](i,0,k).AcetateAgar*h2/dt/DiffAgarAce + 0.5*AgaErr[l][0](i,0,k-1).AcetateAgar + 0.5*AgaErr[l][0](i,0,k+1).AcetateAgar + 0.5*AgaErr[l][0](i-1,0,k).AcetateAgar + 0.5*AgaErr[l][0](i+1,0,k).AcetateAgar)/(h2/dt/DiffAgarAce+2);
    }
}

void Multigrid::residueRestrictionAga(int l)
{
    int d = ipow(2,l);
    int d0 = ipow(2,agarLevels);
    int nx = min(BoxX, BoxX*d0/d/2);
    int nz = min(BoxZAgar, BoxZAgar*d0/d/2);
    int nxp = min(BoxX, BoxX*d0/d);
    int nzp = min(BoxZAgar, BoxZAgar*d0/d);
    double h2 = BoxLength*BoxLength*d*d;
    int i, k;
#pragma omp parallel for default(shared) private(i,k)
    for ( i=0; i<nx+1; i++)
    {
        for ( k=0;k<nz+1;k++)
        {
            restrictionAgaElementwise(l,i,0,k,nx,BoxY,nz,nxp,BoxY,nzp,h2,d);
        }
    }
}

void Multigrid::restrictionAgaElementwise(int l, int i, int j, int k, int nx, int ny, int nz, int nxp, int nyp, int nzp, double h2, int d)
{
    int i0 = min(max((i-nx/2)*d+BoxX/2,0),BoxX);
    int ip = 2*i-nx+nxp/2;
    int kp = 2*k;
    if (i==0 || i==nx || k==nz)
    {
        if (l>=agarLevels)
        {
            if (k==nz && (i>0) && (i<nx))
            {
                AgaRes[l][0](i,0,k).CarbonAgar = (AgaRes[l-1][0](ip,0,kp).CarbonAgar+0.5*(AgaErr[l-1][0](ip+1,0,kp).CarbonAgar+AgaErr[l-1][0](ip-1,0,kp).CarbonAgar+2*AgaErr[l-1][0](ip,0,kp-1).CarbonAgar-4*AgaErr[l-1][0](ip,0,kp).CarbonAgar)*DiffAgarCarbon*dt*4/h2-AgaErr[l-1][0](ip,0,kp).CarbonAgar)/2;
                AgaRes[l][0](i,0,k).OxygenAgar = (AgaRes[l-1][0](ip,0,kp).OxygenAgar+0.5*(AgaErr[l-1][0](2*i+1,0,2*k).OxygenAgar+AgaErr[l-1][0](2*i-1,0,2*k).OxygenAgar+2*AgaErr[l-1][0](2*i,0,2*k-1).OxygenAgar-4*AgaErr[l-1][0](2*i,0,2*k).OxygenAgar)*DiffAgarO2*dt*4/h2-AgaErr[l-1][0](2*i,0,2*k).OxygenAgar)/2;
                AgaRes[l][0](i,0,k).AcetateAgar = (AgaRes[l-1][0](ip,0,kp).AcetateAgar+0.5*(AgaErr[l-1][0](2*i+1,0,2*k).AcetateAgar+AgaErr[l-1][0](2*i-1,0,2*k).AcetateAgar+2*AgaErr[l-1][0](2*i,0,2*k-1).AcetateAgar-4*AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)*DiffAgarAce*dt*4/h2-AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)/2;
            }
            else
            {
                if (i==0)
                {
                    if (k==0)
                    {
                        AgaRes[l][0](i,0,k).CarbonAgar = (AgaRes[l-1][0](ip,0,kp).CarbonAgar+0.5*(2*AgaErr[l-1][0](ip+1,0,kp).CarbonAgar+2*AgaErr[l-1][0](ip,0,kp+1).CarbonAgar-4*AgaErr[l-1][0](ip,0,kp).CarbonAgar)*DiffAgarCarbon*dt*4/h2-AgaErr[l-1][0](ip,0,kp).CarbonAgar)/2;
                        AgaRes[l][0](i,0,k).OxygenAgar = AgaRes[l-1][0](ip,0,kp).OxygenAgar-AgaErr[l-1][0](ip,0,kp).OxygenAgar;
                        AgaRes[l][0](i,0,k).AcetateAgar = (AgaRes[l-1][0](ip,0,kp).AcetateAgar+0.5*(2*AgaErr[l-1][0](2*i+1,0,2*k).AcetateAgar+2*AgaErr[l-1][0](2*i,0,2*k+1).AcetateAgar-4*AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)*DiffAgarAce*dt*4/h2-AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)/2;
                    }
                    else if (k==nz)
                    {
                        AgaRes[l][0](i,0,k).CarbonAgar = (AgaRes[l-1][0](ip,0,kp).CarbonAgar+0.5*(2*AgaErr[l-1][0](ip+1,0,kp).CarbonAgar+2*AgaErr[l-1][0](ip,0,kp-1).CarbonAgar-4*AgaErr[l-1][0](ip,0,kp).CarbonAgar)*DiffAgarCarbon*dt*4/h2-AgaErr[l-1][0](ip,0,kp).CarbonAgar)/2;
                        AgaRes[l][0](i,0,k).OxygenAgar = AgaRes[l-1][0](ip,0,kp).OxygenAgar-AgaErr[l-1][0](ip,0,kp).OxygenAgar;
                        AgaRes[l][0](i,0,k).AcetateAgar = (AgaRes[l-1][0](ip,0,kp).AcetateAgar+0.5*(2*AgaErr[l-1][0](2*i+1,0,2*k).AcetateAgar+2*AgaErr[l-1][0](2*i,0,2*k-1).AcetateAgar-4*AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)*DiffAgarAce*dt*4/h2-AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)/2;
                    }
                    else
                    {
                        AgaRes[l][0](i,0,k).CarbonAgar = (AgaRes[l-1][0](ip,0,kp).CarbonAgar+0.5*(2*AgaErr[l-1][0](ip+1,0,kp).CarbonAgar+AgaErr[l-1][0](ip,0,kp+1).CarbonAgar+AgaErr[l-1][0](ip,0,kp-1).CarbonAgar-4*AgaErr[l-1][0](ip,0,kp).CarbonAgar)*DiffAgarCarbon*dt*4/h2-AgaErr[l-1][0](ip,0,kp).CarbonAgar)/2;
                        AgaRes[l][0](i,0,k).OxygenAgar = AgaRes[l-1][0](ip,0,kp).OxygenAgar-AgaErr[l-1][0](ip,0,kp).OxygenAgar;
                        AgaRes[l][0](i,0,k).AcetateAgar = (AgaRes[l-1][0](ip,0,kp).AcetateAgar+0.5*(2*AgaErr[l-1][0](2*i+1,0,2*k).AcetateAgar+AgaErr[l-1][0](2*i,0,2*k+1).AcetateAgar+AgaErr[l-1][0](2*i,0,2*k-1).AcetateAgar-4*AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)*DiffAgarAce*dt*4/h2-AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)/2;
                    }
                }
                else
                {
                    if (k==0)
                    {
                        AgaRes[l][0](i,0,k).CarbonAgar = (AgaRes[l-1][0](ip,0,kp).CarbonAgar+0.5*(2*AgaErr[l-1][0](ip-1,0,kp).CarbonAgar+2*AgaErr[l-1][0](ip,0,kp+1).CarbonAgar-4*AgaErr[l-1][0](ip,0,kp).CarbonAgar)*DiffAgarCarbon*dt*4/h2-AgaErr[l-1][0](ip,0,kp).CarbonAgar)/2;
                        AgaRes[l][0](i,0,k).OxygenAgar = AgaRes[l-1][0](ip,0,kp).OxygenAgar-AgaErr[l-1][0](ip,0,kp).OxygenAgar;
                        AgaRes[l][0](i,0,k).AcetateAgar = (AgaRes[l-1][0](ip,0,kp).AcetateAgar+0.5*(2*AgaErr[l-1][0](2*i-1,0,2*k).AcetateAgar+2*AgaErr[l-1][0](2*i,0,2*k+1).AcetateAgar-4*AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)*DiffAgarAce*dt*4/h2-AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)/2;
                    }
                    else if (k==nz)
                    {
                        AgaRes[l][0](i,0,k).CarbonAgar = (AgaRes[l-1][0](ip,0,kp).CarbonAgar+0.5*(2*AgaErr[l-1][0](ip-1,0,kp).CarbonAgar+2*AgaErr[l-1][0](ip,0,kp-1).CarbonAgar-4*AgaErr[l-1][0](ip,0,kp).CarbonAgar)*DiffAgarCarbon*dt*4/h2-AgaErr[l-1][0](ip,0,kp).CarbonAgar)/2;
                        AgaRes[l][0](i,0,k).OxygenAgar = AgaRes[l-1][0](ip,0,kp).OxygenAgar-AgaErr[l-1][0](ip,0,kp).OxygenAgar;
                        AgaRes[l][0](i,0,k).AcetateAgar = (AgaRes[l-1][0](ip,0,kp).AcetateAgar+0.5*(2*AgaErr[l-1][0](2*i-1,0,2*k).AcetateAgar+2*AgaErr[l-1][0](2*i,0,2*k-1).AcetateAgar-4*AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)*DiffAgarAce*dt*4/h2-AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)/2;
                    }
                    else
                    {
                        AgaRes[l][0](i,0,k).CarbonAgar = (AgaRes[l-1][0](ip,0,kp).CarbonAgar+0.5*(2*AgaErr[l-1][0](ip-1,0,kp).CarbonAgar+AgaErr[l-1][0](ip,0,kp+1).CarbonAgar+AgaErr[l-1][0](ip,0,kp-1).CarbonAgar-4*AgaErr[l-1][0](ip,0,kp).CarbonAgar)*DiffAgarCarbon*dt*4/h2-AgaErr[l-1][0](ip,0,kp).CarbonAgar)/2;
                        AgaRes[l][0](i,0,k).OxygenAgar = AgaRes[l-1][0](ip,0,kp).OxygenAgar-AgaErr[l-1][0](ip,0,kp).OxygenAgar;
                        AgaRes[l][0](i,0,k).AcetateAgar = (AgaRes[l-1][0](ip,0,kp).AcetateAgar+0.5*(2*AgaErr[l-1][0](2*i-1,0,2*k).AcetateAgar+AgaErr[l-1][0](2*i,0,2*k+1).AcetateAgar+AgaErr[l-1][0](2*i,0,2*k-1).AcetateAgar-4*AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)*DiffAgarAce*dt*4/h2-AgaErr[l-1][0](2*i,0,2*k).AcetateAgar)/2;
                    }
                }
            }
        }
    }
    else if (k==0)
    {
        if (insideColonyDen[0](i0,0,0)>0)
        {
            AgaRes[l][0](i,0,0).CarbonAgar = 0;
            AgaRes[l][0](i,0,0).OxygenAgar = 0;
            AgaRes[l][0](i,0,0).AcetateAgar = 0;
        }
        else if (ip>0 && ip<nxp)
        {
            AgaRes[l][0](i,0,0).CarbonAgar = (AgaRes[l-1][0](ip,0,0).CarbonAgar+0.5*(AgaErr[l-1][0](ip+1,0,0).CarbonAgar+AgaErr[l-1][0](ip-1,0,0).CarbonAgar+2*AgaErr[l-1][0](ip,0,1).CarbonAgar-4*AgaErr[l-1][0](ip,0,0).CarbonAgar)*DiffAgarCarbon*dt*4/h2-AgaErr[l-1][0](ip,0,0).CarbonAgar)*0.5;
            AgaRes[l][0](i,0,0).AcetateAgar = (AgaRes[l-1][0](ip,0,0).AcetateAgar+0.5*(AgaErr[l-1][0](ip+1,0,0).AcetateAgar+AgaErr[l-1][0](ip-1,0,0).AcetateAgar+2*AgaErr[l-1][0](ip,0,1).AcetateAgar-4*AgaErr[l-1][0](ip,0,0).AcetateAgar)*DiffAgarAce*dt*4/h2-AgaErr[l-1][0](ip,0,0).AcetateAgar)*0.5;
            AgaRes[l][0](i,0,0).OxygenAgar = (AgaRes[l-1][0](ip,0,0).OxygenAgar-AgaErr[l-1][0](ip,0,0).OxygenAgar)*0.5;
        }
        else if ((l<agarLevels) && (ip>=0 && ip<=nxp) && ((ip==0) || (ip==nxp)))
        {
            if (ip==0)
            {
                AgaRes[l][0](i,0,0).CarbonAgar = (Aga[l-1][0](ip+2,0,0).CarbonAgar+0.5*AgaErr[l-1][0](ip+2,0,0).CarbonAgar+Aga[l][0](i-1,0,0).CarbonAgar+2*Aga[l][0](i,0,1).CarbonAgar-4*Aga[l][0](i,0,0).CarbonAgar)*DiffAgarCarbon*dt/h2;
                AgaRes[l][0](i,0,0).AcetateAgar = (Aga[l-1][0](ip+2,0,0).AcetateAgar+0.5*AgaErr[l-1][0](ip+2,0,0).AcetateAgar+Aga[l][0](i-1,0,0).AcetateAgar+2*Aga[l][0](i,0,1).AcetateAgar-4*Aga[l][0](i,0,0).AcetateAgar)*DiffAgarAce*dt/h2;
                AgaRes[l][0](i,0,0).OxygenAgar = maxO2_Agar-Aga[l][0](i,0,0).OxygenAgar;
            }
            else if (ip==nxp)
            {
                AgaRes[l][0](i,0,0).CarbonAgar = (Aga[l][0](i+1,0,0).CarbonAgar+0.5*AgaErr[l-1][0](ip-2,0,0).CarbonAgar+Aga[l-1][0](ip-2,0,0).CarbonAgar+2*Aga[l][0](i,0,1).CarbonAgar-4*Aga[l][0](i,0,0).CarbonAgar)*dt*DiffAgarCarbon/h2;
                AgaRes[l][0](i,0,0).AcetateAgar = (Aga[l][0](i+1,0,0).AcetateAgar+0.5*AgaErr[l-1][0](ip-2,0,0).AcetateAgar+Aga[l-1][0](ip-2,0,0).AcetateAgar+2*Aga[l][0](i,0,1).AcetateAgar-4*Aga[l][0](i,0,0).AcetateAgar)*dt*DiffAgarAce/h2;
                AgaRes[l][0](i,0,0).OxygenAgar = maxO2_Agar-Aga[l][0](i,0,0).OxygenAgar;
            }
        }
    }
    else if (l<agarLevels && ip>=0 && ip<=nxp && kp<=nzp) // continue to work from here.
    {
        if (kp==nzp)
        {
            if (ip==0 || ip==nxp)
            {
                AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l][0](i+1,0,k).CarbonAgar+Aga[l][0](i-1,0,k).CarbonAgar+Aga[l][0](i,0,k+1).CarbonAgar+Aga[l][0](i,0,k-1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*DiffAgarCarbon*dt/h2;
                AgaRes[l][0](i,0,k).OxygenAgar = (Aga[l][0](i+1,0,k).OxygenAgar+Aga[l][0](i-1,0,k).OxygenAgar+Aga[l][0](i,0,k+1).OxygenAgar+Aga[l][0](i,0,k-1).OxygenAgar-4*Aga[l][0](i,0,k).OxygenAgar)*DiffAgarO2*dt/h2;
                AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l][0](i+1,0,k).AcetateAgar+Aga[l][0](i-1,0,k).AcetateAgar+Aga[l][0](i,0,k+1).AcetateAgar+Aga[l][0](i,0,k-1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*DiffAgarAce*dt/h2;
            }
            else
            {
                AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l][0](i+1,0,k).CarbonAgar+Aga[l][0](i-1,0,k).CarbonAgar+Aga[l][0](i,0,k+1).CarbonAgar+Aga[l-1][0](ip,0,kp-2).CarbonAgar+0.5*AgaErr[l-1][0](ip,0,kp-2).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*DiffAgarCarbon*dt/h2;
                AgaRes[l][0](i,0,k).OxygenAgar = (Aga[l][0](i+1,0,k).OxygenAgar+Aga[l][0](i-1,0,k).OxygenAgar+Aga[l][0](i,0,k+1).OxygenAgar+Aga[l-1][0](ip,0,kp-2).OxygenAgar+0.5*AgaErr[l-1][0](ip,0,kp-2).OxygenAgar-4*Aga[l][0](i,0,k).OxygenAgar)*DiffAgarO2*dt/h2;
                AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l][0](i+1,0,k).AcetateAgar+Aga[l][0](i-1,0,k).AcetateAgar+Aga[l][0](i,0,k+1).AcetateAgar+Aga[l-1][0](ip,0,kp-2).AcetateAgar+0.5*AgaErr[l-1][0](ip,0,kp-2).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*DiffAgarAce*dt/h2;
            }
        }
        else if (ip==0)
        {
            AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l-1][0](ip+2,0,kp).CarbonAgar+0.5*AgaErr[l-1][0](ip+2,0,kp).CarbonAgar+Aga[l][0](i-1,0,k).CarbonAgar+Aga[l][0](i,0,k+1).CarbonAgar+Aga[l][0](i,0,k-1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*DiffAgarCarbon*dt/h2;
            AgaRes[l][0](i,0,k).OxygenAgar = (Aga[l-1][0](ip+2,0,kp).OxygenAgar+0.5*AgaErr[l-1][0](ip+2,0,kp).OxygenAgar+Aga[l][0](i-1,0,k).OxygenAgar+Aga[l][0](i,0,k+1).OxygenAgar+Aga[l][0](i,0,k-1).OxygenAgar-4*Aga[l][0](i,0,k).OxygenAgar)*DiffAgarO2*dt/h2;
            AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l-1][0](ip+2,0,kp).AcetateAgar+0.5*AgaErr[l-1][0](ip+2,0,kp).AcetateAgar+Aga[l][0](i-1,0,k).AcetateAgar+Aga[l][0](i,0,k+1).AcetateAgar+Aga[l][0](i,0,k-1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*DiffAgarAce*dt/h2;
        }
        else if (ip==nxp)
        {
            AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l][0](i+1,0,k).CarbonAgar+0.5*AgaErr[l-1][0](ip-2,0,kp).CarbonAgar+Aga[l-1][0](ip-2,0,kp).CarbonAgar+Aga[l][0](i,0,k+1).CarbonAgar+Aga[l][0](i,0,k-1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*DiffAgarCarbon*dt/h2;
            AgaRes[l][0](i,0,k).OxygenAgar = (Aga[l][0](i+1,0,k).OxygenAgar+0.5*AgaErr[l-1][0](ip-2,0,kp).OxygenAgar+Aga[l-1][0](ip-2,0,kp).OxygenAgar+Aga[l][0](i,0,k+1).OxygenAgar+Aga[l][0](i,0,k-1).OxygenAgar-4*Aga[l][0](i,0,k).OxygenAgar)*DiffAgarO2*dt/h2;
            AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l][0](i+1,0,k).AcetateAgar+0.5*AgaErr[l-1][0](ip-2,0,kp).AcetateAgar+Aga[l-1][0](ip-2,0,kp).AcetateAgar+Aga[l][0](i,0,k+1).AcetateAgar+Aga[l][0](i,0,k-1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*DiffAgarAce*dt/h2;
        }
        else
        {
            AgaRes[l][0](i,0,k).CarbonAgar = (AgaRes[l-1][0](ip,0,kp).CarbonAgar+0.5*(AgaErr[l-1][0](ip+1,0,kp).CarbonAgar+AgaErr[l-1][0](ip-1,0,kp).CarbonAgar+AgaErr[l-1][0](ip,0,kp+1).CarbonAgar+AgaErr[l-1][0](ip,0,kp-1).CarbonAgar-4*AgaErr[l-1][0](ip,0,kp).CarbonAgar)*DiffAgarCarbon*dt/h2*4-AgaErr[l-1][0](ip,0,kp).CarbonAgar)*0.5;
            AgaRes[l][0](i,0,k).OxygenAgar = (AgaRes[l-1][0](ip,0,kp).OxygenAgar+0.5*(AgaErr[l-1][0](ip+1,0,kp).OxygenAgar+AgaErr[l-1][0](ip-1,0,kp).OxygenAgar+AgaErr[l-1][0](ip,0,kp+1).OxygenAgar+AgaErr[l-1][0](ip,0,kp-1).OxygenAgar-4*AgaErr[l-1][0](ip,0,kp).OxygenAgar)*DiffAgarO2*dt/h2*4-AgaErr[l-1][0](ip,0,kp).OxygenAgar)*0.5;
            AgaRes[l][0](i,0,k).AcetateAgar = (AgaRes[l-1][0](ip,0,kp).AcetateAgar+0.5*(AgaErr[l-1][0](ip+1,0,kp).AcetateAgar+AgaErr[l-1][0](ip-1,0,kp).AcetateAgar+AgaErr[l-1][0](ip,0,kp+1).AcetateAgar+AgaErr[l-1][0](ip,0,kp-1).AcetateAgar-4*AgaErr[l-1][0](ip,0,kp).AcetateAgar)*DiffAgarAce*dt/h2*4-AgaErr[l-1][0](ip,0,kp).AcetateAgar)*0.5;
        }
    }
    else if (l>=agarLevels)
    {
        AgaRes[l][0](i,0,k).CarbonAgar = (AgaRes[l-1][0](ip,0,kp).CarbonAgar+0.5*(AgaErr[l-1][0](ip+1,0,kp).CarbonAgar+AgaErr[l-1][0](ip-1,0,kp).CarbonAgar+AgaErr[l-1][0](ip,0,kp+1).CarbonAgar+AgaErr[l-1][0](ip,0,kp-1).CarbonAgar-4*AgaErr[l-1][0](ip,0,kp).CarbonAgar)*DiffAgarCarbon*dt/h2*4-AgaErr[l-1][0](ip,0,kp).CarbonAgar)*0.5;
        AgaRes[l][0](i,0,k).OxygenAgar = (AgaRes[l-1][0](ip,0,kp).OxygenAgar+0.5*(AgaErr[l-1][0](ip+1,0,kp).OxygenAgar+AgaErr[l-1][0](ip-1,0,kp).OxygenAgar+AgaErr[l-1][0](ip,0,kp+1).OxygenAgar+AgaErr[l-1][0](ip,0,kp-1).OxygenAgar-4*AgaErr[l-1][0](ip,0,kp).OxygenAgar)*DiffAgarO2*dt/h2*4-AgaErr[l-1][0](ip,0,kp).OxygenAgar)*0.5;
        AgaRes[l][0](i,0,k).AcetateAgar = (AgaRes[l-1][0](ip,0,kp).AcetateAgar+0.5*(AgaErr[l-1][0](ip+1,0,kp).AcetateAgar+AgaErr[l-1][0](ip-1,0,kp).AcetateAgar+AgaErr[l-1][0](ip,0,kp+1).AcetateAgar+AgaErr[l-1][0](ip,0,kp-1).AcetateAgar-4*AgaErr[l-1][0](ip,0,kp).AcetateAgar)*DiffAgarAce*dt/h2*4-AgaErr[l-1][0](ip,0,kp).AcetateAgar)*0.5;
    }
}

void Multigrid::interpolationAga(int l)
{
    int d = ipow(2,l);
    int d0 = ipow(2,agarLevels);
    int nx = min(BoxX, BoxX*d0/d/2);
    int nz = min(BoxZAgar, BoxZAgar*d0/d/2);
    int nxp = min(BoxX, BoxX*d0/d/4);
    int nzp = min(BoxZAgar, BoxZAgar*d0/d/4);
    int i, k;
#pragma omp parallel for default(shared) private(i,k)
    for ( i=0; i<nx+1; i++)
    {
        for ( k=0;k<nz+1;k++)
        {
                interpolationAgaElementwise(l,i,0,k,nx,BoxY,nz,nxp,BoxY,nzp);
        }
    }
}

void Multigrid::interpolationAgaElementwise(int l, int i, int j, int k, int nx, int ny, int nz, int nxp, int nyp, int nzp)
{
    int i0 = i%2;
    int k0 = k%2;
    
    if (i0==0)
    {
        int ip = nxp/2-nx/4+i/2;
        if (k0==0)
        {
            int kp = k/2;
            AgaErr[l][0](i,0,k).CarbonAgar += AgaErr[l+1][0](ip,0,kp).CarbonAgar;
            AgaErr[l][0](i,0,k).OxygenAgar += AgaErr[l+1][0](ip,0,kp).OxygenAgar;
            AgaErr[l][0](i,0,k).AcetateAgar += AgaErr[l+1][0](ip,0,kp).AcetateAgar;
        }
        else
        {
            int kp = (k-1)/2;
            AgaErr[l][0](i,0,k).CarbonAgar += 0.5*(AgaErr[l+1][0](ip,0,kp).CarbonAgar+AgaErr[l+1][0](ip,0,kp+1).CarbonAgar);
            AgaErr[l][0](i,0,k).OxygenAgar += 0.5*(AgaErr[l+1][0](ip,0,kp).OxygenAgar+AgaErr[l+1][0](ip,0,kp+1).OxygenAgar);
            AgaErr[l][0](i,0,k).AcetateAgar += 0.5*(AgaErr[l+1][0](ip,0,kp).AcetateAgar+AgaErr[l+1][0](ip,0,kp+1).AcetateAgar);
        }
    }
    else
    {
        int ip = nxp/2-nx/4+(i-1)/2;
        if (k0==0)
        {
            int kp = k/2;
            AgaErr[l][0](i,0,k).CarbonAgar += 0.5*(AgaErr[l+1][0](ip,0,kp).CarbonAgar+AgaErr[l+1][0](ip+1,0,kp).CarbonAgar);
            AgaErr[l][0](i,0,k).OxygenAgar += 0.5*(AgaErr[l+1][0](ip,0,kp).OxygenAgar+AgaErr[l+1][0](ip+1,0,kp).OxygenAgar);
            AgaErr[l][0](i,0,k).AcetateAgar += 0.5*(AgaErr[l+1][0](ip,0,kp).AcetateAgar+AgaErr[l+1][0](ip+1,0,kp).AcetateAgar);
        }
        else
        {
            int kp = (k-1)/2;
            AgaErr[l][0](i,0,k).CarbonAgar += 0.25*(AgaErr[l+1][0](ip,0,kp).CarbonAgar+AgaErr[l+1][0](ip,0,kp+1).CarbonAgar+AgaErr[l+1][0](ip+1,0,kp).CarbonAgar+AgaErr[l+1][0](ip+1,0,kp+1).CarbonAgar);
            AgaErr[l][0](i,0,k).OxygenAgar += 0.25*(AgaErr[l+1][0](ip,0,kp).OxygenAgar+AgaErr[l+1][0](ip,0,kp+1).OxygenAgar+AgaErr[l+1][0](ip+1,0,kp).OxygenAgar+AgaErr[l+1][0](ip+1,0,kp+1).OxygenAgar);
            AgaErr[l][0](i,0,k).AcetateAgar += 0.25*(AgaErr[l+1][0](ip,0,kp).AcetateAgar+AgaErr[l+1][0](ip,0,kp+1).AcetateAgar+AgaErr[l+1][0](ip+1,0,kp).AcetateAgar+AgaErr[l+1][0](ip+1,0,kp+1).AcetateAgar);
            }
    }
}

void Multigrid::rewriteSolBdryAga(int l)
{
    int d = ipow(2,l);
    int d0 = ipow(2,agarLevels);
    int nx = min(BoxX, BoxX*d0/d/2);
    int nz = min(BoxZAgar, BoxZAgar*d0/d/2);
    int nx2 = min(BoxX, BoxX*d0/d/4);
    int i, k, i2, k2, i2p1, k2p1;
    k  = nz;
    k2 = nz/2;
#pragma omp parallel for default(shared) private(i,i2, i2p1)
    for (i=0; i<nx+1; i++)
    {
        i2 = nx2/2+i/2-nx/4;
        i2p1 = nx2/2+(i+1)/2-nx/4;
        Aga[l][0](i,0,k).CarbonAgar = (Aga[l+1][0](i2,0,k2).CarbonAgar+Aga[l+1][0](i2p1,0,k2).CarbonAgar)/2;
        Aga[l][0](i,0,k).OxygenAgar = (Aga[l+1][0](i2,0,k2).OxygenAgar+Aga[l+1][0](i2p1,0,k2).OxygenAgar)/2;
        Aga[l][0](i,0,k).AcetateAgar = (Aga[l+1][0](i2,0,k2).AcetateAgar+Aga[l+1][0](i2p1,0,k2).AcetateAgar)/2;
    }
#pragma omp parallel for default(shared) private(i,k,i2,k2, k2p1)
    for (k=0; k<nz; k++)
    {
        i = 0;
        i2 = nx2/2+i/2-nx/4;
        k2 = k/2;
        k2p1 = (k+1)/2;
        Aga[l][0](i,0,k).CarbonAgar = (Aga[l+1][0](i2,0,k2).CarbonAgar+Aga[l+1][0](i2,0,k2p1).CarbonAgar)/2;
        Aga[l][0](i,0,k).OxygenAgar = (Aga[l+1][0](i2,0,k2).OxygenAgar+Aga[l+1][0](i2,0,k2p1).OxygenAgar)/2;
        Aga[l][0](i,0,k).AcetateAgar = (Aga[l+1][0](i2,0,k2).AcetateAgar+Aga[l+1][0](i2,0,k2p1).AcetateAgar)/2;
        i = nx;
        i2 = nx2/2+i/2-nx/4;
        Aga[l][0](i,0,k).CarbonAgar = (Aga[l+1][0](i2,0,k2).CarbonAgar+Aga[l+1][0](i2,0,k2p1).CarbonAgar)/2;
        Aga[l][0](i,0,k).OxygenAgar = (Aga[l+1][0](i2,0,k2).OxygenAgar+Aga[l+1][0](i2,0,k2p1).OxygenAgar)/2;
        Aga[l][0](i,0,k).AcetateAgar = (Aga[l+1][0](i2,0,k2).AcetateAgar+Aga[l+1][0](i2,0,k2p1).AcetateAgar)/2;
    }
}

void Multigrid::updateSolAga(int l)
{
    int d = ipow(2,l);
    int d0 = ipow(2,agarLevels);
    int nx = min(BoxX, BoxX*d0/d/2);
    int nz = min(BoxZAgar, BoxZAgar*d0/d/2);
    int nxp = min(BoxX, BoxX*d0/d);
    int nzp = min(BoxZAgar, BoxZAgar*d0/d);
    double h2 = BoxLength*BoxLength*d*d;
    int i, k;
    
#pragma omp parallel for default(shared) private(i,k)
    for ( i=0; i<nx+1; i++)
    {
        for ( k=0;k<nz+1;k++)
        {
            updateSolAgaElementwise(l,i,0,k,nx,BoxY,nz,nxp,BoxY,nzp,h2);
        }
    }
}

void Multigrid::updateResAga(int l)
{
    int d = ipow(2,l);
    int d0 = ipow(2,agarLevels);
    int nx = min(BoxX, BoxX*d0/d/2);
    int nz = min(BoxZAgar, BoxZAgar*d0/d/2);
    int nxp = min(BoxX, BoxX*d0/d);
    int nzp = min(BoxZAgar, BoxZAgar*d0/d);
    double h2 = BoxLength*BoxLength*d*d;
    int i, k;
#pragma omp parallel for default(shared) private(i,k)
    for ( i=0; i<nx+1; i++)
    {
        for ( k=0;k<nz+1;k++)
        {
            updateResAgaElementwise(l,i,0,k,nx,BoxY,nz,nxp,BoxY,nzp,h2);
        }
    }
}

void Multigrid::updateSolAgaElementwise(int l, int i, int j, int k, int nx, int ny, int nz, int nxp, int nyp, int nzp, double h2)
{
    AgaRes[l][0](i,0,k).CarbonAgar = 0;
    AgaRes[l][0](i,0,k).OxygenAgar = 0;
    AgaRes[l][0](i,0,k).AcetateAgar = 0;
    int ip = 2*i - nx+nxp/2;
    int kp = 2*k;
    if (l<agarLevels)
    {
        if (l==0 || ip<=0 || ip>=nxp || kp>=nzp)
        {
            Aga[l][0](i,0,k).CarbonAgar += AgaErr[l][0](i,0,k).CarbonAgar;
            Aga[l][0](i,0,k).OxygenAgar += AgaErr[l][0](i,0,k).OxygenAgar;
            Aga[l][0](i,0,k).AcetateAgar += AgaErr[l][0](i,0,k).AcetateAgar;
            Aga[l][0](i,0,k).CarbonAgar = max(min(Aga[l][0](i,0,k).CarbonAgar,maxCarbon),0.0);
            Aga[l][0](i,0,k).OxygenAgar = max(min(Aga[l][0](i,0,k).OxygenAgar,maxO2_Agar),0.0);
            Aga[l][0](i,0,k).AcetateAgar = max(Aga[l][0](i,0,k).AcetateAgar,0.0);
        }
    }
    AgaErr[l][0](i,0,k).CarbonAgar = 0;
    AgaErr[l][0](i,0,k).OxygenAgar = 0;
    AgaErr[l][0](i,0,k).AcetateAgar = 0;
}

void Multigrid::updateResAgaElementwise(int l, int i, int j, int k, int nx, int ny, int nz, int nxp, int nyp, int nzp, double h2)
{
    int ip = 2*i - nx+nxp/2;
    int kp = 2*k;

    if (l<agarLevels)
    {
        if (i==0 || i==nx)
        {
            if (l==agarLevels-1)
            {
                if (i==0)
                {
                    if (k==0)
                    {
                        AgaRes[l][0](i,0,k).CarbonAgar = (2*Aga[l][0](i+1,0,k).CarbonAgar+2*Aga[l][0](i,0,k+1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*DiffAgarCarbon*dt/h2;
                        AgaRes[l][0](i,0,k).OxygenAgar = 0;//maxO2-Aga[l][0](i,0,k).OxygenAgar;
                        AgaRes[l][0](i,0,k).AcetateAgar = (2*Aga[l][0](i+1,0,k).AcetateAgar+2*Aga[l][0](i,0,k+1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*DiffAgarAce*dt/h2;
                    }
                    else if (k==nz)
                    {
                        AgaRes[l][0](i,0,k).CarbonAgar = (2*Aga[l][0](i+1,0,k).CarbonAgar+2*Aga[l][0](i,0,k-1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*DiffAgarCarbon*dt/h2;
                        AgaRes[l][0](i,0,k).OxygenAgar = 0;//maxO2-Aga[l][0](i,0,k).OxygenAgar;
                        AgaRes[l][0](i,0,k).AcetateAgar = (2*Aga[l][0](i+1,0,k).AcetateAgar+2*Aga[l][0](i,0,k-1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*DiffAgarAce*dt/h2;
                    }
                    else
                    {
                        AgaRes[l][0](i,0,k).CarbonAgar = (2*Aga[l][0](i+1,0,k).CarbonAgar+Aga[l][0](i,0,k-1).CarbonAgar+Aga[l][0](i,0,k+1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*DiffAgarCarbon*dt/h2;
                        AgaRes[l][0](i,0,k).OxygenAgar = 0;//maxO2-Aga[l][0](i,0,k).OxygenAgar;
                        AgaRes[l][0](i,0,k).AcetateAgar = (2*Aga[l][0](i+1,0,k).AcetateAgar+Aga[l][0](i,0,k-1).AcetateAgar+Aga[l][0](i,0,k+1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*DiffAgarAce*dt/h2;
                    }
                }
                else
                {
                    if (k==0)
                    {
                        AgaRes[l][0](i,0,k).CarbonAgar = (2*Aga[l][0](i-1,0,k).CarbonAgar+2*Aga[l][0](i,0,k+1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*DiffAgarCarbon*dt/h2;
                        AgaRes[l][0](i,0,k).OxygenAgar = 0;//maxO2-Aga[l][0](i,0,k).OxygenAgar;
                        AgaRes[l][0](i,0,k).AcetateAgar = (2*Aga[l][0](i-1,0,k).AcetateAgar+2*Aga[l][0](i,0,k+1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*DiffAgarAce*dt/h2;
                    }
                    else if (k==nz)
                    {
                        AgaRes[l][0](i,0,k).CarbonAgar = (2*Aga[l][0](i-1,0,k).CarbonAgar+2*Aga[l][0](i,0,k-1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*DiffAgarCarbon*dt/h2;
                        AgaRes[l][0](i,0,k).OxygenAgar = 0;//maxO2-Aga[l][0](i,0,k).OxygenAgar;
                        AgaRes[l][0](i,0,k).AcetateAgar = (2*Aga[l][0](i-1,0,k).AcetateAgar+2*Aga[l][0](i,0,k-1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*DiffAgarAce*dt/h2;
                    }
                    else
                    {
                        AgaRes[l][0](i,0,k).CarbonAgar = (2*Aga[l][0](i-1,0,k).CarbonAgar+Aga[l][0](i,0,k-1).CarbonAgar+Aga[l][0](i,0,k+1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*DiffAgarCarbon*dt/h2;
                        AgaRes[l][0](i,0,k).OxygenAgar = 0;//maxO2-Aga[l][0](i,0,k).OxygenAgar;
                        AgaRes[l][0](i,0,k).AcetateAgar = (2*Aga[l][0](i-1,0,k).AcetateAgar+Aga[l][0](i,0,k-1).AcetateAgar+Aga[l][0](i,0,k+1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*DiffAgarAce*dt/h2;
                    }
                }
            }
        }
        else if (k==nz)
        {
            if (l==agarLevels-1)
            {
                AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l][0](i+1,0,k).CarbonAgar+Aga[l][0](i-1,0,k).CarbonAgar+2*Aga[l][0](i,0,k-1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*dt*DiffAgarCarbon/h2;
                AgaRes[l][0](i,0,k).OxygenAgar = (Aga[l][0](i+1,0,k).OxygenAgar+Aga[l][0](i-1,0,k).OxygenAgar+2*Aga[l][0](i,0,k-1).OxygenAgar-4*Aga[l][0](i,0,k).OxygenAgar)*DiffAgarO2/h2;
                AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l][0](i+1,0,k).AcetateAgar+Aga[l][0](i-1,0,k).AcetateAgar+2*Aga[l][0](i,0,k-1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*dt*DiffAgarAce/h2;
            }
        }
        else if (k==0)
        {
            if (l==0)
            {
                if (insideColonyDen[0](i,0,0)>0)
                {
                    AgaRes[l][0](i,0,k).CarbonAgar = 0;
                    AgaRes[l][0](i,0,k).OxygenAgar = 0;
                    AgaRes[l][0](i,0,k).AcetateAgar = 0;
                }
                else
                {
                    AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l][0](i+1,0,k).CarbonAgar+Aga[l][0](i-1,0,k).CarbonAgar+2*Aga[l][0](i,0,k+1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*dt*DiffAgarCarbon/h2;
                    AgaRes[l][0](i,0,k).OxygenAgar = maxO2_Agar-Aga[l][0](i,0,k).OxygenAgar;
                    AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l][0](i+1,0,k).AcetateAgar+Aga[l][0](i-1,0,k).AcetateAgar+2*Aga[l][0](i,0,k+1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*dt*DiffAgarAce/h2;
                }
            }
            else if ((ip>0 && ip<nxp)!=1)
            {
                if (ip==0)
                {
                    AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l-1][0](ip+2,0,kp).CarbonAgar+Aga[l][0](i-1,0,k).CarbonAgar+2*Aga[l][0](i,0,k+1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*dt*DiffAgarCarbon/h2;
                    AgaRes[l][0](i,0,k).OxygenAgar = maxO2_Agar-Aga[l][0](i,0,k).OxygenAgar;
                    AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l-1][0](ip+2,0,kp).AcetateAgar+Aga[l][0](i-1,0,k).AcetateAgar+2*Aga[l][0](i,0,k+1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*dt*DiffAgarAce/h2;
                }
                else if (ip==nxp)
                {
                    AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l-1][0](ip-2,0,kp).CarbonAgar+Aga[l][0](i+1,0,k).CarbonAgar+2*Aga[l][0](i,0,k+1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*dt*DiffAgarCarbon/h2;
                    AgaRes[l][0](i,0,k).OxygenAgar = maxO2_Agar-Aga[l][0](i,0,k).OxygenAgar;
                    AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l-1][0](ip-2,0,kp).AcetateAgar+Aga[l][0](i+1,0,k).AcetateAgar+2*Aga[l][0](i,0,k+1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*dt*DiffAgarAce/h2;
                }
                else
                {
                    AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l][0](i+1,0,k).CarbonAgar+Aga[l][0](i-1,0,k).CarbonAgar+2*Aga[l][0](i,0,k+1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*dt*DiffAgarCarbon/h2;
                    AgaRes[l][0](i,0,k).OxygenAgar = maxO2_Agar-Aga[l][0](i,0,k).OxygenAgar;
                    AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l][0](i+1,0,k).AcetateAgar+Aga[l][0](i-1,0,k).AcetateAgar+2*Aga[l][0](i,0,k+1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*dt*DiffAgarAce/h2;
                }
            }
        }
        else
        {
            if (l>0)
            {
                if (kp==nzp && ip>0 && ip<nxp)
                {
                    AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l][0](i+1,0,k).CarbonAgar+Aga[l][0](i-1,0,k).CarbonAgar+Aga[l-1][0](ip,0,kp-2).CarbonAgar+Aga[l][0](i,0,k+1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*dt*DiffAgarCarbon/h2;
                    AgaRes[l][0](i,0,k).OxygenAgar = (Aga[l][0](i+1,0,k).OxygenAgar+Aga[l][0](i-1,0,k).OxygenAgar+Aga[l-1][0](ip,0,kp-2).OxygenAgar+Aga[l][0](i,0,k+1).OxygenAgar-4*Aga[l][0](i,0,k).OxygenAgar)*dt*DiffAgarO2/h2;
                    AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l][0](i+1,0,k).AcetateAgar+Aga[l][0](i-1,0,k).AcetateAgar+Aga[l-1][0](ip,0,kp-2).AcetateAgar+Aga[l][0](i,0,k+1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*dt*DiffAgarAce/h2;
                
                }
                else if (ip==0 && kp<nzp)
                {
                    AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l-1][0](ip+2,0,kp).CarbonAgar+Aga[l][0](i-1,0,k).CarbonAgar+Aga[l][0](i,0,k-1).CarbonAgar+Aga[l][0](i,0,k+1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*dt*DiffAgarCarbon/h2;
                    AgaRes[l][0](i,0,k).OxygenAgar = (Aga[l-1][0](ip+2,0,kp).OxygenAgar+Aga[l][0](i-1,0,k).OxygenAgar+Aga[l][0](i,0,k-1).OxygenAgar+Aga[l][0](i,0,k+1).OxygenAgar-4*Aga[l][0](i,0,k).OxygenAgar)*dt*DiffAgarO2/h2;
                    AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l-1][0](ip+2,0,kp).AcetateAgar+Aga[l][0](i-1,0,k).AcetateAgar+Aga[l][0](i,0,k-1).AcetateAgar+Aga[l][0](i,0,k+1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*dt*DiffAgarAce/h2;
                }
                else if (ip==nxp && kp<nzp)
                {
                    AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l][0](i+1,0,k).CarbonAgar+Aga[l-1][0](ip-2,0,kp).CarbonAgar+Aga[l][0](i,0,k-1).CarbonAgar+Aga[l][0](i,0,k+1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*dt*DiffAgarCarbon/h2;
                    AgaRes[l][0](i,0,k).OxygenAgar = (Aga[l][0](i+1,0,k).OxygenAgar+Aga[l-1][0](ip-2,0,kp).OxygenAgar+Aga[l][0](i,0,k-1).OxygenAgar+Aga[l][0](i,0,k+1).OxygenAgar-4*Aga[l][0](i,0,k).OxygenAgar)*dt*DiffAgarO2/h2;
                    AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l][0](i+1,0,k).AcetateAgar+Aga[l-1][0](ip-2,0,kp).AcetateAgar+Aga[l][0](i,0,k-1).AcetateAgar+Aga[l][0](i,0,k+1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*dt*DiffAgarAce/h2;
                }
                else if (ip<=0 || ip>=nxp || kp>=nzp)
                {
                    AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l][0](i+1,0,k).CarbonAgar+Aga[l][0](i-1,0,k).CarbonAgar+Aga[l][0](i,0,k-1).CarbonAgar+Aga[l][0](i,0,k+1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*dt*DiffAgarCarbon/h2;
                    AgaRes[l][0](i,0,k).OxygenAgar = (Aga[l][0](i+1,0,k).OxygenAgar+Aga[l][0](i-1,0,k).OxygenAgar+Aga[l][0](i,0,k-1).OxygenAgar+Aga[l][0](i,0,k+1).OxygenAgar-4*Aga[l][0](i,0,k).OxygenAgar)*dt*DiffAgarO2/h2;
                    AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l][0](i+1,0,k).AcetateAgar+Aga[l][0](i-1,0,k).AcetateAgar+Aga[l][0](i,0,k-1).AcetateAgar+Aga[l][0](i,0,k+1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*dt*DiffAgarAce/h2;
                }
            }
            else
            {
                AgaRes[l][0](i,0,k).CarbonAgar = (Aga[l][0](i+1,0,k).CarbonAgar+Aga[l][0](i-1,0,k).CarbonAgar+Aga[l][0](i,0,k-1).CarbonAgar+Aga[l][0](i,0,k+1).CarbonAgar-4*Aga[l][0](i,0,k).CarbonAgar)*dt*DiffAgarCarbon/h2;
                AgaRes[l][0](i,0,k).OxygenAgar = (Aga[l][0](i+1,0,k).OxygenAgar+Aga[l][0](i-1,0,k).OxygenAgar+Aga[l][0](i,0,k-1).OxygenAgar+Aga[l][0](i,0,k+1).OxygenAgar-4*Aga[l][0](i,0,k).OxygenAgar)*dt*DiffAgarO2/h2;
                AgaRes[l][0](i,0,k).AcetateAgar = (Aga[l][0](i+1,0,k).AcetateAgar+Aga[l][0](i-1,0,k).AcetateAgar+Aga[l][0](i,0,k-1).AcetateAgar+Aga[l][0](i,0,k+1).AcetateAgar-4*Aga[l][0](i,0,k).AcetateAgar)*dt*DiffAgarAce/h2;
            }
        }
        
    }
}

void Multigrid::NestSolAga(int l)
{
    int d = ipow(2,l);
    int d0 = ipow(2,agarLevels);
    int nx = min(BoxX, BoxX*d0/d/2);
    int nz = min(BoxZAgar, BoxZAgar*d0/d/2);
    int nxp = min(BoxX, BoxX*d0/d);
    int nzp = min(BoxZAgar, BoxZAgar*d0/d);
    double h2 = BoxLength*BoxLength*d*d;
    int i,k;
#pragma omp parallel for default(shared) private(i,k)
    for ( i=0; i<nx+1; i++)
    {
        for ( k=0;k<nz+1;k++)
        {
            NestSolAgaElementwise(l,i,0,k,nx,BoxY,nz,nxp,BoxY,nzp,h2);
        }
    }
}

void Multigrid::NestSolAgaElementwise(int l, int i, int j, int k, int nx, int ny, int nz, int nxp, int nyp, int nzp, double h2)
{
    int ip = 2*i - nx+nxp/2;
    int kp = 2*k;
    if (l<agarLevels)
    {
        if (ip>0 && ip<nxp && kp<nzp)
        {
            Aga[l][0](i,0,k).CarbonAgar = Aga[l-1][0](ip,0,kp).CarbonAgar;
            Aga[l][0](i,0,k).OxygenAgar = Aga[l-1][0](ip,0,kp).OxygenAgar;
            Aga[l][0](i,0,k).AcetateAgar = Aga[l-1][0](ip,0,kp).AcetateAgar;
        }
    }
}

double Multigrid::evalUpdateErr()
{
    double er2Carbon = 0;
    double erinfCarbon = 0;
    double er2Oxygen = 0;
    double erinfOxygen = 0;
    double er2Acetate = 0;
    double erinfAcetate = 0;
    int i,k;
    for (int l=0; l<agarLevels; l++)
    {
        int d = ipow(2,l);
        int d0 = ipow(2,agarLevels);
        int nx = min(BoxX, BoxX*d0/d/2);
        int nz = min(BoxZAgar, BoxZAgar*d0/d/2);
        double h2 = BoxLength*BoxLength*d*d;
#pragma omp parallel for default(shared) private(i,k) reduction(+:er2Carbon,er2Oxygen,er2Acetate) reduction(max : erinfCarbon, erinfOxygen, erinfAcetate)
        for ( i=0; i<=nx; i++)
        {
            for ( k=0; k<nz; k++)
            {
                er2Carbon += max(abs(AgaErr[0][0](i,0,k).CarbonAgar),0.0)*h2;
                erinfCarbon = max(erinfCarbon, max(abs(AgaErr[0][0](i,0,k).CarbonAgar),0.0));
                er2Oxygen += max(abs(AgaErr[0][0](i,0,k).OxygenAgar),0.0)*h2;
                erinfOxygen = max(erinfOxygen, max(abs(AgaErr[0][0](i,0,k).OxygenAgar),0.0));
                er2Acetate += max(abs(AgaErr[0][0](i,0,k).AcetateAgar),0.0)*h2;
                erinfAcetate = max(erinfAcetate, max(abs(AgaErr[0][0](i,0,k).AcetateAgar),0.0));
            }
        }
    }
//    cout<<" "<<erinfCarbon<<" "<<" "<<erinfOxygen<<" "<<erinfAcetate<<" ";
    return max(max(erinfCarbon,erinfAcetate),erinfOxygen);
}

double Multigrid::evalUpdateRes()
{
    double res2Carbon = 0;
    double resinfCarbon = 0;
    double res2Oxygen = 0;
    double resinfOxygen = 0;
    double res2Acetate = 0;
    double resinfAcetate = 0;
    int i,k;
    for (int l=0; l<1; l++)
    {
        int d = ipow(2,l);
        int d0 = ipow(2,agarLevels);
        int nx = min(BoxX, BoxX*d0/d/2);
        int nz = min(BoxZAgar, BoxZAgar*d0/d/2);
        double h2 = BoxLength*BoxLength*d*d;
#pragma omp parallel for default(shared) private(i,k) reduction(+:res2Carbon, res2Oxygen, res2Acetate) reduction(max : resinfCarbon, resinfOxygen, resinfAcetate)
        for ( i=0; i<=nx; i++)
        {
            for ( k=0; k<nz; k++)
            {
                res2Carbon += max(abs(AgaRes[l][0](i,0,k).CarbonAgar),0.0)*h2;
                resinfCarbon = max(resinfCarbon, max(abs(AgaRes[l][0](i,0,k).CarbonAgar),0.0));
                res2Oxygen += max(abs(AgaRes[l][0](i,0,k).OxygenAgar),0.0)*h2;
                resinfOxygen = max(resinfOxygen, max(abs(AgaRes[l][0](i,0,k).OxygenAgar),0.0));
                res2Acetate += max(abs(AgaRes[l][0](i,0,k).AcetateAgar),0.0)*h2;
                resinfAcetate = max(resinfAcetate, max(abs(AgaRes[l][0](i,0,k).AcetateAgar),0.0));
            }
        }
    }
//    cout<<"     "<<resinfCarbon<<" "<<resinfOxygen<<" "<<resinfAcetate<<"     "<<endl;
    return max(max(resinfCarbon,resinfOxygen),resinfAcetate);
}

void Multigrid::setup()
{
    resetZero();
    for (int l=0; l<=agarLevels-1; l++)
        updateSolAga(l);
    for (int l=0; l<=agarLevels-1; l++)
        updateResAga(l);
}

void Multigrid::reset()
{
    for (int l=0; l<=agarLevels-1; l++)
        updateSolAga(l);
    for (int l=0; l<agarLevels-1; l++)
        rewriteSolBdryAga(l);
    for (int l=0; l<=agarLevels-1; l++)
        updateResAga(l);
    evalUpdateRes();
    
    resetZero();
    
    for (int l=1; l<=agarLevels; l++)
        NestSolAga(l);
}

void Multigrid::Vcycle(int mu1, int mu2)
{    
    for (int l=0; l<agarMGLevels-1; l++)
    {
        for (int j=0;j<mu1; j++)
            relaxationAgarAtLevel(l);
        residueRestrictionAga(l+1);
    }
    
    for (int j=0;j<mu2; j++)
        relaxationAgarAtLevel(agarMGLevels-1);
    
    for (int l=agarMGLevels-2; l>=0; l--)
    {
        interpolationAga(l);
        for (int j=0;j<mu1; j++)
            relaxationAgarAtLevel(l);
    }
    
}

void Multigrid::solve(Array3D<LocalAga>** CurrentAgar, Array2D<LocalAga>** CurrentWall, Array3D<LocalAga>** PreviousAgar, Array2D<LocalAga>** PreviousWall, int mu1, int mu2, double Dt)
{
    dt = Dt;
    copyToMG(CurrentAgar, CurrentWall);
    setup();
    for (int i = 0; i<3; i++)
    {
        Vcycle(mu1, mu2);
	/*FILE* f1, *f2;
	char f1_name[500];
	char f2_name[500];
	sprintf(f1_name,"output/sol%d.txt",i);
	sprintf(f2_name,"output/err%d.txt",i);
	f1 = fopen(f1_name, "w");
	f2 = fopen(f2_name, "w");
	for (int level=0;level<3;level++)
            {
                Aga[level]->Append(f1,1);
                AgaErr[level]->Append(f2,1);
            }
	fclose(f1);
	fclose(f2);*/
    }
    reset();
    copyFromMG(CurrentAgar, CurrentWall, PreviousAgar, PreviousWall);
}

void Multigrid::copyToMG(Array3D<LocalAga>** CurrentAgar, Array2D<LocalAga>** CurrentWall)
{
    for (int i=1; i<=BoxX-1; i++)
    {
        if (insideColonyDen[0](i,0,0)>0)
        {
            Aga[0][0](i,0,0).CarbonAgar = CurrentWall[0][0](i,0).CarbonAgar;
            Aga[0][0](i,0,0).OxygenAgar = CurrentWall[0][0](i,0).OxygenAgar;
            Aga[0][0](i,0,0).AcetateAgar = CurrentWall[0][0](i,0).AcetateAgar;
        }
    }
    for (int k=0; k<=BoxZAgar; k++)
    {
        Aga[agarLevels-1][0](BoxX,0,k).CarbonAgar = Aga[agarLevels-1][0](0,0,k).CarbonAgar;
        Aga[agarLevels-1][0](BoxX,0,k).OxygenAgar = Aga[agarLevels-1][0](0,0,k).OxygenAgar;
        Aga[agarLevels-1][0](BoxX,0,k).AcetateAgar =Aga[agarLevels-1][0](0,0,k).AcetateAgar;
    }
    /*for (int l=0; l<agarLevels; l++)
    {
        for (int i=1; i<=BoxX-1; i++)
        {
            Aga[l][0](i,0,0).CarbonAgar = CurrentWall[l][0](i,0).CarbonAgar;
            Aga[l][0](i,0,0).OxygenAgar = CurrentWall[l][0](i,0).OxygenAgar;
            Aga[l][0](i,0,0).AcetateAgar = CurrentWall[l][0](i,0).AcetateAgar;
            for (int k=0; k<=BoxZAgar-1; k++)
            {
                Aga[l][0](i,0,k+1).CarbonAgar = CurrentAgar[l][0](i,0,k).CarbonAgar;
                Aga[l][0](i,0,k+1).OxygenAgar = CurrentAgar[l][0](i,0,k).OxygenAgar;
                Aga[l][0](i,0,k+1).AcetateAgar = CurrentAgar[l][0](i,0,k).AcetateAgar;
            }
        }
        if (l<agarLevels-1)
        {
            Aga[l][0](BoxX,0,0).CarbonAgar = CurrentWall[l+1][0](BoxX*3/4,0).CarbonAgar;
            Aga[l][0](BoxX,0,0).OxygenAgar = CurrentWall[l+1][0](BoxX*3/4,0).OxygenAgar;
            Aga[l][0](BoxX,0,0).AcetateAgar = CurrentWall[l+1][0](BoxX*3/4,0).AcetateAgar;
            for (int k=1; k<=BoxZAgar; k++)
            {
                if (k%2==0)
                {
                    Aga[l][0](BoxX,0,k).CarbonAgar = CurrentAgar[l+1][0](BoxX*3/4,0,k/2).CarbonAgar;
                    Aga[l][0](BoxX,0,k).OxygenAgar = CurrentAgar[l+1][0](BoxX*3/4,0,k/2).OxygenAgar;
                    Aga[l][0](BoxX,0,k).AcetateAgar = CurrentAgar[l+1][0](BoxX*3/4,0,k/2).AcetateAgar;
                }
                else
                {
                    Aga[l][0](BoxX,0,k).CarbonAgar = 0.5*(CurrentAgar[l+1][0](BoxX*3/4,0,k/2).CarbonAgar+CurrentAgar[l+1][0](BoxX*3/4,0,k/2+1).CarbonAgar);
                    Aga[l][0](BoxX,0,k).OxygenAgar = 0.5*(CurrentAgar[l+1][0](BoxX*3/4,0,k/2).OxygenAgar+CurrentAgar[l+1][0](BoxX*3/4,0,k/2+1).OxygenAgar);
                    Aga[l][0](BoxX,0,k).AcetateAgar = 0.5*(CurrentAgar[l+1][0](BoxX*3/4,0,k/2).AcetateAgar+CurrentAgar[l+1][0](BoxX*3/4,0,k/2+1).AcetateAgar);
                }
            }
        }
        else
        {
            Aga[l][0](BoxX,0,0).CarbonAgar = 2*CurrentWall[l][0](BoxX-1,0).CarbonAgar-CurrentWall[l][0](BoxX-2,0).CarbonAgar;
            Aga[l][0](BoxX,0,0).OxygenAgar = 2*CurrentWall[l][0](BoxX-1,0).OxygenAgar-CurrentWall[l][0](BoxX-2,0).OxygenAgar;
            Aga[l][0](BoxX,0,0).AcetateAgar = 2*CurrentWall[l][0](BoxX-1,0).AcetateAgar-CurrentWall[l][0](BoxX-2,0).AcetateAgar;
            for (int k=0; k<=BoxZAgar-1; k++)
            {
                Aga[l][0](BoxX,0,k+1).CarbonAgar = 2*CurrentAgar[l][0](BoxX-1,0,k).CarbonAgar-CurrentAgar[l][0](BoxX-2,0,k).CarbonAgar;
                Aga[l][0](BoxX,0,k+1).OxygenAgar = 2*CurrentAgar[l][0](BoxX-1,0,k).OxygenAgar-CurrentAgar[l][0](BoxX-2,0,k).OxygenAgar;
                Aga[l][0](BoxX,0,k+1).AcetateAgar = 2*CurrentAgar[l][0](BoxX-1,0,k).AcetateAgar-CurrentAgar[l][0](BoxX-2,0,k).AcetateAgar;
            }
        }
    }*/
}

void Multigrid::copyFromMG(Array3D<LocalAga>** CurrentAgar, Array2D<LocalAga>** CurrentWall, Array3D<LocalAga>** PreviousAgar, Array2D<LocalAga>** PreviousWall)
{
    for (int l=0; l<agarLevels; l++)
    {
        for (int i=0; i<=BoxX-1; i++)
        {
            CurrentWall[l][0](i,0).CarbonAgar = Aga[l][0](i,0,0).CarbonAgar;
            CurrentWall[l][0](i,0).OxygenAgar = Aga[l][0](i,0,0).OxygenAgar;
            CurrentWall[l][0](i,0).AcetateAgar = Aga[l][0](i,0,0).AcetateAgar;
            PreviousWall[l][0](i,0).CarbonAgar = Aga[l][0](i,0,0).CarbonAgar;
            PreviousWall[l][0](i,0).OxygenAgar = Aga[l][0](i,0,0).OxygenAgar;
            PreviousWall[l][0](i,0).AcetateAgar = Aga[l][0](i,0,0).AcetateAgar;
            for (int k=0; k<=BoxZAgar-1; k++)
            {
                CurrentAgar[l][0](i,0,k).CarbonAgar = Aga[l][0](i,0,k+1).CarbonAgar;
                CurrentAgar[l][0](i,0,k).OxygenAgar = Aga[l][0](i,0,k+1).OxygenAgar;
                CurrentAgar[l][0](i,0,k).AcetateAgar = Aga[l][0](i,0,k+1).AcetateAgar;
                PreviousAgar[l][0](i,0,k).CarbonAgar = Aga[l][0](i,0,k+1).CarbonAgar;
                PreviousAgar[l][0](i,0,k).OxygenAgar = Aga[l][0](i,0,k+1).OxygenAgar;
                PreviousAgar[l][0](i,0,k).AcetateAgar = Aga[l][0](i,0,k+1).AcetateAgar;
            }
        }
    }
}


int UpdateEnvArrayMG_OxygenAcetate(Array3D<LocalEnv>* CurrentColony, Array3D<LocalEnv>* PreviousColony, Array3D<LocalAga>** CurrentAgar, Array3D<LocalAga>** PreviousAgar, Array2D<LocalAga>** CurrentWall, Array2D<LocalAga>** PreviousWall, Array3D<double>& DensityShiftP, Array3D<double>& Density1ShiftP, Array3D<double>& Density2ShiftP, Array2D<double>& WallDensityShiftP, Array2D<double>& WallDensity1ShiftP, Array2D<double>& WallDensity2ShiftP, int minX, int maxX, int minY, int maxY, int maxH, Array2D<double>& Height, Array3D<double>& insideColonyDen, Multigrid& mg)
{
    // loop through nutrient array until steady state
    int count = 0;
    IntCoord positionColony; // position in the colony
    IntCoord positionAgar; // position in the agar
    IntCoord2D positionWall; // position on the wall
    int conv = 1;
    int nbk = 1;
    int zmin, zmax, dir;
    int zi, yi, xi;
    double dt = 100;//min(DiffAgarAce, DiffAgarCarbon)*BoxLength*BoxLength*ipow(4,mg.getLevels());
    int n = 3600*UpdateTime*740/dt;
    double previous_carbon, previous_waste, previous_oxygen, previous_acetate;
//    while (count<maxIter && nbk)
    for (int i=0; i<n; i++)
    {
        
        mg.solve(CurrentAgar, CurrentWall, PreviousAgar, PreviousWall, 3, 100, dt);
        /*FILE* f1, *f2;
        char f1_name[500];
        char f2_name[500];
        sprintf(f1_name,"output2/sol%d.txt",i);
        sprintf(f2_name,"output2/err%d.txt",i);
        f1 = fopen(f1_name, "w");
        f2 = fopen(f2_name, "w");
        for (int level=0;level<3;level++)
        {
            mg.Aga[level]->Append(f1,1);
            mg.AgaErr[level]->Append(f2,1);
        }
        fclose(f1);
        fclose(f2);*/
        for (int j=0; j<200; j++)
        {
            if (count%2==0) {
                zmin = 0;
                zmax = maxH+1;
                dir = 1;
            } else {
                zmin = -maxH-1;
                zmax = 0;
                dir = -1;
            }
            
#pragma omp parallel for default(shared) private(zi,yi,xi,positionColony,positionAgar,positionWall,previous_carbon,previous_waste,previous_acetate,previous_oxygen) reduction(+:conv)
        for (zi=zmin;zi<=zmax;zi++)
        {
            // Colony
            if (zi*dir>0)
            {
                for (yi=minY; yi<=maxY; yi++)
                {
                    for (xi=minX; xi<=maxX; xi++)
                    {
                        
                        positionColony.x = xi;
                        positionColony.y = yi;
                        positionColony.z = zi*dir-1;
                        
                        previous_carbon = PreviousColony->Get(positionColony).Carbon;
                        previous_waste  = PreviousColony->Get(positionColony).Waste;
                        previous_acetate = PreviousColony->Get(positionColony).Acetate;
                        previous_oxygen = PreviousColony->Get(positionColony).Oxygen;
                        
                        if (DensityShiftP.Get(positionColony) > DBL_EPSILON)
                        {
                            UpdateEnvironment_OxygenAcetate(CurrentColony, PreviousColony, PreviousWall[0], DensityShiftP, Density1ShiftP, Density2ShiftP,WallDensityShiftP,positionColony,insideColonyDen);
                        }
                        
                        //if (xi==50 && yi ==50 && zi==0)
                        //printf("ColonyCarbon: %6f \n",CurrentColony->Get(positionColony).Carbon);
                        
                        // check convergence criteria for this grid point
                        conv += int( fabs(CurrentColony->Get(positionColony).Carbon - previous_carbon) > ConvCrit*max(0.1,fabs(CurrentColony->Get(positionColony).Carbon)));
                        //conv += int( fabs(CurrentColony->Get(positionColony).Waste - previous_waste) > ConvCrit*fabs(CurrentColony->Get(positionColony).Waste) );
                        conv += int( fabs(CurrentColony->Get(positionColony).Acetate - previous_acetate) > ConvCrit*max(0.1,fabs(CurrentColony->Get(positionColony).Acetate)));
                        conv += int( fabs(CurrentColony->Get(positionColony).Oxygen - previous_oxygen) > ConvCrit*max(0.1,fabs(CurrentColony->Get(positionColony).Oxygen)));
                        
                    }
                }
            }
            // Colony and Agar interaction surface
            else if (zi==0)
            {
                for (int level=0; level<maxLevelsMG; level++)
                {
                    for (yi=0; yi<BoxY; yi++)
                    {
                        for (xi=0; xi<BoxX; xi++)
                        {
                            positionWall.x=xi;
                            positionWall.y=yi;
                            
                            previous_carbon = PreviousWall[level]->Get(positionWall).CarbonAgar;
                            previous_waste  = PreviousWall[level]->Get(positionWall).WasteAgar;
                            previous_oxygen = PreviousWall[level]->Get(positionWall).OxygenAgar;
                            previous_acetate = PreviousWall[level]->Get(positionWall).AcetateAgar;
                            
                            UpdateWallMultigrid_OxygenAcetate(CurrentColony, PreviousAgar, CurrentWall, PreviousWall, WallDensityShiftP, WallDensity1ShiftP, WallDensity2ShiftP, positionWall, level, Height, insideColonyDen);
                            
                            // check convergence criteria for this grid point
                            
                            conv += int( fabs(CurrentWall[level]->Get(positionWall).CarbonAgar - previous_carbon) > ConvCrit*max(0.1,fabs(CurrentWall[level]->Get(positionWall).CarbonAgar)));
                            // conv += int( fabs(CurrentWall[level]->Get(positionWall).WasteAgar - previous_waste) > ConvCrit*fabs(CurrentWall[level]->Get(positionWall).WasteAgar) );
                            conv += int( fabs(CurrentWall[level]->Get(positionWall).OxygenAgar - previous_oxygen) > ConvCrit*max(0.1,fabs(CurrentWall[level]->Get(positionWall).OxygenAgar)));
                            conv += int( fabs(CurrentWall[level]->Get(positionWall).AcetateAgar - previous_acetate) > ConvCrit*max(0.1,fabs(CurrentWall[level]->Get(positionWall).AcetateAgar)));
                        }
                    }
                }
            }
        }
        
        
        
        if (NutrientGSI==0)
        {
            Array3D<LocalEnv>* swapColony = PreviousColony;
            PreviousColony = CurrentColony;
            CurrentColony = swapColony;
            
            Array3D<LocalAga>** swapAgar = PreviousAgar;
            PreviousAgar = CurrentAgar;
            CurrentAgar = swapAgar;
            
            Array2D<LocalAga>** swapWall = PreviousWall;
            PreviousWall = CurrentWall;
            CurrentWall = swapWall;
        }
        
        //printf("count: %d   conv: %d \n",count, conv);
//        if ((count>minIter)&&(conv==0))
//            nbk=0;
        
        count++;
            
        }

    }
    return count;
    
}

