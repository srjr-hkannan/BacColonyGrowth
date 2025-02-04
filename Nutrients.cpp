#include <float.h>
#include <math.h>
#include "Nutrients.h"
#include "Array.h"
#include "Constants.h"
#include <omp.h>

/*************************************************************************************
 * This solves for the concentration of nutrients (C = carbon, O = oxygen, A = acetate) within the colony.
 * Diffusion in the agar is NOT computed, but the agar-colony interface is treated like a boundary condition.
 * UpdateEnvArray loops through the entire environment array and updates. UpdateEnvironment does the reaction-diffusion calculation for
 * each individual cell in the environment array.
 */

int levelM1x(int x)
{
    return x*2-BoxX/2;
}

int levelP1x(int x)
{
    return x/2+BoxX/4;
}

int levelM1y(int y)
{
    return y*2-BoxY/2;
}

int levelP1y(int y)
{
    return y/2+BoxY/4;
}

int levelM1z(int z)
{
    return z*2+1;
}

int levelP1z(int z)
{
    return (z-1)/2;
}

int UpdateEnvArray_OxygenAcetate(Array3D<LocalEnv>* CurrentColony, Array3D<LocalEnv>* PreviousColony, Array3D<LocalAga>** CurrentAgar, Array3D<LocalAga>** PreviousAgar, Array2D<LocalAga>** CurrentWall, Array2D<LocalAga>** PreviousWall, Array3D<double>& DensityShiftP, Array3D<double>& Density1ShiftP, Array3D<double>& Density2ShiftP, Array2D<double>& WallDensityShiftP, Array2D<double>& WallDensity1ShiftP, Array2D<double>& WallDensity2ShiftP, int minX, int maxX, int minY, int maxY, int maxH, Array2D<double>& Height, Array3D<double>& insideColonyDen)
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

    double previous_carbon, previous_waste, previous_oxygen, previous_acetate;
    
    while (count<maxIter && nbk)
    {
        conv = 0;
        if (count%2==0) {
            zmin = -BoxZAgar;
            zmax = maxH+1;
            dir = 1;
        } else {
            zmin = -maxH-1;
            zmax = BoxZAgar;
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
            // Agar Region
            else
            {
                for (int level=0; level<maxLevelsMG; level++)
                {
                    for (yi=0; yi<BoxY; yi++)
                    {
                        for (xi=0; xi<BoxX; xi++)
                        {
                            positionAgar.x = xi;
                            positionAgar.y = yi;
                            positionAgar.z = -(zi*dir+1);

                            previous_carbon = PreviousAgar[level]->Get(positionAgar).CarbonAgar;
                            previous_waste  = PreviousAgar[level]->Get(positionAgar).WasteAgar;
                            previous_oxygen = PreviousAgar[level]->Get(positionAgar).OxygenAgar;
                            previous_acetate  = PreviousAgar[level]->Get(positionAgar).AcetateAgar;
                            
                            UpdateAgarMultigrid_OxygenAcetate(CurrentAgar, PreviousAgar, CurrentWall, positionAgar, level);

                            // check convergence criteria for this grid point
                            
                            conv += int( fabs(CurrentAgar[level]->Get(positionAgar).CarbonAgar - previous_carbon) > ConvCrit*max(0.1,fabs(CurrentAgar[level]->Get(positionAgar).CarbonAgar)));
                            //conv += int( fabs(CurrentAgar[level]->Get(positionAgar).WasteAgar - previous_waste) > ConvCrit*fabs(CurrentAgar[level]->Get(positionAgar).WasteAgar) );
                            conv += int( fabs(CurrentAgar[level]->Get(positionAgar).OxygenAgar - previous_oxygen) > ConvCrit*max(0.1,fabs(CurrentAgar[level]->Get(positionAgar).OxygenAgar)));
                            conv += int( fabs(CurrentAgar[level]->Get(positionAgar).AcetateAgar - previous_acetate) > ConvCrit*max(0.1,fabs(CurrentAgar[level]->Get(positionAgar).AcetateAgar)));
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
        if ((count>minIter)&&(conv==0))
            nbk=0;

        count++;


    }

    return count;

}

// Colony area
void UpdateEnvironment_OxygenAcetate(Array3D<LocalEnv>* Env, Array3D<LocalEnv>* prevEnv, Array2D<LocalAga>* prevWal, Array3D<double>& DensityShiftP, Array3D<double>& Density1ShiftP,Array3D<double>& Density2ShiftP,Array2D<double>& WallDensityShiftP,IntCoord position, Array3D<double>& insideColonyDen)
{
    IntCoord2D positionwall;
    positionwall.x=position.x;
    positionwall.y=position.y;
    // Inside Colony
    if (insideColonyDen.Get(position)>0.5)//((position.z==0 && WallDensityShiftP.Get(positionwall) > 0.1) || DensityShiftP.Get(position)>0.05)	// if density is zero, do not update
    {
       
        double C = prevEnv->Get(position).Carbon;
        double W = prevEnv->Get(position).Waste;
		double A = prevEnv->Get(position).Acetate;
		double O = prevEnv->Get(position).Oxygen;
        
        double HA = acetateBufferSolve(A, prevEnv, position);

        // calculate growth rate and consumption rates based on current concentrations
        double fC = C/(C+KC_aer);	
        double fW = W/(W+KW);	
	    double fA = A/(A+KA);	
	    double fO = O/(O+KO);	
	    //double fO = 26/27;	
	    double fI = 0;//sqrt(HA/KI); // acetic acid toxic effect

        double Cgr=0; // redistribute C, W and A to vertices, in order to compute growth rate
        double Wgr=0;
		double Agr=0;
        double Ogr=0;

	double mgrCarAers, mgrCarAnas, mgrAceAers, g1s, g2s, g3s;
        // calculate constants for maintenance
        if (maxGrowthRateCarAer != 0.0)
	{
		mgrCarAers = qCarAernt/qCarAer;
		g1s = KC_aer*(mgrCarAers/maxGrowthRateCarAer);	
	}
	else
	{
		mgrCarAers = 0;
		g1s = 1e6;
	}

	if (maxGrowthRateCarAna != 0.0)
        {
                mgrCarAnas = qCarAnant/qCarAna;
                g2s = KC_ana*(mgrCarAnas/maxGrowthRateCarAna);
        }
        else
        {
                mgrCarAnas = 0;
                g2s = 1e6;
        }

	if (maxGrowthRateAceAer != 0.0)
        {
                mgrAceAers = qAceAernt/qAceAer;
                g3s = KA*(mgrAceAers/maxGrowthRateAceAer);
        }
        else
        {
                mgrAceAers = 0;
                g3s = 1e6;
        }
	
	  
        double factorO2 = 1;
        double factorCar = 1;

        if (position.x!=0)
        {		
            if (position.z==0)
            {
                Cgr=(prevWal->Get(position.x-1,position.y).CarbonAgar+prevEnv->Get(position.x-1,position.y,position.z).Carbon
                     +prevWal->Get(position.x,position.y).CarbonAgar+prevEnv->Get(position.x,position.y,position.z).Carbon)/4.0;
                Wgr=(prevWal->Get(position.x-1,position.y).WasteAgar+prevEnv->Get(position.x-1,position.y,position.z).Waste
                     +prevWal->Get(position.x,position.y).WasteAgar+prevEnv->Get(position.x,position.y,position.z).Waste)/4.0;
                Ogr=(prevWal->Get(position.x-1,position.y).OxygenAgar+prevEnv->Get(position.x-1,position.y,position.z).Oxygen
                     +prevWal->Get(position.x,position.y).OxygenAgar+prevEnv->Get(position.x,position.y,position.z).Oxygen)/4.0;
                Agr=(prevWal->Get(position.x-1,position.y).AcetateAgar+prevEnv->Get(position.x-1,position.y,position.z).Acetate
                     +prevWal->Get(position.x,position.y).AcetateAgar+prevEnv->Get(position.x,position.y,position.z).Acetate)/4.0;
            }
            else
            {
                Cgr=(prevEnv->Get(position.x-1,position.y,position.z-1).Carbon+prevEnv->Get(position.x-1,position.y,position.z).Carbon
                     +prevEnv->Get(position.x,position.y,position.z-1).Carbon+prevEnv->Get(position.x,position.y,position.z).Carbon)/4.0;
                Wgr=(prevEnv->Get(position.x-1,position.y,position.z-1).Waste+prevEnv->Get(position.x-1,position.y,position.z).Waste
                     +prevEnv->Get(position.x,position.y,position.z-1).Waste+prevEnv->Get(position.x,position.y,position.z).Waste)/4.0;
                Ogr=(prevEnv->Get(position.x-1,position.y,position.z-1).Oxygen+prevEnv->Get(position.x-1,position.y,position.z).Oxygen
                     +prevEnv->Get(position.x,position.y,position.z-1).Oxygen+prevEnv->Get(position.x,position.y,position.z).Oxygen)/4.0;
                Agr=(prevEnv->Get(position.x-1,position.y,position.z-1).Acetate+prevEnv->Get(position.x-1,position.y,position.z).Acetate
                     +prevEnv->Get(position.x,position.y,position.z-1).Acetate+prevEnv->Get(position.x,position.y,position.z).Acetate)/4.0;
            }
        }
        
        double HAgr = acetateBufferSolve(Agr, prevEnv, position);

        double fCgr_aer = Cgr/(Cgr+KC_aer);
        double fCgr_ana = Cgr/(Cgr+KC_ana);
        double fWgr = Wgr/(Wgr+KW);
	    double fAgr = Agr/(Agr+KA);
        double fOgr = Ogr/(Ogr+KO);
        //double fOgr = 1;
	    //double fIgr = Agr/(Agr+KI);
        double fIgr = 0;//sqrt(HAgr/KI);

        double qC, GR1, pW, qW, GR2, qA, pA, qO, l1, l2, l3;
        
        if (crossFeeding)
        {
            qC = C_rate*fC*Density1ShiftP.Get(position);	// carbon consumption rate
            GR1 = maxGrowthRate1*max(0.0,fCgr_aer-Maintenance_rate/C_rate);	// Type 1 growth rate
            pW = Carbon_to_Waste*qC; // waste production rate

            qW = W_rate*fW*Density2ShiftP.Get(position);	// waste consumption rate
            GR2 = maxGrowthRate2*max(0.0,fWgr-Maintenance_rate/W_rate);	// Type 2 growth rate
        }
        else
        {
            l1 = max((maxGrowthRateCarAer+mgrCarAers)*fCgr_aer - mgrCarAers,0.0);
            l2 = max((maxGrowthRateCarAna+mgrCarAnas)*fCgr_ana - mgrCarAnas,0.0);
            l3 = max((maxGrowthRateAceAer+mgrAceAers)*fAgr - mgrAceAers,0.0);

            qO =(qO2Car*l1*fCgr_aer + qO2Ace*l3*(1-fCgr_aer))*fOgr*max(1-fIgr,0.0) + qO2Carnt*min((Cgr/g1s),1.0)*fCgr_aer*fOgr + qO2Acent*min((Agr/g3s),1.0)*(1-fCgr_aer)*fOgr;
            qC =  (qCarAer*l1*fCgr_aer*fOgr + qCarAna*l2*(1-fOgr))*max(1-fIgr,0.0) + qCarAernt*min((Cgr/g1s),1.0)*fCgr_aer*fOgr + qCarAnant*min((Cgr/g2s),1.0)*(1-fOgr);
            qA = (qAceAer*l3*(1-fCgr_aer)*fOgr)*max(1-fIgr,0.0) + qAceAernt*min((Agr/g3s),1.0)*(1-fCgr_aer)*fOgr;
            pA = (pAceAer*l1*fCgr_aer*fOgr + pAceAna*l2*(1-fOgr))*max(1-fIgr,0.0) + pAceAnant*min((Cgr/g2s),1.0)*(1-fOgr);
            
            	    
            Env->At(position).qCarmnt = qCarAernt*min((Cgr/g1s),1.0)*fCgr_aer*fOgr + qCarAnant*min((Cgr/g2s),1.0)*(1-fOgr);
            Env->At(position).qAcemnt = qAceAernt*min((Agr/g3s),1.0)*(1-fCgr_aer)*fOgr;
            Env->At(position).pAcemnt = pAceAnant*min((Cgr/g2s),1.0)*(1-fOgr);
            Env->At(position).qOxymnt = qO2Carnt*min((Cgr/g1s),1.0)*fCgr_aer*fOgr + qO2Acent*min((Agr/g3s),1.0)*(1-fCgr_aer)*fOgr;

            pW = Carbon_to_Waste*qC; // waste production rate
            qW = 0.0;	// waste consumption rate
            GR2 = maxGrowthRate2*max(0.0,fCgr_aer-Maintenance_rate/C_rate);	// Type 2 growth rate
        }

        l1 = max((maxGrowthRateCarAer+mgrCarAers)*fCgr_aer - mgrCarAers,0.0);
        l2 = max((maxGrowthRateCarAna+mgrCarAnas)*fCgr_ana - mgrCarAnas,0.0);
        l3 = max((maxGrowthRateAceAer+mgrAceAers)*fAgr - mgrAceAers,0.0);

        Env->At(position).GrCarAer = l1*fCgr_aer*fOgr*max(1-fIgr,0.0);
        Env->At(position).GrCarAna = l2*(1-fOgr)*max(1-fIgr,0.0);
        Env->At(position).GrAceAer = l3*(1-fCgr_aer)*fOgr*max(1-fIgr,0.0);
        Env->At(position).GrowthRate1 = Env->Get(position).GrCarAer + Env->Get(position).GrCarAna + Env->Get(position).GrAceAer;
        Env->At(position).GrowthRate2 = Env->At(position).GrowthRate1;

        // gradient of density
        double Dx, Dy, Dz;

        // gradient of C
        double Cx, Cy, Cz;

        // gradient of W
        double Wx, Wy, Wz;
        
        // second derivatives of C, O, and A
        double Cxx, Czz, Wxx, Wzz, Axx, Azz, Oxx, Ozz;
        
        MyAssert(position.x>2,"Need more boxes");

        // if we're at the bottom boundary
        if (position.z==0)
        {
            Dz = 0.0;
            //Dz =(DensityShiftP.Get(position.x,position.y,position.z+1) - WallDensityShiftP.Get(position.x,position.y))/(2.0*BoxLength);

            Wz = (prevEnv->Get(position.x,position.y,position.z+1).Waste - prevWal->Get(position.x,position.y).WasteAgar)/(2.0*BoxLength);
            Wzz = (prevEnv->Get(position.x,position.y,position.z+1).Waste + prevWal->Get(position.x,position.y).WasteAgar- 2.0*prevEnv->Get(position.x,position.y,position.z).Waste)/(BoxLength*BoxLength);

            Cz = (prevEnv->Get(position.x,position.y,position.z+1).Carbon - prevWal->Get(position.x,position.y).CarbonAgar)/(2.0*BoxLength); // BC = 1
            Czz = (prevEnv->Get(position.x,position.y,position.z+1).Carbon + prevWal->Get(position.x,position.y).CarbonAgar - 2.0*prevEnv->Get(position.x,position.y,position.z).Carbon)/(BoxLength*BoxLength);

            Azz = (prevEnv->Get(position.x,position.y,position.z+1).Acetate + prevWal->Get(position.x,position.y).AcetateAgar - 2.0*prevEnv->Get(position.x,position.y,position.z).Acetate)/(BoxLength*BoxLength);
            
            Ozz = (prevEnv->Get(position.x,position.y,position.z+1).Oxygen + prevWal->Get(position.x,position.y).OxygenAgar - 2.0*prevEnv->Get(position.x,position.y,position.z).Oxygen)/(BoxLength*BoxLength);

        }
        else
        {
            Dz = (DensityShiftP.Get(position.x,position.y,position.z+1) - DensityShiftP.Get(position.x,position.y,position.z-1))/(2.0*BoxLength);

            Wz = (prevEnv->Get(position.x,position.y,position.z+1).Waste - prevEnv->Get(position.x,position.y,position.z-1).Waste)/(2.0*BoxLength);
            Wzz = (prevEnv->Get(position.x,position.y,position.z+1).Waste + prevEnv->Get(position.x,position.y,position.z-1).Waste - 2.0*prevEnv->Get(position.x,position.y,position.z).Waste)/(BoxLength*BoxLength);

            Cz = (prevEnv->Get(position.x,position.y,position.z+1).Carbon - prevEnv->Get(position.x,position.y,position.z-1).Carbon)/(2.0*BoxLength);
            Czz = (prevEnv->Get(position.x,position.y,position.z+1).Carbon + prevEnv->Get(position.x,position.y,position.z-1).Carbon - 2.0*prevEnv->Get(position.x,position.y,position.z).Carbon)/(BoxLength*BoxLength);

			Azz = (prevEnv->Get(position.x,position.y,position.z+1).Acetate + prevEnv->Get(position.x,position.y,position.z-1).Acetate - 2.0*prevEnv->Get(position.x,position.y,position.z).Acetate)/(BoxLength*BoxLength);
            
            Ozz = (prevEnv->Get(position.x,position.y,position.z+1).Oxygen + prevEnv->Get(position.x,position.y,position.z-1).Oxygen - 2.0*prevEnv->Get(position.x,position.y,position.z).Oxygen)/(BoxLength*BoxLength);

        }

        Dx = (DensityShiftP.Get(position.x+1,position.y,position.z) - DensityShiftP.Get(position.x-1,position.y,position.z))/(2.0*BoxLength);
        Cx = (prevEnv->Get(position.x+1,position.y,position.z).Carbon - prevEnv->Get(position.x-1,position.y,position.z).Carbon)/(2.0*BoxLength);
        Cxx = (prevEnv->Get(position.x+1,position.y,position.z).Carbon + prevEnv->Get(position.x-1,position.y,position.z).Carbon - 2.0*prevEnv->Get(position.x,position.y,position.z).Carbon)/(BoxLength*BoxLength);

        Wx = (prevEnv->Get(position.x+1,position.y,position.z).Waste - prevEnv->Get(position.x-1,position.y,position.z).Waste)/(2.0*BoxLength);
        Wxx = (prevEnv->Get(position.x+1,position.y,position.z).Waste + prevEnv->Get(position.x-1,position.y,position.z).Waste - 2.0*prevEnv->Get(position.x,position.y,position.z).Waste)/(BoxLength*BoxLength);

		Axx = (prevEnv->Get(position.x+1,position.y,position.z).Acetate + prevEnv->Get(position.x-1,position.y,position.z).Acetate - 2.0*prevEnv->Get(position.x,position.y,position.z).Acetate)/(BoxLength*BoxLength);        
        
        Oxx = (prevEnv->Get(position.x+1,position.y,position.z).Oxygen + prevEnv->Get(position.x-1,position.y,position.z).Oxygen - 2.0*prevEnv->Get(position.x,position.y,position.z).Oxygen)/(BoxLength*BoxLength);

        // Calculate new C concentration
        double Cnew, CAcenew, Onew;
        if (crossFeeding)
        {
            Cnew = prevEnv->Get(position).Carbon + (DiffColonyCarbon*(Cxx + Czz ) - qC*0.25*DensityShiftP.Get(position))*Cdt;
        }
        else
        {
            Cnew = prevEnv->Get(position).Carbon + (DiffColonyCarbon*(Cxx + Czz ) - qC*0.25*DensityShiftP.Get(position))*Cdt;
	    CAcenew = prevEnv->Get(position).Acetate + (DiffColonyAce*(Axx + Azz) + pA*0.25*DensityShiftP.Get(position) - qA*0.25*DensityShiftP.Get(position))*Cdt;
            Onew = prevEnv->Get(position).Oxygen + (DiffColonyO2*(Oxx + Ozz) - qO*0.25*DensityShiftP.Get(position))*Cdt;
        }
        Env->At(position).Carbon = max(0.0,min(Cnew,maxCarbon));	// for stability
        Env->At(position).Oxygen = max(0.0,min(Onew,maxO2));
        if (NutrientGSI==1)
        {
            prevEnv->At(position).Carbon = max(0.0,min(Cnew,maxCarbon));	// for stability
            prevEnv->At(position).Oxygen = max(0.0,min(Onew,maxO2));
        }
		Env->At(position).Acetate = max(0.0, CAcenew);	// for stability
		if (NutrientGSI == 1)
		{
			prevEnv->At(position).Acetate = max(0.0, CAcenew);	// for stability
		}
        
        double Wnew = prevEnv->Get(position).Waste + (DiffColony*(Wxx + Wzz ) + pW*DensityShiftP.Get(position) - qW*DensityShiftP.Get(position))*Cdt;
        Env->At(position).Waste = max(0.0,Wnew);	// for stability
        if (NutrientGSI==1)
        {
            prevEnv->At(position).Waste = max(0.0,Wnew);// for stability
        }

    }
    // Exterior Colony Surface
    else if (insideColonyDen.Get(position)>0.3 && position.z>0)
    {
        int ix = position.x;
        int iy = position.y;
        int iz = position.z;

        double Cnew = (prevEnv->Get(ix+1,iy,iz).Carbon*(insideColonyDen.Get(ix+1,iy,iz)>0.5)+prevEnv->Get(ix-1,iy,iz).Carbon*(insideColonyDen.Get(ix-1,iy,iz)>0.5)+prevEnv->Get(ix,iy,iz+1).Carbon*(insideColonyDen.Get(ix,iy,iz+1)>0.5)+prevEnv->Get(ix,iy,iz-1).Carbon*(insideColonyDen.Get(ix,iy,iz-1)>0.5))/((insideColonyDen.Get(ix+1,iy,iz)>0.5)+(insideColonyDen.Get(ix-1,iy,iz)>0.5)+(insideColonyDen.Get(ix,iy,iz+1)>0.5)+(insideColonyDen.Get(ix,iy,iz-1)>0.5));

        double Wnew = (prevEnv->Get(ix+1,iy,iz).Waste*(insideColonyDen.Get(ix+1,iy,iz)>0.5)+prevEnv->Get(ix-1,iy,iz).Waste*(insideColonyDen.Get(ix-1,iy,iz)>0.5)+prevEnv->Get(ix,iy,iz+1).Waste*(insideColonyDen.Get(ix,iy,iz+1)>0.5)+prevEnv->Get(ix,iy,iz-1).Waste*(insideColonyDen.Get(ix,iy,iz-1)>0.5))/((insideColonyDen.Get(ix+1,iy,iz)>0.5)+(insideColonyDen.Get(ix-1,iy,iz)>0.5)+(insideColonyDen.Get(ix,iy,iz+1)>0.5)+(insideColonyDen.Get(ix,iy,iz-1)>0.5));

		double CAcenew = (prevEnv->Get(ix+1,iy,iz).Acetate*(insideColonyDen.Get(ix+1,iy,iz)>0.5)+prevEnv->Get(ix-1,iy,iz).Acetate*(insideColonyDen.Get(ix-1,iy,iz)>0.5)+prevEnv->Get(ix,iy,iz+1).Acetate*(insideColonyDen.Get(ix,iy,iz+1)>0.5)+prevEnv->Get(ix,iy,iz-1).Acetate*(insideColonyDen.Get(ix,iy,iz-1)>0.5))/((insideColonyDen.Get(ix+1,iy,iz)>0.5)+(insideColonyDen.Get(ix-1,iy,iz)>0.5)+(insideColonyDen.Get(ix,iy,iz+1)>0.5)+(insideColonyDen.Get(ix,iy,iz-1)>0.5));;
        
        double Onew = maxO2;

        Env->At(position).Carbon = max(0.0,min(Cnew, maxCarbon));
        Env->At(position).Oxygen = max(0.0,min(Onew, maxO2));
        if (NutrientGSI==1)
        {
            prevEnv->At(position).Carbon = max(0.0,min(Cnew,maxCarbon));	// for stability
            prevEnv->At(position).Oxygen = max(0.0,min(Onew,maxO2));
        }
		Env->At(position).Acetate = max(0.0, CAcenew);   // for stability
		if (NutrientGSI == 1)
		{
			prevEnv->At(position).Acetate = max(0.0, CAcenew);  // for stability
		}
        Env->At(position).Waste = max(0.0,Wnew);	// for stability
        if (NutrientGSI==1)
        {
            prevEnv->At(position).Waste = max(0.0,Wnew);// for stability
        }
    }
    // Colony Wall Interaction
    else if (insideColonyDen.Get(position)>0.3 && position.z==0)
    {
        int ix = position.x;
        int iy = position.y;
        int iz = position.z;

        double Cnew = (prevEnv->Get(ix+1,iy,iz).Carbon*(insideColonyDen.Get(ix+1,iy,iz)>0.5)+prevEnv->Get(ix-1,iy,iz).Carbon*(insideColonyDen.Get(ix-1,iy,iz)>0.5)+prevEnv->Get(ix,iy,iz+1).Carbon*(insideColonyDen.Get(ix,iy,iz+1)>0.5))/((insideColonyDen.Get(ix+1,iy,iz)>0.5)+(insideColonyDen.Get(ix-1,iy,iz)>0.5)+(insideColonyDen.Get(ix,iy,iz+1)>0.5));

        double Wnew = (prevEnv->Get(ix+1,iy,iz).Waste*(insideColonyDen.Get(ix+1,iy,iz)>0.5)+prevEnv->Get(ix-1,iy,iz).Waste*(insideColonyDen.Get(ix-1,iy,iz)>0.5)+prevEnv->Get(ix,iy,iz+1).Waste*(insideColonyDen.Get(ix,iy,iz+1)>0.5))/((insideColonyDen.Get(ix+1,iy,iz)>0.5)+(insideColonyDen.Get(ix-1,iy,iz)>0.5)+(insideColonyDen.Get(ix,iy,iz+1)>0.5));

		double CAcenew = (prevEnv->Get(ix+1,iy,iz).Acetate*(insideColonyDen.Get(ix+1,iy,iz)>0.5)+prevEnv->Get(ix-1,iy,iz).Acetate*(insideColonyDen.Get(ix-1,iy,iz)>0.5)+prevEnv->Get(ix,iy,iz+1).Acetate*(insideColonyDen.Get(ix,iy,iz+1)>0.5))/((insideColonyDen.Get(ix+1,iy,iz)>0.5)+(insideColonyDen.Get(ix-1,iy,iz)>0.5)+(insideColonyDen.Get(ix,iy,iz+1)>0.5));
        
        Env->At(position).Carbon = max(0.0,min(Cnew, maxCarbon));    // for stability
        
        double Onew = maxO2;

        Env->At(position).Oxygen = max(0.0,min(Onew, maxO2));
        if (NutrientGSI==1)
        {
            prevEnv->At(position).Carbon = max(0.0,min(Cnew,maxCarbon));	// for stability
            prevEnv->At(position).Oxygen = max(0.0,min(Onew,maxO2));
        }
		Env->At(position).Acetate = max(0.0, CAcenew);    // for stability
		if (NutrientGSI == 1)
		{
			prevEnv->At(position).Acetate = max(0.0, CAcenew);	// for stability
		}
        Env->At(position).Waste = max(0.0,Wnew);	// for stability
        if (NutrientGSI==1)
        {
            prevEnv->At(position).Waste = max(0.0,Wnew);// for stability
        }
    }
}

// Colony Agar interaction surface
void UpdateWallMultigrid_OxygenAcetate(Array3D<LocalEnv>* Env, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wal, Array2D<LocalAga>** prevWal, Array2D<double>& WallDensityShiftP, Array2D<double>& WallDensity1ShiftP, Array2D<double>& WallDensity2ShiftP, IntCoord2D position, int level, Array2D<double>& Height, Array3D<double>& insideColonyDen)
{
    //if ((position.x>=minX) && (position.x<=maxX) && (position.y>=minY) && (position.y<=maxY))
    IntCoord positionWall3D;
    positionWall3D.x = position.x;
    positionWall3D.y = position.y;
    positionWall3D.z = 0;
    // Level 0 and within the Colony
    if (level==0 && insideColonyDen.Get(positionWall3D)>0.5) //(level==0 && WallDensityShiftP.Get(position) > 0.8 && Height.At(position) > cellRadius+1e-7)
    {
        
        double C = prevWal[level]->Get(position).CarbonAgar;
        double W = prevWal[level]->Get(position).WasteAgar;
		double O = prevWal[level]->Get(position).OxygenAgar;
		double A = prevWal[level]->Get(position).AcetateAgar;
        
        
        // calculate growth rate and consumption rates based on current concentrations
        double fC = C/(C+KC_aer);	// carbon transporter rate
        double fW = W/(W+KW);	// waste transporter rate
	    double fO = O / (O + KO);	// oxygen transporter rate
	    double fA = A / (A + KA);	// acetate transporter rate
        double fI = 0;
              
        switch (InterfaceCondition)
        {
            case 1:
            {
                double Cnew = (DiffColonyCarbon*Env->Get(positionWall3D).Carbon + DiffAgarCarbon*prevAga[level]->Get(positionWall3D).CarbonAgar)/(DiffColonyCarbon+DiffAgarCarbon);
                double Onew = (DiffColonyO2*Env->Get(positionWall3D).Oxygen + DiffAgarO2*prevAga[level]->Get(positionWall3D).OxygenAgar)/(DiffColonyO2+DiffAgarO2);
                
                Wal[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon)); // for stability
                Wal[level]->At(position).OxygenAgar = max(0.0,min(Onew,maxO2));
                if (NutrientGSI==1)
                {
                    prevWal[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon)); // for stability
                    prevWal[level]->At(position).OxygenAgar = max(0.0,min(Onew,maxO2));
                }
				double CAcenew = (DiffColonyAce*Env->Get(positionWall3D).Acetate + DiffAgarAce * prevAga[level]->Get(positionWall3D).AcetateAgar) / (DiffColonyAce + DiffAgarAce);
				Wal[level]->At(position).AcetateAgar = max(0.0, CAcenew); // for stability
				if (NutrientGSI == 1)
				{
					prevWal[level]->At(position).AcetateAgar = max(0.0, CAcenew); // for stability
				}
				double Wnew = (DiffColony*Env->Get(positionWall3D).Waste + DiffAgar*prevAga[level]->Get(positionWall3D).WasteAgar)/(DiffColony+DiffAgar);
                Wal[level]->At(position).WasteAgar = max(0.0,Wnew); // for stability
                if (NutrientGSI==1)
                {
                    prevWal[level]->At(position).WasteAgar = max(0.0,Wnew); // for stability
                }
                break;
            }            
            default:
            {
                printf("Wrong switch on interface condition!");
            }
        }
    }
    // Level 1, Level 2, Level 3, Level 0 Outside Colony
    else
    {
        double Cxx, Czz, Wxx, Wzz, Axx, Azz;
        if (level>0 && (position.x>=BoxX/4-1) && position.x<=BoxX*3/4)
        {
            if (position.x==BoxX/4-1)
            {
                Wxx = (prevWal[level]->Get(position.x-1,position.y).WasteAgar+prevWal[level-1]->Get(0,levelM1y(position.y)).WasteAgar-2*prevWal[level]->Get(position.x,position.y).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
                Cxx = (prevWal[level]->Get(position.x-1,position.y).CarbonAgar+prevWal[level-1]->Get(0,levelM1y(position.y)).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
				Axx = (prevWal[level]->Get(position.x-1,position.y).AcetateAgar+prevWal[level-1]->Get(0,levelM1y(position.y)).AcetateAgar-2*prevWal[level]->Get(position.x,position.y).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
            else if (position.x==BoxX*3/4)
            {
                Wxx = (prevWal[level]->Get(position.x+1,position.y).WasteAgar+prevWal[level-1]->Get(BoxX-2,levelM1y(position.y)).WasteAgar-2*prevWal[level]->Get(position.x,position.y).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
                Cxx = (prevWal[level]->Get(position.x+1,position.y).CarbonAgar+prevWal[level-1]->Get(BoxX-2,levelM1y(position.y)).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
				Axx = (prevWal[level]->Get(position.x+1,position.y).AcetateAgar+prevWal[level-1]->Get(BoxX-2,levelM1y(position.y)).AcetateAgar-2*prevWal[level]->Get(position.x,position.y).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
        }
        else if (level<maxLevelsMG-1 && (position.x==0 || position.x==BoxX-1))
        {
            int ly=position.y%2;
            if (position.x==0)
            {
                Wxx = (prevWal[level]->Get(position.x+2,position.y).WasteAgar+1.0/(ly+1)*prevWal[level+1]->Get(BoxX/4-1,levelP1y(position.y)).WasteAgar+ly*1.0/(ly+1)*prevWal[level+1]->Get(BoxX/4-1,(levelP1y(position.y)+1)%BoxY).WasteAgar-2*prevWal[level]->Get(position.x,position.y).WasteAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
                Cxx = (prevWal[level]->Get(position.x+2,position.y).CarbonAgar+1.0/(ly+1)*prevWal[level+1]->Get(BoxX/4-1,levelP1y(position.y)).CarbonAgar+ly*1.0/(ly+1)*prevWal[level+1]->Get(BoxX/4-1,(levelP1y(position.y)+1)%BoxY).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
				Axx = (prevWal[level]->Get(position.x+2,position.y).AcetateAgar+1.0/(ly+1)*prevWal[level+1]->Get(BoxX/4-1,levelP1y(position.y)).AcetateAgar+ly*1.0/(ly+1)*prevWal[level+1]->Get(BoxX/4-1,(levelP1y(position.y)+1)%BoxY).AcetateAgar-2*prevWal[level]->Get(position.x,position.y).AcetateAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
            }
            else
            {
                Wxx = (prevWal[level]->Get(position.x-1,position.y).WasteAgar+1.0/(ly+1)*prevWal[level+1]->Get(BoxX*3/4,levelP1y(position.y)).WasteAgar+ly*1.0/(ly+1)*prevWal[level+1]->Get(BoxX*3/4,(levelP1y(position.y)+1)%BoxY).WasteAgar-2*prevWal[level]->Get(position.x,position.y).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
                Cxx = (prevWal[level]->Get(position.x-1,position.y).CarbonAgar+1.0/(ly+1)*prevWal[level+1]->Get(BoxX*3/4,levelP1y(position.y)).CarbonAgar+ly*1.0/(ly+1)*prevWal[level+1]->Get(BoxX*3/4,(levelP1y(position.y)+1)%BoxY).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
				Axx = (prevWal[level]->Get(position.x-1,position.y).AcetateAgar+1.0/(ly+1)*prevWal[level+1]->Get(BoxX*3/4,levelP1y(position.y)).AcetateAgar+ly*1.0/(ly+1)*prevWal[level+1]->Get(BoxX*3/4,(levelP1y(position.y)+1)%BoxY).AcetateAgar-2*prevWal[level]->Get(position.x,position.y).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
        }
        //Level 3 Sides (Delta S)
        else if (level==maxLevelsMG-1 && (position.x==0 || position.x==BoxX-1))
        {
            if (position.x==0)
            {
                Wxx = (prevWal[level]->Get(position.x+1,position.y).WasteAgar-1.0*prevWal[level]->Get(position.x,position.y).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
                Cxx = (prevWal[level]->Get(position.x+1,position.y).CarbonAgar-1.0*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
				Axx = (prevWal[level]->Get(position.x+1,position.y).AcetateAgar-1.0*prevWal[level]->Get(position.x,position.y).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
            else if (position.x==(BoxX-1))
            {
                Wxx = (prevWal[level]->Get(position.x-1,position.y).WasteAgar-1.0*prevWal[level]->Get(position.x,position.y).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
                Cxx = (prevWal[level]->Get(position.x-1,position.y).CarbonAgar-1.0*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
				Axx = (prevWal[level]->Get(position.x-1,position.y).AcetateAgar-1.0*prevWal[level]->Get(position.x,position.y).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
            }
        }
        else
        {
            Wxx = (prevWal[level]->Get(position.x-1,position.y).WasteAgar+prevWal[level]->Get(position.x+1,position.y).WasteAgar-2*prevWal[level]->Get(position.x,position.y).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Cxx = (prevWal[level]->Get(position.x-1,position.y).CarbonAgar+prevWal[level]->Get(position.x+1,position.y).CarbonAgar-2*prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
			Axx = (prevWal[level]->Get(position.x-1,position.y).AcetateAgar+prevWal[level]->Get(position.x+1,position.y).AcetateAgar-2*prevWal[level]->Get(position.x,position.y).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }


        ///////////////////////////////////////////


        Wzz = 2.0*(prevAga[level]->Get(position.x,position.y,0).WasteAgar-prevWal[level]->Get(position.x,position.y).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Czz = 2.0*(prevAga[level]->Get(position.x,position.y,0).CarbonAgar-prevWal[level]->Get(position.x,position.y).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
		Azz = 2.0*(prevAga[level]->Get(position.x,position.y,0).AcetateAgar-prevWal[level]->Get(position.x,position.y).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
        
        double Cnew = prevWal[level]->Get(position).CarbonAgar + DiffAgarCarbon*(Cxx + Czz)*Cdt;
        double Onew = maxO2;
        Wal[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon));	// for stability
        Wal[level]->At(position).OxygenAgar = max(0.0,min(Onew,maxO2));
        if (NutrientGSI==1)
        {
            prevWal[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon));	// for stability
            prevWal[level]->At(position).OxygenAgar = max(0.0,min(Onew,maxO2));
        }
        
	double CAcenew = prevWal[level]->Get(position).AcetateAgar +DiffAgarAce * (Axx + Azz)*Cdt;
	Wal[level]->At(position).AcetateAgar = max(0.0, CAcenew);	// for stability
	if (NutrientGSI == 1)
	{
		prevWal[level]->At(position).AcetateAgar = max(0.0, CAcenew);	// for stability
	}

        double Wnew = prevWal[level]->Get(position).WasteAgar + DiffAgar*(Wxx + Wzz)*Cdt;
        Wal[level]->At(position).WasteAgar = max(0.0,Wnew);	// for stability
        if (NutrientGSI==1)
        {
            prevWal[level]->At(position).WasteAgar = max(0.0,Wnew);	// for stability
        }

    }

}
// Agar area
void UpdateAgarMultigrid_OxygenAcetate(Array3D<LocalAga>** Aga, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wal, IntCoord position, int level)
{

    // second derivatives of C, W, A, O
    double Cxx, Czz, Wxx, Wzz, Axx, Azz, Oxx, Ozz;
    
    
    // if we're at the upper boundary
    if (position.z==0)
    {

        Wzz = (prevAga[level]->Get(position.x,position.y,position.z+1).WasteAgar+ Wal[level]->Get(position.x,position.y).WasteAgar- 2.0*prevAga[level]->Get(position.x,position.y,position.z).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Czz = (prevAga[level]->Get(position.x,position.y,position.z+1).CarbonAgar + Wal[level]->Get(position.x,position.y).CarbonAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
		Azz = (prevAga[level]->Get(position.x,position.y,position.z+1).AcetateAgar + Wal[level]->Get(position.x,position.y).AcetateAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Ozz = (prevAga[level]->Get(position.x,position.y,position.z+1).OxygenAgar + Wal[level]->Get(position.x,position.y).OxygenAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).OxygenAgar)/(BoxLength*BoxLength*pow(4.0,level));

    }
    // Level 1 - maxLevelsMG at the boundary between Levels
    else if (level>0 && position.x>=BoxX/4 && position.x<3*BoxX/4 && position.z<=BoxZAgar/2)
    {
        if (position.z==BoxZAgar/2)
        {
            Wzz = (prevAga[level]->Get(position.x,position.y,position.z+1).WasteAgar+ prevAga[level-1]->Get(levelM1x(position.x),levelM1y(position.y),BoxZAgar-1).WasteAgar- 2.0*prevAga[level]->Get(position.x,position.y,position.z).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Czz = (prevAga[level]->Get(position.x,position.y,position.z+1).CarbonAgar + prevAga[level-1]->Get(levelM1x(position.x),levelM1y(position.y),BoxZAgar-1).CarbonAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
			Azz = (prevAga[level]->Get(position.x,position.y,position.z+1).AcetateAgar + prevAga[level-1]->Get(levelM1x(position.x),levelM1y(position.y),BoxZAgar-1).AcetateAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Ozz = (prevAga[level]->Get(position.x,position.y,position.z+1).OxygenAgar + prevAga[level-1]->Get(levelM1x(position.x),levelM1y(position.y),BoxZAgar-1).OxygenAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).OxygenAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
    }
    // Bottom boundary of Level 0 - (maxLevelsMG - 1)
    else if (level<maxLevelsMG-1 && position.z==BoxZAgar-1)
    {
        int lx=position.x%2;
        Wzz = (prevAga[level]->Get(position.x,position.y,position.z-2).WasteAgar+ 1.0/(lx+1)*prevAga[level+1]->Get(levelP1x(position.x),levelP1y(position.y),BoxZAgar/2).WasteAgar+ lx*1.0/(lx+1)*prevAga[level+1]->Get((levelP1x(position.x)+1)%BoxX,levelP1y(position.y),BoxZAgar/2).WasteAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).WasteAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
        Czz = (prevAga[level]->Get(position.x,position.y,position.z-2).CarbonAgar + 1.0/(lx+1)*prevAga[level+1]->Get(levelP1x(position.x),levelP1y(position.y),BoxZAgar/2).CarbonAgar + lx*1.0/(lx+1)*prevAga[level+1]->Get((levelP1x(position.x)+1)%BoxX,levelP1y(position.y),BoxZAgar/2).CarbonAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
		Azz = (prevAga[level]->Get(position.x,position.y,position.z-2).AcetateAgar + 1.0/(lx+1)*prevAga[level+1]->Get(levelP1x(position.x),levelP1y(position.y),BoxZAgar/2).AcetateAgar + lx*1.0/(lx+1)*prevAga[level+1]->Get((levelP1x(position.x)+1)%BoxX,levelP1y(position.y),BoxZAgar/2).AcetateAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
        Ozz = (prevAga[level]->Get(position.x,position.y,position.z-2).OxygenAgar + 1.0/(lx+1)*prevAga[level+1]->Get(levelP1x(position.x),levelP1y(position.y),BoxZAgar/2).OxygenAgar + lx*1.0/(lx+1)*prevAga[level+1]->Get((levelP1x(position.x)+1)%BoxX,levelP1y(position.y),BoxZAgar/2).OxygenAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).OxygenAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
    }
    // Bottom boundary (Delta b) of Level maxLevelsMG
    else if (level==maxLevelsMG-1 && position.z==BoxZAgar-1)
    {
        Wzz = 2.0*(prevAga[level]->Get(position.x,position.y,position.z-1).WasteAgar-prevAga[level]->Get(position.x,position.y,position.z).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Czz = 2.0*(prevAga[level]->Get(position.x,position.y,position.z-1).CarbonAgar-prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
		Azz = 2.0*(prevAga[level]->Get(position.x,position.y,position.z-1).AcetateAgar-prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Ozz = 2.0*(prevAga[level]->Get(position.x,position.y,position.z-1).OxygenAgar-prevAga[level]->Get(position.x,position.y,position.z).OxygenAgar)/(BoxLength*BoxLength*pow(4.0,level));
    }
    else
    {
        Wzz = (prevAga[level]->Get(position.x,position.y,position.z+1).WasteAgar + prevAga[level]->Get(position.x,position.y,position.z-1).WasteAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Czz = (prevAga[level]->Get(position.x,position.y,position.z+1).CarbonAgar + prevAga[level]->Get(position.x,position.y,position.z-1).CarbonAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
		Azz = (prevAga[level]->Get(position.x,position.y,position.z+1).AcetateAgar + prevAga[level]->Get(position.x,position.y,position.z-1).AcetateAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Ozz = (prevAga[level]->Get(position.x,position.y,position.z+1).OxygenAgar + prevAga[level]->Get(position.x,position.y,position.z-1).OxygenAgar - 2.0*prevAga[level]->Get(position.x,position.y,position.z).OxygenAgar)/(BoxLength*BoxLength*pow(4.0,level));
    }

    ////////////////////////////////////////////////

    if (level>0 && (position.x>=BoxX/4-1) && position.x<=BoxX*3/4 && position.z<=BoxZAgar/2-1)
    {
        if (position.x==BoxX/4-1)
        {
            Wxx = (prevAga[level]->Get(position.x-1,position.y,position.z).WasteAgar+prevAga[level-1]->Get(0,levelM1y(position.y),levelM1z(position.z)).WasteAgar-2*prevAga[level]->Get(position.x,position.y,position.z).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Cxx = (prevAga[level]->Get(position.x-1,position.y,position.z).CarbonAgar+prevAga[level-1]->Get(0,levelM1y(position.y),levelM1z(position.z)).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
			Axx = (prevAga[level]->Get(position.x-1,position.y,position.z).AcetateAgar+prevAga[level-1]->Get(0,levelM1y(position.y),levelM1z(position.z)).AcetateAgar-2*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Oxx = (prevAga[level]->Get(position.x-1,position.y,position.z).OxygenAgar+prevAga[level-1]->Get(0,levelM1y(position.y),levelM1z(position.z)).OxygenAgar-2*prevAga[level]->Get(position.x,position.y,position.z).OxygenAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
        else if (position.x==BoxX*3/4)
        {
            Wxx = (prevAga[level]->Get(position.x+1,position.y,position.z).WasteAgar+prevAga[level-1]->Get(BoxX-2,levelM1y(position.y),levelM1z(position.z)).WasteAgar-2*prevAga[level]->Get(position.x,position.y,position.z).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Cxx = (prevAga[level]->Get(position.x+1,position.y,position.z).CarbonAgar+prevAga[level-1]->Get(BoxX-2,levelM1y(position.y),levelM1z(position.z)).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
			Axx = (prevAga[level]->Get(position.x+1,position.y,position.z).AcetateAgar+prevAga[level-1]->Get(BoxX-2,levelM1y(position.y),levelM1z(position.z)).AcetateAgar-2*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Oxx = (prevAga[level]->Get(position.x+1,position.y,position.z).OxygenAgar+prevAga[level-1]->Get(BoxX-2,levelM1y(position.y),levelM1z(position.z)).OxygenAgar-2*prevAga[level]->Get(position.x,position.y,position.z).OxygenAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
    }
    else if (level<maxLevelsMG-1 && (position.x==0 || position.x==BoxX-1))
    {
        int lz=(position.z+1)%2;
        if (position.x==0)
        {
            Wxx = (prevAga[level]->Get(position.x+2,position.y,position.z).WasteAgar+1.0/(lz+1)*prevAga[level+1]->Get(BoxX/4-1,levelP1y(position.y),levelP1z(position.z)).WasteAgar+lz*1.0/(lz+1)*prevAga[level+1]->Get(BoxX/4-1,levelP1y(position.y),(levelP1z(position.z)+1)%BoxZAgar).WasteAgar-2*prevAga[level]->Get(position.x,position.y,position.z).WasteAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
            Cxx = (prevAga[level]->Get(position.x+2,position.y,position.z).CarbonAgar+1.0/(lz+1)*prevAga[level+1]->Get(BoxX/4-1,levelP1y(position.y),levelP1z(position.z)).CarbonAgar+lz*1.0/(lz+1)*prevAga[level+1]->Get(BoxX/4-1,levelP1y(position.y),(levelP1z(position.z)+1)%BoxZAgar).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
			Axx = (prevAga[level]->Get(position.x+2,position.y,position.z).AcetateAgar+1.0/(lz+1)*prevAga[level+1]->Get(BoxX/4-1,levelP1y(position.y),levelP1z(position.z)).AcetateAgar+lz*1.0/(lz+1)*prevAga[level+1]->Get(BoxX/4-1,levelP1y(position.y),(levelP1z(position.z)+1)%BoxZAgar).AcetateAgar-2*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
            Oxx = (prevAga[level]->Get(position.x+2,position.y,position.z).OxygenAgar+1.0/(lz+1)*prevAga[level+1]->Get(BoxX/4-1,levelP1y(position.y),levelP1z(position.z)).OxygenAgar+lz*1.0/(lz+1)*prevAga[level+1]->Get(BoxX/4-1,levelP1y(position.y),(levelP1z(position.z)+1)%BoxZAgar).OxygenAgar-2*prevAga[level]->Get(position.x,position.y,position.z).OxygenAgar)/4/(BoxLength*BoxLength*pow(4.0,level));
        }
        else
        {
            Wxx = (prevAga[level]->Get(position.x-1,position.y,position.z).WasteAgar+1.0/(lz+1)*prevAga[level+1]->Get(BoxX*3/4,levelP1y(position.y),levelP1z(position.z)).WasteAgar+lz*1.0/(lz+1)*prevAga[level+1]->Get(BoxX*3/4,levelP1y(position.y),(levelP1z(position.z)+1)%BoxZAgar).WasteAgar-2*prevAga[level]->Get(position.x,position.y,position.z).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Cxx = (prevAga[level]->Get(position.x-1,position.y,position.z).CarbonAgar+1.0/(lz+1)*prevAga[level+1]->Get(BoxX*3/4,levelP1y(position.y),levelP1z(position.z)).CarbonAgar+lz*1.0/(lz+1)*prevAga[level+1]->Get(BoxX*3/4,levelP1y(position.y),(levelP1z(position.z)+1)%BoxZAgar).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
			Axx = (prevAga[level]->Get(position.x-1,position.y,position.z).AcetateAgar+1.0/(lz+1)*prevAga[level+1]->Get(BoxX*3/4,levelP1y(position.y),levelP1z(position.z)).AcetateAgar+lz*1.0/(lz+1)*prevAga[level+1]->Get(BoxX*3/4,levelP1y(position.y),(levelP1z(position.z)+1)%BoxZAgar).AcetateAgar-2*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Oxx = (prevAga[level]->Get(position.x-1,position.y,position.z).OxygenAgar+1.0/(lz+1)*prevAga[level+1]->Get(BoxX*3/4,levelP1y(position.y),levelP1z(position.z)).OxygenAgar+lz*1.0/(lz+1)*prevAga[level+1]->Get(BoxX*3/4,levelP1y(position.y),(levelP1z(position.z)+1)%BoxZAgar).OxygenAgar-2*prevAga[level]->Get(position.x,position.y,position.z).OxygenAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
    }
    else if (level==maxLevelsMG-1 && (position.x==0 || position.x==BoxX-1))
    {
        if (position.x==0)
        {
            Wxx = (prevAga[level]->Get(position.x+1,position.y,position.z).WasteAgar-1*prevAga[level]->Get(position.x,position.y,position.z).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Cxx = (prevAga[level]->Get(position.x+1,position.y,position.z).CarbonAgar-1.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
			Axx = (prevAga[level]->Get(position.x+1,position.y,position.z).AcetateAgar-1.0*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Oxx = (prevAga[level]->Get(position.x+1,position.y,position.z).OxygenAgar+maxO2-2*prevAga[level]->Get(position.x,position.y,position.z).OxygenAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
        else if (position.x==(BoxX-1))
        {
            Wxx = (prevAga[level]->Get(position.x-1,position.y,position.z).WasteAgar-1*prevAga[level]->Get(position.x,position.y,position.z).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Cxx = (prevAga[level]->Get(position.x-1,position.y,position.z).CarbonAgar-1.0*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
			Axx = (prevAga[level]->Get(position.x-1,position.y,position.z).AcetateAgar-1.0*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Oxx = (maxO2+prevAga[level]->Get(position.x-1,position.y,position.z).OxygenAgar-2*prevAga[level]->Get(position.x,position.y,position.z).OxygenAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
    }
    else
    {
        Wxx = (prevAga[level]->Get(position.x-1,position.y,position.z).WasteAgar+prevAga[level]->Get(position.x+1,position.y,position.z).WasteAgar-2*prevAga[level]->Get(position.x,position.y,position.z).WasteAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Cxx = (prevAga[level]->Get(position.x-1,position.y,position.z).CarbonAgar+prevAga[level]->Get(position.x+1,position.y,position.z).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
		Axx = (prevAga[level]->Get(position.x-1,position.y,position.z).AcetateAgar+prevAga[level]->Get(position.x+1,position.y,position.z).AcetateAgar-2*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Oxx = (prevAga[level]->Get(position.x-1,position.y,position.z).OxygenAgar+prevAga[level]->Get(position.x+1,position.y,position.z).OxygenAgar-2*prevAga[level]->Get(position.x,position.y,position.z).OxygenAgar)/(BoxLength*BoxLength*pow(4.0,level));
    }


    ///////////////////////////////////////////

    // Calculate new concentrations for C, O, A, W
    double Cnew = prevAga[level]->Get(position).CarbonAgar + DiffAgarCarbon*(Cxx + Czz)*Cdt;
    double Onew = prevAga[level]->Get(position).OxygenAgar + DiffAgarO2*(Oxx + Ozz)*Cdt;
    Aga[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon));	// for stability
    Aga[level]->At(position).OxygenAgar = max(0.0,min(Onew,maxO2));
    if (NutrientGSI==1)
    {
        prevAga[level]->At(position).CarbonAgar = max(0.0,min(Cnew,maxCarbon));	// for stability
        prevAga[level]->At(position).OxygenAgar = max(0.0,min(Onew,maxO2));
    }

	double CAcenew = prevAga[level]->Get(position).AcetateAgar + DiffAgarAce * (Axx + Azz)*Cdt;
	Aga[level]->At(position).AcetateAgar = max(0.0, CAcenew);	// for stability
	if (NutrientGSI == 1)
	{
		prevAga[level]->At(position).AcetateAgar = max(0.0, CAcenew);	// for stability
	}

    double Wnew = prevAga[level]->Get(position).WasteAgar + DiffAgar*(Wxx + Wzz)*Cdt;
    Aga[level]->At(position).WasteAgar = max(0.0,Wnew);	// for stability
    if (NutrientGSI==1)
    {
        prevAga[level]->At(position).WasteAgar = max(0.0,Wnew);	// for stability
    }

}

double acetateBufferSolve(double Cac, Array3D<LocalEnv>* prevEnv, IntCoord position)
{
    double aH = -(Cb*Kb)/(pow(10,-4)+Kb);
    // Initial Conditions x = OH-, y = A-, z = B-
    double x = 1e-7, y = Cac, z = Cb;
    double s[3];

    for (int it = 0; it < 100; it++)
    {
        double J[][4] = {{2*x + y + z + (Kw/c1) + aH, x, x, -(pow(x,2) + x*(y + z + (Kw/c1) + aH) - Kw)}, 
                         {y, 2*y + x + z + Ka + aH, y, -(pow(y,2) + y*(x + z + (Ka) + aH) - Ka*Cac)}, 
                         {z, z, 2*z + x + y + Kb + aH, -(pow(z,2) + z*(x + y + (Kb) + aH) - Kb*Cb)}};
        
        double m1 = J[1][0]/J[0][0], m2 = J[2][0]/J[0][0];
        for (int i=0;i<4;i++)
        {
            J[1][i] -= m1*J[0][i];
            J[2][i] -= m2*J[0][i];
        }
        double m3 = J[2][1]/J[1][1];
        for (int i=0;i<4;i++)
        {
            J[2][i] -= m3*J[1][i];
        }
        
        s[2] = J[2][3]/J[2][2];
        s[1] = J[1][3]/J[1][1] - (J[1][2]/J[1][1])*s[2];
        s[0] = J[0][3]/J[0][0] - (J[0][2]/J[0][0])*s[2] - (J[0][1]/J[0][0])*s[1];

        double max = 0;
        for (int i = 0;i<3;i++)
        {
            if (fabs(s[i]) > max)
                max = fabs(s[i]);
        }
        
        x += s[0];
        y += s[1];
        z += s[2];

        if (max < pow(10,-7))
            break;
    }
    prevEnv->At(position).pH = -log((x + y + z + aH)/1000)/log(10);
    prevEnv->At(position).Hydrox = x;
    prevEnv->At(position).Acetic = y;
    prevEnv->At(position).Buffer = z;

    return (Cac - y);
}

int nutrientDepletionAgar_OxygenAcetate(Array3D<LocalAga>** CurrentAgar, Array3D<LocalAga>** PreviousAgar, Array2D<LocalAga>** CurrentWall, Array2D<LocalAga>** PreviousWall, Array3D<double>& insideColonyDen)
{
    int count = 0;
    IntCoord positionAgar;
    int zi, yi, xi;
    int level = maxLevelsMG-1;
    double h2 = BoxLength*BoxLength*ipow(4,level);
    double dt = min(DiffAgarAce, DiffAgarCarbon)/10*h2;
    int n = 3600*UpdateTime*740/dt;

    restrictionAgar(CurrentAgar, PreviousAgar);

    while (count<=n)
    {
        yi = 0;
#pragma omp parallel for default(shared) private(zi,xi,positionAgar)
        for (zi=0;zi<BoxZAgar;zi++)
        {
            for (xi=0; xi<BoxX; xi++)
            {
                positionAgar.x = xi;
                positionAgar.y = yi;
                positionAgar.z = zi;

                UpdateDepletionAgar(CurrentAgar, PreviousAgar, CurrentWall, positionAgar, insideColonyDen, level, dt);
            }
        }
        
        
#pragma omp parallel for default(shared) private(zi,xi,positionAgar)
        for (zi=0;zi<BoxZAgar;zi++)
        {
            for (xi=0; xi<BoxX; xi++)
            {
                positionAgar.x = xi;
                positionAgar.y = yi;
                positionAgar.z = zi;
                PreviousAgar[level]->At(positionAgar).CarbonAgar = CurrentAgar[level]->Get(positionAgar.x,positionAgar.y,positionAgar.z).CarbonAgar;
                PreviousAgar[level]->At(positionAgar).AcetateAgar = CurrentAgar[level]->Get(positionAgar.x,positionAgar.y,positionAgar.z).AcetateAgar;
            }
        }
        
        count++;

    }
    
    interpolation(CurrentAgar, PreviousAgar, CurrentWall, PreviousWall);
    
    return count;
}

void restrictionAgar(Array3D<LocalAga>** Aga, Array3D<LocalAga>** prevAga)
{
    int zi, yi=0, xi;
    IntCoord pos1, pos2;
#pragma omp parallel for default(shared) private(zi,xi,pos1,pos2)
    for (zi=0; zi<BoxZAgar;zi++)
    {
        for (xi=0; xi<BoxX; xi++)
        {
            pos1.x = xi;
            pos1.y = yi;
            pos1.z = zi;
            int l = max(int(ceil(log((1+zi)*ipow(2,maxLevelsMG-1)*1.0/BoxZAgar)/log(2))),0);
            int lx = 0;
            if (xi<BoxX/2)
                lx = ceil(log(1-xi*2.0/BoxX)/log(2)+maxLevelsMG-1);
            else if (xi>BoxX/2)
                lx = floor(log(xi*2.0/BoxX-1)/log(2)+maxLevelsMG-1)+1;
            if (lx>l)
                l = lx;
            pos2.x = (xi-BoxX/2)*ipow(2,maxLevelsMG-1-l)+BoxX/2;
            pos2.y = 0;
            pos2.z = (zi+1)*ipow(2,maxLevelsMG-1-l)-1;
            if (l<maxLevelsMG-1)
            {
                Aga[maxLevelsMG-1]->At(pos1).CarbonAgar = Aga[l]->Get(pos2).CarbonAgar;
                Aga[maxLevelsMG-1]->At(pos1).AcetateAgar = Aga[l]->Get(pos2).AcetateAgar;
                prevAga[maxLevelsMG-1]->At(pos1).CarbonAgar = prevAga[l]->Get(pos2).CarbonAgar;
                prevAga[maxLevelsMG-1]->At(pos1).AcetateAgar = prevAga[l]->Get(pos2).AcetateAgar;
            }
        }
    }
}

void interpolation(Array3D<LocalAga>** Aga, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wall, Array2D<LocalAga>** prevWall)
{
    int zi, yi=0, xi;
    IntCoord pos1, pos2;
    double a1, a2, c1, c2;
    for (int l = maxLevelsMG-1; l>0; l--)
    {
#pragma omp parallel for default(shared) private(zi,xi,pos1,pos2, a1, a2, c1, c2)
        for (zi=0; zi<BoxZAgar;zi++)
        {
            for (xi=0; xi<BoxX; xi++)
            {
                a1 = 0;
                a2 = 0;
                c1 = 0;
                c2 = 0;
                pos1.x = xi;
                pos1.y = 0;
                pos1.z = zi;
                if (zi==0)
                {
                    for (int i1 = 0; i1<2; i1++)
                    {
                        pos2.x = (xi+i1)/2+BoxX/4;
                        pos2.y = 0;
                        pos2.z = 0;
                        c1 += Aga[l]->Get(pos2).CarbonAgar/4;
                        c2 += prevAga[l]->Get(pos2).CarbonAgar/4;
                        a1 += Aga[l]->Get(pos2).AcetateAgar/4;
                        a2 += prevAga[l]->Get(pos2).AcetateAgar/4;
                    }
                    int lx = 0;
                    if (xi<BoxX/2)
                        lx = ceil(log(1-xi*2.0/BoxX)/log(2)+l);
                    else if (xi>BoxX/2)
                        lx = floor(log(xi*2.0/BoxX-1)/log(2)+l)+1;
                    if (lx<0)
                        lx = 0;
                    pos2.x = (xi-BoxX/2)*ipow(2,l-lx)+BoxX/2;
                    c1 += Wall[lx]->Get(pos2.x,0).CarbonAgar/2;
                    c2 += prevWall[lx]->Get(pos2.x,0).CarbonAgar/2;
                    a1 += Wall[lx]->Get(pos2.x,0).AcetateAgar/2;
                    a2 += prevWall[lx]->Get(pos2.x,0).AcetateAgar/2;
                }
                else
                {
                    for (int i1 = 0; i1<2; i1++)
                    {
                        for (int i2 = 0; i2<2; i2++)
                        {
                            pos2.x = (xi+i1)/2+BoxX/4;
                            pos2.y = 0;
                            pos2.z = (zi-i2)/2;
                            c1 += Aga[l]->Get(pos2).CarbonAgar/4;
                            c2 += prevAga[l]->Get(pos2).CarbonAgar/4;
                            a1 += Aga[l]->Get(pos2).AcetateAgar/4;
                            a2 += prevAga[l]->Get(pos2).AcetateAgar/4;
                        }
                    }
                }
                Aga[l-1]->At(pos1).CarbonAgar = c1;
                prevAga[l-1]->At(pos1).CarbonAgar = c2;
                Aga[l-1]->At(pos1).AcetateAgar = a1;
                prevAga[l-1]->At(pos1).AcetateAgar = a2;
            }
        }
    }
}

void UpdateDepletionAgar(Array3D<LocalAga>** Aga, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wal, IntCoord position, Array3D<double>& insideColonyDen, int level, double dt)
{
    double Cxx, Axx, Czz, Azz;
    if (position.x==0)
    {
        Cxx = (2*prevAga[level]->Get(position.x+1,position.y,position.z).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Axx = (2*prevAga[level]->Get(position.x+1,position.y,position.z).AcetateAgar-2*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
    }
    else if (position.x==BoxX-1)
    {
        Cxx = (prevAga[level]->Get(position.x-1,position.y,position.z).CarbonAgar-prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Axx = (prevAga[level]->Get(position.x-1,position.y,position.z).AcetateAgar-prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
    }
    else
    {
        Cxx = (prevAga[level]->Get(position.x+1,position.y,position.z).CarbonAgar+prevAga[level]->Get(position.x-1,position.y,position.z).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Axx = (prevAga[level]->Get(position.x+1,position.y,position.z).AcetateAgar+prevAga[level]->Get(position.x-1,position.y,position.z).AcetateAgar-2*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
    }
    if (position.z==0)
    {
        int lx = 0;
        if (position.x<BoxX/2)
            lx = ceil(log(1-position.x*2.0/BoxX)/log(2)+maxLevelsMG-1);
        else if (position.x>BoxX/2)
            lx = floor(log(position.x*2.0/BoxX-1)/log(2)+maxLevelsMG-1)+1;
        if (lx<0)
            lx=0;
        int x = (position.x-BoxX/2)*ipow(2,maxLevelsMG-1-lx)+BoxX/2;
        IntCoord positionWall3D;
        positionWall3D.x = x;
        positionWall3D.y = 0;
        positionWall3D.z = 0;
        double Qc = 0, Qa = 0;
        if (lx==0 && insideColonyDen.Get(positionWall3D)>0.5)
        {
//            Qc = (Wal[0]->Get(x, 0).CarbonAgar - prevAga[0]->Get(x, 0, 0).CarbonAgar);
//            Qa = (Wal[0]->Get(x, 0).AcetateAgar - prevAga[0]->Get(x, 0, 0).AcetateAgar);
            Qc = (Wal[0]->Get(x, 0).CarbonAgar - prevAga[level]->Get(position.x, 0, 0).CarbonAgar)/pow(2.0,level);
            Qa = (Wal[0]->Get(x, 0).AcetateAgar - prevAga[level]->Get(position.x, 0, 0).AcetateAgar)/pow(2.0,level);
            Czz = 0;
            Azz = 0;
        }
        else
        {
            Czz = (prevAga[level]->Get(position.x,position.y,position.z+1).CarbonAgar+Qc*pow(2.0,level)-prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
            Azz = (prevAga[level]->Get(position.x,position.y,position.z+1).AcetateAgar+Qa*pow(2.0,level)-prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
        }
    }
    else if (position.z==BoxZAgar-1)
    {
        Czz = (2*prevAga[level]->Get(position.x,position.y,position.z-1).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Azz = (2*prevAga[level]->Get(position.x,position.y,position.z-1).AcetateAgar-2*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
    }
    else
    {
        Czz = (prevAga[level]->Get(position.x,position.y,position.z+1).CarbonAgar+prevAga[level]->Get(position.x,position.y,position.z-1).CarbonAgar-2*prevAga[level]->Get(position.x,position.y,position.z).CarbonAgar)/(BoxLength*BoxLength*pow(4.0,level));
        Azz = (prevAga[level]->Get(position.x,position.y,position.z+1).AcetateAgar+prevAga[level]->Get(position.x,position.y,position.z-1).AcetateAgar-2*prevAga[level]->Get(position.x,position.y,position.z).AcetateAgar)/(BoxLength*BoxLength*pow(4.0,level));
    }
    
    Aga[level]->At(position).CarbonAgar += dt*DiffAgarCarbon*(Cxx+Czz);
//    Aga[level]->At(position).CarbonAgar = min(max(Aga[level]->At(position).CarbonAgar,0.0),maxCarbon);
    Aga[level]->At(position).AcetateAgar += dt*DiffAgarAce*(Axx+Azz);
    
}

