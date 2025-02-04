#ifndef NUTRIENTS_H_
#define NUTRIENTS_H_

#include "tools.h"
#include "Constants.h"
//#include "Array.h"

template<typename T>
class Array3D;

struct LocalEnv
{
	LocalEnv()
	{
		Carbon = 0.0;
		GrowthRate1 = 0.0;
		GrowthRate2 = 0.0;
		Waste = 0.0;
        	Oxygen = maxO2;
       		Acetate = 0.0;

		pH = 7.0;
		Buffer = 0.0;
		Acetic = 0.0;
		Hydrox = 0.0;
        
		// Acetate consumption/production
		qCarmnt = 0.0;
		qAcemnt = 0.0;
		pAcemnt = 0.0;
		qOxymnt = 0.0;

        // Component of the growth rate
        GrCarAer = 0.0;
        GrCarAna = 0.0;
        GrAceAer = 0.0;
        
    }

	LocalEnv(double _carbon, double _waste, double _growthrate1, double _growthrate2, double _oxygen, double _acetate): Carbon(_carbon), Waste(_waste), GrowthRate1(_growthrate1), GrowthRate2(_growthrate2), Oxygen(_oxygen), Acetate(_acetate) {}

	double Carbon;
	double Waste;
	double GrowthRate1;
	double GrowthRate2;
    double Oxygen;
    double Acetate;
	
	double pH;
	double Buffer;
	double Acetic;
	double Hydrox;

	// Acetate consumption/production
	double qCarmnt;
	double qAcemnt;
	double pAcemnt;
	double qOxymnt;
    
    // Component of the growth rate
    double GrCarAer;
    double GrCarAna;
    double GrAceAer;
};

struct LocalAga
{
	LocalAga()
	{
		CarbonAgar = maxCarbon;
		WasteAgar = 0.0;
        	OxygenAgar = maxO2;
        	AcetateAgar = 0.0;
	}

	LocalAga(double _carbonAgar, double _wasteAgar, double _growthrateAgar1, double _growthrateAgar2, double _oxygenAgar, double _acetateAgar): CarbonAgar(_carbonAgar), WasteAgar(_wasteAgar), OxygenAgar(_oxygenAgar), AcetateAgar(_acetateAgar) {}

	double CarbonAgar;
	double WasteAgar;
    double OxygenAgar;
    double AcetateAgar;
};

template<typename T>
class Array2D;





int UpdateEnvArray(Array3D<LocalEnv>* CurrentColony, Array3D<LocalEnv>* PreviousColony, Array3D<LocalAga>** CurrentAgar, Array3D<LocalAga>** PreviousAgar, Array2D<LocalAga>** CurrentWall, Array2D<LocalAga>** PreviousWall, Array3D<double>& DensityShiftP, Array3D<double>& Density1ShiftP, Array3D<double>& Density2ShiftP, Array2D<double>& WallDensityShiftP,Array2D<double>& WallDensity1ShiftP,Array2D<double>& WallDensity2ShiftP, int minX, int maxX, int minY, int maxY, int maxH, Array2D<double>& Height, Array3D<double>& insideColonyDen);
void UpdateEnvironment(Array3D<LocalEnv>* Env, Array3D<LocalEnv>* prevEnv, Array2D<LocalAga>* prevWal, Array3D<double>& DensityShiftP, Array3D<double>& Density1ShiftP, Array3D<double>& Density2ShiftP, Array2D<double>& WallDensityShiftP,IntCoord position, Array3D<double>& insideColonyDen);
void UpdateAgar(Array3D<LocalAga>* Aga, Array3D<LocalAga>* prevAga, Array2D<LocalAga>* Wal, IntCoord position);
void UpdateAgarMultigrid(Array3D<LocalAga>** Aga, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wal, IntCoord position, int level);
void UpdateWall(Array3D<LocalEnv>* Env, Array3D<LocalAga>* prevAga, Array2D<LocalAga>* Wal, Array2D<LocalAga>* prevWal, Array2D<double>& WallDensityShiftP, Array2D<double>& WallDensity1ShiftP, Array2D<double>& WallDensity2ShiftP,IntCoord2D position, Array2D<double>& Height, Array3D<double>& insideColonyDen);
void UpdateWallMultigrid(Array3D<LocalEnv>* Env, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wal, Array2D<LocalAga>** prevWal, Array2D<double>& WallDensityShiftP, Array2D<double>& WallDensity1ShiftP, Array2D<double>& WallDensity2ShiftP, IntCoord2D position, int level, Array2D<double>& Height, Array3D<double>& insideColonyDen);
void ShiftGrowthrate(Array3D<LocalEnv>* CurrentColony, Array2D<LocalAga>* CurrentWall, int BoxX, int BoxY, int BoxZ);

int UpdateEnvArray_OxygenAcetate(Array3D<LocalEnv>* CurrentColony, Array3D<LocalEnv>* PreviousColony, Array3D<LocalAga>** CurrentAgar, Array3D<LocalAga>** PreviousAgar, Array2D<LocalAga>** CurrentWall, Array2D<LocalAga>** PreviousWall, Array3D<double>& DensityShiftP, Array3D<double>& Density1ShiftP, Array3D<double>& Density2ShiftP, Array2D<double>& WallDensityShiftP,Array2D<double>& WallDensity1ShiftP,Array2D<double>& WallDensity2ShiftP, int minX, int maxX, int minY, int maxY, int maxH, Array2D<double>& Height, Array3D<double>& insideColonyDen);
void UpdateEnvironment_OxygenAcetate(Array3D<LocalEnv>* Env, Array3D<LocalEnv>* prevEnv, Array2D<LocalAga>* prevWal, Array3D<double>& DensityShiftP, Array3D<double>& Density1ShiftP, Array3D<double>& Density2ShiftP, Array2D<double>& WallDensityShiftP,IntCoord position, Array3D<double>& insideColonyDen);
void UpdateAgarMultigrid_OxygenAcetate(Array3D<LocalAga>** Aga, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wal, IntCoord position, int level);
void UpdateWallMultigrid_OxygenAcetate(Array3D<LocalEnv>* Env, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wal, Array2D<LocalAga>** prevWal, Array2D<double>& WallDensityShiftP, Array2D<double>& WallDensity1ShiftP, Array2D<double>& WallDensity2ShiftP, IntCoord2D position, int level, Array2D<double>& Height, Array3D<double>& insideColonyDen);
double acetateBufferSolve(double Cac, Array3D<LocalEnv>* prevEnv, IntCoord position);
double acetateBufferSolve(double Cac);
int nutrientDepletionAgar_OxygenAcetate(Array3D<LocalAga>** CurrentAgar, Array3D<LocalAga>** PreviousAgar, Array2D<LocalAga>** CurrentWall, Array2D<LocalAga>** PreviousWall, Array3D<double>& insideColonyDen);
void restrictionAgar(Array3D<LocalAga>** Aga, Array3D<LocalAga>** prevAga);
void interpolation(Array3D<LocalAga>** Aga, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wall, Array2D<LocalAga>** prevWall);
void UpdateDepletionAgar(Array3D<LocalAga>** Aga, Array3D<LocalAga>** prevAga, Array2D<LocalAga>** Wal, IntCoord position, Array3D<double>& insideColonyDen, int level, double dt);
#endif
