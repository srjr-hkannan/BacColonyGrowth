#ifndef CONSTANTS_H_
#define CONSTANTS_H_

// default values for all variables

// physical variables
extern double cellRadius;
extern double L_divide;		
extern double k_cc;			
extern double k_wc;		
extern double varL;
extern double varAngle;
extern double var_pos;
extern double viscosity;		
extern double wall_rough;
extern double gamma_t;
extern double cell_mu;
extern double wall_mu;
extern double density_threshold;
extern double tension;
extern double DH;

// time
extern double t_max;
extern double initial_dt;
extern double OutputTime;
extern double t0;
extern double UpdateTime;
extern double maxO2_Agar;
// non-physical constants
extern int BoxX;
extern int BoxY;
extern int BoxZ;
extern int BoxZAgar;
extern int maxLevelsMG;

extern double BoxLength;
extern int FilterLen;

// directory name
extern char DirName[500];

// Restart information
extern int restart;
extern int restartIndex;
extern char restartDir[500];

// number max cells to simulate (for allocating memory)
extern int maxCells;	// maximum number of cells in the simulation

// nutrient constants
extern double Tortuosity;
extern double KC_aer;
extern double KC_ana;

extern double KW;
extern double KO;
extern double KA;
extern double KI;
extern double Carbon_to_Waste;
extern double C_rate;
extern double W_rate;
extern double qCarAer;
extern double qCarAna;
extern double qAceAer;
extern double qO2Car;
extern double qO2Ace;
extern double pAceAer;
extern double pAceAna;
extern double DiffColony;
extern double DiffAgar;
extern double DiffColonyCarbon;
extern double DiffAgarCarbon;
extern double DiffColonyO2;
extern double DiffAgarO2;
extern double DiffColonyAce;
extern double DiffAgarAce;
extern double maxCarbon;
extern double maxO2;
extern double maxGrowthRate1;
extern double maxGrowthRate2;
extern double maxGrowthRateCarAer;
extern double maxGrowthRateCarAna;
extern double maxGrowthRateAceAer;
extern double Maintenance_rate;
extern double Cdt;
extern double ConvCrit;
extern int minIter;
extern int maxIter;
extern int InterfaceCondition;
extern bool NutrientGSI;
extern double Rc;

// maintenance constants
extern double qCarAernt;
extern double qCarAnant;
extern double qAceAernt;
extern double qO2Carnt;
extern double qO2Acent;
extern double pAceAnant;
extern double qEthaAnant;
extern double qForAnant;

// pH calculation constants
extern double Kw;
extern double Ka;
extern double Kb;
extern double Cb;
extern double c1;

// colony constants
extern bool crossFeeding;
extern int refinementGridHeight;


#endif /* CONSTANTS_H_ */
