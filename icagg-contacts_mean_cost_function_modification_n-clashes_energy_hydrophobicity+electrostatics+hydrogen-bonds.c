#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <time.h>
//=============Parameters you might want to change===============


#define MCSTEPS			(2000)
#define PROTSTEP        	(100)

#define STRUCTS	                (4005)
#define PAIRS			(1000)
#define CONTACT_PAIR    	(2)
#define SCALE			(1.)		//d(CM1-CM2)=SCALE*(Rg1+Rg2)
#define SCALE_1                 (1.2)           //scale for adjusting CM-distance in the initial pose/orientation sampling
#define SCALE_2                 (0.9)		//scale for adjusting CM-distance in the initial pose/orientation sampling	

#define SCALE_CLASHES   	(0.5)

#define HARDCORE_R		(1.)
#define CONTACT_R		(2.33)

#define Tm			(0.15)
 
//Tm = temperature at which structures were sampled
//ensemble sizes: NAT(WT)pH5.2=4005
//ensemble sizes: NAT(WT)pH6.2=4005
//ensemble sizes: NAT(WT)pH7.2=4005
//ensemble sizes: NAT(D76N)pH5.2=4005
//ensemble sizes: NAT(D76N)pH6.2=4005
//ensemble sizes: NAT(D76N)pH7.2=4005
//ensemble sizes: I1(D76N)pH5.2=4005
//ensemble sizes: I1(D76N)pH6.2=4005
//ensemble sizes: I1(D76N)pH7.2=4005
//ensemble sizes: I2(D76N)pH5.2=4005
//ensemble sizes: I2(D76N)pH6.2=4005
//ensemble sizes: I2(D76N)pH7.2=4005

//end to parameters you might want to change=====================

#define MAXLENLINE 	(256)
#define LENCHARGESFILE  (287)
#define MAXRESIDUE      (99)
#define MAXATOMNUM	(1100)
#define MAXCONTACT  	(50000)
#define HIS_TYPE_HIS	(0)
#define HIS_TYPE_HSD	(1)
#define HIS_TYPE_HSC	(2)
#define NUM_AA_TYPE	(26)
#define NUM_Atom_TYPE	(32)
#define PI		(3.14159265358979323846)


//==============Directories you might want to change========================================================================================
char workingdir[]="/home/rjloureiro/DOCKING/B2M-D76N-I1-PH5p2-newcf/";
char writingdir[]="/home/rjloureiro/DOCKING/DOCKING_B2M_D76N-I1-PH5p2_newcf_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new/";
char rootdir[]="/home/rjloureiro";
char ss[]="DOCKING_B2M_D76N-I1-PH5p2_newcf_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new";
//==========================================================================================================================================


int ResNum[2];
int nAtom[2];

typedef struct	{

	int Res_Type;
	char AtmName[5];
	char Atm_Type[5];
	float Atm_Charge;
	
}MYCHARGES,*PMYCHARGES;

typedef struct	{

	double x,y,z;
	double r0;
	
	int NodeIndex;
	int Res_Type;
	int AtomIndex;
	int Atom_Type;
	int ASP_Type;
	char AtmName[5];
	char AtmHBlabel[5];
	
}MYATOM,*PMYATOM;

typedef struct	{

	double Rg;
	double RMSD;
	int contEn;
	
}MYCONF, *PMYCONF;

typedef struct	{

	double r01;
	double r02;
	double rsq;
	
}MYDISTANCES, *PMYDISTANCES;

double radii[]={1.61, 1.76, 1.88, 1.88, 1.88, 1.64, 1.64, 1.64, 1.64, 1.42, 1.46, 1.77, 1.77, 1.05, 0.58};

char AA_Name[NUM_AA_TYPE][5]={"METN","ILE","ILEN","VAL","LEU","PHE","CYS","MET","METC","ALA","GLY","THR","SER","TRP","TYR","PRO","HIS","HISH","GLN",
"ASN","GLU","GLUH","ASP","ASPH","LYS","ARG"};/*check*/

float atomic_solvation_parameters[6]={18.0, -7.0, 18.0, -20.0, -34.0, 0.0};

char atom_name[NUM_Atom_TYPE][4]={"CB","CG","CG1","CG2","CD","CD1","CD2","CE","CE1","CE2","CE3","CH2","CZ","CZ2","CZ3","ND1","ND2","NE","NE1","NE2",
"NH1","NH2","NZ","OD1","OD2","OE1","OE2","OG","OG1","OH","SG","SD"};

float atomic_SASA[NUM_AA_TYPE][NUM_Atom_TYPE]={
{24.0,  29.6, 0.0,  0.0,  0.0,  0.0,  0.0,  74.8, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  36.4},
{11.0,  0.0,  37.0, 58.4, 0.0,  43.7, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{11.0,  0.0,  37.0, 58.4, 0.0,  43.7, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{11.1,  0.0,  57.4, 59.9, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{23.1,  9.6, 0.0,  0.0,  0.0,  63.4, 61.7, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{29.3,  0.6,  0.0,  0.0,  0.0,  22.1, 21.9, 0.0,  36.7, 36.1, 0.0,  0.0,  37.7, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{38.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  65.0, 0.0},
{24.0,  29.6, 0.0,  0.0,  0.0,  0.0,  0.0,  74.8, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  36.4},
{24.0,  29.6, 0.0,  0.0,  0.0,  0.0,  0.0,  74.8, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  36.4},
{71.9, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{18.0,  0.0,  0.0,  63.1, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  33.5, 0.0,  0.0,  0.0},
{44.3, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  41.5, 0.0,  0.0,  0.0,  0.0},
{30.1,  1.3,  0.0,  0.0,  0.0,  30.2, 2.2,  0.0,  0.0,  3.3,  21.5, 37.9, 0.0,  38.5, 34.1, 0.0,  0.0,  0.0,  29.8, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{29.2,  0.4,  0.0,  0.0,  0.0,  22.6, 21.5, 0.0,  34.8, 34.1, 0.0,  0.0,  2.6, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  52.9, 0.0,  0.0},
{39.1,   44.3,  0.0,  0.0,  27.6,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{34.0,  1.6,  0.0,  0.0,  0.0,  0.0,  33.0, 0.0,  51.6, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  11.4, 0.0,  0.0,  0.0,  30.5, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{34.0,  1.6,  0.0,  0.0,  0.0,  0.0,  33.0, 0.0,  51.6, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  11.4, 0.0,  0.0,  0.0,  30.5, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{26.6,  30.9, 0.0,  0.0,  3.7, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  60.0, 0.0,  0.0,  0.0,  0.0,  0.0,  34.2, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{32.1,  3.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  59.1, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  31.1, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{26.8,  33.2, 0.0,  0.0,  5.5, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  41.4, 41.5, 0.0,  0.0,  0.0,  0.0,  0.0},
{26.8,  33.2, 0.0,  0.0,  5.5, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  41.4, 41.5, 0.0,  0.0,  0.0,  0.0,  0.0},
{33.5,  5.1, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  38.7, 40.9, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{33.5,  5.1, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  38.7, 40.9, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{25.9,  22.6, 0.0,  0.0,  28.0, 0.0,  0.0,  36.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  74.6, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{25.9,  23.9, 0.0,  0.0,  29.6, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  2.2, 0.0,  0.0,  0.0,  0.0,  13.3, 0.0,  0.0,  58.6, 63.4, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
};/*atomic solvent accessible surface area*/

/*float hydropathic_inter_residue_potentials[NUM_AA_TYPE][NUM_AA_TYPE]={
{-1.0000, -1.0000, -0.9700, -0.9200, -0.8100, -0.7800, -0.7100, -0.7100, -0.7000, +0.0000, +0.0000, -0.0000, +0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000},
{-1.0000, -1.0000, -0.9700, -0.9200, -0.8100, -0.7800, -0.7100, -0.7100, -0.7000, +0.0000, +0.0000, -0.0000, +0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000},
{-0.9700, -0.9700, -0.9300, -0.8900, -0.7800, -0.7500, -0.6800, -0.6800, -0.6700, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000},
{-0.9200, -0.9200, -0.8900, -0.8400, -0.7300, -0.7000, -0.6300, -0.6300, -0.6200, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000},
{-0.8100, -0.8100, -0.7800, -0.7300, -0.6200, -0.5900, -0.5200, -0.5200, -0.5100, +0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000},
{-0.7800, -0.7800, -0.7500, -0.7000, -0.5900, -0.5600, -0.4900, -0.4900, -0.4800, +0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000},
{-0.7100, -0.7100, -0.6800, -0.6300, -0.5200, -0.4900, -0.4200, -0.4200, -0.4100, +0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000},
{-0.7100, -0.7100, -0.6800, -0.6300, -0.5200, -0.4900, -0.4200, -0.4200, -0.4100, +0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000},
{-0.7000, -0.7000, -0.6700, -0.6200, -0.5100, -0.4800, -0.4100, -0.4100, -0.4000, -0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, -0.0000, -0.0000, +0.0000, -0.0000},
{+0.0000, +0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, -0.0000, -0.0000, +0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000},
{+0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, +0.0000, -0.0000, -0.0000, -0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000},
{-0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, +0.0000, -0.0000, +0.0000, +0.0000, -0.0000, -0.0000, +0.0000, -0.0000},
{+0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, +0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000},
{-0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, +0.0000, +0.0000, -0.0000, +0.0000},
{+0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, +0.0000, -0.0000, -0.0000, -0.0000, +0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000},
{+0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, +0.0000, -0.0000, -0.0000, -0.0000, -0.0000, +0.7100, +0.7100, +0.7500, +0.7500, +0.7500, +0.7500, +0.7500, +0.7500, +0.7900, +0.8600},
{+0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, +0.0000, -0.0000, -0.0000, -0.0000, -0.0000, +0.7100, +0.7100, +0.7500, +0.7500, +0.7500, +0.7500, +0.7500, +0.7500, +0.7900, +0.8600},
{+0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, +0.7500, +0.7500, +0.7800, +0.7800, +0.7800, +0.7800, +0.7800, +0.7800, +0.8300, +0.8900},
{+0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, -0.0000, -0.0000, +0.7500, +0.7500, +0.7800, +0.7800, +0.7800, +0.7800, +0.7800, +0.7800, +0.8300, +0.8900},
{-0.0000, -0.0000, +0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, +0.0000, +0.0000, -0.0000, -0.0000, -0.0000, +0.7500, +0.7500, +0.7800, +0.7800, -0.0000, -0.0000, -0.0000, -0.0000, +0.0000, -0.0000},
{-0.0000, -0.0000, +0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, +0.0000, +0.0000, -0.0000, -0.0000, -0.0000, +0.7500, +0.7500, +0.7800, +0.7800, -0.0000, -0.0000, -0.0000, -0.0000, +0.0000, -0.0000},
{+0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, +0.0000, -0.0000, -0.0000, +0.0000, -0.0000, +0.7500, +0.7500, +0.7800, +0.7800, -0.0000, -0.0000, +0.0000, +0.0000, +0.0000, -0.0000},
{+0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, +0.0000, -0.0000, -0.0000, +0.0000, -0.0000, +0.7500, +0.7500, +0.7800, +0.7800, -0.0000, -0.0000, +0.0000, +0.0000, +0.0000, -0.0000},
{+0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, +0.0000, +0.0000, -0.0000, -0.0000, -0.0000, +0.7900, +0.7900, +0.8300, +0.8300, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000},
{+0.0000, +0.0000, -0.0000, +0.0000, +0.0000, +0.0000, +0.0000, +0.0000, -0.0000, -0.0000, +0.0000, -0.0000, -0.0000, +0.0000, -0.0000, +0.8600, +0.8600, +0.8900, +0.8900, -0.0000, -0.0000, -0.0000, -0.0000, +0.0000, +0.0000},
};/*hydropathic inter-residue contact potentials between the 20 amino acids*/

/*double SASA_AA[MAXRESIDUE];*/

MYCHARGES Charges[LENCHARGESFILE];

MYATOM Atoms[MAXATOMNUM][2];

MYATOM Moved[MAXATOMNUM];

MYATOM RotMoved[MAXATOMNUM];

MYATOM stack1[MAXATOMNUM];

MYATOM stack2[MAXATOMNUM];

MYCONF Conf[2];
MYDISTANCES dist;

double MEAN_R[3][2];

float q0[MAXCONTACT];
float q1[MAXCONTACT];
/*double sasa_aa[MAXCONTACT][CONTACT_PAIR];*/
int it[MAXCONTACT];
double rd[MAXCONTACT];
int ct0[MAXCONTACT];
int ct1[MAXCONTACT];
int at0[MAXCONTACT];
int at1[MAXCONTACT];
int ri0[MAXCONTACT];
int ri1[MAXCONTACT];
double ro0[MAXCONTACT];
double ro1[MAXCONTACT];
int asp0[MAXCONTACT];
int asp1[MAXCONTACT];
double sum_vdW[MAXCONTACT];
double hidropathy_0[MAXCONTACT];
double hidropathy_1[MAXCONTACT];
double hp[MAXCONTACT];
double ep[MAXCONTACT];
double hbp[MAXCONTACT];
double ip[MAXCONTACT];


//CODE PARTS
int main(void);
void READER(int a, int b);
void CENTRALIZER(double a);
void INITIALORIENTATOR(void);
int TYPECONTACTS(int ct0[MAXCONTACT], int ct1[MAXCONTACT], int at0[MAXCONTACT], int at1[MAXCONTACT], int ri0[MAXCONTACT], int ri1[MAXCONTACT], double ro0[MAXCONTACT], double ro1[MAXCONTACT], double rd[MAXCONTACT], int it[MAXCONTACT], float q0[MAXCONTACT], float q1[MAXCONTACT], int asp0[MAXCONTACT], int asp1[MAXCONTACT]);
int CONTACTS(void);
int CLASHES(void);
int TYPECONTACTSMOVED(int ct0[MAXCONTACT], int ct1[MAXCONTACT], int at0[MAXCONTACT], int at1[MAXCONTACT], int ri0[MAXCONTACT], int ri1[MAXCONTACT], double ro0[MAXCONTACT], double ro1[MAXCONTACT], double rd[MAXCONTACT], int it[MAXCONTACT], float q0[MAXCONTACT], float q1[MAXCONTACT], int asp0[MAXCONTACT], int asp1[MAXCONTACT]);
int CONTACTSMOVED(void);
int CLASHESMOVED(void);
void MCMOVE(double a);
void READCHARGES(void);
/*void READPOTENTIALS(void);
/*double RESIDUECHARGES(int residues_index, int contacts);*/
/*void READSASA(void);*/


//PDB HANDLING
double abbs(double a);
int IsANumber(char c);
int Get_AA_Type(char SzAAName[]);
int Get_Atom_Type(char SzAtomName[]);
void AssignAtomRadii(void);
int QueryAtomType(char *s, char *res);
void ExportSnapshot(int Index1, int Index2);


//RANDOM NUMBER GENERATOR HEADER
void RandomInitialise(int ij,int kl);
double RandomUniform(void);
double RandomGaussian(double mean,double stddev);
int RandomInt(int lower,int upper);
double RandomDouble(double lower,double upper);


double restdist=0;
double currentdist=0;
double olddist=0;
double mcinternal1=0.025;			//perturbation amplitude for translations
double mcinternal2=0.025;			//perturbation amplitude for rotations
double TOL=0.2;					//tolerance of CM-diffusion in percent of d(CM1-CM2)
double TOL0=0.2;
double x = 1.3;


int DIR;		

//	COMMENT on DIR
//	direction of docking: 
//	+1 ... 	positive x-axis
//	-1 ... 	negative x-axis
//	+2 ...	positive y-axis
//	-2 ...	negative y-axis
//	+3 ...	positive z-axis
//	-3 ...	negative z-axis

int ORIENTATION;

//      COMMENT on ORIENTATION
//      initial orientation of docking:
//      0  ...  no rotation
//	+1 ... 	+90º rotation arround x-axis
//	-1 ... 	-90º rotation arround x-axis
//      +2 ...  +90º rotation arround y-axis
//      -2 ...  -90º rotation arround y-axis
//      +3 ...  +90º rotation arround z-axis
//      -3 ...  -90º rotation arround z-axis

int main(void)
{

	FILE *output1;
	FILE *output2;
	/*FILE *output3;
	FILE *output4;
	FILE *output5;
	FILE *output6;
	FILE *output7;
	FILE *output8;*/
	
	int i,j,k,l,m,o,p,e,mc,helper1,helper2,index,ncontacts,N,contact_flag,restdistchange1,restdistchange2;
	int newrestdist1,newrestdist2;	
	int n = 2;
	int cc,clashinit[6][7][3],contactinit[6][7][3],hitwall;
	int hydropathic_flag;
	int residue_index_0;
	int residue_index_1;
	int accepted, rejected;	
	int initial_counter;
	int final_counter;
	int counter_before_mcmove;
	int counter_after_mcmove;
	double SASA_AA[MAXRESIDUE];
	char oo[256], pp[256], qq[256], rr[256], s[256], tt[256], uu[256], vv[256];
	/*double diffcontactsenergy[PAIRS][MCSTEPS],diffcontactsenergy_sum[PAIRS],diffcontactsenergy_mean[PAIRS];*/
	double shapeinit[6][7][3];
	double Etransl, Econtnew, Etranslold, Econtold,Eallold,Eallnew;
	double contactnew, contactold;
	double initclash,clashnew,clashold;
	double Pacc,Pdecl;
	double dummie;
	double lastTOL;
	int CONTACTS_SUM = 0;
	double CONTACTS_MEAN;
	double contacts_diff_sum = 0.0;
	double VARIANCE_CONTACTS;
	double VARIANCE_CONTACTS_MEAN;
	int CLASHES_SUM = 0;
	double CLASHES_MEAN;
	double ENERGY_SUM = 0.0;
	double ENERGY_MEAN;
	double energy_diff_sum = 0.0;
	double VARIANCE_ENERGY;
	double VARIANCE_ENERGY_MEAN;
	double initial_energy_sum = 0.0;
	double initial_hp_energy_sum = 0.0;
	double initial_ep_energy_sum = 0.0;
	double initial_hbp_energy_sum = 0.0;
	double energy_sum_before_mcmove = 0.0;
	double energy_sum_after_mcmove = 0.0;
	double final_energy_sum = 0.0;
	double final_hp_energy_sum = 0.0;
	double final_ep_energy_sum = 0.0;
	double final_hbp_energy_sum = 0.0;
	
	
	int CCONTACTS[PAIRS];
	int CCLASHES[PAIRS];
	double ENERGY[PAIRS];

	/*float * ptq0; 
	float * ptq1;
	int * ptct0;
	int * ptct1;
	int * ptri0;
	int * ptri1;
	double * ptro0;
	double * ptro1;
	double * ptrd;
	int * ptit;
	double * ptsvdw;
	double * ptep;
	double * pthp;
	double * pthbp;
	double * ptip;*/

		
	mkdir(ss,457);
	
	cc=time(NULL);
	RandomInitialise(cc%20000,cc%10000);
	//READPOTENTIALS();
	READCHARGES();
	
	accepted=0;
	rejected=0;
	
	chdir(writingdir);
	

	sprintf(oo,"B2M_D76N-I1-PH5p2_2000_T-0.15_newcf_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new.dat");
	output1=fopen(oo,"a");
	
	sprintf(pp,"B2M_D76N-I1-PH5p2_contacts_energy_mean_2000_T-0.15_newcf_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new.dat");
	output2=fopen(pp,"a");
	
	/*sprintf(qq,"B2M_D76N_I2_pH6p2_contacts_energy_mean_10000_T-0.7_n-clashes_energy_hidropathy+electrostatics_5_test_potentials.dat");
	output3=fopen(qq,"a");
	
	sprintf(qq,"B2M_D76N_I2_pH6p2_contacts_energy_mean_10000_T-0.7_n-clashes_energy_hidropathy+electrostatics_5_residues_distance-mean.dat");
	output3=fopen(qq,"a");
	
	sprintf(rr,"B2M_D76N_I2_pH6p2_contacts_energy_mean_10000_T-0.7_n-clashes_energy_hidropathy+electrostatics_5_sum_vdW-mean.dat");
	output4=fopen(rr,"a");
	
	sprintf(s,"B2M_D76N_I2_pH6p2_contacts_energy_mean_10000_T-0.7_n-clashes_energy_hidropathy+electrostatics_5_1.4_sum_vdW-mean.dat");
	output5=fopen(s,"a");

	sprintf(tt,"B2M_D76N_I2_pH6p2_contacts_energy_mean_10000_T-0.7_n-clashes_energy_hidropathy+electrostatics_5_2.3_sum_vdW-mean.dat");
	output6=fopen(tt,"a");
	
	sprintf(uu,"B2M_D76N_I2_pH6p2_contacts_energy_mean_10000_T-0.7_n-clashes_energy_hidropathy+electrostatics_5_clashes_cutoff-mean.dat");
	output7=fopen(uu,"a");

	sprintf(vv,"B2M_D76N_I2_pH6p2_contacts_energy_mean_10000_T-0.7_n-clashes_energy_hidropathy+electrostatics_5_contacts_cutoff-mean.dat");
	output8=fopen(vv,"a");*/

	/*sprintf(uu,"Total_energy_mc-step_10000_T-0.7_n-clashes_energy_hidropathy+electrostatics_5.dat");
	output3=fopen(uu,"a");
	
	/*sprintf(vv,"B2M_D76N_I2_pH6p2_diffcontactsenergy_mean_10000_T-0.7_n-clashes_energy_hidropathy+electrostatics_5.dat");
	output8=fopen(vv,"a");*/
	
	
	/*for (e=0; e < PAIRS; e++)
	{
		diffcontactsenergy_sum[e] = 0.0;
	}*/
	
	chdir(workingdir);
	
	/*READSASA();*/
		
	for (k=1;k<=PAIRS;k++)
	{
		
			chdir(workingdir);
			
			//PICKING TWO STRUCTURE INDICES
			i=RandomInt(1,STRUCTS);
			j=RandomInt(1,STRUCTS);
							
			READER(i,j);
					
			chdir(writingdir);
			
			index=0;
				
			//restdist=SCALE*(Conf[0].Rg + Conf[1].Rg);
		
			restdist=31.95;	//SCALE*(2*Rg(D76N-I1))
			
			hitwall=0;
			
			//BEGINNING of INITIAL DIRECTION/ORIENTATION sampling
		
			for (DIR=-3;DIR<=-1;DIR++)
			{
				for(ORIENTATION=-3;ORIENTATION<=-1;ORIENTATION++)
				{
					CENTRALIZER(restdist);
					clashinit[DIR+3][ORIENTATION+3][0] = CLASHES();
					contactinit[DIR+3][ORIENTATION+3][0] = -CONTACTS();
					clashinit[DIR+3][ORIENTATION+3][1]=DIR;
					contactinit[DIR+3][ORIENTATION+3][1]=DIR;
					clashinit[DIR+3][ORIENTATION+3][2]=ORIENTATION;
					contactinit[DIR+3][ORIENTATION+3][2]=ORIENTATION;
					chdir(workingdir);
					READER(i,j);
					chdir(writingdir);										
				}
				for(ORIENTATION=1;ORIENTATION<=4;ORIENTATION++) 
				{
					CENTRALIZER(restdist);
					clashinit[DIR+3][ORIENTATION+2][0] = CLASHES();
					contactinit[DIR+3][ORIENTATION+2][0] = -CONTACTS();
					clashinit[DIR+3][ORIENTATION+2][1]=DIR;
					contactinit[DIR+3][ORIENTATION+2][1]=DIR;
					clashinit[DIR+3][ORIENTATION+2][2]=ORIENTATION;
					contactinit[DIR+3][ORIENTATION+2][2]=ORIENTATION;
					chdir(workingdir);
					READER(i,j);
					chdir(writingdir);															
				}
			}
			for (DIR=1;DIR<=3;DIR++)
			{
				for(ORIENTATION=-3;ORIENTATION<=-1;ORIENTATION++)
				{
					CENTRALIZER(restdist);
					clashinit[DIR+2][ORIENTATION+3][0] = CLASHES();
					contactinit[DIR+2][ORIENTATION+3][0] = -CONTACTS();
					clashinit[DIR+2][ORIENTATION+3][1]=DIR;
					contactinit[DIR+2][ORIENTATION+3][1]=DIR;
					clashinit[DIR+2][ORIENTATION+3][2]=ORIENTATION;
					contactinit[DIR+2][ORIENTATION+3][2]=ORIENTATION;
					chdir(workingdir);
					READER(i,j);
					chdir(writingdir);																				
				}
				for(ORIENTATION=1;ORIENTATION<=4;ORIENTATION++) 
				{
					CENTRALIZER(restdist);
					clashinit[DIR+2][ORIENTATION+2][0] = CLASHES();
					contactinit[DIR+2][ORIENTATION+2][0] = -CONTACTS();
					clashinit[DIR+2][ORIENTATION+2][1]=DIR;
					contactinit[DIR+2][ORIENTATION+2][1]=DIR;
					clashinit[DIR+2][ORIENTATION+2][2]=ORIENTATION;
					contactinit[DIR+2][ORIENTATION+2][2]=ORIENTATION;
					chdir(workingdir);
					READER(i,j);
					chdir(writingdir);																									
				}
			}
		
			DIR=0;
			ORIENTATION=0;
			dummie=0;
			for (o=0;o<6;o++)
			{
				for (p=0;p<7;p++)
				{
					restdistchange1=0;
					restdistchange2=0;
					
					//INITIAL OPTIMIZATION OF CLASHES AND CONTACTS
					
					if (clashinit[o][p][0]>=200)
					{
						DIR=contactinit[o][p][1];
						ORIENTATION=contactinit[o][p][2];
						restdist=SCALE_1*restdist;
						CENTRALIZER(restdist);
						restdistchange1=1;
						contactinit[o][p][0] = -CONTACTS();
						clashinit[o][p][0] = CLASHES();
						chdir(workingdir);
						READER(i,j);
						chdir(writingdir);
						restdist=32.59;						
					}
					else if (contactinit[o][p][0]>=-3000)
					{
						DIR=contactinit[o][p][1];
						ORIENTATION=contactinit[o][p][2];
						restdist=SCALE_2*restdist;
						CENTRALIZER(restdist);
						restdistchange2=1;
						contactinit[o][p][0] = -CONTACTS();
						clashinit[o][p][0] = CLASHES();
						chdir(workingdir);
						READER(i,j);
						chdir(writingdir);
						restdist=32.59;						
					}

					printf("%d	%d	%d\n",contactinit[o][p][1],contactinit[o][p][2],contactinit[o][p][0]);
										
					/*if(clashinit[o][p][0]==0)
					{
						shapeinit[o][p][0]=0;
						continue;
					}
					
					else
					{
						shapeinit[o][p][0]=(0.99*contactinit[o][p][0])/(0.01*clashinit[o][p][0]);
					
						if (shapeinit[o][p][0]<dummie)
						{
							dummie=shapeinit[o][p][0];
						}
					}*/
					
					if (contactinit[o][p][0]<dummie)
					{
						dummie=contactinit[o][p][0];
						if(restdistchange1)
						{
							newrestdist1=1;
						}
						else
						{
							newrestdist1=0;
						}
						if(restdistchange2)
						{
							newrestdist2=1;
						}
						else
						{
							newrestdist2=0;
						}
					}
				}
			}
		
			printf("\n\n");
			for (o=0;o<6;o++)
			{
				for (p=0;p<7;p++)
				{
					if (dummie==contactinit[o][p][0])
					{
						DIR=contactinit[o][p][1];
						ORIENTATION=contactinit[o][p][2];
						break;
					}
				}
			}
			
			if(newrestdist1==1)
			{
				restdist=SCALE_1*restdist;
			}
			if(newrestdist2==1)
			{
				restdist=SCALE_2*restdist;
			}			
		
			CENTRALIZER(restdist);
			currentdist=restdist;
			//ExportSnapshot(i,j);
			//sprintf(oo,"MCRUN_I_H_BETA2M_%d_%d_DIR=%d.dat",i,j,(int)DIR);
			//output1=fopen(oo,"w");
			
			//END of INITIAL DIRECTION/ORIENTATION sampling
			
			//BEGINNING of the computation of the INITIAL nº of CONTACTS, CLASHES and INTERMOLECULAR ENERGY
						
			ncontacts = CONTACTS();
				
			if (ncontacts > MAXCONTACT)
			{
				contact_flag=1;
				continue;
			}

			for (N=0; N<ncontacts; N++)
			{
				it[N]=0;
			}
				
			for (N=0; N<MAXATOMNUM; N++)
			{
				strcpy(Atoms[N][0].AtmHBlabel,"NoHB");
				strcpy(Atoms[N][1].AtmHBlabel,"NoHB");
			}
			
			/*ptq0 = (float *) malloc(nmalloc * sizeof(float));
			ptq1 = (float *) malloc(nmalloc * sizeof(float));
			ptct0 = (int *) malloc(nmalloc * sizeof(int));
			ptct1 = (int *) malloc(nmalloc * sizeof(int));
			ptri0 = (int *) malloc(nmalloc * sizeof(int));
			ptri1 = (int *) malloc(nmalloc * sizeof(int));
			ptro0 = (double *) malloc(nmalloc * sizeof(double));
			ptro1 = (double *) malloc(nmalloc * sizeof(double));
			ptrd = (double *) malloc(nmalloc * sizeof(double));
			ptit = (int *) malloc(nmalloc * sizeof(int));
			ptsvdw = (double *) malloc(nmalloc * sizeof(double));
			ptep = (double *) malloc(nmalloc * sizeof(double));
			pthp = (double *) malloc(nmalloc * sizeof(double));
			pthbp = (double *) malloc(nmalloc * sizeof(double));
			ptip = (double *) malloc(nmalloc * sizeof(double));*/
						
			initial_counter = TYPECONTACTS(ct0, ct1, at0, at1, ri0, ri1, ro0, ro1, rd, it, q0, q1, asp0, asp1);
						
			initial_energy_sum = 0.0;
			initial_hp_energy_sum = 0.0;
			initial_ep_energy_sum = 0.0;
			initial_hbp_energy_sum = 0.0;

			hydropathic_flag = 0;	
			for (m=0; m<initial_counter; m++)
			{
				sum_vdW[m] = HARDCORE_R*(ro0[m] + ro1[m]);
				
				/*ct0 = ptct0;
				ct1 = ptct1;
				ri0 = ptri0;
				ri1 = ptri1;
								
				/*q[m][0] = RESIDUECHARGES(ri0, ct0);
				q[m][1] = RESIDUECHARGES(ri1, ct1);
				
				/*sasa_aa[m][0] = SASA_AA[ri0-1];
				sasa_aa[m][1] = SASA_AA[ri1-1];
				
				/*fprintf(output3, "%lf\n", residues_distance[m]);
				fprintf(output4, "%lf\n", sum_vdW[m]);
				fprintf(output5, "%lf\n", 1.4*sum_vdW[m]);
				fprintf(output6, "%lf\n", 2.3*sum_vdW[m]);
				fprintf(output7, "%lf\n", clashes_cutoff[m]);
				fprintf(output8, "%lf\n", contacts_cutoff[m]);*/
				
				if (it[m] == 1 || it[m] == 3 || it[m] == 4) //TYPE OF INTERACTION: 1=electro;3=electro+hidro;4=electro+HB.
				{
					if (rd[m] > sum_vdW[m] & rd[m] < 1.4*sum_vdW[m]) //1st well of the electrostatic potential
					{
						if (q0[m] < 0.0 & q1[m] < 0.0)
						{	
							ep[m] = 0.4;
						}
						else if (q0[m] > 0.0 & q1[m] > 0.0)
						{
							ep[m] = 0.4;
						}
						else if (q0[m] < 0.0 & q1[m] > 0.0)
						{
							ep[m] = -0.4;
						}
						else if (q0[m] > 0.0 & q1[m] < 0.0)
						{
							ep[m] = -0.4;
						}
						else
						{
							ep[m] = 0.0;
						}
					}
					else if (rd[m] > 1.4*sum_vdW[m] & rd[m] < 2.33*sum_vdW[m]) //2nd well of the electrostatic potential
					{
						if (q0[m] < 0.0 & q1[m] < 0.0)
						{
							ep[m] = 0.3*0.4;
						}
						else if (q0[m] > 0.0 & q1[m] > 0.0)
						{
							ep[m] = 0.3*0.4;
						}
						else if (q0[m] < 0.0 & q1[m] > 0.0)
						{
							ep[m] = -0.3*0.4;
						}
						else if (q0[m] > 0.0 & q1[m] < 0.0)
						{
							ep[m] = -0.3*0.4;
						}
						else
						{
							ep[m] = 0.0;
						}
					}
					else
					{
						ep[m] = 0.0;
					}
					hp[m] = 0.0;
					hbp[m] = 0.0;
					if ((it[m] == 3) && (rd[m] < 1.6*sum_vdW[m]))   //TYPE OF INTERACTION: 3=electro+hidro;
					{
						hidropathy_0[m] = -((atomic_solvation_parameters[asp0[m]]*atomic_SASA[ct0[m]][at0[m]])/n)/3000;
						hidropathy_1[m] = -((atomic_solvation_parameters[asp1[m]]*atomic_SASA[ct1[m]][at1[m]])/n)/3000;
						
						if (hidropathy_0[m] <= -0.0539995)
						{
							hidropathy_0[m] = -0.1;
						}
						else if ((hidropathy_0[m] >= -0.0332997) && (hidropathy_0[m] <= 0.0484162))
						{
							hidropathy_0[m] = 0.0;
						}
						else if ((hidropathy_0[m] >= 0.0617161) && (hidropathy_0[m] <= 0.138332))
						{
							hidropathy_0[m] = 0.1;	
						}
						else if (hidropathy_0[m] >= 0.332063)
						{
							hidropathy_0[m] = 0.4;
						}
						
						if (hidropathy_1[m] <= -0.0539995)
						{
							hidropathy_1[m] = -0.1;
						}
						else if ((hidropathy_1[m] >= -0.0332997) && (hidropathy_1[m] <= 0.0484162))
						{
							hidropathy_1[m] = 0.0;
						}
						else if ((hidropathy_1[m] >= 0.0617161) && (hidropathy_1[m] <= 0.138332))
						{
							hidropathy_1[m] = 0.1;	
						}
						else if (hidropathy_1[m] >= 0.332063)
						{
							hidropathy_1[m] = 0.4;
						}
						
						hp[m] = hidropathy_0[m]+hidropathy_1[m];	
					}
						
					if (it[m] == 4)   //TYPE OF INTERACTION: 4=electro+HB;
					{
						hbp[m] = -1.0;
					}
					else
					{
						hbp[m] = 0.0;
					}
				}
				else if ((it[m] == 2) && (rd[m] < 1.6*sum_vdW[m])) //TYPE OF INTERACTION: 2=hidro.
				{
					hidropathy_0[m] = -((atomic_solvation_parameters[asp0[m]]*atomic_SASA[ct0[m]][at0[m]])/n)/3000;
					hidropathy_1[m] = -((atomic_solvation_parameters[asp1[m]]*atomic_SASA[ct1[m]][at1[m]])/n)/3000;
					
					if (hidropathy_0[m] <= -0.0539995)
					{
						hidropathy_0[m] = -0.1;
					}
					else if ((hidropathy_0[m] >= -0.0332997) && (hidropathy_0[m] <= 0.0484162))
					{
						hidropathy_0[m] = 0.0;
					}
					else if ((hidropathy_0[m] >= 0.0617161) && (hidropathy_0[m] <= 0.138332))
					{
						hidropathy_0[m] = 0.1;	
					}
					else if (hidropathy_0[m] >= 0.332063)
					{
						hidropathy_0[m] = 0.4;
					}
			
					if (hidropathy_1[m] <= -0.0539995)
					{
						hidropathy_1[m] = -0.1;
					}
					else if ((hidropathy_1[m] >= -0.0332997) && (hidropathy_1[m] <= 0.0484162))
					{
						hidropathy_1[m] = 0.0;
					}
					else if ((hidropathy_1[m] >= 0.0617161) && (hidropathy_1[m] <= 0.138332))
					{
						hidropathy_1[m] = 0.1;	
					}
					else if (hidropathy_1[m] >= 0.332063)
					{
						hidropathy_1[m] = 0.4;
					}
					
					hp[m] = hidropathy_0[m]+hidropathy_1[m];	
					ep[m] = 0.0;
					hbp[m] = 0.0;
				}
				/*else if ((it[m] == 2) && ((ri0[m] != residue_index_0) || (ri1[m] != residue_index_1)))
				{
					hp[m] = 0.3*hydropathic_inter_residue_potentials[ct0[m]][ct1[m]];
					ep[m] = 0.0;
					hbp[m] = 0.0;
					residue_index_0 = ri0[m];
					residue_index_1 = ri1[m];
					hydropathic_flag = 1;
				}*/
				else
				{
					ep[m] = 0.0;
					hp[m] = 0.0;
					hbp[m] = 0.0;
				}					
				/*fprintf(output3,"electrostatic potential: %lf\n", electrostatic_potential[m]);
				fprintf(output3,"hydropathic inter-residue potential: %lf\n", hydropathic_inter_residue_potentials[contact_type[m][0]][contact_type[m][1]]);*/
				ip[m] = ep[m] + hp[m] + hbp[m];
				/*fprintf(output3,"interaction potential: %lf\n", interaction_potential[m]);*/
				/*f[m][0] = tanh(5*tan(PI*sasa_aa[m][0]/2));
				f[m][1] = tanh(5*tan(PI*sasa_aa[m][1]/2));*/
				initial_energy_sum = initial_energy_sum + /*f[m][0]*f[m][1]**/ip[m];
				initial_hp_energy_sum = initial_hp_energy_sum + hp[m];
				initial_ep_energy_sum = initial_ep_energy_sum + ep[m];
				initial_hbp_energy_sum = initial_hbp_energy_sum + hbp[m];
				/*fprintf(output8,"%d    %d    %d-%d     %lf     %lf\n",m+1, initial_counter, contact_type_initial[m][0], contact_type_initial[m][1],hydropathic_inter_residue_potentials[contact_type_initial[m][0]][contact_type_initial[m][1]], initial_energy_sum); 
				/*printf("%d    %d    %d-%d     %lf     %lf\n",m+1, initial_counter, contact_type_initial[m][0],contact_type_initial[m][1],hydropathic_inter_residue_potentials[contact_type_initial[m][0]][contact_type_initial[m][1]], initial_energy_sum);*/
			}
			/*free(ptq0);
			free(ptq1);
			free(ptct0);
			free(ptct1);
			free(ptri0);
			free(ptri1);
			free(ptro0);
			free(ptro1);
			free(ptrd);
			free(ptit);
			free(ptsvdw);
			free(ptep);
			free(pthp);
			free(pthbp);
			free(ptip);*/								
			printf("%d	%d	%d	%d	%4.3lf	%4.3lf	%4.3lf	%4.3lf	%4.3lf	%4.3lf	%4.3lf	%d	%d       %d	  %d\n",k,i,j,0,(1-TOL)*restdist,currentdist,(1+TOL0)*restdist,initial_energy_sum,initial_hp_energy_sum,initial_ep_energy_sum,initial_hbp_energy_sum,CONTACTS(),initial_counter,CLASHES(),DIR);
			/*fprintf(output7,"%d	%d	%d	%d	%4.3lf	%4.3lf	%4.3lf	%4.3lf	%d       %d	  %d\n",k,i,j,0,(1-TOL)*restdist,currentdist,(1+TOL0)*restdist,initial_energy_sum,CONTACTS(),CLASHES(),DIR);*/
			//fprintf(output1,"%d	%d	%d	%4.3lf	%4.3lf	%4.3lf	    %d      %d      %4.3lf\n",i,j,index,currentdist,restdist,initial_energy_sum,CONTACTS(),CLASHES(),Pacc);
			//fflush(output1);

			//END of the computation of the INITIAL nº of CONTACTS, CLASHES and INTERMOLECULAR ENERGY
			
			//BEGINNING OF MONTE CARLO 
			
			initclash=(double)CLASHES();
			
			accepted=0;
			rejected=0;
			Pacc=0;
			Pdecl=0;
			for (index=1;index<=MCSTEPS;index++)
			{
		
			
				if (currentdist<=(1+0.01)*(1-TOL)*restdist)
				{	
					lastTOL=TOL;
					TOL=TOL+0.01;
					hitwall++;
					if ((1-TOL)*restdist <= 0)
						TOL=lastTOL;
				}
						
			
				//copying information of structure 2 to a stack variable which will be
				//subjected to a MC move
			
				for (m=0;m<nAtom[1];m++)
				{
					Moved[m]=Atoms[m][1];
				}
		
				olddist=currentdist;
			
				//total energy for l-th Markov state
			
				
						
				//=== translating + rotating ====================================
			
				dummie=RandomUniform();
				MCMOVE(dummie);
			
				//===============================================================
			
				//total energy for possible (l+1)-th Markov state
				
				//BEGINNING of the computation of the nº of CONTACTS, CLASHES and INTERMOLECULAR ENERGY BEFORE THE MOVE
				
				ncontacts = CONTACTS();

				if (ncontacts > MAXCONTACT)
				{
					contact_flag=1;
					break;
				}
				
				for (N=0; N<ncontacts; N++)
				{
					it[N]=0;
				}
				
				for (N=0; N<MAXATOMNUM; N++)
				{
					strcpy(Atoms[N][0].AtmHBlabel,"NoHB");
					strcpy(Atoms[N][1].AtmHBlabel,"NoHB");
				}
			
				/*ptq0 = (float *) malloc(nmalloc * sizeof(float));
				ptq1 = (float *) malloc(nmalloc * sizeof(float));
				ptct0 = (int *) malloc(nmalloc * sizeof(int));
				ptct1 = (int *) malloc(nmalloc * sizeof(int));
				ptri0 = (int *) malloc(nmalloc * sizeof(int));
				ptri1 = (int *) malloc(nmalloc * sizeof(int));
				ptro0 = (double *) malloc(nmalloc * sizeof(double));
				ptro1 = (double *) malloc(nmalloc * sizeof(double));
				ptrd = (double *) malloc(nmalloc * sizeof(double));
				ptit = (int *) malloc(nmalloc * sizeof(int));
				ptsvdw = (double *) malloc(nmalloc * sizeof(double));
				ptep = (double *) malloc(nmalloc * sizeof(double));
				pthp = (double *) malloc(nmalloc * sizeof(double));
				pthbp = (double *) malloc(nmalloc * sizeof(double));
				ptip = (double *) malloc(nmalloc * sizeof(double));*/
				
				counter_before_mcmove = TYPECONTACTS(ct0, ct1, at0, at1, ri0, ri1, ro0, ro1, rd, it, q0, q1, asp0, asp1);
				
				energy_sum_before_mcmove = 0.0;
				hydropathic_flag = 0;
				for (m=0; m<counter_before_mcmove; m++)
				{
					sum_vdW[m] = HARDCORE_R*(ro0[m] + ro1[m]);
				
					/*ct0 = ptct0;
					ct1 = ptct1;
					ri0 = ptri0;
					ri1 = ptri1;
								
					/*q[m][0] = RESIDUECHARGES(ri0, ct0);
					q[m][1] = RESIDUECHARGES(ri1, ct1);
				
					/*sasa_aa[m][0] = SASA_AA[ri0-1];
					sasa_aa[m][1] = SASA_AA[ri1-1];
				
					/*fprintf(output3, "%lf\n", residues_distance[m]);
					fprintf(output4, "%lf\n", sum_vdW[m]);
					fprintf(output5, "%lf\n", 1.4*sum_vdW[m]);
					fprintf(output6, "%lf\n", 2.3*sum_vdW[m]);
					fprintf(output7, "%lf\n", clashes_cutoff[m]);
					fprintf(output8, "%lf\n", contacts_cutoff[m]);*/
				
					if (it[m] == 1 || it[m] == 3 || it[m] == 4) //TYPE OF INTERACTION: 1=electro;3=electro+hidro;4=electro+HB.
					{
						if (rd[m] > sum_vdW[m] & rd[m] < 1.4*sum_vdW[m])   //1st well of the electrostatic potential
						{
							if (q0[m] < 0.0 & q1[m] < 0.0)
							{	
								ep[m] = 0.4;
							}
							else if (q0[m] > 0.0 & q1[m] > 0.0)
							{
								ep[m] = 0.4;
							}
							else if (q0[m] < 0.0 & q1[m] > 0.0)
							{
								ep[m] = -0.4;
							}
							else if (q0[m] > 0.0 & q1[m] < 0.0)
							{
								ep[m] = -0.4;
							}
							else
							{
								ep[m] = 0.0;
							}
						}
						else if (rd[m] > 1.4*sum_vdW[m] & rd[m] < 2.33*sum_vdW[m])  //2nd well of the electrostatic potential
						{
							if (q0[m] < 0.0 & q1[m] < 0.0)
							{
								ep[m] = 0.3*0.4;
							}
							else if (q0[m] > 0.0 & q1[m] > 0.0)
							{
								ep[m] = 0.3*0.4;
							}
							else if (q0[m] < 0.0 & q1[m] > 0.0)
							{
								ep[m] = -0.3*0.4;
							}
							else if (q0[m] > 0.0 & q1[m] < 0.0)
							{
								ep[m] = -0.3*0.4;
							}
							else
							{
								ep[m] = 0.0;
							}
						}
						else
						{
							ep[m] = 0.0;
						}
						hp[m] = 0.0;
						hbp[m] = 0.0;
						if ((it[m] == 3) && (rd[m] < 1.6*sum_vdW[m]))  //TYPE OF INTERACTION: 3=electro+hidro;
						{
							hidropathy_0[m] = -((atomic_solvation_parameters[asp0[m]]*atomic_SASA[ct0[m]][at0[m]])/n)/3000;
							hidropathy_1[m] = -((atomic_solvation_parameters[asp1[m]]*atomic_SASA[ct1[m]][at1[m]])/n)/3000;
							
							if (hidropathy_0[m] <= -0.0539995)
							{
								hidropathy_0[m] = -0.1;
							}
							else if ((hidropathy_0[m] >= -0.0332997) && (hidropathy_0[m] <= 0.0484162))
							{
								hidropathy_0[m] = 0.0;
							}
							else if ((hidropathy_0[m] >= 0.0617161) && (hidropathy_0[m] <= 0.138332))
							{
								hidropathy_0[m] = 0.1;	
							}
							else if (hidropathy_0[m] >= 0.332063)
							{
								hidropathy_0[m] = 0.4;
							}
			
							if (hidropathy_1[m] <= -0.0539995)
							{
								hidropathy_1[m] = -0.1;
							}
							else if ((hidropathy_1[m] >= -0.0332997) && (hidropathy_1[m] <= 0.0484162))
							{
								hidropathy_1[m] = 0.0;
							}
							else if ((hidropathy_1[m] >= 0.0617161) && (hidropathy_1[m] <= 0.138332))
							{
								hidropathy_1[m] = 0.1;	
							}
							else if (hidropathy_1[m] >= 0.332063)
							{
								hidropathy_1[m] = 0.4;
							}
							
							hp[m] = hidropathy_0[m]+hidropathy_1[m];	
						}
						
						if (it[m] == 4)   //TYPE OF INTERACTION: 4=electro+HB;
						{
							hbp[m] = -1.0;
						}
						else
						{
							hbp[m] = 0.0;
						}
					}
					else if ((it[m] == 2) && (rd[m] < 1.6*sum_vdW[m])) //TYPE OF INTERACTION: 2=hidro.
					{
						hidropathy_0[m] = -((atomic_solvation_parameters[asp0[m]]*atomic_SASA[ct0[m]][at0[m]])/n)/3000;
						hidropathy_1[m] = -((atomic_solvation_parameters[asp1[m]]*atomic_SASA[ct1[m]][at1[m]])/n)/3000;
							
						if (hidropathy_0[m] <= -0.0539995)
						{
							hidropathy_0[m] = -0.1;
						}
						else if ((hidropathy_0[m] >= -0.0332997) && (hidropathy_0[m] <= 0.0484162))
						{
							hidropathy_0[m] = 0.0;
						}
						else if ((hidropathy_0[m] >= 0.0617161) && (hidropathy_0[m] <= 0.138332))
						{
							hidropathy_0[m] = 0.1;	
						}
						else if (hidropathy_0[m] >= 0.332063)
						{
							hidropathy_0[m] = 0.4;
						}
			
						if (hidropathy_1[m] <= -0.0539995)
						{
							hidropathy_1[m] = -0.1;
						}
						else if ((hidropathy_1[m] >= -0.0332997) && (hidropathy_1[m] <= 0.0484162))
						{
							hidropathy_1[m] = 0.0;
						}
						else if ((hidropathy_1[m] >= 0.0617161) && (hidropathy_1[m] <= 0.138332))
						{
							hidropathy_1[m] = 0.1;	
						}
						else if (hidropathy_1[m] >= 0.332063)
						{
							hidropathy_1[m] = 0.4;
						}
							
						hp[m] = hidropathy_0[m]+hidropathy_1[m];	
						ep[m] = 0.0;
						hbp[m] = 0.0;
					}
					/*else if ((it[m] == 2) && ((ri0[m] != residue_index_0) || (ri1[m] != residue_index_1)))
					{
						hp[m] = 0.3*hydropathic_inter_residue_potentials[ct0[m]][ct1[m]];
						ep[m] = 0.0;
						hbp[m] = 0.0;
						residue_index_0 = ri0[m];
						residue_index_1 = ri1[m];
						hydropathic_flag = 1;
					}*/
					else
					{
						ep[m] = 0.0;
						hp[m] = 0.0;
						hbp[m] = 0.0;
					}					
					/*fprintf(output3,"electrostatic potential before mcmove: %lf\n", electrostatic_potential_before_mcmove[m]);
					fprintf(output3,"hydropathic inter-residue potential before mcmove: %lf\n", hydropathic_inter_residue_potentials[contact_type_before_mcmove[m][0]][contact_type_before_mcmove[m][1]]);*/
					ip[m] = ep[m] + hp[m] + hbp[m];
					/*fprintf(output3,"interaction potential before mcmove: %lf\n", interaction_potential_before_mcmove[m]);*/
					/*f[m][0] = tanh(5*tan(PI*sasa_aa[m][0]/2));
					f[m][1] = tanh(5*tan(PI*sasa_aa[m][1]/2));*/
					energy_sum_before_mcmove = energy_sum_before_mcmove + /*f[m][0]*f[m][1]**/ip[m];
				}
				/*free(ptq0);
				free(ptq1);
				free(ptct0);
				free(ptct1);
				free(ptri0);
				free(ptri1);
				free(ptro0);
				free(ptro1);
				free(ptrd);
				free(ptit);
				free(ptsvdw);
				free(ptep);
				free(pthp);
				free(pthbp);
				free(ptip);*/								
				//END of the computation of the nº of CONTACTS, CLASHES and INTERMOLECULAR ENERGY BEFORE THE MOVE
				
				//BEGINNING of the computation of the nº of CONTACTS, CLASHES and INTERMOLECULAR ENERGY AFTER THE MOVE

				ncontacts = CONTACTSMOVED();

				if (ncontacts > MAXCONTACT)
				{
					contact_flag=1;
					break;
				}
				
				for (N=0; N<ncontacts; N++)
				{
					it[N]=0;
				}
				
				for (N=0; N<MAXATOMNUM; N++)
				{
					strcpy(Atoms[N][0].AtmHBlabel,"NoHB");
					strcpy(Moved[N].AtmHBlabel,"NoHB");
				}
							
				/*ptq0 = (float *) malloc(nmalloc * sizeof(float));
				ptq1 = (float *) malloc(nmalloc * sizeof(float));
				ptct0 = (int *) malloc(nmalloc * sizeof(int));
				ptct1 = (int *) malloc(nmalloc * sizeof(int));
				ptri0 = (int *) malloc(nmalloc * sizeof(int));
				ptri1 = (int *) malloc(nmalloc * sizeof(int));
				ptro0 = (double *) malloc(nmalloc * sizeof(double));
				ptro1 = (double *) malloc(nmalloc * sizeof(double));
				ptrd = (double *) malloc(nmalloc * sizeof(double));
				ptit = (int *) malloc(nmalloc * sizeof(int));
				ptsvdw = (double *) malloc(nmalloc * sizeof(double));
				ptep = (double *) malloc(nmalloc * sizeof(double));
				pthp = (double *) malloc(nmalloc * sizeof(double));
				pthbp = (double *) malloc(nmalloc * sizeof(double));
				ptip = (double *) malloc(nmalloc * sizeof(double));*/

				counter_after_mcmove = TYPECONTACTSMOVED(ct0, ct1, at0, at1, ri0, ri1, ro0, ro1, rd, it, q0, q1, asp0, asp1);
				
				energy_sum_after_mcmove = 0.0;
				hydropathic_flag = 0;
				for (m=0; m<counter_after_mcmove; m++)
				{
					sum_vdW[m] = HARDCORE_R*(ro0[m] + ro1[m]);
				
					/*ct0 = ptct0;
					ct1 = ptct1;
					ri0 = ptri0;
					ri1 = ptri1;
								
					/*q[m][0] = RESIDUECHARGES(ri0, ct0);
					q[m][1] = RESIDUECHARGES(ri1, ct1);
				
					/*sasa_aa[m][0] = SASA_AA[ri0-1];
					sasa_aa[m][1] = SASA_AA[ri1-1];
				
					/*fprintf(output3, "%lf\n", residues_distance[m]);
					fprintf(output4, "%lf\n", sum_vdW[m]);
					fprintf(output5, "%lf\n", 1.4*sum_vdW[m]);
					fprintf(output6, "%lf\n", 2.3*sum_vdW[m]);
					fprintf(output7, "%lf\n", clashes_cutoff[m]);
					fprintf(output8, "%lf\n", contacts_cutoff[m]);*/
				
					if (it[m] == 1 || it[m] == 3 || it[m] == 4)
					{
						if (rd[m] > sum_vdW[m] & rd[m] < 1.4*sum_vdW[m])
						{
							if (q0[m] < 0.0 & q1[m] < 0.0)
							{	
								ep[m] = 0.4;
							}
							else if (q0[m] > 0.0 & q1[m] > 0.0)
							{
								ep[m] = 0.4;
							}
							else if (q0[m] < 0.0 & q1[m] > 0.0)
							{
								ep[m] = -0.4;
							}
							else if (q0[m] > 0.0 & q1[m] < 0.0)
							{
								ep[m] = -0.4;
							}
							else
							{
								ep[m] = 0.0;
							}
						}
						else if (rd[m] > 1.4*sum_vdW[m] & rd[m] < 2.33*sum_vdW[m])
						{
							if (q0[m] < 0.0 & q1[m] < 0.0)
							{
								ep[m] = 0.3*0.4;
							}
							else if (q0[m] > 0.0 & q1[m] > 0.0)
							{
								ep[m] = 0.3*0.4;
							}
							else if (q0[m] < 0.0 & q1[m] > 0.0)
							{
								ep[m] = -0.3*0.4;
							}
							else if (q0[m] > 0.0 & q1[m] < 0.0)
							{
								ep[m] = -0.3*0.4;
							}
							else
							{
								ep[m] = 0.0;
							}
						}
						else
						{
							ep[m] = 0.0;
						}
						hp[m] = 0.0;
						hbp[m] = 0.0;
						if ((it[m] == 3) && (rd[m] < 1.6*sum_vdW[m]))
						{
							hidropathy_0[m] = -((atomic_solvation_parameters[asp0[m]]*atomic_SASA[ct0[m]][at0[m]])/n)/3000;
							hidropathy_1[m] = -((atomic_solvation_parameters[asp1[m]]*atomic_SASA[ct1[m]][at1[m]])/n)/3000;
							
							if (hidropathy_0[m] <= -0.0539995)
							{
								hidropathy_0[m] = -0.1;
							}
							else if ((hidropathy_0[m] >= -0.0332997) && (hidropathy_0[m] <= 0.0484162))
							{
								hidropathy_0[m] = 0.0;
							}
							else if ((hidropathy_0[m] >= 0.0617161) && (hidropathy_0[m] <= 0.138332))
							{
								hidropathy_0[m] = 0.1;	
							}
							else if (hidropathy_0[m] >= 0.332063)
							{
								hidropathy_0[m] = 0.4;
							}
			
							if (hidropathy_1[m] <= -0.0539995)
							{
								hidropathy_1[m] = -0.1;
							}
							else if ((hidropathy_1[m] >= -0.0332997) && (hidropathy_1[m] <= 0.0484162))
							{
								hidropathy_1[m] = 0.0;
							}
							else if ((hidropathy_1[m] >= 0.0617161) && (hidropathy_1[m] <= 0.138332))
							{
								hidropathy_1[m] = 0.1;	
							}
							else if (hidropathy_1[m] >= 0.332063)
							{
								hidropathy_1[m] = 0.4;
							}
							
							hp[m] = hidropathy_0[m]+hidropathy_1[m];	
						}
						
						if (it[m] == 4)
						{
							hbp[m] = -1.0;
						}
						else
						{
							hbp[m] = 0.0;
						}
					}
					else if ((it[m] == 2) && (rd[m] < 1.6*sum_vdW[m]))
					{
						hidropathy_0[m] = -((atomic_solvation_parameters[asp0[m]]*atomic_SASA[ct0[m]][at0[m]])/n)/3000;
						hidropathy_1[m] = -((atomic_solvation_parameters[asp1[m]]*atomic_SASA[ct1[m]][at1[m]])/n)/3000;
							
						if (hidropathy_0[m] <= -0.0539995)
						{
							hidropathy_0[m] = -0.1;
						}
						else if ((hidropathy_0[m] >= -0.0332997) && (hidropathy_0[m] <= 0.0484162))
						{
							hidropathy_0[m] = 0.0;
						}
						else if ((hidropathy_0[m] >= 0.0617161) && (hidropathy_0[m] <= 0.138332))
						{
							hidropathy_0[m] = 0.1;	
						}
						else if (hidropathy_0[m] >= 0.332063)
						{
							hidropathy_0[m] = 0.4;
						}
			
						if (hidropathy_1[m] <= -0.0539995)
						{
							hidropathy_1[m] = -0.1;
						}
						else if ((hidropathy_1[m] >= -0.0332997) && (hidropathy_1[m] <= 0.0484162))
						{
							hidropathy_1[m] = 0.0;
						}
						else if ((hidropathy_1[m] >= 0.0617161) && (hidropathy_1[m] <= 0.138332))
						{
							hidropathy_1[m] = 0.1;	
						}
						else if (hidropathy_1[m] >= 0.332063)
						{
							hidropathy_1[m] = 0.4;
						}
							
						hp[m] = hidropathy_0[m]+hidropathy_1[m];	
						ep[m] = 0.0;
						hbp[m] = 0.0;
					}
					/*else if ((it[m] == 2) && ((ri0[m] != residue_index_0) || (ri1[m] != residue_index_1)))
					{
						hp[m] = 0.3*hydropathic_inter_residue_potentials[ct0[m]][ct1[m]];
						ep[m] = 0.0;
						hbp[m] = 0.0;
						residue_index_0 = ri0[m];
						residue_index_1 = ri1[m];
						hydropathic_flag = 1;
					}*/
					else
					{
						ep[m] = 0.0;
						hp[m] = 0.0;
						hbp[m] = 0.0;
					}					
					/*fprintf(output3,"electrostatic potential after mcmove: %lf\n", electrostatic_potential_after_mcmove[m]);
					fprintf(output3,"hydropathic inter-residue potential after mcmove: %lf\n", hydropathic_inter_residue_potentials[contact_type_after_mcmove[m][0]][contact_type_after_mcmove[m][1]]);*/
					ip[m] = ep[m] + hp[m] + hbp[m];
					/*fprintf(output3,"interaction potential after mcmove: %lf\n", interaction_potential_after_mcmove[m]);*/
					/*f[m][0] = tanh(5*tan(PI*sasa_aa[m][0]/2));
					f[m][1] = tanh(5*tan(PI*sasa_aa[m][1]/2));*/
					energy_sum_after_mcmove = energy_sum_after_mcmove + /*f[m][0]*f[m][1]**/ip[m];
				}
				/*free(ptq0);
				free(ptq1);
				free(ptct0);
				free(ptct1);
				free(ptri0);
				free(ptri1);
				free(ptro0);
				free(ptro1);
				free(ptrd);
				free(ptit);
				free(ptsvdw);
				free(ptep);
				free(pthp);
				free(pthbp);
				free(ptip);*/								
				//END of the computation of the nº of CONTACTS, CLASHES and INTERMOLECULAR ENERGY AFTER THE MOVE
								
				Econtold = energy_sum_before_mcmove;
				Econtnew = energy_sum_after_mcmove;
				contactnew = (double)CONTACTSMOVED();
				contactold = (double)CONTACTS(); 
				clashnew = (double)CLASHESMOVED();
				clashold = (double)CLASHES();
				
						
				//METROPOLIS CRITERION
			
				if (SCALE_CLASHES*(clashnew - clashold) < 0.)
				{
					/*double Energydiffmetropolis;
					Energydiffmetropolis = Econtnew - Econtold;
					printf("energydiff is %lf\n", Energydiffmetropolis);*/
					if ((1-SCALE_CLASHES)*(Econtnew - Econtold) < 0.)
					{
						accepted++;
						for (m = 0; m<nAtom[1]; m++)
						{
							Atoms[m][1] = Moved[m];
						}
						e = k-1;
						mc = index-1;
						/*diffcontactsenergy[e][mc] = Econtnew - Econtold;*/
					}
					else
					{
						if (exp(-((1-SCALE_CLASHES)*(Econtnew - Econtold)) / Tm) >= RandomUniform())
						{
							accepted++;
							for (m = 0; m<nAtom[1]; m++)
							{
								Atoms[m][1] = Moved[m];
							}
							e = k-1;
							mc = index-1;
							/*diffcontactsenergy[e][mc] = Econtnew - Econtold;*/
						}
						else
						{
							rejected++;
							e = k-1;
							mc = index-1;
							/*diffcontactsenergy[e][mc] = Econtnew - Econtold;*/
						}
					}
				}
				else
				{
					/*double Clashesdiffmetropolis;
					Clashesdiffmetropolis = clashnew - clashold;
					printf("clashesdiff is %lf\n", Clashesdiffmetropolis);*/
					if ((exp(-(SCALE_CLASHES*(clashnew - clashold)) / Tm)) >= RandomUniform())
					{
						if ((1-SCALE_CLASHES)*(Econtnew - Econtold) < 0)
						{
							accepted++;
							for (m = 0; m<nAtom[1]; m++)
							{
								Atoms[m][1] = Moved[m];
							}
							e = k-1;
							mc = index-1;
							/*diffcontactsenergy[e][mc] = Econtnew - Econtold;*/
						}
						else
						{
							if (exp(-((1-SCALE_CLASHES)*(Econtnew - Econtold)) / Tm) >= RandomUniform())
							{
								accepted++;
								for (m = 0; m<nAtom[1]; m++)
								{
									Atoms[m][1] = Moved[m];
								}
								e = k-1;
								mc = index-1;
								/*diffcontactsenergy[e][mc] = Econtnew - Econtold;*/
							}
							else
							{
								rejected++;
								e = k-1;
								mc = index-1;
								/*diffcontactsenergy[e][mc] = Econtnew - Econtold;*/
							}
						}
					}
					else
					{
						rejected++;
						e = k-1;
						mc = index-1;
						/*diffcontactsenergy[e][mc] = Econtnew - Econtold;*/
					}
				}
				/*printf("accepted\n");
				printf("%d\n", accepted);
				printf("rejected\n");
				printf("%d\n", rejected);*/
				Pacc=(double)accepted/((double)index);
				Pdecl=(double)rejected/((double)index);
			
				//METROPOLIS CONVERGENCE OPTIMIZATION
			
				if (Pacc>0.55)
				{
					mcinternal1=mcinternal1+0.00001*mcinternal1;
					mcinternal2=mcinternal2+0.00001*mcinternal2;
				}
				else
				{
					if (Pacc<0.45)
					{
						mcinternal1=mcinternal1-0.00001*mcinternal1;
						mcinternal2=mcinternal2-0.00001*mcinternal2;
					}
				}
				
				if (mcinternal1>=1.0 & mcinternal2>=1.0)
				{
					mcinternal1=0.025;
					mcinternal2=0.025;
				}
				
				if (mcinternal1<=0.01 & mcinternal2<=0.01)
				{
					mcinternal1=0.025;
					mcinternal2=0.025;
				}								
				
				//PROTOCOL PART
				
				//COMPUTATION OF CONTACTS, CLASHES AND ENERGIES AFTER METROPOLIS CRITERION
				
				ncontacts = CONTACTS();

				if (ncontacts > MAXCONTACT)
				{
					contact_flag=1;
					break;
				}
				
				for (N=0; N<ncontacts; N++)
				{
					it[N]=0;
				}
				
				for (N=0; N<MAXATOMNUM; N++)
				{
					strcpy(Atoms[N][0].AtmHBlabel,"NoHB");
					strcpy(Atoms[N][1].AtmHBlabel,"NoHB");
				}
							
				/*ptq0 = (float *) malloc(nmalloc * sizeof(float));
				ptq1 = (float *) malloc(nmalloc * sizeof(float));
				ptct0 = (int *) malloc(nmalloc * sizeof(int));
				ptct1 = (int *) malloc(nmalloc * sizeof(int));
				ptri0 = (int *) malloc(nmalloc * sizeof(int));
				ptri1 = (int *) malloc(nmalloc * sizeof(int));
				ptro0 = (double *) malloc(nmalloc * sizeof(double));
				ptro1 = (double *) malloc(nmalloc * sizeof(double));
				ptrd = (double *) malloc(nmalloc * sizeof(double));
				ptit = (int *) malloc(nmalloc * sizeof(int));
				ptsvdw = (double *) malloc(nmalloc * sizeof(double));
				ptep = (double *) malloc(nmalloc * sizeof(double));
				pthp = (double *) malloc(nmalloc * sizeof(double));
				pthbp = (double *) malloc(nmalloc * sizeof(double));
				ptip = (double *) malloc(nmalloc * sizeof(double));*/

				final_counter = TYPECONTACTS(ct0, ct1, at0, at1, ri0, ri1, ro0, ro1, rd, it, q0, q1, asp0, asp1);
				
				final_energy_sum = 0.0;
				final_hp_energy_sum = 0.0;
				final_ep_energy_sum = 0.0;
				final_hbp_energy_sum = 0.0;
				
				hydropathic_flag = 0;
				for (m=0; m<final_counter; m++)
				{
					sum_vdW[m] = HARDCORE_R*(ro0[m] + ro1[m]);
				
					/*ct0 = ptct0;
					ct1 = ptct1;
					ri0 = ptri0;
					ri1 = ptri1;
								
					/*q[m][0] = RESIDUECHARGES(ri0, ct0);
					q[m][1] = RESIDUECHARGES(ri1, ct1);
				
					/*sasa_aa[m][0] = SASA_AA[ri0-1];
					sasa_aa[m][1] = SASA_AA[ri1-1];
				
					/*fprintf(output3, "%lf\n", residues_distance[m]);
					fprintf(output4, "%lf\n", sum_vdW[m]);
					fprintf(output5, "%lf\n", 1.4*sum_vdW[m]);
					fprintf(output6, "%lf\n", 2.3*sum_vdW[m]);
					fprintf(output7, "%lf\n", clashes_cutoff[m]);
					fprintf(output8, "%lf\n", contacts_cutoff[m]);*/
				
					if (it[m] == 1 || it[m] == 3 || it[m] == 4)
					{
						if (rd[m] > sum_vdW[m] & rd[m] < 1.4*sum_vdW[m])
						{
							if (q0[m] < 0.0 & q1[m] < 0.0)
							{	
								ep[m] = 0.4;
							}
							else if (q0[m] > 0.0 & q1[m] > 0.0)
							{
								ep[m] = 0.4;
							}
							else if (q0[m] < 0.0 & q1[m] > 0.0)
							{
								ep[m] = -0.4;
							}
							else if (q0[m] > 0.0 & q1[m] < 0.0)
							{
								ep[m] = -0.4;
							}
							else
							{
								ep[m] = 0.0;
							}
						}
						else if (rd[m] > 1.4*sum_vdW[m] & rd[m] < 2.33*sum_vdW[m])
						{
							if (q0[m] < 0.0 & q1[m] < 0.0)
							{
								ep[m] = 0.3*0.4;
							}
							else if (q0[m] > 0.0 & q1[m] > 0.0)
							{
								ep[m] = 0.3*0.4;
							}
							else if (q0[m] < 0.0 & q1[m] > 0.0)
							{
								ep[m] = -0.3*0.4;
							}
							else if (q0[m] > 0.0 & q1[m] < 0.0)
							{
								ep[m] = -0.3*0.4;
							}
							else
							{
								ep[m] = 0.0;
							}
						}
						else
						{
							ep[m] = 0.0;
						}
						hp[m] = 0.0;
						hbp[m] = 0.0;
						if ((it[m] == 3) && (rd[m] < 1.6*sum_vdW[m]))
						{
							hidropathy_0[m] = -((atomic_solvation_parameters[asp0[m]]*atomic_SASA[ct0[m]][at0[m]])/n)/3000;
							hidropathy_1[m] = -((atomic_solvation_parameters[asp1[m]]*atomic_SASA[ct1[m]][at1[m]])/n)/3000;
							
							if (hidropathy_0[m] <= -0.0539995)
							{
								hidropathy_0[m] = -0.1;
							}
							else if ((hidropathy_0[m] >= -0.0332997) && (hidropathy_0[m] <= 0.0484162))
							{
								hidropathy_0[m] = 0.0;
							}
							else if ((hidropathy_0[m] >= 0.0617161) && (hidropathy_0[m] <= 0.138332))
							{
								hidropathy_0[m] = 0.1;	
							}
							else if (hidropathy_0[m] >= 0.332063)
							{
								hidropathy_0[m] = 0.4;
							}
			
							if (hidropathy_1[m] <= -0.0539995)
							{
								hidropathy_1[m] = -0.1;
							}
							else if ((hidropathy_1[m] >= -0.0332997) && (hidropathy_1[m] <= 0.0484162))
							{
								hidropathy_1[m] = 0.0;
							}
							else if ((hidropathy_1[m] >= 0.0617161) && (hidropathy_1[m] <= 0.138332))
							{
								hidropathy_1[m] = 0.1;	
							}
							else if (hidropathy_1[m] >= 0.332063)
							{
								hidropathy_1[m] = 0.4;
							}
							
							hp[m] = hidropathy_0[m]+hidropathy_1[m];	
						}
						
						if (it[m] == 4)
						{
							hbp[m] = -1.0;
						}
						else
						{
							hbp[m] = 0.0;
						}
					}
					else if ((it[m] == 2) && (rd[m] < 1.6*sum_vdW[m]))
					{
						hidropathy_0[m] = -((atomic_solvation_parameters[asp0[m]]*atomic_SASA[ct0[m]][at0[m]])/n)/3000;
						hidropathy_1[m] = -((atomic_solvation_parameters[asp1[m]]*atomic_SASA[ct1[m]][at1[m]])/n)/3000;
							
						if (hidropathy_0[m] <= -0.0539995)
						{
							hidropathy_0[m] = -0.1;
						}
						else if ((hidropathy_0[m] >= -0.0332997) && (hidropathy_0[m] <= 0.0484162))
						{
							hidropathy_0[m] = 0.0;
						}
						else if ((hidropathy_0[m] >= 0.0617161) && (hidropathy_0[m] <= 0.138332))
						{
							hidropathy_0[m] = 0.1;	
						}
						else if (hidropathy_0[m] >= 0.332063)
						{
							hidropathy_0[m] = 0.4;
						}
		
						if (hidropathy_1[m] <= -0.0539995)
						{
							hidropathy_1[m] = -0.1;
						}
						else if ((hidropathy_1[m] >= -0.0332997) && (hidropathy_1[m] <= 0.0484162))
						{
							hidropathy_1[m] = 0.0;
						}
						else if ((hidropathy_1[m] >= 0.0617161) && (hidropathy_1[m] <= 0.138332))
						{
							hidropathy_1[m] = 0.1;	
						}
						else if (hidropathy_1[m] >= 0.332063)
						{
							hidropathy_1[m] = 0.4;
						}
							
						hp[m] = hidropathy_0[m]+hidropathy_1[m];	
						ep[m] = 0.0;
						hbp[m] = 0.0;
					}
					/*else if ((it[m] == 2) && ((ri0[m] != residue_index_0) || (ri1[m] != residue_index_1)))
					{
						hp[m] = 0.3*hydropathic_inter_residue_potentials[ct0[m]][ct1[m]];
						ep[m] = 0.0;
						hbp[m] = 0.0;
						residue_index_0 = ri0[m];
						residue_index_1 = ri1[m];
						hydropathic_flag = 1;
					}*/
					else
					{
						ep[m] = 0.0;
						hp[m] = 0.0;
						hbp[m] = 0.0;
					}					
					/*fprintf(output3,"electrostatic potential final: %lf\n", electrostatic_potential_final[m]);
					fprintf(output3,"hydropathic inter-residue potential final: %lf\n", hydropathic_inter_residue_potentials[contact_type_final[m][0]][contact_type_final[m][1]]);*/
					ip[m] = ep[m] + hp[m] + hbp[m];
					/*fprintf(output3,"interaction potential final: %lf\n", interaction_potential_final[m]);*/
					/*f[m][0] = tanh(5*tan(PI*sasa_aa[m][0]/2));
					f[m][1] = tanh(5*tan(PI*sasa_aa[m][1]/2));*/
					final_energy_sum = final_energy_sum + /*f[m][0]*f[m][1]**/ip[m];
					final_hp_energy_sum = final_hp_energy_sum + hp[m];
					final_ep_energy_sum = final_ep_energy_sum + ep[m];
					final_hbp_energy_sum = final_hbp_energy_sum + hbp[m];
				}
				/*free(ptq0);
				free(ptq1);
				free(ptct0);
				free(ptct1);
				free(ptri0);
				free(ptri1);
				free(ptro0);
				free(ptro1);
				free(ptrd);
				free(ptit);
				free(ptsvdw);
				free(ptep);
				free(pthp);
				free(pthbp);
				free(ptip);*/								

				/*if (final_energy_sum <= -10.0)
				{
					fprintf(output1,"%d	%d	%d	%4.3lf	%4.3lf	%4.3lf	%d      %d      %4.3lf\n",i,j,index,currentdist,restdist,final_energy_sum,CONTACTS(),CLASHES(),Pacc);
					fflush(output1);
					ExportSnapshot(i,j);
					break;
				}*/
				
				//COMPUTATION OF CONTACTS, CLASHES AND ENERGIES AT PROTSTEP MC STEP				
				if ((index % PROTSTEP == 0)&&(index>0))
				{
					ncontacts = CONTACTS();

					if (ncontacts > MAXCONTACT)
					{
						contact_flag=1;
						break;
					}
				
					for (N=0; N<ncontacts; N++)
					{
						it[N]=0;
					}
				
					for (N=0; N<MAXATOMNUM; N++)
					{
						strcpy(Atoms[N][0].AtmHBlabel,"NoHB");
						strcpy(Atoms[N][1].AtmHBlabel,"NoHB");
					}

					/*ptq0 = (float *) malloc(nmalloc * sizeof(float));
					ptq1 = (float *) malloc(nmalloc * sizeof(float));
					ptct0 = (int *) malloc(nmalloc * sizeof(int));
					ptct1 = (int *) malloc(nmalloc * sizeof(int));
					ptri0 = (int *) malloc(nmalloc * sizeof(int));
					ptri1 = (int *) malloc(nmalloc * sizeof(int));
					ptro0 = (double *) malloc(nmalloc * sizeof(double));
					ptro1 = (double *) malloc(nmalloc * sizeof(double));
					ptrd = (double *) malloc(nmalloc * sizeof(double));
					ptit = (int *) malloc(nmalloc * sizeof(int));
					ptsvdw = (double *) malloc(nmalloc * sizeof(double));
					ptep = (double *) malloc(nmalloc * sizeof(double));
					pthp = (double *) malloc(nmalloc * sizeof(double));
					pthbp = (double *) malloc(nmalloc * sizeof(double));
					ptip = (double *) malloc(nmalloc * sizeof(double));*/

					final_counter = TYPECONTACTS(ct0, ct1, at0, at1, ri0, ri1, ro0, ro1, rd, it, q0, q1, asp0, asp1);
				
					final_energy_sum = 0.0;
					final_hp_energy_sum = 0.0;
					final_ep_energy_sum = 0.0;
					final_hbp_energy_sum = 0.0;

					hydropathic_flag = 0;
					for (m=0; m<final_counter; m++)
					{
						sum_vdW[m] = HARDCORE_R*(ro0[m] + ro1[m]);
				
						/*ct0 = ptct0;
						ct1 = ptct1;
						ri0 = ptri0;
						ri1 = ptri1;
								
						/*q[m][0] = RESIDUECHARGES(ri0, ct0);
						q[m][1] = RESIDUECHARGES(ri1, ct1);
				
						/*sasa_aa[m][0] = SASA_AA[ri0-1];
						sasa_aa[m][1] = SASA_AA[ri1-1];
				
						/*fprintf(output3, "%lf\n", residues_distance[m]);
						fprintf(output4, "%lf\n", sum_vdW[m]);
						fprintf(output5, "%lf\n", 1.4*sum_vdW[m]);
						fprintf(output6, "%lf\n", 2.3*sum_vdW[m]);
						fprintf(output7, "%lf\n", clashes_cutoff[m]);
						fprintf(output8, "%lf\n", contacts_cutoff[m]);*/
				
						if (it[m] == 1 || it[m] == 3 || it[m] == 4)
						{
							if (rd[m] > sum_vdW[m] & rd[m] < 1.4*sum_vdW[m])
							{
								if (q0[m] < 0.0 & q1[m] < 0.0)
								{	
									ep[m] = 0.4;
								}
								else if (q0[m] > 0.0 & q1[m] > 0.0)
								{
									ep[m] = 0.4;
								}
								else if (q0[m] < 0.0 & q1[m] > 0.0)
								{
									ep[m] = -0.4;
								}
								else if (q0[m] > 0.0 & q1[m] < 0.0)
								{
									ep[m] = -0.4;
								}
								else
								{
									ep[m] = 0.0;
								}
							}
							else if (rd[m] > 1.4*sum_vdW[m] & rd[m] < 2.33*sum_vdW[m])
							{
								if (q0[m] < 0.0 & q1[m] < 0.0)
								{
									ep[m] = 0.3*0.4;
								}
								else if (q0[m] > 0.0 & q1[m] > 0.0)
								{
									ep[m] = 0.3*0.4;
								}
								else if (q0[m] < 0.0 & q1[m] > 0.0)
								{
									ep[m] = -0.3*0.4;
								}
								else if (q0[m] > 0.0 & q1[m] < 0.0)
								{
									ep[m] = -0.3*0.4;
								}
								else
								{
									ep[m] = 0.0;
								}
							}
							else
							{
								ep[m] = 0.0;
							}
							hp[m] = 0.0;
							hbp[m] = 0.0;
							if ((it[m] == 3) && (rd[m] < 1.6*sum_vdW[m]))
							{
								hidropathy_0[m] = -((atomic_solvation_parameters[asp0[m]]*atomic_SASA[ct0[m]][at0[m]])/n)/3000;
								hidropathy_1[m] = -((atomic_solvation_parameters[asp1[m]]*atomic_SASA[ct1[m]][at1[m]])/n)/3000;
							
								if (hidropathy_0[m] <= -0.0539995)
								{
									hidropathy_0[m] = -0.1;
								}
								else if ((hidropathy_0[m] >= -0.0332997) && (hidropathy_0[m] <= 0.0484162))
								{
									hidropathy_0[m] = 0.0;
								}
								else if ((hidropathy_0[m] >= 0.0617161) && (hidropathy_0[m] <= 0.138332))
								{
									hidropathy_0[m] = 0.1;	
								}
								else if (hidropathy_0[m] >= 0.332063)
								{
									hidropathy_0[m] = 0.4;
								}
			
								if (hidropathy_1[m] <= -0.0539995)
								{
									hidropathy_1[m] = -0.1;
								}
								else if ((hidropathy_1[m] >= -0.0332997) && (hidropathy_1[m] <= 0.0484162))
								{
									hidropathy_1[m] = 0.0;
								}
								else if ((hidropathy_1[m] >= 0.0617161) && (hidropathy_1[m] <= 0.138332))
								{
									hidropathy_1[m] = 0.1;	
								}
								else if (hidropathy_1[m] >= 0.332063)
								{
									hidropathy_1[m] = 0.4;
								}
							
								hp[m] = hidropathy_0[m]+hidropathy_1[m];	
							}
						
							if (it[m] == 4)
							{
								hbp[m] = -1.0;
							}
							else
							{
								hbp[m] = 0.0;
							}
						}
						else if ((it[m] == 2) && (rd[m] < 1.6*sum_vdW[m]))
						{
							hidropathy_0[m] = -((atomic_solvation_parameters[asp0[m]]*atomic_SASA[ct0[m]][at0[m]])/n)/3000;
							hidropathy_1[m] = -((atomic_solvation_parameters[asp1[m]]*atomic_SASA[ct1[m]][at1[m]])/n)/3000;
							
							if (hidropathy_0[m] <= -0.0539995)
							{
								hidropathy_0[m] = -0.1;
							}
							else if ((hidropathy_0[m] >= -0.0332997) && (hidropathy_0[m] <= 0.0484162))
							{
								hidropathy_0[m] = 0.0;
							}
							else if ((hidropathy_0[m] >= 0.0617161) && (hidropathy_0[m] <= 0.138332))
							{
								hidropathy_0[m] = 0.1;	
							}
							else if (hidropathy_0[m] >= 0.332063)
							{
								hidropathy_0[m] = 0.4;
							}
			
							if (hidropathy_1[m] <= -0.0539995)
							{
								hidropathy_1[m] = -0.1;
							}
							else if ((hidropathy_1[m] >= -0.0332997) && (hidropathy_1[m] <= 0.0484162))
							{
								hidropathy_1[m] = 0.0;
							}
							else if ((hidropathy_1[m] >= 0.0617161) && (hidropathy_1[m] <= 0.138332))
							{
								hidropathy_1[m] = 0.1;	
							}
							else if (hidropathy_1[m] >= 0.332063)
							{
								hidropathy_1[m] = 0.4;
							}
							
							hp[m] = hidropathy_0[m]+hidropathy_1[m];	
							ep[m] = 0.0;
							hbp[m] = 0.0;
						}
						/*else if ((it[m] == 2) && ((ri0[m] != residue_index_0) || (ri1[m] != residue_index_1)))
						{
							hp[m] = 0.3*hydropathic_inter_residue_potentials[ct0[m]][ct1[m]];
							ep[m] = 0.0;
							hbp[m] = 0.0;
							residue_index_0 = ri0[m];
							residue_index_1 = ri1[m];
							hydropathic_flag = 1;
						}*/
						else
						{
							ep[m] = 0.0;
							hp[m] = 0.0;
							hbp[m] = 0.0;
						}					
						/*fprintf(output3,"electrostatic potential final: %lf\n", electrostatic_potential_final[m]);
						fprintf(output3,"hydropathic inter-residue potential final: %lf\n", hydropathic_inter_residue_potentials[contact_type_final[m][0]][contact_type_final[m][1]]);*/
						ip[m] = ep[m] + hp[m] + hbp[m];
						/*fprintf(output3,"interaction potential final: %lf\n", interaction_potential_final[m]);*/
						/*f[m][0] = tanh(5*tan(PI*sasa_aa[m][0]/2));
						f[m][1] = tanh(5*tan(PI*sasa_aa[m][1]/2));*/
						final_energy_sum = final_energy_sum + /*f[m][0]*f[m][1]**/ip[m];
						final_hp_energy_sum = final_hp_energy_sum + hp[m];
						final_ep_energy_sum = final_ep_energy_sum + ep[m];
						final_hbp_energy_sum = final_hbp_energy_sum + hbp[m];
					}
					/*free(ptq0);
					free(ptq1);
					free(ptct0);
					free(ptct1);
					free(ptri0);
					free(ptri1);
					free(ptro0);
					free(ptro1);
					free(ptrd);
					free(ptit);
					free(ptsvdw);
					free(ptep);
					free(pthp);
					free(pthbp);
					free(ptip);*/								
					/*fprintf(output7,"%d    %d    %d     %d   %4.3lf  %4.3lf  %4.3lf  %4.3lf  %d    %d   %4.3lf   %7.6lf    %7.6lf  %d\n",k,i,j,index,(1-TOL)*restdist,currentdist,(1+TOL0)*restdist,final_energy_sum,CONTACTS(),CLASHES(),Pacc, mcinternal1,mcinternal2,hitwall);*/
					//fprintf(output1,"%d	%d	%d	%4.3lf	%4.3lf	%4.3lf	%d  %d %4.3lf\n",i,j,index,currentdist,restdist,final_energy_sum,CONTACTS(),CLASHES(),Pacc);
					//fflush(output1);
					printf("%d    %d    %d     %d   %4.3lf  %4.3lf  %4.3lf  %4.3lf  %4.3lf	%4.3lf	%4.3lf	%d	%d    %d   %4.3lf   %7.6lf    %7.6lf  %d\n",k,i,j,index,(1-TOL)*restdist,currentdist,(1+TOL0)*restdist,final_energy_sum,final_hp_energy_sum,final_ep_energy_sum,final_hbp_energy_sum,CONTACTS(),final_counter,CLASHES(),Pacc, mcinternal1,mcinternal2,hitwall);
				}
				//COMPUTATION OF CONTACTS, CLASHES AND ENERGIES AT LAST MC STEP
				if (index == MCSTEPS)
				{
					ncontacts = CONTACTS();

					if (ncontacts > MAXCONTACT)
					{
						contact_flag=1;
						break;
					}
				
					for (N=0; N<ncontacts; N++)
					{
						it[N]=0;
					}
				
					for (N=0; N<MAXATOMNUM; N++)
					{
						strcpy(Atoms[N][0].AtmHBlabel,"NoHB");
						strcpy(Atoms[N][1].AtmHBlabel,"NoHB");
					}

					/*ptq0 = (float *) malloc(nmalloc * sizeof(float));
					ptq1 = (float *) malloc(nmalloc * sizeof(float));
					ptct0 = (int *) malloc(nmalloc * sizeof(int));
					ptct1 = (int *) malloc(nmalloc * sizeof(int));
					ptri0 = (int *) malloc(nmalloc * sizeof(int));
					ptri1 = (int *) malloc(nmalloc * sizeof(int));
					ptro0 = (double *) malloc(nmalloc * sizeof(double));
					ptro1 = (double *) malloc(nmalloc * sizeof(double));
					ptrd = (double *) malloc(nmalloc * sizeof(double));
					ptit = (int *) malloc(nmalloc * sizeof(int));
					ptsvdw = (double *) malloc(nmalloc * sizeof(double));
					ptep = (double *) malloc(nmalloc * sizeof(double));
					pthp = (double *) malloc(nmalloc * sizeof(double));
					pthbp = (double *) malloc(nmalloc * sizeof(double));
					ptip = (double *) malloc(nmalloc * sizeof(double));*/
					
					final_counter = TYPECONTACTS(ct0, ct1, at0, at1, ri0, ri1, ro0, ro1, rd, it, q0, q1, asp0, asp1);
				
					final_energy_sum = 0.0;
					final_hp_energy_sum = 0.0;
					final_ep_energy_sum = 0.0;
					final_hbp_energy_sum = 0.0;

					hydropathic_flag = 0;
					for (m=0; m<final_counter; m++)
					{
						sum_vdW[m] = HARDCORE_R*(ro0[m] + ro1[m]);
				
						/*ct0 = ptct0;
						ct1 = ptct1;
						ri0 = ptri0;
						ri1 = ptri1;
								
						/*q[m][0] = RESIDUECHARGES(ri0, ct0);
						q[m][1] = RESIDUECHARGES(ri1, ct1);
				
						/*sasa_aa[m][0] = SASA_AA[ri0-1];
						sasa_aa[m][1] = SASA_AA[ri1-1];
				
						/*fprintf(output3, "%lf\n", residues_distance[m]);
						fprintf(output4, "%lf\n", sum_vdW[m]);
						fprintf(output5, "%lf\n", 1.4*sum_vdW[m]);
						fprintf(output6, "%lf\n", 2.3*sum_vdW[m]);
						fprintf(output7, "%lf\n", clashes_cutoff[m]);
						fprintf(output8, "%lf\n", contacts_cutoff[m]);*/
				
						if (it[m] == 1 || it[m] == 3 || it[m] == 4)
						{
							if (rd[m] > sum_vdW[m] & rd[m] < 1.4*sum_vdW[m])
							{
								if (q0[m] < 0.0 & q1[m] < 0.0)
								{	
									ep[m] = 0.4;
								}
								else if (q0[m] > 0.0 & q1[m] > 0.0)
								{
									ep[m] = 0.4;
								}
								else if (q0[m] < 0.0 & q1[m] > 0.0)
								{
									ep[m] = -0.4;
								}
								else if (q0[m] > 0.0 & q1[m] < 0.0)
								{
									ep[m] = -0.4;
								}
								else
								{
									ep[m] = 0.0;
								}
							}
							else if (rd[m] > 1.4*sum_vdW[m] & rd[m] < 2.33*sum_vdW[m])
							{
								if (q0[m] < 0.0 & q1[m] < 0.0)
								{
									ep[m] = 0.3*0.4;
								}
								else if (q0[m] > 0.0 & q1[m] > 0.0)
								{
									ep[m] = 0.3*0.4;
								}
								else if (q0[m] < 0.0 & q1[m] > 0.0)
								{
									ep[m] = -0.3*0.4;
								}
								else if (q0[m] > 0.0 & q1[m] < 0.0)
								{
									ep[m] = -0.3*0.4;
								}
								else
								{
									ep[m] = 0.0;
								}
							}
							else
							{
								ep[m] = 0.0;
							}
							hp[m] = 0.0;
							hbp[m] = 0.0;
							if ((it[m] == 3) && (rd[m] < 1.6*sum_vdW[m]))
							{
								hidropathy_0[m] = -((atomic_solvation_parameters[asp0[m]]*atomic_SASA[ct0[m]][at0[m]])/n)/3000;
								hidropathy_1[m] = -((atomic_solvation_parameters[asp1[m]]*atomic_SASA[ct1[m]][at1[m]])/n)/3000;								
								
								if (hidropathy_0[m] <= -0.0539995)
								{
									hidropathy_0[m] = -0.1;
								}
								else if ((hidropathy_0[m] >= -0.0332997) && (hidropathy_0[m] <= 0.0484162))
								{
									hidropathy_0[m] = 0.0;
								}
								else if ((hidropathy_0[m] >= 0.0617161) && (hidropathy_0[m] <= 0.138332))
								{
									hidropathy_0[m] = 0.1;	
								}
								else if (hidropathy_0[m] >= 0.332063)
								{
									hidropathy_0[m] = 0.4;
								}
			
								if (hidropathy_1[m] <= -0.0539995)
								{
									hidropathy_1[m] = -0.1;
								}
								else if ((hidropathy_1[m] >= -0.0332997) && (hidropathy_1[m] <= 0.0484162))
								{
									hidropathy_1[m] = 0.0;
								}
								else if ((hidropathy_1[m] >= 0.0617161) && (hidropathy_1[m] <= 0.138332))
								{
									hidropathy_1[m] = 0.1;	
								}
								else if (hidropathy_1[m] >= 0.332063)
								{
									hidropathy_1[m] = 0.4;
								}
							
								hp[m] = hidropathy_0[m]+hidropathy_1[m];
							}
						
							if (it[m] == 4)
							{
								hbp[m] = -1.0;
							}
							else
							{
								hbp[m] = 0.0;
							}
						}
						else if ((it[m] == 2) && (rd[m] < 1.6*sum_vdW[m]))
						{
							hidropathy_0[m] = -((atomic_solvation_parameters[asp0[m]]*atomic_SASA[ct0[m]][at0[m]])/n)/3000;
							hidropathy_1[m] = -((atomic_solvation_parameters[asp1[m]]*atomic_SASA[ct1[m]][at1[m]])/n)/3000;
							
							if (hidropathy_0[m] <= -0.0539995)
							{
								hidropathy_0[m] = -0.1;
							}
							else if ((hidropathy_0[m] >= -0.0332997) && (hidropathy_0[m] <= 0.0484162))
							{
								hidropathy_0[m] = 0.0;
							}
							else if ((hidropathy_0[m] >= 0.0617161) && (hidropathy_0[m] <= 0.138332))
							{
								hidropathy_0[m] = 0.1;	
							}
							else if (hidropathy_0[m] >= 0.332063)
							{
								hidropathy_0[m] = 0.4;
							}
			
							if (hidropathy_1[m] <= -0.0539995)
							{
								hidropathy_1[m] = -0.1;
							}
							else if ((hidropathy_1[m] >= -0.0332997) && (hidropathy_1[m] <= 0.0484162))
							{
								hidropathy_1[m] = 0.0;
							}
							else if ((hidropathy_1[m] >= 0.0617161) && (hidropathy_1[m] <= 0.138332))
							{
								hidropathy_1[m] = 0.1;	
							}
							else if (hidropathy_1[m] >= 0.332063)
							{
								hidropathy_1[m] = 0.4;
							}
							
							hp[m] = hidropathy_0[m]+hidropathy_1[m];	
							ep[m] = 0.0;
							hbp[m] = 0.0;
						}
						/*else if ((it[m] == 2) && ((ri0[m] != residue_index_0) || (ri1[m] != residue_index_1)))
						{
							hp[m] = 0.3*hydropathic_inter_residue_potentials[ct0[m]][ct1[m]];
							ep[m] = 0.0;
							hbp[m] = 0.0;
							residue_index_0 = ri0[m];
							residue_index_1 = ri1[m];
							hydropathic_flag = 1;
						}*/
						else
						{
							ep[m] = 0.0;
							hp[m] = 0.0;
							hbp[m] = 0.0;
						}					
						/*fprintf(output3,"electrostatic potential final: %lf\n", electrostatic_potential_final[m]);
						fprintf(output3,"hydropathic inter-residue potential final: %lf\n", hydropathic_inter_residue_potentials[contact_type_final[m][0]][contact_type_final[m][1]]);*/
						ip[m] = ep[m] + hp[m] + hbp[m];
						/*fprintf(output3,"interaction potential final: %lf\n", interaction_potential_final[m]);*/
						/*f[m][0] = tanh(5*tan(PI*sasa_aa[m][0]/2));
						f[m][1] = tanh(5*tan(PI*sasa_aa[m][1]/2));*/
						final_energy_sum = final_energy_sum + /*f[m][0]*f[m][1]**/ip[m];
						final_hp_energy_sum = final_hp_energy_sum + hp[m];
						final_ep_energy_sum = final_ep_energy_sum + ep[m];
						final_hbp_energy_sum = final_hbp_energy_sum + hbp[m];							
					}
					/*free(ptq0);
					free(ptq1);
					free(ptct0);
					free(ptct1);
					free(ptri0);
					free(ptri1);
					free(ptro0);
					free(ptro1);
					free(ptrd);
					free(ptit);
					free(ptsvdw);
					free(ptep);
					free(pthp);
					free(pthbp);
					free(ptip);*/								
					fprintf(output1,"%d	%d	%d	%4.3lf	%4.3lf	%4.3lf	%4.3lf	%4.3lf	%4.3lf	%d 	%d     %d      %4.3lf\n",i,j,index,currentdist,restdist,final_energy_sum,final_hp_energy_sum,final_ep_energy_sum,final_hbp_energy_sum,CONTACTS(),final_counter,CLASHES(),Pacc);
					fflush(output1);
					
					ExportSnapshot(i,j);
				}
		
			}
			//END OF MONTE CARLO
			/*fprintf(output7,"\n");*/
			
			e = k-1;
			
			CCONTACTS[e] = CONTACTS();
			CCLASHES[e] = CLASHES();
			ENERGY[e] = final_energy_sum;
			CONTACTS_SUM = CONTACTS_SUM + CCONTACTS[e];
			CLASHES_SUM = CLASHES_SUM + CCLASHES[e];
			ENERGY_SUM = ENERGY_SUM + ENERGY[e];
			
			if (k % 100 == 0)
			{
				CONTACTS_MEAN = (CONTACTS_SUM/(double)k);
	
				CLASHES_MEAN = (CLASHES_SUM/(double)k);
				
				ENERGY_MEAN = (ENERGY_SUM/(double)k);
				
				for (e=0; e < k; e++)
				{
					contacts_diff_sum = contacts_diff_sum + (CCONTACTS[e] - CONTACTS_MEAN)*(CCONTACTS[e] - CONTACTS_MEAN);
				}
				
				VARIANCE_CONTACTS = contacts_diff_sum/(double)k;
	
				VARIANCE_CONTACTS_MEAN = VARIANCE_CONTACTS/(double)k;
				
				for (e=0; e < k; e++)
				{
					energy_diff_sum = energy_diff_sum + (ENERGY[e] - ENERGY_MEAN)*(ENERGY[e] - ENERGY_MEAN);
				}
				
				VARIANCE_ENERGY = energy_diff_sum/(double)k;
	
				VARIANCE_ENERGY_MEAN = VARIANCE_ENERGY/(double)k;
				
				fprintf(output2,"%lf     %lf     %lf     %lf     %lf\n", CONTACTS_MEAN, CLASHES_MEAN, ENERGY_MEAN, VARIANCE_CONTACTS_MEAN, VARIANCE_ENERGY_MEAN);
	
			}
										
	}
	
	/*for (e=0; e < PAIRS; e++)
	{
		for (mc=0; mc < MCSTEPS; mc++)
		{
			diffcontactsenergy_sum[e] = diffcontactsenergy_sum[e] + diffcontactsenergy[e][mc];
			fprintf(output8,"%d     %d     %lf     \n", e+1, mc+1, diffcontactsenergy[e][mc]);
		}
		diffcontactsenergy_mean[e] = diffcontactsenergy_sum[e]/(double)MCSTEPS;
		fprintf(output8,"%d     %lf     \n", e+1, diffcontactsenergy_mean[e]);
	}*/
	
	
	fclose(output1);
	
	fclose(output2);
	
	/*fclose(output3);
	
	fclose(output4);

	fclose(output5);
	
	fclose(output6);
	
	fclose(output7);
	
	fclose(output8);*/
	
		
	return 0;
	
}


void MCMOVE(double a)
{
	int i,j,k;
	double translrand;
	double rotrand1;
	double rotrand2;
	long double helper1,helper2,helper3,helper4;
	double rotvector[3];
	double rotmat[3][3];
	double oldmoved_x[MAXATOMNUM];
	double oldmoved_y[MAXATOMNUM];
	double oldmoved_z[MAXATOMNUM];		
	
	long double tester1,tester2;
	
	if (a<0.5)
	{
	
		//translational MC-submove along x-axis
		
		while (1)
		{
			translrand=RandomUniform()-0.5;
			/*printf("%lf\n", mcinternal1); */
			currentdist=currentdist+mcinternal1*translrand;
			if ((currentdist<=(1+TOL0)*restdist)&&(currentdist>=(1-TOL)*restdist))
				break;
		}
	
		if (abbs(DIR)==1)
		{
			for (j=0;j<nAtom[1];j++)
			{
				Moved[j].x = Moved[j].x + mcinternal1*translrand;
			}
			
		}
		else
		{
			if (abbs(DIR)==2)
			{
				for (j=0;j<nAtom[1];j++)
				{
					Moved[j].y = Moved[j].y + mcinternal1*translrand;
				}
				
			}
			else
			{
				if (abbs(DIR)==3)
				{
					for (j=0;j<nAtom[1];j++)
					{
						Moved[j].z = Moved[j].z + mcinternal1*translrand;
					}
					
				}
			}
		}
	}
	else
	{
	
		/*for (i=0;i<nAtom[0];i++)
		{
			controldist[i][0]=sqrt((Moved[i].x - CM2[0])*(Moved[i].x - CM2[0])+(Moved[i].y - CM2[1])*(Moved[i].y - CM2[1])+(Moved[i].z - CM2[2])*(Moved[i].z - CM2[2]));
		}*/
	
	
		//rotational MC-submove
	
		while (1)
		{
			rotrand1=PI*RandomUniform();			//theta on unit sphere in [0,pi]
			rotrand2=2*PI*RandomUniform();			//phi on unit sphere in [0,2pi]
	
			helper1=sin(rotrand1);
			helper2=cos(rotrand1);
			helper3=sin(rotrand2);
			helper4=cos(rotrand2);
	
			rotvector[0]=helper1*helper4;					//random unit vector in spherical coord. defining rotation axis
			rotvector[1]=helper1*helper3;
			rotvector[2]=helper4;
	
			tester1=rotvector[0]*rotvector[0]+rotvector[1]*rotvector[1]+rotvector[2]*rotvector[2];
	
			rotvector[0]=rotvector[0]/tester1;
			rotvector[1]=rotvector[1]/tester1;
			rotvector[2]=rotvector[2]/tester1;
	
			tester2=rotvector[0]*rotvector[0]+rotvector[1]*rotvector[1]+rotvector[2]*rotvector[2];
	
			if (tester2==1.)
			{
				break;
			}	
		}
		/*printf("%lf\n", mcinternal2);*/
	
		helper1=2*PI*mcinternal2*RandomUniform(); 		//alpha
	
		helper2=sin(helper1);
		helper3=cos(helper1);
		helper4=1-helper3;
	
		//first matrix index: line
		//second matrix index: column
	
		//setting up rotation matrix
		rotmat[0][0] = rotvector[0]*rotvector[0]*helper4+helper3;
		rotmat[0][1] = rotvector[0]*rotvector[1]*helper4-rotvector[2]*helper2;
		rotmat[0][2] = rotvector[0]*rotvector[2]*helper4+rotvector[1]*helper2;
		rotmat[1][0] = rotvector[0]*rotvector[1]*helper4+rotvector[2]*helper2;
		rotmat[1][1] = rotvector[1]*rotvector[1]*helper4+helper3;
		rotmat[1][2] = rotvector[1]*rotvector[2]*helper4-rotvector[0]*helper2;
		rotmat[2][0] = rotvector[0]*rotvector[2]*helper4-rotvector[1]*helper2;
		rotmat[2][1] = rotvector[1]*rotvector[2]*helper4+rotvector[0]*helper2;
		rotmat[2][2] = rotvector[2]*rotvector[2]*helper4+helper3;
		
		//recalculating the coordinates of the structure stored in Moved in relation to its centre of mass
	
		for (j=0;j<3;j++)
		{
			MEAN_R[j][1]=0;
		}
		
		for (j=0;j<nAtom[1];j++)
		{
			MEAN_R[0][1]=MEAN_R[0][1]+Moved[j].x;
			MEAN_R[1][1]=MEAN_R[1][1]+Moved[j].y;
			MEAN_R[2][1]=MEAN_R[2][1]+Moved[j].z;
		}
		
		for (j=0;j<3;j++)
		{
			MEAN_R[j][1]=MEAN_R[j][1]/nAtom[1];
		}
		
		for (i=0;i<nAtom[1];i++)
		{
			RotMoved[i].x = Moved[i].x - MEAN_R[0][1];
			RotMoved[i].y = Moved[i].y - MEAN_R[1][1];
			RotMoved[i].z = Moved[i].z - MEAN_R[2][1];
		}

		//performing rotation on structure stored in RotMoved via matrix multiplication
	
		for (i=0;i<nAtom[1];i++)
		{
			oldmoved_x[i] = RotMoved[i].x;
			oldmoved_y[i] = RotMoved[i].y;
			oldmoved_z[i] = RotMoved[i].z;
			RotMoved[i].x = rotmat[0][0]*oldmoved_x[i] + rotmat[0][1]*oldmoved_y[i] + rotmat[0][2]*oldmoved_z[i];
			RotMoved[i].y = rotmat[1][0]*oldmoved_x[i] + rotmat[1][1]*oldmoved_y[i] + rotmat[1][2]*oldmoved_z[i];
			RotMoved[i].z = rotmat[2][0]*oldmoved_x[i] + rotmat[2][1]*oldmoved_y[i] + rotmat[2][2]*oldmoved_z[i];
		}
		
		//recalculating the coordinates of the structure stored in RotMoved in relation to the origin of the cartesian axes
		
		for (i=0;i<nAtom[1];i++)
		{
			Moved[i].x = RotMoved[i].x + MEAN_R[0][1];
			Moved[i].y = RotMoved[i].y + MEAN_R[1][1];
			Moved[i].z = RotMoved[i].z + MEAN_R[2][1];
		}
	}
}


int TYPECONTACTSMOVED(int ct0[MAXCONTACT], int ct1[MAXCONTACT], int at0[MAXCONTACT], int at1[MAXCONTACT], int ri0[MAXCONTACT], int ri1[MAXCONTACT], double ro0[MAXCONTACT], double ro1[MAXCONTACT], double rd[MAXCONTACT], int it[MAXCONTACT], float q0[MAXCONTACT], float q1[MAXCONTACT], int asp0[MAXCONTACT], int asp1[MAXCONTACT])
{
	FILE *test;
	
	int i, k, l, n, o, p, r, maxcontactflag, counter, len1, len2, residue_index_0, residue_index_1, s, t, ss, tt;
	int AtomIndex_1, AtomIndex_2;
	int electrostatic;
	float dDA, dc1, dc2, dc3;
	/*float sum_inter_residue_potential = 0.0;
	float inter_residue_potential[MAXCONTACT];
	char s[256];
	
	/*sprintf(s, "test_pairs-energies_pos-mcmove.dat");
	test = fopen(s,"a");*/
	
	/*for (i = 0; i < MAXCONTACT; i++)
	{
		inter_residue_potential[i] = 0.0;
	}*/
	
	for (k = 0; k<nAtom[0]; k++)
	{
		strcpy(Atoms[k][0].AtmHBlabel, "NoHB");
	}
	
	for (l = 0; l<nAtom[1]; l++)
	{
		strcpy(Moved[l].AtmHBlabel, "NoHB");
	}
	
	counter = 0;
	for (i = Atoms[0][0].NodeIndex; i <= Atoms[nAtom[0] - 1][0].NodeIndex; i++)
	{
		//picking coordinates of residue with index i.
		
		/*if (contactflag)
		{
			break;
		}*/

		len1 = 0;
		s = 0;
		for (k = 0; k<nAtom[0]; k++)
		{
			if (Atoms[k][0].NodeIndex == i)
			{
				stack1[len1].x = Atoms[k][0].x;
				stack1[len1].y = Atoms[k][0].y;
				stack1[len1].z = Atoms[k][0].z;
				stack1[len1].r0 = Atoms[k][0].r0;
				stack1[len1].Res_Type = Atoms[k][0].Res_Type;
				stack1[len1].Atom_Type = Atoms[k][0].Atom_Type;
				stack1[len1].NodeIndex = Atoms[k][0].NodeIndex;
				stack1[len1].AtomIndex = Atoms[k][0].AtomIndex;
				strcpy(stack1[len1].AtmName, Atoms[k][0].AtmName);
				len1++;	
			}
		}
		for (k = Moved[0].NodeIndex; k <= Moved[nAtom[1] - 1].NodeIndex; k++)
		{
			//picking coordinates of residue with index k 
			
			/*if (contactflag)
			{
				break;
			}*/

			len2 = 0;
			t = 0;
			for (l = 0; l<nAtom[1]; l++)
			{
				if (Moved[l].NodeIndex == k)
				{
					stack2[len2].x = Moved[l].x;
					stack2[len2].y = Moved[l].y;
					stack2[len2].z = Moved[l].z;
					stack2[len2].r0 = Moved[l].r0;
					stack2[len2].Res_Type = Moved[l].Res_Type;
					stack2[len2].Atom_Type = Moved[l].Atom_Type;
					stack2[len2].NodeIndex = Moved[l].NodeIndex;
					stack2[len2].AtomIndex = Moved[l].AtomIndex;
					strcpy(stack2[len2].AtmName, Moved[l].AtmName);
					len2++;		
				}
			}
					
			for (n=0;n<len1;n++)
			{
				/*if (contactflag)
				{
					break;
				}*/
				for (o=0;o<len2;o++)
				{
					dist.r01 = HARDCORE_R*(stack1[n].r0 + stack2[o].r0);
					dist.r02 = CONTACT_R*(stack1[n].r0 + stack2[o].r0);
					dist.rsq = (stack2[o].x - stack1[n].x)*(stack2[o].x - stack1[n].x) + (stack2[o].y - stack1[n].y)*(stack2[o].y - stack1[n].y) + (stack2[o].z - stack1[n].z)*(stack2[o].z - stack1[n].z);
					
					
					if ((dist.rsq>=dist.r01*dist.r01)&&(dist.rsq<=dist.r02*dist.r02))
					{
						for (p = 0; p < LENCHARGESFILE; p++)
						{
							if ((strcmp(stack1[n].AtmName,Charges[p].AtmName) == 0) && (stack1[n].Res_Type == Charges[p].Res_Type))
							{
								for (r = 0; r < LENCHARGESFILE; r++)
								{
									if ((strcmp(stack2[o].AtmName,Charges[r].AtmName) == 0) && (stack2[o].Res_Type == Charges[r].Res_Type))
									{
										q0[counter] = Charges[p].Atm_Charge; //GROMOS54A7 point charges of each interacting atom
										q1[counter] = Charges[r].Atm_Charge; //GROMOS54A7 point charges of each interacting atom
									}
								}
							}							
						}
						if (q0[counter] < 0.0 & q1[counter] < 0.0)
						{	
							electrostatic++;
							it[counter] = 1; /*possible electrostatic interaction*/
							ct0[counter] = stack1[n].Res_Type;
							ct1[counter] = stack2[o].Res_Type;
							ro0[counter] = stack1[n].r0;
							ro1[counter] = stack2[o].r0;
							ri0[counter] = stack1[n].NodeIndex;
							ri1[counter] = stack2[o].NodeIndex;
							rd[counter] = sqrt(dist.rsq);
							
							//HYDROPHOBIC interactions between non-hydrogen side-chain atoms
							
							if (((strcmp(stack1[n].AtmName,"N") != 0) && (strcmp(stack1[n].AtmName,"H") != 0) && (strcmp(stack1[n].AtmName,"H1") != 0) &&
							(strcmp(stack1[n].AtmName,"H2") != 0) && (strcmp(stack1[n].AtmName,"H3") != 0) && (strcmp(stack1[n].AtmName,"CA") != 0) &&
							(strcmp(stack1[n].AtmName,"C") != 0) && (strcmp(stack1[n].AtmName,"O") != 0) && (strcmp(stack1[n].AtmName,"O1") != 0) &&
							(strcmp(stack1[n].AtmName,"O2") != 0)) && ((strcmp(stack2[o].AtmName,"N") != 0) && (strcmp(stack2[o].AtmName,"H") != 0) && (strcmp(stack2[o].AtmName,"H1") != 0) &&
							(strcmp(stack2[o].AtmName,"H2") != 0) && (strcmp(stack2[o].AtmName,"H3") != 0) && (strcmp(stack2[o].AtmName,"CA") != 0) &&
							(strcmp(stack2[o].AtmName,"C") != 0) && (strcmp(stack2[o].AtmName,"O") != 0) && (strcmp(stack2[o].AtmName,"O1") != 0) && (strcmp(stack2[o].AtmName,"O2") != 0)))
							{
								if (strncmp(stack1[n].AtmName,"H", 1) != 0 && strncmp(stack2[o].AtmName,"H", 1) != 0)
								{
									it[counter] = 3; /*possible electrostatic and hydropathic interaction*/
									ct0[counter] = stack1[n].Res_Type;
									ct1[counter] = stack2[o].Res_Type;
									at0[counter] = stack1[n].Atom_Type;
									at1[counter] = stack2[o].Atom_Type;
									
									//Assigning the ATOMIC SOLVATION PARAMETERS(Cummings, 1995) to the different atom types
									
									if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") == 0) || (strcmp(stack1[n].AtmName,"NZ") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 4;
									}
									else if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") != 0) && (strcmp(stack1[n].AtmName,"NZ") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") == 0) || (strcmp(stack1[n].AtmName,"OE2") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 3;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") != 0) && (strcmp(stack1[n].AtmName,"OE2") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if (strncmp(stack1[n].AtmName,"C", 1) == 0)
									{
										stack1[n].ASP_Type = 0;
									}
									else if (strncmp(stack1[n].AtmName,"S", 1) == 0)
									{
										stack1[n].ASP_Type = 2;
									}
									else
									{
										stack1[n].ASP_Type = 5;
									}
									
									if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") == 0) || (strcmp(stack2[o].AtmName,"NZ") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 4;
									}
									else if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") != 0) && (strcmp(stack2[o].AtmName,"NZ") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") == 0) || (strcmp(stack2[o].AtmName,"OE2") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 3;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") != 0) && (strcmp(stack2[o].AtmName,"OE2") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if (strncmp(stack2[o].AtmName,"C", 1) == 0)
									{
										stack2[o].ASP_Type = 0;
									}
									else if (strncmp(stack2[o].AtmName,"S", 1) == 0)
									{
										stack2[o].ASP_Type = 2;
									}									
									else
									{
										stack2[o].ASP_Type = 5;
									}
									asp0[counter] = stack1[n].ASP_Type;
									asp1[counter] = stack2[o].ASP_Type;
									ro0[counter] = stack1[n].r0;
									ro1[counter] = stack2[o].r0;
									ri0[counter] = stack1[n].NodeIndex;
									ri1[counter] = stack2[o].NodeIndex;
									rd[counter] = sqrt(dist.rsq);
								}
							}
							counter++;																
						}
						else if (q0[counter] > 0.0 & q1[counter] > 0.0)
						{
							electrostatic++;
							it[counter] = 1; /*possible electrostatic interaction*/
							ct0[counter] = stack1[n].Res_Type;
							ct1[counter] = stack2[o].Res_Type;
							ro0[counter] = stack1[n].r0;
							ro1[counter] = stack2[o].r0;
							ri0[counter] = stack1[n].NodeIndex;
							ri1[counter] = stack2[o].NodeIndex;
							rd[counter] = sqrt(dist.rsq);
							
							//HYDROPHOBIC interactions between non-hydrogen side-chain atoms
							
							if (((strcmp(stack1[n].AtmName,"N") != 0) && (strcmp(stack1[n].AtmName,"H") != 0) && (strcmp(stack1[n].AtmName,"H1") != 0) &&
							(strcmp(stack1[n].AtmName,"H2") != 0) && (strcmp(stack1[n].AtmName,"H3") != 0) && (strcmp(stack1[n].AtmName,"CA") != 0) &&
							(strcmp(stack1[n].AtmName,"C") != 0) && (strcmp(stack1[n].AtmName,"O") != 0) && (strcmp(stack1[n].AtmName,"O1") != 0) &&
							(strcmp(stack1[n].AtmName,"O2") != 0)) && ((strcmp(stack2[o].AtmName,"N") != 0) && (strcmp(stack2[o].AtmName,"H") != 0) && (strcmp(stack2[o].AtmName,"H1") != 0) &&
							(strcmp(stack2[o].AtmName,"H2") != 0) && (strcmp(stack2[o].AtmName,"H3") != 0) && (strcmp(stack2[o].AtmName,"CA") != 0) &&
							(strcmp(stack2[o].AtmName,"C") != 0) && (strcmp(stack2[o].AtmName,"O") != 0) && (strcmp(stack2[o].AtmName,"O1") != 0) && (strcmp(stack2[o].AtmName,"O2") != 0)))
							{
								if (strncmp(stack1[n].AtmName,"H", 1) != 0 && strncmp(stack2[o].AtmName,"H", 1) != 0)
								{
									it[counter] = 3; /*possible electrostatic and hydropathic interaction*/
									ct0[counter] = stack1[n].Res_Type;
									ct1[counter] = stack2[o].Res_Type;
									at0[counter] = stack1[n].Atom_Type;
									at1[counter] = stack2[o].Atom_Type;
									
									//Assigning the ATOMIC SOLVATION PARAMETERS(Cummings, 1995) to the different atom types
									
									if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") == 0) || (strcmp(stack1[n].AtmName,"NZ") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 4;
									}
									else if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") != 0) && (strcmp(stack1[n].AtmName,"NZ") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") == 0) || (strcmp(stack1[n].AtmName,"OE2") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 3;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") != 0) && (strcmp(stack1[n].AtmName,"OE2") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if (strncmp(stack1[n].AtmName,"C", 1) == 0)
									{
										stack1[n].ASP_Type = 0;
									}
									else if (strncmp(stack1[n].AtmName,"S", 1) == 0)
									{
										stack1[n].ASP_Type = 2;
									}
									else
									{
										stack1[n].ASP_Type = 5;
									}
									
									if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") == 0) || (strcmp(stack2[o].AtmName,"NZ") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 4;
									}
									else if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") != 0) && (strcmp(stack2[o].AtmName,"NZ") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") == 0) || (strcmp(stack2[o].AtmName,"OE2") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 3;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") != 0) && (strcmp(stack2[o].AtmName,"OE2") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if (strncmp(stack2[o].AtmName,"C", 1) == 0)
									{
										stack2[o].ASP_Type = 0;
									}
									else if (strncmp(stack2[o].AtmName,"S", 1) == 0)
									{
										stack2[o].ASP_Type = 2;
									}									
									else
									{
										stack2[o].ASP_Type = 5;
									}
									asp0[counter] = stack1[n].ASP_Type;
									asp1[counter] = stack2[o].ASP_Type;
									ro0[counter] = stack1[n].r0;
									ro1[counter] = stack2[o].r0;
									ri0[counter] = stack1[n].NodeIndex;
									ri1[counter] = stack2[o].NodeIndex;
									rd[counter] = sqrt(dist.rsq);
								}
							}
							counter++;																
						}
						else if (q0[counter] < 0.0 & q1[counter] > 0.0)
						{
							electrostatic++;
							it[counter] = 1; /*possible electrostatic interaction*/
							ct0[counter] = stack1[n].Res_Type;
							ct1[counter] = stack2[o].Res_Type;
							ro0[counter] = stack1[n].r0;
							ro1[counter] = stack2[o].r0;
							ri0[counter] = stack1[n].NodeIndex;
							ri1[counter] = stack2[o].NodeIndex;
							rd[counter] = sqrt(dist.rsq);
							
							//ASSESSING the presence of HYDROGEN BONDS: CHECKING if the DONOR-ACCEPTOR DISTANCE IS <=3.2 Angstrom
							
							if (((strncmp(stack1[n].AtmName,"O",1) == 0) || (strncmp(stack1[n].AtmName,"N",1) == 0) || (strncmp(stack1[n].AtmName,"S",1) == 0)) && (strncmp(stack2[o].AtmName,"H",1) == 0))
							{
								if ((strcmp(stack2[o].AtmName,"H2") == 0) || (strcmp(stack2[o].AtmName,"HH12") == 0) || (strcmp(stack2[o].AtmName,"HH22") == 0) || (strcmp(stack2[o].AtmName,"HD22") == 0) || (strcmp(stack2[o].AtmName,"HE22") == 0) || (strcmp(stack2[o].AtmName,"HZ2") == 0))
								{
									if (strncmp(stack2[o-2].AtmName,"N",1) == 0)
									{	
										AtomIndex_1 = stack1[n].AtomIndex;
										AtomIndex_2 = stack2[o-2].AtomIndex;
										if ((strcmp(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Moved[AtomIndex_2-1].AtmHBlabel,"HB") != 0))
										{
											dDA = sqrt((stack2[o-2].x - stack1[n].x)*(stack2[o-2].x - stack1[n].x) + (stack2[o-2].y -
											stack1[n].y)*(stack2[o-2].y - stack1[n].y) + (stack2[o-2].z - stack1[n].z)*(stack2[o-2].z -
											stack1[n].z));					
											if (dDA<=3.2)
											{
												it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
												strcpy(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB");
												strcpy(Moved[AtomIndex_2-1].AtmHBlabel,"HB");
											}
										}	
									}
								}
								else if ((strcmp(stack2[o].AtmName,"H3") == 0) || (strcmp(stack2[o].AtmName,"HZ3") == 0))
								{
									if (strncmp(stack2[o-3].AtmName,"N",1) == 0)
									{
										AtomIndex_1 = stack1[n].AtomIndex;
										AtomIndex_2 = stack2[o-3].AtomIndex;
										if ((strcmp(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Moved[AtomIndex_2-1].AtmHBlabel,"HB") != 0))
										{
											dDA = sqrt((stack2[o-3].x - stack1[n].x)*(stack2[o-3].x - stack1[n].x) + (stack2[o-3].y -
											stack1[n].y)*(stack2[o-3].y - stack1[n].y) + (stack2[o-3].z - stack1[n].z)*(stack2[o-3].z -
											stack1[n].z));					
											if (dDA<=3.2)
											{
												it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
												strcpy(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB");
												strcpy(Moved[AtomIndex_2-1].AtmHBlabel,"HB"); 
											}
										}
									}	
								}
								else
								{
									if ((strncmp(stack2[o-1].AtmName,"N",1) == 0) || (strncmp(stack2[o-1].AtmName,"O",1) == 0))
									{
										AtomIndex_1 = stack1[n].AtomIndex;
										AtomIndex_2 = stack2[o-1].AtomIndex;
										if ((strcmp(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Moved[AtomIndex_2-1].AtmHBlabel,"HB") != 0))
										{
											dDA = sqrt((stack2[o-1].x - stack1[n].x)*(stack2[o-1].x - stack1[n].x) + (stack2[o-1].y -
											stack1[n].y)*(stack2[o-1].y - stack1[n].y) + (stack2[o-1].z - stack1[n].z)*(stack2[o-1].z -
											stack1[n].z));					
											if (dDA<=3.2)
											{
												it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
												strcpy(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB");
												strcpy(Moved[AtomIndex_2-1].AtmHBlabel,"HB"); 
											}
										}
									}
								}									
							}
							
							//HYDROPHOBIC interactions between non-hydrogen side-chain atoms
							
							if (((strcmp(stack1[n].AtmName,"N") != 0) && (strcmp(stack1[n].AtmName,"H") != 0) && (strcmp(stack1[n].AtmName,"H1") != 0) &&
							(strcmp(stack1[n].AtmName,"H2") != 0) && (strcmp(stack1[n].AtmName,"H3") != 0) && (strcmp(stack1[n].AtmName,"CA") != 0) &&
							(strcmp(stack1[n].AtmName,"C") != 0) && (strcmp(stack1[n].AtmName,"O") != 0) && (strcmp(stack1[n].AtmName,"O1") != 0) &&
							(strcmp(stack1[n].AtmName,"O2") != 0)) && ((strcmp(stack2[o].AtmName,"N") != 0) && (strcmp(stack2[o].AtmName,"H") != 0) && (strcmp(stack2[o].AtmName,"H1") != 0) &&
							(strcmp(stack2[o].AtmName,"H2") != 0) && (strcmp(stack2[o].AtmName,"H3") != 0) && (strcmp(stack2[o].AtmName,"CA") != 0) &&
							(strcmp(stack2[o].AtmName,"C") != 0) && (strcmp(stack2[o].AtmName,"O") != 0) && (strcmp(stack2[o].AtmName,"O1") != 0) && (strcmp(stack2[o].AtmName,"O2") != 0)))
							{
								if (strncmp(stack1[n].AtmName,"H", 1) != 0 && strncmp(stack2[o].AtmName,"H", 1) != 0)
								{
									it[counter] = 3; /*possible electrostatic and hydropathic interaction*/
									ct0[counter] = stack1[n].Res_Type;
									ct1[counter] = stack2[o].Res_Type;
									at0[counter] = stack1[n].Atom_Type;
									at1[counter] = stack2[o].Atom_Type;
									
									//Assigning the ATOMIC SOLVATION PARAMETERS(Cummings, 1995) to the different atom types
									
									if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") == 0) || (strcmp(stack1[n].AtmName,"NZ") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 4;
									}
									else if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") != 0) && (strcmp(stack1[n].AtmName,"NZ") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") == 0) || (strcmp(stack1[n].AtmName,"OE2") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 3;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") != 0) && (strcmp(stack1[n].AtmName,"OE2") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if (strncmp(stack1[n].AtmName,"C", 1) == 0)
									{
										stack1[n].ASP_Type = 0;
									}
									else if (strncmp(stack1[n].AtmName,"S", 1) == 0)
									{
										stack1[n].ASP_Type = 2;
									}
									else
									{
										stack1[n].ASP_Type = 5;
									}
									
									if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") == 0) || (strcmp(stack2[o].AtmName,"NZ") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 4;
									}
									else if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") != 0) && (strcmp(stack2[o].AtmName,"NZ") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") == 0) || (strcmp(stack2[o].AtmName,"OE2") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 3;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") != 0) && (strcmp(stack2[o].AtmName,"OE2") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if (strncmp(stack2[o].AtmName,"C", 1) == 0)
									{
										stack2[o].ASP_Type = 0;
									}
									else if (strncmp(stack2[o].AtmName,"S", 1) == 0)
									{
										stack2[o].ASP_Type = 2;
									}									
									else
									{
										stack2[o].ASP_Type = 5;
									}
									asp0[counter] = stack1[n].ASP_Type;
									asp1[counter] = stack2[o].ASP_Type;
									ro0[counter] = stack1[n].r0;
									ro1[counter] = stack2[o].r0;
									ri0[counter] = stack1[n].NodeIndex;
									ri1[counter] = stack2[o].NodeIndex;
									rd[counter] = sqrt(dist.rsq);	
								}
							}
							counter++;																
						}
						else if (q0[counter] > 0.0 & q1[counter] < 0.0)
						{
							electrostatic++;
							it[counter] = 1; /*possible electrostatic interaction*/
							ct0[counter] = stack1[n].Res_Type;
							ct1[counter] = stack2[o].Res_Type;
							ro0[counter] = stack1[n].r0;
							ro1[counter] = stack2[o].r0;
							ri0[counter] = stack1[n].NodeIndex;
							ri1[counter] = stack2[o].NodeIndex;
							rd[counter] = sqrt(dist.rsq);
							
							//ASSESSING the presence of HYDROGEN BONDS: CHECKING if the DONOR-ACCEPTOR DISTANCE IS <=3.2 Angstrom
							
							if ((strncmp(stack1[n].AtmName,"H",1) == 0) && ((strncmp(stack2[o].AtmName,"O",1) == 0) || (strncmp(stack2[o].AtmName,"N",1) == 0) || (strncmp(stack2[o].AtmName,"S",1) == 0))) 
							{
								if ((strcmp(stack1[n].AtmName,"H2") == 0) || (strcmp(stack1[n].AtmName,"HH12") == 0) || (strcmp(stack1[n].AtmName,"HH22") == 0) || (strcmp(stack1[n].AtmName,"HD22") == 0) || (strcmp(stack1[n].AtmName,"HE22") == 0) || (strcmp(stack1[n].AtmName,"HZ2") == 0))
								{
									if (strncmp(stack1[n-2].AtmName,"N",1) == 0)
									{	
										AtomIndex_1 = stack1[n-2].AtomIndex;
										AtomIndex_2 = stack2[o].AtomIndex;
										if ((strcmp(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Moved[AtomIndex_2-1].AtmHBlabel,"HB") != 0))
										{
											dDA = sqrt((stack2[o].x - stack1[n-2].x)*(stack2[o].x - stack1[n-2].x) + (stack2[o].y -
											stack1[n-2].y)*(stack2[o].y - stack1[n-2].y) + (stack2[o].z - stack1[n-2].z)*(stack2[o].z -
											stack1[n-2].z));
											if (dDA<=3.2)
											{
												it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
												strcpy(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB");
												strcpy(Moved[AtomIndex_2-1].AtmHBlabel,"HB");
											}
										}
									}
								}
								else if ((strcmp(stack1[n].AtmName,"H3") == 0) || (strcmp(stack1[n].AtmName,"HZ3") == 0))
								{	
									if (strncmp(stack1[n-3].AtmName,"N",1) == 0)
									{	
										AtomIndex_1 = stack1[n-3].AtomIndex;
										AtomIndex_2 = stack2[o].AtomIndex;
										if ((strcmp(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Moved[AtomIndex_2-1].AtmHBlabel,"HB") != 0))
										{
											dDA = sqrt((stack2[o].x - stack1[n-3].x)*(stack2[o].x - stack1[n-3].x) + (stack2[o].y -
											stack1[n-3].y)*(stack2[o].y - stack1[n-3].y) + (stack2[o].z - stack1[n-3].z)*(stack2[o].z -
											stack1[n-3].z));
											if (dDA<=3.2)
											{
												it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
												strcpy(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB");
												strcpy(Moved[AtomIndex_2-1].AtmHBlabel,"HB");
											}
										}
									}
								}
								else
								{
									if ((strncmp(stack1[n-1].AtmName,"N",1) == 0) || (strncmp(stack1[n-1].AtmName,"O",1) == 0))
									{
										AtomIndex_1 = stack1[n-1].AtomIndex;
										AtomIndex_2 = stack2[o].AtomIndex;
										if ((strcmp(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Moved[AtomIndex_2-1].AtmHBlabel,"HB") != 0))
										{
											dDA = sqrt((stack2[o].x - stack1[n-1].x)*(stack2[o].x - stack1[n-1].x) + (stack2[o].y -
											stack1[n-1].y)*(stack2[o].y - stack1[n-1].y) + (stack2[o].z - stack1[n-1].z)*(stack2[o].z -
											stack1[n-1].z));
											if (dDA<=3.2)
											{
												it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
												strcpy(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB");
												strcpy(Moved[AtomIndex_2-1].AtmHBlabel,"HB"); 
											}
										}
									}
								}		
							}
							
							//HYDROPHOBIC interactions between non-hydrogen side-chain atoms
							
							if (((strcmp(stack1[n].AtmName,"N") != 0) && (strcmp(stack1[n].AtmName,"H") != 0) && (strcmp(stack1[n].AtmName,"H1") != 0) &&
							(strcmp(stack1[n].AtmName,"H2") != 0) && (strcmp(stack1[n].AtmName,"H3") != 0) && (strcmp(stack1[n].AtmName,"CA") != 0) &&
							(strcmp(stack1[n].AtmName,"C") != 0) && (strcmp(stack1[n].AtmName,"O") != 0) && (strcmp(stack1[n].AtmName,"O1") != 0) &&
							(strcmp(stack1[n].AtmName,"O2") != 0)) && ((strcmp(stack2[o].AtmName,"N") != 0) && (strcmp(stack2[o].AtmName,"H") != 0) && (strcmp(stack2[o].AtmName,"H1") != 0) &&
							(strcmp(stack2[o].AtmName,"H2") != 0) && (strcmp(stack2[o].AtmName,"H3") != 0) && (strcmp(stack2[o].AtmName,"CA") != 0) &&
							(strcmp(stack2[o].AtmName,"C") != 0) && (strcmp(stack2[o].AtmName,"O") != 0) && (strcmp(stack2[o].AtmName,"O1") != 0) && (strcmp(stack2[o].AtmName,"O2") != 0)))
							{
								if (strncmp(stack1[n].AtmName,"H", 1) != 0 && strncmp(stack2[o].AtmName,"H", 1) != 0)
								{
									it[counter] = 3; /*possible electrostatic and hydropathic interaction*/
									ct0[counter] = stack1[n].Res_Type;
									ct1[counter] = stack2[o].Res_Type;
									at0[counter] = stack1[n].Atom_Type;
									at1[counter] = stack2[o].Atom_Type;
									
									//Assigning the ATOMIC SOLVATION PARAMETERS(Cummings, 1995) to the different atom types
									
									if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") == 0) || (strcmp(stack1[n].AtmName,"NZ") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 4;
									}
									else if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") != 0) && (strcmp(stack1[n].AtmName,"NZ") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") == 0) || (strcmp(stack1[n].AtmName,"OE2") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 3;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") != 0) && (strcmp(stack1[n].AtmName,"OE2") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if (strncmp(stack1[n].AtmName,"C", 1) == 0)
									{
										stack1[n].ASP_Type = 0;
									}
									else if (strncmp(stack1[n].AtmName,"S", 1) == 0)
									{
										stack1[n].ASP_Type = 2;
									}
									else
									{
										stack1[n].ASP_Type = 5;
									}
									
									if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") == 0) || (strcmp(stack2[o].AtmName,"NZ") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 4;
									}
									else if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") != 0) && (strcmp(stack2[o].AtmName,"NZ") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") == 0) || (strcmp(stack2[o].AtmName,"OE2") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 3;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") != 0) && (strcmp(stack2[o].AtmName,"OE2") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if (strncmp(stack2[o].AtmName,"C", 1) == 0)
									{
										stack2[o].ASP_Type = 0;
									}
									else if (strncmp(stack2[o].AtmName,"S", 1) == 0)
									{
										stack2[o].ASP_Type = 2;
									}									
									else
									{
										stack2[o].ASP_Type = 5;
									}
									asp0[counter] = stack1[n].ASP_Type;
									asp1[counter] = stack2[o].ASP_Type;
									ro0[counter] = stack1[n].r0;
									ro1[counter] = stack2[o].r0;
									ri0[counter] = stack1[n].NodeIndex;
									ri1[counter] = stack2[o].NodeIndex;
									rd[counter] = sqrt(dist.rsq);	
								}
							}
							counter++;																
						}
						
						//HYDROPHOBIC interactions between non-hydrogen side-chain atoms
																										
						else if (((strcmp(stack1[n].AtmName,"N") != 0) && (strcmp(stack1[n].AtmName,"H") != 0) && (strcmp(stack1[n].AtmName,"H1") != 0) &&
						(strcmp(stack1[n].AtmName,"H2") != 0) && (strcmp(stack1[n].AtmName,"H3") != 0) && (strcmp(stack1[n].AtmName,"CA") != 0) &&
						(strcmp(stack1[n].AtmName,"C") != 0) && (strcmp(stack1[n].AtmName,"O") != 0) && (strcmp(stack1[n].AtmName,"O1") != 0) &&
						(strcmp(stack1[n].AtmName,"O2") != 0)) && ((strcmp(stack2[o].AtmName,"N") != 0) && (strcmp(stack2[o].AtmName,"H") != 0) && (strcmp(stack2[o].AtmName,"H1") != 0) &&
						(strcmp(stack2[o].AtmName,"H2") != 0) && (strcmp(stack2[o].AtmName,"H3") != 0) && (strcmp(stack2[o].AtmName,"CA") != 0) &&
						(strcmp(stack2[o].AtmName,"C") != 0) && (strcmp(stack2[o].AtmName,"O") != 0) && (strcmp(stack2[o].AtmName,"O1") != 0) && (strcmp(stack2[o].AtmName,"O2") != 0)))
						{
							if (strncmp(stack1[n].AtmName,"H", 1) != 0 && strncmp(stack2[o].AtmName,"H", 1) != 0)
							{
								it[counter] = 2; /*possible hydropathic interaction*/
								ct0[counter] = stack1[n].Res_Type;
								ct1[counter] = stack2[o].Res_Type;
								at0[counter] = stack1[n].Atom_Type;
								at1[counter] = stack2[o].Atom_Type;
								
								//Assigning the ATOMIC SOLVATION PARAMETERS(Cummings, 1995) to the different atom types
								
								if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") == 0) || (strcmp(stack1[n].AtmName,"NZ") == 0))) //fully charged group
								{
									stack1[n].ASP_Type = 4;
								}
								else if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") != 0) && (strcmp(stack1[n].AtmName,"NZ") != 0)))
								{
									stack1[n].ASP_Type = 1;
								}
								else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") == 0) || (strcmp(stack1[n].AtmName,"OE2") == 0))) //fully charged group
								{
									stack1[n].ASP_Type = 3;
								}
								else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") != 0) && (strcmp(stack1[n].AtmName,"OE2") != 0)))
								{
									stack1[n].ASP_Type = 1;
								}
								else if (strncmp(stack1[n].AtmName,"C", 1) == 0)
								{
									stack1[n].ASP_Type = 0;
								}
								else if (strncmp(stack1[n].AtmName,"S", 1) == 0)
								{
									stack1[n].ASP_Type = 2;
								}
								else
								{
									stack1[n].ASP_Type = 5;
								}
									
								if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") == 0) || (strcmp(stack2[o].AtmName,"NZ") == 0))) //fully charged group
								{
									stack2[o].ASP_Type = 4;
								}
								else if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") != 0) && (strcmp(stack2[o].AtmName,"NZ") != 0)))
								{
									stack2[o].ASP_Type = 1;
								}
								else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") == 0) || (strcmp(stack2[o].AtmName,"OE2") == 0))) //fully charged group
								{
									stack2[o].ASP_Type = 3;
								}
								else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") != 0) && (strcmp(stack2[o].AtmName,"OE2") != 0)))
								{
									stack2[o].ASP_Type = 1;
								}
								else if (strncmp(stack2[o].AtmName,"C", 1) == 0)
								{
									stack2[o].ASP_Type = 0;
								}
								else if (strncmp(stack2[o].AtmName,"S", 1) == 0)
								{
									stack2[o].ASP_Type = 2;
								}									
								else
								{
									stack2[o].ASP_Type = 5;
								}
								asp0[counter] = stack1[n].ASP_Type;
								asp1[counter] = stack2[o].ASP_Type;
								ro0[counter] = stack1[n].r0;
								ro1[counter] = stack2[o].r0;
								ri0[counter] = stack1[n].NodeIndex;
								ri1[counter] = stack2[o].NodeIndex;
								rd[counter] = sqrt(dist.rsq);
								counter++;
							}
						}
						/*inter_residue_potential[counter] = hydropathic_inter_residue_potentials[stack1[n].Res_Type][stack2[o].Res_Type];*/
						/*if (counter == MAXCONTACT - 1)
						{
							break;
						}
						/*clashes_cutoff[counter] = dist.r01;
						contacts_cutoff[counter] = dist.r02;*/
					}
				}
			
			}
		}
							
	}
	
	/*for (m = 0; m < counter; m++)
	{
		for (k = 0; k <= nAtom[0]; k++)
		{
			n = residues_index[counter][0];
			if (n == Atoms[k][0].NodeIndex & strcmp(Atoms[k][0].AtmName,"CB") == 0)
			{
				x[counter][0] = Atoms[k][0].x;
				y[counter][0] = Atoms[k][0].y;
				z[counter][0] = Atoms[k][0].z;
				ro[counter][0] = Atoms[k][0].r0;
			}
		}
		for (l = 0; l <= nAtom[1]; l++)
		{
			o = residues_index[counter][1];	
			if (o == Moved[l].NodeIndex & strcmp(Moved[l].AtmName,"CB") == 0)
			{
				x[counter][1] = Moved[l].x;
				y[counter][1] = Moved[l].y;
				z[counter][1] = Moved[l].z;
				ro[counter][1] = Moved[l].r0;
			}
		}
		dist.rsq = (x[counter][0] - x[counter][1])*(x[counter][0] - x[counter][1]) + (y[counter][0] - y[counter][1])*(y[counter][0] - y[counter][1]) + (z[counter][0] - z[counter][1])*(z[counter][0] - z[counter][1]);
		residues_distance[counter] = sqrt(dist.rsq);
	}

	/*for (counter = 0; counter < MAXCONTACT; counter++)
	{
		sum_inter_residue_potential = sum_inter_residue_potential + inter_residue_potential[counter];
	}*/
	
	return counter;
	fclose(test);
}


int CONTACTSMOVED(void)
{
	int n,o,counter;
	
	counter=0;
	for (n=0;n<nAtom[0];n++)
	{
		for (o=0;o<nAtom[1];o++)
		{
			dist.r01 = HARDCORE_R*(Atoms[n][0].r0+Moved[o].r0);
			dist.r02 = CONTACT_R*(Atoms[n][0].r0+Moved[o].r0);
			dist.rsq = (Moved[o].x-Atoms[n][0].x)*(Moved[o].x-Atoms[n][0].x)+(Moved[o].y-Atoms[n][0].y)*(Moved[o].y-Atoms[n][0].y)+(Moved[o].z-Atoms[n][0].z)*(Moved[o].z-Atoms[n][0].z);
				 				
			if ((dist.rsq>=dist.r01*dist.r01)&&(dist.rsq<=dist.r02*dist.r02))
			{
				counter++;
			}
			
		}
							
	}
	return counter;
}


int CLASHESMOVED(void)
{
	int n,o,counter;
	
	counter=0;
	for (n=0;n<nAtom[0];n++)
	{
		for (o=0;o<nAtom[1];o++)
		{
			dist.r01 = HARDCORE_R*(Atoms[n][0].r0+Moved[o].r0);
			dist.r02 = CONTACT_R*(Atoms[n][0].r0+Moved[o].r0);
			dist.rsq = (Moved[o].x-Atoms[n][0].x)*(Moved[o].x-Atoms[n][0].x)+(Moved[o].y-Atoms[n][0].y)*(Moved[o].y-Atoms[n][0].y)+(Moved[o].z-Atoms[n][0].z)*(Moved[o].z-Atoms[n][0].z);
				 				
			if (dist.rsq<(dist.r01*dist.r01))
			{
				counter++;
			}
											
		}
							
	}
	return counter;
}


int TYPECONTACTS(int ct0[MAXCONTACT], int ct1[MAXCONTACT], int at0[MAXCONTACT], int at1[MAXCONTACT], int ri0[MAXCONTACT], int ri1[MAXCONTACT], double ro0[MAXCONTACT], double ro1[MAXCONTACT], double rd[MAXCONTACT], int it[MAXCONTACT], float q0[MAXCONTACT], float q1[MAXCONTACT], int asp0[MAXCONTACT], int asp1[MAXCONTACT])
{
	FILE *test;
	
	int i, k, l, n, o, p, r, maxcontactflag, counter, len1, len2, residue_index_0, residue_index_1, s, t, ss, tt;
	int AtomIndex_1, AtomIndex_2;
	int electrostatic = 0;
	float dDA, dc1, dc2, dc3;
	/*float sum_inter_residue_potential = 0.0;
	float inter_residue_potential[MAXCONTACT];*/
	
	/*for (i = 0; i < MAXCONTACT; i++)
	{
		inter_residue_potential[i] = 0.0;
	}*/
	
	for (k = 0; k<nAtom[0]; k++)
	{
		strcpy(Atoms[k][0].AtmHBlabel, "NoHB");
	}
	
	for (l = 0; l<nAtom[1]; l++)
	{
		strcpy(Atoms[l][1].AtmHBlabel, "NoHB");
	}
	
	
	counter = 0;
	s = 0;
	t = 0;
	for (i = Atoms[0][0].NodeIndex; i <= Atoms[nAtom[0] - 1][0].NodeIndex; i++)
	{
		//picking coordinates of residue with index i.
		
		/*if (contactflag)
		{
			break;
		}*/

		len1 = 0;
		s = 0;
		for (k = 0; k<nAtom[0]; k++)
		{
			if (Atoms[k][0].NodeIndex == i)
			{
				stack1[len1].x = Atoms[k][0].x;
				stack1[len1].y = Atoms[k][0].y;
				stack1[len1].z = Atoms[k][0].z;
				stack1[len1].r0 = Atoms[k][0].r0;
				stack1[len1].Res_Type = Atoms[k][0].Res_Type;
				stack1[len1].Atom_Type = Atoms[k][0].Atom_Type;
				stack1[len1].NodeIndex = Atoms[k][0].NodeIndex;
				stack1[len1].AtomIndex = Atoms[k][0].AtomIndex;
				strcpy(stack1[len1].AtmName, Atoms[k][0].AtmName);
				len1++;	
			}
		}
		for (k = Atoms[0][1].NodeIndex; k <= Atoms[nAtom[1] - 1][1].NodeIndex; k++)
		{
			//picking coordinates of residue with index k 
			
			/*if (maxcontactflag == 1)
			{
				break;
			}*/

			len2 = 0;
			t = 0;
			for (l = 0; l<nAtom[1]; l++)
			{
				if (Atoms[l][1].NodeIndex == k)
				{
					stack2[len2].x = Atoms[l][1].x;
					stack2[len2].y = Atoms[l][1].y;
					stack2[len2].z = Atoms[l][1].z;
					stack2[len2].r0 = Atoms[l][1].r0;
					stack2[len2].Res_Type = Atoms[l][1].Res_Type;
					stack2[len2].Atom_Type = Atoms[l][1].Atom_Type;
					stack2[len2].NodeIndex = Atoms[l][1].NodeIndex;
					stack2[len2].AtomIndex = Atoms[l][1].AtomIndex;
					strcpy(stack2[len2].AtmName, Atoms[l][1].AtmName);
					len2++;		
				}
			}		
			for (n=0;n<len1;n++)
			{
				/*if (contactflag)
				{
					break;
				}*/
				for (o=0;o<len2;o++)
				{
					dist.r01 = HARDCORE_R*(stack1[n].r0 + stack2[o].r0);
					dist.r02 = CONTACT_R*(stack1[n].r0 + stack2[o].r0);
					dist.rsq = (stack2[o].x - stack1[n].x)*(stack2[o].x - stack1[n].x) + (stack2[o].y - stack1[n].y)*(stack2[o].y - stack1[n].y) + (stack2[o].z - stack1[n].z)*(stack2[o].z - stack1[n].z);
					
					if ((dist.rsq>=dist.r01*dist.r01)&&(dist.rsq<=dist.r02*dist.r02))
					{
						for (p = 0; p < LENCHARGESFILE; p++)
						{
							if ((strcmp(stack1[n].AtmName,Charges[p].AtmName) == 0) && (stack1[n].Res_Type == Charges[p].Res_Type))
							{
								for (r = 0; r < LENCHARGESFILE; r++)
								{
									if ((strcmp(stack2[o].AtmName,Charges[r].AtmName) == 0) && (stack2[o].Res_Type == Charges[r].Res_Type))
									{
										q0[counter] = Charges[p].Atm_Charge; //GROMOS54A7 point charges of each interacting atom
										q1[counter] = Charges[r].Atm_Charge; //GROMOS54A7 point charges of each interacting atom
									}
								}	
							}							
						}
						if (q0[counter] < 0.0 & q1[counter] < 0.0)
						{	
							electrostatic++;
							it[counter] = 1; /*possible electrostatic interaction*/
							ct0[counter] = stack1[n].Res_Type;
							ct1[counter] = stack2[o].Res_Type;
							ro0[counter] = stack1[n].r0;
							ro1[counter] = stack2[o].r0;
							ri0[counter] = stack1[n].NodeIndex;
							ri1[counter] = stack2[o].NodeIndex;
							rd[counter] = sqrt(dist.rsq);
							
							//HYDROPHOBIC interactions between non-hydrogen side-chain atoms
							
							if (((strcmp(stack1[n].AtmName,"N") != 0) && (strcmp(stack1[n].AtmName,"H") != 0) && (strcmp(stack1[n].AtmName,"H1") != 0) &&
							(strcmp(stack1[n].AtmName,"H2") != 0) && (strcmp(stack1[n].AtmName,"H3") != 0) && (strcmp(stack1[n].AtmName,"CA") != 0) &&
							(strcmp(stack1[n].AtmName,"C") != 0) && (strcmp(stack1[n].AtmName,"O") != 0) && (strcmp(stack1[n].AtmName,"O1") != 0) &&
							(strcmp(stack1[n].AtmName,"O2") != 0)) && ((strcmp(stack2[o].AtmName,"N") != 0) && (strcmp(stack2[o].AtmName,"H") != 0) && (strcmp(stack2[o].AtmName,"H1") != 0) &&
							(strcmp(stack2[o].AtmName,"H2") != 0) && (strcmp(stack2[o].AtmName,"H3") != 0) && (strcmp(stack2[o].AtmName,"CA") != 0) &&
							(strcmp(stack2[o].AtmName,"C") != 0) && (strcmp(stack2[o].AtmName,"O") != 0) && (strcmp(stack2[o].AtmName,"O1") != 0) && (strcmp(stack2[o].AtmName,"O2") != 0)))
							{
								if (strncmp(stack1[n].AtmName,"H", 1) != 0 && strncmp(stack2[o].AtmName,"H", 1) != 0)
								{
									it[counter] = 3; /*possible electrostatic and hydropathic interaction*/
									ct0[counter] = stack1[n].Res_Type;
									ct1[counter] = stack2[o].Res_Type;
									at0[counter] = stack1[n].Atom_Type;
									at1[counter] = stack2[o].Atom_Type;
									
									//Assigning the ATOMIC SOLVATION PARAMETERS(Cummings, 1995) to the different atom types
									
									if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") == 0) || (strcmp(stack1[n].AtmName,"NZ") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 4;
									}
									else if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") != 0) && (strcmp(stack1[n].AtmName,"NZ") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") == 0) || (strcmp(stack1[n].AtmName,"OE2") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 3;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") != 0) && (strcmp(stack1[n].AtmName,"OE2") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if (strncmp(stack1[n].AtmName,"C", 1) == 0)
									{
										stack1[n].ASP_Type = 0;
									}
									else if (strncmp(stack1[n].AtmName,"S", 1) == 0)
									{
										stack1[n].ASP_Type = 2;
									}
									else
									{
										stack1[n].ASP_Type = 5;
									}
									
									if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") == 0) || (strcmp(stack2[o].AtmName,"NZ") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 4;
									}
									else if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") != 0) && (strcmp(stack2[o].AtmName,"NZ") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") == 0) || (strcmp(stack2[o].AtmName,"OE2") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 3;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") != 0) && (strcmp(stack2[o].AtmName,"OE2") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if (strncmp(stack2[o].AtmName,"C", 1) == 0)
									{
										stack2[o].ASP_Type = 0;
									}
									else if (strncmp(stack2[o].AtmName,"S", 1) == 0)
									{
										stack2[o].ASP_Type = 2;
									}									
									else
									{
										stack2[o].ASP_Type = 5;
									}
									asp0[counter] = stack1[n].ASP_Type;
									asp1[counter] = stack2[o].ASP_Type;
									ro0[counter] = stack1[n].r0;
									ro1[counter] = stack2[o].r0;
									ri0[counter] = stack1[n].NodeIndex;
									ri1[counter] = stack2[o].NodeIndex;
									rd[counter] = sqrt(dist.rsq);
								}
							}
							counter++;																
						}
						else if (q0[counter] > 0.0 & q1[counter] > 0.0)
						{
							electrostatic++;
							it[counter] = 1; /*possible electrostatic interaction*/
							ct0[counter] = stack1[n].Res_Type;
							ct1[counter] = stack2[o].Res_Type;
							ro0[counter] = stack1[n].r0;
							ro1[counter] = stack2[o].r0;
							ri0[counter] = stack1[n].NodeIndex;
							ri1[counter] = stack2[o].NodeIndex;
							rd[counter] = sqrt(dist.rsq);
							
							//HYDROPHOBIC interactions between non-hydrogen side-chain atoms
							
							if (((strcmp(stack1[n].AtmName,"N") != 0) && (strcmp(stack1[n].AtmName,"H") != 0) && (strcmp(stack1[n].AtmName,"H1") != 0) &&
							(strcmp(stack1[n].AtmName,"H2") != 0) && (strcmp(stack1[n].AtmName,"H3") != 0) && (strcmp(stack1[n].AtmName,"CA") != 0) &&
							(strcmp(stack1[n].AtmName,"C") != 0) && (strcmp(stack1[n].AtmName,"O") != 0) && (strcmp(stack1[n].AtmName,"O1") != 0) &&
							(strcmp(stack1[n].AtmName,"O2") != 0)) && ((strcmp(stack2[o].AtmName,"N") != 0) && (strcmp(stack2[o].AtmName,"H") != 0) && (strcmp(stack2[o].AtmName,"H1") != 0) &&
							(strcmp(stack2[o].AtmName,"H2") != 0) && (strcmp(stack2[o].AtmName,"H3") != 0) && (strcmp(stack2[o].AtmName,"CA") != 0) &&
							(strcmp(stack2[o].AtmName,"C") != 0) && (strcmp(stack2[o].AtmName,"O") != 0) && (strcmp(stack2[o].AtmName,"O1") != 0) && (strcmp(stack2[o].AtmName,"O2") != 0)))
							{
								if (strncmp(stack1[n].AtmName,"H", 1) != 0 && strncmp(stack2[o].AtmName,"H", 1) != 0)
								{
									it[counter] = 3; /*possible electrostatic and hydropathic interaction*/
									ct0[counter] = stack1[n].Res_Type;
									ct1[counter] = stack2[o].Res_Type;
									at0[counter] = stack1[n].Atom_Type;
									at1[counter] = stack2[o].Atom_Type;
									
									//Assigning the ATOMIC SOLVATION PARAMETERS(Cummings, 1995) to the different atom types
									
									if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") == 0) || (strcmp(stack1[n].AtmName,"NZ") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 4;
									}
									else if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") != 0) && (strcmp(stack1[n].AtmName,"NZ") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") == 0) || (strcmp(stack1[n].AtmName,"OE2") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 3;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") != 0) && (strcmp(stack1[n].AtmName,"OE2") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if (strncmp(stack1[n].AtmName,"C", 1) == 0)
									{
										stack1[n].ASP_Type = 0;
									}
									else if (strncmp(stack1[n].AtmName,"S", 1) == 0)
									{
										stack1[n].ASP_Type = 2;
									}
									else
									{
										stack1[n].ASP_Type = 5;
									}
									
									if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") == 0) || (strcmp(stack2[o].AtmName,"NZ") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 4;
									}
									else if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") != 0) && (strcmp(stack2[o].AtmName,"NZ") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") == 0) || (strcmp(stack2[o].AtmName,"OE2") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 3;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") != 0) && (strcmp(stack2[o].AtmName,"OE2") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if (strncmp(stack2[o].AtmName,"C", 1) == 0)
									{
										stack2[o].ASP_Type = 0;
									}
									else if (strncmp(stack2[o].AtmName,"S", 1) == 0)
									{
										stack2[o].ASP_Type = 2;
									}									
									else
									{
										stack2[o].ASP_Type = 5;
									}
									asp0[counter] = stack1[n].ASP_Type;
									asp1[counter] = stack2[o].ASP_Type;
									ro0[counter] = stack1[n].r0;
									ro1[counter] = stack2[o].r0;
									ri0[counter] = stack1[n].NodeIndex;
									ri1[counter] = stack2[o].NodeIndex;
									rd[counter] = sqrt(dist.rsq);
								}
							}
							counter++;																
						}
						else if (q0[counter] < 0.0 & q1[counter] > 0.0)
						{
							electrostatic++;
							it[counter] = 1; /*possible electrostatic interaction*/
							ct0[counter] = stack1[n].Res_Type;
							ct1[counter] = stack2[o].Res_Type;
							ro0[counter] = stack1[n].r0;
							ro1[counter] = stack2[o].r0;
							ri0[counter] = stack1[n].NodeIndex;
							ri1[counter] = stack2[o].NodeIndex;
							rd[counter] = sqrt(dist.rsq);
							
							//ASSESSING the presence of HYDROGEN BONDS: CHECKING if the DONOR-ACCEPTOR DISTANCE IS <=3.2 Angstrom
							if (((strncmp(stack1[n].AtmName,"O",1) == 0) || (strncmp(stack1[n].AtmName,"N",1) == 0) || (strncmp(stack1[n].AtmName,"S",1) == 0)) && (strncmp(stack2[o].AtmName,"H",1) == 0))
							{
								if ((strcmp(stack2[o].AtmName,"H2") == 0) || (strcmp(stack2[o].AtmName,"HH12") == 0) || (strcmp(stack2[o].AtmName,"HH22") == 0) || (strcmp(stack2[o].AtmName,"HD22") == 0) || (strcmp(stack2[o].AtmName,"HE22") == 0) || (strcmp(stack2[o].AtmName,"HZ2") == 0))
								{
									if (strncmp(stack2[o-2].AtmName,"N",1) == 0)
									{	
										AtomIndex_1 = stack1[n].AtomIndex;
										AtomIndex_2 = stack2[o-2].AtomIndex;
										if ((strcmp(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Atoms[AtomIndex_2-1][1].AtmHBlabel,"HB") != 0))
										{
											dDA = sqrt((stack2[o-2].x - stack1[n].x)*(stack2[o-2].x - stack1[n].x) + (stack2[o-2].y -
											stack1[n].y)*(stack2[o-2].y - stack1[n].y) + (stack2[o-2].z - stack1[n].z)*(stack2[o-2].z -
											stack1[n].z));					
											if (dDA<=3.2)
											{
												it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
												strcpy(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB");
												strcpy(Atoms[AtomIndex_2-1][1].AtmHBlabel,"HB");
											}
										}	
									}
								}
								else if ((strcmp(stack2[o].AtmName,"H3") == 0) || (strcmp(stack2[o].AtmName,"HZ3") == 0))
								{
									if (strncmp(stack2[o-3].AtmName,"N",1) == 0)
									{
										AtomIndex_1 = stack1[n].AtomIndex;
										AtomIndex_2 = stack2[o-3].AtomIndex;
										if ((strcmp(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Atoms[AtomIndex_2-1][1].AtmHBlabel,"HB") != 0))
										{
											dDA = sqrt((stack2[o-3].x - stack1[n].x)*(stack2[o-3].x - stack1[n].x) + (stack2[o-3].y -
											stack1[n].y)*(stack2[o-3].y - stack1[n].y) + (stack2[o-3].z - stack1[n].z)*(stack2[o-3].z -
											stack1[n].z));					
											if (dDA<=3.2)
											{
												it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
												strcpy(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB");
												strcpy(Atoms[AtomIndex_2-1][1].AtmHBlabel,"HB");
											}
										}
									}	
								}
								else
								{
									if ((strncmp(stack2[o-1].AtmName,"N",1) == 0) || (strncmp(stack2[o-1].AtmName,"O",1) == 0))
									{
										AtomIndex_1 = stack1[n].AtomIndex;
										AtomIndex_2 = stack2[o-1].AtomIndex;
										if ((strcmp(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Atoms[AtomIndex_2-1][1].AtmHBlabel,"HB") != 0))
										{
											dDA = sqrt((stack2[o-1].x - stack1[n].x)*(stack2[o-1].x - stack1[n].x) + (stack2[o-1].y -
											stack1[n].y)*(stack2[o-1].y - stack1[n].y) + (stack2[o-1].z - stack1[n].z)*(stack2[o-1].z -
											stack1[n].z));					
											if (dDA<=3.2)
											{
												it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
												strcpy(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB");
												strcpy(Atoms[AtomIndex_2-1][1].AtmHBlabel,"HB"); 
											}
										}
									}
								}									
							}
							
							//HYDROPHOBIC interactions between non-hydrogen side-chain atoms
							
							if (((strcmp(stack1[n].AtmName,"N") != 0) && (strcmp(stack1[n].AtmName,"H") != 0) && (strcmp(stack1[n].AtmName,"H1") != 0) &&
							(strcmp(stack1[n].AtmName,"H2") != 0) && (strcmp(stack1[n].AtmName,"H3") != 0) && (strcmp(stack1[n].AtmName,"CA") != 0) &&
							(strcmp(stack1[n].AtmName,"C") != 0) && (strcmp(stack1[n].AtmName,"O") != 0) && (strcmp(stack1[n].AtmName,"O1") != 0) &&
							(strcmp(stack1[n].AtmName,"O2") != 0)) && ((strcmp(stack2[o].AtmName,"N") != 0) && (strcmp(stack2[o].AtmName,"H") != 0) && (strcmp(stack2[o].AtmName,"H1") != 0) &&
							(strcmp(stack2[o].AtmName,"H2") != 0) && (strcmp(stack2[o].AtmName,"H3") != 0) && (strcmp(stack2[o].AtmName,"CA") != 0) &&
							(strcmp(stack2[o].AtmName,"C") != 0) && (strcmp(stack2[o].AtmName,"O") != 0) && (strcmp(stack2[o].AtmName,"O1") != 0) && (strcmp(stack2[o].AtmName,"O2") != 0)))
							{
								if (strncmp(stack1[n].AtmName,"H", 1) != 0 && strncmp(stack2[o].AtmName,"H", 1) != 0)
								{
									it[counter] = 3; /*possible electrostatic and hydropathic interaction*/
									ct0[counter] = stack1[n].Res_Type;
									ct1[counter] = stack2[o].Res_Type;
									at0[counter] = stack1[n].Atom_Type;
									at1[counter] = stack2[o].Atom_Type;
									
									//Assigning the ATOMIC SOLVATION PARAMETERS(Cummings, 1995) to the different atom types
									
									if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") == 0) || (strcmp(stack1[n].AtmName,"NZ") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 4;
									}
									else if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") != 0) && (strcmp(stack1[n].AtmName,"NZ") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") == 0) || (strcmp(stack1[n].AtmName,"OE2") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 3;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") != 0) && (strcmp(stack1[n].AtmName,"OE2") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if (strncmp(stack1[n].AtmName,"C", 1) == 0)
									{
										stack1[n].ASP_Type = 0;
									}
									else if (strncmp(stack1[n].AtmName,"S", 1) == 0)
									{
										stack1[n].ASP_Type = 2;
									}
									else
									{
										stack1[n].ASP_Type = 5;
									}
									
									if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") == 0) || (strcmp(stack2[o].AtmName,"NZ") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 4;
									}
									else if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") != 0) && (strcmp(stack2[o].AtmName,"NZ") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") == 0) || (strcmp(stack2[o].AtmName,"OE2") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 3;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") != 0) && (strcmp(stack2[o].AtmName,"OE2") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if (strncmp(stack2[o].AtmName,"C", 1) == 0)
									{
										stack2[o].ASP_Type = 0;
									}
									else if (strncmp(stack2[o].AtmName,"S", 1) == 0)
									{
										stack2[o].ASP_Type = 2;
									}									
									else
									{
										stack2[o].ASP_Type = 5;
									}
									asp0[counter] = stack1[n].ASP_Type;
									asp1[counter] = stack2[o].ASP_Type;
									ro0[counter] = stack1[n].r0;
									ro1[counter] = stack2[o].r0;
									ri0[counter] = stack1[n].NodeIndex;
									ri1[counter] = stack2[o].NodeIndex;
									rd[counter] = sqrt(dist.rsq);
								}
							}
							counter++;																
						}
						else if (q0[counter] > 0.0 & q1[counter] < 0.0)
						{
							electrostatic++;
							it[counter] = 1; /*possible electrostatic interaction*/
							ct0[counter] = stack1[n].Res_Type;
							ct1[counter] = stack2[o].Res_Type;
							at0[counter] = stack1[n].Atom_Type;
							at1[counter] = stack2[o].Atom_Type;
							ro0[counter] = stack1[n].r0;
							ro1[counter] = stack2[o].r0;
							ri0[counter] = stack1[n].NodeIndex;
							ri1[counter] = stack2[o].NodeIndex;
							rd[counter] = sqrt(dist.rsq);
							
							//ASSESSING the presence of HYDROGEN BONDS: CHECKING if the DONOR-ACCEPTOR DISTANCE IS <=3.2 Angstrom
							
							if ((strncmp(stack1[n].AtmName,"H",1) == 0) && ((strncmp(stack2[o].AtmName,"O",1) == 0) || (strncmp(stack2[o].AtmName,"N",1) == 0) || (strncmp(stack2[o].AtmName,"S",1) == 0))) 
							{
								if ((strcmp(stack1[n].AtmName,"H2") == 0) || (strcmp(stack1[n].AtmName,"HH12") == 0) || (strcmp(stack1[n].AtmName,"HH22") == 0) || (strcmp(stack1[n].AtmName,"HD22") == 0) || (strcmp(stack1[n].AtmName,"HE22") == 0) || (strcmp(stack1[n].AtmName,"HZ2") == 0))
								{
									if (strncmp(stack1[n-2].AtmName,"N",1) == 0)
									{	
										AtomIndex_1 = stack1[n-2].AtomIndex;
										AtomIndex_2 = stack2[o].AtomIndex;
										if ((strcmp(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Atoms[AtomIndex_2-1][1].AtmHBlabel,"HB") != 0))
										{
											dDA = sqrt((stack2[o].x - stack1[n-2].x)*(stack2[o].x - stack1[n-2].x) + (stack2[o].y -
											stack1[n-2].y)*(stack2[o].y - stack1[n-2].y) + (stack2[o].z - stack1[n-2].z)*(stack2[o].z -
											stack1[n-2].z));
											if (dDA<=3.2)
											{
												it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
												strcpy(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB");
												strcpy(Atoms[AtomIndex_2-1][1].AtmHBlabel,"HB"); 
											}
										}
									}
								}
								else if ((strcmp(stack1[n].AtmName,"H3") == 0) || (strcmp(stack1[n].AtmName,"HZ3") == 0))
								{	
									if (strncmp(stack1[n-3].AtmName,"N",1) == 0)
									{	
										AtomIndex_1 = stack1[n-3].AtomIndex;
										AtomIndex_2 = stack2[o].AtomIndex;
										if ((strcmp(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Atoms[AtomIndex_2-1][1].AtmHBlabel,"HB") != 0))
										{
											dDA = sqrt((stack2[o].x - stack1[n-3].x)*(stack2[o].x - stack1[n-3].x) + (stack2[o].y -
											stack1[n-3].y)*(stack2[o].y - stack1[n-3].y) + (stack2[o].z - stack1[n-3].z)*(stack2[o].z -
											stack1[n-3].z));
											if (dDA<=3.2)
											{
												it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
												strcpy(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB");
												strcpy(Atoms[AtomIndex_2-1][1].AtmHBlabel,"HB");
											}
										}
									}
								}
								else
								{
									if ((strncmp(stack1[n-1].AtmName,"N",1) == 0) || (strncmp(stack1[n-1].AtmName,"O",1) == 0))
									{
										AtomIndex_1 = stack1[n-1].AtomIndex;
										AtomIndex_2 = stack2[o].AtomIndex;
										if ((strcmp(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Atoms[AtomIndex_2-1][1].AtmHBlabel,"HB") != 0))
										{
											dDA = sqrt((stack2[o].x - stack1[n-1].x)*(stack2[o].x - stack1[n-1].x) + (stack2[o].y -
											stack1[n-1].y)*(stack2[o].y - stack1[n-1].y) + (stack2[o].z - stack1[n-1].z)*(stack2[o].z -
											stack1[n-1].z));
											if (dDA<=3.2)
											{
												it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
												strcpy(Atoms[AtomIndex_1-1][0].AtmHBlabel,"HB");
												strcpy(Atoms[AtomIndex_2-1][1].AtmHBlabel,"HB");
											}
										}
									}
								}		
							}
							
							//HYDROPHOBIC interactions between non-hydrogen side-chain atoms
							
							if (((strcmp(stack1[n].AtmName,"N") != 0) && (strcmp(stack1[n].AtmName,"H") != 0) && (strcmp(stack1[n].AtmName,"H1") != 0) &&
							(strcmp(stack1[n].AtmName,"H2") != 0) && (strcmp(stack1[n].AtmName,"H3") != 0) && (strcmp(stack1[n].AtmName,"CA") != 0) &&
							(strcmp(stack1[n].AtmName,"C") != 0) && (strcmp(stack1[n].AtmName,"O") != 0) && (strcmp(stack1[n].AtmName,"O1") != 0) &&
							(strcmp(stack1[n].AtmName,"O2") != 0)) && ((strcmp(stack2[o].AtmName,"N") != 0) && (strcmp(stack2[o].AtmName,"H") != 0) && (strcmp(stack2[o].AtmName,"H1") != 0) &&
							(strcmp(stack2[o].AtmName,"H2") != 0) && (strcmp(stack2[o].AtmName,"H3") != 0) && (strcmp(stack2[o].AtmName,"CA") != 0) &&
							(strcmp(stack2[o].AtmName,"C") != 0) && (strcmp(stack2[o].AtmName,"O") != 0) && (strcmp(stack2[o].AtmName,"O1") != 0) && (strcmp(stack2[o].AtmName,"O2") != 0)))
							{
								if (strncmp(stack1[n].AtmName,"H", 1) != 0 && strncmp(stack2[o].AtmName,"H", 1) != 0)
								{
									it[counter] = 3; /*possible electrostatic and hydropathic interaction*/
									ct0[counter] = stack1[n].Res_Type;
									ct1[counter] = stack2[o].Res_Type;
									at0[counter] = stack1[n].Atom_Type;
									at1[counter] = stack2[o].Atom_Type;
									
									//Assigning the ATOMIC SOLVATION PARAMETERS(Cummings, 1995) to the different atom types
									
									if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") == 0) || (strcmp(stack1[n].AtmName,"NZ") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 4;
									}
									else if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") != 0) && (strcmp(stack1[n].AtmName,"NZ") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") == 0) || (strcmp(stack1[n].AtmName,"OE2") == 0))) //fully charged group
									{
										stack1[n].ASP_Type = 3;
									}
									else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") != 0) && (strcmp(stack1[n].AtmName,"OE2") != 0)))
									{
										stack1[n].ASP_Type = 1;
									}
									else if (strncmp(stack1[n].AtmName,"C", 1) == 0)
									{
										stack1[n].ASP_Type = 0;
									}
									else if (strncmp(stack1[n].AtmName,"S", 1) == 0)
									{
										stack1[n].ASP_Type = 2;
									}
									else
									{
										stack1[n].ASP_Type = 5;
									}
									
									if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") == 0) || (strcmp(stack2[o].AtmName,"NZ") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 4;
									}
									else if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") != 0) && (strcmp(stack2[o].AtmName,"NZ") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") == 0) || (strcmp(stack2[o].AtmName,"OE2") == 0))) //fully charged group
									{
										stack2[o].ASP_Type = 3;
									}
									else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") != 0) && (strcmp(stack2[o].AtmName,"OE2") != 0)))
									{
										stack2[o].ASP_Type = 1;
									}
									else if (strncmp(stack2[o].AtmName,"C", 1) == 0)
									{
										stack2[o].ASP_Type = 0;
									}
									else if (strncmp(stack2[o].AtmName,"S", 1) == 0)
									{
										stack2[o].ASP_Type = 2;
									}									
									else
									{
										stack2[o].ASP_Type = 5;
									}
									asp0[counter] = stack1[n].ASP_Type;
									asp1[counter] = stack2[o].ASP_Type;
									ro0[counter] = stack1[n].r0;
									ro1[counter] = stack2[o].r0;
									ri0[counter] = stack1[n].NodeIndex;
									ri1[counter] = stack2[o].NodeIndex;
									rd[counter] = sqrt(dist.rsq);
								}
							}
							counter++;																
						}
						
						//HYDROPHOBIC interactions between non-hydrogen side-chain atoms
																										
						else if (((strcmp(stack1[n].AtmName,"N") != 0) && (strcmp(stack1[n].AtmName,"H") != 0) && (strcmp(stack1[n].AtmName,"H1") != 0) &&
						(strcmp(stack1[n].AtmName,"H2") != 0) && (strcmp(stack1[n].AtmName,"H3") != 0) && (strcmp(stack1[n].AtmName,"CA") != 0) &&
						(strcmp(stack1[n].AtmName,"C") != 0) && (strcmp(stack1[n].AtmName,"O") != 0) && (strcmp(stack1[n].AtmName,"O1") != 0) &&
						(strcmp(stack1[n].AtmName,"O2") != 0)) && ((strcmp(stack2[o].AtmName,"N") != 0) && (strcmp(stack2[o].AtmName,"H") != 0) && (strcmp(stack2[o].AtmName,"H1") != 0) &&
						(strcmp(stack2[o].AtmName,"H2") != 0) && (strcmp(stack2[o].AtmName,"H3") != 0) && (strcmp(stack2[o].AtmName,"CA") != 0) &&
						(strcmp(stack2[o].AtmName,"C") != 0) && (strcmp(stack2[o].AtmName,"O") != 0) && (strcmp(stack2[o].AtmName,"O1") != 0) && (strcmp(stack2[o].AtmName,"O2") != 0)))
						{
							if (strncmp(stack1[n].AtmName,"H", 1) != 0 && strncmp(stack2[o].AtmName,"H", 1) != 0)
							{
								it[counter] = 2; /*possible hydropathic interaction*/
								ct0[counter] = stack1[n].Res_Type;
								ct1[counter] = stack2[o].Res_Type;
								at0[counter] = stack1[n].Atom_Type;
								at1[counter] = stack2[o].Atom_Type;
								
								//Assigning the ATOMIC SOLVATION PARAMETERS(Cummings, 1995) to the different atom types
								
								if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") == 0) || (strcmp(stack1[n].AtmName,"NZ") == 0))) //fully charged group
								{
									stack1[n].ASP_Type = 4;
								}
								else if ((strncmp(stack1[n].AtmName,"N", 1) == 0) && ((strcmp(stack1[n].AtmName,"NH2") != 0) && (strcmp(stack1[n].AtmName,"NZ") != 0)))
								{
									stack1[n].ASP_Type = 1;
								}
								else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") == 0) || (strcmp(stack1[n].AtmName,"OE2") == 0))) //fully charged group
								{
									stack1[n].ASP_Type = 3;
								}
								else if ((strncmp(stack1[n].AtmName,"O", 1) == 0) && ((strcmp(stack1[n].AtmName,"OD2") != 0) && (strcmp(stack1[n].AtmName,"OE2") != 0)))
								{
									stack1[n].ASP_Type = 1;
								}
								else if (strncmp(stack1[n].AtmName,"C", 1) == 0)
								{
									stack1[n].ASP_Type = 0;
								}
								else if (strncmp(stack1[n].AtmName,"S", 1) == 0)
								{
									stack1[n].ASP_Type = 2;
								}
								else
								{
									stack1[n].ASP_Type = 5;
								}
									
								if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") == 0) || (strcmp(stack2[o].AtmName,"NZ") == 0))) //fully charged group
								{
									stack2[o].ASP_Type = 4;
								}
								else if ((strncmp(stack2[o].AtmName,"N", 1) == 0) && ((strcmp(stack2[o].AtmName,"NH2") != 0) && (strcmp(stack2[o].AtmName,"NZ") != 0)))
								{
									stack2[o].ASP_Type = 1;
								}
								else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") == 0) || (strcmp(stack2[o].AtmName,"OE2") == 0))) //fully charged group
								{
									stack2[o].ASP_Type = 3;
								}
								else if ((strncmp(stack2[o].AtmName,"O", 1) == 0) && ((strcmp(stack2[o].AtmName,"OD2") != 0) && (strcmp(stack2[o].AtmName,"OE2") != 0)))
								{
									stack2[o].ASP_Type = 1;
								}
								else if (strncmp(stack2[o].AtmName,"C", 1) == 0)
								{
									stack2[o].ASP_Type = 0;
								}
								else if (strncmp(stack2[o].AtmName,"S", 1) == 0)
								{
									stack2[o].ASP_Type = 2;
								}									
								else
								{
									stack2[o].ASP_Type = 5;
								}
								asp0[counter] = stack1[n].ASP_Type;
								asp1[counter] = stack2[o].ASP_Type;								
								ro0[counter] = stack1[n].r0;
								ro1[counter] = stack2[o].r0;
								ri0[counter] = stack1[n].NodeIndex;
								ri1[counter] = stack2[o].NodeIndex;
								rd[counter] = sqrt(dist.rsq);
								counter++;
							}
						}
						/*inter_residue_potential[counter] = hydropathic_inter_residue_potentials[stack1[n].Res_Type][stack2[o].Res_Type];*/
						
						/*clashes_cutoff[counter] = dist.r01;
						contacts_cutoff[counter] = dist.r02;/*
						/*if (counter == MAXCONTACT - 1)
						{
							break;
							contactflag = 1;
						}*/
					}
				}
			
			}
		}		
	}
	/*for (counter = 0; counter < MAXCONTACT; counter++)
	{
		sum_inter_residue_potential = sum_inter_residue_potential + inter_residue_potential[counter];
	}*/
	
	return counter;
	//fclose(test);
}


int CONTACTS(void)
{
	int n,o,counter;
	
	counter=0;
	for (n=0;n<nAtom[0];n++)
	{
		for (o=0;o<nAtom[1];o++)
		{
							
			dist.r01 = HARDCORE_R*(Atoms[n][0].r0+Atoms[o][1].r0);
			dist.r02 = CONTACT_R*(Atoms[n][0].r0+Atoms[o][1].r0);
			dist.rsq = (Atoms[o][1].x-Atoms[n][0].x)*(Atoms[o][1].x-Atoms[n][0].x)+(Atoms[o][1].y-Atoms[n][0].y)*(Atoms[o][1].y-Atoms[n][0].y)+(Atoms[o][1].z-Atoms[n][0].z)*(Atoms[o][1].z-Atoms[n][0].z);
				 				
			if ((dist.rsq>=dist.r01*dist.r01)&&(dist.rsq<=dist.r02*dist.r02))
			{
				counter++;
			}
			
		}
							
	}
	return counter;
}


int CLASHES(void)
{
	int n,o,counter;
	
	counter=0;
	for (n=0;n<nAtom[0];n++)
	{
		for (o=0;o<nAtom[1];o++)
		{
							
			dist.r01 = HARDCORE_R*(Atoms[n][0].r0+Atoms[o][1].r0);
			dist.r02 = CONTACT_R*(Atoms[n][0].r0+Atoms[o][1].r0);
			dist.rsq = (Atoms[o][1].x-Atoms[n][0].x)*(Atoms[o][1].x-Atoms[n][0].x)+(Atoms[o][1].y-Atoms[n][0].y)*(Atoms[o][1].y-Atoms[n][0].y)+(Atoms[o][1].z-Atoms[n][0].z)*(Atoms[o][1].z-Atoms[n][0].z);
				 				
			if (dist.rsq<(dist.r01*dist.r01))
			{
				counter++;
			}
											
		}
							
	}
	return counter;
}

void READCHARGES(void)
{
	FILE *input;
	
	int i = 0;
	
	char s[256], szLine[MAXLENLINE], szAtomName[8], szResName[8], szAtomType[8], buff[256],*ReadStatus;
	
	float AtomCharge;
	
	sprintf(s, "charges.dat"); //file containing the GROMOS54A7 atomic point charges
	
	input = fopen(s,"r");
	
	if(input == NULL)	
	{
		printf("Fail to open charges file: %s\nQuit.\n",s);
		exit(0);
	}

	while (!feof(input))
	{	
		fscanf(input, "%s   %s   %s   %f\n", &szResName, &szAtomName, &szAtomType, &AtomCharge);				
		if(strncmp(szResName, "CYS2", 4)==0)
		{		
			strcpy(szResName, "CYS");
		}
		Charges[i].Res_Type = Get_AA_Type(szResName);
		strcpy(Charges[i].AtmName, szAtomName);
		strcpy(Charges[i].Atm_Type, szAtomType);
		Charges[i].Atm_Charge = AtomCharge;	
		i++;	
	}
	fclose(input);
}

/*double RESIDUECHARGES(int residues_index, int contacts)
{
	double q;
					
	if (residues_index == 34) 
	{
		q = -1.0;
	}
	else if (residues_index == 38)
	{
		q = -1.0;
	}
	else if (residues_index == 53) 
	{
		q = -1.0;
	}
	else if (residues_index == 59)
	{
		q = -1.0;
	}
	else if (residues_index == 96)
	{
		q = -1.0;
	}
	else if (residues_index == 98)
	{
		q = -1.0;
	}
	else if (residues_index == 99)
	{
		q = -1.0;
	}
	else if (residues_index == 16)
	{
		q = -1.0;
	}
	else if (residues_index == 36)
	{
		q = -1.0;
	}
	else if (residues_index == 44)
	{
		q = -1.0;
	}
	else if (residues_index == 47)
	{
		q = -1.0;
	}
	else if (residues_index == 50)
	{
		q = -1.0;
	}
	else if (residues_index == 69)
	{
		q = -1.0;
	}
	else if (residues_index == 74)
	{
		q = -1.0;
	}
	else if (residues_index == 77)
	{
		q = -1.0;
	}
	else if (residues_index == 13)
	{
		q = +0.05;
	}
	else if (residues_index == 31)
	{
		q = +0.10;
	}
	else if (residues_index == 51)
	{
		q = +0.17;
	}
	else if (residues_index == 84)
	{
		q = +0.02;
	}
	else if (residues_index == 1)
	{
		q = +0.31;
	}
	else if (strcmp(AA_Name[contacts],"LYS") == 0) 
	{
		q = +1.0;
	}
	else if (strcmp(AA_Name[contacts],"ARG") == 0)
	{
		q = +1.0;
	}
	else
	{
		q = 0.0;
	}			
	
	return q;
}*/

/*void READSASA(void)
{
	FILE *input;
	
	int i;
	
	double sasa_aa;
	
	char s[256];
	
	sprintf(s, "B2M_D76N_I2_relative_SASA.dat");
	
	input = fopen(s,"r");
	
	if(input == NULL)	
	{
		printf("Fail to open SASA file: %s\nQuit.\n",s);
		exit(0);
	}
	
	while (!feof(input))
	{	
		/*printf("tag1sasa\n");*/
		/*for (i = 0; i<MAXRESIDUE; i++)
		{
			/*printf("tag2sasa\n");*/
			/*fscanf(input, "%lf	  ", &sasa_aa);
			SASA_AA[i] = sasa_aa;
			/*printf("%.3lf\n", SASA_AA[i]);*/
		/*}
	}
	
	fclose(input);
}*/

void CENTRALIZER(double a)
{
	int i,j,k;

	//protein centralizer
	for (i=0;i<2;i++)
	{
		for (j=0;j<3;j++)
		{
			MEAN_R[j][i]=0;
		}
	}
	
	
	for (i=0;i<2;i++)
	{
		for (j=0;j<nAtom[i];j++)
		{
			MEAN_R[0][i]=MEAN_R[0][i]+Atoms[j][i].x;
			MEAN_R[1][i]=MEAN_R[1][i]+Atoms[j][i].y;
			MEAN_R[2][i]=MEAN_R[2][i]+Atoms[j][i].z;
		}
		for (j=0;j<3;j++)
		{
			MEAN_R[j][i]=MEAN_R[j][i]/nAtom[i];
		}
	}
	
	//ROUNDING ERRORS DUE TO PDB FILE FORMAT
	//for (i=0;i<2;i++)
	//{
	//	printf("%d	%lf	%lf	%lf\n",i,MEAN_R[0][i],MEAN_R[1][i],MEAN_R[2][i]);
	//}
	//printf("\n\n");
	
	
	for (i=0;i<2;i++)
	{
		for (j=0;j<nAtom[i];j++)
		{
			Atoms[j][i].x = Atoms[j][i].x - MEAN_R[0][i];
			Atoms[j][i].y = Atoms[j][i].y - MEAN_R[1][i];
			Atoms[j][i].z = Atoms[j][i].z - MEAN_R[2][i];
		}
	}
	for (i=0;i<2;i++)
	{
		for (j=0;j<3;j++)
		{
			MEAN_R[j][i]=0;
		}
	}
	for (i=0;i<2;i++)
	{
		for (j=0;j<nAtom[i];j++)
		{
			MEAN_R[0][i]=MEAN_R[0][i]+Atoms[j][i].x;
			MEAN_R[1][i]=MEAN_R[1][i]+Atoms[j][i].y;
			MEAN_R[2][i]=MEAN_R[2][i]+Atoms[j][i].z;
		}
		for (j=0;j<3;j++)
		{
			MEAN_R[j][i]=MEAN_R[j][i]/nAtom[i];
		}
	}
		
	//CONTROL OUTPUT: CENTRALIZED PROTEINS with CMs at (0,0,0):
	//for (i=0;i<2;i++)
	//{
	//	printf("%d	%lf	%lf	%lf\n",i,MEAN_R[0][i],MEAN_R[1][i],MEAN_R[2][i]);
	//}
	
	
	//MOVE STRUCTURE 2 along X/Y/Z-axis for a=Rg1+Rg2 in DIR-direction //
	/////////////////////////////////////////////////////////////////
	if (DIR==1)
	{
		for (j=0;j<nAtom[1];j++)
		{
			Atoms[j][1].x = Atoms[j][1].x + a;
		}
		INITIALORIENTATOR();
	}
	else
	{
		if (DIR==-1)
		{
			for (j=0;j<nAtom[1];j++)
			{
				Atoms[j][1].x = Atoms[j][1].x - a;
			}
			INITIALORIENTATOR();
		}
		else
		{
			if (DIR == 2)
			{
				for (j=0;j<nAtom[1];j++)
				{
					Atoms[j][1].y = Atoms[j][1].y + a;
				}
				INITIALORIENTATOR();
			}
			else
			{
				if (DIR == -2)
				{
					for (j=0;j<nAtom[1];j++)
					{
						Atoms[j][1].y = Atoms[j][1].y - a;
					}
					INITIALORIENTATOR();
				}
				else
				{
					if (DIR == 3)
					{
						for (j=0;j<nAtom[1];j++)
						{
							Atoms[j][1].z = Atoms[j][1].z + a;
						}
						INITIALORIENTATOR();	
					}
					else
					{
						if (DIR == -3)
						{
							for (j=0;j<nAtom[1];j++)
							{
								Atoms[j][1].z = Atoms[j][1].z - a;
							}
							INITIALORIENTATOR();
						}
					}
				
				}
			}
		
		
		}
	}	
	
	
	//CONTROL OUTPUT: MOVED PROTEINS - CMs
	//PROTEIN 1: (0,0,0)
	//PROTEIN 2: (restdist,0,0)
	
	//for (i=0;i<2;i++)
	//{
	//	printf("%d	%lf	%lf	%lf\n",i,MEAN_R[0][i],MEAN_R[1][i],MEAN_R[2][i]);
	//}
	
}


void INITIALORIENTATOR(void)
{
	int i,j,k;
	double oldmoved_x[MAXATOMNUM];
	double oldmoved_y[MAXATOMNUM];
	double oldmoved_z[MAXATOMNUM];		

	//recalculating the coordinates of the structure stored in Atoms in relation to its centre of mass
	
	for (i=0;i<2;i++)
	{
		for (j=0;j<3;j++)
		{
			MEAN_R[j][i]=0;
		}
	}
		
	for (i=0;i<2;i++)
	{
		for (j=0;j<nAtom[i];j++)
		{
			MEAN_R[0][i]=MEAN_R[0][i]+Atoms[j][i].x;
			MEAN_R[1][i]=MEAN_R[1][i]+Atoms[j][i].y;
			MEAN_R[2][i]=MEAN_R[2][i]+Atoms[j][i].z;
		}
		for (j=0;j<3;j++)
		{
			MEAN_R[j][i]=MEAN_R[j][i]/nAtom[i];
		}
	}
		
	for (i=0;i<nAtom[1];i++)
	{
		RotMoved[i].x = Atoms[i][1].x - MEAN_R[0][1];
		RotMoved[i].y = Atoms[i][1].y - MEAN_R[1][1];
		RotMoved[i].z = Atoms[i][1].z - MEAN_R[2][1];
	}

	//performing -/+90º rotation over each axis on structure stored in RotMoved via matrix multiplication
	
	
	if (ORIENTATION == 1)
	{	
		for (i=0;i<nAtom[1];i++)
		{
			oldmoved_x[i] = RotMoved[i].x;
			oldmoved_y[i] = RotMoved[i].y;
			oldmoved_z[i] = RotMoved[i].z;
			RotMoved[i].x = oldmoved_x[i];
			RotMoved[i].y = -oldmoved_z[i];
			RotMoved[i].z = oldmoved_y[i];
		}
	}
	else if (ORIENTATION == -1)
	{
		for (i=0;i<nAtom[1];i++)
		{
			oldmoved_x[i] = RotMoved[i].x;
			oldmoved_y[i] = RotMoved[i].y;
			oldmoved_z[i] = RotMoved[i].z;
			RotMoved[i].x = oldmoved_x[i];
			RotMoved[i].y = oldmoved_z[i];
			RotMoved[i].z = -oldmoved_y[i];
		}
	}
	else if (ORIENTATION == 2)
	{
		for (i=0;i<nAtom[1];i++)
		{
			oldmoved_x[i] = RotMoved[i].x;
			oldmoved_y[i] = RotMoved[i].y;
			oldmoved_z[i] = RotMoved[i].z;
			RotMoved[i].x = oldmoved_z[i];
			RotMoved[i].y = oldmoved_y[i];
			RotMoved[i].z = -oldmoved_x[i];
		}
	}
	else if (ORIENTATION == -2)
	{
		for (i=0;i<nAtom[1];i++)
		{
			oldmoved_x[i] = RotMoved[i].x;
			oldmoved_y[i] = RotMoved[i].y;
			oldmoved_z[i] = RotMoved[i].z;
			RotMoved[i].x = -oldmoved_z[i];
			RotMoved[i].y = oldmoved_y[i];
			RotMoved[i].z = oldmoved_x[i];
		}
	}
	else if (ORIENTATION == 3)
	{
		for (i=0;i<nAtom[1];i++)
		{
			oldmoved_x[i] = RotMoved[i].x;
			oldmoved_y[i] = RotMoved[i].y;
			oldmoved_z[i] = RotMoved[i].z;
			RotMoved[i].x = oldmoved_y[i];
			RotMoved[i].y = -oldmoved_x[i];
			RotMoved[i].z = oldmoved_z[i];
		}
	}
	else if (ORIENTATION == -3)
	{
		for (i=0;i<nAtom[1];i++)
		{
			oldmoved_x[i] = RotMoved[i].x;
			oldmoved_y[i] = RotMoved[i].y;
			oldmoved_z[i] = RotMoved[i].z;
			RotMoved[i].x = -oldmoved_y[i];
			RotMoved[i].y = oldmoved_x[i];
			RotMoved[i].z = oldmoved_z[i];
		}
	}
	else
	{
		if (ORIENTATION == 4)
		{
			for (i=0;i<nAtom[1];i++)
			{
				oldmoved_x[i] = RotMoved[i].x;
				oldmoved_y[i] = RotMoved[i].y;
				oldmoved_z[i] = RotMoved[i].z;
				RotMoved[i].x = oldmoved_x[i];
				RotMoved[i].y = oldmoved_y[i];
				RotMoved[i].z = oldmoved_z[i];
			}
		}
	}
				
		
	//recalculating the coordinates of the structure stored in RotMoved in relation to the origin of the cartesian axes
		
	for (i=0;i<nAtom[1];i++)
	{
		Atoms[i][1].x = RotMoved[i].x + MEAN_R[0][1];
		Atoms[i][1].y = RotMoved[i].y + MEAN_R[1][1];
		Atoms[i][1].z = RotMoved[i].z + MEAN_R[2][1];
	}

	////////////////////////////////////////////////////////////////
	for (i=0;i<2;i++)
	{
		for (j=0;j<3;j++)
		{
			MEAN_R[j][i]=0;
		}
	}
	for (i=0;i<2;i++)
	{
		for (j=0;j<nAtom[i];j++)
		{
			MEAN_R[0][i]=MEAN_R[0][i]+Atoms[j][i].x;
			MEAN_R[1][i]=MEAN_R[1][i]+Atoms[j][i].y;
			MEAN_R[2][i]=MEAN_R[2][i]+Atoms[j][i].z;
		}
		for (j=0;j<3;j++)
		{
			MEAN_R[j][i]=MEAN_R[j][i]/nAtom[i];
		}
	}
}	

	
void READER(int a, int b)
{
	FILE *input;
	
	int i,j,k;
	int ReadCoord, HIS_Type;
	
	int cter;

	char s[256], szLine[MAXLENLINE], szAtomName[8], szResName[8], buff[256],*ReadStatus;
	
	for (i=0;i<2;i++)
	{
		input=NULL;
		Conf[i].contEn=0;
		Conf[i].Rg=0;
		Conf[i].RMSD=0;
		ResNum[i]=0;
		nAtom[i]=0;
		
		if (i==0)
		{
			j=a;
		}
		else
		{
			j=b;
		}
		
		if(j < 10)	
		{
			sprintf(s, "D76N-I1-pH5p2-chain-newcf-000%d.pdb",j);
		}
		else 
			if(j < 100)	
			{
				sprintf(s, "D76N-I1-pH5p2-chain-newcf-00%d.pdb",j);
			}
			else 
				if(j < 1000)	
				{
					sprintf(s, "D76N-I1-pH5p2-chain-newcf-0%d.pdb",j);
				}
				else	
				{
					sprintf(s, "D76N-I1-pH5p2-chain-newcf-%d.pdb",j);
				}
	
		
		
		
		//printf("Filename: [%s]\n", s); fflush(stdout);
		
		///////////////////////////////////////////////////////////////////////////////
		
			
		input = fopen(s,"r");
		
			
		if(input == NULL)	
		{
			printf("Fail to open PDB file: %s\nQuit.\n",s);
			exit(0);
		}
	
		k=0;
		while(1)
		{
			
			ReadStatus = fgets(szLine, MAXLENLINE, input);
			
			if(ReadStatus == NULL)	
			{
				break;
			}
			
			if (k==1)
			{
				ReadCoord=sscanf(szLine+10, "%d", &(Conf[i].contEn));
				ReadCoord=sscanf(szLine+24, "%lf", &(Conf[i].Rg));
				ReadCoord=sscanf(szLine+43, "%lf", &(Conf[i].RMSD));
			}
								
			if (strncmp(szLine, "ATOM", 4)==0)
			{
				sscanf(szLine+5, "%d", &(Atoms[nAtom[i]][i].AtomIndex));  
				sscanf(szLine+12, "%s", szAtomName);
				sscanf(szLine+26, "%d", &(Atoms[nAtom[i]][i].NodeIndex));
											
				if(szLine[17] == ' ')	
				{
					sscanf(szLine+17, "%s", szResName);
				}
				else if(szLine[17] == 'A')	//only read atom in chain A
				{	
					sscanf(szLine+18, "%s", szResName);
				}
				else
				{
					continue;
				}
				
				if(strncmp(szResName, "ARGH", 4)==0)
				{
					strcpy(szResName, "ARG");
				}
				if(strncmp(szResName, "LYSH", 4)==0)
				{
					strcpy(szResName, "LYS");
				}												
				if((strncmp(szResName, "ILE", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 1)) //N-terminal residue: 
				{
					strcpy(szResName, "ILEN");
				}
				if((strncmp(szResName, "MET", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 99)) //C-terminal residue:
				{
					strcpy(szResName, "METC");
				}				
				strcpy(Atoms[nAtom[i]][i].AtmName, szAtomName);
				
				//ASSIGNING the atom type for all side-chain atoms
				if((strcmp(szAtomName, "N")!=0) && (strncmp(szAtomName, "H", 1)!=0) && (strncmp(szAtomName, "1", 1)!=0) && (strncmp(szAtomName, "2", 1)!=0) && (strcmp(szAtomName, "CA")!=0) && (strcmp(szAtomName, "C")!=0) && (strcmp(szAtomName, "O")!=0) && (strcmp(szAtomName, "O1")!=0) && (strcmp(szAtomName, "O2")!=0))
				{
					Atoms[nAtom[i]][i].Atom_Type = Get_Atom_Type(szAtomName);
				}
				Atoms[nAtom[i]][i].Res_Type = Get_AA_Type(szResName);
					
				ReadCoord = sscanf(szLine+30, "%lf %lf %lf", &(Atoms[nAtom[i]][i].x), &(Atoms[nAtom[i]][i].y), &(Atoms[nAtom[i]][i].z));
			
				if(ReadCoord != 3)	
				{
					printf("Fail to read coordinates from PDB file: %s, Atom %d AtomName = %s ResName = %s\nQuit\n", 
					s, nAtom[i]+1, szAtomName, szResName);
					fclose(input);
					exit(1);
				}
			
				nAtom[i]++;
					
			}
					
			ResNum[i] = Atoms[nAtom[i]-1][i].NodeIndex - Atoms[0][i].NodeIndex + 1;
			k++;									
		}
		fclose(input);
		AssignAtomRadii();
		
	}
	
	//printf("READING COMPLETE\n\n");
	

}

/*void READPOTENTIALS(void)
{
	FILE *input;

	int i, j;
	double k;

	input = fopen(matrix_MJ_inter-residues_potentials.dat, "r");

	if (input == NULL)
	{
		printf("Fail to open dat file: matrix_MJ_inter-residues_potentials.dat\nQuit.\n");
		exit(0);
	}

	while (!feof(input))
	{
		for (i = 0; i < NUM_AA_TYPE; i++)
		{
			for (j = 0; j < NUM_AA_TYPE; j++)
			{
				scanf(input, "%lf", & k);
				MJ_inter-residue_potentials[i][j] = k;
			}
		}
	}
}
*/
void ExportSnapshot(int Index1, int Index2)	//export the current protein structure. It can be used for resumed simulation from unexpected termination
{
	FILE *fOut1;
	FILE *fOut2;
	
	int i;
	char szName1[256];
	char szName2[256];
	char AA_Name_0[5];
	char AA_Name_1[5];
		
	
	sprintf(szName1, "D76N_I1_PAIRS_ID1=%d_ID2=%d_DIR=%d_newcf__a.pdb", Index1,Index2,(int)DIR);
	sprintf(szName2, "D76N_I1_PAIRS_ID1=%d_ID2=%d_DIR=%d_newcf__b.pdb", Index1,Index2,(int)DIR);

	
	fOut1 = fopen(szName1, "w");
	fOut2 = fopen(szName2, "w");
	
	fprintf(fOut1, "REMARK PDB file generated by GO potential code written by Lei Huang.\n");
	//fprintf(fOut1, "REMARK E = %7.3lf Rg = %7.3lf  CaRMSD = %7.3lf \n", (double)Conf[0].contEn, Conf[0].Rg, Conf[0].RMSD);
	for(i=0; i<nAtom[0]; i++)	
	{
		strcpy(AA_Name_0, AA_Name[Atoms[i][0].Res_Type]);
		
		if(strncmp(AA_Name_0, "ARGH", 4)==0)
		{
			strcpy(AA_Name_0, "ARG");
		}
		if(strncmp(AA_Name_0, "LYSH", 4)==0)
		{
			strcpy(AA_Name_0, "LYS");
		}
		if(strncmp(AA_Name_0, "HISH", 4)==0)
		{
			strcpy(AA_Name_0, "HIS");
		}		
		if(strncmp(AA_Name_0, "METC", 4)==0)
		{
			strcpy(AA_Name_0, "MET");
		}
		if(strncmp(AA_Name_0, "ILEN", 4)==0)
		{
			strcpy(AA_Name_0, "ILE");
		}
		
		if((strncmp(Atoms[i][0].AtmName, "HD11", 4)==0) || (strncmp(Atoms[i][0].AtmName, "HD12", 4)==0) || (strncmp(Atoms[i][0].AtmName, "HD21", 4)==0) || (strncmp(Atoms[i][0].AtmName, "HD22", 4)==0) || (strncmp(Atoms[i][0].AtmName, "HE11", 4)==0) || (strncmp(Atoms[i][0].AtmName, "HE12", 4)==0) || (strncmp(Atoms[i][0].AtmName, "HE21", 4)==0) || (strncmp(Atoms[i][0].AtmName, "HE22", 4)==0) || (strncmp(Atoms[i][0].AtmName, "HH11", 4)==0) || (strncmp(Atoms[i][0].AtmName, "HH12", 4)==0) || (strncmp(Atoms[i][0].AtmName, "HH21", 4)==0) || (strncmp(Atoms[i][0].AtmName, "HH22", 4)==0))
		{
			fprintf(fOut1, "ATOM%7d %-4s %-3s A%4d    %8.3lf%8.3lf%8.3lf  1.00  0.00\n", 
			i+1, Atoms[i][0].AtmName, AA_Name_0, Atoms[i][0].NodeIndex, Atoms[i][0].x,Atoms[i][0].y,Atoms[i][0].z);
		}
		else
		{
			fprintf(fOut1, "ATOM%7d  %-3s %-3s A%4d    %8.3lf%8.3lf%8.3lf  1.00  0.00\n", 
			i+1, Atoms[i][0].AtmName, AA_Name_0, Atoms[i][0].NodeIndex, Atoms[i][0].x,Atoms[i][0].y,Atoms[i][0].z);
		}		
				
	}
	
	fprintf(fOut2, "REMARK PDB file generated by GO potential code written by Lei Huang.\n");
	//fprintf(fOut2, "REMARK E = %7.3lf Rg = %7.3lf  CaRMSD = %7.3lf \n", (double)Conf[1].contEn, Conf[1].Rg, Conf[1].RMSD);
	for(i=0; i<nAtom[1]; i++)	
	{
		strcpy(AA_Name_1, AA_Name[Atoms[i][1].Res_Type]);
		
		if(strncmp(AA_Name_1, "ARGH", 4)==0)
		{
			strcpy(AA_Name_1, "ARG");
		}
		if(strncmp(AA_Name_1, "LYSH", 4)==0)
		{
			strcpy(AA_Name_1, "LYS");
		}
		if(strncmp(AA_Name_1, "HISH", 4)==0)
		{
			strcpy(AA_Name_1, "HIS");
		}				
		if(strncmp(AA_Name_1, "METC", 4)==0)
		{
			strcpy(AA_Name_1, "MET");
		}
		if(strncmp(AA_Name_1, "ILEN", 4)==0)
		{
			strcpy(AA_Name_1, "ILE");
		}
		
		if((strncmp(Atoms[i][1].AtmName, "HD11", 4)==0) || (strncmp(Atoms[i][1].AtmName, "HD12", 4)==0) || (strncmp(Atoms[i][1].AtmName, "HD21", 4)==0) || (strncmp(Atoms[i][1].AtmName, "HD22", 4)==0) || (strncmp(Atoms[i][1].AtmName, "HE11", 4)==0) || (strncmp(Atoms[i][1].AtmName, "HE12", 4)==0) || (strncmp(Atoms[i][1].AtmName, "HE21", 4)==0) || (strncmp(Atoms[i][1].AtmName, "HE22", 4)==0) || (strncmp(Atoms[i][1].AtmName, "HH11", 4)==0) || (strncmp(Atoms[i][1].AtmName, "HH12", 4)==0) || (strncmp(Atoms[i][1].AtmName, "HH21", 4)==0) || (strncmp(Atoms[i][1].AtmName, "HH22", 4)==0))
		{
			fprintf(fOut2, "ATOM%7d %-4s %-3s A%4d    %8.3lf%8.3lf%8.3lf  1.00  0.00\n", 
			i+1, Atoms[i][1].AtmName, AA_Name_1, Atoms[i][1].NodeIndex, Atoms[i][1].x,Atoms[i][1].y,Atoms[i][1].z);
		}
		else
		{
			fprintf(fOut2, "ATOM%7d  %-3s %-3s A%4d    %8.3lf%8.3lf%8.3lf  1.00  0.00\n", 
			i+1, Atoms[i][1].AtmName, AA_Name_1, Atoms[i][1].NodeIndex, Atoms[i][1].x,Atoms[i][1].y,Atoms[i][1].z);
		}		
				
	}
	
	fclose(fOut1);
	fclose(fOut2);
	
}

//pdb-compatible screen output
//for (i=0;i<nAtom;i++)
//{
//	printf("ATOM%7d	%-3s	%-3s	A%4d	%8.3lf	%8.3lf	%8.3lf	1.00 0.00\n",
//	i+1, Atoms[i].AtmName, AA_Name[Atoms[i].Res_Type],Atoms[i].NodeIndex, Atoms[i].x,Atoms[i].y,Atoms[i].z);
//}

//output of r0-assignment, tested, works!
//for (i=0;i<2;i++)
//{
//	for (j=0;j<nAtom[i];j++)
//	{
//		printf("%d  ATOM%7d	%-3s	%-3s	A%4d	%8.3lf	%8.3lf	%8.3lf	1.00 %8.3lf\n",i,j+1, Atoms[j][i].AtmName,
//		AA_Name[Atoms[j][i].Res_Type],Atoms[j][i].NodeIndex, Atoms[j][i].x, Atoms[j][i].y, Atoms[j][i].z, Atoms[j][i].r0);
//	}
//}
	
//control output of protein folding observables
//printf("%d	%d	%d\n",Conf[0].contEn,Conf[1].contEn,Conf[0].contEn+Conf[1].contEn);
//printf("%lf	%lf	%lf\n",Conf[0].Rg,Conf[1].Rg,Conf[0].Rg+Conf[1].Rg);
//printf("%lf	%lf	%lf\n\n",Conf[0].RMSD,Conf[1].RMSD,Conf[0].RMSD+Conf[1].RMSD);

//TEST OUTPUT for stack structure
/*for (j=0;j<nAtom[1];j++)
{
	printf("ATOM%7d	%-3s	%-3s	A%4d	%8.3lf	%8.3lf	%8.3lf	1.00 %8.3lf\n",j+1, Moved[j].AtmName,
		AA_Name[Moved[j].Res_Type],Moved[j].NodeIndex, Moved[j].x, Moved[j].y, Moved[j].z, Moved[j].r0);
}*/

double abbs(double a)
{
	if (a < 0)
		return -a;
	else
		return a;
}
int IsANumber(char c)			//to check whether c is a number [0-9] or not
{
	if((c >= '0') && (c<='9'))	
	{
		return 1;
	}
	else	
	{
		return 0;
	}
}
int Get_AA_Type(char SzAAName[])
{
	int AA_Type=-1;
	int i;
	
	for (i=0;i<NUM_AA_TYPE; i++)
	{
		if (strcmp(SzAAName, AA_Name[i])==0)
		{
			AA_Type=i;
			break;
		}	
	}
	
	if (AA_Type<0)
	{
		printf("Fatal error in AA matching!");
		exit(1);
	}
	return AA_Type;
}
int Get_Atom_Type(char SzAtomName[])
{
	int Atom_Type=-1;
	int i;
	
	for (i=0;i<NUM_Atom_TYPE; i++)
	{
		if (strcmp(SzAtomName, atom_name[i])==0)
		{
			Atom_Type=i;
			break;
		}	
	}
	
	if (Atom_Type<0)
	{
		printf("Fatal error in Atom matching!");
		exit(1);
	}
	return Atom_Type;
}
void AssignAtomRadii(void)
{
	int n,i, Type;
	
	for (i=0;i<2;i++)
	{
		for (n=0;n<nAtom[i];n++)
		{
			Type = QueryAtomType(Atoms[n][i].AtmName, AA_Name[Atoms[n][i].Res_Type]);
			Atoms[n][i].r0=radii[Type];
		}
	}


}
int QueryAtomType(char *s, char *res)	//from the atom name, s, and residue name, res, to the type of atoms. The type is used to assign the radii. This procedure is same as the REMC code 
{
	if (!strncmp(s,"C",1)) {
		if (!strcmp(s,"C") || (!strcmp(s,"CG") && (!strcmp(res,"ASN") || !strcmp(res,"ASP") || !strcmp(res,"ASPH") || !strcmp(res,"HIS") || !strcmp(res,"HISH") || !strcmp(res,"PHE") || !strcmp(res,"TYR") ||
		!strcmp(res,"TRP"))) || (!strcmp(s,"CD") && (!strcmp(res,"GLN") || !strcmp(res,"GLU") || !strcmp(res,"GLUH"))) || (!strcmp(s,"CZ") && (!strcmp(res,"ARG") || !strcmp(res,"ARGH") || !strcmp(res,"TYR"))) ||  (!strcmp(s,"CD2") && !strcmp(res,"TRP")) || (!strcmp(s,"CE2") && !strcmp(res,"TRP")))
			return 0;
		else if ((!strcmp(s,"CD1") && (!strcmp(res,"PHE") || !strcmp(res,"TYR"))) || (!strcmp(s,"CD2") && (!strcmp(res,"HIS") || !strcmp(res,"HISH") || !strcmp(res,"PHE") || !strcmp(res,"TYR"))) || (!strcmp(s,"CZ") && !strcmp(res,"PHE")) || (!strcmp(res,"TRP") && (!strncmp(s,"CH",2) || !strncmp(s,"CZ",2) || !strcmp(s,"CE3") || !strcmp(s,"CD1"))))
			return 1;
		else if ((!strcmp(s,"CA") && strcmp(res,"GLY")) || (!strcmp(s,"CB") && (!strcmp(res,"ILE") || !strcmp(res,"ILEN") || !strcmp(res,"THR") || !strcmp(res,"VAL"))) || (!strcmp(s,"CG") && !strcmp(res,"LEU")))
			return 2;
		else if ((!strcmp(s,"CB") && !strcmp(res,"ALA")) || (!strcmp(s,"CD1") && (!strcmp(res,"ILE") || !strcmp(res,"ILEN") || !strcmp(res,"LEU"))) || (!strcmp(s,"CD2") &&
		!strcmp(res,"LEU")) || (!strcmp(s,"CG1") && !strcmp(res,"VAL")) || (!strcmp(s,"CG2") && (!strcmp(res,"ILE") || !strcmp(res,"ILEN") || !strcmp(res,"THR") || !strcmp(res,"VAL"))))
			return 4;
		else
			return 3;
	}
	else if (!strncmp(s,"N",1)) {
		if ((!strcmp(s,"N") && !strcmp(res,"PRO")) || (!strcmp(s,"NE2") && (!strcmp(res,"HIS") || !strcmp(res,"HISH"))))
			return 5;
		else if (!strcmp(s,"NZ"))
			return 8;
		else if (!strcmp(s,"N") || (!strcmp(s,"ND1") && (!strcmp(res,"HIS") || !strcmp(res,"HISH"))) || (!strcmp(s,"NE") && (!strcmp(res,"ARG") || !strcmp(res,"ARGH"))) || (!strcmp(s,"NE1") && !strcmp(res,"TRP")))
			return 6;
		else
			return 7;
	}
	else if (!strncmp(s,"O",1)) {
		if (!strcmp(s,"OH") || !strcmp(s,"OG") || !strcmp(s,"OG1"))
			return 10;
		else
			return 9;
	}
	else if (!strncmp(s,"S",1)) {
		if (!strcmp(res,"MET") || !strcmp(res,"METC") || !strcmp(res,"METN"))
			return 11;
		else
			return 12;
	}
	else if (!strncmp(s,"H",1)) 
	{
		if ((!strcmp(s,"HD1") || !strcmp(s,"HD2") || !strcmp(s,"HE1") || !strcmp(s,"HE2") || !strcmp(s,"HZ")) && !strcmp(res,"PHE"))
		{
			return 14;
		}
		else if ((!strcmp(s,"HD1") || !strcmp(s,"HE3") || !strcmp(s,"HZ2") || !strcmp(s,"HZ3") || !strcmp(s,"HH2")) && !strcmp(res,"TRP"))
		{
			return 14;
		}
		else if ((!strcmp(s,"HD1") || !strcmp(s,"HD2") || !strcmp(s,"HE1") || !strcmp(s,"HE2")) && !strcmp(res,"TYR"))
		{
			return 14;
		}
		else if ((!strcmp(s,"HD2") || !strcmp(s,"HE1")) && !strcmp(res,"HIS"))
		{
			return 14;
		}
		else if ((!strcmp(s,"HD2") || !strcmp(s,"HE1")) && !strcmp(res,"HISH"))		
		{
			return 14;
		}
		else
		{
			return 13;
		}
	}
	else
	{
		return -500;
	}	
}


//RANDOM NUMBER GENERATOR FUNCTIONS
#define FALSE 0
#define TRUE 1

/*
   This Random Number Generator is based on the algorithm in a FORTRAN
   version published by George Marsaglia and Arif Zaman, Florida State
   University; ref.: see original comments below.
   At the fhw (Fachhochschule Wiesbaden, W.Germany), Dept. of Computer
   Science, we have written sources in further languages (C, Modula-2
   Turbo-Pascal(3.0, 5.0), Basic and Ada) to get exactly the same test
   results compared with the original FORTRAN version.
   April 1989
   Karl-L. Noell <NOELL@DWIFH1.BITNET>
      and  Helmut  Weber <WEBER@DWIFH1.BITNET>

   This random number generator originally appeared in "Toward a Universal
   Random Number Generator" by George Marsaglia and Arif Zaman.
   Florida State University Report: FSU-SCRI-87-50 (1987)
   It was later modified by F. James and published in "A Review of Pseudo-
   random Number Generators"
   THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
   (However, a newly discovered technique can yield
   a period of 10^600. But that is still in the development stage.)
   It passes ALL of the tests for random number generators and has a period
   of 2^144, is completely portable (gives bit identical results on all
   machines with at least 24-bit mantissas in the floating point
   representation).
   The algorithm is a combination of a Fibonacci sequence (with lags of 97
   and 33, and operation "subtraction plus one, modulo one") and an
   "arithmetic sequence" (using subtraction).

   Use IJ = 1802 & KL = 9373 to test the random number generator. The
   subroutine RANMAR should be used to generate 20000 random numbers.
   Then display the next six random numbers generated multiplied by 4096*4096
   If the random number generator is working properly, the random numbers
   should be:
           6533892.0  14220222.0  7275067.0
           6172232.0  8354498.0   10633180.0
*/

/* Globals */
double u[97],c,cd,cm;
int i97,j97;
int test = FALSE;

/*
   This is the initialization routine for the random number generator.
   NOTE: The seed variables can have values between:    0 <= IJ <= 31328
                                                        0 <= KL <= 30081
   The random number sequences created by these two seeds are of sufficient
   length to complete an entire calculation with. For example, if sveral
   different groups are working on different parts of the same calculation,
   each group could be assigned its own IJ seed. This would leave each group
   with 30000 choices for the second seed. That is to say, this random
   number generator can create 900 million different subsequences -- with
   each subsequence having a length of approximately 10^30.
*/
void RandomInitialise(int ij,int kl)
{
   double s,t;
   int ii,i,j,k,l,jj,m;

   /*
      Handle the seed range errors
         First random number seed must be between 0 and 31328
         Second seed must have a value between 0 and 30081
   */
   if (ij < 0 || ij > 31328 || kl < 0 || kl > 30081) {
		ij = 1802;
		kl = 9373;
   }

   i = (ij / 177) % 177 + 2;
   j = (ij % 177)       + 2;
   k = (kl / 169) % 178 + 1;
   l = (kl % 169);

   for (ii=0; ii<97; ii++) {
      s = 0.0;
      t = 0.5;
      for (jj=0; jj<24; jj++) {
         m = (((i * j) % 179) * k) % 179;
         i = j;
         j = k;
         k = m;
         l = (53 * l + 1) % 169;
         if (((l * m % 64)) >= 32)
            s += t;
         t *= 0.5;
      }
      u[ii] = s;
   }

   c    = 362436.0 / 16777216.0;
   cd   = 7654321.0 / 16777216.0;
   cm   = 16777213.0 / 16777216.0;
   i97  = 97;
   j97  = 33;
   test = TRUE;
}

/* 
   This is the random number generator proposed by George Marsaglia in
   Florida State University Report: FSU-SCRI-87-50
*/
double RandomUniform(void)
{
   double uni;

   /* Make sure the initialisation routine has been called */
   if (!test) 
   	RandomInitialise(1802,9373);

   uni = u[i97-1] - u[j97-1];
   if (uni <= 0.0)
      uni++;
   u[i97-1] = uni;
   i97--;
   if (i97 == 0)
      i97 = 97;
   j97--;
   if (j97 == 0)
      j97 = 97;
   c -= cd;
   if (c < 0.0)
      c += cm;
   uni -= c;
   if (uni < 0.0)
      uni++;

   return(uni);
}

/*
  ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
  THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
  The function returns a normally distributed pseudo-random number
  with a given mean and standard devaiation.  Calls are made to a
  function subprogram which must return independent random
  numbers uniform in the interval (0,1).
  The algorithm uses the ratio of uniforms method of A.J. Kinderman
  and J.F. Monahan augmented with quadratic bounding curves.
*/
double RandomGaussian(double mean,double stddev)
{
   double  q,u,v,x,y;

	/*  
		Generate P = (u,v) uniform in rect. enclosing acceptance region 
      Make sure that any random numbers <= 0 are rejected, since
      gaussian() requires uniforms > 0, but RandomUniform() delivers >= 0.
	*/
   do {
      u = RandomUniform();
      v = RandomUniform();
   	if (u <= 0.0 || v <= 0.0) {
       	u = 1.0;
       	v = 1.0;
   	}
      v = 1.7156 * (v - 0.5);

      /*  Evaluate the quadratic form */
      x = u - 0.449871;
   	y = fabs(v) + 0.386595;
      q = x * x + y * (0.19600 * y - 0.25472 * x);

      /* Accept P if inside inner ellipse */
      if (q < 0.27597)
			break;

      /*  Reject P if outside outer ellipse, or outside acceptance region */
    } while ((q > 0.27846) || (v * v > -4.0 * log(u) * u * u));

    /*  Return ratio of P's coordinates as the normal deviate */
    return (mean + stddev * v / u);
}

/*
   Return random integer within a range, lower -> upper INCLUSIVE
*/
int RandomInt(int lower,int upper)
{
   return((int)(RandomUniform() * (upper - lower + 1)) + lower);
}

/*
   Return random float within a range, lower -> upper
*/
double RandomDouble(double lower,double upper)
{
   return((upper - lower) * RandomUniform() + lower);
}



