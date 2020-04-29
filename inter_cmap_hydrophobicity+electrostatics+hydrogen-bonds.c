#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
//=============Parameters you might want to change===============

//ENSSTRUCTS = number of conformations in intermediate ensembles
//B2M-WT-NAT(pH5.2)_new: 4005
//B2M-WT-NAT(pH6.2)_new: 4005
//B2M-WT-NAT(pH7.2)_new: 4005
//NAT(D76N)pH5.2=4005
//NAT(D76N)pH6.2=4005
//NAT(D76N)pH7.2=4005
//I1(D76N)pH5.2=4005
//I1(D76N)pH6.2=4005
//I1(D76N)pH7.2=4005
//I2(D76N)pH5.2=4005
//I2(D76N)pH6.2=4005
//I2(D76N)pH7.2=4005
//D76N-NAT(DMD Go model) = 1399

#define ALLPAIRS        (1000)

#define HARDCORE_R	     (1.)
#define CONTACT_R	     (2.33)
#define CONTACT_R_hydropathy (1.60)
#define LOCAL_RANGE	     (0)

//end to parameters you might want to change=====================

#define MAXLENLINE 	(256)
#define LENCHARGESFILE  (287)
#define MAXATOMNUM	(1100)
#define MAXRESIDUE      (100)
#define MAXCONTACT  	(50000)
#define HIS_TYPE_HIS	(0)
#define HIS_TYPE_HSD	(1)
#define HIS_TYPE_HSC	(2)
#define NUM_AA_TYPE	(26)
#define NUM_Atom_TYPE	(32)

//==============Directories you might want to change=============
char workingdir[]="/home/rjloureiro/DOCKING/DOCKING_B2M-D76N-I1-PH7p2---B2M-D76N-I2-PH7p2_newcf_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new";
char writingdir[]="/home/rjloureiro/DOCKING/INTERMAPS";
char rootdir[]="/home/rjloureiro";
char ss[]="INTERMAPS";
//===============================================================

int ResNum[2];
int nAtom[2];
int PAIRTRACK[2];

double cmap1[MAXATOMNUM][MAXATOMNUM];		//atomic resolution contact map
double cmap2[MAXRESIDUE][MAXRESIDUE];		//C-alpha resolution contact map

double Cmap1[MAXATOMNUM][MAXATOMNUM];		//atomic resolution contact map
double Cmap2[MAXRESIDUE][MAXRESIDUE];		//C-alpha resolution contact map

double cMAP1[MAXATOMNUM][MAXATOMNUM];		//atomic resolution contact map
double cMAP2[MAXRESIDUE][MAXRESIDUE];		//C-alpha resolution contact map

double CMAP1[MAXATOMNUM][MAXATOMNUM];		//atomic resolution contact map
double CMAP2[MAXRESIDUE][MAXRESIDUE];		//C-alpha resolution contact map



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

	double xx;
	double yy;
	double zz;
	double r0;
	
}MYCONF, *PMYCONF;

typedef struct	{

	double r01;
	double r02;
	double r03;
	double rsq1;
	
}MYDISTANCES, *PMYDISTANCES;

double radii[]={1.61, 1.76, 1.88, 1.88, 1.88, 1.64, 1.64, 1.64, 1.64, 1.42, 1.46, 1.77, 1.77, 1.05, 0.58};

char AA_Name[NUM_AA_TYPE][5]={"METN","ILE","ILEN","VAL","LEU","PHE","CYS","MET","METC","ALA","GLY","THR","SER","TRP","TYR","PRO","HIS","HISH","GLN",
"ASN","GLU","GLUH","ASP","ASPH","LYS","ARG"};/*check*/

float atomic_solvation_parameters[6]={12.0, -116.0, -18.0, -175.0, -186.0, 0.0};

char atom_name[NUM_Atom_TYPE][4]={"CB","CG","CG1","CG2","CD","CD1","CD2","CE","CE1","CE2","CE3","CH2","CZ","CZ2","CZ3","ND1","ND2","NE","NE1","NE2",
"NH1","NH2","NZ","OD1","OD2","OE1","OE2","OG","OG1","OH","SG","SD"};

float atomic_SASA[NUM_AA_TYPE][NUM_Atom_TYPE]={ 
{73.0,  42.0, 0.0,  0.0,  0.0,  0.0,  0.0,  76.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  40.0},
{41.0,  0.0,  41.0, 72.0, 0.0,  72.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{41.0,  0.0,  41.0, 72.0, 0.0,  72.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{44.0,  0.0,  76.0, 76.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{68.0,  17.0, 0.0,  0.0,  0.0,  68.0, 68.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{73.0,  5.0,  0.0,  0.0,  0.0,  32.0, 33.0, 0.0,  39.0, 39.0, 0.0,  0.0,  39.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{91.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  79.0, 0.0},
{73.0,  42.0, 0.0,  0.0,  0.0,  0.0,  0.0,  76.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  40.0},
{73.0,  42.0, 0.0,  0.0,  0.0,  0.0,  0.0,  76.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  40.0},
{137.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{58.0,  0.0,  0.0,  82.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  40.0, 0.0,  0.0,  0.0},
{105.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  43.0, 0.0,  0.0,  0.0,  0.0},
{70.0,  6.0,  0.0,  0.0,  0.0,  43.0, 6.0,  0.0,  0.0,  8.0,  31.0, 39.0, 0.0,  38.0, 38.0, 0.0,  0.0,  0.0,  26.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{73.0,  6.0,  0.0,  0.0,  0.0,  32.0, 33.0, 0.0,  38.0, 38.0, 0.0,  0.0,  13.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  38.0, 0.0,  0.0},
{0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{77.0,  8.0,  0.0,  0.0,  0.0,  0.0,  44.0, 0.0,  53.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  25.0, 0.0,  0.0,  0.0,  28.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{77.0,  8.0,  0.0,  0.0,  0.0,  0.0,  44.0, 0.0,  53.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  25.0, 0.0,  0.0,  0.0,  28.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{74.0,  40.0, 0.0,  0.0,  16.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  56.0, 0.0,  0.0,  0.0,  0.0,  0.0,  35.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{79.0,  19.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  59.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  37.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{75.0,  41.0, 0.0,  0.0,  27.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  32.0, 37.0, 0.0,  0.0,  0.0,  0.0,  0.0},
{75.0,  41.0, 0.0,  0.0,  27.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  32.0, 37.0, 0.0,  0.0,  0.0,  0.0,  0.0},
{82.0,  29.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  37.0, 38.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{82.0,  29.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  37.0, 38.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{76.0,  36.0, 0.0,  0.0,  31.0, 0.0,  0.0,  42.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  60.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
{74.0,  36.0, 0.0,  0.0,  32.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  12.0, 0.0,  0.0,  0.0,  0.0,  15.0, 0.0,  0.0,  49.0, 59.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
};/*atomic solvent accessible surface area*/

MYATOM Atoms[MAXATOMNUM][2];

MYCHARGES Charges[LENCHARGESFILE];

MYDISTANCES dist;

int IDs[ALLPAIRS][3];

int it[MAXCONTACT];
float q0[MAXCONTACT];
float q1[MAXCONTACT];
int ct0[MAXCONTACT];
int ct1[MAXCONTACT];
int at0[MAXCONTACT];
int at1[MAXCONTACT];
int asp0[MAXCONTACT];
int asp1[MAXCONTACT];
double hidropathy_0[MAXCONTACT];
double hidropathy_1[MAXCONTACT];
double hp[MAXCONTACT];


//FUNCTIONS
double abbs(double a);
int IsANumber(char c);
int Get_AA_Type(char SzAAName[]);
int Get_Atom_Type(char SzAtomName[]);
void AssignAtomRadii(void);
int QueryAtomType(char *s, char *res);
void READPDBIDS(void);
void READPDBs(int a, int b, int c);
void READCHARGES(void);
void READCHARGES(void);


//int ALLPAIRS = 0;


int main(void)
{
	FILE *input1;
	FILE *input2;
	FILE *ccmap1;

	int i,j,k,l,m,n,o,p,r,len1,len2;
	int ReadCoord, HIS_Type;
	int filedummie=0;
	int N=2;
	
	int contactflag=0;
	int readerflag=0;
	int counter=0;
	int contactcounter,contactcounter1,contactcounter2,contactcounter3,contactcounter4;
	int pairtrackcounter=0;
	int dummie1, ID1, ID2;
	float dummie2, dummie3, dummie4;
	
	double dDA;
	double r_ik_sq,r01,r02;
	
	double total=0.0;
			
	contactcounter=0;
	contactcounter1=0;
	contactcounter2=0;
	contactcounter3=0;
	contactcounter4=0;
		
	
	READCHARGES();

	for (i=0;i<MAXATOMNUM;i++)
	{
		for (j=0;j<MAXATOMNUM;j++)
		{
			cmap1[i][j]=0.0;

		}
	}
	
	for (i=0;i<MAXRESIDUE;i++)
	{
		for (j=0;j<MAXRESIDUE;j++)
		{
			cmap2[i][j]=0.0;

		}
	}
		

	for (i=0;i<MAXATOMNUM;i++)
	{
		for (j=0;j<MAXATOMNUM;j++)
		{
			Cmap1[i][j]=0.0;

		}
	}
	
	for (i=0;i<MAXRESIDUE;i++)
	{
		for (j=0;j<MAXRESIDUE;j++)
		{
			Cmap2[i][j]=0.0;

		}
	}

	for (i=0;i<MAXATOMNUM;i++)
	{
		for (j=0;j<MAXATOMNUM;j++)
		{
			cMAP1[i][j]=0.0;

		}
	}
	
	for (i=0;i<MAXRESIDUE;i++)
	{
		for (j=0;j<MAXRESIDUE;j++)
		{
			cMAP2[i][j]=0.0;

		}
	}
		

	for (i=0;i<MAXATOMNUM;i++)
	{
		for (j=0;j<MAXATOMNUM;j++)
		{
			CMAP1[i][j]=0.0;

		}
	}
	
	for (i=0;i<MAXRESIDUE;i++)
	{
		for (j=0;j<MAXRESIDUE;j++)
		{
			CMAP2[i][j]=0.0;

		}
	}
		
	mkdir(ss,457);
	
	
	chdir(workingdir);
	
	READPDBIDS();

	//ALLPAIRS=1;

	for (i=0;i<ALLPAIRS;i++)
	{
		printf("%d	%d	%d\n",i+1,IDs[i][0],IDs[i][1]);
	}
	
	for (j=0;j<ALLPAIRS;j++)
	{
		counter=0;
				
		READPDBs(IDs[j][0],IDs[j][1],j+1);
		
		AssignAtomRadii();
		
		for (N=0; N<MAXCONTACT; N++)
		{
			it[N]=0;
		}

		for (N=0; N<MAXATOMNUM; N++)
		{
			strcpy(Atoms[N][0].AtmHBlabel,"NoHB");
			strcpy(Atoms[N][1].AtmHBlabel,"NoHB");
		}
			
		//taking the atom number i and check contact with the atom k
		
		//i-axis = vertical axis in matrix = structure 1 (aa.pdb)
		//k-axis = horizontal axis in matrix = structure 2 (bb.pdb)
				
		for (i=Atoms[0][0].AtomIndex;i<=Atoms[nAtom[0]-1][0].AtomIndex;i++)
		{					
			for (k=Atoms[0][1].AtomIndex;k<=Atoms[nAtom[1]-1][1].AtomIndex;k++)
			{
				dist.rsq1=0;

				dist.r01 = HARDCORE_R*(Atoms[i-1][0].r0+Atoms[k-1][1].r0);
				dist.r02 = CONTACT_R*(Atoms[i-1][0].r0+Atoms[k-1][1].r0);
				dist.r03 = CONTACT_R_hydropathy*(Atoms[i-1][0].r0+Atoms[k-1][1].r0);

				dist.rsq1 = (Atoms[k-1][1].x-Atoms[i-1][0].x)*(Atoms[k-1][1].x-Atoms[i-1][0].x)+(Atoms[k-1][1].y-Atoms[i-1][0].y)*(Atoms[k-1][1].y-Atoms[i-1][0].y)+(Atoms[k-1][1].z-Atoms[i-1][0].z)*(Atoms[k-1][1].z-Atoms[i-1][0].z);
				//dist.rsq2 = (Atoms[k-1][1].x-Atoms[i-1][0].x)*(Atoms[k-1][1].x-Atoms[i-1][0].x)+(Atoms[k-1][1].y-Atoms[i-1][0].y)*(Atoms[k-1][1].y-Atoms[i-1][0].y)+(Atoms[k-1][1].z-Atoms[i-1][0].z)*(Atoms[k-1][1].z-Atoms[i-1][0].z);
							
				//contactflag++;

				/*if ((j==0)&&(i==19)&&(k==72))
				{
					printf("%d	%lf	%lf	%lf\n",contactflag,(stack1[n].xx-stack2[o].xx)*(stack1[n].xx-stack2[o].xx),(stack1[n].yy-stack2[o].yy)*(stack1[n].yy-stack2[o].yy),stack1[n].zz-stack2[o].zz)*(stack1[n].zz-stack2[o].zz);
				}*/
							
				/*if ((j==0)&&(i==72)&&(k==19))
				{
					printf("%lf	%lf	%lf\n",dist.r01,sqrt(dist.rsq1),dist.r02);
				}*/
				
				if ((dist.rsq1>=dist.r01*dist.r01)&&(dist.rsq1<=dist.r02*dist.r02))
				{
					for (p = 0; p < LENCHARGESFILE; p++)
					{
						if ((strcmp(Atoms[i-1][0].AtmName,Charges[p].AtmName) == 0) && (Atoms[i-1][0].Res_Type == Charges[p].Res_Type))
						{
							for (r = 0; r < LENCHARGESFILE; r++)
							{
								if ((strcmp(Atoms[k-1][1].AtmName,Charges[r].AtmName) == 0) && (Atoms[k-1][1].Res_Type == Charges[r].Res_Type))
								{
									q0[counter] = Charges[p].Atm_Charge;
									q1[counter] = Charges[r].Atm_Charge;
								}
							}
						}
					}
					if (q0[counter] < 0.0 & q1[counter] < 0.0)
					{	
						it[counter] = 1; /*possible electrostatic interaction*/
						if (((strcmp(Atoms[i-1][0].AtmName,"N") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H1") != 0) &&
						(strcmp(Atoms[i-1][0].AtmName,"H2") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H3") != 0) && (strcmp(Atoms[i-1][0].AtmName,"CA") != 0) &&
						(strcmp(Atoms[i-1][0].AtmName,"C") != 0) && (strcmp(Atoms[i-1][0].AtmName,"O") != 0) && (strcmp(Atoms[i-1][0].AtmName,"O1") != 0) &&
						(strcmp(Atoms[i-1][0].AtmName,"O2") != 0)) && ((strcmp(Atoms[k-1][1].AtmName,"N") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H1") != 0) &&
						(strcmp(Atoms[k-1][1].AtmName,"H2") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H3") != 0) && (strcmp(Atoms[k-1][1].AtmName,"CA") != 0) &&
						(strcmp(Atoms[k-1][1].AtmName,"C") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O1") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O2") != 0)))
						{
							if ((strncmp(Atoms[i-1][0].AtmName,"H", 1) != 0) && (strncmp(Atoms[k-1][1].AtmName,"H", 1) != 0) && (dist.rsq1<=dist.r03*dist.r03))
							{
								it[counter] = 3;
							}
						}
					}
					else if (q0[counter] > 0.0 & q1[counter] > 0.0)
					{
						it[counter] = 1; /*possible electrostatic interaction*/
						if (((strcmp(Atoms[i-1][0].AtmName,"N") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H1") != 0) &&
						(strcmp(Atoms[i-1][0].AtmName,"H2") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H3") != 0) && (strcmp(Atoms[i-1][0].AtmName,"CA") != 0) &&
						(strcmp(Atoms[i-1][0].AtmName,"C") != 0) && (strcmp(Atoms[i-1][0].AtmName,"O") != 0) && (strcmp(Atoms[i-1][0].AtmName,"O1") != 0) &&
						(strcmp(Atoms[i-1][0].AtmName,"O2") != 0)) && ((strcmp(Atoms[k-1][1].AtmName,"N") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H1") != 0) &&
						(strcmp(Atoms[k-1][1].AtmName,"H2") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H3") != 0) && (strcmp(Atoms[k-1][1].AtmName,"CA") != 0) &&
						(strcmp(Atoms[k-1][1].AtmName,"C") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O1") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O2") != 0)))
						{
							if ((strncmp(Atoms[i-1][0].AtmName,"H", 1) != 0) && (strncmp(Atoms[k-1][1].AtmName,"H", 1) != 0) && (dist.rsq1<=dist.r03*dist.r03))
							{
								it[counter] = 3;
							}
						}						
					}												
					else if (q0[counter] < 0.0 & q1[counter] > 0.0)
					{
						it[counter] = 1; /*possible electrostatic interaction*/
						if (((strncmp(Atoms[i-1][0].AtmName,"O",1) == 0) || (strncmp(Atoms[i-1][0].AtmName,"N",1) == 0) || (strncmp(Atoms[i-1][0].AtmName,"S",1) == 0)) && (strncmp(Atoms[k-1][1].AtmName,"H",1) == 0))
						{
							if ((strcmp(Atoms[k-1][1].AtmName,"H2") == 0) || (strcmp(Atoms[k-1][1].AtmName,"HH12") == 0) || (strcmp(Atoms[k-1][1].AtmName,"HH22") == 0) || (strcmp(Atoms[k-1][1].AtmName,"HD22") == 0) || (strcmp(Atoms[k-1][1].AtmName,"HE22") == 0) || (strcmp(Atoms[k-1][1].AtmName,"HZ2") == 0))
							{
								if (strncmp(Atoms[k-3][1].AtmName,"N",1) == 0)
								{	
									if ((strcmp(Atoms[i-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Atoms[k-3][1].AtmHBlabel,"HB") != 0))
									{
										dDA = sqrt((Atoms[k-3][1].x - Atoms[i-1][0].x)*(Atoms[k-3][1].x - Atoms[i-1][0].x) + (Atoms[k-3][1].y -
										Atoms[i-1][0].y)*(Atoms[k-3][1].y - Atoms[i-1][0].y) + (Atoms[k-3][1].z - Atoms[i-1][0].z)*(Atoms[k-3][1].z -
										Atoms[i-1][0].z));					
										if (dDA<=3.2)
										{
											it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
											strcpy(Atoms[i-1][0].AtmHBlabel,"HB");
											strcpy(Atoms[k-3][1].AtmHBlabel,"HB");
										}
									}	
								}
							}
							else if ((strcmp(Atoms[k-1][1].AtmName,"H3") == 0) || (strcmp(Atoms[k-1][1].AtmName,"HZ3") == 0))
							{
								if (strncmp(Atoms[k-4][1].AtmName,"N",1) == 0)
								{
									if ((strcmp(Atoms[i-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Atoms[k-4][1].AtmHBlabel,"HB") != 0))
									{
										dDA = sqrt((Atoms[k-4][1].x - Atoms[i-1][0].x)*(Atoms[k-4][1].x - Atoms[i-1][0].x) + (Atoms[k-4][1].y -
										Atoms[i-1][0].y)*(Atoms[k-4][1].y - Atoms[i-1][0].y) + (Atoms[k-4][1].z - Atoms[i-1][0].z)*(Atoms[k-4][1].z -
										Atoms[i-1][0].z));					
										if (dDA<=3.2)
										{
											it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
											strcpy(Atoms[i-1][0].AtmHBlabel,"HB");
											strcpy(Atoms[k-4][1].AtmHBlabel,"HB");
										}
									}
								}	
							}
							else
							{
								if ((strncmp(Atoms[k-2][1].AtmName,"N",1) == 0) || (strncmp(Atoms[k-2][1].AtmName,"O",1) == 0))
								{
									if ((strcmp(Atoms[i-1][0].AtmHBlabel,"HB") != 0) && (strcmp(Atoms[k-2][1].AtmHBlabel,"HB") != 0))
									{
										dDA = sqrt((Atoms[k-2][1].x - Atoms[i-1][0].x)*(Atoms[k-2][1].x - Atoms[i-1][0].x) + (Atoms[k-2][1].y -
										Atoms[i-1][0].y)*(Atoms[k-2][1].y - Atoms[i-1][0].y) + (Atoms[k-2][1].z - Atoms[i-1][0].z)*(Atoms[k-2][1].z -
										Atoms[i-1][0].z));					
										if (dDA<=3.2)
										{
											it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
											strcpy(Atoms[i-1][0].AtmHBlabel,"HB");
											strcpy(Atoms[k-2][1].AtmHBlabel,"HB");
										}
									}
								}
							}									
						}
						if (((strcmp(Atoms[i-1][0].AtmName,"N") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H1") != 0) &&
						(strcmp(Atoms[i-1][0].AtmName,"H2") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H3") != 0) && (strcmp(Atoms[i-1][0].AtmName,"CA") != 0) &&
						(strcmp(Atoms[i-1][0].AtmName,"C") != 0) && (strcmp(Atoms[i-1][0].AtmName,"O") != 0) && (strcmp(Atoms[i-1][0].AtmName,"O1") != 0) &&
						(strcmp(Atoms[i-1][0].AtmName,"O2") != 0)) && ((strcmp(Atoms[k-1][1].AtmName,"N") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H1") != 0) &&
						(strcmp(Atoms[k-1][1].AtmName,"H2") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H3") != 0) && (strcmp(Atoms[k-1][1].AtmName,"CA") != 0) &&
						(strcmp(Atoms[k-1][1].AtmName,"C") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O1") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O2") != 0)))
						{
							if ((strncmp(Atoms[i-1][0].AtmName,"H", 1) != 0) && (strncmp(Atoms[k-1][1].AtmName,"H", 1) != 0) && (dist.rsq1<=dist.r03*dist.r03))
							{
								it[counter] = 3;
							}
						}												
					}
					else if (q0[counter] > 0.0 & q1[counter] < 0.0)
					{
						it[counter] = 1; /*possible electrostatic interaction*/
						if ((strncmp(Atoms[i-1][0].AtmName,"H",1) == 0) && ((strncmp(Atoms[k-1][1].AtmName,"O",1) == 0) || (strncmp(Atoms[k-1][1].AtmName,"N",1) == 0) || (strncmp(Atoms[k-1][1].AtmName,"S",1) == 0))) 
						{
							if ((strcmp(Atoms[i-1][0].AtmName,"H2") == 0) || (strcmp(Atoms[i-1][0].AtmName,"HH12") == 0) || (strcmp(Atoms[i-1][0].AtmName,"HH22") == 0) || (strcmp(Atoms[i-1][0].AtmName,"HD22") == 0) || (strcmp(Atoms[i-1][0].AtmName,"HE22") == 0) || (strcmp(Atoms[i-1][0].AtmName,"HZ2") == 0))
							{
								if (strncmp(Atoms[i-3][0].AtmName,"N",1) == 0)
								{	
									if ((strcmp(Atoms[i-3][0].AtmHBlabel,"HB") != 0) && (strcmp(Atoms[k-1][1].AtmHBlabel,"HB") != 0))
									{
										dDA = sqrt((Atoms[k-1][1].x - Atoms[i-3][0].x)*(Atoms[k-1][1].x - Atoms[i-3][0].x) + (Atoms[k-1][1].y -
										Atoms[i-3][0].y)*(Atoms[k-1][1].y - Atoms[i-3][0].y) + (Atoms[k-1][1].z - Atoms[i-3][0].z)*(Atoms[k-1][1].z -
										Atoms[i-3][0].z));
										if (dDA<=3.2)
										{
											it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
											strcpy(Atoms[i-3][0].AtmHBlabel,"HB");
											strcpy(Atoms[k-1][1].AtmHBlabel,"HB");
										}
									}
								}
							}
							else if ((strcmp(Atoms[i-1][0].AtmName,"H3") == 0) || (strcmp(Atoms[i-1][0].AtmName,"HZ3") == 0))
							{	
								if (strncmp(Atoms[i-4][0].AtmName,"N",1) == 0)
								{	
									if ((strcmp(Atoms[i-4][0].AtmHBlabel,"HB") != 0) && (strcmp(Atoms[k-1][1].AtmHBlabel,"HB") != 0))
									{
										dDA = sqrt((Atoms[k-1][1].x - Atoms[i-4][0].x)*(Atoms[k-1][1].x - Atoms[i-4][0].x) + (Atoms[k-1][1].y -
										Atoms[i-4][0].y)*(Atoms[k-1][1].y - Atoms[i-4][0].y) + (Atoms[k-1][1].z - Atoms[i-4][0].z)*(Atoms[k-1][1].z -
										Atoms[i-4][0].z));
										if (dDA<=3.2)
										{
											it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
											strcpy(Atoms[i-4][0].AtmHBlabel,"HB");
											strcpy(Atoms[k-1][1].AtmHBlabel,"HB");
										}
									}
								}
							}
							else
							{
								if ((strncmp(Atoms[i-2][0].AtmName,"N",1) == 0) || (strncmp(Atoms[i-2][0].AtmName,"O",1) == 0))
								{
									if ((strcmp(Atoms[i-2][0].AtmHBlabel,"HB") != 0) && (strcmp(Atoms[k-1][1].AtmHBlabel,"HB") != 0))
									{
										dDA = sqrt((Atoms[k-1][1].x - Atoms[i-2][0].x)*(Atoms[k-1][1].x - Atoms[i-2][0].x) + (Atoms[k-1][1].y -
										Atoms[i-2][0].y)*(Atoms[k-1][1].y - Atoms[i-2][0].y) + (Atoms[k-1][1].z - Atoms[i-2][0].z)*(Atoms[k-1][1].z -
										Atoms[i-2][0].z));
										if (dDA<=3.2)
										{
											it[counter] = 4; /*possible electrostatic and hydrogen bond interaction*/	
											strcpy(Atoms[i-2][0].AtmHBlabel,"HB");
											strcpy(Atoms[k-1][1].AtmHBlabel,"HB"); 
										}
									}
								}
							}		
						}
						if (((strcmp(Atoms[i-1][0].AtmName,"N") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H1") != 0) &&
						(strcmp(Atoms[i-1][0].AtmName,"H2") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H3") != 0) && (strcmp(Atoms[i-1][0].AtmName,"CA") != 0) &&
						(strcmp(Atoms[i-1][0].AtmName,"C") != 0) && (strcmp(Atoms[i-1][0].AtmName,"O") != 0) && (strcmp(Atoms[i-1][0].AtmName,"O1") != 0) &&
						(strcmp(Atoms[i-1][0].AtmName,"O2") != 0)) && ((strcmp(Atoms[k-1][1].AtmName,"N") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H1") != 0) &&
						(strcmp(Atoms[k-1][1].AtmName,"H2") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H3") != 0) && (strcmp(Atoms[k-1][1].AtmName,"CA") != 0) &&
						(strcmp(Atoms[k-1][1].AtmName,"C") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O1") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O2") != 0)))
						{
							if ((strncmp(Atoms[i-1][0].AtmName,"H", 1) != 0) && (strncmp(Atoms[k-1][1].AtmName,"H", 1) != 0) && (dist.rsq1<=dist.r03*dist.r03))
							{
								it[counter] = 3;
							}
						}												
					}
					else if (((strcmp(Atoms[i-1][0].AtmName,"N") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H1") != 0) &&
					(strcmp(Atoms[i-1][0].AtmName,"H2") != 0) && (strcmp(Atoms[i-1][0].AtmName,"H3") != 0) && (strcmp(Atoms[i-1][0].AtmName,"CA") != 0) &&
					(strcmp(Atoms[i-1][0].AtmName,"C") != 0) && (strcmp(Atoms[i-1][0].AtmName,"O") != 0) && (strcmp(Atoms[i-1][0].AtmName,"O1") != 0) &&
					(strcmp(Atoms[i-1][0].AtmName,"O2") != 0)) && ((strcmp(Atoms[k-1][1].AtmName,"N") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H1") != 0) &&
					(strcmp(Atoms[k-1][1].AtmName,"H2") != 0) && (strcmp(Atoms[k-1][1].AtmName,"H3") != 0) && (strcmp(Atoms[k-1][1].AtmName,"CA") != 0) &&
					(strcmp(Atoms[k-1][1].AtmName,"C") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O1") != 0) && (strcmp(Atoms[k-1][1].AtmName,"O2") != 0)))
					{
						if ((strncmp(Atoms[i-1][0].AtmName,"H", 1) != 0) && (strncmp(Atoms[k-1][1].AtmName,"H", 1) != 0) && (dist.rsq1<=dist.r03*dist.r03))
						{
							it[counter] = 2;
						}
					}												
					
					if ((it[counter] == 1) || (it[counter] == 3) || (it[counter] == 4))
					{
						cmap1[i][k]=cmap1[i][k]+1.0;
						//contactflag=1;
						//break;

						contactcounter++;	
					}
					else if (it[counter] == 2)
					{
						cMAP1[i][k]=cMAP1[i][k]+1.0;
						//contactflag=1;
						//break;

						contactcounter++;	
					}
					/*else if (it[counter] == 3)
					{
						cmap1[i][k]=cmap1[i][k]+1.0;
						cMAP1[i][k]=cMAP1[i][k]+1.0;
						CMAP1[i][k]=CMAP1[i][k]+1.0;
						//contactflag=1;
						//break;

						contactcounter1++;
						contactcounter2++;
						contactcounter4++;	
					}					
					else if (it[counter] == 4)
					{
						cmap1[i][k]=cmap1[i][k]+1.0;
						Cmap1[i][k]=Cmap1[i][k]+1.0;
						//contactflag=1;
						//break;

						contactcounter1++;
						contactcounter3++;	
					}*/
					counter++;					
				}
				/*if ((j==0)&&(i==72)&&(k==19))
				{
					printf("%d	%lf	%lf	%lf	%lf	%d	%d\n",contactflag,dist.r01,sqrt(dist.rsq1),sqrt(dist.rsq2),dist.r02,len1,len2);
				}
				if ((j==0)&&(i==19)&&(k==72))
				{
					printf("%d	%lf	%lf	%lf	%lf	%d	%d\n",contactflag,dist.r01,sqrt(dist.rsq1),sqrt(dist.rsq2),dist.r02,len1,len2);
				}*/

				//contactflag=0;		
			}		
		}			
					
	}
	
	for (i=Atoms[0][0].AtomIndex;i<=Atoms[nAtom[0]-1][0].AtomIndex;i++)
	{
		for (k=Atoms[0][1].AtomIndex;k<=Atoms[nAtom[1]-1][1].AtomIndex;k++)
		{
			cmap1[i][k]=cmap1[i][k]/(double)contactcounter;
			cMAP1[i][k]=cMAP1[i][k]/(double)contactcounter;
			//Cmap1[i][k]=Cmap1[i][k]/(double)contactcounter3;
			//CMAP1[i][k]=CMAP1[i][k]/(double)contactcounter4;
		}
	}

	for (i=Atoms[0][0].NodeIndex;i<=Atoms[nAtom[0]-1][0].NodeIndex;i++)
	{
		for (j=Atoms[0][1].NodeIndex;j<=Atoms[nAtom[1]-1][1].NodeIndex;j++)
		{
			for (k=Atoms[0][0].AtomIndex;k<=Atoms[nAtom[0]-1][0].AtomIndex;k++)
			{
				for (l=Atoms[0][1].AtomIndex;l<=Atoms[nAtom[1]-1][1].AtomIndex;l++)
				{
					if (i==Atoms[k-1][0].NodeIndex && j==Atoms[l-1][1].NodeIndex)
					{
						cmap2[i][j]=cmap2[i][j]+cmap1[k][l];
						cMAP2[i][j]=cMAP2[i][j]+cMAP1[k][l];
						//Cmap2[i][j]=Cmap2[i][j]+Cmap1[k][l];
						//CMAP2[i][j]=CMAP2[i][j]+CMAP1[k][l];
					}
				}
			}
		}
	}
					
					
	chdir(writingdir);
	ccmap1=fopen("INTER-CMAP-B2M_D76N-I1-pH7p2---B2M-D76N-I2-pH7p2_2000_T-0.15_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new.dat","w");

	for (i=Atoms[0][0].NodeIndex;i<=Atoms[nAtom[0]-1][0].NodeIndex;i++)
	{
		for (j=Atoms[0][1].NodeIndex;j<=Atoms[nAtom[1]-1][1].NodeIndex;j++)
		{
			fprintf(ccmap1,"%lf	",(double)cmap2[i][j]+cMAP2[i][j]);
			total=total+cmap2[i][j]+cMAP2[i][j];
		}
		fprintf(ccmap1,"\n");
	}
	printf("%lf", total);
	//printf("%d	%d\n",Atoms[0][0].NodeIndex, Atoms[nAtom[0]-1][0].NodeIndex);
	
	
	fclose(ccmap1);
		
	return 0;
	
}

void READPDBs(int a, int b, int c)
{
	FILE *input;
	
	int i,j,k,m;
	int ReadCoord, HIS_Type;
	
	int cter;

	char s1[256], s2[256], szLine[MAXLENLINE], szAtomName[3], szResName[4], buff[256],*ReadStatus;
	
	for (i=0;i<2;i++)
	{
		input=NULL;

		ResNum[i]=0;
		nAtom[i]=0;
		
		if (i==0)
		{
			
			for (m=1;m<=3;m++)
			{
				sprintf(s1, "D76N_I1_D76N_I2_PAIRS_ID1=%d_ID2=%d_DIR=%d_newcf__a.pdb", a,b,m);
				sprintf(s2, "D76N_I1_D76N_I2_PAIRS_ID1=%d_ID2=%d_DIR=%d_newcf__a.pdb", a,b,-m);
			
				input=fopen(s1,"r");
				
				if (input==NULL)
				{
					input=fopen(s2,"r");
					if(input != NULL)	
					{
						printf("%d	Filename: [%s]\n", c, s2); fflush(stdout);
						break;
					}
					
				}
				else
				{
					printf("%d	Filename: [%s]\n", c, s1); fflush(stdout);
					break;
				}
			}
		}
		else
		{
			for (m=1;m<=3;m++)
			{
				sprintf(s1, "D76N_I1_D76N_I2_PAIRS_ID1=%d_ID2=%d_DIR=%d_newcf__b.pdb", a,b,m);
				sprintf(s2, "D76N_I1_D76N_I2_PAIRS_ID1=%d_ID2=%d_DIR=%d_newcf__b.pdb", a,b,-m);
				
				input=fopen(s1,"r");
				
				if (input==NULL)
				{
					input=fopen(s2,"r");
					if(input != NULL)	
					{
						printf("%d	Filename: [%s]\n", c, s2); fflush(stdout);
						break;
					}
				}
				else
				{
					printf("%d	Filename: [%s]\n", c, s1); fflush(stdout);
					break;
				}
			}
		}
		
		
		
		if (input==NULL)
		{
			printf("Failed to open PDB pair %d + %d. Shutting down", a,b);
			exit(0);
		}	
		
		//printf("Filename: [%s]\n", s); fflush(stdout);
		
		///////////////////////////////////////////////////////////////////////////////
		
	
		
	
		k=0;
		while(1)
		{
			
			ReadStatus = fgets(szLine, MAXLENLINE, input);
			
			if(ReadStatus == NULL)	
			{
				break;
			}
			
										
			if (strncmp(szLine, "ATOM", 4)==0)
			{
				sscanf(szLine+5, "%d", &(Atoms[nAtom[i]][i].AtomIndex));
				sscanf(szLine+12, "%s", szAtomName);
				sscanf(szLine+22, "%d", &(Atoms[nAtom[i]][i].NodeIndex));
			
				if(szLine[16] == ' ')	
				{
					sscanf(szLine+16, "%s", szResName);
				}
				else if(szLine[16] == 'A')	//only read atom in chain A
				{	
					sscanf(szLine+17, "%s", szResName);
				}
				else
				{
					continue;
				}
						
			
				//start	to handle special cases
					
				if(strncmp(szResName, "HIS", 3)==0)
				{
					HIS_Type = HIS_TYPE_HIS;
				}
				if(strncmp(szResName, "HSD", 3)==0)
				{
					HIS_Type = HIS_TYPE_HSD;
					strcpy(szResName, "HIS");
				}
				if(strncmp(szResName, "HSC", 3)==0)
				{
					HIS_Type = HIS_TYPE_HSC;
					strcpy(szResName, "HIS");
				}
				
				//IMPORTANT:Consider the PROTONATION states of histidines at each pH/system and change histidines name accordingly
				
				if (i==0)
				{
					/*if((strncmp(szResName, "HIS", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 13))
					{
						strcpy(szResName, "HISH");
					}
					if((strncmp(szResName, "HIS", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 31))
					{
						strcpy(szResName, "HISH");
					}
					if((strncmp(szResName, "HIS", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 51))
					{
						strcpy(szResName, "HISH");
					}
					if((strncmp(szResName, "HIS", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 84))
					{
						strcpy(szResName, "HISH");
					}*/
				}
				else
				{
					/*if((strncmp(szResName, "HIS", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 13))
					{
						strcpy(szResName, "HISH");
					}
					if((strncmp(szResName, "HIS", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 31))
					{
						strcpy(szResName, "HISH");
					}
					if((strncmp(szResName, "HIS", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 51))
					{
						strcpy(szResName, "HISH");
					}
					if((strncmp(szResName, "HIS", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 84))
					{
						strcpy(szResName, "HISH");
					}*/
				}
				if((strncmp(szResName, "ILE", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 1))
				{
					strcpy(szResName, "ILEN");
				}
				if((strncmp(szResName, "MET", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 1))
				{
					strcpy(szResName, "METN");
				}				
				if((strncmp(szResName, "MET", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 94))
				{
					strcpy(szResName, "METC");
				}
				if((strncmp(szResName, "MET", 3)==0) && (Atoms[nAtom[i]][i].NodeIndex == 99))
				{
					strcpy(szResName, "METC");
				}								
					
				if((strncmp(szAtomName, "CD", 2)==0) && (strncmp(szResName, "ILE", 3)==0))
				{
					strcpy(szAtomName, "CD1");
				}
					//end	to handle special cases
				
				
				strcpy(Atoms[nAtom[i]][i].AtmName, szAtomName);
				if((strcmp(szAtomName, "N")!=0) && (strncmp(szAtomName, "H", 1)!=0) && (strcmp(szAtomName, "CA")!=0) && (strcmp(szAtomName, "C")!=0) && (strcmp(szAtomName, "O")!=0) && (strcmp(szAtomName, "O1")!=0) && (strcmp(szAtomName, "O2")!=0))
				{
					Atoms[nAtom[i]][i].Atom_Type = Get_Atom_Type(szAtomName);
				}
				sscanf(szLine+22, "%d", &(Atoms[nAtom[i]][i].NodeIndex));
				Atoms[nAtom[i]][i].Res_Type = Get_AA_Type(szResName);
					
				ReadCoord = sscanf(szLine+30, "%lf %lf %lf", &(Atoms[nAtom[i]][i].x), &(Atoms[nAtom[i]][i].y), &(Atoms[nAtom[i]][i].z));
			
				if(ReadCoord != 3)	
				{
					printf("Fail to read coordinates from PDB file:ID pair %d / %d Atom %d AtomName = %s ResName = %s\nQuit\n", 
					a,b, nAtom+1, szAtomName, szResName);
					fclose(input);
					exit(1);
				}
			
				nAtom[i]++;
					
			}
					
			ResNum[i] = Atoms[nAtom[i]-1][i].NodeIndex - Atoms[0][i].NodeIndex + 1;
			k++;
												
		}
		fclose(input);
	}
	
	//printf("READING COMPLETE\n\n");
	

}

//pdb-compatible screen output
//for (i=0;i<nAtom;i++)
//{
//	printf("ATOM%7d	%-3s	%-3s	A%4d	%8.3lf	%8.3lf	%8.3lf	1.00 0.00\n",
//	i+1, Atoms[i].AtmName, AA_Name[Atoms[i].Res_Type],Atoms[i].NodeIndex, Atoms[i].x,Atoms[i].y,Atoms[i].z);
//}

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
	int n, Type;
	
	for (n=0;n<nAtom[0];n++)
	{
		Type = QueryAtomType(Atoms[n][0].AtmName, AA_Name[Atoms[n][0].Res_Type]);
		Atoms[n][0].r0=radii[Type];
	}
	
	for (n=0;n<nAtom[1];n++)
	{
		Type = QueryAtomType(Atoms[n][1].AtmName, AA_Name[Atoms[n][1].Res_Type]);
		Atoms[n][1].r0=radii[Type];
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

void READCHARGES(void)
{
	FILE *input;
	
	int i = 0;
	
	char s[256], szLine[MAXLENLINE], szAtomName[8], szResName[8], szAtomType[8], buff[256],*ReadStatus;
	
	float AtomCharge;
	
	sprintf(s, "charges.dat");
	
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

void READPDBIDS(void)
{
	FILE *handler;
	
	int i,dummie1,ID1,ID2;
	float dummie2,dummie3,dummie4,dummie5,dummie6,dummie7,dummie8,dummie9,dummie10,dummie11,dummie12,dummie13,dummie14;
	
	//dummie4 = contacts	
	
	i=0;
	handler=fopen("B2M_D76N-I1-pH7p2---B2M-D76N-I2-pH7p2_2000_T-0.15_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new.dat","r");
	
	while (!feof(handler)) 
	{
		if (fscanf(handler, "%d	%d	%d	%f	%f	%f	%f	%f	%f	%f	%f      %f      %f", &ID1, &ID2, &dummie1, &dummie2, &dummie3, &dummie4, &dummie5, &dummie6, &dummie7, &dummie8, &dummie9, &dummie10, &dummie11) > 0)
		{
    			if (1)
			{
				IDs[i][0]=ID1;
				IDs[i][1]=ID2;
				IDs[i][2]=(int)dummie8;
				i++;
				//ALLPAIRS++
			}
    		}
	};
	fclose(handler);
	
}



