#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

//=============Parameters you might want to change===============

//ENSSTRUCTS = number of conformations in intermediate ensembles
//D76N-NAT: 25498
//D76N-I1: 514
//D76N-I2: 2496
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

#define ALLPAIRS	(50)
#define PAIR            (2)

#define HARDCORE_R	(1.)
#define CONTACT_R	(2.33)
#define LOCAL_RANGE	(0)

//end to parameters you might want to change=====================

#define MAXLENLINE 	 (256)
#define MAXATOMNUM	 (1085)
#define MAXRESIDUE	 (99)
#define N_MAX_CMAP       (9801)
#define NUMBER_INTERVALS (10000)
#define INTERVAL_WIDTH   (10)
#define HIS_TYPE_HIS	 (0)
#define HIS_TYPE_HSD	 (1)
#define HIS_TYPE_HSC	 (2)
#define NUM_AA_TYPE	 (20)

//==============Directories you might want to change=============
char workingdir[] = "/home/rjloureiro/DOCKING/hotspots";
char writingdir[] = "/home/rjloureiro/DOCKING/hotspots";
char rootdir[] = "/home/rjloureiro";
char ss[] = "INTERMAPS";
//===============================================================


int top_50_contacts[ALLPAIRS][PAIR];

int hotspots_counter[MAXRESIDUE];

double hotspots_counter_normalized[MAXRESIDUE];

void READTOP50CONTACTS(void);


int main(void)
{
	FILE *output;
	
	int i, j, k;
	
	chdir(workingdir);
	
	READTOP50CONTACTS();
	
	for (k = 0; k < MAXRESIDUE; k++)
	{
		hotspots_counter[k] = 0;
	}
	
	
	for (k = 0; k < MAXRESIDUE; k++)
	{
		for (i = 0; i < ALLPAIRS; i++)
		{
			for (j = 0; j < PAIR; j++)
			{
				if(top_50_contacts[i][j] == (k+1))
				{
					hotspots_counter[k]++;
				}
			}
		}
	}
	
	for (k = 0; k < MAXRESIDUE; k++)
	{
		hotspots_counter_normalized[k] = hotspots_counter[k]/(double)ALLPAIRS;
	}
	
	chdir(writingdir);
	
	output = fopen("hotspots_top_50_contacts_B2M_D76N-I1-pH7p2---B2M-D76N-I2-pH7p2_2000_T-0.15_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new.dat", "w");
	
	for (k = 0; k < MAXRESIDUE; k++)
	{
		fprintf(output, "%d       %lf        \n", k+1, hotspots_counter_normalized[k]);
	}
	
	fclose(output);
	
	return 0;
	
}


void READTOP50CONTACTS(void)
{
	FILE *handler;
	
	int i, j, aa;
	
	handler = fopen("top_50_contactsonly_B2M_D76N-I1-pH7p2---B2M-D76N-I2-pH7p2_2000_T-0.15_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new.dat", "r");
	
	while (!feof(handler))
	{
		for (i = 0; i < ALLPAIRS; i++)
		{
			for (j = 0; j < PAIR; j++)
			{
				fscanf(handler, "%d ", & aa);
				if (aa < 0)
				{
					aa = -aa;
				}
				top_50_contacts[i][j] = aa;
				printf ("%d\n", top_50_contacts[i][j]);
			}
		}
	}
}		
	
	
	
	
	
	
			
		
