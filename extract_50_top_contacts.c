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
char workingdir[] = "/home/rjloureiro/DOCKING/INTERMAPS";
char writingdir[] = "/home/rjloureiro/DOCKING/hotspots";
char rootdir[] = "/home/rjloureiro";
char ss[] = "INTERMAPS";
//===============================================================


double cmap[MAXRESIDUE][MAXRESIDUE];

double max_cmap[N_MAX_CMAP];

int index_i[N_MAX_CMAP];

int index_j[N_MAX_CMAP];

double max_cmap_50[ALLPAIRS];

int aa_pairs_a[ALLPAIRS];

int aa_pairs_b[ALLPAIRS];

void READCMAP(void);


int main(void)
{
	FILE *output;
	FILE *contactsonly;
	
	int i, j, a, b, k, n, o, p;
	int counterflag;
	double max;
	
	chdir(workingdir);
	
	READCMAP();
	
	for (k = 0; k < N_MAX_CMAP; k++)
	{
		max_cmap[k] = 0.0;
	}
	
	for (k = 0; k < N_MAX_CMAP; k++)
	{
		index_i[k] = 0;
	}
	
	for (k = 0; k < N_MAX_CMAP; k++)
	{
		index_j[k] = 0;
	}
	
	for (k = 0; k < ALLPAIRS; k++)
	{
		max_cmap_50[k] = 0.0;
	}
	
	for (k = 0; k < ALLPAIRS; k++)
	{
		aa_pairs_a[k] = 0;
	}
	
	for (k = 0; k < ALLPAIRS; k++)
	{
		aa_pairs_b[k] = 0;
	}
	
	
	for (p = 0; p < N_MAX_CMAP; p++)
	{
		max = -1.0;
		for (i = 0; i < MAXRESIDUE; i++)
		{
			for (j = 0; j < MAXRESIDUE; j++)
			{	
				if (p == 0)
				{
					if (cmap[i][j] >= max)
					{
						max = cmap[i][j];
						max_cmap[p] = cmap[i][j];
						index_i[p] = i;
						index_j[p] = j;
					}
				}
				else
				{
					if ((cmap[i][j] >= max) & (cmap[i][j] <= max_cmap[p - 1]))
					{
						counterflag = 0;
						for (n = 1; n < ALLPAIRS; n++)
						{
							if (p-n >= 0)
							{
								if ((i == index_i[p-n]) & (j == index_j[p-n]))
								{
									counterflag = 1;
									break;
								}
							}	
						}
						if (counterflag == 0)
						{
							max = cmap[i][j];
							max_cmap[p] = cmap[i][j];
							index_i[p] = i;
							index_j[p] = j;	
						}
					}
				}
			}	
		}
	}
	
	
	o = 0;
	a = 0;
	b = 0;
	for (p = 0; p < ALLPAIRS; p++)
	{
		if (o < ALLPAIRS)
		{
			max_cmap_50[p] = max_cmap[o];
			aa_pairs_a[a] = index_i[o];
			aa_pairs_b[b] = index_j[o];
			o++;
			a++;
			b++;
		}
	}
	
	
	chdir(writingdir);
	
	output = fopen("top_50_contacts_B2M_D76N-I1-pH7p2---B2M-D76N-I2-pH7p2_2000_T-0.15_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new.dat", "w");
	
	//IMPORTANT: CHANGE + 1 to +6 from D76N to DN6 FOR CORRECT NUMBERING of the residues//
	
	for (p = 0; p < ALLPAIRS; p++)
	{
		fprintf(output, "%d-%d: %lf      \n", aa_pairs_a[p] + 1, aa_pairs_b[p] + 1, max_cmap_50[p]);
	}
	
	contactsonly = fopen("top_50_contactsonly_B2M_D76N-I1-pH7p2---B2M-D76N-I2-pH7p2_2000_T-0.15_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new.dat", "w");
	
	for (p = 0; p < ALLPAIRS; p++)
	{
		fprintf(contactsonly, "%d-%d\n", aa_pairs_a[p] + 1, aa_pairs_b[p] + 1);
	}
	
	fclose(output);
	
	fclose(contactsonly);
	
	return 0;
		
}


void READCMAP(void)
{
	FILE *handler;
	
	int i, j;
	double contact_p;
	
	handler = fopen("INTER-CMAP-B2M_D76N-I1-pH7p2---B2M-D76N-I2-pH7p2_2000_T-0.15_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new.dat","r");
	
	while (!feof(handler))
	{
		for (i = 0; i<MAXRESIDUE; i++)
		{
			for (j = 0; j<MAXRESIDUE; j++)
			{
	
				fscanf(handler, "%lf	  ", &contact_p);
				cmap[i][j] = contact_p;	
				printf ("%lf\n", cmap[i][j]);	
			}
		}
	}
}


		
			

