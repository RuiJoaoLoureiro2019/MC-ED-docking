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
//D76N-NAT(DMD Go model): 1399

#define ALLPAIRS	(2000)

#define HARDCORE_R	(0.80)
#define CONTACT_R	(2.33)
#define LOCAL_RANGE	(0)

//end to parameters you might want to change=====================

#define MAXLENLINE 	 (256)
#define MAXATOMNUM	 (1034)
#define MAXRESIDUE	 (99)
#define NUMBER_INTERVALS (10000)
#define INTERVAL_WIDTH   (1.0)
#define HIS_TYPE_HIS	 (0)
#define HIS_TYPE_HSD	 (1)
#define HIS_TYPE_HSC	 (2)
#define NUM_AA_TYPE	 (25)

//==============Directories you might want to change=============
char workingdir[] = "/home/rjloureiro/DOCKING/DOCKING_B2M_D76N_I2_PH7p2_newcf_n-clashes_0.1_energy_electrostatics_atticus";
char writingdir[] = "/home/rjloureiro/DOCKING";
char rootdir[] = "/home/rjloureiro";
char ss[] = "INTERMAPS";
//===============================================================


double ContactsEnergy[ALLPAIRS];

double contactsenergyinterval[NUMBER_INTERVALS];

int n_contactsenergyperinterval[NUMBER_INTERVALS];

double contactsenergyperinterval_pdf[NUMBER_INTERVALS];

void READCONTACTSENERGY(void);


int main(void)
{
	
	FILE *output_pdf;
	FILE *output_pdf_for_plot;
	FILE *output_n_pdf_for_plot;
	FILE *output_energy_contacts;
	
	int i, m, n;
	double j;
	double k;
	double max, min, l;
	
	
	for (i = 0; i < NUMBER_INTERVALS; i++)
	{
		contactsenergyinterval[i] = 0.0;
	}
	
	for (i = 0; i < NUMBER_INTERVALS; i++)
	{
		n_contactsenergyperinterval[i] = 0.0;
	}
	
	for (i = 0; i < NUMBER_INTERVALS; i++)
	{
		contactsenergyperinterval_pdf[i] = 0.0;
	}
	
	chdir(workingdir);
	
	READCONTACTSENERGY();
	
	for (i = 0; i<ALLPAIRS; i++)
	{
		printf("%d	%f\n", i+1, ContactsEnergy[i]);
	}
	
	j = -200.0;
	for (i = 0; i < ALLPAIRS; i++)
	{
		if (ContactsEnergy[i] > j)
		{
			j = ContactsEnergy[i];
		}
	}				
	
	max = j;
	printf("%f\n", max);
	
	k = 200.0;
	for (i = 0; i < ALLPAIRS; i++)
	{
		if (ContactsEnergy[i] < k)
		{
			k = ContactsEnergy[i];
		}
	}
	
	min = k; 
	
	l = min;
	printf("%f\n", l);
	
	for (i = 0; i < NUMBER_INTERVALS;)
	{
		if ((l >= min) && (l < max))
		{
			if (i == 0)
			{
				contactsenergyinterval[i] = l;
				i++;
			}
			else
			{
				if ((l + INTERVAL_WIDTH) <  max)
				{	
					contactsenergyinterval[i] = l + INTERVAL_WIDTH;
					l = l + INTERVAL_WIDTH;
					i++;
				}
				else
				{
					break;
				}
			}
		}
		
		else
		{
			break;
		}
		
	}
	
	n = i;
	
	
	for (i = 0; i < n; i++)
	{
		for (m = 0; m < ALLPAIRS; m++)
		{
			if ((ContactsEnergy[m] >= contactsenergyinterval[i]) && (ContactsEnergy[m] < contactsenergyinterval[i+1]))
			{
				n_contactsenergyperinterval[i]++;
			}
			else 
			{
				if ((ContactsEnergy[m] >= contactsenergyinterval[i]) && ((i+1) == n))
				{
					if ((ContactsEnergy[m] >= contactsenergyinterval[i]) && (ContactsEnergy[m] <= max))
					{
						n_contactsenergyperinterval[i]++;
					}
					
					else
					{
						continue;
					}
				}
			}
				
		}
	}
	
	for (i = 0; i < n; i++)
	{
		contactsenergyperinterval_pdf[i] = n_contactsenergyperinterval[i]/((double)ALLPAIRS*INTERVAL_WIDTH);
	}	
	
	chdir(writingdir);
	
	output_pdf = fopen("pdf_1_energy_B2M-D76N-I2-PH7p2_2000mc_T-0.21_newcf_n-clashes_0.1_energy_electrostatics.dat", "w");
	
	output_pdf_for_plot = fopen("pdf_1_energy_B2M-D76N-I2-PH7p2_2000mc_T-0.21_newcf_n-clashes_0.1_energy_electrostatics_plot.dat", "w");
	
	output_n_pdf_for_plot = fopen("pdf_1_energy_B2M-D76N-I2-PH7p2_2000mc_T-0.21_newcf_n-clashes_0.1_energy_electrostatics_n_plot.dat", "w");
	
	output_energy_contacts = fopen("contacts_energy_B2M-D76N-I2-PH7p2_2000mc_T-0.21_newcf_n-clashes_0.1_energy_electrostatics.dat", "w");	
	
	for (i = 0; i < n; i++)
	{
		fprintf(output_pdf, "%f          %d          %f\n", contactsenergyinterval[i], n_contactsenergyperinterval[i], contactsenergyperinterval_pdf[i]);
		fprintf(output_pdf_for_plot, "%f         %f\n", contactsenergyinterval[i], contactsenergyperinterval_pdf[i]);
		fprintf(output_n_pdf_for_plot, "%f       %d\n", contactsenergyinterval[i], n_contactsenergyperinterval[i]);
	}
	
	for (m = 0; m < ALLPAIRS; m++)
	{
		fprintf(output_energy_contacts, "%f\n", ContactsEnergy[m]);
	}
	
	fclose(output_pdf);
	
	fclose(output_pdf_for_plot);
	
	fclose(output_n_pdf_for_plot);
	
	fclose(output_energy_contacts);
	
	return 0;		
}


void READCONTACTSENERGY(void)
{
	FILE *handler;

	int i,dummie1,ID1,ID2,Contacts,TypeContacts,Clashes;
	float dummie2,dummie3,dummie5, EnergyContacts, EnergyHidro, EnergyElectro, EnergyHB;

	i=0;

	handler=fopen("B2M_D76N_I2_pH7p2_2000_T-0.21_n-clashes_0.1_energy_electrostatics_atticus.dat","r");

	while (!feof(handler))
	{
			if (fscanf(handler, "%d	%d	%d	%f	%f	%f 	%f	%f	%f	%d	%d	%d	%f", &ID1, &ID2, &dummie1, &dummie2, &dummie3, &EnergyContacts, &EnergyHidro, &EnergyElectro, &EnergyHB, &Contacts, &TypeContacts, &Clashes, &dummie5) > 0)
			{
				if (1)
				{
					ContactsEnergy[i]=EnergyContacts;
					//printf("%d	%d	%d	%f	%f	%f\n",ID1,ID2,dummie1,dummie2,dummie3,EnergyContacts);
					i++;
				}
			}
	};
	fclose(handler);

}	
	
