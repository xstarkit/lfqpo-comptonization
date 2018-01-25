// Include librarries
#include<stdio.h> 
#include<stdlib.h>
#include<string.h>

void main ()
{
  int       runnum;
  FILE     *fp;
  char     out_file[20];
  char     num[10];
  
  int i;

for(i=0; i<=63; i++)
{

  runnum = i;
  strcpy(out_file,"mtc__inp.000000");
  sprintf(num,"%06i",runnum);
  memcpy(out_file+9,num,6);


  fp = fopen(out_file,"w");

    fprintf(fp,"%s","------temp-------energy?---\n");
/*  'energy' for which the cumul. distr. of ro1' will be calculated  */
   fprintf(fp,"%i %i ", 100, 1);
   fprintf(fp, "\n");
   fprintf(fp,"%s", "--spectrum---rad_temp--alpha---emin------emax----nop_spec-spec_file_name--\n");
   fprintf(fp,"%i %lf %lf %lf %i %i %s", 0, 0.1, 0.6, 0.001, 10, 200, "photar.dat");
   fprintf(fp, "\n");   
   fprintf(fp,"%s", "---e_min_out--e_max_out--num_bin_inp--num_bin_out--nop_grid--grid_fname-\n");
   fprintf(fp,"%lf %i %i %i %i %s", 0.001, 1000, 10, 10, 128, "ginga_ear.dat");
   /*   printf("%.2e %.2e %ld %ld \n",emin_o,emax_o,num_bin_i,num_bin_o); */
   fprintf(fp, "\n");
   fprintf(fp,"%s", "--number_of_MC_steps-\n");
   fprintf(fp,"%i", 100000000);
   fprintf(fp, "\n");
   fprintf(fp,"%s","-spatial_distr--r_emission----irradiation--theta_irrad--nangle---\n");
   fprintf(fp,"%i %i %i %i %i", 6, -1, 0, 60, 10);
   fprintf(fp,"\n");
   fprintf(fp,"%s", "--torus:-dh/dr-theta0---prectheta----precphi---prec_option-(also:Rin)\n");
   fprintf(fp,"%i %i %i %lf %i", 1, 15, 15, i*360./64., 2);
   fprintf(fp,"\n");
   fprintf(fp,"%s", "--taumax----Rin[R_g]------dens_alpha---rho0/Rout----\n");
   fprintf(fp,"%lf %lf %lf %lf", 1.76, 5.0, 0.01, 30.0);
   fprintf(fp,"\n");
   fprintf(fp,"%s", "---cold_refl----b_bulk--\n");
   fprintf(fp,"%i %lf", 1, 0.); 
   fprintf(fp,"\n");
   fprintf(fp,"%s", "--dist_max-----dist_points----\n");
   fprintf(fp,"%i %i", 30, 256);
   fprintf(fp,"\n");
   fprintf(fp,"%s", "---rand_seed---xo--b--\n");
   fprintf(fp,"%i %i %i", -47157, 1, 1);
   fprintf(fp,"\n");
   fprintf(fp,"%s", "---num_en_dist--\n");
   fprintf(fp,"%i", 0);
   fprintf(fp,"\n");
   fprintf(fp,"%s", "--max_scsp----\n");
   fprintf(fp,"%i", 5);
   fprintf(fp,"\n");
   fprintf(fp,"%i %i", 0, 0);
   fprintf(fp,"\n");
   fprintf(fp,"%i %i", 1, 1);
   fprintf(fp,"\n");
   fprintf(fp,"%i %i", 2, 2);
   fprintf(fp,"\n");
   fprintf(fp,"%i %i", 3, 3);
   fprintf(fp,"\n");
   fprintf(fp,"%i %i", 4, 100);
   fprintf(fp,"\n");
   fprintf(fp,"%s", "--bh_mass--bh_spin--bh_mdot--rotation--\n");
   fprintf(fp,"%lf %lf %lf %i", 10.0, 0.3, 0.01, 0);
   fprintf(fp,"\n");   
   fclose(fp);
}
}

