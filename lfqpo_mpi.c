/*
Author: Bei You (youbeiyb@gmail.com) 

  Compton scattering in a uniform, hot  plasma cloud.
  Based on Gorecki & Wilczewski and Pozdnyakov, Sobol & Sunyaev

  This code is the combination of Piotr's (Comptonization) and Michal's (GR effect) codes.
  
  Currently, it only concentrates on the wedge geometry for Lense-Thirring precession project (LFQPO). 
  It will be extented to other geometries after this project.

-----------------------------------------------------------------------------
*/

#include "mtc_incl_def.c"
#include "mtc_incl_code.c"
#include "quadrat.c"
#include "sim5lib.c"
#include "mpi.h"

/* 链表节点定义 */
typedef struct tag_LIST_NODE
{
    struct tag_LIST_NODE *pNext;
    double pos[4], mom[4], pp, ee;
}LIST_NODE;

/* 链表定义 */
LIST_NODE *pHead = NULL;
LIST_NODE *pTail = NULL;

#define vect_copy(v1,v2) {v2[0]=v1[0]; v2[1]=v1[1]; v2[2]=v1[2]; v2[3]=v1[3];}

double disk_photons_pdf(double x);
double disk_gfactor_K(double r, double a, double l, double q, double rms);
sim5distrib dist;

long int Ntot=0, Ntot_all, 
		 Ntot_tmp=0, Ntot_tmp_all, Ntot_disk=0, Ntot_torus=0, Ntot_refl=0, 
		 Nlimit0=1e8, Nlimit;
int refl_rept = 200, line_rept=10, comp_rept=5;
int num_bin_line = 200; 
double elmin_o = 1.0/511.0, elmax_o = 10.0/511.0;
//--------------------------------------------------------------- 
// notice: the grid here has to be consistent with the grid defined in another code 'photon_table.c' 
float ra[nrd+1], fa[nfd+1], tha[nth+1], pha[nph+1];
float inf_ph0[nrd+1][nfd+1][nth+1][nph+1], inf_ph1[nrd+1][nfd+1][nth+1][nph+1],
inf_ph2[nrd+1][nfd+1][nth+1][nph+1], inf_ph3[nrd+1][nfd+1][nth+1][nph+1],
inf_ph4[nrd+1][nfd+1][nth+1][nph+1], inf_ph5[nrd+1][nfd+1][nth+1][nph+1],
inf_ph6[nrd+1][nfd+1][nth+1][nph+1], inf_ph7[nrd+1][nfd+1][nth+1][nph+1];
//--------------------------------------------------------------- 

int main(int argc, char *argv[])
{
    long int il, istep;

    double    position[4], momentum[4], momentum_torus[4], weight;
    double    position_esc[4], momentum_esc[4], weight_esc;
    double 	  gi, gf; // dotprod(p,U)
    int       runnum;
    int 	  i, indx, iindx, findx;
	double 	  E_loc, E_inf;

  runnum = 0;
  if (argc == 2)   {
    runnum = atoi(argv[1]);
  }

int j,k;
  for(i=0; i<1000; i++) {
    for(j=0; j<10; j++) {
      for(k=0; k<4; k++){
		Ndisk[i][j][k] = 0;
        Ndisk_all[i][j][k] = 0;
	    Ncomp[i][j][k] = 0;
	    Ncomp_all[i][j][k] = 0;
        Nrefl[i][j][k] = 0;
        Nrefl_all[i][j][k] = 0;
        Nline[i][j][k] = 0;
        Nline_all[i][j][k] = 0;        
      }
    }
  }

/* initialize the geometry of:  1. the torus
                                2. the disk
                                3. the axis of the torus */
    input_data(runnum); // read input parameters and set up the geometry of the torus
  if (do_reflection) {
    read_opacity();
  }    
       
	double r_min = Rdisk;          // inner radius of the disk
    double r_max = Rmax;         // outer radius of the disk 
/* Note: if r_min and r_rmax are changed, check the following rmin and rmax as well */ 

    // setup NT disk
    disk_nt_setup(bh_mass, bh_spin, bh_mdot, 0.1, 0);

    // setup photon distribution
    distrib_init(&dist, &disk_photons_pdf, log(r_min), log(r_max), 100);    

/* convert the axis of the torus (0, 0, 1) in torus frame to Cartesian coordiante of BH, according to (prectheta, precphi), 
in order to check if photon is inside of the torus by means of Dot Product */
    if (precession_option<2) {
      torus_axis1(); // precession model (1)
    } else {
      torus_axis2(); // precession model (2)
    } 
//	printf("axis_vec[0]=%e, axis_vec[1]=%e, axis_vec[2]=%e\n", axis_vec[0], axis_vec[1], axis_vec[2]);
//	getchar();

//--------------------------------------------------------------- 
// prepare the grid and read the data file "photon_table.dat" 
// notice: the grid here has to be consistent with the grid defined in another code 'photon_table.c'      
double dr = (log10(r_max)-log10(r_min))/nrd,
       df = (f_max-f_min)/nfd, 
       dt = (t_max-t_min)/nth, // cos(th0)
       dp = (p_max-p_min)/nph; // cos(ph0) 
int ix, iy, iz, it, ii;  
double r, phi, rnd_th, rnd_ph;

//FILE *data = fopen("photon_table.dat", "r"); // opens file for reading 
FILE *data; // opens file for reading 

char table_file[20];
char num[10];
strcpy(table_file,"table000000.dat");
sprintf(num,"%06i",runnum);
memcpy(table_file+5,num,6);		
data = fopen(table_file,"r");         /* output observed torus Comptonization spectrum */
  
if(data)
{
printf("Read photon_table from file\n");
for (ix=0; ix<=nrd; ix++) {	
	for (iy=0; iy<=nfd; iy++) {
		for(iz=0; iz<=nth; iz++){
			for(it=0; it<=nph; it++){
        r = pow(10.0,(log10(r_min)+ix*dr));
        ra[ix] = r;
		phi = f_min + iy*df; 
		fa[iy] = phi;
		rnd_th = t_min + iz*dt;
		tha[iz] = rnd_th;
		rnd_ph = p_min + it*dp;
		pha[it] = rnd_ph;
				fscanf(data, "%f %f %f %f %f %f %f %f", &inf_ph0[ix][iy][iz][it], &inf_ph1[ix][iy][iz][it],	
														&inf_ph2[ix][iy][iz][it], &inf_ph3[ix][iy][iz][it],	
														&inf_ph4[ix][iy][iz][it], &inf_ph5[ix][iy][iz][it],	
														&inf_ph6[ix][iy][iz][it], &inf_ph7[ix][iy][iz][it]);																																						
			}
		}	
	}
} 
}
fclose(data);
printf("reading disk_table is done\n");
//getchar();

read_sphere();
printf("reading sphere_table is done\n");

/*int ij, ik;
for (ii=0; ii<=npsp; ii++) printf("cosa = %e\n", cosa[ii]);	
	for (ij=0; ij<=ntsp; ij++) printf("kta = %e\n", kta[ij]);	
		for(ik=0; ik<=nfsp; ik++) printf("kfa = %e\n", kfa[ik]);	
getchar();
*/ 
//---------------------------------------------------------------           
        
int not_scattering, i_path;

//---------------------------------------------------------------------
int myid, numprocs;

    MPI_Init(&argc,&argv);

    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    Nlimit=Nlimit0/numprocs;
//---------------------------------------------------------------------
int ij=0;
double position0[4], momentum0[4], E_inf0;
idum += (100*myid);  
for (il=myid; il<=mc_steps; il+=numprocs) {// generate photons	
	if(il%1000000==0){ 
      /*        output_results(runnum,il); */
      writelog(runnum,il); 
    }
	
	generate_photon_gr(position, momentum, &E_inf);  	
//	continue;  
    vect_copy(position, position0);
    vect_copy(momentum, momentum0);
    E_inf0 = E_inf;
    for(ij=0;ij<comp_rept;ij++){
			int icomp=0;
			do{ // loop for multi-scattering in the torus
				icomp++;
				if(icomp>100){printf("too many Comp scattering in torus\n"); getchar();}
				
				/* determine the scattering position (position, momentum) and escaping position (position_esc, momentum_esc) */
				not_scattering = raytrace_in_torus(position, momentum, E_inf); 

				if(not_scattering==-1) break; /* failed raytrace, skipping this photon, Normally raytrace inside torus should always succeed */
			
				// this is the not-scattering part; photon will go to infinity, or disk (then reflection)
				i_path = 0;
				if(not_scattering==1){		
					i_path = raytrace_esc_torus(position, momentum, E_inf);	
//					break;
//					printf("here2\n"); 				
					if(i_path==1){ 
//						printf("here to torus\n");
						continue; // photon enter into torus again, so directly move to "raytrace_in_torus"						
					}else{
						break;
					}
				}	
				
				if(not_scattering==0){ 																		     																	     
					prod_pU(position, momentum, &gi); 
					E_loc = E_inf/gi; /* calculate the photon energy in the position where it will be scattered */	
					// this is the scattering part;
					scattering_gr(position, momentum, &E_loc); // this is the scattering; global variable "E_loc" is already changed. 
					prod_pU(position, momentum, &gf); // Since photon momentum is changed after scattering, new "gf" is required for calculating photon energy at infinity					
					E_inf = E_loc*gf;
				}	
					
			}while((E_inf > emin_o)&&(E_inf < emax_o));
    vect_copy(position0, position);
    vect_copy(momentum0, momentum);
    E_inf = E_inf0;    
    }                                    
} // end of all photons 
  printf("thread %d done\n", myid); 
                  
//---------------------------------------------------------------------
  int nos = num_bin_o*log10(emax_o/emin_o);
  for (i=0; i<nos; i++){
      for(j=0; j<nangle; j++) {
		for(k=0; k<MAXphi; k++) {
	MPI_Reduce(&Ndisk[i][j][k],&Ndisk_all[i][j][k], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Ncomp[i][j][k],&Ncomp_all[i][j][k], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Nrefl[i][j][k],&Nrefl_all[i][j][k], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);    
		}
      }
   }
  int nos1 = num_bin_line; 
  for (i=0; i<nos1; i++){
      for(j=0; j<nangle; j++) {
		for(k=0; k<MAXphi; k++) {
    MPI_Reduce(&Nline[i][j][k],&Nline_all[i][j][k], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		}
      }
   }   	
   MPI_Reduce(&Ntot,&Ntot_all, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&Ntot_tmp,&Ntot_tmp_all, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);   
   
if (myid == 0){
  for (i=0; i<nos; i++){
      for(j=0; j<nangle; j++) {
        for(k=0; k<MAXphi; k++) { 
				out_disk[i][j][k] = Ndisk_all[i][j][k];
				out_comp[i][j][k] = Ncomp_all[i][j][k];
				out_refl[i][j][k] = Nrefl_all[i][j][k];				
		}
      }	
   }
  for (i=0; i<nos1; i++){
      for(j=0; j<nangle; j++) {
        for(k=0; k<MAXphi; k++) { 
				out_line[i][j][k] = Nline_all[i][j][k];
		}
      }	
   }   
  distrib_done(&dist);
  output_spec(runnum,Ntot_all); /* Output spectrum: 
					                 Disk reflection spectrum;
					                 Torus Comptonization spectrum;
								*/
}   
MPI_Finalize();

	double speed;						
	speed = (1.0*mc_steps)/(time(NULL)-start_time);
    printf(" Speed of computations: %.1f photons/sec \n",speed);		
    					
    return 0;
} // end of main program 

/*  ***********************************************************************************************************************  */

writelog(runnum,noph)
   long int noph;
        int runnum;
{
    FILE   *fl;
    double speed;
    char   log_file[20],num[10];

    strcpy(log_file,"mtc__log.000000");
    sprintf(num,"%06i",runnum);
    memcpy(log_file+9,num,6);

    if (noph == 0) {
      fl = fopen(log_file,"w");
      fprintf(fl,"Uniform density run.  %s \n",
	      asctime(localtime(&start_time)));
      fprintf(fl,
	      "Spatial distribution: %d (%.1f),  Rin: %.2f,  rndseed: %ld \n",
	      spat_distr,r_emission,Rin,idum0);
      fprintf(fl,"Tau: %.2f,  plasma temp (keV): %.1f,  T0 (keV): %.3f \n",
	      tau_max,511./med_temp,rad_temp*511.);
    } else {
      fl = fopen(log_file,"a");
      speed = (1.0*noph)/(time(NULL)-start_time);
      printf("%6li  %.1f \n",noph,speed);
      fprintf(fl,"%6li  %.1f \n",noph,speed);
    }
    fclose(fl);
    
    return 0;
}


int write_spec(position, momentum, energy)
double position[], momentum[], energy;
{
/* -------------------- The block for saving photon energy and position ---------------------------------- */

	    // get final theta-position at infinity
        double theta_inf = position[2];	// cos(theta)
        
        // get final phi-position at infinity
        double phi_inf = position[3];
        phi_inf = reduce_angle_2pi(phi_inf)*180.0/M_PI;
        
		/* bins in \phi: 0-20, 90-110, 180-200, 270-290 degs */
	    int findx = -1; 
	    int i, indx, iindx;
	    for(i=0; i<MAXphi; i++) {
	      if (fabs(phi_inf-(i*90+10)) <= 10) {  
	       findx = phi_inf/90.;
	      }
	    }
	    
	if(theta_inf>=0.0){ // only deal with photon energy in the interesting range 	
		indx = num_bin_o*log10(energy/emin_o); // indx corresponding to photon energy at infinity 
			    
	  if(findx>=0){
	    iindx = theta_inf*nangle;
		Ncomp[indx][iindx][findx] += 1.0;		 
	  }	
	}                
/* -------------------- The block for saving photon energy and position ---------------------------------- */
return 0;	
}


output_spec(runnum,noph)
  long int   noph;
  int        runnum;
{
  FILE   *fp1, *fp2, *fp3, *fp4;
  int     i,nos,j,k; 
  double  e,st,si,sr,er,d,de;

  char   out_file[20],
         inp_spec[20],
         reflspec[20],
         linespec[20];
  char    num[10];

  if (noph == 0) return 0;

  strcpy(out_file,"mcomp000000.dat");
  strcpy(inp_spec,"inpsp000000.dat");
  strcpy(reflspec,"refls000000.dat");
  strcpy(linespec,"lines000000.dat");  


  sprintf(num,"%06i",runnum);
  memcpy(out_file+5,num,6);
  memcpy(inp_spec+5,num,6);
  memcpy(reflspec+5,num,6);
  memcpy(linespec+5,num,6);  
		
fp1 = fopen(inp_spec,"w");         /* output observed Novikov-Thorne disk spectrum */
fp2 = fopen(reflspec,"w");         /* output observed reflection spectrum from the disk */
fp3 = fopen(out_file,"w");         /* output observed torus Comptonization spectrum */
fp4 = fopen(linespec,"w");         /* output observed Fe Ka line */

double er_disk, er_torus; 
  nos = num_bin_o*log10(emax_o/emin_o);
  for (i=0; i<nos; i++){
    e = emin_o*pow(10.,(i+0.5)/num_bin_o );
    de = emin_o*(pow(10.,(i+1.0)/num_bin_o)-pow(10.,(1.0*i)/num_bin_o));
    er = e/de/noph;
    for(j=0; j<nangle; j++) {
		for(k=0; k<MAXphi; k++) {// 基本思想是, 当实验次数足够多时, 各个成分之间的比例一定; 比如对于disk成分, 不依赖noph
			fprintf(fp1,"%e %.3e %.3e ",e*511., out_disk[i][j][k], out_disk[i][j][k]*(e/de/Ntot_tmp_all)*nangle);
			fprintf(fp2,"%e %.3e %.3e ",e*511., out_refl[i][j][k], out_refl[i][j][k]*er*nangle/refl_rept/comp_rept);			
			fprintf(fp3,"%e %.3e %.3e ",e*511., out_comp[i][j][k], out_comp[i][j][k]*er*nangle/comp_rept); // 在 noph*comp_rept 的基础上
		}
		fprintf(fp1,"  ");
		fprintf(fp2,"  ");		
		fprintf(fp3,"  ");
    }
    fprintf(fp1,"\n");
    fprintf(fp2,"\n");    
    fprintf(fp3,"\n");
  }
  fclose(fp1);
  fclose(fp2);  
  fclose(fp3);
//--------------------------------------------------------------------
  int nos1 = num_bin_line;
  de = (elmax_o-elmin_o)/num_bin_line;
  for (i=0; i<nos1; i++){
    e = elmin_o + (i+0.5)*de;
    er = e/de/noph;
    for(j=0; j<nangle; j++) {
		for(k=0; k<MAXphi; k++) {
			fprintf(fp4,"%e %.3e %.3e ",e*511., out_line[i][j][k], out_line[i][j][k]*er*nangle/refl_rept/line_rept/comp_rept);
		}
		fprintf(fp4,"  ");
    }
    fprintf(fp4,"\n");
  }
  fclose(fp4);

return 0;
}

int generate_photon_gr(position, momentum, E_inf)
double position[], momentum[], *E_inf; 
{

double r, T, energy;
double th0, ph0, momentum_on[4];
double dl;
double Omega;
//raytrace_data rtd;
sim5tetrad t;
sim5metric m;	
geodesic gd;
int status;
double U[4], gf;
double cos_inc, phi_inf, phi, dphi, gd_P;
int findx; 
int i, indx, iindx;
int iter = 0;
int	inside_torus = 0;
do{		iter++;
        // draw an emission radius from distribution        
        do {
        r = exp(distrib_hit(&dist));
        T = disk_nt_temp(r); 
//		printf("r = %.2e, T = %.2e\n", r, 2.7*T*8.617328149741e-8);
		} while ((2.7*T*8.617328149741e-8/511.) < emin); // convert temperature to keV (normalized to 511 keV), to set up minimum temperature (maximum disk outer radius)

/* Generate energy from a planckian distribution of a corresponding temperature */
		do{
		  energy = blackbody_photon_energy_random(T)/511.; /* photon energy in disk-rest frame */
		}while ( (energy < emin) || (energy > emax));	

     	// generate pos[4] in B-L spherical coordinate
	    position[0] = 0.0;
	    position[1] = r;
	    position[2] = 0.0;
	    position[3] = urand*M_PI*2.0;
	    phi = position[3];
    
        // get metric and tetrad
        kerr_metric(bh_spin, r, 0.0, &m);   
		Omega = OmegaK(r,bh_spin);
        tetrad_azimuthal(&m, Omega, &t); 
        fourvelocity_azimuthal(Omega, &m, U);    
     
        // pick direction in disk rotating frame
        //double th0 = urand*M_PI/2.;
        th0 = acos(sqrt(1.-sim5urand())); // iCDF that corresponds to PDF(theta)=sin(theta)*cos(theta)	
        if(rnd()<0.5) th0 = M_PI-th0;	
        ph0 = urand*M_PI*2.;   	            
        momentum_on[0] = 1.0;
        momentum_on[1] = sin(th0)*cos(ph0);
        momentum_on[2] = cos(th0);
        momentum_on[3] = sin(th0)*sin(ph0);
        on2bl(momentum_on, momentum, &t); // transfer the direction from disk frame to B-L 

//	if(fabs(momentum[1])<1e-7) continue; 
        
        gf = (momentum[0]*m.g00 + momentum[3]*m.g03) / dotprod(momentum, U, &m);
		*E_inf = energy*gf;   // photon energy measured by observer at infinity, namely E_inf  

// case 1:		
//--------------------------------------------------------------------------------------------------
// do interpolation:
//-----------------------------------------------------
double rnd_th = cos(th0);
double rnd_ph = ph0;
double tmp_ph0, tmp_ph1,tmp_ph2,tmp_ph3,tmp_ph4,tmp_ph5,tmp_ph6,tmp_ph7;
int interp = 0;
interp = interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph1,&tmp_ph1);      
if(interp==0){
	if(tmp_ph1<0.0) continue; // photon travels to infinity but with cosi < 0
    if(tmp_ph1>=1e4){// photon goes to infinity, assume 1e4 larger enough than torus maximum radius 
			Ntot++; // !!! first accumulate photon number
			if(Ntot_tmp<Nlimit){ // If the disk photon is enough for good statistic, don't have to save
	interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph0,&tmp_ph0);
	interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph2,&tmp_ph2);
	interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph3,&tmp_ph3);  				 
    interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph4,&tmp_ph4);
    interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph5,&tmp_ph5);
    interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph6,&tmp_ph6);
    interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph7,&tmp_ph7);
				// return theta_infinity angle and g-factor
				cos_inc = tmp_ph2;
				phi_inf = rad2deg(reduce_angle_2pi(tmp_ph3));
				//	printf("inc = %e, phi = %e, g = %e\n",cos_inc, phi_inf, g_inf);
				
						/* bins in \phi: 0-20, 90-110, 180-200, 270-290 degs */
						findx = -1;
						for(i=0; i<MAXphi; i++) {
							if (fabs(phi_inf-(i*90+10)) <= 10) {  
								findx = phi_inf/90.;
							}
						}							
						indx = num_bin_o*log10((energy*gf)/emin_o); // indx corresponding to photon energy at infinity 	
						iindx = cos_inc*nangle;	 								
					if((cos_inc>0.0) && (findx>=0)){ // only deal with photon energy in the interesting range								
								Ndisk[indx][iindx][findx] += 1;
					} 
			Ntot_tmp=Ntot; // before reaching Nlimit, counting and saving disk photon	
			}
			continue; // continue to next photon						
	}else{	
	// then photon goes to torus; return position[], momentum[] and E_inf as initial parameters for Comptonization		
		interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph0,&tmp_ph0);
		interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph2,&tmp_ph2);
		interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph3,&tmp_ph3);
		double interp_pos[4];
		interp_pos[0]=tmp_ph0; interp_pos[1]=tmp_ph1; interp_pos[2]=tmp_ph2; interp_pos[3]=tmp_ph3;		

	  if(inside_dotprod(interp_pos)){// nin some cases, near the boundaries of the torus, the interpolated position is just out of the torus.
									 // if so, apply raytrace()
		position[0]=tmp_ph0; position[1]=tmp_ph1; position[2]=tmp_ph2; position[3]=tmp_ph3;		  		
		interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph4,&tmp_ph4);
		interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph5,&tmp_ph5);
		interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph6,&tmp_ph6);
		interp4d(r,phi,rnd_th,rnd_ph,ra,fa,tha,pha,inf_ph7,&tmp_ph7); 				
		momentum[0]=tmp_ph4; momentum[1]=tmp_ph5; momentum[2]=tmp_ph6; momentum[3]=tmp_ph7;
		norm_vect(position,momentum); // normalize the k-vector
		return 0; // end of this subroutine 'generate_photon_gr'	
      }
	}
	
}
// the interpolation is not applicable for this case
// do direct raytracing for this photon                                                                				
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 	

		double pos[4], k[4];
		vect_copy(position, pos);
		vect_copy(momentum, k);

        // get geodesic
        // this geodesic can also be used in case 3
		geodesic_init_src(bh_spin, r, 0.0, momentum, momentum[1]<0?1:0, &gd, &status);
		
				if ((!status) || (r < gd.rp)) {
				fprintf(stderr, "WRN: geodesic_init_src failed (status=%d, r=%.2e, rp=%.2e)\n", status, r, gd.rp);
				continue;
				}

				if(gd.m2p==1.0)	continue;

/*				// sometime a trajectory is found, which is correct, but unphysical
				// (photon would have pass radial turning point, which lies bellow horizon)
				if ((creal(gd.r1) < r_bh(bh_spin)) && (momentum[1] < 0.0)) {
				// unphysical solution - geodesic goes to observer through BH
				continue;
				}
*/	
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<	
	
// case 2:
			if((gd.rp>Rout) && (gd.cos_i>0.0)){	
				gd_P = geodesic_P_int(&gd, r, momentum[1]<0?1:0);
                if(((gd_P > gd.Rpa) && (gd.nrr == 2)) || (gd.m2p==1.0)) continue;                
				dphi = geodesic_position_azm(&gd, r, 0.0, gd_P);
				if(isnan(dphi)) continue;
									
				// return theta_infinity angle and g-factor
				cos_inc = gd.cos_i;
				phi_inf = rad2deg(reduce_angle_2pi(phi + dphi));
				//	printf("inc = %e, phi = %e, g = %e\n",cos_inc, phi_inf, g_inf);
				
						/* bins in \phi: 0-20, 90-110, 180-200, 270-290 degs */
						findx = -1;
						for(i=0; i<MAXphi; i++) {
							if (fabs(phi_inf-(i*90+10)) <= 10) {  
								findx = phi_inf/90.;
							}
						}							
						indx = num_bin_o*log10((energy*gf)/emin_o); // indx corresponding to photon energy at infinity 	
						iindx = cos_inc*nangle;	 								
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 					
				// if the radial turning point is outside torus extent (meanwhile photon moves inwards BH), 
				// photon will never enter into torus; Instead, going to infinity
				Ntot++; // !!! first accumulate photon number	
					if(Ntot_tmp<Nlimit){ // If the disk photon is enough for good statistic, don't have to save 					
						if(findx>=0){ // only deal with photon energy in the interesting range	  
								Ndisk[indx][iindx][findx] += 1;
						}
					Ntot_tmp=Ntot; // before reaching Nlimit, counting and saving disk photon	
					}			
				continue; // !!! must go to the new loop
			} 										
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<	
// case 3:                                                                     
        /* ==================== The block for ray-tracing ================================= */
        raytrace_data rtd;
	    raytrace_prepare(bh_spin, position, momentum, NULL, 1.0, 0, &rtd);
     		
    int iloop = 0;   		
 	do{ iloop++;
	    dl = 1e6; // use maximal step     
	    raytrace(position , momentum, NULL, &dl, &rtd);
	    // after raytrace(), posi and momi already changed
	    
	    		// stop condition: if applied, stop raytrace
				if (position[1] < r_bh(bh_spin)) {
//				printf("raytrace(): photon goes to bh event horizon\n");
				break;  
				}
				// also stop if relative error this step is too large
				if (rtd.error>1e-2) {
//				printf("raytrace(): aborted due to large error\n");
				break;
				}

		if(position[1]<=(1.1*Rout)){		
        inside_torus = inside_dotprod(position);
		} 
				
		if((position[1]>(1.1*Rout))&&(momentum[1]>0.0)){ // if this applies, use geodesic function instead of raytrace()
			if(gd.cos_i>0.0){
				gd_P = geodesic_P_int(&gd, r, k[1]<0?1:0);
                if(((gd_P > gd.Rpa) && (gd.nrr == 2)) || (gd.m2p==1.0)) break;
				dphi = geodesic_position_azm(&gd, r, 0.0, gd_P);
				if(isnan(dphi)) break;				
									
				// return theta_infinity angle and g-factor
				cos_inc = gd.cos_i;
				phi_inf = rad2deg(reduce_angle_2pi(phi + dphi));
					
						/* bins in \phi: 0-20, 90-110, 180-200, 270-290 degs */
						findx = -1;
						for(i=0; i<MAXphi; i++) {
							if (fabs(phi_inf-(i*90+10)) <= 10) {  
								findx = phi_inf/90.;
							}
						}							
						indx = num_bin_o*log10((energy*gf)/emin_o); // indx corresponding to photon energy at infinity 	
						iindx = cos_inc*nangle;										
			
		// !!! very important condition, so that to reduce the spending time to find the entering-torus photon 	
			Ntot++; // !!! first accumulate photon number	
					if(Ntot_tmp<Nlimit){ // If the disk photon is enough for good statistic, don't have to save 					
						if(findx>=0){ // only deal with photon energy in the interesting range	  
								Ndisk[indx][iindx][findx] += 1;
						}
					Ntot_tmp=Ntot; // before reaching Nlimit, counting and saving disk photon	
					}
			}							
			break;
		}
				 
		if(iloop>10000) {printf("too much step \n");break;}   
			 
	}while(inside_torus!=1);
		
}while((inside_torus!=1) && (iter < 10000));     	

if(iter>=10000){
	printf("iter: %d\n", iter);
//	getchar();
}    	
                 				                                                        
return 0;
}


int raytrace_out_torus(position, momentum, i_path) 
/* a function of raytracing photon which is outside torus */
double position[], momentum[];
int *i_path; /* = 3, failed trajactory;
				= 0, to infinity
				= 1, to disk
				= 2, to torus
			 */	  	
{
double dl, costi, costf;
double posi[4], posf[4], ki[4], kf[4];
int i, inside_torus = 0;;
*i_path = 3; // default value 

//return 0;	
if(position[1]<Rout){
// case 1#: first check if photon enter the torus again within "R = Rout" ------------------------------------------------------
        /* ==================== The block for ray-tracing ================================= */
	    raytrace_data rtd;
	    raytrace_prepare(bh_spin, position, momentum, NULL, 10.0, 0, &rtd);
int icase1=0;		
	do{
icase1++;
if(icase1>100000) {printf("too many raytrace in case1 which is for photon escaping torus\n"); getchar();}		
	    dl = 1e6; // use maximal step   	  
//	    if(position[1]>0.9*Rout){
			vect_copy(position, posi);
			vect_copy(momentum, ki);
//		}
	    raytrace(position , momentum, NULL, &dl, &rtd);
//	    if(position[1]>0.9*Rout){
			vect_copy(position, posf);
			vect_copy(momentum, kf);
//		}
	    // after raytrace(), posi and momi already changed
	    
	    		// stop condition: if applied, stop raytrace
				if ((position[1] < r_bh(bh_spin))) {
//				printf("out torus1, raytrace(): photon goes to bh event horizon\n");
				*i_path = 3;  
				return 0;  
				}
				// also stop if relative error this step is too large
				if (rtd.error>1e-2) {
//				printf("out torus2, raytrace(): aborted due to large error (generally it is due to photon going to BH)\n");
				*i_path = 3; // failure step; as photon goes to BH event horizon, 
				return 0;
				}

				if(isnan(position[1]) || isnan(position[2]) || isnan(position[3])){ // failed geodesic	
				*i_path = 3;  
				return 0;  
				}								
		
        int inside_torus = inside_dotprod(position);
		if(inside_torus==1){
			 *i_path = 2;
			 return 0; 
		}
			
		if(position[1]>=Rout){ // escape out of sphere of radius = Rout		
				for(i=0;i<4;i++){
					position[i] = posi[i] + (Rout-posi[1])/(posf[1]-posi[1])*(posf[i]-posi[i]);	
					momentum[i] = ki[i] + (Rout-posi[1])/(posf[1]-posi[1])*(kf[i]-ki[i]);					
				}	
				position[1] = Rout; 				
				norm_vect(position, momentum); // normalize k-vector	
//				printf("just out of sphere\n");				
				break; // go to case 2 
		}  	
	}while(1);		
		/* ==================== End of the block for ray-tracing ================================= */ 
}		
//return 0;

// case 2#: if photon escape the south sphere, ignore it and then go to case 3# ----------------------------------------
if(position[2]<=0.0){
	*i_path = 3;
	return 0;
}	
if(momentum[1]<=0.0){
printf("stop here\n");
printf("k = %e\n", momentum[1]);
getchar();
}
/*printf("out of sphere:\n");	
printf("pos=%e, %e, %e, %e\n", position[0], position[1], position[2], position[3] );
printf("k=%e, %e, %e, %e\n", momentum[0], momentum[1], momentum[2], momentum[3]);
*/	
// case 3#: out of Rout, check if photon goes to infinity or disk ------------------------------------------------------	

/* independent compare the interpolation with the raytrace */
//double pos[4], k[4];	
//raytrace_check(position, momentum, pos, k);	
//printf("raytrace check:\n");	
//printf("pos=%e, %e, %e, %e\n", pos[0], pos[1], pos[2], pos[3]);
//printf("k=%e, %e, %e, %e\n", k[0], k[1], k[2], k[3]);	

//--------------------------------------------------------------------------------------------------
// do interpolation:
//-----------------------------------------------------
double r = position[1];
double m = position[2];
double phi_ori = position[3];
double momentum_on[4], k_loc[4], kk;
sim5tetrad t; 
sim5metric metric;
                
// get metric and tetrad
kerr_metric(bh_spin, position[1], position[2], &metric);
                   
tetrad_zamo(&metric, &t);
bl2on(momentum, momentum_on, &t);
vect_copy(momentum_on, k_loc);
k_loc[0] = 1.0; k_loc[1] /= k_loc[0]; k_loc[2] /= k_loc[0]; k_loc[3] /= k_loc[0]; 
kk = sqrt(sqr(k_loc[1]) + sqr(k_loc[2]) + sqr(k_loc[3]));
k_loc[1] /= kk; k_loc[2] /= kk; k_loc[3] /= kk;  
  
double th0=k_loc[2], ph0=acos(k_loc[1]/sin(acos(th0)));
if(sin(acos(th0))*sin(ph0)/k_loc[3]<0.0){
 ph0 = M_2PI-ph0;	
}
//printf("cost0=%e, ph0=%e\n", th0, ph0);
double tmp_ph0, tmp_ph1,tmp_ph2,tmp_ph3,tmp_ph4,tmp_ph5,tmp_ph6,tmp_ph7;
int interp = 0;
interp = interp3d(m,th0,ph0,cosa,kta,kfa,sp_ph1,&tmp_ph1);
//printf("interp = %d\n", interp);  
if(interp==0){icc1++;
//			  return 0; 
interp3d(m,th0,ph0,cosa,kta,kfa,sp_ph0,&tmp_ph0);  
interp3d(m,th0,ph0,cosa,kta,kfa,sp_ph2,&tmp_ph2);  
interp3d(m,th0,ph0,cosa,kta,kfa,sp_ph3,&tmp_ph3);  
interp3d(m,th0,ph0,cosa,kta,kfa,sp_ph4,&tmp_ph4);  
interp3d(m,th0,ph0,cosa,kta,kfa,sp_ph5,&tmp_ph5);  
interp3d(m,th0,ph0,cosa,kta,kfa,sp_ph6,&tmp_ph6);  
interp3d(m,th0,ph0,cosa,kta,kfa,sp_ph7,&tmp_ph7);  

position[0] = tmp_ph0, position[1] = tmp_ph1, position[2] = tmp_ph2, position[3] = tmp_ph3 + phi_ori; // remember to shift the azimuthal angle
momentum[0] = tmp_ph4, momentum[1] = tmp_ph5, momentum[2] = tmp_ph6, momentum[3] = tmp_ph7;
/*printf("interp:\n");	
printf("pos=%e, %e, %e, %e\n", position[0], position[1], position[2], position[3] );
printf("k=%e, %e, %e, %e\n", momentum[0], momentum[1], momentum[2], momentum[3]);
getchar();
*/ 
/*if(fabs((position[3]-pos[3])/position[3])>1e-1)
{
printf("big error1\n");
printf("position=%e, %e, %e, %e\n", position[0], position[1], position[2], position[3] );
printf("pos=%e, %e, %e, %e\n", pos[0], pos[1], pos[2], pos[3] );	
//	getchar();
}
*/ 
	if(tmp_ph1<=Rout){
		printf("should not\n");	
		getchar();
	}
	if(position[1]>=1e6){
		if(position[2]>=0.0) {
			*i_path = 0; // to infinity
		}else{
			*i_path = 3; // failed geodesic
		}	
			return 0;
	}else{
			norm_vect(position, momentum);	// normalize k-vector
			*i_path = 1; // to disk
			return 0;	
	}	

}else{ 
	if(interp==-1) icc2++;
	if(interp==-2) icc3++;
	if(interp==-3) icc4++;	
//	   return 0;
	  // the interpolation is not applicable for this case
	  // do direct raytracing for this photon                                                                				
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<	
		                                                                                                 
        /* ==================== The block for ray-tracing ================================= */
	    raytrace_data rtd;
	    raytrace_prepare(bh_spin, position, momentum, NULL, 10.0, 0, &rtd);
int icase2=0;		
	do{
icase2++;
if(icase2>100000) {printf("too many raytrace in case1 which is for photon escaping meridian\n"); getchar();}
	    dl = 1e9; // use maximal step   
			vect_copy(position, posi);
			vect_copy(momentum, ki);	  
	    raytrace(position , momentum, NULL, &dl, &rtd);
			vect_copy(position, posf);
			vect_copy(momentum, kf);
	    // after raytrace(), posi and momi already changed
	    
	    		// stop condition: if applied, stop raytrace
				if ((position[1] < r_bh(bh_spin))) {
//				printf("out torus1, raytrace(): photon goes to bh event horizon\n");
				*i_path = 3;  
				return 0;  
				}
				// also stop if relative error this step is too large
				if (rtd.error>1e-2) {
//				printf("out torus2, raytrace(): aborted due to large error (generally it is due to photon going to BH)\n");
				*i_path = 3; // failure step; as photon goes to BH event horizon, 
				return 0;
				}
				
				if(isnan(position[1]) || isnan(position[2]) || isnan(position[3])){ // failed geodesic	
				*i_path = 3;  
				return 0;  
				}				
      
			if((position[1]<=1e6)&&((posi[2]>0.0)&&(posf[2]<0.0))){ // chech if photon go across the disk; outer radius is 1e6 	
//if(fabs((position[3]-pos[3])/position[3])>1e-1)
//{				     							  
			    *i_path = 1; 
				for(i=0;i<4;i++){
					position[i] = posi[i] + (-posi[2])/(posf[2]-posi[2])*(posf[i]-posi[i]);	
					momentum[i] = ki[i] + (-posi[2])/(posf[2]-posi[2])*(kf[i]-ki[i]);					
				}
				norm_vect(position, momentum); // normalize k-vector	
						
//printf("big error2\n");
//printf("position=%e, %e, %e, %e\n", position[0], position[1], position[2], position[3] );
//printf("pos=%e, %e, %e, %e\n", pos[0], pos[1], pos[2], pos[3] );	
//	getchar();
//}
			    return 0;
			}
			
			if(position[1]>1e9){
//if(fabs((position[3]-pos[3])/position[3])>1e-1)
//{
//printf("big error3\n");
//printf("position=%e, %e, %e, %e\n", position[0], position[1], position[2], position[3] );
//printf("pos=%e, %e, %e, %e\n", pos[0], pos[1], pos[2], pos[3] );	
//	getchar();
//}			
				if(position[2]>=0.0){
					 *i_path = 0; // to infinity; if pos[2]<=0, ignore this photon
				}else{
					 *i_path = 3;
				}	 	 			    
			    return 0;
			}					  	
	}while(1);		
		/* ==================== End of the block for ray-tracing ================================= */  		
}		
return 0;
}

// probability distribution function of emitted photons from disk
double disk_photons_pdf(double x)
{
    double hardf = 1.0;
    double cos_mu = -1.0;
    double r = exp(x);
    double gtt = -1. + 2./r;
    double gtf = -2.*bh_spin/r;
    double gff = sqr(r) + sqr(bh_spin) + 2.*sqr(bh_spin)/r;
    double Omega = 1./(bh_spin + pow(r,1.5));
    double U_t = sqrt(-1.0/(gtt + 2.*Omega*gtf + sqr(Omega)*gff)) * (gtt + Omega*gtf);

    double T = sqrt4((-U_t)*disk_nt_flux(r)/sb_sigma);
    
    return M_PI * blackbody_photons_total(T, hardf, cos_mu) * r*r;
}


// g-factor for keplerian disk
double disk_gfactor_K(double r, double a, double l, double q, double rms)
//*********************************************************
// g-factor at given point
// relativistic correction to radiation output that includes
// Doppler effect and gravitational redshift
{
	if (r<=rms) return(0.0);
	double OmegaK = 1./(a + pow(r,1.5));
	double r2 = sqr(r);
	double a2 = sqr(a);
    double E = sqrt(1. - 2./r * sqr(1.-a*OmegaK) - (r2+a2)*sqr(OmegaK));
    return E / (1. - OmegaK*l); 
}

/* convert the axis of the torus (0, 0, 1) in torus frame to Cartesian coordiante of BH, 
   according to (prectheta, precphi) in order to check if photon is inside of the torus 
   by means of Dot Product between unit vector and photon position vector*/

    torus_axis1()
{
	double prec_theta, prec_phi;
	double posi[3], posf[3];
	double R[3][3];
	
	prec_theta = prectheta;
	prec_phi = precphi;

    /* initial unit vector of torus axis in the torus frame */
    posi[0] = 0.;
    posi[1] = 0.;
    posi[2] = 1.; 

    invrotmatrix(R,prec_phi,prec_theta);
    multmatvec(posf,R,posi);

	axis_vec[0] = posf[0];    // x-component of the unit vector of torus axis
	axis_vec[1] = posf[1];    // y-component of the unit vector of torus axis
	axis_vec[2] = posf[2];    // z-component of the unit vector of torus axis
	
    return 0;
	
}

/* convert the axis of the torus (0, 0, 1) in torus frame to Cartesian coordiante of BH, 
   according to (prectheta, precphi) in order to check if photon is inside of the torus 
   by means of Dot Product between unit vector and photon position vector*/
    torus_axis2()
{
	double prec_theta, prec_phi;
	double posi[3], posf[3], posm[3];
	double R[3][3];
	
	prec_theta = prectheta;
	prec_phi = precphi;

    /* initial unit vector of torus axis in the torus frame */
    posi[0] = 0.;
    posi[1] = 0.;
    posi[2] = 1.; 

    invrotmatrix(R,prec_phi,prec_theta);
    multmatvec(posm,R,posi);

    invrotmatrix(R,0.0,prec_theta);
    multmatvec(posf,R,posm);

	axis_vec[0] = posf[0];    // x-component of the unit vector of torus axis
	axis_vec[1] = posf[1];    // y-component of the unit vector of torus axis
	axis_vec[2] = posf[2];    // z-component of the unit vector of torus axis
	
    return 0;
	
}


/* convert photon position from B-L spherical coordinate to cartesian coordinate */
	bl2cart(position_bl, position_cart)
	double position_bl[], position_cart[];
{
	position_cart[0] = position_bl[0];
	position_cart[1] = position_bl[1]*fabs(sqrt(1.-sqr(position_bl[2])))*cos(position_bl[3]);   // x
	position_cart[2] = position_bl[1]*fabs(sqrt(1.-sqr(position_bl[2])))*sin(position_bl[3]);   // y
	position_cart[3] = position_bl[1]*position_bl[2];                                           // z
	
	return 0;
	
}

/* check if photon is inside torus  */
int inside_dotprod(position_bl)
    double position_bl[];
{
	double  position_cart[4];
	int inside_dot, inside_r, inside_theta;
	
	inside_dot = 0;

	bl2cart(position_bl, position_cart);

	double dot_prod = axis_vec[0]*position_cart[1] + axis_vec[1]*position_cart[2] + axis_vec[2]*position_cart[3];
	
	inside_r = ((position_bl[1]>Rin) && (position_bl[1]<Rout)) ; // to chick if photon is inside torus in r-direction
	inside_theta = (fabs(dot_prod/position_bl[1]) <= fabs(sin(torus_theta0))); // to chick if photon is inside torus in theta-direction
	inside_dot = inside_r*inside_theta;
		
	return inside_dot; // 1: inside of torus; 0: outside of torus

}


int prod_pU(position, momentum, gf) /* dot product between photon momentum and medium 4-velocity */
double position[], momentum[], *gf;
{
		double U[4],Omega,ell_torus;

        sim5metric m;
                
        // get metric and tetrad
        kerr_metric(bh_spin, position[1], position[2], &m);
                  
		ell_torus = ellK(position[1], bh_spin); // torus specific angular momentum
            
		Omega = Omega_from_ell(ell_torus, &m);
//		double Omegak = OmegaK(position[1],bh_spin); // for accretion disk
                      
        // dot product of photon direction and medium motion: p[4]*U[4]
        fourvelocity_azimuthal(Omega, &m, U);        
		// gravitational redshift factor with respect to infinity
		*gf = (momentum[0]*m.g00 + momentum[3]*m.g03)/dotprod(momentum, U, &m);
//        printf("pU = %e\n", dotprod(momentum, U, &m));

	return 0;
}        

int bl2torus(position_bl, momentum_bl, momentum_torus) 
/* convert photons momentum from B-L to torus-rest frame */
double position_bl[], momentum_bl[], momentum_torus[];
{
//	double momentum_loc[4], diri[3], dirf[3], dirm[3];
	    
	/* first convert momentum from B-L to torus-rest frame */
	
        sim5tetrad t; 
        sim5metric m;
                
        // get metric and tetrad
        kerr_metric(bh_spin, position_bl[1], position_bl[2], &m);
                   
		double ell_torus = ellK(position_bl[1], bh_spin); // torus specific angular momentum
        
        double Omega;                  
		Omega = Omega_from_ell(ell_torus, &m);
        //double Omega = OmegaK(r,bh_spin); // for accretion disk
        
        tetrad_azimuthal(&m, Omega, &t); 
//        printf("===bl2torus:\n");
//        printf("%e, %e, %e, %e\n", momentum_bl[0], momentum_bl[1], momentum_bl[2], momentum_bl[3]);
//        printf("just see null %e\n", dotprod(momentum_bl, momentum_bl, &m));
		bl2on(momentum_bl, momentum_torus, &t);
//		printf("%e, %e, %e, %e\n", momentum_torus[0], momentum_torus[1], momentum_torus[2], momentum_torus[3]);
//		printf("just see null %e\n", -sqr(momentum_torus[0])+sqr(momentum_torus[1])+sqr(momentum_torus[2])+sqr(momentum_torus[3]));
//		getchar();
	return 0;
}

int torus2bl(position_bl, momentum_torus, momentum_bl) 
/* convert photons momentum from torus frame to B-L */
double position_bl[], momentum_bl[], momentum_torus[];
{
//	double momentum_loc[4], diri[3], dirf[3], dirm[3];
		
	/* second, convert momentum from torus-rest frame to B-L */
        sim5tetrad t; 
        sim5metric m;
                
        // get metric and tetrad
        kerr_metric(bh_spin, position_bl[1], position_bl[2], &m);
                    
		double ell_torus = ellK(position_bl[1], bh_spin); // torus specific angular momentum
        
        double Omega;            
		Omega = Omega_from_ell(ell_torus, &m);
        //double Omega = OmegaK(r,bh_spin); // for accretion disk
        
        tetrad_azimuthal(&m, Omega, &t);              
//		on2bl(momentum_torus, momentum_bl, &t);	
//		printf("===torus2bl:\n");
//		printf("%e, %e, %e, %e\n", momentum_torus[0], momentum_torus[1], momentum_torus[2], momentum_torus[3]);
//		printf("just see null %e\n", -sqr(momentum_torus[0])+sqr(momentum_torus[1])+sqr(momentum_torus[2])+sqr(momentum_torus[3]));
		on2bl(momentum_torus, momentum_bl, &t);	
//		printf("%e, %e, %e, %e\n", momentum_bl[0], momentum_bl[1], momentum_bl[2], momentum_bl[3]);
//        printf("just see null %e\n", dotprod(momentum_bl, momentum_bl, &m));
//		getchar();
		
	return 0;
}
   
double ray_tau(double l_i, double energy) // calculate the depth of each step in ray-trace
{
	double dtau, sigm_mean;
    static int ins=1, inc=1;


	/* calculate depth */
     c_hunt(sig_e,nop_se,energy,&ins);
     sigm_mean = sigmm[ins]+(energy-sig_e[ins])*(sigmm[ins+1]-sigmm[ins])/
                                                (sig_e[ins+1]-sig_e[ins]);
   /* l_i is now in Rg !! */
     dtau = l_i*tTh1Rg*sigm_mean;  
//   printf("li = %e, tTh1Rg = %e, sigm_mean = %e\n", l_i, tTh1Rg, sigm_mean);
//   printf("energy = %e, tau = %e\n", energy, l_i*tTh1Rg*sigm_mean);
//   getchar();

   return dtau;
}

int raytrace_in_torus(position,momentum,energy)  
/* assuming photon is inside torus, given position[] and momentum[], to determine:
1. accumulated escaping probalibity along geodesic to boundary
2. if return -1, raytrace failure	
   if return 1, escaping
   if return 0, scattering 											
*/																                                                     
double position[],momentum[];
double energy;
{
double  posi_orig[4], momi_orig[4], momi_torus[4];
double  dl, p_sc, p_esc, dtau, tau;
int inside_torus;
double energy1, energy2, energy_m;
double gf; // gravitational redshift factor
                
/* =============== repeat raytrace to determine the scattering position, if photon is scattered =================*/
    double posmax[4],posmin[4],
		   mommax[4],mommin[4],
		   pmin, pmax;
	int i;
	double x;	
	x = urand; 
	if(x==0.0) x+=0.0001;
	if(x==1.0) x-=0.0001; // ensure that x = (0,1), not [0,1]

        /* copy vectors to another temporary ones which are used in raytrace */
        vect_copy(position, posi_orig); // posi -> posi_orig
        vect_copy(momentum, momi_orig); // momi -> momi_orig
                                 
        /* ==================== The block for ray-tracing ================================= */    
	    raytrace_data rtd;           
	    raytrace_prepare(bh_spin, posi_orig, momi_orig, NULL, 1.0, 0, &rtd);

		tau  = 0.0;
		dtau = 0.0;
		p_sc = 0.0; // accumulated scattering probalibity: p_sc = 1 - exp(-tau), where tau += sigma*dl 
		p_esc = 0.0;
		
		prod_pU(posi_orig, momi_orig, &gf); 
		energy1 = energy/gf; /* calculate the photon energy once it enters into the torus 
								'energy' is not changed in this subroutine */		
		energy2 = energy1; /* energy in torus-rest frame */
int iray=0;
		do{ // loop to save each step along geodesic						
iray++;
if(iray>100000){printf("too many raytrace in torus\n"); getchar();}
			dl = 1e6; // have to use maximal step; rather than 0						
									
			vect_copy(posi_orig, posmin);
			vect_copy(momi_orig, mommin);
			pmin = p_sc;
												
			raytrace(posi_orig , momi_orig, NULL, &dl, &rtd); 
				// stop if relative error this step is too large
				if (rtd.error>1e-2) {
				printf("raytrace(): aborted due to large error\n");
				return -1;
				}			
			
			vect_copy(posi_orig, posmax);
			vect_copy(momi_orig, mommax);

			inside_torus = inside_dotprod(posi_orig);	
		  if(inside_torus){	

			prod_pU(posi_orig, momi_orig, &gf);
			energy2 = energy/gf; // photon energy (in torus-rest frame) changed after each raytrace step
			
			energy_m = 0.5*(energy1 + energy2);
			dtau = ray_tau(dl, energy_m); 
			energy1 = energy2;	// energy1 as the initial energy for the next step	  
			
			tau += dtau; // accumulate optical depth before escaping torus
			p_sc = 1.0 - exp(-tau);
			pmax = p_sc;
			if(p_sc>x) break; // condition to stop raytrace
		  }else{
			if(iray==1) return -1; // just pass the torus  			  
			if(((posmin[1]-Rout)*(posmax[1]-Rout))<0.0){// escaping the torus from the outer boundary 								
				for (i=0; i<4; i++){
					position[i] = posmin[i] + (Rout-posmin[1])*(posmax[i]-posmin[i])/
												   (posmax[1]-posmin[1]); 
					momentum[i] = mommin[i] + (Rout-posmin[1])*(mommax[i]-mommin[i])/
												   (posmax[1]-posmin[1]); 	
				}				  	
				position[1] = Rout;							   					    
				return 1;				
			}else{	  
				vect_copy(posi_orig, position);
				vect_copy(momi_orig, momentum);				
				return 1;
			}
		  }
		}while(1);	
		/* ==================== The block for ray-tracing ================================= */
					/* do interpolation to determine new pos[], mom[], energy and weight, corresponding to x */
//					printf("interpolation:\n");
					for (i=0; i<4; i++){
					position[i] = posmin[i] + (x-pmin)*(posmax[i]-posmin[i])/
												   (pmax-pmin); 
					momentum[i] = mommin[i] + (x-pmin)*(mommax[i]-mommin[i])/
												   (pmax-pmin); 						    
					}
/* =============== end of determining the scattering position, if photon is scattered =================*/

return 0;
}
//==============================================================================================

/* 链表初始化 */
int list_init(void)
{
    pHead = malloc(sizeof(LIST_NODE));
    pHead->pNext = NULL;
    pTail = pHead;
    return 0;
}

/* 链表添加节点 */
int list_node_add(position, momentum, prob, el)
double position[], momentum[], prob, el; // photon position, momentum, energy and accumulated probability
{
	int i; 
	
    for (i=0; i<4; i++) pTail->pos[i] = position[i];	
    for (i=0; i<4; i++) pTail->mom[i] = momentum[i];
    pTail->pp = prob;
    pTail->ee = el;
    
    pTail->pNext = malloc(sizeof(LIST_NODE));
    pTail = pTail->pNext;
    pTail->pNext = NULL;
    return 0;
}

int list_free(void)  
{  
    LIST_NODE *pt = NULL;  
  
    while (pHead != NULL)  
    {   
            pt = pHead;  
            pHead = pt->pNext;
            free(pt);
    }
    return 0;
}

int scattering_gr(position, momentum, E_loc) // scattering subroutine 
double position[], momentum[], *E_loc;
{
	double momentum_torus[4];
					
        sim5tetrad t; 
        sim5metric m;
                
        // get metric and tetrad
        kerr_metric(bh_spin, position[1], position[2], &m);                   
		double ell_torus = ellK(position[1], bh_spin); // torus specific angular momentum                      
		double Omega = Omega_from_ell(ell_torus, &m);
        //double Omega = OmegaK(r,bh_spin); // for accretion disk        
        tetrad_azimuthal(&m, Omega, &t); 
        
        // ------------------------------------------------	
        /* convert photons momentum from B-L to torus-rest frame */	    
		bl2on(momentum, momentum_torus, &t);
			
		double w_cs, energy;	
		energy = (*E_loc);	
		momentum_norm_to_null(momentum_torus); // normalize the photon momentum so that k = [1, kx, ky, kz]
	   	double mom[3];
	   	mom[0] = momentum_torus[1]; mom[1] = momentum_torus[3]; mom[2] = momentum_torus[2];
	   	// have to rearrange the orders of these three Cartesian components; 
	   	// As the definitions of Piotr's and Michal's are different

	   	scattering(position, mom, &energy, w_cs); // Actually position[] is not used in the scattering; 
	   	momentum_torus[1] = mom[0]; momentum_torus[3] = mom[1]; momentum_torus[2] = mom[2];
	   	// rearrange back

		*E_loc = energy;
		momentum_norm_to_null(momentum_torus);
		// normalization to ensure that new momentum is null	

		on2bl(momentum_torus, momentum, &t); // convert photon momentum back to B-L frame
	
	return 0;
}

int null_norm_to_momentum(double k[], double energy)
{
	double kx = k[1]/k[0];
	double ky = k[2]/k[0];
	double kz = k[3]/k[0];
	double kk = sqrt(kx*kx + ky*ky + kz*kz);
	kx /= kk; ky /= kk; kz /= kk; 		
	k[1] = energy*(kx);
	k[2] = energy*(ky);
	k[3] = energy*(kz);
	k[0] = energy;
	
	return 0;
}

int momentum_norm_to_null(double k[])
{
	double kx = k[1]/k[0];
	double ky = k[2]/k[0];
	double kz = k[3]/k[0];
	double kk = sqrt(kx*kx + ky*ky + kz*kz);
	kx /= kk; ky /= kk; kz /= kk; 	
	k[1] = kx;
	k[2] = ky;
	k[3] = kz;
	k[0] = 1.0;
	return 0;
}

int raytrace_esc_torus(position_esc, momentum_esc, E_inf)
double position_esc[], momentum_esc[], E_inf;
{
	int i_path;/* photon trajectory by raytrace:
				0: to infinity; 1: to disk (for reflection); 2: to torus; 3: to BH */
	double gf; // gravitational redshift factor
	double E_loc, weight_esc = 1.0;
	int ir;
//printf("here5\n");	
//return 0;
	raytrace_out_torus(position_esc, momentum_esc, &i_path); /* determine whether escaping photon goes to infinity, or reflection */
//printf("here6: i_path = %d\n", i_path);	
//return 0;
	if(i_path==3) return -1;	// failed raytrace, going on to next scattering of this photon	

	if(i_path==0){// from torus to infinity
			write_spec(position_esc, momentum_esc, E_inf);	
	} 			                                                                                                

	if(i_path==1){ // this is the reflection
			prod_pU(position_esc, momentum_esc, &gf); // After doing raytrace, calculating "g_f" in order to determine the photon energy at the next position
			E_loc = E_inf/gf;
		for(ir=0; ir<refl_rept; ir++)	
			disk_reflection(position_esc, momentum_esc, E_loc, weight_esc);
	}	
	if(i_path==2) return 1; // photon enter into the torus again 
return 0;
}

int disk_reflection(position, momentum, ee, ww)
double position[], momentum[], ee, ww;
{
  double mom_i[3],mom_f[3],pos_l[3],
         energy,weight,weight_min,mu_prime,sine2,term,bound;	 
  int i, j;
  double pos[4], mom_r[4];

  vect_copy(position, pos);
  
  energy = ee;
  weight = ww;
    
  // both postion and direction have to be in cartesian coordinate !!! 
  bl2torus(pos, momentum, mom_r); // convert the B-L direction to the disk-rest frame 
//  printf("mom[] = %e %e %e %e\n", mom_r[0], mom_r[1], mom_r[2], mom_r[3]);
/* NOTE: rerange the order of the x,y,z direction in order to be consistent with mom_i[] which is differently defined by Piotr */
  mom_i[0] = mom_r[1]/mom_r[0]; 
  mom_i[1] = mom_r[3]/mom_r[0];
  mom_i[2] = mom_r[2]/mom_r[0];
//  for (i=0; i<3; i++) mom_i[i] = mom_r[i+1]/mom_r[0]; // normalization, make sure that kx*kx + ky*ky + kz*kz = 1

  pos_l[0] = pos[1]*fabs(sqrt(1.0-sqr(pos[2])))*cos(pos[3]);
  pos_l[1] = pos[1]*fabs(sqrt(1.0-sqr(pos[2])))*sin(pos[3]);
  pos_l[2] = pos[1]*pos[2];	
  
//  printf("input:\n");
//  printf("energy = %e\n", energy);			
//  printf("x = %e, y = %e, z = %e\n", pos_l[0], pos_l[1], pos_l[2]);
//  printf("kx = %e, ky = %e, kz = %e\n", mom_i[0], mom_i[1], mom_i[2]);
//  printf("k.k = %e\n", sqr(mom_i[0]) + sqr(mom_i[1]) + sqr(mom_i[2]));
//  getchar();
	
  double r0 = r_bh(bh_spin);
  weight_min = weight*1e-5;
//  if (weight_min < wghtmin) weight_min = wghtmin;

double pos_l0[3], mom_i0[3], energy0;
// individual photon
  int i_esc=0, iscat = 0, iline = 0;
  do{
	iscat++;
	if(iscat>100){printf("too many scattering in the disk\n"); break;}  
	
	i_esc = raytrace_in_disk(pos, pos_l, mom_i, energy, &weight, 0);
    if(i_esc==0){   
		do{
		mu_prime = 1-(exp(log(1+2*energy)*rnd())-1)/energy;
		sine2 = 1-mu_prime*mu_prime;
		term  = 1/(1+energy*(1-mu_prime));
		bound = term*(term+1/term-sine2)/2;
		} while (rnd() > bound);
		transform(mom_i,mom_f,sqrt(sine2),mu_prime);
		for (i=0; i<3; i++) mom_i[i] = mom_f[i];        
		energy *= term;
	}
	if(i_esc==-1) break; // photon is absorbed in the disk
	if(i_esc==1) break;	// photon escape from the disk 
	if(i_esc==2){ // fluorescent photon is emitted
		energy = 6.4/511.0; // Fe Ka line energy
		energy0 = energy;
		vect_copy(pos_l, pos_l0); vect_copy(mom_i, mom_i0);  
		iline = 1;
		break;
	}	
  }while(1);

if(iline==1){  
//printf("here\n"); getchar();	
for(j=0; j<line_rept; j++){
	
  energy = energy0;
  vect_copy(pos_l0, pos_l); vect_copy(mom_i0, mom_i);  
			
  i_esc=0, iscat = 0;
  do{
	iscat++;
	if(iscat>100){printf("too many scattering in the disk\n"); break;}  
	
	i_esc = raytrace_in_disk(pos, pos_l, mom_i, energy, &weight, 1);
    if(i_esc==0){   
		do{
		mu_prime = 1-(exp(log(1+2*energy)*rnd())-1)/energy;
		sine2 = 1-mu_prime*mu_prime;
		term  = 1/(1+energy*(1-mu_prime));
		bound = term*(term+1/term-sine2)/2;
		} while (rnd() > bound);
		transform(mom_i,mom_f,sqrt(sine2),mu_prime);
		for (i=0; i<3; i++) mom_i[i] = mom_f[i];        
		energy *= term;
	}
	if(i_esc==-1) break; // photon is absorbed in the disk
	if(i_esc==1) break;	// photon escape from the disk 
	if(i_esc==2){ // fluorescent photon is emitted
		energy = 6.4/511.0; // Fe Ka line energy 
	}	
  }while(1);
//printf("here\n"); getchar();   
} 
//printf("here\n"); getchar(); 
}  

return 0;
}

int raytrace_in_disk(posi, pos, mom, energy, weight, status)
// pos[], mom[] are 3-vector in cartesian coordinate, in disk-rest frame
double posi[], pos[], mom[], energy, *weight;
int status; // 0 is for continuum photon; 1 is for fluorescent photon
{
  double sigma_es,sigma_abs,sigma_tot,p_esc,wesc,lambda,d,fi;
  int    indx,iindx,i,findx;
	
  sigma_es = sigma_KN(energy);
//if(fabs(energy-6.4/511.0)<1e-5){
//  sigma_abs = 2.53*sigma_es;
//}else{	  
  sigma_abs= sigma_bf(energy);
//}  
  sigma_tot= sigma_es+sigma_abs;

  double iref = rnd();          
  if (mom[2] <= 0) { /* downwards - no escape */
    p_esc = 0;
  } else {
    d = fabs(pos[2]/mom[2]);
    p_esc = exp(-d*sigma_tot); 
    
    if(iref<=p_esc) { // photon escape from the disk
		wesc = (*weight);
		int feline = status;
		raytrace_infinity(posi, mom, energy, wesc, feline); // assume that the photon enters and is reflected 
													// at the same position on the disk: posi[], rather than pos[]	
		return 1; 													
	}																						
  }
  
  lambda = -log(1-(1-p_esc)*iref)/sigma_tot;
  
  for (i=0; i<3; i++)
    pos[i] += lambda*mom[i];  /* position of the next scattering */  
						
double abs_scat = rnd();
if(abs_scat>(sigma_es/sigma_tot)){ // photon is absorbed 
	if(energy>(7.1/511.)){ // photon energy is larger than k edge energy E_k = 7.1
		double sigma_fe = 3.8e-20*3.3e-5*1.21/6.65e-25; // sig_fe is criss-section per atom for absorption by neutral iron  !!! NOT per hydrogen
														// sigma_fe = sig_fe*Az, cross-section per hydrogen atom.
														// So need to multiply it by 1.21, and in units of sigma(Thomson)
		double p_fe = sigma_fe/sigma_tot;
//		printf("P_fe = %e\n", p_fe);
//		getchar();		
		double fe_rnd = rnd(); // the probability that an incident photon is absorbed by an iron atom
		if(fe_rnd < p_fe){ // yes, absorbed by an iron atom
			double p_rnd = rnd(); // probability to compare with fluorescent yield = 0.34
			if(p_rnd < 0.34){ // fluorescent photon is emitted
				return 2; // fluorescent photon proceeds in the same way as the continuum photon 
			}else{
			    return -1; // fluorescent photon is not emitted
			}		
		}else{
			return -1; // not absorbed by an iron atom
		}															
	}else{ // smaller than E_k
		return -1;
	}	 
}						
          
return 0;
}

int raytrace_infinity(pos, dir, energy, weight, iph)
// pos[] are 4-vector in spherical coordinate, dir[] are 3-vector in cartesian coordinate, in disk-rest frame
double pos[], dir[], energy, weight;
int iph; // 0 is for continuum photon; 1 is for fluorescent photon
{
	double kx=dir[0], ky=dir[1], kz=dir[2];
    // make sure [kx,ky,kz] is a unit vector
    double kk = sqrt(kx*kx + ky*ky + kz*kz);
    kx /= kk; ky /= kk; kz /= kk;

    double r = pos[1];
    double m = 0.0;
    double phi = pos[3];

    // get metric coefficients
    sim5metric metric;
    kerr_metric(bh_spin, r, m, &metric);
    
    /* Keplerian four-velocity */
    double u[4];
    double OmK = OmegaK(r,bh_spin);
	fourvelocity_azimuthal(OmK, &metric, u);
	
	double k[4], k_loc[4];
/* NOTE: rerange the order of the x,y,z direction in order to be consistent with k_loc[] which is differently defined by Michal */
	k_loc[0]=1.0;
	k_loc[1]=kx;
	k_loc[2]=kz;
	k_loc[3]=ky;
	sim5tetrad t;
    tetrad_azimuthal(&metric, OmK, &t);
    on2bl(k_loc, k, &t);
//    printf("k_loc = %e %e %e %e\n", k_loc[0], k_loc[1],k_loc[2],k_loc[3]);
//    printf("pos = %e %e %e\n", r, m, phi);
//    printf("xyz = %e %e %e\n", x, y, z);
//--------------------------------------------------------------------------------------------------
    // get photon constants of motion
    double q, l;
    photon_motion_constants(bh_spin, r, m, k, &l, &q);

    // given the position and motion constants, find the geodesic parameters
    int status;
    geodesic gd;
    geodesic_init_src(bh_spin, r, m, k, k[1]<0?1:0, &gd, &status);
    if ((!status) || (r < gd.rp)) {
        fprintf(stderr, "WRN: geodesic_init_src failed (status=%d, r=%.2e, rp=%.2e)\n", status, r, gd.rp);
	printf("2\n");
        return -1;
    }
    
    double gd_P = geodesic_P_int(&gd, r, k[1]<0?1:0);

    // sometime a trajectory is found, which is correct, but unphysical
    // (photon would have pass radial turning point, which lies bellow horizon)
    if ((creal(gd.r1) < r_bh(bh_spin)) && (gd_P > gd.Rpa)) {
        // unphysical solution - geodesic goes to observer through BH
        return -1;
    }
    if(((gd_P > gd.Rpa) && (gd.nrr == 2)) || (gd.m2p==1.0)) return -1;
    double dphi = geodesic_position_azm(&gd, r, m, gd_P);
    if(isnan(dphi)) return -1;    

    // return theta_infinity angle and g-factor
    double cos_inc, phi_inf, g_inf;
    cos_inc = gd.cos_i;
    phi_inf = rad2deg(reduce_angle_2pi(phi + dphi));
    g_inf = (k[0]*metric.g00 + k[3]*metric.g03) / dotprod(k, u, &metric);
    energy*=g_inf;
//	printf("inc = %e, phi = %e, g = %e\n",cos_inc, phi_inf, g_inf);
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
		/* bins in \phi: 0-20, 90-110, 180-200, 270-290 degs */
	    int findx = -1; 
	    int i, indx, iindx;
	    for(i=0; i<MAXphi; i++) {
	      if (fabs(phi_inf-(i*90+10)) <= 10) {  
	       findx = phi_inf/90.;
	      }
	    }
	    
	if((cos_inc>0.0)&&(findx>=0)){ // only deal with photon energy in the interesting range 
		iindx = cos_inc*nangle; 	
		if(iph==0){
		indx = num_bin_o*log10(energy/emin_o); // indx corresponding to continuum photon energy at infinity 
			    Nrefl[indx][iindx][findx] += weight;	
		}
		if(iph==1){
			if((energy<elmin_o)||(energy>elmax_o)) {printf("out of energy range\n"); getchar();}
		indx = num_bin_line*(energy-elmin_o)/(elmax_o-elmin_o); // indx corresponding to line photon energy at infinity 	
				Nline[indx][iindx][findx] += weight;		
		}	               
	}    
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return 0;
}

sphe2cart(position, v_sphe, v_cart)
double position[], v_sphe[], v_cart[];
{
			/* x-component */
			v_cart[0] = cos(position[3])*fabs(sqrt(1.-sqr(position[2])))*v_sphe[0] 
						+ position[1]*cos(position[3])*position[2]*v_sphe[1]
						- position[1]*sin(position[3])*fabs(sqrt(1.-sqr(position[2])))*v_sphe[2];
			
			/* y-component */			
			v_cart[1] = sin(position[3])*fabs(sqrt(1.-sqr(position[2])))*v_sphe[0] 
						+ position[1]*sin(position[3])*position[2]*v_sphe[1]
						+ position[1]*cos(position[3])*fabs(sqrt(1.-sqr(position[2])))*v_sphe[2];
			
			/* z-component */			
			v_cart[2] = position[2]*v_sphe[0] - position[1]*fabs(sqrt(1.-sqr(position[2])))*v_sphe[1];			
  
return 0;	
}

cart2sphe(position, v_cart, v_sphe)
double position[], v_cart[], v_sphe[];
{
		double r = sqrt(sqr(position[1])+sqr(position[2])+sqr(position[3]));
		double xy = sqr(position[1])+sqr(position[2]);
			/* r-component */
			v_sphe[0] = (position[1]*v_cart[0] + position[2]*v_cart[1] + position[3]*v_cart[2])/r;
			
			/* theta-component */		
			v_sphe[1] = (position[3]*(position[1]*v_cart[0] + position[2]*v_cart[1]) - xy*v_cart[2])/sqr(r)/sqrt(xy); 		
			
			/* phi-component */	
			v_sphe[2] = (position[1]*v_cart[1] - position[2]*v_cart[0])/xy;					
  
return 0;	
}

transform(omega_i,omega_f,st,ct)
    double omega_i[],omega_f[],
           st,ct;
{
    double x,y,r2,o3,dd,cph,sph,norm;
    int    i;

    do {
      x = rnd()-0.5;
      y = rnd()-0.5;
      r2 = x*x+y*y;
    }
    while (r2 > 0.25);
    o3 = 1 - omega_i[2]*omega_i[2];
    dd = sqrt(o3*r2);
    cph = x/dd;
    sph = y/dd;
    omega_f[0] = st*( omega_i[1]*sph-omega_i[0]*omega_i[2]*cph)+omega_i[0]*ct;
    omega_f[1] = st*(-omega_i[0]*sph-omega_i[1]*omega_i[2]*cph)+omega_i[1]*ct;
    omega_f[2] = st*cph*o3+omega_i[2]*ct;

    norm = 1/sqrt(omega_f[0]*omega_f[0]+omega_f[1]*omega_f[1]+omega_f[2]*omega_f
[2]);
    if (fabs(norm-1.0) > 0.001)
    {
      for (i=0; i<3; i++)
         omega_f[i] *= norm;
    }

    return 0;
}


double sigma_bf(e)
   double e;
{
   double  sigma;

   if (e >= elop) {
     sigma = coop/(e*e*e);
   } else {
     if (e > enop[0]) {
       c_hunt(enop,MAXA-1,e,&jlo_op);
       sigma = opac[jlo_op]+(e-enop[jlo_op])*(opac[jlo_op+1]-opac[jlo_op])/
	                                     (enop[jlo_op+1]-enop[jlo_op]);
     } else {
       sigma = opac[0];
     }
   }

   return(sigma);
}

double sigma_KN(e)
     double e;
{
  double  sigma,z;
/*
   Klein-Nishina cross section in units of sigma(Thomson),
   e is energy in units of the electron rest mass
   Per hydrogen atom! It means multiplying it by 1.21, which is the
   number of electrons per H atom for cosmic abundances

*/  
  z = 1+2*e;
  if (e < 0.02)
    sigma = (1+e*(2+e*(1.2-0.5*e)))/z/z;  /*  Hubbell 1980  */
  else
    sigma = 0.375/e*((1-2/e-2/e/e)*log(z)+0.5+4/e-0.5/z/z);

  sigma *= 1.21; 
  return(sigma);
}

read_opacity()
{
   double ee,dum,op;
   int    i;
   FILE  *fp,*fopen();

   if ((fp=fopen("absorp.dat","r")) != NULL)
   {
      for (i=0; i<MAXA; i++)
      {
        fscanf(fp,"%lf %lf %lf",&ee,&op,&dum);
  	enop[i] = ee/511;
        opac[i] = op;
      }
      fclose(fp);

      elop = enop[MAXA-1];
      coop = opac[MAXA-1]*pow(elop,3.0);
      /*      printf("opac: %.3e  %.4e \n",elop*511,coop); */
   }
   else
   {
      printf("Cannot find file with opacity data \n");
      exit(13);
   }

   return 0;
}

int interp4d(x,y,z,t,xa,ya,za,ta,f,interp)
double x,y,z,t;
float xa[],ya[],za[],ta[],f[nrd+1][nfd+1][nth+1][nph+1];
double *interp;	   
{
	
if((x<=xa[0])||(x>=xa[nrd])) return -1;	
if((y<=ya[0])||(y>=ya[nfd])) return -1;	
if((z<=za[0])||(z>=za[nth])) return -1;
if((t<=ta[0])||(t>=ta[nph])) return -1;
	int ix,iy,iz,it;
	c_table_hunt(xa, nrd, x, &ix);
	c_table_hunt(ya, nfd, y, &iy);
	c_table_hunt(za, nth, z, &iz);
	c_table_hunt(ta, nph, t, &it);
/*if((ix==0)||(ix==nr)) return -1;	
if((iy==0)||(iy==nf)) return -1;
if((iz==0)||(iz==nth)) return -1;
if((it==0)||(it==nph)) return -1;
*/ 
	
//	printf("%d %d %d %d\n", ix,iy,iz,it);
//	getchar();	
int i, j, k,l;
double dd;

	double dx = (x-xa[ix])/(xa[ix+1]-xa[ix]);
	double dy = (y-ya[iy])/(ya[iy+1]-ya[iy]);
	double dz = (z-za[iz])/(za[iz+1]-za[iz]);
	double dt = (t-ta[it])/(ta[it+1]-ta[it]);			

	double c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15;
	c0=f[ix][iy][iz][it];
	c1=f[ix+1][iy][iz][it]-f[ix][iy][iz][it];
	c2=f[ix][iy+1][iz][it]-f[ix][iy][iz][it];
	c3=f[ix][iy][iz+1][it]-f[ix][iy][iz][it];
	c4=f[ix][iy][iz][it+1]-f[ix][iy][iz][it];
	
	c5=f[ix+1][iy+1][iz][it]-f[ix][iy+1][iz][it]-f[ix+1][iy][iz][it]+f[ix][iy][iz][it];
	c6=f[ix+1][iy][iz+1][it]-f[ix][iy][iz+1][it]-f[ix+1][iy][iz][it]+f[ix][iy][iz][it];
	c7=f[ix+1][iy][iz][it+1]-f[ix][iy][iz][it+1]-f[ix+1][iy][iz][it]+f[ix][iy][iz][it];
	c8=f[ix][iy+1][iz+1][it]-f[ix][iy][iz+1][it]-f[ix][iy+1][iz][it]+f[ix][iy][iz][it];
	c9=f[ix][iy+1][iz][it+1]-f[ix][iy][iz][it+1]-f[ix][iy+1][iz][it]+f[ix][iy][iz][it];
   c10=f[ix][iy][iz+1][it+1]-f[ix][iy][iz][it+1]-f[ix][iy][iz+1][it]+f[ix][iy][iz][it];
   
   c11=f[ix+1][iy+1][iz+1][it]-f[ix][iy+1][iz+1][it]-f[ix+1][iy][iz+1][it]-f[ix+1][iy+1][iz][it]
      +f[ix][iy][iz+1][it]+f[ix+1][iy][iz][it]+f[ix][iy+1][iz][it]-f[ix][iy][iz][it]; 
        				
   c12=f[ix][iy+1][iz+1][it+1]-f[ix][iy][iz+1][it+1]-f[ix][iy+1][iz][it+1]-f[ix][iy+1][iz+1][it]
      +f[ix][iy][iz][it+1]+f[ix][iy+1][iz][it]+f[ix][iy][iz+1][it]-f[ix][iy][iz][it]; 
      
   c13=f[ix+1][iy][iz+1][it+1]-f[ix][iy][iz+1][it+1]-f[ix+1][iy][iz][it+1]-f[ix+1][iy][iz+1][it]
      +f[ix][iy][iz][it+1]+f[ix+1][iy][iz][it]+f[ix][iy][iz+1][it]-f[ix][iy][iz][it]; 
          
   c14=f[ix+1][iy+1][iz][it+1]-f[ix][iy+1][iz][it+1]-f[ix+1][iy][iz][it+1]-f[ix+1][iy+1][iz][it]
      +f[ix][iy][iz][it+1]+f[ix+1][iy][iz][it]+f[ix][iy+1][iz][it]-f[ix][iy][iz][it];   
      
c15=f[ix+1][iy+1][iz+1][it+1]-f[ix][iy+1][iz+1][it+1]-f[ix+1][iy][iz+1][it+1]-f[ix+1][iy+1][iz][it+1]-f[ix+1][iy][iz][it]
      +f[ix][iy][iz+1][it+1]+f[ix][iy+1][iz][it+1]+f[ix][iy+1][iz+1][it]+f[ix+1][iy][iz][it+1]+f[ix+1][iy][iz+1][it]+f[ix+1][iy+1][iz][it]
      -f[ix][iy][iz][it+1]-f[ix+1][iy][iz][it]-f[ix][iy+1][iz][it]-f[ix][iy][iz+1][it]
      +f[ix][iy][iz][it];
      
      *interp = c0 + c1*dx + c2*dy + c3*dz + c4*dt
			  + c5*dx*dy + c6*dx*dz + c7*dx*dt + c8*dy*dz + c9*dy*dt + c10*dz*dt
			  + c11*dx*dy*dz + c12*dy*dz*dt + c13*dz*dt*dx + c14*dx*dy*dt
			  + c15*dx*dy*dz*dt;	
			
for(i=0; i<=1; i++){
	for(j=0; j<=1; j++){
		for(k=0; k<=1; k++){
			for(l=0; l<=1; l++){
	dd = fabs(f[ix+i][iy+j][iz+k][it+l]-f[ix][iy][iz][it]);
	if((dd>=1e4)||(f[ix+i][iy+j][iz+k][it+l]==-1.0)){ 
//	printf("dd= %e\n", dd);
//	getchar();
	return -2; // this means that the grid around f[ix][iy][iz][it] is not successive; 
						   // interpolation is not valid.
						   // should do independent raytracing for this point
	}			
			}
		}
	}
}					  
			  		  					
	return 0;
}


c_table_hunt(float xx[], int n, double x, int *jlo)
{
  int long jm,jhi,inc;
  int ascnd;

  ascnd=(xx[n] > xx[1]);
  if (*jlo <= 0 || *jlo > n) {
    *jlo=0;
    jhi=n+1;
  } else {
    inc=1;
    if (x >= xx[*jlo] == ascnd) {
      if (*jlo == n) return (*jlo);
      jhi=(*jlo)+1;
      while (x >= xx[jhi] == ascnd) {
	*jlo=jhi;
	inc += inc;
	jhi=(*jlo)+inc;
	if (jhi > n) {
	  jhi=n+1;
	  break;
	}
      }
    } else {
      if (*jlo == 1) {
	*jlo=0;
	return (*jlo);
      }
      jhi=(*jlo)--;
      while (x < xx[*jlo] == ascnd) {
	jhi=(*jlo);
	inc <<= 1;
	if (inc >= jhi) {
	  *jlo=0;
	  break;
	}
	else *jlo=jhi-inc;
      }
    }
  }
  while (jhi-(*jlo) != 1) {
    jm=(jhi+(*jlo)) >> 1;
    if (x > xx[jm] == ascnd)
      *jlo=jm;
    else
      jhi=jm;
  }

  return (*jlo);
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */

int read_sphere()
{
	double r_min = Rout;    	 // inner radius of the meridian
    double r_max = Rmax;         // outer radius of the meridian    
//---------------------------------------------------------------------------------
double r = r_min; // define the radius of the ring; the longitude of the meridian is assumed to be 0.0

int ii, ij, ik, it;
double dcos=(cossp_max-cossp_min)/npsp,
	   dkt=(ktsp_max-ktsp_min)/ntsp,
	   dkf=(kfsp_max-kfsp_min)/nfsp; 

FILE *data = fopen("sphere_table.dat", "r"); // opens new file for writing
if(data)
{
printf("Read sphere_table from file\n");
for (ii=0; ii<=npsp; ii++) {	
	for (ij=0; ij<=ntsp; ij++) {
		for(ik=0; ik<=nfsp; ik++){

			cosa[ii] = cossp_min+ii*dcos;
			kta[ij] = ktsp_min+ij*dkt; // cos(theta)
			kfa[ik] = kfsp_min+ik*dkf;
			for (it = 0; it <= 7; it++)
				fscanf(data, "%le", &sp_ph[ii][ij][ik][it]);						
		}	
	}
} 
}
fclose(data);
for (ii=0; ii<=npsp; ii++) {	
	for (ij=0; ij<=ntsp; ij++) {
		for(ik=0; ik<=nfsp; ik++){				
				sp_ph0[ii][ij][ik] = sp_ph[ii][ij][ik][0];
				sp_ph1[ii][ij][ik] = sp_ph[ii][ij][ik][1];
				sp_ph2[ii][ij][ik] = sp_ph[ii][ij][ik][2];
				sp_ph3[ii][ij][ik] = sp_ph[ii][ij][ik][3];
				sp_ph4[ii][ij][ik] = sp_ph[ii][ij][ik][4];
				sp_ph5[ii][ij][ik] = sp_ph[ii][ij][ik][5];
				sp_ph6[ii][ij][ik] = sp_ph[ii][ij][ik][6];
				sp_ph7[ii][ij][ik] = sp_ph[ii][ij][ik][7];			  	
		}	
	}
}
//printf("reading data is done\n");
//getchar();
//--------------------------------------------------------------- 
return 0;
}	


int interp3d(x,y,z,xa,ya,za,f,interp)
double x,y,z,xa[],ya[],za[],f[npsp+1][ntsp+1][nfsp+1];
double *interp;	   
{
	
if((x<=xa[0])||(x>=xa[npsp])) return -1;	
if((y<=ya[0])||(y>=ya[ntsp])) return -1;	// be aware of the max and min value
if((z<=za[0])||(z>=za[nfsp])) return -1;

	int ix,iy,iz;
	c_hunt(xa, npsp, x, &ix);
	c_hunt(ya, ntsp, y, &iy);
	c_hunt(za, nfsp, z, &iz);

/*if((ix==0)||(ix==nr)) return -1;	
if((iy==0)||(iy==nf)) return -1;
if((iz==0)||(iz==nth)) return -1;
if((it==0)||(it==nph)) return -1;
*/ 
	
//	printf("%d %d %d\n", ix,iy,iz);
//	printf("%e %e %e\n", x, y, z);	
//	printf("%e %e %e\n", xa[ix], ya[iy], za[iz]);
//	getchar();	
int i, j, k;
double dd;
		
	double dx = (x-xa[ix])/(xa[ix+1]-xa[ix]);
	double dy = (y-ya[iy])/(ya[iy+1]-ya[iy]);
	double dz = (z-za[iz])/(za[iz+1]-za[iz]);			

	double c0,c1,c2,c3,c4,c5,c6,c7;
	c0=f[ix][iy][iz];
	c1=f[ix+1][iy][iz]-f[ix][iy][iz];
	c2=f[ix][iy+1][iz]-f[ix][iy][iz];
	c3=f[ix][iy][iz+1]-f[ix][iy][iz];
	
	c4=f[ix+1][iy+1][iz]-f[ix][iy+1][iz]-f[ix+1][iy][iz]+f[ix][iy][iz];
	c5=f[ix][iy+1][iz+1]-f[ix][iy][iz+1]-f[ix][iy+1][iz]+f[ix][iy][iz];
	c6=f[ix+1][iy][iz+1]-f[ix][iy][iz+1]-f[ix+1][iy][iz]+f[ix][iy][iz];

    c7=f[ix+1][iy+1][iz+1]-f[ix][iy+1][iz+1]-f[ix+1][iy][iz+1]-f[ix+1][iy+1][iz]
      +f[ix+1][iy][iz]+f[ix][iy][iz+1]+f[ix][iy+1][iz]-f[ix][iy][iz];  
      
      *interp = c0 + c1*dx + c2*dy + c3*dz
			  + c4*dx*dy + c5*dy*dz + c6*dx*dz + c7*dx*dy*dz;
			  
for(i=0; i<=1; i++){
	for(j=0; j<=1; j++){
		for(k=0; k<=1; k++){
//			printf("%e, %e\n", f[ix][iy][iz], f[ix+i][iy+j][iz+k]);
	dd = fabs(f[ix+i][iy+j][iz+k]-f[ix][iy][iz]);
	if(f[ix+i][iy+j][iz+k]==-1.0){
/*	printf("dd= %e, %e\n", f[ix+i][iy+j][iz+k], f[ix][iy][iz]);		
	printf("ix=%d, iy=%d, iz=%d\n", ix, iy, iz);
	printf("i=%d, j=%d, k=%d\n", i, j, k);
	printf("xa=%e, ya=%e, za=%e\n", xa[ix+i], ya[iy+j], za[iz+k]);			 	
*/
	return -2; // this means that the grid around f[ix][iy][iz][it] is not successive; 
						   // interpolation is not valid.
						   // should do independent raytracing for this point
	}	
	if(dd>=1e9){
/*	printf("ix=%d, iy=%d, iz=%d\n", ix, iy, iz);
	printf("i=%d, j=%d, k=%d\n", i, j, k);
	printf("xa=%e, ya=%e, za=%e\n", xa[ix+i], ya[iy+j], za[iz+k]);
	printf("xa=%e, ya=%e, za=%e\n", xa[ix], ya[iy], za[iz]);			
		 printf("dd= %e, %e\n", f[ix+i][iy+j][iz+k], f[ix][iy][iz]);
		 getchar();
*/ 
		 return -3;			
	}
		}
	}
}			  
			  
	return 0;
}

int raytrace_check(position, momentum, pos, k)
double position[], momentum[], pos[], k[];
{		
vect_copy(position,pos);
vect_copy(momentum,k);
        /* ==================== The block for ray-tracing ================================= */
        raytrace_data rtd;
	    raytrace_prepare(bh_spin, pos, k, NULL, 1.0, 0, &rtd);
int inside_torus=0; 
	double p0i=0.0, p1i=0.0, p2i=0.0, p3i=0.0,
		   p0f=0.0, p1f=0.0, p2f=0.0, p3f=0.0,
		   k0i=0.0, k1i=0.0, k2i=0.0, k3i=0.0,
		   k0f=0.0, k1f=0.0, k2f=0.0, k3f=0.0;
	double costi, costf;	       		  		
 	do{
	    double dl = 1e9; // use maximal step  
	    
        costi = pos[2];
		p0i=pos[0], p1i=pos[1], p2i=pos[2], p3i=pos[3],
		k0i=k[0], k1i=k[1], k2i=k[2], k3i=k[3];	    
	       
	    raytrace(pos, k, NULL, &dl, &rtd);
	    // after raytrace(), posi and momi already changed
	    
        costf = pos[2];
		p0f=pos[0], p1f=pos[1], p2f=pos[2], p3f=pos[3],
		k0f=k[0], k1f=k[1], k2f=k[2], k3f=k[3];   
	    
	    		// stop condition: if applying, stop raytrace
				if (pos[1] <= r_bh(bh_spin)) {
				printf("raytrace(): going to BH\n");			
				break;
				}
								
				// also stop if relative error this step is too large
				if (rtd.error>1e-2) {
				printf("raytrace(): aborted due to large error\n");;			
				break;
				}
				
		if (pos[1] > 1e9) {		
//			printf("inf: pos = %e, %e, %e\n", pos[1], pos[2], pos[3]);	
			break;
		}
				
/*		if(pos[1]<=(1.1*Rout)){		
			inside_torus = inside_dotprod(pos);
		}  
		if(inside_torus){
			printf("torus: pos = %e, %e, %e\n", pos[1], pos[2], pos[3]);		
			break;
		}
*/		
        if((pos[1] <= 1e6)&&(costi>0.0)&&(costf<0.0)){	
				pos[0] = p0i + (-costi)/(costf-costi)*(p0f-p0i);
			    pos[1] = p1i + (-costi)/(costf-costi)*(p1f-p1i);
			    pos[2] = p2i + (-costi)/(costf-costi)*(p2f-p2i);
			    pos[3] = p3i + (-costi)/(costf-costi)*(p3f-p3i);
			    k[0] = k0i + (-costi)/(costf-costi)*(k0f-k0i);	
			    k[1] = k1i + (-costi)/(costf-costi)*(k1f-k1i);
			    k[2] = k2i + (-costi)/(costf-costi)*(k2f-k2i);
			    k[3] = k3i + (-costi)/(costf-costi)*(k3f-k3i);					
//			printf("disk: pos = %e, %e, %e\n", pos[1], pos[2], pos[3]);	
			break;
		}		
				   
	}while(1);
return 0;	
}

int norm_vect(position, momentum)
double position[], momentum[];
{
double momentum_on[4], k_loc[4], kk;
sim5tetrad t; 
sim5metric metric;
                
// get metric and tetrad
kerr_metric(bh_spin, position[1], position[2], &metric);
                   
tetrad_zamo(&metric, &t);
bl2on(momentum, momentum_on, &t);
vect_copy(momentum_on, k_loc);
k_loc[0] = 1.0; k_loc[1] /= k_loc[0]; k_loc[2] /= k_loc[0]; k_loc[3] /= k_loc[0]; 
kk = sqrt(sqr(k_loc[1]) + sqr(k_loc[2]) + sqr(k_loc[3]));
k_loc[1] /= kk; k_loc[2] /= kk; k_loc[3] /= kk; 
on2bl(k_loc,momentum, &t); 

return 0;	
}	
