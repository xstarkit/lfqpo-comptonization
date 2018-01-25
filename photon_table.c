/*
Author: Bei You (youbeiyb@gmail.com)

This program is to tabulate the terminus of individual photon from the disk by use of raytracing code.
The grid of the table model is determined by position(r, phi) and momentum(cosi, cosf).
The grid of the table model contains the following information: position[4] and momentum[4];
If the photon goes to either infinity or the torus surface.
Therefore, for any given photon with (r0, phi0, cosi0, cosf0), we can derive the terminus by use of 4-D interpolation, no need for raytracing,
which can significantly reduce the computing time.
Sometime, photon from the disk travels just along the surface, and then goes to infinity; In this case (in the code, that means 'status=-1'), interpolation is nor applicable, should do
the raytracing.
   
*/

#include "mtc_incl_def.c"
#include "mtc_incl_code.c"
#include "quadrat.c"
#include "sim5lib.c"

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

long int Ntot = 0, Ntot1=0, Ntot2=0;

int main(int argc, char *argv[])
{
double    position[4], momentum[4];
int runnum;

  runnum = 0;
  if (argc == 2)   {
    runnum = atoi(argv[1]);
  }

/* initialize the geometry of:  1. the torus
                                2. the disk
                                3. the axis of the torus */
    input_data(runnum); // read input parameters and set up the geometry of the torus
      
	double r_min = Rdisk;     // inner radius of the disk
    double r_max = Rmax;         // outer radius of the disk    

/* convert the axis of the torus (0, 0, 1) in torus frame to Cartesian coordiante of BH, according to (prectheta, precphi), 
in order to check if photon is inside of the torus by means of Dot Product */
    if (precession_option<2) {
      torus_axis1(); // precession model (1)
    } else {
      torus_axis2(); // precession model (2)
    } 
//	printf("axis_vec[0]=%e, axis_vec[1]=%e, axis_vec[2]=%e\n", axis_vec[0], axis_vec[1], axis_vec[2]);
//	getchar();

//---------------------------------------------------------------------------------
double r, th0, ph0, momentum_on[4], dl;
double Omega;	
int status;
double U[4], cos_inc, phi_inf, phi, dphi, gd_P;

float inf_ph[nrd+1][nfd+1][nth+1][nph+1][8];

double dr = (log10(r_max)-log10(r_min))/nrd,
       df = (f_max-f_min)/nfd, 
       dt = (t_max-t_min)/nth, // cos(th0)
       dp = (p_max-p_min)/nph; // cos(ph0)
int ix, iy, iz, it, ii, iter;
double rnd_th, rnd_ph;

FILE *data; // opens new file for writing

char table_file[20];
char num[10];
strcpy(table_file,"table000000.dat");
sprintf(num,"%06i",runnum);
memcpy(table_file+5,num,6);		
data = fopen(table_file,"w");         


for (ix=0; ix<=nrd; ix++) {	
	for (iy=0; iy<=nfd; iy++) {
		for(iz=0; iz<=nth; iz++){
			for(it=0; it<=nph; it++){
				for(ii=0; ii<=7; ii++){
					inf_ph[ix][iy][iz][it][ii] = 0.0;
				}
			}
		}	
	}
}

for (ix=0; ix<=nrd; ix++) {	
	for (iy=0; iy<=nfd; iy++) {
		for(iz=0; iz<=nth; iz++){
			for(it=0; it<=nph; it++){
				iter++;
				if(iter%100000==0) printf("iter=%d\n", iter);
        sim5tetrad t;
        sim5metric m;
        geodesic gd;

        // draw an emission radius from distribution        
        r = pow(10.0,(log10(r_min)+ix*dr));
		phi = f_min + iy*df; 
		rnd_th = t_min + iz*dt;
		rnd_ph = p_min + it*dp;

     	// generate pos[4] in B-L spherical coordinate
	    position[0] = 0.0;
	    position[1] = r;
	    position[2] = 0.0;
	    position[3] = phi;
    
        // get metric and tetrad
        kerr_metric(bh_spin, r, 0.0, &m);   
		Omega = OmegaK(r,bh_spin);
        tetrad_azimuthal(&m, Omega, &t); 
        fourvelocity_azimuthal(Omega, &m, U);    
     
        // pick direction in disk rotating frame
        //double th0 = urand*M_PI/2.;
        th0 = rnd_th;
        ph0 = rnd_ph;
        momentum_on[0] = 1.0;
        momentum_on[1] = sqrt(1.0-rnd_th*rnd_th)*cos(ph0);
        momentum_on[2] = rnd_th;
        momentum_on[3] = sqrt(1.0-rnd_th*rnd_th)*sin(ph0);
//        momentum_on[1] = sin(th0)*cos(ph0);
//        momentum_on[2] = cos(th0);
//        momentum_on[3] = sin(th0)*sin(ph0);
        on2bl(momentum_on, momentum, &t); // transfer the direction from disk frame to B-L 

                // get geodesic
				geodesic_init_src(bh_spin, r, 0.0, momentum, momentum[1]<0?1:0, &gd, &status);
				if ((!status) || (r < gd.rp)) {
				//~ printf("geodesic_init_src failed (status=%d, r=%.2e, rp=%.2e)\n", status, r, gd.rp);
				inf_ph[ix][iy][iz][it][1] = -1.0; // radius being negative for failed raytracing
				for(ii=0; ii<=7; ii++){
					fprintf(data, "%f ", inf_ph[ix][iy][iz][it][ii]);
				}
					fprintf(data, "\n");				
				getchar();
				continue;
				}

				if(gd.m2p==1.0){
				//~ printf("gd.m2p=1.0\n");
				inf_ph[ix][iy][iz][it][1] = -1.0; // radius being negative for wrong raytracing
				for(ii=0; ii<=7; ii++){
					fprintf(data, "%f ", inf_ph[ix][iy][iz][it][ii]);
				}
					fprintf(data, "\n");				
				getchar();					
				continue;
				}

				// sometime a trajectory is found, which is correct, but unphysical
				// (photon would have pass radial turning point, which lies bellow horizon)
if ((creal(gd.r1) < r_bh(bh_spin)) && (momentum[1] < 0.0)) {
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 	                                                                  
        /* ==================== The block for ray-tracing ================================= */
        raytrace_data rtd;
	    raytrace_prepare(bh_spin, position, momentum, NULL, 0.1, 0, &rtd);
     		
 int iloop = 0, inside_torus = 0;   		
 	do{ iloop++;
	    dl = 1e9; // use maximal step     
	    raytrace(position , momentum, NULL, &dl, &rtd);
	    // after raytrace(), posi and momi already changed
	    
	    		// stop condition: if applying, stop raytrace
				if (position[1] <= r_bh(bh_spin)) {
				//~ printf("raytrace(): going to BH\n");
				inf_ph[ix][iy][iz][it][1] = -1.0;
				for(ii=0; ii<=7; ii++){
					fprintf(data, "%f ", inf_ph[ix][iy][iz][it][ii]);
				}
					fprintf(data, "\n");				
				getchar();
				break;
				}
								
				// also stop if relative error this step is too large
				if (rtd.error>1e-2) {
				//~ printf("raytrace(): aborted due to large error\n");;
				inf_ph[ix][iy][iz][it][1] = -1.0; // radius being negative for wrong raytracing due to large error
				for(ii=0; ii<=7; ii++){
					fprintf(data, "%f ", inf_ph[ix][iy][iz][it][ii]);
				}
					fprintf(data, "\n");				
				getchar();
				break;
				}
				
		if (position[1] > 1e6) {
			inf_ph[ix][iy][iz][it][1] = -1.0;			
			printf("if here, check the code\n");
			getchar();
				for(ii=0; ii<=7; ii++){
					fprintf(data, "%f ", inf_ph[ix][iy][iz][it][ii]);
				}
					fprintf(data, "\n");			
			break;
		}
				
		if(position[1]<=(1.1*Rout)){		
        inside_torus = inside_dotprod(position);
		}  
		if(inside_torus){
			Ntot2++;
			    inf_ph[ix][iy][iz][it][0] = position[0];
			    inf_ph[ix][iy][iz][it][1] = position[1];
			    inf_ph[ix][iy][iz][it][2] = position[2];
			    inf_ph[ix][iy][iz][it][3] = position[3];
			    inf_ph[ix][iy][iz][it][4] = momentum[0];	
			    inf_ph[ix][iy][iz][it][5] = momentum[1];
			    inf_ph[ix][iy][iz][it][6] = momentum[2];
			    inf_ph[ix][iy][iz][it][7] = momentum[3];			    			    		    			    			    
				for(ii=0; ii<=7; ii++){
					fprintf(data, "%f ", inf_ph[ix][iy][iz][it][ii]);
				}
					fprintf(data, "\n");
		}	
		if(iloop>10000) {
			printf("too much step \n");
			getchar();			
			break;
		}    
	}while(inside_torus!=1);			
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
continue; // go to next grid
}

				gd_P = geodesic_P_int(&gd, r, momentum[1]<0?1:0);
				dphi = geodesic_position_azm(&gd, r, 0.0, gd_P);
				// only save theta angle and azimuthal angle for photons at infinity
				inf_ph[ix][iy][iz][it][1] = 1e6;				
				inf_ph[ix][iy][iz][it][2] = gd.cos_i;
				inf_ph[ix][iy][iz][it][3] = reduce_angle_2pi(phi + dphi);				

				// if the radial turning point is outside torus extent (meanwhile photon moves inwards BH), 
				// photon will never enter into torus; Instead, going to infinity				
				if(gd.rp>Rout){// write out to the file 	
					if(gd.cos_i>0.0){ 												
						Ntot1++;
						inf_ph[ix][iy][iz][it][1] = 1e6;					
					}else{ // ignore the photon which travel below the disk and finally goes to infinity
						inf_ph[ix][iy][iz][it][1] = -1e6;							
					}
					for(ii=0; ii<=7; ii++){
						fprintf(data, "%f ", inf_ph[ix][iy][iz][it][ii]);
					}
						fprintf(data, "\n");						
					continue; // !!! go to next grid
				}			
							    			     
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 	
// !!!!!!!!!!!!! only photons moving inwards are possible to enter into torus                                                                    
        /* ==================== The block for ray-tracing ================================= */
        raytrace_data rtd;
	    raytrace_prepare(bh_spin, position, momentum, NULL, 0.1, 0, &rtd);
     		
 int iloop = 0, inside_torus = 0;   		
 	do{ iloop++;
	    dl = 1e9; // use maximal step     
	    raytrace(position , momentum, NULL, &dl, &rtd);
	    // after raytrace(), posi and momi already changed
	    
	    		// stop condition: if applying, stop raytrace
				if (position[1] <= r_bh(bh_spin)) {
				printf("raytrace(): going to BH; But should not\n");
				inf_ph[ix][iy][iz][it][1] = -1.0;
				for(ii=0; ii<=7; ii++){
					fprintf(data, "%f ", inf_ph[ix][iy][iz][it][ii]);
				}
					fprintf(data, "\n");				
				getchar();
				break;  
				}
			
				// also stop if relative error this step is too large
				if (rtd.error>1e-2) {
				printf("raytrace(): aborted due to large error; But should not\n");;
				inf_ph[ix][iy][iz][it][1] = -1.0; // radius being negative for wrong raytracing due to large error
				for(ii=0; ii<=7; ii++){
					fprintf(data, "%f ", inf_ph[ix][iy][iz][it][ii]);
				}
					fprintf(data, "\n");				
				getchar();
				break;
				}
				
		if((position[1]>1.1*Rout)&&(momentum[1]>0.0)){ // if this applies, it means that photon do not hit torus, and go to infinity
		// !!! very important condition, so that to reduce the spending time to find the entering-torus photon 	
		
				if(gd.cos_i>0.0){ 												
					Ntot1++;
					inf_ph[ix][iy][iz][it][1] = 1e6;					
				}else{ // ignore the photon which travel below the disk and finally goes to infinity
					inf_ph[ix][iy][iz][it][1] = -1e6;							
				}
				for(ii=0; ii<=7; ii++){
					fprintf(data, "%f ", inf_ph[ix][iy][iz][it][ii]);
				}
					fprintf(data, "\n");						
		
				break;
		}
				
		if(position[1]<=(1.1*Rout)){		
			inside_torus = inside_dotprod(position);
		}  
		if(inside_torus){
			Ntot2++;
			    inf_ph[ix][iy][iz][it][0] = position[0];
			    inf_ph[ix][iy][iz][it][1] = position[1];
			    inf_ph[ix][iy][iz][it][2] = position[2];
			    inf_ph[ix][iy][iz][it][3] = position[3];
			    inf_ph[ix][iy][iz][it][4] = momentum[0];	
			    inf_ph[ix][iy][iz][it][5] = momentum[1];
			    inf_ph[ix][iy][iz][it][6] = momentum[2];
			    inf_ph[ix][iy][iz][it][7] = momentum[3];			    			    		    			    			    
				for(ii=0; ii<=7; ii++){
					fprintf(data, "%f ", inf_ph[ix][iy][iz][it][ii]);
				}
					fprintf(data, "\n");
		}	
		if(iloop>10000) {
			printf("too much step \n");
			getchar();			
			break;
		}    
	}while(inside_torus!=1);			
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			}
		}	
	}
}              
printf("N = %ld, %ld, %ld\n", Ntot1, Ntot2, Ntot1/Ntot2);                        			 								 
	double speed;						
	speed = (nrd+1)*(nfd+1)*(nth+1)*(nph+1)/(time(NULL)-start_time);
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
		out_comp[indx][iindx][findx] += 1.0;		 
	  }	
	}                
/* -------------------- The block for saving photon energy and position ---------------------------------- */
return 0;	
}


output_spec(runnum,noph)
  long int   noph;
  int        runnum;
{
  FILE   *fp1, *fp2, *fp3;
  int     i,nos,j,k; 
  double  e,st,si,sr,er,d,de;

  char   out_file[20],
         inp_spec[20],
         reflspec[20];
  char    num[10];

  if (noph == 0) return 0;

  strcpy(out_file,"mcomp000000.dat");
  strcpy(inp_spec,"inpsp000000.dat");
  strcpy(reflspec,"refls000000.dat");


  sprintf(num,"%06i",runnum);
  memcpy(out_file+5,num,6);
  memcpy(inp_spec+5,num,6);
  memcpy(reflspec+5,num,6);
		
fp1 = fopen(inp_spec,"w");         /* output observed Novikov-Thorne disk spectrum */
fp2 = fopen(reflspec,"w");         /* output observed reflection spectrum from the disk */
fp3 = fopen(out_file,"w");         /* output observed torus Comptonization spectrum */

double er_disk, er_torus; 
  nos = num_bin_o*log10(emax_o/emin_o);
  for (i=0; i<nos; i++){
    e = emin_o*pow(10.,(i+0.5)/num_bin_o );
    de = emin_o*(pow(10.,(i+1.0)/num_bin_o)-pow(10.,(1.0*i)/num_bin_o));
    er = e/de/noph;
    for(j=0; j<nangle; j++) {
		for(k=0; k<MAXphi; k++) {
			fprintf(fp1,"%e %.3e ",e*511., out_disk[i][j][k]*er*nangle);
			fprintf(fp2,"%e %.3e ",e*511., out_refl[i][j][k]*er*nangle);
			fprintf(fp3,"%e %.3e ",e*511., out_comp[i][j][k]*er*nangle);
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
        ph0 = urand*M_PI*2.;
        momentum_on[0] = 1.0;
        momentum_on[1] = sin(th0)*cos(ph0);
        momentum_on[2] = cos(th0);
        momentum_on[3] = sin(th0)*sin(ph0);
        on2bl(momentum_on, momentum, &t); // transfer the direction from disk frame to B-L 

	if(fabs(momentum[1])<1e-7) continue; 
        
        gf = (momentum[0]*m.g00 + momentum[3]*m.g03) / dotprod(momentum, U, &m);
		*E_inf = energy*gf;   // photon energy measured by observer at infinity, namely E_inf  
                // get geodesic
				geodesic_init_src(bh_spin, r, 0.0, momentum, momentum[1]<0?1:0, &gd, &status);
				if ((!status) || (r < gd.rp)) {
				fprintf(stderr, "WRN: geodesic_init_src failed (status=%d, r=%.2e, rp=%.2e)\n", status, r, gd.rp);
				continue;
				}

				if(gd.m2p==1.0)	continue;

				// sometime a trajectory is found, which is correct, but unphysical
				// (photon would have pass radial turning point, which lies bellow horizon)
				if ((creal(gd.r1) < r_bh(bh_spin)) && (momentum[1] < 0.0)) {
				// unphysical solution - geodesic goes to observer through BH
				continue;
				}
				
				gd_P = geodesic_P_int(&gd, r, momentum[1]<0?1:0);
				dphi = geodesic_position_azm(&gd, r, 0.0, gd_P);
									
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
				if((momentum[1]<0.0)&&(gd.rp>=Rout)){
					if((cos_inc>0.0) && (findx>=0)){ // only deal with photon energy in the interesting range	  
								out_disk[indx][iindx][findx] += 1;
					}	
					Ntot++; // !!! accumulate photon number		
				continue; // !!! must go to the new loop
				} 			
				if(momentum[1]>=0.0){ 
					if((cos_inc>0.0) && (findx>=0)){ // only deal with photon energy in the interesting range					  
								out_disk[indx][iindx][findx] += 1;	
					}	
					Ntot++; // !!! accumulate photon number			
				continue;					    			    
				}				
//	break;				
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<	
// !!!!!!!!!!!!! only photons moving inwards are possible to enter into torus 
                                                                     
        /* ==================== The block for ray-tracing ================================= */
        raytrace_data rtd;
	    raytrace_prepare(bh_spin, position, momentum, NULL, 1.0, 0, &rtd);
     		
 int iloop = 0;   		
 	do{ iloop++;
	    dl = 1e6; // use maximal step     
	    raytrace(position , momentum, NULL, &dl, &rtd);
	    // after raytrace(), posi and momi already changed
	    
	    		// stop condition: if applied, stop raytrace
				if ((position[1] < r_bh(bh_spin)) || (position[1] > 1e6)) {
//				printf("raytrace(): photon goes to bh event horizon\n");
				break;  
				}
				// also stop if relative error this step is too large
				if (rtd.error>1e-2) {
//				printf("raytrace(): aborted due to large error\n");
				break;
				}
				
		if((position[1]>Rout)&&(momentum[1]>0.0)){
		// !!! very important condition, so that to reduce the spending time to find the entering-torus photon 
					if((cos_inc>0.0) && (findx>=0)){ // only deal with photon energy in the interesting range		
								out_disk[indx][iindx][findx] += 1;		
					}	
			Ntot++; // !!! accumulate photon number					
			break;
		}
				
		if(position[1]<=(1.1*Rout)){		
        inside_torus = inside_dotprod(position);
		}  
		if(iloop>10000) {printf("too much step \n");break;}    
	}while(inside_torus!=1);
//	getchar();			
	
}while((inside_torus!=1) && (iter < 10000));     	

		Ntot++; // !!! accumulate photon number	

if(iter>=10000){
	printf("iter: %d\n", iter);
//	getchar();
}    	
                 				                                                        
return 0;
}


int raytrace_out_torus(position, momentum, i_path) 
/* a function of raytracing photon which is outside torus */
double position[], momentum[];
int *i_path; 
{
		double dl, costi, costf;
		                                                                                                  
        /* ==================== The block for ray-tracing ================================= */
	    raytrace_data rtd;
	    raytrace_prepare(bh_spin, position, momentum, NULL, 10.0, 0, &rtd);
		
	*i_path = 3; // default value 
	do{
	    dl = 1e4; // use maximal step   
	    costi = position[2];	  
	    raytrace(position , momentum, NULL, &dl, &rtd);
	    costf = position[2]; 
	    // after raytrace(), posi and momi already changed
	    
	    		// stop condition: if applied, stop raytrace
				if ((position[1] < r_bh(bh_spin))) {
				printf("out torus1, raytrace(): photon goes to bh event horizon\n");
				*i_path = 3;  
				return 0;  
				}
				// also stop if relative error this step is too large
				if (rtd.error>1e-2) {
				printf("out torus2, raytrace(): aborted due to large error\n");
				*i_path = 3; // failure step; as photon goes to BH event horizon, 
				return 0;;
				}
       
		if(position[1]>Rdisk){
			if((costi>0.0)&&(costf<0.0)){ // chech if photon go across the disk; 		     							  
			    *i_path = 1;
			    return 0;
			}else if((costi<0.0)&&(costf>0.0)){ // hit the disk from the bottom; then ignore it
				*i_path = 3;
				return 0;
			}
		}  
//	printf("pos = %e, %e\n", acos(position[2])/3.1415*180, position[3]);	
	}while(position[1]<1e4);	
//	getchar();
	*i_path = 0;		
		/* ==================== End of the block for ray-tracing ================================= */  		
		
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

		do{ // loop to save each step along geodesic						
			dl = 0.1/tTh1Rg; // have to use maximal step; rather than 0		
			inside_torus = inside_dotprod(posi_orig);	

		  if(inside_torus){					
									
			vect_copy(posi_orig, posmin);
			vect_copy(momi_orig, mommin);
			pmin = p_sc;									
			raytrace(posi_orig , momi_orig, NULL, &dl, &rtd); 
			vect_copy(posi_orig, posmax);
			vect_copy(momi_orig, mommax);

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
			vect_copy(posi_orig, position);
			vect_copy(momi_orig, momentum);				
			return 1;
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
	
	raytrace_out_torus(position_esc, momentum_esc, &i_path); /* determine whether escaping photon goes to infinity, or reflection */

	if(i_path==3) return -1;	// failed raytrace, going on to next scattering of this photon	

	if(i_path==0){// from torus to infinity
			write_spec(position_esc, momentum_esc, E_inf);	
	} 			                                                                                                

	if(i_path==1){ // this is the reflection
			prod_pU(position_esc, momentum_esc, &gf); // After doing raytrace, calculating "g_f" in order to determine the photon energy at the next position
			E_loc = E_inf/gf;
			disk_reflection(position_esc, momentum_esc, E_loc, weight_esc);
	}	
 
return 0;
}

int disk_reflection(position, momentum, ee, ww)
double position[], momentum[], ee, ww;
{
  double mom_i[3],mom_f[3],pos_l[3],
         energy,weight,weight_min,mu_prime,sine2,term,bound;	 
  int i;
  double mom_r[4];
  
  energy = ee;
  weight = ww;
    
  // both postion and direction have to be in cartesian coordinate !!! 
  bl2torus(position, momentum, mom_r); // convert the B-L direction to the disk-rest frame 
//  printf("mom[] = %e %e %e %e\n", mom_r[0], mom_r[1], mom_r[2], mom_r[3]);
/* NOTE: rerange the order of the x,y,z direction in order to be consistent with mom_i[] which is differently defined by Piotr */
  mom_i[0] = mom_r[1]/mom_r[0]; 
  mom_i[1] = mom_r[3]/mom_r[0];
  mom_i[2] = mom_r[2]/mom_r[0];
//  for (i=0; i<3; i++) mom_i[i] = mom_r[i+1]/mom_r[0]; // normalization, make sure that kx*kx + ky*ky + kz*kz = 1

  pos_l[0] = position[1]*fabs(sqrt(1.0-sqr(position[2])))*cos(position[3]);
  pos_l[1] = position[1]*fabs(sqrt(1.0-sqr(position[2])))*sin(position[3]);
  pos_l[2] = position[1]*position[2];				
//  printf("x = %e, y = %e, z = %e\n", pos_l[0], pos_l[1], pos_l[2]);
//  printf("kx = %e, ky = %e, kz = %e\n", mom_i[0], mom_i[1], mom_i[2]);
//  printf("k.k = %e\n", sqr(mom_i[0]) + sqr(mom_i[1]) + sqr(mom_i[2]));
//  getchar();
	
  double r0 = r_bh(bh_spin);
  weight_min = weight*1e-5;
//  if (weight_min < wghtmin) weight_min = wghtmin;
  
  int i_esc;
  i_esc = raytrace_in_disk(position, pos_l, mom_i, energy, &weight);	
//  if(i_esc==1) return 0;

  do {
    do {
      mu_prime = 1-(exp(log(1+2*energy)*rnd())-1)/energy;
      sine2 = 1-mu_prime*mu_prime;
      term  = 1/(1+energy*(1-mu_prime));
      bound = term*(term+1/term-sine2)/2;
    } while (rnd() > bound);
    transform(mom_i,mom_f,sqrt(sine2),mu_prime);
    for (i=0; i<3; i++)
      mom_i[i] = mom_f[i];
    
    energy *= term;
    i_esc = raytrace_in_disk(position, pos_l, mom_i, energy, &weight);
//    if(i_esc==1) return 0;
//    if((pos_l[0]<=r0) || (pos_l[1]<=r0) || (pos_l[2]<=r0))	return 0;
  }
  while ( (weight > weight_min) && (energy > emin_o) );  	
return 0;		
}

int raytrace_in_disk(posi, pos, mom, energy, weight)
// pos[], mom[] are 3-vector in cartesian coordinate, in disk-rest frame
double posi[], pos[], mom[], energy, *weight;
{
  double sigma_es,sigma_abs,sigma_tot,p_esc,wesc,lambda,d,fi;
  int    indx,iindx,i,findx;

  sigma_es = sigma_KN(energy);
  sigma_abs= sigma_bf(energy);
  sigma_tot= sigma_es+sigma_abs;
          
  if (mom[2] <= 0) { /* downwards - no escape */
    p_esc = 0;
  } else {
    d = fabs(pos[2]/mom[2]);
    p_esc = exp(-d*sigma_tot);
    wesc = p_esc*(*weight);
			raytrace_infinity(posi, mom, energy, wesc); // assume that the photon enters and is reflected 
														// at the same position on the disk: posi[], rather than pos[]
  }
  lambda = -log(1-(1-p_esc)*rnd())/sigma_tot;
  
  for (i=0; i<3; i++)
    pos[i] += lambda*mom[i];                 /* position of the next    
						scattering               */  
  (*weight) *= (1-p_esc)*sigma_es/sigma_tot;
          
return 0;
}

int raytrace_infinity(pos, dir, energy, weight)
// pos[] are 4-vector in spherical coordinate, dir[] are 3-vector in cartesian coordinate, in disk-rest frame
double pos[], dir[], energy, weight;
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
        return -1;
    }

    double gd_P = geodesic_P_int(&gd, r, k[1]<0?1:0);

    // sometime a trajectory is found, which is correct, but unphysical
    // (photon would have pass radial turning point, which lies bellow horizon)
    if ((creal(gd.r1) < r_bh(bh_spin)) && (gd_P > gd.Rpa)) {
        // unphysical solution - geodesic goes to observer through BH
        return -1;
    }

    double dphi = geodesic_position_azm(&gd, r, m, gd_P);

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
	    
	if(cos_inc>0.0){ // only deal with photon energy in the interesting range 	
		indx = num_bin_o*log10(energy/emin_o); // indx corresponding to photon energy at infinity 
			    
	  if(findx>=0){
	    iindx = cos_inc*nangle;    
	    out_refl[indx][iindx][findx] += weight;		
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
