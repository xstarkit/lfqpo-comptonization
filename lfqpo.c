/*
Author: Bei You 

  Compton scattering in a uniform, hot  plasma cloud.
  Based on Gorecki & Wilczewski and Pozdnyakov, Sobol & Sunyaev

  This code is the combination of Piotr's (Comptonization) and Michal's (GR effect) codes.
  
  Currently, it only concentrates on the wedge geometry for Lense-Thirring precession project (LFQPO). 
  It will be extented to other geometries after this project.

*/

#include "mtc_incl_def.c"
#include "mtc_incl_code.c"
#include "quadrat.c"

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

int main(int argc, char *argv[])
{
    long int il, istep;

    double    position[4], momentum[4], momentum_torus[4], weight;
    double    position_esc[4], momentum_esc[4], weight_esc;
    double 	  gi, gf; // dotprod(p,U)
    int       runnum;
    int 	  i, indx, iindx, findx;
	double 	  theta_inf, phi_inf;
	double 	  E_loc, E_inf;

  runnum = 0;
  if (argc == 2)   {
    runnum = atoi(argv[1]);
  }


/* initialize the geometry of:  1. the torus
                                2. the disk
                                3. the axis of the torus */
    input_data(runnum); // read input parameters and set up the geometry of the torus
    
	double r_min = Rdisk;          // inner radius of the disk
    double r_max = 5000.0;         // outer radius of the disk 
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
        
int i_path; /* photon trajectory by raytrace:
				0: to infinity; 1: to disk (for reflection); 2: to torus; 3: to BH */

int not_scattering, if_failure;

for (il=0; il<=mc_steps; il++) { // generate photons
	
    if ((il % 200) == 0) { 
      /*        output_results(runnum,il); */
      writelog(runnum,il); 
    }
	
	generate_photon_gr(position, momentum, &E_loc, &gf);
			
	weight = 1;
	wghtmin = min_weight*weight;                                                  	
		 
			do{ // loop for multi-scattering in the torus
				if_failure = 0;
				/* determine the scattering position (position, momentum) and escaping position (position_esc, momentum_esc) */
				not_scattering = raytrace_in_torus(position, momentum, E_loc, &weight, position_esc, momentum_esc, &weight_esc); 
				if(not_scattering==-1) break; // failed raytrace, skipping this photon; Normally raytrace inside torus should always succeed
				prod_pU(position_esc, momentum_esc, &gi); 
				E_inf = E_loc*(gi/gf); /* calculate the photon energy in the position where it escapes from the torus */
									   /* Notice that, E_inf has to be calculated before E_loc is changed */
				prod_pU(position, momentum, &gi); 
				E_loc = E_loc*(gi/gf); /* calculate the photon energy in the position where it will be scattered next time */					
// ----------------------------------------------------------------------------------------------------------------------------------------			
				// this is the not-scattering part; photon will go to infinity, or disk (then reflection)
				if_failure = raytrace_esc_torus(position_esc, momentum_esc, E_inf, weight_esc); 												     											   
				if(if_failure==-1) break; // failed raytrace, skipping this photon
				
				// this is the scattering part;
				if(not_scattering==1) break; // only escaping part is done
				scattering_gr(position, momentum, &E_loc); // this is the scattering; global variable "E_loc" is already changed. 
// ----------------------------------------------------------------------------------------------------------------------------------------						
				prod_pU(position, momentum, &gf); // Since photon momentum is changed after scattering, new "gf" is required for calculating photon energy 	
			}while((weight > wghtmin)&&(E_loc > emin_o)&&(E_loc < emax_o));
								    
} // end of all photons 
	output_spec(runnum,il); /* output spectrum: 
					                   Disk reflection spectrum;
					                   Torus Comptonization spectrum;
							*/
	double speed;						
	speed = (1.0*mc_steps)/(time(NULL)-start_time);
    printf(" Speed of computations: %.1f photons/sec \n",speed);

    distrib_done(&dist);		
    					
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


int write_spec(position, momentum, energy, i_spec, weight)
double position[], momentum[], energy, weight;
int i_spec;
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
	    
	    switch (i_spec){
	    case 0: {out_disk[indx][iindx][findx] += weight;
				 break;}
		case 1: {out_refl[indx][iindx][findx] += weight;
				 break;}
		case 2: {out_comp[indx][iindx][findx] += weight;			
				 break;}	 
		}	
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
		
//fp1 = fopen(inp_spec,"w");         /* output observed Novikov-Thorne disk spectrum */
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
//			fprintf(fp1,"%e %.3e ",e*511., out_disk[i][j][k]*er*nangle);
			fprintf(fp2,"%e %.3e ",e*511., out_refl[i][j][k]*er*nangle);
			fprintf(fp3,"%e %.3e ",e*511., out_comp[i][j][k]*er*nangle);
		}
//		fprintf(fp1,"  ");
		fprintf(fp2,"  ");
		fprintf(fp3,"  ");
    }
//    fprintf(fp1,"\n");
    fprintf(fp2,"\n");
    fprintf(fp3,"\n");
  }
//  fclose(fp1);
  fclose(fp2);
  fclose(fp3);

return 0;
}


int generate_photon_gr(position, momentum, energy, gf)
double position[], momentum[], *energy, *gf; 
{

double r, T;
double g, th0, ph0, momentum_on[4];
double dl;
int inside_torus;
double Omega;
raytrace_data rtd;
sim5tetrad t;
sim5metric m;	
geodesic gd;
int status;

do{		
        // draw an emission radius from distribution        
        do {
        r = exp(distrib_hit(&dist));
        T = disk_nt_temp(r); 
//		printf("r = %.2e, T = %.2e\n", r, 2.7*T*8.617328149741e-8);
		} while ((2.7*T*8.617328149741e-8/511.) < emin); // convert temperature to keV (normalized to 511 keV), to set up minimum temperature (maximum disk outer radius)

/* Generate energy from a planckian distribution of a corresponding temperature */
		  do{
		  *energy = blackbody_photon_energy_random(T)/511.; /* photon energy in disk-rest frame */
		  }while ( ((*energy) < emin) || ((*energy) > emax));		  

     	// generate pos[4] in B-L spherical coordinate
	    position[0] = 0.0;
	    position[1] = r;
	    position[2] = 0.0;
	    position[3] = urand*M_PI*2.0;
    
        // get metric and tetrad
        kerr_metric(bh_spin, r, 0.0, &m);   
		Omega = OmegaK(r,bh_spin);
        tetrad_azimuthal(&m, Omega, &t);    
     
	  do{
        // pick direction in disk rotating frame
        //double th0 = urand*M_PI/2.;
        th0 = acos(sqrt(1.-sim5urand())); // iCDF that corresponds to PDF(theta)=sin(theta)*cos(theta)
        ph0 = urand*M_PI*2.;
        momentum_on[0] = 1.0;
        momentum_on[1] = sin(th0)*cos(ph0);
        momentum_on[2] = cos(th0);
        momentum_on[3] = sin(th0)*sin(ph0);
        on2bl(momentum_on, momentum, &t); // transfer the direction from disk frame to B-L 
                // get geodesic
				geodesic_init_src(bh_spin, r, 0.0, momentum, momentum[1]<0?1:0, &gd, &status);
				if ((!status) || (r < gd.rp)) {
				fprintf(stderr, "WRN: geodesic_init_src failed (status=%d, r=%.2e, rp=%.2e)\n", status, r, gd.rp);
				continue;
				}

				// sometime a trajectory is found, which is correct, but unphysical
				// (photon would have pass radial turning point, which lies bellow horizon)
				if ((creal(gd.r1) < r_bh(bh_spin)) && (momentum[1] < 0.0)) {
				// unphysical solution - geodesic goes to observer through BH
				continue;
				}
				// if the radial turning point is outside torus extent (meanwhile photon moves inwards BH), 
				// photon will never enter into torus; Instead, going to infinity
				if((momentum[1]<0.0)&&(gd.rp>Rout)) continue; 
				
       }while(momentum[1]>0.0); // only photons moving inwards are possible to enter into torus 
                                                                                
        /* ==================== The block for ray-tracing ================================= */
	    raytrace_prepare(bh_spin, position, momentum, NULL, 1.0, 0, &rtd);
		
	inside_torus = 0;
	do{
	    dl = 1e9; // use maximal step       
	    raytrace(position , momentum, NULL, &dl, &rtd);
	    // after raytrace(), posi and momi already changed
	    
	    		// stop condition: if applied, stop raytrace
				if ((position[1] < r_bh(bh_spin))) {
//				printf("raytrace(): photon goes to bh event horizon\n");
				break;  
				}
				// also stop if relative error this step is too large
				if (rtd.error>1e-2) {
//				printf("raytrace(): aborted due to large error\n");
				break;
				}
				
		if((position[1]>Rout)&&(momentum[1]>0.0)) break; 
		// !!! very important condition, so that to reduce the spending time to find the entering-torus photon //
				
		if(position[1]<=(Rout+0.1)){		
        inside_torus = inside_dotprod(position);
		}        
		
	}while(inside_torus!=1);			
	
}while(inside_torus!=1);      
			
		prod_pU(position, momentum, gf); // (*gf) is negative
		*energy = -(*energy)*(*gf); // photon energy measured in torus frame, 
									// which is redshifted with relative to the initial energy emitted from the disk
//		printf("E_ini = %e\n", *energy);
                                                       
return 0;
}


int raytrace_out_torus(position, momentum, i_path, type) 
/* a function of raytracing photon which is outside torus */
double position[], momentum[];
int *i_path, type; /* if type = 3, torus, disk is taken into account; 
					  if type = 2, only disk is taken into account;
					  if type = 1, photon travel in vacuum;
					*/
{
		double dl, costi, costf;
		int inside_torus;
		                                                                                                  
        /* ==================== The block for ray-tracing ================================= */
	    raytrace_data rtd;
	    raytrace_prepare(bh_spin, position, momentum, NULL, 1.0, 0, &rtd);
		
int with_torus = 1, with_disk = 1;	    
switch(type){
	case 3:{
		with_torus = 1;
		with_disk = 1;
		break;
	}   
	case 2:{
		with_torus = 0;
		with_disk = 1;
		break;
	}
	case 1:{
		with_torus = 0;
		with_disk = 0;
		break;
	}
}
	*i_path = 3; // default value 
	do{
	    dl = 1e9; // use maximal step   
	    costi = position[2];	  
	    raytrace(position , momentum, NULL, &dl, &rtd);
	    costf = position[2]; 
	    // after raytrace(), posi and momi already changed
	    
	    		// stop condition: if applied, stop raytrace
				if ((position[1] < r_bh(bh_spin))) {
				printf("out torus1, raytrace(): photon goes to bh event horizon\n");
				*i_path = 3;  
				break;  
				}
				// also stop if relative error this step is too large
				if (rtd.error>1e-2) {
				printf("out torus2, raytrace(): aborted due to large error\n");
				*i_path = 3; // failure step; as photon goes to BH event horizon, 
				break;
				}
        inside_torus = inside_dotprod(position);        
        if(position[1]>1e9) { // Is this value large enough to stop raytrace ? 
			*i_path = 0;
			return 0;
		} 
		if((with_disk==1)&&(position[1]>Rdisk)&&(costi>0.0)&&(costf<0.0)){ // chech if photon go across the disk; 
															  // if photon go across the disk from the bottom, keep the default of i_path = 3 
			*i_path = 1;
			return 0;
		} 
		if((with_torus==1)&&(inside_torus)){
			*i_path = 2;
			return 0;
		}			
	}while(1);			
		/* ==================== End of the block for ray-tracing ================================= */  		
		/* weight is not changed during ray-tracing */
		
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


int prod_pU(position, momentum, pU) /* dot product between photon momentum and medium 4-velocity */
double position[], momentum[], *pU;
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
//        printf("Omega = %e, Omegak = %e\n", Omega, Omegak);
		*pU = dotprod(momentum, U, &m);
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

int raytrace_in_torus(position,momentum,energy,weight,position_esc,momentum_esc,weight_esc)  
/* assuming photon is inside torus, given position[] and momentum[], to determine:
1. accumulated escaping probalibity along geodesic to boundary
2. if return -1, raytrace failure	
   if return 1, no scattering, only have the escaping photon
   if return 0, both scattering and escaping													
*/																                                                     
double position[],momentum[],*weight;
double position_esc[],momentum_esc[],*weight_esc;
double energy;
{
double  posi_orig[4], momi_orig[4], momi_torus[4];
double  dl, p_sc, p_esc, dtau, tau;
int inside_torus;
double energy1, energy2, energy_m;
double gf_i, gf_f; // gravitational redshift factor
                
        /* copy vectors to another temporary ones which are used in raytrace */
        vect_copy(position, posi_orig); // posi -> posi_orig
        vect_copy(momentum, momi_orig); // momi -> momi_orig
                         
        list_init(); // initialize list_node
        
        /* ==================== The block for ray-tracing ================================= */       
	    raytrace_data rtd;
	    raytrace_prepare(bh_spin, posi_orig, momi_orig, NULL, 1.0, 0, &rtd);

		tau  = 0.0;
		dtau = 0.0;
		p_sc = 0.0; // accumulated scattering probalibity: p_sc = 1 - exp(-tau), where tau += sigma*dl 
		p_esc = 0.0;
		
		energy2 = energy; /* 1. energy in torus-rest frame; 
							 2. E_loc should not be changed after this subroutine */
		energy1 = energy2;

		do{ // loop to save each step along geodesic						
			dl = 1e9; // have to use maximal step; rather than 0		
			inside_torus = inside_dotprod(posi_orig);	

			if(inside_torus){					
			tau += dtau; // accumulate optical depth before escaping torus
			p_sc = 1.0 - exp(-tau);
			p_esc = exp(-tau);
			
            list_node_add(posi_orig , momi_orig, p_sc, energy2);  /* starting from the first step, save data into link-node; 
														   If, after the first step, photon is outside of torus, meaning no scattering */															  											   
									   						   
			}else{	
					// if no scattering, get weight after escaping
					/* copy temorary vectors to final ones as output */
					vect_copy(posi_orig, position_esc); // posi_orig -> position
					vect_copy(momi_orig, momentum_esc); // momi_orig -> momentum
					(*weight_esc) = p_esc*(*weight);
					if(p_sc==0.0){
						return 1; /* if here, that means that after the first step photon is just out of the torus, 
									then save the position and momention as escape one, and no scattering */
					}else{
						break; // once outside of torus, jump out of raytrace loop	
					}
			}	
			prod_pU(posi_orig, momi_orig, &gf_i);	
								
			raytrace(posi_orig , momi_orig, NULL, &dl, &rtd); 
			
			prod_pU(posi_orig, momi_orig, &gf_f);
			energy2 = energy2*(gf_f/gf_i); // photon energy (in torus-rest frame) changed after each raytrace step
			
			energy_m = 0.5*(energy1 + energy2);
			dtau = ray_tau(dl, energy_m); 
//			printf("dtau = %e\n", dtau);
			energy1 = energy2;	// energy1 as the initial energy for the next step	  
	 
			// stop condition: if applied, stop raytrace
			if ((posi_orig[1] < r_bh(bh_spin))) {
			printf("in torus, raytrace(): photon goes to bh event horizon\n");
			return -1;    
			}
	
			// also stop if relative error this step is too large
			if (rtd.error>1e-2) {
	        printf("in torus, raytrace(): aborted due to large error\n");
	        return -1;
			}
		}while(1);	
		/* ==================== The block for ray-tracing ================================= */
/* =============== to determine the scattering position, if photon is scattered =================*/
    double posmax[4],posmin[4],
		   mommax[4],mommin[4],
		   pmin, pmax, elmax, elmin;
	int i;
	double x;
			
					x = urand; 
					if(x==0.0) x+=0.001;
					if(x==1.0) x-=0.001; // ensure that x = (0,1), not [0,1]
					x = x*p_sc; // critical random value (0,1) to determine the position where scattering happens 

					LIST_NODE *pb = pHead;		
					while(pb->pNext) // loop to find the index so that p(i) < x < p(i+1)
					{
						/* 如果record值大于value，返回深度index */
						if(pb->pp > x){
						for (i=0; i<4; i++) posmax[i] = pb->pos[i];
						for (i=0; i<4; i++) mommax[i] = pb->mom[i];
						pmax = pb->pp;
						elmax = pb->ee;						
						break;
						}else{
						for (i=0; i<4; i++) posmin[i] = pb->pos[i];
						for (i=0; i<4; i++) mommin[i] = pb->mom[i];
						pmin = pb->pp;
						elmin = pb->ee;
						}	
					/* 流程走到这里，说明需要判断下一个节点 */
					pb = pb->pNext;
					/* 深度index更新 */
					}
					list_free();
					/* do interpolation to determine new pos[], mom[], energy and weight, corresponding to x */
//					printf("interpolation:\n");
					for (i=0; i<4; i++){
					position[i] = posmin[i] + (x-pmin)*(posmax[i]-posmin[i])/
												   (pmax-pmin); 
					momentum[i] = mommin[i] + (x-pmin)*(mommax[i]-mommin[i])/
												   (pmax-pmin); 
//					printf("position = %e, momentum = %e\n", position[i], momentum[i]);							    
					}
					(*weight) *= p_sc; // new weight 
//					printf("weight = %e\n", *weight);			
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

int raytrace_esc_torus(double position_esc[], double momentum_esc[], double E_inf, double weight_esc)
{
	int i_path;
	int i_spec; 	/* 0: accumulate photon weight as direct disk spectrum 
					   1:                             reflection spectrum
					   2:                             compton spectrum */
	double gf_i, gf_f; // gravitational redshift factor
	prod_pU(position_esc, momentum_esc, &gf_i); // Once photon escapes from the torus, calculating "g_o" 

	raytrace_out_torus(position_esc, momentum_esc, &i_path, 2); /* determine whether escaping photon goes to infinity, or reflection */
	if(i_path==3) return -1;	// failed raytrace, going on to next scattering of this photon	
		
	prod_pU(position_esc, momentum_esc, &gf_f); // After doing raytrace, calculating "g_f" in order to determine the photon energy at the next position
	E_inf = E_inf*(gf_f/gf_i);
	if(i_path==0){
			i_spec = 2; // from torus to infinity
			write_spec(position_esc, momentum_esc, E_inf, i_spec, weight_esc);	
	}			                                                                                                
	if(i_path==1){ // this is the reflection
			disk_reflection(position_esc, momentum_esc, E_inf, weight_esc);
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
  if (weight_min < wghtmin) weight_min = wghtmin;
  
  raytrace_in_disk(pos_l, mom_i, energy, &weight);	
  
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
    
    raytrace_in_disk(pos_l, mom_i, energy, &weight);
    if((pos_l[0]<=r0) || (pos_l[1]<=r0) || (pos_l[2]<=r0))	return 0;
  }
  while ( (weight > weight_min) && (energy > emin_o) );  	
		
return 0;		
}

int raytrace_in_disk(pos, mom, energy, weight)
// pos[], mom[] are 3-vector in cartesian coordinate, in disk-rest frame
double pos[], mom[], energy, *weight;
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
    
	raytrace_infinity(pos, mom, energy, wesc);    
    
  }
  lambda = -log(1-(1-p_esc)*rnd())/sigma_tot;
//  printf("lambda = %e\n", lambda);
//  getchar();
  for (i=0; i<3; i++)
    pos[i] += lambda*mom[i];                 /* position of the next    
						scattering               */  
  (*weight) *= (1-p_esc)*sigma_es/sigma_tot;
    
    return 0;
}

int raytrace_infinity(pos, dir, energy, weight)
// pos[], dir[] are 3-vector in cartesian coordinate, in disk-rest frame
double pos[], dir[], energy, weight;
{
	double kx=dir[0], ky=dir[1], kz=dir[2];
    // make sure [kx,ky,kz] is a unit vector
    double kk = sqrt(kx*kx + ky*ky + kz*kz);
    kx /= kk; ky /= kk; kz /= kk;

    // convert cartesian coordinates to spherical; using m=cos(theta)
    double x=pos[0], y=pos[1], z=pos[2];
    double r = sqrt(x*x + y*y + z*z);
    double m = z/r;
    double phi = atan2(y,x);

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

/*    #ifdef DEBUG
    // debug reporting
    printf("gd_P=%e  (Prp=%e)\n", gd_P, gd.Rpa);
    printf("cos_inc=%e (%.2f)\n", gd.cos_i, acos(gd.cos_i)/3.1415*180);
    printf("phi=%f + %f = %f\n", phi/3.1415*180, dphi/3.1415*180, reduce_angle_2pi(phi + dphi)/M_PI*180);
    #endif
*/
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
/*
    // independent check by integrating geodecis equation
    raytrace_data rtd;
    double position[4];
    position[0] = 0.0;position[1] = r;position[2] = m;position[3] = phi;
    raytrace_prepare(bh_spin, position, k, NULL, 0.001, 0, &rtd);
    while (1) {
        double dl = 1e9; // use maximal step
        raytrace(position, k, NULL, &dl, &rtd);
        // stop condition:
        if ((position[1] < r_bh(bh_spin)) || (position[1] > 1e9)) break;
        // also stop if relative error this step is too large
        if (rtd.error>1e-2) {
            fprintf(stderr, "raytrace(): aborted due to large error\n");
            break;
        }
    }
    double k_inf[4];
    bl2torus(position, k, k_inf);
    if (position[1] > 1e9) {
        fprintf(stderr, "ge-inc=%e  ge-phi=%e, g = %e\n", position[2], rad2deg(reduce_angle_2pi(position[3])), k_inf[0]);
    }
*/	

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
