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

int i_spec = 0; /* 0: accumulate photon weight as direct disk spectrum 
                   1:                             reflection spectrum
                   2:                             compton spectrum */

int main(int argc, char *argv[])
{
    long int il, istep;

    double    position[4], momentum[4], energy, weight;
    int       runnum;
    int 	  indx, iindx, findx;

  runnum = 0;
  if (argc == 2)   {
    runnum = atoi(argv[1]);
  }


/* initialize the geometry of:  1. the torus
                                2. the disk
                                3. the axis of the torus */

    input_data(runnum); // read input parameters and set up the geometry of the torus
    // double r_min = r_ms(bh_spin);
        
//  istep = mc_steps/5;

int i_path; /* photon path by raytrace:
				0: to infinity; 1: to disk (for reflection); 2: to torus; 3: to BH */
int not_scattering, i_weight = 1;


for (il=0; il<=mc_steps; il++) { // generate photons
	

	i_path = 0;
	not_scattering = 0;
	
	generate_photon_gr(position, momentum, energy);
	weight = 1;
	wghtmin = min_weight*weight;

	do{ // large loop: raytrace -> reflection (or scattering) -> raytrace		
		
		do{// one single run of raytrace-> infinity (or reflection, scattering)
			raytrace_out_torus(position, momentum, energy, i_path); /* determine whether photon from the disk goes to infinity,
			                                                                                                or reflection, or scattering, or BH */
			if(i_path=3)break; // to BH
			if(i_path=0)break; // to infinity
			if(i_path=1){ // this is the reflection
			     //reflection()	// do reflection; Will be completed soon
			     i_spec = 1; 	
			}
			if(i_path=2)break;
			
		}while(1);
	
		if(i_path=3){
					 weight = 0.0;
					 break;
					 } // to BH
		if(i_path=0)break; // to infinity
		if(i_path=2){ // to torus
			do{ // loop for multi-scattering in the torus
				not_scattering = raytrace_in_torus(position, momentum, energy, weight);  
				if(not_scattering){
					break;
				}else{
					i_spec = 2;
					if((weight < wghtmin)) { /* if photon weight is too small, ignore loop of scattering */
					i_weight = 0;
					break;	
					}else{continue;
					}
				}

			}while(1);	
			if(i_weight=0) {weight = 0.0;
						break; // photon is highly absorpted, ignore 
			}	
		}						
	}while(1); // only when i_path = 0 or 3, exit raytrace loop to next stage: store photon energy and weight into arrays 
	
/* -------------------- The block for saving photon energy and position ---------------------------------- */
		int i;
	if((energy > emin_o) && (energy < emax_o)){ // only deal with photon energy in the interesting range 	
		int indx = num_bin_o*log10(energy/emin_o); // indx corresponding to photon energy at infinity 
		
	    // get final theta-position at infinity
        double theta_inf = position[2];	// cos(theta)
        
        // get final phi-position at infinity
        double phi_inf = position[3];
        
		/* bins in \phi: 0-20, 90-110, 180-200, 270-290 degs */
	    findx = -1;
	    for(i=0; i<MAXphi; i++) {
	      if (fabs(phi_inf-(i*90+10)) <= 10) {  
	       findx = phi_inf/90.;
	      }
	    }
	    
	    int iindx = theta_inf*nangle;
	    
	    switch (i_spec){
	    case 0: {out_disk[indx][iindx][findx] += weight;
				 break;}
		case 1: {out_refl[indx][iindx][findx] += weight;
				 break;}
		case 2: {out_comp[indx][iindx][findx] += weight;
				 break;}	 
		}                
	}
/* -------------------- The block for saving photon energy and position ---------------------------------- */

	} // end of all photons 

	output_spec(runnum,il); /* output spectrum: Novikov-Thorne disk spectrum;
					                   Disk reflection spectrum;
					                   Torus Comptonization spectrum;
					*/
    return 0;
} // end of main program 

output_spec(runnum,noph)
  long int   noph;
  int        runnum;
{
  FILE   *fp;
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
		
fp = fopen(inp_spec,"w");         /* output observed Novikov-Thorne disk spectrum */
  nos = num_bin_o*log10(emax_o/emin_o);
  for (i=0; i<nos; i++){
    e = emin_o*pow(10.,(i+0.5)/num_bin_o );
    de = emin_o*(pow(10.,(i+1.0)/num_bin_o)-pow(10.,(1.0*i)/num_bin_o));
    er = e/de/noph;
    for(j=0; j<nangle; j++) {
		for(k=0; k<MAXphi; k++) {
			fprintf(fp,"%e %.3e ",e*511, out_disk[i][j][k]*er*nangle);
		}
		fprintf(fp,"  ");
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

fp = fopen(reflspec,"w");         /* output observed reflection spectrum from the disk */
  nos = num_bin_o*log10(emax_o/emin_o);
  for (i=0; i<nos; i++){
    e = emin_o*pow(10.,(i+0.5)/num_bin_o );
    de = emin_o*(pow(10.,(i+1.0)/num_bin_o)-pow(10.,(1.0*i)/num_bin_o));
    er = e/de/noph;
    for(j=0; j<nangle; j++) {
		for(k=0; k<MAXphi; k++) {
			fprintf(fp,"%e %.3e ",e*511, out_refl[i][j][k]*er*nangle);
		}
		fprintf(fp,"  ");
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

fp = fopen(out_file,"w");         /* output observed torus Comptonization spectrum */
  nos = num_bin_o*log10(emax_o/emin_o);
  for (i=0; i<nos; i++){
    e = emin_o*pow(10.,(i+0.5)/num_bin_o );
    de = emin_o*(pow(10.,(i+1.0)/num_bin_o)-pow(10.,(1.0*i)/num_bin_o));
    er = e/de/noph;
    for(j=0; j<nangle; j++) {
		for(k=0; k<MAXphi; k++) {
			fprintf(fp,"%e %.3e ",e*511, out_comp[i][j][k]*er*nangle);
		}
		fprintf(fp,"  ");
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

return 0;
}


int generate_photon_gr(position, momentum, E_loc)
double position[], momentum[], E_loc; 
{
	
	double r_min = Rdisk;          // inner radius of the disk
    double r_max = 3000.0;         // outer radius of the disk 
    // setup NT disk
    disk_nt_setup(bh_mass, bh_spin, bh_mdot, 0.1, 0);

    /* convert the axis of the torus (0, 0, 1) in torus frame to Cartesian coordiante of BH, according to (prectheta, precphi), 
       in order to check if photon is inside of the torus by means of Dot Product */
    if (precession_option<2) {
      torus_axis1(); // precession model (1)
    } else {
      torus_axis2(); // precession model (2)
    } 
	printf("%e %e %e\n", axis_vec[0], axis_vec[1], axis_vec[2]);


    // setup photon distribution
    sim5distrib dist;
    distrib_init(&dist, &disk_photons_pdf, log(r_min), log(r_max), 100);

        sim5tetrad t;
        sim5metric m;

/* end of set-up of the torus-disk geometry */
	
		double r, T;
        // draw an emission radius from distribution        
        do {
        r = exp(distrib_hit(&dist));
        T = disk_nt_temp(r); 
		} while (2.7*T < (emin*511.)); 


     	// generate pos[4] in B-L spherical coordinate
	    position[0] = 0.0;
	    position[1] = r;
	    position[2] = 0.0;
	    position[3] = urand*M_PI*2.0;
    
        // get metric and tetrad
        kerr_metric(bh_spin, r, 0.0, &m);
        double Omega = OmegaK(r,bh_spin);
        tetrad_azimuthal(&m, Omega, &t);
        //fprintf(stderr, "Omega=%.2e\n", Omega);

        // pick direction in disk rotating frame
        //double th0 = urand*M_PI/2.;
        double th0 = acos(sqrt(1.-sim5urand())); // iCDF that corresponds to PDF(theta)=sin(theta)*cos(theta)
        double ph0 = urand*M_PI*2.;
        double momentum_on[4];
        momentum_on[0] = 1.0;
        momentum_on[1] = sin(th0)*cos(ph0);
        momentum_on[2] = cos(th0);
        momentum_on[3] = sin(th0)*sin(ph0);
				
        on2bl(momentum_on, momentum, &t); // transfer the direction from disk frame to B-L
               
        /*
		Generate energy from a planckian distribution of a corresponding 
		temperature
		*/
		do {
		  E_loc = blackbody_photon_energy_random(T)/511.; /* photon energy in disk-rest frame */
		} while ( (E_loc < emin) || (E_loc > emax));

    distrib_done(&dist);

return 0;
}


int raytrace_out_torus(position, momentum, energy, i_path) 
/* a function of raytracing photon which is outside torus */
double position[], momentum[], energy;
int i_path;
{
		double posi[4], momi[4], posf[4], momf[4];
		double dl, pu_i, pu_f;
		int inside_torus;
		           
        pu_i = prod_pU(position, momentum); /* dot product between photon momentum and medium 4-velocity, at the first step */
                       
        /* ==================== The block for ray-tracing ================================= */
	    raytrace_data rtd;
	    raytrace_prepare(bh_spin, position, momentum, NULL, 0.001, 0, &rtd);

		do{
	    dl = 1e9; // use maximal step
	    
	    /* copy input position[] and momentum[] as initial ones at each step */
	    vect_copy(position, posi); // position -> posi 
        vect_copy(momentum, momi); // momentum -> momi     
	    
	    raytrace(position , momentum, NULL, &dl, &rtd); 
	    // after raytrace(), posi and momi already changed
	    
	    		// stop condition: if applied, stop raytrace
				if ((position[1] < r_bh(bh_spin))) {
				printf("raytrace(): photon goes to bh event horizon\n");
				i_path = 3;  
				break;  
				}
				// also stop if relative error this step is too large
				if (rtd.error>1e-2) {
				printf("raytrace(): aborted due to large error\n");
				i_path = 3; // failure step; In order to ignore this photon, regarding such photon as one going to BH event horizon, 
				break;
				}
	    
	    /* position and momentum after each step */
	    vect_copy(position, posf); // posi -> posf
        vect_copy(momentum, momf); // momi -> momf
        
        if(posf[1]>1e9) { // Is this value large enough to stop raytrace ? 
			i_path = 0;
			break;
		}	
        
        if((posi[2]>0.0)&&(posf[2]<0.0)){ // chech if photon go across the disk
			i_path = 1;
			break;
		}	 
        
        inside_torus = inside_dotprod(posf);	
		if(inside_torus){
			i_path = 2;
			break;
		}	
        		
		}while(1);			
		/* ==================== End of the block for ray-tracing ================================= */
		
		pu_f = prod_pU(posf, momf); /* dot product between photon momentum and medium 4-velocity, at the last step */   		
		energy *= pu_f/pu_i;	
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
	
	inside_r = position_bl[1]>Rin && position_bl[1]<Rout ; // to chick if photon is inside torus in r-direction
	inside_theta = fabs(dot_prod/position_bl[1]) <= fabs(sin(torus_theta0)); // to chick if photon is inside torus in theta-direction
	inside_dot = inside_r*inside_theta;
		
	return inside_dot; // 1: inside of torus; 0: outside of torus

}


  prod_pU(position, momentum) /* dot product between photon momentum and medium 4-velocity */
double position[], momentum[];
{
		double U[4],Omega,ell_torus, pU;

        sim5tetrad t; /* Question: will this t, m, gd conflict with main program ? */
        sim5metric m;
                
        // get metric and tetrad
        kerr_metric(bh_spin, position[1], position[2], &m);
        ell_torus = ellK(position[1], bh_spin); // torus specific angular momentum

		Omega = Omega_from_ell(ell_torus, &m);
        //double Omega = OmegaK(r,bh_spin); // for accretion disk
        
        tetrad_azimuthal(&m, Omega, &t);
              
        // dot product of photon direction and medium motion: p[4]*U[4]
        fourvelocity_azimuthal(Omega, &m, U[4]);
        pU = dotprod(momentum, U, &m);

	return pU;
}        

int bl2torus(position_bl, momentum_bl, momentum_torus) 
/* convert photons momentum from B-L to torus-rest frame */
double position_bl[], momentum_bl[], momentum_torus[];
{
	double momentum_loc[4], diri[3], dirf[3], dirm[3];
	
	/* first convert momentum from B-L to torus-rest frame */
        sim5tetrad t; 
        sim5metric m;
                
        // get metric and tetrad
        kerr_metric(bh_spin, position_bl[1], position_bl[2], &m);
        double ell_torus = ellK(position_bl[1], bh_spin); // torus specific angular momentum
         
		double Omega = Omega_from_ell(ell_torus, &m);
        //double Omega = OmegaK(r,bh_spin); // for accretion disk
        
        tetrad_azimuthal(&m, Omega, &t);              
		bl2on(momentum_bl, momentum_loc, &t);
	
	/* secondly, in torus-rest frame, taking precession into account */
	    diri[0] = momentum_loc[1];
	    diri[1] = momentum_loc[2];
	    diri[2] = momentum_loc[3];

		double R[3][3];
		double prec_theta = prectheta;
		double prec_phi = precphi;	    
		
	        if (precession_option<2) {
			    rotmatrix(R,prec_phi,prec_theta);
				multmatvec(dirf,R,diri);
				
				// Question: don't know how to deal with t-component (momentum_loc[0]), in the context of precession
				momentum_torus[0] = dirf[0];
				momentum_torus[1] = dirf[1];
				momentum_torus[2] = dirf[2];
			
			} else {
			    rotmatrix(R,0.0,prec_theta);
				multmatvec(dirm,R,diri);
				
				rotmatrix(R,prec_phi,prec_theta);
				multmatvec(dirf,R,dirm);
				
				// Question: don't know how to deal with t-component (momentum_loc[0]), in the context of precession
				momentum_torus[0] = dirf[0];
				momentum_torus[1] = dirf[1];
				momentum_torus[2] = dirf[2];
			}

	return 0;
}

int torus2bl(position_bl, momentum_torus, momentum_bl) 
/* convert photons momentum from torus frame to B-L */
double position_bl[], momentum_bl[], momentum_torus[];
{
	double momentum_loc[4], diri[3], dirf[3], dirm[3];
	
	/* first convert momentum from torus-rest frame to B-L */
        sim5tetrad t; 
        sim5metric m;
                
        // get metric and tetrad
        kerr_metric(bh_spin, position_bl[1], position_bl[2], &m);
        double ell_torus = ellK(position_bl[1], bh_spin); // torus specific angular momentum
         
		double Omega = Omega_from_ell(ell_torus, &m);
        //double Omega = OmegaK(r,bh_spin); // for accretion disk
        
        tetrad_azimuthal(&m, Omega, &t);              
		on2bl(momentum_torus, momentum_loc, &t);
	
	/* secondly, in B-L frame, taking precession into account */
	    diri[0] = momentum_loc[1];
	    diri[1] = momentum_loc[2];
	    diri[2] = momentum_loc[3];

		double R[3][3];
		double prec_theta = prectheta;
		double prec_phi = precphi;	    
		
	        if (precession_option<2) {
			    invrotmatrix(R,prec_phi,prec_theta);
				multmatvec(dirf,R,diri);
				
				// Question: don't know how to deal with t-component (momentum_loc[0]), in the context of precession
				momentum_bl[0] = dirf[0];
				momentum_bl[1] = dirf[1];
				momentum_bl[2] = dirf[2];
			
			} else {
			    invrotmatrix(R,0.0,prec_theta);
				multmatvec(dirm,R,diri);
				
				invrotmatrix(R,prec_phi,prec_theta);
				multmatvec(dirf,R,dirm);
				
				// Question: don't know how to deal with t-component (momentum_loc[0]), in the context of precession
				momentum_bl[0] = dirf[0];
				momentum_bl[1] = dirf[1];
				momentum_bl[2] = dirf[2];
			}

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

   return dtau;
}

int raytrace_in_torus(position,momentum,energy, weight)  
/* assuming photon is inside torus, given position[] and momentum[], to determine:
1. accumulated escaping probalibity along geodesic to boundary
2. if not_scattering = 1 (escaping),   posf[], momf[]: position and momentum at the surface of the torus;   																     
3. if not_scattering = 0 (scattering), posf[], momf[]: the momentum after scattering, at scattering postion 
4. if not_scattering = -1, raytrace failure													     
5. in any case, new energy and weight */																                                                     
double position[],momentum[],energy,weight;
{
double  posi_orig[4], momi_orig[4];
double  dl, p_sc, p_esc, dtau, tau;
int inside_torus;
double energy1, energy2, energy_m;
double pU_i, pU_f;
                  
        /* for the sake of not changing posi and momi, copy vectors to another temporary ones which are used in raytrace */
        vect_copy(position, posi_orig); // posi -> posi_orig
        vect_copy(momentum, momi_orig); // momi -> momi_orig
                
        list_init(); // initialize list_node
        
        /* ==================== The block for ray-tracing ================================= */       
	    raytrace_data rtd;
	    raytrace_prepare(bh_spin, posi_orig, momi_orig, NULL, 0.001, 0, &rtd);

		dl   = 0.0; // traveling distance of each step
		tau  = 0.0;
		dtau = 0.0;
		p_sc = 0.0; // accumulated scattering probalibity: p_sc = 1 - exp(-tau), where tau += sigma*dl 
		p_esc = 0.0;
		
		energy1 = energy; /* energy in torus-rest frame */
        energy2 = energy;
        energy_m = energy;

		while(1){ // loop to save each step along geodesic						
			
			inside_torus = inside_dotprod(posi_orig);	
			if(inside_torus){					
			tau += dtau; // accumulate traveling distance before escaping torus
			p_sc = 1.0 - exp(-tau);
			p_esc = exp(-tau);
            list_node_add(posi_orig , momi_orig, p_sc, energy2);  /* starting from the first step, save data into link-node; 
														   If, after the first step, photon is outside of torus, assuming no scattering */							   
			}else{		
			break; // once outside of torus, jump out of raytrace loop	
			}	
			
			pU_i = prod_pU(posi_orig , momi_orig);
			raytrace(posi_orig , momi_orig, NULL, &dl, &rtd); 
			// after each run of raytrace(), posi and momi already changed, although keep the same name	
			pU_f = prod_pU(posi_orig , momi_orig);
			
			energy2 = energy1*(pU_f/pU_i);
			energy_m = 0.5*(energy1 + energy2);
			dtau = ray_tau(dl, energy_m); 
			energy1 = energy2;	// energy1 as the initial energy for the next step	  
	    	   
			// stop condition: if applied, stop raytrace
			if ((posi_orig[1] < r_bh(bh_spin))) {
			printf("raytrace(): photon goes to bh event horizon\n");
			return -1;    
			}
			// also stop if relative error this step is too large
			if (rtd.error>1e-2) {
	        printf("raytrace(): aborted due to large error\n");
	        return -1;
			}
					
		}	
		/* ==================== The block for ray-tracing ================================= */

        
/* =============== to determine the scattering position, if photon is scattered =================*/

    double posmax[4],posmin[4],
		   mommax[4],mommin[4],
		   pmin, pmax, elmax, elmin;
	int i;
	int not_scattering = 0;
	double x;

		double rnd_sc = rnd();
		if(p_sc > rnd_sc){ // scattering part					
					x = rnd()*p_sc; // critical random value to determine the position where scattering happens 

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
					/* do interpolation to determine new pos[], mom[], energy and weight, corresponding to x */
					for (i=0; i<4; i++){
					position[i] = posmin[i] + (x-pmin)*(posmax[i]-posmin[i])/
												   (pmax-pmin); 
					momentum[i] = mommin[i] + (x-pmin)*(mommax[i]-mommin[i])/
												   (pmax-pmin);  
					}
					energy = elmin + (x-pmin)*(elmax-elmin)/(pmax-pmin); 							
					scattering_gr(position, momentum, energy); // this is the scattering 				
					weight *= p_sc;
		}
		else{ // if no scattering, get the shifted energy and weight after escaping
					/* copy temorary vectors to final ones as output*/
					not_scattering = 1;
					vect_copy(posi_orig, position); // posi_orig -> position
					vect_copy(momi_orig, momentum); // momi_orig -> momentum
					energy = energy1;
					weight *= p_esc;
		}		
/* =============== end of determining the scattering position, if photon is scattered =================*/

return not_scattering; /* if return 0, stop ray-trace;
						  if return 1, photon escape outside torus without scattering */
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

int scattering_gr(position, momentum, energy) // scattering subroutine 
double position[], momentum[], energy;
{
	double v,mu,x,Omega_new[3],vCart[3],momf_torus[3],momentum_bl[3];
    double enr_p;
	
	bl2torus(position, momentum, momf_torus); // convert photon direction from B-L to torus-rest frame
	enr_p = energy; /* enr_p and energy are photon energy in torus-rest frame */
					
	generate_velocity(&v,&mu,&x,enr_p); /* electron's velocity           */

	sphere2cart(v,mu,momf_torus,vCart);   /*  Cartesian components of the electron's
														 velocity in coordinate system with 
														 z-axis along photon's momentum      */

	new_Omega(v,mu,x,vCart,momf_torus,Omega_new,&enr_p);  /* this is scattering -
																	    new direction & energy of photon      */
   
	energy = enr_p; // in torus-rest frame 										
	torus2bl(position, Omega_new, momentum_bl); // convert from torus-rest frame to B-L
	momentum[1] = momentum_bl[1]; // momentum[0] not changed, because don't know how to deal with it
	momentum[2] = momentum_bl[2];
	momentum[3] = momentum_bl[3];
	
	return 0;
}
