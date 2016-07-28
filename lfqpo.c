/*

  Compton scattering in a uniform, hot  plasma cloud.
  Based on Gorecki & Wilczewski and Pozdnyakov, Sobol & Sunyaev

  This code is the combination of Piotr's (Comptonization) and Michal's (GR effect) codes.
  
  Currently, it only concentrates on the wedge geometry for Lense-Thirring precession project (LFQPO). 
  It will be extented to other geometries after this project.

*/

#include "mtc_incl_def.c"
#include "mtc_incl_code.c"
#include "quadrat.c"


double disk_photons_pdf(double x);
double disk_gfactor_K(double r, double a, double l, double q, double rms);




int main(int argc, char *argv[])
{
    long int il, istep;

    double    position[4], momentum_on[4],momentum_bl[4],
              r;
    int       runnum;
    int 	  i_inf =0, i_torus = 0;


  runnum = 0;
  if (argc == 2)   {
    runnum = atoi(argv[1]);
  }


/* initialize the geometry of:  1. the torus
                                2. the disk
                                3. the axis of the torus */

    input_data(runnum); // read input parameters and set up the geometry of the torus
    // double r_min = r_ms(bh_spin);
    double r_min = Rdisk;          // inner radius of the disk
    double r_max = 3000.0;         // outer radius of the disk 
    // setup NT disk
    disk_nt_setup(bh_mass, bh_spin, bh_mdot, 0.1, 0);

    // setup the axis of the torus with given (prectheta, precphi) in order to check if photon is inside of the torus be means of Dot Product
    axis_torus();
    //printf("vx=%.2e  vy=%.2e  vz=%.2e\n", axis_vec[0], axis_vec[1], axis_vec[2]);

    // setup photon distribution
    sim5distrib dist;
    distrib_init(&dist, &disk_photons_pdf, log(r_min), log(r_max), 100);

/* end of set-up of the torus-disk geometry */


  istep = mc_steps/5;

  for (il=0; il<=mc_steps; il++) { // generate photons

        sim5tetrad t;
        sim5metric m;
        geodesic gd;
        int status;

        // draw an emission radius from distribution        
        r = exp(distrib_hit(&dist));
        //fprintf(stderr, "r=%.1f\n", r);

     	// generate pos[4]
	    position[0] = 0.0;
	    position[1] = r;
	    position[2] = 0.0;
	    position[3] = urand*M_PI*2.0;

        // get metric and tetrad
        kerr_metric(bh_spin, r, 0.0, &m);
        double Omega = OmegaK(r,bh_spin);
        tetrad_azimuthal(&m, Omega, &t);
        //fprintf(stderr, "Omega=%.2e\n", Omega);

        // pick direction
        double th0 = urand*M_PI/2.;
        double ph0 = urand*M_PI*2.;
        momentum_on[0] = 1.0;
        momentum_on[1] = sin(th0)*cos(ph0);
        momentum_on[2] = cos(th0);
        momentum_on[3] = sin(th0)*sin(ph0);
        on2bl(momentum_on, momentum_bl, &t);
        //fprintf(stderr, "k[0]=%.2e  k[1]=%.2e  k[2]=%.2e  k[3]=%.2e\n", k[0], k[1], k[2], k[3]);
        //fprintf(stdout, "k[0]=%.2e  k[1]=%.2e  k[2]=%.2e  k[3]=%.2e\n", k[0], k[1], k[2], k[3]);
        //fprintf(stdout, "p[0]=%.2e  p[1]=%.2e  p[2]=%.2e  p[3]=%.2e\n", position[0], position[1], position[2], position[3]);
        
        // get geodesic
        // **************************************************************
        // IMPORTANT:
        // geodesic has to be initialized first before ray-tracing
        // **************************************************************
        geodesic_init_src(bh_spin, r, 0.0, momentum_bl, momentum_bl[1]<0?1:0, &gd, &status); // check if the geodesic of photon from disk to infinity can be found, 
                                                                                             // ignoring the existence of the torus
        if ((!status) || (r < gd.rp)) {
            fprintf(stderr, "WRN: geodesic_init_src failed (status=%d, r=%.2e, rp=%.2e)\n", status, r, gd.rp);
            continue;
        }

        // get temperature at r, g-factor and pick photon energy in the local disk
        double T = disk_nt_temp(r);
        double g_disk_inf = disk_gfactor_K(r, bh_spin, gd.l, gd.q, r_min); 
        double E_loc = blackbody_photon_energy_random(T);                  


/* ==================== The block for ray-tracing ================================= */
	raytrace_data rtd;
	raytrace_prepare(bh_spin, position, momentum_bl, NULL, 0.001, 0, &rtd);
	int outside_torus = 1; // first initialize outside_torus =1
	while (1) {
	    double dl = 1e9; // use maximal step
	    raytrace(position, momentum_bl, NULL, &dl, &rtd);
	    	   
	    // stop condition:
	    if ((position[1] < r_bh(bh_spin)) || (position[1] > 1e9)) break;
	    // also stop if relative error this step is too large
	    if (rtd.error>1e-2) {
	        fprintf(stderr, "raytrace(): aborted due to large error\n");
	        break;
	    }
	    
		int inside_torus = inside_dotprod(position); 
		if(inside_torus) { // if inside_torus = 1, then photon is in the torus
		i_torus++;	
	    printf("Inside of torus: %d\n", i_torus);
	    outside_torus = 0; 
	    break; // stop ray-tracing of photon which goes into torus
	    }	
	} // end of ray-tracing of each photon
	
/* ==================== The block for ray-tracing ================================= */

/* -------------------- The block for saving photon energy and position ---------------------------------- */
	if(outside_torus == 1) { // if after ray-tracing, outside_torus keep to be 1, this means that photon doesn't go into torus
		i_inf++; 	
		double E_inf = g_disk_inf*E_loc;
	    // get final theta-position at infinity
        double theta_inf = acos(gd.cos_i);	
	} else { // when outside_torus = 0, photon goes into torus. we should start the Comptonization
	    
        continue;
	}
	
	printf("photon going to infinity: %d\n", i_inf);

/* -------------------- The block for saving photon energy and position ---------------------------------- */


  } // end of all photons 
/* ***************************************************************************************************************** */

    distrib_done(&dist);


    return 0;
} // end of main program 


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
    
    return blackbody_photons_total(T, hardf, cos_mu) * r*r;
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

/* setup the axis of the torus with given (prectheta, precphi) in order to check if photon is inside of the torus by means of Dot Product of unit vector and photon position vector*/
/* suppose the precession starting from positive x-axis, to x=0, to negative x-axis */
    axis_torus()
{
	double prec_theta, prec_phi;
	prec_theta = prectheta;
	prec_phi = precphi;
	axis_vec[0] = cos(prec_theta);                  // z-component of the unit vector of torus axis
	axis_vec[1] = sin(prec_theta)*cos(prec_phi);    // x-component of the unit vector of torus axis
	axis_vec[2] = sin(prec_theta)*sin(prec_phi);    // y-component of the unit vector of torus axis
	
    return 0;
	
}


/* convert photon postion from B-L to cartesian coordinate */
	bl2cart(position_bl, position_cart)
	double position_bl[], position_cart[];
{
	position_cart[0] = position_bl[0];
	position_cart[1] = position_bl[1]*fabs(sqrt(1.-sqr(position_bl[2])))*sin(position_bl[3]);   // x
	position_cart[2] = position_bl[1]*fabs(sqrt(1.-sqr(position_bl[2])))*cos(position_bl[3]);   // y
	position_cart[3] = position_bl[1]*position_bl[2];                                            // z
	
	return 0;
	
}


/* check if photon is inside torus  */
   inside_dotprod(position_bl)
    double position_bl[];
{
	double  position_cart[4];
	int inside_dot, inside_r, inside_theta;

	bl2cart(position_bl, position_cart);
	double dot_prod = axis_vec[0]*position_cart[1] + axis_vec[1]*position_cart[2] + axis_vec[2]*position_cart[3];
	
	inside_r = position_bl[1]>Rin && position_bl[1]<Rout ; // to chick if photon is inside torus in r-direction
	inside_theta = fabs(dot_prod/position_bl[1]) <= fabs(sin(torus_theta0)); // to chick if photon is inside torus in theta-direction
	inside_dot = inside_r*inside_theta;
		
	return inside_dot;

}
