/*
Author: Bei You (youbeiyb@gmail.com)

1. This program is to tabulate the terminus of individual photon from a half-meridian along the sphere by use of raytracing code.

2. For gicen radius r and spin, the grid of the table model for the half-ring (assuming the longitude is 0.0) is determined by position (m = costheta) and momentum(th0, ph0).

3. The grid of the table model contains the following information: position[4] and momentum[4];

4. If the photon goes to either infinity or the disk surface.

5. Therefore, for any given photon with (pos[], mom[]) escaping from the sphere with radius r, we can derive the terminus by use of 3-D interpolation, no need for raytracing,
   which can significantly reduce the computing time.
--------------------------------------------------------------------
*/

#include "mtc_incl_def.c"
#include "mtc_incl_code.c"
#include "quadrat.c"
#include "sim5lib.c"

#define vect_copy(v1,v2) {v2[0]=v1[0]; v2[1]=v1[1]; v2[2]=v1[2]; v2[3]=v1[3];}

long int Ntot = 0, Ntot1=0, Ntot2=0, Ntot3=0, Ntot4=0;
int i_bh=0, i_step=0, iter;

int main(int argc, char *argv[])
{
int runnum;

  runnum = 0;
  if (argc == 2)   {
    runnum = atoi(argv[1]);
  }

/* initialize the geometry of:  1. the torus
                                2. the disk
                                3. the axis of the torus */
    input_data(runnum); // read input parameters and set up the geometry of the torus
      
    sphere_raytrace(); // calculate the table model   
      
printf("N = %ld, %ld, %ld, %ld, %ld, %d\n", Ntot, Ntot1, Ntot2, Ntot3, Ntot4, iter);  
printf("i_bh = %d, i_step=%d\n", i_bh, i_step);                      			 								 
	double speed;						
	speed = (npsp+1)*(ntsp+1)*(nfsp+1)/(time(NULL)-start_time);
    printf(" Speed of computations: %.1f photons/sec \n",speed);
		    					
    return 0;
} // end of main program 

/*  ***********************************************************************************************************************  */
int sphere_raytrace()
{
//	double r_min = ir*Rdisk;     // inner radius of the disk
	double r_min = Rout;     // inner radius of the disk
    double r_max = Rmax;         // outer radius of the disk    
//---------------------------------------------------------------------------------
double r_sp = r_min; // define the radius of the sphere; the longitude of the meridian is assumed to be 0.0
double r = r_sp;

int ii, ij, ik, it;
double dcos=(cossp_max-cossp_min)/npsp,
	   dkt=(ktsp_max-ktsp_min)/ntsp,
	   dkf=(kfsp_max-kfsp_min)/nfsp;  

FILE *data = fopen("sphere_table.dat", "w"); // opens new file for writing

for (ii=0; ii<=npsp; ii++) {	
	for (ij=0; ij<=ntsp; ij++) {
		for(ik=0; ik<=nfsp; ik++){
			for(it=0; it<=7; it++){
					sp_ph[ii][ij][ik][it] = 0.0; 
			}
		}	
	}
}

double m,th0,ph0;

for (ii=0; ii<=npsp; ii++) {	
	for (ij=0; ij<=ntsp; ij++) {
		for(ik=0; ik<=nfsp; ik++){
	iter++;
	if(iter%10000==0) printf("iter=%d\n", iter);
	m=cossp_min+ii*dcos;
	th0=acos(ktsp_min+ij*dkt);
	ph0=kfsp_min+ik*dkf;
		
    // get metric coefficients
    sim5metric metric;
    kerr_metric(bh_spin, r, m, &metric);

    // create ZAMO frame and renormalize the k vector in it,
    // then convert it back to B-L frame
    double k_loc[4], k[4], kk;
    sim5tetrad t;
    tetrad_zamo(&metric, &t);
	k_loc[1] = sin(th0)*cos(ph0);
    k_loc[2] = cos(th0);
    k_loc[3] = sin(th0)*sin(ph0);  
    kk = sqrt(sqr(k_loc[1]) + sqr(k_loc[2]) + sqr(k_loc[3]));
    k_loc[0] = 1.0; k_loc[1] /= kk; k_loc[2] /= kk; k_loc[3] /= kk;    
    on2bl(k_loc, k, &t);
  
	Ntot++; //count how many photons's r-direction is position, which corresponds that photon goes to infinity	     
    raytrace_data rtd;
    double pos[4], costi, costf, dl;
	pos[0]=0.0, pos[1]=r, pos[2]=m, pos[3]=0.0;
    raytrace_prepare(bh_spin, pos, k, NULL, 0.1, 0, &rtd);
	double p0i=0.0, p1i=0.0, p2i=0.0, p3i=0.0,
		   p0f=0.0, p1f=0.0, p2f=0.0, p3f=0.0,
		   k0i=0.0, k1i=0.0, k2i=0.0, k3i=0.0,
		   k0f=0.0, k1f=0.0, k2f=0.0, k3f=0.0;
	do{ //printf("pos=%e, %e, %e\n", pos[1], pos[2], pos[3]);
        dl = 1e9; // use maximal step
        costi = pos[2];
		p0i=pos[0], p1i=pos[1], p2i=pos[2], p3i=pos[3],
		k0i=k[0], k1i=k[1], k2i=k[2], k3i=k[3];
		
        raytrace(pos, k, NULL, &dl, &rtd);        

        costf = pos[2];
		p0f=pos[0], p1f=pos[1], p2f=pos[2], p3f=pos[3],
		k0f=k[0], k1f=k[1], k2f=k[2], k3f=k[3];

		if(isnan(pos[1]) || isnan(pos[2]) || isnan(pos[3])){ // failed geodesic
			Ntot2++;	
//			printf("wrong geodesic\n");
				sp_ph[ii][ij][ik][1] = -1.0;
				for(it=0; it<=7; it++){
					fprintf(data, "%e ", sp_ph[ii][ij][ik][it]);
				}
					fprintf(data, "\n");			
			break;
		}
//printf("pos = %e, %e, %e, %e\n", pos[0], pos[1], pos[2], pos[3]);
//printf("k = %e, %e, %e, %e\n", k[0], k[1], k[2], k[3]);		
/*        if((k[1]<=0.0)||(pos[1]<=r_min)){ // from here, we learn that once photon's r-velocity is position, it will continue to infinity, 
										  // never turn to negative ones and travel into the sphere
			printf("should not\n");
			getchar();
		}	
*/        
        // stop condition:
        if (pos[1] < r_bh(bh_spin)){ 
			i_bh++;
//			printf("go to BH\n");
//			printf("should not\n");
//			getchar();
				sp_ph[ii][ij][ik][1] = -1.0;
				for(it=0; it<=7; it++){
					fprintf(data, "%e ", sp_ph[ii][ij][ik][it]);
				}
					fprintf(data, "\n");
			break;		
		}
        // also stop if relative error this step is too large
        if (rtd.error>1e-2) {
//			printf("r=%e, m=%e, th0=%e, ph0=%e\n", r, m, th0, ph0);
//			getchar();
//			printf("raytrace(): aborted due to large error\n");	
				sp_ph[ii][ij][ik][1] = -1.0;
				for(it=0; it<=7; it++){
					fprintf(data, "%e ", sp_ph[ii][ij][ik][it]);
				}
					fprintf(data, "\n");		
			i_step++;		
            break;
        }
/* 
the following formula is used to estimate pos[] and k[] on the disk for reflection 
y0 = y1 + [(0-cost1)/(cost2-cost1)]*(y2-y1)
*/        
        if((pos[1]<=1e6)&&((costi>0.0)&&(costf<0.0))){ // the ourer radius of the disk for reflection is 1e6		
				Ntot3++;
				sp_ph[ii][ij][ik][0] = p0i + (-costi)/(costf-costi)*(p0f-p0i);
			    sp_ph[ii][ij][ik][1] = p1i + (-costi)/(costf-costi)*(p1f-p1i);
			    sp_ph[ii][ij][ik][2] = p2i + (-costi)/(costf-costi)*(p2f-p2i);
			    sp_ph[ii][ij][ik][3] = p3i + (-costi)/(costf-costi)*(p3f-p3i);
			    sp_ph[ii][ij][ik][4] = k0i + (-costi)/(costf-costi)*(k0f-k0i);	
			    sp_ph[ii][ij][ik][5] = k1i + (-costi)/(costf-costi)*(k1f-k1i);
			    sp_ph[ii][ij][ik][6] = k2i + (-costi)/(costf-costi)*(k2f-k2i);
			    sp_ph[ii][ij][ik][7] = k3i + (-costi)/(costf-costi)*(k3f-k3i);		  		    		    			    			    
				for(it=0; it<=7; it++){
					fprintf(data, "%e ", sp_ph[ii][ij][ik][it]);
				}
					fprintf(data, "\n");
//		printf("m=%e, th0=%e, ph0=%e\n", m, th0, ph0);
//		getchar();			
			break; 	
		}
			
		if(pos[1]>1e9){ // define the infinity to be 1e9
//			if(pos[2]>0.0){
				Ntot1++;
				sp_ph[ii][ij][ik][0] = pos[0];
			    sp_ph[ii][ij][ik][1] = pos[1];
			    sp_ph[ii][ij][ik][2] = pos[2];
			    sp_ph[ii][ij][ik][3] = pos[3];
			    sp_ph[ii][ij][ik][4] = k[0];	
			    sp_ph[ii][ij][ik][5] = k[1];
			    sp_ph[ii][ij][ik][6] = k[2];
			    sp_ph[ii][ij][ik][7] = k[3];
//			}
//			if(pos[2]<=0.0){
//				 Ntot4++;
//				 sp_ph[ii][ij][ik][1] = -1.0;    					
//			}
			for(it=0; it<=7; it++){
				fprintf(data, "%e ", sp_ph[ii][ij][ik][it]);
			}
				fprintf(data, "\n");			
			break; // stop raytracing for photons which go to infinity			       
		}
    }while(1);
//printf("pos=%e, %e, %e, %e\n", pos[0], pos[1], pos[2], pos[3] );
//printf("k=%e, %e, %e, %e\n", k[0], k[1], k[2], k[3]);   
//getchar();	 
		}	
	}
}  
//---------------------------------------------------------------------------------
return 0;
}	
