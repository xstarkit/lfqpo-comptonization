/* November 1998 
          - tests of F(x) and Inverse_F(x); seems that analytical approximation
	    gives inaccurate results, especially in optically thick limit

 */

scattering(pos,mom,energy,weight)
    double pos[],mom[],*energy,weight;

{
    double v,mu,x,
           Omega_new[3],
	   vCart[3];
    double olde;
    int    i,ii;

    double mom_p[3],bn,enr_p;

    if (nscatt < MAXsc1) nscatt ++;
    if (nscatt > nscatt_glob) nscatt_glob = nscatt;

    olde = *energy;

    /* bulk velocity effect */

    bn = 1-b_bulk*mom[2];
    mom_p[0] =  mom[0]/g_bulk /bn;
    mom_p[1] =  mom[1]/g_bulk /bn;
    mom_p[2] = (mom[2]-b_bulk)/bn;
    /* energy in the comoving frame */
    enr_p = (*energy)*g_bulk*(1-b_bulk*mom[2]);

    /*
    if (nscatt==1) {
      ii = num_bin_i*log10(enr_p/emin);
      inp_s[ii][1] += weight;
    }
    */

    generate_velocity(&v,&mu,&x,enr_p); /* electron's velocity           */

    sphere2cart(v,mu,mom_p,vCart);   /*  Cartesian components of the electron's
                                       velocity in coordinate system with 
				       z-axis along photon's momentum      */

    new_Omega(v,mu,x,vCart,mom_p,Omega_new,&enr_p);  /* this is scattering -
                                     new direction & energy of photon      */
   
    bn = 1+b_bulk*Omega_new[2];
    mom[0] = Omega_new[0]/g_bulk/bn;
    mom[1] = Omega_new[1]/g_bulk/bn;
    mom[2] = (Omega_new[2]+b_bulk)/bn;
    *energy = enr_p*g_bulk*(1+b_bulk*Omega_new[2]);

   /*  for (i=0; i<3; i++)
       mom[i]=Omega_new[i]; */     /* new direction of photon's mom. */
    energ2 = 2*(*energy);

    /*    printf("vel: %.3f, olde: %.2f   newe:  %.2f \n",v,olde*511,(*energy)*511); */

    return 0;
}

generate_velocity(v,mu,x,energy)
   double *v,*mu,*x,energy;
/*
   input:  energy
   output: v,mu,x
   generates electron velocity according to the distribution ro1',
   direction of its motion (mu=cos(theta)) and value of 'x' parameter

   Global variables: arrays x_cum,y_cum, nop_rho1[]
*/
{
     double t,x1,x2,Fx1,Fx2,Fx,e1,v1,v2;
     int    ix,ix1,ix2,ie;

/*  generate value of velocity from the cumulative distribution       */
     t = rnd();
     
     /*  for the photon energy */
     if (energy<rho1_ens[0]) {
       ix = nop_rho1[0]*t;
       c_hunt(x_cum[0],nop_rho1[0],t,&ix);
       *v = y_cum[0][ix]+(t-x_cum[0][ix])*(y_cum[0][ix+1]-y_cum[0][ix])/
           	                          (x_cum[0][ix+1]-x_cum[0][ix]);
     } else {
       if (energy>rho1_ens[max_roe]) {
	 ix = nop_rho1[max_roe]*t;
	 c_hunt(x_cum[max_roe],nop_rho1[max_roe],t,&ix);
	 *v = y_cum[max_roe][ix]+(t-x_cum[max_roe][ix])*
	   (y_cum[max_roe][ix+1]-y_cum[max_roe][ix])/
	   (x_cum[max_roe][ix+1]-x_cum[max_roe][ix]);
       } else {
	 c_hunt(rho1_ens,max_roe,energy,&ie);
	 e1 = (energy-rho1_ens[ie])/(rho1_ens[ie+1]-rho1_ens[ie]);
	 
	 ix1 = nop_rho1[ie]*t;
	 c_hunt(x_cum[ie],nop_rho1[ie],t,&ix1);
	 v1 = y_cum[ie][ix1]+(t-x_cum[ie][ix1])*
	   (y_cum[ie][ix1+1]-y_cum[ie][ix1])/
	   (x_cum[ie][ix1+1]-x_cum[ie][ix1]);

	 ix2 = nop_rho1[ie+1]*t;
	 c_hunt(x_cum[ie+1],nop_rho1[ie+1],t,&ix2);
	 v2 = y_cum[ie+1][ix2]+(t-x_cum[ie+1][ix2])*
	   (y_cum[ie+1][ix2+1]-y_cum[ie+1][ix2])/
	   (x_cum[ie+1][ix2+1]-x_cum[ie+1][ix2]);

	 *v = v1*(1-e1) + v2*e1;
       }
     }

     if (fabs(*v) >= 1) {
       printf(" v>c ! %.2f  e: %.3f \n",(*v),energy*511);
       if (*v > 0) 
         *v = 0.999999;
       else
         *v = -0.999999;
     }
/*  generate 'x'                         */
     x1x2_fun(*v,&x1,&x2);
     Fx1 = F(x1);
     Fx2 = F(x2);
     Fx = Fx1 + rnd()*(Fx2-Fx1);
     *x = Invers_F(Fx);
     if ( ((*x)-x1)*(x2-x1) < 0 ) /* inverse failed */
     {
       *x = x1 + (x2-x1)/(Fx2-Fx1)*(Fx-Fx1);  /* linear interpolation */
     }
/*  calculate 'mu'   G&W (A16)            */
     *mu = 1/(*v)*(1-0.5*(*x)/gamma(*v)/energy);
     if (fabs(*mu) > 1)
     {
       if (*mu>0) *mu=1-1e-6; else *mu=-1+1e-6;
     }
     return 0;
}

new_Omega(v,mu,x,vCart,Omega_old,Omega_new,energy)
   double v,mu,x,vCart[],Omega_old[],Omega_new[],*energy;
/*
    input :v,mu,x,vCart,Omega_old,energy
    output:Omega_new,energy

    generates the direction of the scattered photon (Appendix B of G&W)
*/
{
   double g,mu_prim,c_scat,y,x_prime,xxp,gama,
          v1[3];
   int    i;
   
   for (i=0; i<3; i++)
      v1[i] = vCart[i]/v;                     /* normalize v        */
   gama = gamma(v);
   do
   {
      g = 2*rnd()-1;
      mu_prim = (v + g)/(1+v*g);              /*   (B5)             */

      sphere2cart(1.,mu_prim,v1,Omega_new);   /*   (B6)             */
      c_scat = 0;               /* comp. cosine of scattering angle */
      for (i=0; i<3; i++)
	c_scat += Omega_old[i]*Omega_new[i];  
      x_prime = x/(1 + (*energy)/gama*(1-c_scat)/(1-v*mu_prim));/*  (B3)   */
      xxp = x_prime/x;
      y = 0.5*X(x,x_prime)*xxp*xxp;
   }
   while (rnd() > y);

   *energy= (*energy)*(1-v*mu)/(1-v*mu_prim+(*energy)/gama*(1-c_scat));  /* (3) */

   return 0;
}

sphere2cart(vec_len,mu,s_axis,cart_vec)
   double vec_len,mu,
          s_axis[],
          cart_vec[];
/* input:  vec_len,mu,
           s_axis   -  z-axis of the spherical coordinate system
   output: cart_vec

   Transforms vector of length 'vec_len' and cos(theta) = 'mu' from spherical
   to cartesian coordinates, assuming random polar angle phi  (A17) of G&W
*/
{
   double st,sph,cph,mu2,Oz2,p1,Ozc,phi;

   phi = two_pi*rnd();

   if (s_axis[2]==1)
   {
     st = sin(acos(mu));
     cart_vec[0] = vec_len*st*sin(phi);
     cart_vec[1] = vec_len*st*cos(phi);
     cart_vec[2] = vec_len*mu;
   }
   else
   {
     mu2 = 1-mu*mu;
     Oz2 = 1-s_axis[2]*s_axis[2];
     p1 = sqrt(mu2/Oz2);
     sph = sin(phi);
     cph = cos(phi);
     Ozc = s_axis[2]*cph;
     cart_vec[0] = vec_len*(s_axis[0]*mu-p1*(s_axis[0]*Ozc+s_axis[1]*sph));
     cart_vec[1] = vec_len*(s_axis[1]*mu-p1*(s_axis[1]*Ozc-s_axis[0]*sph));
     cart_vec[2] = vec_len*(s_axis[2]*mu+sqrt(mu2*Oz2)*cph);
   }

   return 0;
}

generate_photon(pos,mom,energy,weight)
   double pos[],mom[],*energy,*weight; 
/*
   generates new photon
*/
{
  double el,phi,st,ct,  r,r_hit;
  int    indx,ei;

  total_dist = 0;
  nscatt = 0;
  el = rnd()*elog;          /* uniform distribution of energies */
  *energy = emin*exp(el);   /* in LOG scale                     */
  energ2 = (*energy)*2;
  
  if ((spat_distr==3)) {
    generate_disk_with_sphere(pos,mom,energy,weight);
    energ2 = 2*(*energy);
    indx = num_bin_i*log10((*energy)/emin); /* store input spectrum */
    /*    if ((indx>=0) && (indx<MAXP) ) inp_s[indx] += *weight; */
    wghtmin = min_weight*(*weight);
    return 0;
  }

  if ((spat_distr==6)) {
    generate_disk_with_torus(pos,mom,energy,weight);
    energ2 = 2*(*energy);
    indx = num_bin_i*log10((*energy)/emin); /* store input spectrum */
    /*    if ((indx>=0) && (indx<MAXP) ) inp_s[indx] += *weight; */
    wghtmin = min_weight*(*weight);
    return 0;
  }
  
  switch (spectrum) {         /* 'weight' is photon number for this energy */
  case 0 : {
    /*	  *weight = Planck(*energy)/(*energy);  */
    do {
      *energy = Planck_photons(T_bb);
    } while ( (*energy < emin) || (*energy > emax) );
    *weight = 1;
    energ2 = 2*(*energy); 
    break;
  }
  case 1 : { 
    *weight = power_law(*energy)/(*energy);
    break;
  }
  case 2 : { 
    ei= el/elog*nop_spec;   /* linear interpolation        */ 
    *weight = spect[ei] + ((*energy)-energ[ei])/(energ[ei+1]-energ[ei])*(spect[ei+1]-spect[ei]);
    break;
  }
  
  case 3 : { 
    *weight = diskbb(*energy)/(*energy); 
    break;
  }
  case 4 : { /* monochromatic input photons */
    *weight = 1;
    *energy = T_bb;
    energ2 = 2*(*energy);
    break ;
  }
  }

  indx = num_bin_i*log10((*energy)/emin); /* store input spectrum */
#ifdef DEBUG
  if (debug_cond) printf("ini energy: %e  %d \n",(*energy),indx);
#endif
  if ((indx>=0) && (indx<MAXP) ) inp_s[indx][0] += *weight;
  wghtmin = min_weight*(*weight);
  
  switch (spat_distr) {
  case 0 :
  case 1 :
  case 4 : {
    ct = 2*rnd()-1;
    break ;
  }
  case 2 : {   /* slab */
    switch (irradiation) {
    case 0 : {             /* semi-isotropic      */
      ct = rnd();
      break;
    }
    case 1 : {             /* cosine - like   ??? */
      ct = sqrt(rnd());
      break;
    } 
    case 2 : {             /* beamed irradiation  */
      ct = cos(theta_irrad*pi/180);
      break;
    }
    case 3 : {       /* Chandrasekhar-Sobol distribution (equiv. to sqrt?? */
      ct = 0.5*(sqrt(1+8*rnd())-1);
      break;
    }
    }
    break ;
  }
  case 5 : {  /* active region - hemisphere H=R = Rout */
    
    do {
      r = 3*Rout*rnd();          /* distance of emission */

      switch (irradiation) {
      case 0 : { ct = rnd(); break; }
      case 1 : { ct = sqrt(rnd()); break; }
      case 2 : { exit(13); }
      case 3 : { ct = 0.5*(sqrt(1+8*rnd())-1); break; }
      }
      
      r_hit = r - Rout*sqrt(1/(ct*ct)-1.);   /* Rout is height here */
    } while  (fabs(r_hit) > Rout);           /* here Rout is Rout   */
    
    pos[0] = r_hit;
    pos[1] = 0;
    pos[2] = 0;
    
    break;
  }
  case 3 : 
  case 6 : {
    printf("we should not be here !\n");
    exit(13);
  }
  }

  st = sqrt(1-ct*ct);
  phi = two_pi*rnd();
  mom[0] = st*cos(phi);      /* cartesian components of momentum */
  mom[1] = st*sin(phi);
  mom[2] = ct;
  
  /* generate initial position; in case of a shell, position should be
     right at the inner edge of it, in the direction mom[] */
  initial_position(mom,pos);

  return 0;
}

double generate_disk_bb(r_emit)
     double *r_emit;
{
  double t_emit,energy;

  if (r_emission < 0.) { /* distribution of r  - "normal" case */
    do {
      *r_emit = fabs(r_emission)*Rdisk*pow(rnd(),-4.);  /* -4 !!! */ 
      /* in units of the outer radius */
      t_emit = T_bb*pow((*r_emit)/Rdisk,-0.75);
    } while (2.7*t_emit < emin);
  } else {     /* a ring */
    *r_emit = Rdisk*r_emission;
    t_emit = T_bb;
  }
  
  /*
    Generate energy from a planckian distribution of a corresponding 
    temperature
  */
  do {
    energy = Planck_photons(t_emit);
  } while ( (energy < emin) || (energy > emax));
  
  return energy;
}

/* ====== start of the torus section  ============================= */

generate_disk_with_torus(pos,mom,energy,weight)
     double pos[],mom[],*energy,*weight;
     /*
Generates seed photon from a disk external to a TORUS. 
Disk temperature varies as r^(-0.75).
      */
{
  double r_emit,ct,st,phi,rq,t;
  int    i,j,tindx,eindx,iter,hits;
  double posi[3],momi[3],ipos[3];
  int    ins1,ins2,ins3,numsol,wall;
  double generate_disk_bb();

  iter = 0;
  hits = 0;
  do {   /* loop to ensure that the photon hits the source */

    iter ++;

  /* First, generate position in the outer disk from which the photon 
     is emitted. This is either from the r*r^(-3)/r^(-3/4) distribution, 
     in which case make sure it is not too far away, OR this is from a ring. */

    *energy = generate_disk_bb(&r_emit);
    *weight = 1;

    eindx = num_bin_i*log10((*energy)/emin);

  /* generate direction of motion */

  /* 
     there are TWO  phi angles here:
     1. global position of the point of emission (because of precession, 
                                                    it matters !!)
     2. local azimuthal angle of emission
     
  */

  /* this is the global point of emission */
    phi = two_pi*rnd();
    posi[0] = r_emit*cos(phi);
    posi[1] = r_emit*sin(phi);
    posi[2] = 0;

    /* direction of emission wr to the global coordinate system */
    phi = two_pi*rnd();
    switch (irradiation) {
    case 0 : {             /* semi-isotropic      */
      ct = rnd();
      break; }
    case 1 : {             /* cosine - like   ??? */
      ct = sqrt(rnd());
      break; }
    case 2 : {             /* beamed irradiation  */
      printf("We should not be here (1) \n");
      exit(13); }
    case 3 : {       /* Chandrasekhar-Sobolev distribution (equiv. to sqrt?? */
      /* 1 + 2*mu    */
      ct = 0.5*(sqrt(1+8*rnd())-1);
      break; }
    }

    if( rnd() < 0.5) ct = -ct;
    st = sqrt(1-ct*ct);
    
    momi[0] = st*cos(phi);
    momi[1] = st*sin(phi);
    momi[2] = ct;

    if (precession_option<2) {
      prectransfsimple(posi,momi,pos,mom,prectheta,precphi,1); 
    } else {
      prectransf(posi,momi,pos,mom,prectheta,precphi,1); 
    }

    hits = checkhit(pos,mom,ipos,&t);
    
    if (hits) {
      for(j=0; j<3; j++) 
	ipos[j] += mom[j]*(t*1e-6);

      rq = ipos[0]*ipos[0]+ipos[1]*ipos[1]+ipos[2]*ipos[2];
      if (rq>torus_r2) {
	printf("Dziwne ...  %lf  %lf\n",sqrt(rq),sqrt(torus_r2));
      }
      
      for(j=0; j<3; j++) pos[j] = ipos[j];
    }

    /* store all emitted photons */
    inp_s[eindx][0] += *weight;
    
  } while  ( (!hits) && (iter < 10000 ));
  
  /* here store photons intercepted by the comptonizing cloud */
  inp_s[eindx][1] += *weight;

  if (iter >= 10000) {
    printf("iter: %d \n",iter); 
  }
      
  return 0;
}


checkhit(pos,mom,ipos,t)
double pos[],mom[],ipos[],*t;
{
  int    numfar,numtop,numinn,hit,i,j;
  double ipos2[3][3],tt[3],ti[3],tf;
  
  *t = 1e12;
  hit = 0;
  
  faredge(pos,mom,&numfar,ipos,&tf);
  if (numfar>0) {
    hit = 1;
    *t = tf;
  }
  
  topbottom(pos,mom,&numtop,ipos2,tt);
  if (numtop > 0) {
    for(i=1; i<=numtop; i++) {
      if (tt[i] < (*t)) {
        *t = tt[i];
        hit = 1;
        for(j=0; j<3; j++) ipos[j] = ipos2[i][j];
      }
    }
  }
  
  nearedge(pos,mom,&numinn,ipos2,ti);
  if (numinn > 0) {
    for(i=1; i<=numinn; i++) {
      if (ti[i] < (*t)) {
        *t = ti[i];
        hit = 1;
        for(j=0; j<3; j++) ipos[j] = ipos2[i][j];
      }
    }
  }
  
  /*
    if (hit && (numtop>0) && (numfar==0)) {
    printf("==================================================== \n");
    printf("numfar = %d; tfar = %.3lf \n",numfar,tf);
    printf("numtop = %d; tt = %.3lf  %.3lf \n",numtop,tt[1],tt[2]);
    printf("numinn = %d; ti = %.3lf  %.3lf \n",numinn,ti[1],ti[2]);
    printf("t=%.3lf, ipos: %.2lf %.2lf %.2lf \n",*t,ipos[0],ipos[1],ipos[2]);
  }
  */

  return hit;
}


prectransf(posi,diri,posf,dirf,theta,phi,trans)
double posi[],posf[],diri[],dirf[],theta,phi;
int    trans;
{
  double R[3][3],posm[3],dirm[3];

  if (trans == 1) {
    rotmatrix(R,0.174533,theta);   /* that's 10 degs */
    multmatvec(posm,R,posi);
    multmatvec(dirm,R,diri);
    
    rotmatrix(R,phi,theta);
    multmatvec(posf,R,posm);
    multmatvec(dirf,R,dirm);
  } else {   /* inverese */
    invrotmatrix(R,phi,theta);
    multmatvec(posm,R,posi);
    multmatvec(dirm,R,diri);
    
    invrotmatrix(R,0.174533,theta);
    multmatvec(posf,R,posm);
    multmatvec(dirf,R,dirm);
  }
  return 0;
}

prectransfsimple(posi,diri,posf,dirf,theta,phi,trans)
double posi[],posf[],diri[],dirf[],theta,phi;
int    trans;
{
  double R[3][3];

  if (trans==1) {
    rotmatrix(R,phi,theta);
    multmatvec(posf,R,posi);
    multmatvec(dirf,R,diri);
  } else {
    invrotmatrix(R,phi,theta);
    multmatvec(posf,R,posi);
    multmatvec(dirf,R,diri);
  }

  return 0;
}


rotmatrix(R,phi,theta)    /* Byron-Fuller convention, B&F, p. 15 */
double R[3][3],phi,theta; /* but it works if anty-clockwise revolution */
{                         /* of theta, rather than clockwise     */
  double sf,cf,st,ct;

  sf = sin(phi);
  cf = cos(phi);
  st = sin(theta);
  ct = cos(theta);

  /* row col  */
  R[0][0] = cf*ct;
  R[0][1] = sf*ct;
  R[0][2] = -st;
  R[1][0] = -sf;
  R[1][1] = cf;
  R[1][2] = 0;
  R[2][0] = cf*st;
  R[2][1] = sf*st;
  R[2][2] = ct;

  return 0;
  
}

invrotmatrix(R,phi,theta)    /* inverse matrix (transposed) */
double R[3][3],phi,theta; 
{                         
  double sf,cf,st,ct;

  sf = sin(phi);
  cf = cos(phi);
  st = sin(theta);
  ct = cos(theta);

  /* row col  */
  R[0][0] = cf*ct;
  R[0][1] = -sf; 
  R[0][2] = cf*st;
  R[1][0] = sf*ct;
  R[1][1] = cf;
  R[1][2] = sf*st;
  R[2][0] = -st;
  R[2][1] = 0;
  R[2][2] = ct;

  return 0;
  
}

multmatvec(A,R,B)
double A[3],B[3],R[3][3];  /* A=R*B */
{
  int i,j;

  for(i=0; i<3; i++) {
    A[i] = 0;
    for(j=0; j<3; j++) {
      A[i] += R[i][j]*B[j];
    }
  }

  return 0;

}
/* ======  end of torus section =============================================*/

generate_disk_with_sphere(pos,mom,energy,weight)
     double pos[],mom[],*energy,*weight;
     /*
Generates seed photon from a disk external to a spherical source. 
Disk temperature varies as r^(-0.75).
      */
{
  double r_emit,rratio,
    ct,st,phi,psi,psicrit,rq,rp,h2,ll,phi_crit,rad;
  int    i,tindx,eindx,iter,hits;
  double generate_disk_bb();

  iter = 0;
  do {   /* loop to ensure that the photons hits the source */

    iter ++;

  /* First, generate position in the outer disk from which the photon 
     is emitted. This is either from the r*r^(-3)/r^(-3/4) distribution, 
     in which case make sure it is not too far away, OR this is from a ring. */

    *energy = generate_disk_bb(&r_emit);

    eindx = num_bin_i*log10((*energy)/emin);

    *weight = 1;

  /* generate direction of motion */

    rratio = Rdisk/r_emit;

    /* Chandrasekhar-Sobolev angular distribution: 1+2*mu, mu from normal  */
    if (rratio>1) {                   /* inside the sphere */
      switch (irradiation) {
      case 0 : {             /* semi-isotropic      */
	ct = rnd();
	break; }
      case 1 : {             /* cosine - like   ??? */
	ct = sqrt(rnd());
	break; }
      case 2 : {             /* beamed irradiation  */
	printf("We should not be here (1) \n");
	exit(13); }
      case 3 : {   /* Chandrasekhar-Sobolev distribution (equiv. to sqrt?? */
	ct = 0.5*(sqrt(1+8*rnd())-1);
	break; }
      }

      st = sqrt(1-ct*ct);
      phi = two_pi*rnd();

      mom[0] = st*cos(phi);      /* cartesian components of momentum */
      mom[1] = st*sin(phi);
      mom[2] = ct;
      
      pos[0] = r_emit;
      pos[1] = 0;
      pos[2] = 0;
      tindx = ct*nangle+2;
      inp_s[eindx][    0] += *weight;
      /*      inp_s[eindx][tindx] += *weight;  */

      hits = 1;
    } else {

      phi = two_pi*rnd();
      switch (irradiation) {
      case 0 : {             /* semi-isotropic      */
	ct = rnd();
	break; }
      case 1 : {             /* cosine - like   ??? */
	ct = sqrt(rnd());
	break; }
      case 2 : {             /* beamed irradiation  */
	printf("We should not be here (1) \n");
	exit(13); }
      case 3 : {       /* Chandrasekhar-Sobol distribution (equiv. to sqrt?? */
	ct = 0.5*(sqrt(1+8*rnd())-1);
	break; }
      }
      st = sqrt(1-ct*ct);

      if (iter==10000) {
	ct = 0;
	st = -pi;
	*weight = 0;
	printf("iter 10000 \n");
      }
      mom[0] = st*cos(phi);
      mom[1] = st*sin(phi);
      mom[2] = ct;

      pos[0] = r_emit;
      pos[1] = 0;
      pos[2] = 0;

      rq = r_emit*r_emit;
      psicrit = -sqrt(rq-Rq)/r_emit;
      psi = mom[0];
      
      tindx = ct*nangle+2;
      inp_s[eindx][    0] += *weight;

      rp = r_emit*psi;
      h2 = sqrt(rq-rp*rp);
      hits = (psi<psicrit) && (h2 < Rout);

      if (!hits)
	inp_s[eindx][tindx] += *weight;

    }
  } while  ( (!hits) && (iter < 10000 ));
  

  /* here store photons intercepted by the comptonizing cloud */
  inp_s[eindx][1] += *weight;

  /*  printf("iter: %ld \n",iter); */
      
  if (rratio <= 1) {
    ll = -rp - sqrt(Rq-h2*h2);
    if (ll < 0) {
      printf("ll<0 ! %.3f  %.3f\n",ll,rp);
    }
    total_dist = ll;
    for (i=0; i<3; i++) 
      pos[i] = pos[i] + ll*mom[i];
    
    rad = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    if (fabs(rad/Rout-1) > 1e-3) {
      printf("Wrong ! %lf \n",rad);
      exit(13);
    }
  }

  return 0;
}

initial_position(mom,pos)
  double mom[],pos[];
/*
   generates initial position of a photon, according to assumed
   spatial distribution of sources
   returns Cartesian coordinates
*/
{
  int     i,ir=1;
  double  x,y,z,r2,ct,phi;

  switch (spat_distr) {
  case 0 : {                  /* point-like central source */
    for (i=0; i<3; i++)
      pos[i] = 1e-9*Rout + Rin*mom[i];
    break;
  }
  case 1 : {                   /* uniform distribution      */
    do  {
      x = rnd();
      y = rnd();
      z = rnd();
      r2 = x*x+y*y+z*z;
    } while ((r2>Rq) || (r2<Rin2));
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    break;
  }
  case 2 : {                   /* slab                  */
    pos[0] = 0;
    pos[1] = 0;
    pos[2] = 0;
    break;
  }
  case 4 : { /* Sunyaev-Titarchuk strange distribution */
    x = rnd();
    c_hunt(ST_distr,nop_ST,x,&ir);
    r2 = ST_rad[ir] + (x-ST_distr[ir])*(ST_rad  [ir+1]-ST_rad  [ir])/
                                       (ST_distr[ir+1]-ST_distr[ir]);
    ct = rnd();
    phi= rnd()*two_pi;
    pos[0] = r2*ct*cos(phi);
    pos[1] = r2*ct*sin(phi);
    pos[2] = r2*sqrt(1-ct*ct);
    break ;
  }
  case 5 : { /* active region - hemisphere */
    /* all done earlier */
    break ; 
  }
  }
  return 0;
}

double distance(pos,mom,nempty,te)
     double pos[],mom[],te[];
     int    *nempty;
/*
    calculates the distance from point 'pos' to the boundary of 
    the region along the direction 'mom' of the photon's motion. (9) in G&W.

    dempty is the distance travelled through the inner void, in the case
    of shell; 'distance' is in this case, the distance travelled through 
    the plasma
*/
{
  double r,rq,rp,psi,d,h2,psicrit,rsq,rtv,vsq,t,intpos[3],dempty;
  int   i,ib;

    *nempty = 0;
    for(i=0; i<5; i++) te[i] = 0;

    switch (spat_distr) {              /* sphere */
    case 0:
    case 1:
    case 3:
    case 4:
      {
	rq = pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];
	r = sqrt(rq);
	psicrit = -sqrt(rq-Rin2)/r;
	psi = (pos[0]*mom[0]+pos[1]*mom[1]+pos[2]*mom[2])/r;
	rp = r*psi;
	h2 = rq-rp*rp;
	d = sqrt(Rq-h2)-rp;
	/* the condition involving psicrit appears because the direction
	   of motion is important !!! */
	if ( (psi<psicrit) && (h2 < Rin2)) {
	  dempty = 2*sqrt(Rin2-h2);
	  d -= dempty;
	}
	if (d < 0.) {
	  printf("d<0 ! %.3f d= %.3f r= %.3f  rp=%.3f\n",dempty,d,r,rp);
	}
	break;
      }
    case 2 :   /* slab */
      {
	psi = mom[2];
	if (psi>0)
	  d = (Rout-pos[2])/psi;
	else
	  d = fabs(pos[2]/psi);
	
	/*	 if (d<0) {
		 printf("dist<0!, pos: %.4f,  mom: %.4f
		 }
	*/
	break;
      }
    case 5 :  /* externally illuminated source */
      {
	rsq = pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];
	rtv = mom[0]*pos[0]+mom[1]*pos[1]+mom[2]*pos[2];
	vsq = mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2];

	d = (sqrt(rtv*rtv-vsq*(rsq-Rq))-rtv)/vsq;

	if (mom[2]<0) 
	  d = d<(-pos[2]/mom[2])?d:(-pos[2]/mom[2]);
	
	break;
      }
    case 6 :  /* externally illuminated torus */
      {

	whichwall(pos,mom,intpos,&ib,&t,nempty,te);

	d = t;
	for(i=1; i<=(*nempty); i++) {
	  d -= (te[2*i]-te[2*i-1]);
	}
	  
	break;
      }


    }

    return d;
}

#include "mtc_incl_wedge.c"

double X(x, xprim)
   double x,xprim;
/*            Formula (2) from G&W                              */
{
   double x1,x2,xx;

   x1 = 1/x;
   x2 = 1/xprim;
   xx = x*x2+xprim*x1+4*(x1-x2)*(1+x1-x2);

   return xx;
}

int jl_fff = 1;
double Invers_F(Fx)
   double Fx;
{
  double xx;

  c_hunt(invsigint,MAXFF1,Fx,&jl_fff);
  xx = invsigval[jl_fff] + (Fx-invsigint[jl_fff])*
    (invsigval[jl_fff+1]-invsigval[jl_fff])/
    (invsigint[jl_fff+1]-invsigint[jl_fff]);
  
  return xx;

}

/*
double Invers_F(Fx)
   double Fx;
{
   double  Ft  = 2.885259743,
           Fx1 = 5.881379e-4,
           Fx2 = 2.366133e-2,
           y,t;

   if (Fx<=Fx1)
      y = sqrt(1.5*Fx);
   else
      if (Fx<=Fx2)
      {
         y = 0.333333*(acos(1-9*Fx)+2*two_pi);
         y = cos(y);
         y += 0.5;
      }
      else
      {
         t = log(Fx) + Ft;
         y = -0.5+t*(0.228 + t*(4.708e-3 + t*(3.114e-4 + t*(-1.921e-5 + t*2.873e-7))));
         y = pow(10.,y);
      }
   return y;
}
*/

double sigma_integr_(v)
   double v;
/*
    integrand in expression for <sigma>
    global variable : energ2
*/
 {
   double g,f,x1,x2,sigm,lor;

   lor = gamma(v);
   g = energ2*lor;
   f = g*v;
   x1 = g-f;
   x2 = g+f;
   sigm = (F(x2)-F(x1))*Relat_v(v)/(g*g*v);

   return 0.375*sigm;
}

double Relat_v(v)

   double v;
/* relativistic distribution in velocities */
{
   double v0,v1,v2,v3,p0;
   double  pf;
   
   if (med_temp < 511.) {
     v0 = v;
     v1 = 1/((1-v0)*(1+v0));
     v2 = v0*v1;
     v3 = sqrt(v1);
     p0 = v2*v2*v3;
     p0 *= exp(-med_temp*v3);
     pf = p0*norm_Rev;   /* global variable, needs to be calculated first */
   } else {
     v0 = v;
     p0 = v0*v0*exp(-med_temp*v0*v0);
     pf = p0*norm_Rev;
   }
    
   return pf;
}

double ro1prim(v)
   double v; 
/*
   distribution ro1' (formula (A9)) from G&W
*/
{
   double t1,t2,x1,x2;

   t1 = Relat_v(v);
   t1 *= 1/v-v;
   x1x2_fun(v,&x1,&x2);
   t2 = F(x2)-F(x1);
   t1 *= t2;
   t1 *= norm_ro1;     /* global variable, needs to be calculated first */
   return t1;
}


int jl_xxx=1;
double F(x)
  double x; 
{
  double xx;

  c_hunt(invsigval,MAXFF1,x,&jl_xxx);
  xx = invsigint[jl_xxx] + (x-invsigval[jl_xxx])*
    (invsigint[jl_xxx+1]-invsigint[jl_xxx])/
    (invsigval[jl_xxx+1]-invsigval[jl_xxx]);

  return xx;
}

double F_old();

double F_old(x)
  double x; 
/*  Analitical approximation to integral of y*sigma(y)  (PSS || GW) */
{
  double x1,y;

  x1 = x+1;
  if (x < 0.5) {
    y = x*x*(0.166666667+0.5/x1+x*(0.047-x*0.03));
  } else {
    if (x < 3.5) {
      y = x1*log(x1)-0.94*x-0.00925;
    } else {
      y = x1*log(x1)-0.5*x-13.16*log(2+0.076*x)+9.214;
    }
  }

  return y;
  /*
  if (x > 3.5) {
    y = x1*log(x1)-0.5*x-13.16*log(2+0.076*x)+9.214;
  } else {
    if (x > 0.5) {
      y = x1*log(x1)-0.94*x-0.00925;
    } else {
      y = x*x*(0.166666667+0.5/x1+x*(0.047-x*0.03));
    }
  }
  return y;
  */
}

x1x2_fun(v,x1,x2)
   double v,*x1,*x2;
/*
   functions x1(v), x2(v2) defined after eq. (A4) in G&W
*/
{
    double g,f;

    g = energ2*gamma(v);
    f = g*v;
    *x1 = g-f;
    *x2 = g+f;

    return 0;
}

double Et0;

double diskbb(en)
     double en;
{
  double dsbbint(),x2,dd,en2;
  
  Et0 = en/rad_temp;
  x2 = log(30./Et0);
  if (x2 < 1) x2 = 1;
  dd = qsimp(dsbbint,0.,x2);
  en2 = en*en;

  return norm_spc*dd*en2*en2;
}

double dsbbint(lx)
     double lx;
{
  double fac,db,x;

  x = exp(lx);
  fac = Et0*exp(0.75*lx);
  if (fac < 1e-6) {
    db = 1./fac;
  } else {
    if (fac < 70) {
      db = 1./(exp(fac)-1);
    } else {
      db = 0;
    }
  }

  db *= x*x;

  return db;
}

double power_law(en)
   double en;
/*  
    power-law incident spectrum - energy per log interval
    global variables: alpha, norm_spc
*/
{
   double x;

   x = pow(en,-alpha+1);
   x *= norm_spc;

   return x;
}

double Planck_photons(T0)
     double T0;
/*
  generates photon's energies according to Planck distribution

*/
{

  double x1,x2,x3,x4,x1x,alpha,x,sj;
  int    j;

  x1 = rnd();
  x2 = rnd();
  x3 = rnd();
  x4 = rnd();
 
  x1x = x1*1.202;

  if (x1x < 1.) {
    alpha = 1.;
  } else {
    sj = 0.;
    j = 0;
    do {
      j ++;
      sj += 1./(j*j*j);
    } while (sj < x1x);
    alpha = j;
  }
  
  x = -T0/alpha*log(x2*x3*x4);

  return x;
}


double Planck(en)
   double en;
/*
    Planck incident spectrum  - energy per log interval
    global variables: rad_temp, norm_spc
*/
{
   double eq,x;

   eq = en*en;
   eq = eq*eq;
   x = eq/(exp(en/rad_temp)-1);
   x *= norm_spc;

   return x;
}

double rnd()
{
  double r;
  r = ran2(&idum);
  return r;
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 2.0e-17
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software . */

find(xa,x,n)
  double xa[],x;
  int    n;
{
  int  k,klo,khi;

  klo = 0;
  khi = n-1;
  while (khi-klo > 1) 
   {
      k = (khi+klo) / 2;
      if (xa[k] > x) 
         khi = k;
      else
         klo = k;
   }
   return klo;
}

c_hunt(double xx[], int n, double x, int *jlo)
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


double gamma(v)   /*   Lorentz factor  */
  double v;
{
  return 1/sqrt((1-v)*(1+v));
}

void derivs(x,y,dydx)
   double x,y[],dydx[];
{
   dydx[0] = ro1prim(x);
}

output_results(runnum,noph)
  long int   noph;
  int        runnum;
{
  FILE   *fp;
  int     i,nos,j,k; 
  double  e,st,si,sr,er,d,de;

  char   out_file[20],
         inp_spec[20],
         distfile[20],
         scatfile[20],
         photfile[20],
         scatspec[20],
         reflspec[20];
  char    num[10];

  if (noph == 0) return 0;

  strcpy(out_file,"mcomp000000.dat");
  strcpy(inp_spec,"inpsp000000.dat");
  strcpy(distfile,"lcurv000000.dat");
  strcpy(scatfile,"nscat000000.dat");
  strcpy(photfile,"photr000000.dat");
  strcpy(scatspec,"sctsp000000.dat");
  strcpy(reflspec,"refls000000.dat");


  sprintf(num,"%06i",runnum);
  memcpy(out_file+5,num,6);
  memcpy(photfile+5,num,6);
  memcpy(inp_spec+5,num,6);
  memcpy(distfile+5,num,6);
  memcpy(scatfile+5,num,6);
  memcpy(scatspec+5,num,6);
  memcpy(reflspec+5,num,6);
  
  fp = fopen(inp_spec,"w");         /* write input spectrum */
  nos = num_bin_i*log10(emax/emin);
  for (i=0; i<nos; i++) {
    e = emin*pow(10.,(i+0.5)/num_bin_i);
    de = emin*(pow(10.,(i+1.0)/num_bin_i)-pow(10.,(1.0*i)/num_bin_i));
    er = e/de/noph;
    si = inp_s[i][0]*er;
    sr = inp_s[i][1]*er;
    total_input += si*de;
    fprintf(fp,"%e %e %e ",e*511,si,sr);
    for (j=2; j<=nangle+1; j++) {
      fprintf(fp,"%.3e ",inp_s[i][j]*er*nangle);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  fp = fopen(out_file,"w");
  nos = num_bin_o*log10(emax_o/emin_o);
  for (i=0; i<nos; i++)    {
    e = emin_o*pow(10.,(i+0.5)/num_bin_o );
    de = emin_o*(pow(10.,(i+1.0)/num_bin_o)-pow(10.,(1.0*i)/num_bin_o));
    er = e/de/noph;
    st = out_s[i]*er;
    sr = out_r[i]*er;
    total_output_up += st*de;
    total_output_down += sr*de;
    fprintf(fp,"%e %e %e ",e*511,st,sr);
    if ((spat_distr == 2) || (spat_distr == 3) || 
	(spat_distr == 5) || (spat_distr == 6)) {
      for(j=0; j<nangle; j++) {
	for(k=0; k<MAXphi; k++) {
	  fprintf(fp,"%.3e ",out_si[i][j][k]*er*nangle);
	}
	fprintf(fp,"  ");
      }
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  /* cold, reflected  spectra */


  if (do_reflection) {
    fp = fopen(reflspec,"w");
    nos = num_bin_o*log10(emax_o/emin_o);
    for (i=0; i<nos; i++)    {
      e = emin_o*pow(10.,(i+0.5)/num_bin_o );
      de = emin_o*(pow(10.,(i+1.0)/num_bin_o)-pow(10.,(1.0*i)/num_bin_o));
      er = e/de/noph;
      total_refl += reflect[i][0][0]*er*de;
      fprintf(fp,"%e %.3e %.3e ",e*511,reflect[i][0][0]*er,incid[i]*er);
      for(j=1; j<=nangle; j++) {
	for (k=1; k<=MAXphi; k++) {
	  fprintf(fp,"%.3e ",reflect[i][j][k]*er*nangle);
	}
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
  
  /* spectra after given number of scatterings */
  
  if (max_spsc > 0) {
  fp = fopen(scatspec,"w");       
  nos = num_bin_o*log10(emax_o/emin_o);
  for (i=0; i<nos; i++) {
    e = emin_o*pow(10.,(i+0.5)/num_bin_o);
    de = emin_o*(pow(10.,(i+1.0)/num_bin_o)-pow(10.,(1.0*i)/num_bin_o));
    er = e/de/noph;
    fprintf(fp,"%.2e ",e*511);
    for (j=0; j<max_spsc; j++) {
      si = scatt_spec[i][j]*er;
      fprintf(fp,"%.2e ",si);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  }

  if (nop_grid) {
    fp = fopen(photfile,"w");
    for (i=0; i<nop_grid; i++)    {
      st = photar[i]/noph;
      fprintf(fp,"%e %e \n",ear[i]*511,st);
    }
    fprintf(fp,"%e \n",ear[nop_grid]*511);
    fclose(fp);
  }

  if (num_dist_en > 0) {
  fp = fopen(distfile,"w");
  for (i=0; i<dist_points; i++) {
    d = total_dist_max/dist_points*(i+0.5);
    fprintf(fp,"%.3e ",d);
    for (j=0; j<num_dist_en; j++) {
      fprintf(fp,"%.3e ",dist_hist[j][i]/noph);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  }

  /* scatterings histogram */

  if (num_dist_en > 0) {
  fp = fopen(scatfile,"w");
  for (i=0; i<nscatt_glob; i++) {
    fprintf(fp,"%5i ",i);
    for (j=0; j<num_dist_en; j++) {
      fprintf(fp,"%.3e ",scatt_hist[j][i]/noph);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  }

  return 0;
}


input_data(runnum)
     int runnum;
/*
   Prepares initial data: parameters of the medium, etc.
*/
{
   FILE    *fp,*fc;
   char     fcum_name[80],
            spec_file_name[80],buff[100],en_grid_file[80];
   double   energy,e1,e2,x,y,
            v1 = 1e-6,v2 = 1 - 1e-6;
   int      readc,i,j,k,n1,n2;

   double   tor_dhdr,ctgt;
   char     inp_file[20];
   char     num[10];

   strcpy(inp_file,"mtc__inp.000000");
   sprintf(num,"%06i",runnum);
   memcpy(inp_file+9,num,6);
   

   if ((fp = fopen(inp_file,"r")) == NULL)
   {
      printf("File with input data not found - %s \n",inp_file);
      exit(13);
   }
   printf("\nReading input data from file: %s \n",inp_file);
      
   fscanf(fp,"%s",buff);
/*  'energy' for which the cumul. distr. of ro1' will be calculated  */
   fscanf(fp,"%lf %lf ",&med_temp,&energy);
   fscanf(fp,"%s",buff);
   fscanf(fp,"%d %lf %lf %lf %lf %d %s",&spectrum,&rad_temp,&alpha,
                              &emin,&emax,&nop_spec,spec_file_name);
   fscanf(fp,"%s",buff);
   fscanf(fp,"%lf %lf %d %d %d %s",&emin_o,&emax_o,&num_bin_i,&num_bin_o,
	  &nop_grid,en_grid_file);
   /*   printf("%.2e %.2e %ld %ld \n",emin_o,emax_o,num_bin_i,num_bin_o); */
   fscanf(fp,"%s",buff);
   fscanf(fp,"%ld",&mc_steps);
   fscanf(fp,"%s",buff);
   fscanf(fp,"%d %lf %d %lf %d",
	  &spat_distr,&r_emission,&irradiation,&theta_irrad,&nangle);
   fscanf(fp,"%s",buff);
   fscanf(fp,"%lf %lf %lf %lf %d",
	  &tor_dhdr,&torus_theta0,&prectheta,&precphi,&precession_option);
   fscanf(fp,"%s",buff);
   fscanf(fp,"%lf %lf %lf %lf",&tau_max,&Rin,&dens_alpha,&rho0);
   printf("spat_distr: %d,  tau: %.2f, Rin: %.2f,   number of photons: %ld\n",
	  spat_distr,tau_max,Rin,mc_steps);
   fscanf(fp,"%s",buff);
   fscanf(fp,"%d %lf",&do_reflection,&b_bulk); 
   fscanf(fp,"%s",buff);
   fscanf(fp,"%lf %d",&total_dist_max,&dist_points);
   fscanf(fp,"%s",buff);
   fscanf(fp,"%ld %li %li",&idum,&xo,&b);
   fscanf(fp,"%s",buff);
   fscanf(fp,"%d",&num_dist_en);

   if (num_dist_en > 0) {
   for(i=0; i<num_dist_en; i++) {
     fscanf(fp,"%lf %lf",&e1,&e2);
     dist_energies[0][i] = e1;
     dist_energies[1][i] = e2;
   }
   }
   fscanf(fp,"%s",buff);
   fscanf(fp,"%d",&max_spsc);
   if (max_spsc > 0) {
   for(i=0; i<max_spsc; i++) {
     fscanf(fp,"%d %d",&n1,&n2);
     spsc[0][i] = n1;
     spsc[1][i] = n2;
   }
   }

/* Bei: set-up of the accretion disk */
   fscanf(fp,"%s",buff);
   fscanf(fp,"%lf %lf %lf",&bh_mass, &bh_spin, &bh_mdot);
/* Bei: set-up of the accretion disk */

   fclose(fp);

   med_temp /= 511;       // temperature of the torus: 100 keV here
   energy /= 511.;        // initial photon energy ~ 1 keV
   rad_temp /= 511.;
   T_bb = rad_temp;       // sed photon temperature ~ 0.1 keV
   emin /= 511.;
   emax /= 511.;
   emin_o /= 511.;
   emax_o /= 511.;


   if (spat_distr == 6) {     /* wedge stuff */
    torus_theta0 *= pi/180.;                 // half-torus opening angle: 15 degree
    ctgt = 1./tan(torus_theta0);         
    wedge_ctgt2 = ctgt*ctgt;
    wedge_sintheta = sin(torus_theta0);
    wedge_costheta = cos(torus_theta0);
    Rout = rho0;                             // outer radius of the torus: 30Rg

    tTh1Rg = tau_max/(Rout-Rin);
    torus_r1 = Rin*Rin;                      // inner radius of the torus
    torus_r2 = Rout*Rout;                    // outer radius of the torus
    Rdisk = 1.0001*Rout;                     // inner radius of the disk
    prectheta *= pi/180; 
    precphi  *= pi/180;
  } 

   return 0;

}

prepare_ST_distr()
{
  double dr,r,ST_fun(),s;
  int    i;
  FILE *fp;

  if (Rin>0.) {
    printf("Sunyaev-Titarchuk spatial distribution AND Rin>0 ?? \n");
    Rin = 0;    
  }

  nop_ST = 100;
  dr = tau_max/nop_ST;
  s = 0;
  ST_distr[0] = 0;
  ST_rad[0] = 0;
  for(i=1; i<=nop_ST; i++) {
    r = dr*(i-0.5);
    s += ST_fun(r)*r*r;
    ST_rad[i] = r;
    ST_distr[i] = s;
  }
  for(i=1;i<=nop_ST; i++) {
    ST_distr[i] /= ST_distr[nop_ST];
  }

  /*  fp = fopen("ST_distr.dat","w");
  for (i=0; i<=nop_ST; i++) {
    fprintf(fp,"%e %e \n",ST_rad[i],ST_distr[i]);
  }
  fclose(fp); */
  return 0;
}

double ST_fun(tau)
     double tau;
{
  double t,st;
  
  if (tau < 1e-5) {
    st = 1;
  } else {
    t = tau*pi/tau_max;
    st = sin(t)/t;
  }

  return st;
}

prepare_cumulative_distr()
/*
   Integrates distribution to obtain the cumulative distribution for ro1'.
   Calculates global variables: norm_Rev  norm_ro1

   Uses global variables: med_temp
   sets global variable   energ2
*/
{ 
   double   v1 = 1e-6,
            v2 = 1 - 1e-6,
            acc = 1e-12,
            hfirst = (v2-v1)/1E4,
            ystart[MAXV];
   int      i,nok,nbad,nop,j;
   FILE      *fp;

   if (emax_o < 1.95) {
     max_roe = MAXRE-3;
   } else { 
     max_roe = MAXRE-1;
   }

   for (j=0; j<=max_roe; j++) {
/*     printf("Calculating cumulative distribution for ro1' for e= %.3f ... ",
	    rho1_ens[j]);*/
     rho1_ens[j] /= 511.;
     energ2 = 2*rho1_ens[j];
   
     norm_ro1 = 1;
     norm_ro1 = 1/qsimp(ro1prim,v1,v2);  /* of the two distributions         */
     ystart[0] = 0;
     odeint(1,derivs,ystart,1,v1,v2,acc,hfirst,1E-30,&nok,&nbad,&nop);
     nop_rho1[j] = nop;
     for (i=0; i<nop; i++) {
       y_cum[j][i+1] = xp[i];       /* inversion of cumulative distribution */
       x_cum[j][i+1] = yp[0][i];
     }
     /*     printf("finished. NoP: %i \n",nop); */
     
     /*     fp = fopen("cumul.dat","w");  
     for(i=0; i<nop_rho1[j]; i++) {
       fprintf(fp,"%e  %e\n",x_cum[j][i],y_cum[j][i]);  
     }
     fclose(fp);  
     getchar(); */
   }

   return 0;
}

prepare_spectrum(spec_file_name)
   char spec_file_name[];
{
    double e,s,sum,qtrap();
    int    n,i;
    FILE   *fp;

    switch (spectrum)
    {
        case 0 :
        {
           norm_spc = 1/qsimp(Planck,emin,emax);
/*
           fp=fopen("test.dat","w");
           n = 100;
           for (i=0; i<n; i++)
           {
              e = exp(log(emin)+i*log(emax/emin)/n);
              s = Planck(e);
              fprintf(fp,"%e %e %e\n",e,s,s/e);
           }
           fclose(fp);
*/
           break;
        }
        case 1 :
        {
           norm_spc = 1/qsimp(power_law,emin,emax);
           break;
        }
        case 2 :  /* read from disk, normalize, etc.     */
        {
           if ((fp=fopen(spec_file_name,"r")) != NULL)
           {
              sum = 0;
              emin = 1e38;
              emax = -1e38;
              for (i=0; (i<MAXP) && (! feof(fp)); i++)
              {
                 fscanf(fp,"%lf %lf",&e,&s);
                 if (e<emin) emin = e; 
                 if (e>emax) emax = e;
                 *(spect+i) = s/e;
                 *(energ+i) = e;
                 sum += *(spect+i);
              }
              fclose(fp);
              nop_spec = i-1;
              sum *= log(emax/emin)/nop_spec;
              for (i=0; i<nop_spec; i++)
                 *(spect+i) /= sum;
           }
           else
           {
              printf("Can not open file %s \n",spec_file_name);
              exit(13);
           }
           break;
        }
        case 3 : {   /* disk black body -> tabulate the spectrum */
	  norm_spc = 1/qtrap(diskbb,emin,emax);
	  nop_spec = num_bin_i*log10(emax/emin);
	  /*	  fp=fopen("test.dat","w"); */
	  for(i=0; i<nop_spec; i++) {
	    e = exp(log(emin)+i*log(emax/emin)/(nop_spec-1.));
	    energ[i] = e;
	    spect[i] = diskbb(e)/e;
	    /*	    fprintf(fp,"%e %e \n",energ[i]*511,spect[i]); */
	  }
	  /*	  fclose(fp); */
	  spectrum = 2; 
	  break;
	}
        case 4 : {
	  break ;
        }
        default:
        {
	  printf("Unknown type of spectrum : %i \n",spectrum);
          exit(13);
        }
    }
    elog = log(emax/emin);

    return 0;
}


prepare_sigma()
{
  double e,e1,e2,v1=1e-6,v2=0.9999999;
  int    i;
  FILE   *fp;

  /*  fp=fopen("sigma.dat","w");  */
  nop_se = 100;
  e1 = 0.5*emin_o;
  e2 = 2*emax_o;
  for(i=1; i<=nop_se; i++) {
    e = exp(log(e1)+(i-1.)*log(e2/e1)/(nop_se-1.));
    energ2 = 2*e;
    sig_e[i] = e;
    sigmm[i] = qsimp(sigma_integr_,v1,v2);
    /*    fprintf(fp,"%e %e \n",e*511,sigmm[i]);   */
  }
  /*  fclose(fp);
      exit(13);  */
  
  return 0;
}


double stheta,ctheta;

prepare_sigma_bulk()
{
  double e1,e2,ct1,ct2,e,ct,ss,
         v1 = 0.000001, v2 = 0.999999;
  double quad3d(double (*func)(double, double, double), double x1, double x2);
  double func(double vel, double alpha, double phi);

  int    j,i;

  FILE   *fp;


  /*  fp = fopen("sigma.dat","w"); */

  nop_se = 50;
  e1 = 0.1/511.;
  e2 = 1.1*emax_o;

  nop_ct = 21;
  ct1 = -1;
  ct2 =  1;
  for(i=1; i<=nop_se; i++) {
    e = exp(log(e1)+(i-1.)*log(e2/e1)/(nop_se-1.));
    energ2 = 2*e;
    sig_e[i] = e;
    /*     fprintf(fp,"%lf ",e*511); */
    for(j=1; j<=nop_ct; j++) {
      if (i == 1) {
	ct = ct1 + (ct2-ct1)/(nop_ct-1.)*(j-1.);
	sig_ct[j] = ct;
      } 
      ctheta = sig_ct[j];
      stheta = sqrt(1-ctheta*ctheta);
      ss = quad3d(func,v1,v2);
      sigma_mean[i][j] = ss;
      /*      fprintf(fp,"%lf ",ss); */
    }
    /*    fprintf(fp,"\n"); */
  }

  /*  fclose(fp); */

  return 0;
}

double func(double vel, double alpha, double phi)
{
  double s1,
    xKN(double x), xxx(double vel,double alpha, double phi);
  
  s1 = 0.375/two_pi*sin(alpha)*Relat_v(vel)/(energ2*gamma(vel))*
    xKN(xxx(vel,alpha,phi));
  
  return s1;
}


double xKN(double x)
{
  double x1,f;

  x1 = x+1;
  if (x < 0.0001) {
    f = x*(1./3.-11./6.*x)+0.5*(1-1/x1/x1);
  } else {
    f = (1-4/x-8/x/x)*log(x1) + 8/x+0.5*(1-1/x1/x1);
  }

  return f;
}


double  xxx(double vel, double alpha, double phi)
{
  double cross(double vel, double alpha, double phi);
  double xx;
   
  xx = energ2*gamma(vel)*cross(vel,alpha,phi);

  return xx;
}

double cross(double vel, double alpha, double phi)
{
  double ca,c1,c2;
  
  ca = cos(alpha);
  c1 = vel*sin(alpha)*cos(phi)*stheta/gamma(vel);
  c2 = ctheta*(vel*ca+b_bulk);

  return 1 - (c1+c2)/(1+b_bulk*vel*ca);
}


double yy1(double x)
{
  return 0;
}

double yy2(double x)
{
  return pi;
}


double z1(double x, double y)
{
  return 0;
}

double z2(double x, double y)
{
  return two_pi;
}

/*
double qgaus(double (*func)(double), double a, double b)
{
  int j;
  double xr,xm,dx,s;
  static double x[]={0.0,0.1488743389,0.4333953941,
		    0.6794095682,0.8650633666,0.9739065285};
  static double w[]={0.0,0.2955242247,0.2692667193,
		    0.2190863625,0.1494513491,0.0666713443};
  
  xm=0.5*(b+a);
  xr=0.5*(b-a);
  s=0;
  for (j=1;j<=5;j++) {
    dx=xr*x[j];
    s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
  }
  return s *= xr;
}
*/
/* (C) Copr. 1986-92 Numerical Recipes Software . */


static double xsav,ysav;
static double (*nrfunc)(double,double,double);

double quad3d(double (*func)(double, double, double), double x1, double x2)
{
  double qgaus(double (*func)(double), double a, double b);
  double f1(double x);
  
  nrfunc=func;
  return qgaus(f1,x1,x2);
}

double f1(double x)
{
  double qgaus(double (*func)(double), double a, double b);
  double f2(double y);
  double yy1(double),yy2(double);
  
  xsav=x;
  return qgaus(f2,yy1(x),yy2(x));
}

double f2(double y)
{
  double qgaus(double (*func)(double), double a, double b);
  double f3(double z);
  double z1(double,double),z2(double,double);
  
  ysav=y;
  return qgaus(f3,z1(xsav,y),z2(xsav,y));
}

double f3(double z)
{
  return (*nrfunc)(xsav,ysav,z);
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */

