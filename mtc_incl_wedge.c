whichwall(curpos,curdir,intpos,ib,t,nempty,te)
double curpos[],curdir[],intpos[],te[];
int    *ib;  /* ib: 1: inner, 2: outer, 3: top */
double *t;
int    *nempty;
{
  int numtop,numfar,numinn;
  int i,inside();
  double xyz2top[3][3],xyz2inn[3][3],xyzfar[3],ti[3],tt[3],to;
  double tf,x,y,z,dempty;

  FILE *fp;
      
  *nempty = 0;
  *t = 0;
  *ib = 0;
  te[1] = 0;
  te[2] = 0;
  te[3] = 0;
  te[4] = 0;
  /*     try to guess the wall which will be intersected */
  topbottom(curpos,curdir,&numtop,xyz2top,tt);

  /*     the answer is anyway ambigiuos, so try inner/outeredge */
  nearedge(curpos,curdir,&numinn,xyz2inn,ti);
   faredge(curpos,curdir,&numfar,xyzfar,&to);

#ifdef DEBUG
   if (debug_cond) {
     printf("which wall: \n");
     printf("  top:  ns %d tt: %.3lf %.3lf \n",numtop,tt[1],tt[2]);
     printf("  near: ns %d tt: %.3lf %.3lf \n",numinn,ti[1],ti[2]);
     printf("  far:  ns %d tt: %.3lf \n",numfar,to);
     fp = fopen("path.dat","w");
     for(i=0; i<=200; i++) {
       tf = to*i/200.;
       x = curpos[0]+tf*curdir[0];
       y = curpos[1]+tf*curdir[1];
       z = curpos[2]+tf*curdir[2];
       fprintf(fp,"%lf %lf %lf %lf %lf %lf\n",tf,x,y,z,sqrt(x*x+y*y),sqrt(x*x+y*y)/sqrt(x*x+y*y+z*z)-wedge_costheta);
     }
     fclose(fp);
     printf("  path saved; enter a number ");
     scanf("%lf",&x);
     printf("----\n");
   }
#endif

  if (numtop == 0) {
    if (numinn == 0) {
      if (numfar == 1) {
	*t = to;
	for(i=0; i<3; i++) { intpos[i] = xyzfar[i];}
	*ib = 2;
      } else {
	printf("numfar != 1  = %d  after numtop=numinn=0 ; nr=%d\n",numfar,debug);
      }
    }
    if (numinn == 1) {
      *t = ti[1];
      for(i=0; i<3; i++) { intpos[i] = xyz2inn[1][i];}
      *ib = 1;
    }
    if (numinn == 2) {
      if (numfar == 1) {
	*nempty = 1;
	te[1] = ti[1];
	te[2] = ti[2];
	if (te[2]-te[1] < 0) {
	  printf("1 %.3lf  %.3lf \n",te[1],te[2]);
	}
	*t = to;
	for(i=0; i<3; i++) { intpos[i] = xyzfar[i];}
	*ib = 2;
      } else {
	printf("numfar != 1  = %d  after numtop=0 numinn=2 \n",numfar);
      }
    }
  }

  if (numtop == 1) {
    if (numinn == 0) {
      *t = tt[1];
      for(i=0; i<3; i++) { intpos[i] = xyz2top[1][i];}
      *ib = 3;
    }
    if (numinn == 1) {
      if (numfar == 1) {   /* inner escape, enter through top, escape far */
	*nempty = 1;
	dempty = tt[1]-ti[1];
	if (dempty>0) {
	  te[1] = ti[1];
	  te[2] = tt[1];
	} else {
	  dempty = - dempty;
	  te[1] = tt[1];
	  te[2] = ti[1];
	}
	if (dempty < 0) {
	  printf("2 %.3lf  %.3lf \n",te[1],te[2]);
	}
	*t = to;
	for(i=0; i<3; i++) { intpos[i] = xyzfar[i];}
	*ib = 2;
      } else {
	printf("numfar != 1  = %d after numtop=1 numinn=1 \n",numfar);
      }
    }
    if (numinn == 2) { /* cross inner, esc. thr. top -  impossible in 2 D */
      *nempty = 1;
      dempty = ti[2]-ti[1];
      te[1] = ti[1];
      te[2] = ti[2];
      if (te[2]<te[1]) {
	printf("4. %.3lf  %.3lf \n",te[1],te[2]);
      }
      *t = tt[1];
      for(i=0; i<3; i++) { intpos[i] = xyz2top[1][i]; }
      *ib = 3;
    }
  }

  if (numtop == 2) {
    if (numfar == 1) {
      if (numinn == 0)  {
	*nempty = 1;
	dempty = tt[2]-tt[1];
	te[1] = tt[1];
	te[2] = tt[2];
	if (dempty < 0) {
	  printf("3 %.3lf  %.3lf \n",te[1],te[2]);
	}
	*t = to;
	for(i=0; i<3; i++) { intpos[i] = xyzfar[i];}
	*ib = 2;
      } else {
	/*	printf("numtop=2; numfar=1 but numinn = %d   %d\n",numinn,debug);
		printf("%.2lf  %.2lf    %.2lf  %.2lf \n",ti[1],ti[2],tt[1],tt[2]); */
	*t = to;
	for(i=0; i<3; i++) { intpos[i] = xyzfar[i];}
	*ib = 2;
	*nempty = 2;
	if (tt[2] > ti[2]) {
	  te[1] = ti[1];
	  te[2] = ti[2];
	  te[3] = tt[1];
	  te[4] = tt[2];
	} else {
	  te[1] = tt[1];
	  te[2] = tt[2];
	  te[3] = ti[1];
	  te[4] = ti[2];
	}
	if (te[3] < te[2]) {
	  printf("problem with te: %.2lf %.2lf %.2lf %.2lf \n",te[1],te[2],te[3],te[4]);
	}
      }
    } else {
      printf("numfar != 1 after numtop=2 \n");
    }
  }

#ifdef DEBUG
  if (debug_cond) {
    printf("t = %.3lf,  nempty = %d, te: %.3lf %.3lf %.3lf %.3lf \n",*t,*nempty,
	   te[1],te[2],te[3],te[4]);
  }
#endif  
  
  return 0;
}

      
faredgestrait(curpos,curdir,numsol,xyzint,t)
double curpos[], curdir[], xyzint[],*t;
int    *numsol; 
{
  double a,b,c,p,q,d,z;

  a = curdir[0]*curdir[0] + curdir[1]*curdir[1];
  b = curpos[0]*curdir[0] + curpos[1]*curdir[1];
  /*  c < 0 */
  c = curpos[0]*curpos[0] + curpos[1]*curpos[1] - torus_r2;
      
  p = b/a;
  q = c/a;
      
  d = p*p-q;
      
  if (d < 0.) {
    /*    printf("d < 0 in faredge: %lf\n",d); */
    *numsol = 0; 
    return 0 ;
  }
      
  *t = -q/(sqrt(d)+p);
      
  z = curpos[2]+(*t)*curdir[2];

  if (fabs(z) <= torus_h2) {
    xyzint[0] = curpos[0] + (*t)*curdir[0];
    xyzint[1] = curpos[1] + (*t)*curdir[1];
    xyzint[2] = z;
    *numsol = 1;
    /*    printf(" far edge: t=%.3lf \n",*t); */
  } else {
    *numsol = 0;
  }
      
  return 0;
}

faredge(curpos,curdir,numsol,xyzint,t)
double curpos[], curdir[], xyzint[],*t;
int    *numsol; 
{
  double a,b,c,p,q,d,z,ds,t1,t2;

  a = curdir[0]*curdir[0] + curdir[1]*curdir[1] + curdir[2]*curdir[2];
  b = curpos[0]*curdir[0] + curpos[1]*curdir[1] + curpos[2]*curdir[2];
  /*  c < 0 */
  c = curpos[0]*curpos[0] + curpos[1]*curpos[1] + curpos[2]*curpos[2]- torus_r2;
      
  p = b/a;
  q = c/a;
      
  d = p*p-q;
      
  if (d < 0.) {
    /*    printf("d < 0 in faredge: %lf\n",d); */
    *numsol = 0; 
    return 0 ;
  }
      
  ds = sqrt(d);

  t1 = -p-ds;
  /*      t2 = -p+ds */
  t2 = -q/(p+ds);

  if (t2<0) {
    /*    printf("both solutions <0 in faredge \n"); */
    *numsol = 0;
    return 0;
  }

  if (t1>0) {
    *t = t1;
  } else {
    *t = t2;
  }

  z = curpos[2]+(*t)*curdir[2];

  if (fabs(z) <= torus_h2) {
    xyzint[0] = curpos[0] + (*t)*curdir[0];
    xyzint[1] = curpos[1] + (*t)*curdir[1];
    xyzint[2] = z;
    *numsol = 1;
  } else {
    *numsol = 0;
  }
      
  return 0;
}


nearedge(curpos,curdir,numsol,xyzint,ti)
double curpos[],curdir[],xyzint[3][3],ti[];
int *numsol;
{
  double a,b,c,p,q,d,ds,tt,z,t[3];
  int  i;
  
  ti[1] = 0;
  ti[2] = 0;
  a = curdir[0]*curdir[0] + curdir[1]*curdir[1];
  b = curpos[0]*curdir[0] + curpos[1]*curdir[1];
  /*    c > 0 */
  c = curpos[0]*curpos[0] + curpos[1]*curpos[1] - torus_r1;
  
  p = b/a;
  q = c/a;
      
  d = p*p-q;

      
  if (d < 0.) { 
    *numsol = 0;
    return 0;
  }
      
  ds = sqrt(d);
  /*      t(1) = -p+ds */
  t[1] = -q/(p+ds);
  t[2] = -p-ds;

  *numsol = 0;
  if (t[1]*t[2] < 0.) { 
    printf("Something very strange happened in nearedge.. nscat=%d\n",nscatt);
    writeposdir(curpos,curdir);
    printf("t1,2: %.3lf %.3lf \n",t[1],t[2]);
#ifdef DEBUG    
    scanf("%lf",&ds);
#endif
    return ;
  }
      
  if (t[1] < 0.) {
    /*    printf("near: t[1]<0  %lf \n",t[1]); */
    *numsol = 0;
    return 0;
  } else {
    if (t[1] > t[2]) {
      tt = t[1];
      t[1] = t[2];
      t[2] = tt;
    }
    for(i=1; i<=2; i++) {
      z = curpos[2]+t[i]*curdir[2];
      if (fabs(z) < torus_h1) { 
	(*numsol) ++;
	xyzint[(*numsol)][0] = curpos[0] + t[i]*curdir[0];
	xyzint[(*numsol)][1] = curpos[1] + t[i]*curdir[1];
	xyzint[(*numsol)][2] = z;
	ti[(*numsol)] = t[i];
      }
    }
  }

  return 0;
}

    
topbottom(pos,dir,numsol,xyz2int,t)
/* 
   for a torus with wedge-like shape

 */

double pos[],dir[],xyz2int[3][3],t[];
int *numsol;
{
  void exit();
  double a,b,c,p,q,d,x,y,z,rq,tt[2],t1,t2;
  int    nt,i;


  t[1] = 0;
  t[2] = 0;

  a = dir[0]*dir[0]+dir[1]*dir[1]-dir[2]*dir[2]*wedge_ctgt2;
  b = pos[0]*dir[0]+pos[1]*dir[1]-pos[2]*dir[2]*wedge_ctgt2;
  c = pos[0]*pos[0]+pos[1]*pos[1]-pos[2]*pos[2]*wedge_ctgt2;

  p = b/a;
  q = c/a;
      
  d = p*p-q;
      
  if (d < 0.) {
    /*    printf("d < 0 in topbottom: %lf\n",d); */
    *numsol = 0; 
    return 0;
  } else {
    if (d == 0.) {
      /*      printf("One degenerate solution \n"); */
      nt = 1;
      tt[1] = -p;
    } else {
      t1 = -p-sqrt(d);
      t2 = -p+sqrt(d);
      if (t1*t2 > 0) {
	if (t1 > 0) {
#ifdef DEBUG
	  if (debug_cond) 
	    printf("Both t's positive,  %.2lf  %.2lf\n",t1,t2); 
#endif
	  nt = 2;
	  tt[1] = t1;
	  tt[2] = t2;
	} else {
	  /*   printf("Both t's negative, no solutions %.2lf %.2lf\n",t1,t2); */
	  nt = 0;
	}
      } else {
	/*	printf("One solution %.2lf,   %.2lf\n",t2,t1); */
	nt = 1;
	tt[1] = t2;
      }
    }
  }

  *numsol = 0;
  for(i=1; i<=nt; i++) {
    x = pos[0]+tt[i]*dir[0];
    y = pos[1]+tt[i]*dir[1];
    z = pos[2]+tt[i]*dir[2];
    rq = x*x + y*y;
#ifdef DEBUG
    if (debug_cond) 
      printf("top: checking: %.2lf, %.2lf %.2lf; rq=%.3lf  %.3lf; %.3lf %.3lf\n",x,y,z,sqrt(rq),sqrt(rq+z*z),sqrt(torus_r1),sqrt(torus_r2)); 
#endif
    if ((rq > torus_r1) && (rq+z*z /* rounded outer edge */ < torus_r2)) {
      (*numsol) ++;
      t[*numsol] = tt[i];
      xyz2int[*numsol][0] = x;
      xyz2int[*numsol][1] = y;
      xyz2int[*numsol][2] = z;
#ifdef DEBUG
      if (debug_cond) 
	printf("     added ... %lf  %d \n",tt[i],(*numsol)); 
#endif
    }

    rq = sqrt(x*x+y*y+z*z);
    if (fabs(fabs(z/rq)/wedge_sintheta-1) > 1e-4) {
      printf("Inaccurate solution in top ???  %lf %lf\n",z/rq,wedge_sintheta);
    }
    
  } 

  
  return 0;
}


writeposdir(pos,dir)
double pos[],dir[];
{
  int i;
  double rq,Rq;

  printf("pos:  ");
  for(i=0; i<3; i++) {
    printf("%9.5lf  ",pos[i]);
  }
  printf("\n");
  printf("dir:  ");
  for(i=0; i<3; i++) {
    printf("%9.5lf  ",dir[i]);
  }
  printf("\n");
  rq = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
  Rq = sqrt(rq*rq+pos[2]*pos[2]);
  printf("rq : %lf,  Rq = %lf, theta=%lf %lf \n--\n",rq,Rq,
	 180/pi*acos(rq/Rq),180/pi*asin(pos[2]/Rq));

  return 0;
}


inside(pos)
double pos[];
{
  double rq;
  int    ins,ins1,ins2,ins3;

  rq = pos[0]*pos[0]+pos[1]*pos[1];
  ins1 = rq>torus_r1;
  ins2 = rq+pos[2]*pos[2] < torus_r2;
  ins3 = sqrt(rq)/sqrt(rq+pos[2]*pos[2]) > wedge_costheta; 

  ins = ins1 && ins2 && ins3;

#ifdef DEBUG
  if (!ins) {
    printf("Inside?\n%d %d %d \n",ins1,ins2,ins3);
    printf("%lf %lf \n",rq,torus_r1);
    printf("%lf %lf \n",rq+pos[2]*pos[2],torus_r2);
    printf("%lf %lf \n---\n",sqrt(rq)/sqrt(rq+pos[2]*pos[2]),wedge_costheta);
  }
#endif

  return ins;
}
