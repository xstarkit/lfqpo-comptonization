//#include <math.h>
//#include <stdio.h>

#ifdef debug
FILE *fd;
#endif
/*
 Procedures for integrating ordinary differential equations (Runge-Kutta
 method, Bulirsch-Stoer method), adopted from "Numerical Receipes".

 if you change MAXV or MAXL here, you must also change them in the main program
*/
#define MAXV 3 
#define MAXL 1000
extern  double xp[MAXL],yp[MAXV][MAXL];

/*        This is main procedure                                    */
odeint(method,derivs,ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,nop)
   void   (*derivs)();
   double ystart[],x1,x2,eps,h1,hmin;
   int    method,nvar,*nok,*nbad,*nop;
/*    
    method:   1 - Runge-Kutta method
              0 - Bulirsch-Stoer method
*/
{
   double  maxstp = 10000,
           two    = 2,
           zero   = 0,
           tiny   = 1E-30;
   double  xsav,x,hnext,hdid,h,
           yscal[MAXV],y[MAXV],dydx[MAXV],
           dxsav;
   int     nstp,i,
           kmax,kount;
   void    exit();

   if (nvar > MAXV)
   {
     printf("Number of variables greater than arrays' dimensions \n");
     exit(13);
   }
#ifdef debug
   fd = fopen("debug_c.inf","w");
   printf("Debug file `debug.inf' created \n");
#endif
   kmax = MAXL;
   dxsav = fabs(x2-x1)/kmax;
   x = x1;
   if (x2 > x1)
     h = fabs(h1);
   else
     h = -fabs(h1);
   *nok = 0;
   *nbad = 0;
   kount = 0;
   for (i=0; i<nvar; i++)
      y[i] = ystart[i];
   xsav = x-dxsav*two;
   nstp = 1;
   for (nstp=1; nstp<maxstp; nstp++)
   {
      derivs(x,y,dydx);
      for (i=0; i<nvar; i++)
         yscal[i] = fabs(y[i])+fabs(dydx[i]*h)+tiny;
      if (kmax > 0)
      {
         if (fabs(x-xsav) > fabs(dxsav))
         {
            if (kount < kmax-1)
            {
               xp[kount] = x;
               for (i=0; i<nvar; i++)
                  yp[i][kount] = y[i];
               xsav = x;
               ++kount;
            }
         }
      }
      if (((x+h-x2)*(x+h-x1)) > zero)
         h = x2-x;
#ifdef debug
      printf("%10f \n",x);
      fprintf(fd,"x=%10f \n",x);
#endif
      if (method % 2)
           rkqc(derivs,y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext);
      else
         bsstep(derivs,y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext);
      if (hdid == h)
         ++(*nok);
      else
         ++(*nbad);
      if (((x-x2)*(x2-x1)) >= zero)
      {
         for (i=0; i<nvar; i++)
            ystart[i] = y[i];
         if (kmax != 0)
         {
            xp[kount] = x;
            for (i=0; i<nvar; i++)
               yp[i][kount] = y[i];
            ++kount;
         }
         goto norm_exit;
      }
      if (fabs(hnext) < hmin)
      {
         printf("pause in routine ODEINT\n");
         printf("stepsize to small \n");
         getchar();
      }
      h = hnext;
   }
   printf("pause in routine ODEINT - too many steps \n");
   getchar();
   norm_exit:
#ifdef debug
   fclose(fd);
#endif
   *nop = kount;
   return 0;
}

rk4(derivs,y,dydx,n,x,h,yout)
  void   (*derivs)();
  double y[],dydx[],yout[];
  double x,h;
  int    n;
{
   int    i;
   double xh,hh,h6,
          dym[MAXV],dyt[MAXV],yt[MAXV];

   hh = h*0.5;
   h6 = h/6.0;
   xh = x+hh;
   for (i=0; i<n; i++)
      yt[i] = y[i]+hh*dydx[i];
   derivs(xh,yt,dyt);
   for (i=0; i<n; i++)
      yt[i] = y[i]+hh*dyt[i];
   derivs(xh,yt,dym);
   for (i=0; i<n; i++)
   {
      yt[i] = y[i]+h*dym[i];
      dym[i] = dyt[i]+dym[i];
   }
   derivs(x+h,yt,dyt);
   for (i=0; i<n; i++)
      yout[i] = y[i]+h6*(dydx[i]+dyt[i]+2*dym[i]);

   return 0;
}

rkqc(derivs,y,dydx,n,x,htry,eps,yscal,hdid,hnext)
   void   (*derivs)();
   double y[],dydx[],*x,htry,eps,yscal[],*hdid,*hnext;
   int    n;
{
   double  pgrow = -0.20,
           pshrnk= -0.25,
           fcor  =  0.06666666,   /* 1/15 */
           one   =  1.0,
           safety=  0.9,
           errcon=  6E-4;
   double  xsav,hh,h,temp,errmax,
           dysav[MAXV],ysav[MAXV],ytemp[MAXV];
   int     i;

   xsav = *x;
   for (i=0; i<n; i++)
   {
      ysav[i] = y[i];
      dysav[i] = dydx[i];
   }
   h = htry;
   do
   {
      hh = 0.5*h;
      rk4(derivs,ysav,dysav,n,xsav,hh,ytemp);
      *x = xsav+hh;
      derivs(*x,ytemp,dydx);
      rk4(derivs,ytemp,dydx,n,*x,hh,y);
      *x = xsav+h;
      if (*x == xsav)
      {
         printf("pause in routine RKQC \n");
         printf("stepsize to small \n");
         getchar();
      }
      rk4(derivs,ysav,dysav,n,xsav,h,ytemp);
      errmax = 0;
      for (i=0; i<n; i++)
      {
         ytemp[i] = y[i]-ytemp[i];
         temp = fabs(ytemp[i]/yscal[i]);
         if (errmax < temp)
           errmax = temp;
      }
      errmax /= eps;
      if (errmax > one)
        h = safety*h*exp(pshrnk*log(errmax));
   }
   while (errmax > one);

   *hdid = h;
   if (errmax > errcon)
      *hnext = safety*h*exp(pgrow*log(errmax));
   else
      *hnext = 4.0*h;

   for (i=0; i<n; i++)
      y[i] += ytemp[i]*fcor;

   return 0;
}


rk_dumb(derivs,vstart,nvar,x1,x2,nstep)
  void    (*derivs)();
  double  vstart[MAXV],x1,x2;
  int     nvar,nstep;
{
   int    kount,k,i;
   double x,h,v[MAXV],vout[MAXV],dv[MAXV];

   for (i=0; i<nvar; i++)
   {
      v[i] = vstart[i];
      yp[i][1] = v[i];
   }
   xp[1] = x1;
   x = x1;
   h = (x2-x1)/nstep;
   for (k=1; k<=nstep; k++)
   {
      derivs(x,v,dv);
      rk4(derivs,v,dv,nvar,x,h,vout);
      if (x+h == x)
      {
         printf("pause in routine RKDUMB \n");
         printf("stepsize to small \n");
         getchar();
      }
      x = x+h;
      xp[k+1] = x;
      for (i=0; i<nvar; i++)
      {
         v[i] = vout[i];
         yp[i][k+1] = v[i];
      }
   }
   return 0;
}

bsstep(derivs,y,dydx,nv,x,htry,eps,yscal,hdid,hnext)
   void    (*derivs)();
   double  y[],dydx[],*x,htry,eps,yscal[],*hdid,*hnext;
   int     nv;
{
   int     imax  = 11,
           nuse  = 7;
   double  one   = 1,
           shrink= 0.95,  
           grow  = 1.2;
   int     nseq[12]; 
   int     i,j,exit_cond,im;
   double  xsav,xest,h,errmax,
           ysav[MAXV],dysav[MAXV],yseq[MAXV],yerr[MAXV];
   
   nseq[0]=1; nseq[1]=2; nseq[2]=4; nseq[3]=6; nseq[4]=8; nseq[5]=12;
   nseq[6]=16; nseq[7]=24; nseq[8]=32; nseq[9]=48; nseq[10]=64; nseq[11]=96;
   h = htry;
   xsav = *x;
   for (i=0; i<nv; i++)
   {
      ysav[i] = y[i];
      dysav[i] = *(dydx+i);
   }
   exit_cond = 0;
   do
   {
      i = 0;
#ifdef debug
         fprintf(fd,"i  errmax ");
#endif
      while ( (i<imax) && (!exit_cond))
      {
         i++;
         mmid(derivs,ysav,dysav,nv,xsav,h,nseq[i],yseq);
         xest = h/nseq[i];
         xest *= xest;
         rzextr(i,xest,yseq,y,yerr,nv,nuse);
         errmax = 0.0;
         for (j=0; j<nv; j++)
         {
            if (errmax < fabs(yerr[j]/yscal[j])) 
               errmax = fabs(yerr[j]/yscal[j]);
         }
         errmax = errmax/eps;
#ifdef debug
	 printf("%4i %6e ",i,errmax);
	 fprintf(fd,"%4i%4i %6e ",i,nseq[i],errmax);
#endif
         if (errmax < one) 
         {
            *x += h;
            *hdid = h;
            if (i == nuse) 
              *hnext = h*shrink;
            else
              if (i == nuse-1)
                 *hnext = h*grow;
              else
                 *hnext = (h*nseq[nuse-1])/nseq[i];
            exit_cond = 1;
         }
      }
#ifdef debug
         fprintf(fd,"\n");
         printf("\n");
#endif
      if (! exit_cond)
      {
	 h = 0.25*h;
	 im = (imax-nuse)/2;
	 if ( im > 0 )
	    for (i=1; i<=im; i++)
              h = h/2;
         if ((*x)+h == *x)
         {
            printf("pause in routine BSSTEP \n");
            printf("step size underflow \n");
            getchar();
         }
      }
   }
   while (! exit_cond);

   return 0;
}

mmid(derivs,y,dydx,nvar,xs,htot,nstep,yout)
   void   (*derivs)();
   double y[],dydx[],xs,htot,yout[];
   int    nvar,nstep;
{
   double x,swap,h2,h,
          ym[MAXV],yn[MAXV];
   int    n,i; 

   h = htot/nstep;
   for (i=0; i<nvar; i++)
   {
      ym[i] = y[i];
      yn[i] = y[i]+h*dydx[i];
   }
   x = xs+h;
   derivs(x,yn,yout);
   h2 = 2*h;
   for (n=2; n<=nstep; n++)
   {
      for (i=0; i<nvar; i++)
      {
          swap = ym[i]+h2*yout[i];
          ym[i] = yn[i];
          yn[i] = swap;
      }
      x = x+h;
      derivs(x,yn,yout);
   }
   for (i=0; i<nvar; i++)
      yout[i] = 0.5*(ym[i]+yn[i]+h*yout[i]);

   return 0;
}

#define ncol     7
#define glimax  11
#define glnmax  10
#define glncol   7

rzextr(iest,xest,yest,yz,dy,nv,nuse)
   double xest,yest[],yz[],dy[];
   int    iest,nv,nuse;
{
   int     m1,k,j,j1;
   double  yy,v,ddy,c,b1,b,fx[ncol+1];

   static  double  glx[glimax+1],
                   gld[glnmax+1][ncol+1];

   glx[iest] = xest;
   if (iest == 1) 
   {
      for (j=0; j<nv; j++)
      {
         yz[j] = yest[j];
         gld[j][1] = yest[j];
         dy[j] = yest[j];
      }
   }
   else
   {
      if (iest < nuse) 
         m1 = iest;
      else
         m1 = nuse;
      for (k=1; k<=m1-1; k++)
         fx[k+1] = glx[iest-k]/xest;
      for (j=0; j<nv; j++)
      {
         yy = yest[j];
         v = gld[j][1];
         c = yy;
         gld[j][1] = yy;
         for (k=2; k<=m1; k++)
         {
            b1 = fx[k]*v;
            b = b1-c;
            if (b != 0.0)
            {
               b = (c-v)/b;
               ddy = c*b;
               c = b1*b;
            }
            else
            {
               ddy = v;
            }
            v = gld[j][k];
            gld[j][k] = ddy;
            yy = yy+ddy;
         }
         dy[j] = ddy;
         yz[j] = yy;
      }
   }

   return 0;
}
#include <math.h>

/*
   Various mathods of numerical quadratures

   functions:

  double Hermite_Quadrature(f)  -  32-point Gauss-Hermite quadrature  (MatLib)
  double Laguerre_Quadrature(f) -  32-point Gauss-Lagueree quadrature (MatLib)

  double qsimp(double (*f)(), double a, double b)  - procedure from Num.Rec

*/

double Hermite_Quadrature(f)
  double (*f)();
/*  Gauss-Hermite quadrature, rank 32 */
{
    double x,y;

    x = 0.71258139098307276E+01;
    y = 0.7310676427384162E-22*(f(x)+f(-x));
    x = 0.64094981492696604E+01;
    y += 0.9231736536518292E-18*(f(x)+f(-x));
    x = 0.58122259495159138E+01;
    y += 0.11973440170928487E-14*(f(x)+f(-x));
    x = 0.52755509865158801E+01;
    y += 0.42150102113264476E-12*(f(x)+f(-x));
    x = 0.47771645035025964E+01;
    y += 0.59332914633966386E-10*(f(x)+f(-x));
    x = 0.43055479533511984E+01;
    y += 0.40988321647708966E-08*(f(x)+f(-x));
    x = 0.38537554854714446E+01;
    y += 0.15741677925455940E-06*(f(x)+f(-x));
    x = 0.34171674928185707E+01;
    y += 0.36505851295623761E-05*(f(x)+f(-x));
    x = 0.29924908250023742E+01;
    y += 0.54165840618199826E-04*(f(x)+f(-x));
    x = 0.25772495377323175E+01;
    y += 0.53626836552797205E-03*(f(x)+f(-x));
    x = 0.21694991836061122E+01;
    y += 0.36548903266544281E-02*(f(x)+f(-x));
    x = 0.17676541094632016E+01;
    y += 0.17553428831573430E-01*(f(x)+f(-x));
    x = 0.13703764109528718E+01;
    y += 0.60458130955912614E-01*(f(x)+f(-x));
    x = 0.9765004635896828E+00;
    y += 0.15126973407664248E+00*(f(x)+f(-x));
    x = 0.58497876543593245E+00;
    y += 0.27745814230252990E+00*(f(x)+f(-x));
    x = 0.19484074156939933E+00;
    y += 0.37523835259280239E+00*(f(x)+f(-x));

    return(y);
}

double Laguerre_Quadrature(f)
  double (*f)();
/*  Gauss-Laguerre quadrature, rank 32 */
{
    double x,y;

      x = 0.11175139809793770E+03;
      y = 0.45105361938989742E-47*f(x);
      x = 0.9882954286828397E+02;
      y += 0.13386169421062563E-41*f(x);
      x = 0.8873534041789240E+02;
      y += 0.26715112192401370E-37*f(x);
      x = 0.8018744697791352E+02;
      y += 0.11922487600982224E-33*f(x);
      x = 0.7268762809066271E+02;
      y += 0.19133754944542243E-30*f(x);
      x = 0.65975377287935053E+02;
      y += 0.14185605454630369E-27*f(x);
      x = 0.59892509162134018E+02;
      y += 0.56612941303973594E-25*f(x);
      x = 0.54333721333396907E+02;
      y += 0.13469825866373952E-22*f(x);
      x = 0.49224394987308639E+02;
      y += 0.20544296737880454E-20*f(x);
      x = 0.44509207995754938E+02;
      y += 0.21197922901636186E-18*f(x);
      x = 0.40145719771539442E+02;
      y += 0.15421338333938234E-16*f(x);
      x = 0.36100494805751974E+02;
      y += 0.8171823443420719E-15*f(x);
      x = 0.32346629153964737E+02;
      y += 0.32378016577292665E-13*f(x);
      x = 0.28862101816323475E+02;
      y += 0.9799379288727094E-12*f(x);
      x = 0.25628636022459248E+02;
      y += 0.23058994918913361E-10*f(x);
      x = 0.22630889013196774E+02;
      y += 0.42813829710409289E-09*f(x);
      x = 0.19855860940336055E+02;
      y += 0.63506022266258067E-08*f(x);
      x = 0.17292454336715315E+02;
      y += 0.7604567879120781E-07*f(x);
      x = 0.14931139755522557E+02;
      y += 0.7416404578667552E-06*f(x);
      x = 0.12763697986742725E+02;
      y += 0.59345416128686329E-05*f(x);
      x = 0.10783018632539972E+02;
      y += 0.39203419679879472E-04*f(x);
      x = 0.8982940924212596E+01;
      y += 0.21486491880136419E-03*f(x);
      x = 0.7358126733186241E+01;
      y += 0.9808033066149551E-03*f(x);
      x = 0.59039585041742439E+01;
      y += 0.37388162946115248E-02*f(x);
      x = 0.46164567697497674E+01;
      y += 0.11918214834838557E-01*f(x);
      x = 0.34922132730219945E+01;
      y += 0.31760912509175070E-01*f(x);
      x = 0.25283367064257949E+01;
      y += 0.70578623865717442E-01*f(x);
      x = 0.17224087764446454E+01;
      y += 0.12998378628607176E+00*f(x);
      x = 0.10724487538178176E+01;
      y += 0.19590333597288104E+00*f(x);
      x = 0.57688462930188643E+00;
      y += 0.23521322966984801E+00*f(x);
      x = 0.23452610951961854E+00;
      y += 0.21044310793881323E+00*f(x);
      x = 0.44489365833267018E-01;
      y += 0.10921834195238497E+00*f(x);

      return(y);
}

trapzd(f,a,b,s,n)
   double a,b,*s,(*f)();
   int    n;
{
   int    j;
   double x,tnm,sum,del;
   static int it;

   if (n == 1) 
   {
      *s = 0.5*(b-a)*(f(a)+f(b));
      it = 1;
   }
   else
   {
      tnm = it;
      del = (b-a)/tnm;
      x = a+0.5*del;
      sum = 0.0;
      for (j=1; j<=it; j++)
      {
         sum += f(x);
         x += del;
      }
      *s = 0.5*(*s+(b-a)*sum/tnm);
      it *= 2;
   }
   return 0;
}
trapzd1(f,a,b,s,n)
   double a,b,*s,(*f)();
   int    n;
{
   int    j;
   double x,tnm,sum,del;
   static int it;

   if (n == 1) 
   {
      *s = 0.5*(b-a)*(f(a)+f(b));
      it = 1;
   }
   else
   {
      tnm = it;
      del = (b-a)/tnm;
      x = a+0.5*del;
      sum = 0.0;
      for (j=1; j<=it; j++)
      {
         sum += f(x);
         x += del;
      }
      *s = 0.5*(*s+(b-a)*sum/tnm);
      it *= 2;
   }
   return 0;
}

double accuracy = 1e-6;

set_accuracy(eps)
   double eps;
{
   accuracy = eps;
   return 0;
}

double qsimp(f,a,b)
   double a,b,(*f)();
{
   int     jmax= 30;
   double  os  = -1E30,
           ost = -1E30,
           st,s;
   int     j,exit_cond;

   j = 0;
   do
   {
      j++;
      trapzd(f,a,b,&st,j);
      s = (4*st-ost)/3;
      exit_cond = (fabs(s-os) < accuracy*fabs(os));
      if (! exit_cond)
      {
        os = s;
        ost = st;
      }
   }
   while ((j<jmax) && (! exit_cond));
   if (! exit_cond)
   {
     printf("pause in QSIMP - too many steps \n"); 
     getchar();
   }
   return s;
}


#define EPS 1.0e-5
#define JMAX 20

double qtrap(double (*func)(double), double a, double b)
{
  /*  double trapzd(double (*func)(double), double a, double b, int n); */
  int j;
  double s,olds;
  
  olds = -1.0e30;
  for (j=1;j<=JMAX;j++) {
    trapzd1(func,a,b,&s,j);
    if (fabs(s-olds) < EPS*fabs(olds)) return s;
    olds=s;
  }
  printf("Too many steps in routine qtrap\n");
  return 0.0;
}
#undef EPS
#undef JMAX
/* (C) Copr. 1986-92 Numerical Recipes Software . */
