#include <stdio.h>
#include <math.h>
#include "salt.h"

double density0(double s,double t) {
  double A,B,C,D,dens0;
  A = 1.001685e-04 + t * ( -1.120083e-06 + t * 6.536332e-09 );
  A = 999.842594 + t * (  6.793952e-02 + t * ( -9.095290e-03 + t * A ) );
  B = 7.6438e-05 + t * ( -8.2467e-07 + t * 5.3875e-09 );
  B = 0.824493 + t * ( -4.0899e-03 + t * B );
  C = -5.72466e-03 + t * ( 1.0227e-04 - t * 1.6546e-06 );
  D = 4.8314e-04;
  dens0 = A + s * (  B + C * sqrt(s) + D * s );
  return dens0;
}

double sal2dens(double s,double t,double p) {
  double d0,d,K,E,F,G,H,I,J,M,N,t2,t3,t4,s1p5,pb;
  t2 = t * t;t3 = t2 * t;t4 = t3 * t;
  d0 = density0(s,t);
  E = 19652.21 + 148.4206 * t - 2.327105 * t2 + 1.360477e-2 * t3 - 5.155288e-5 * t4;
  F = 54.6746 - 0.603459 * t + 1.09987e-2 * t2 - 6.1670e-5 * t3;
  G = 7.944e-2 + 1.6483e-2 * t - 5.3009e-4 * t2;
  H = 3.239908 + 1.43713e-3 * t + 1.16092e-4 * t2 - 5.77905e-7 * t3;
  I = 2.2838e-3 - 1.0981e-5 * t - 1.6078e-6 * t2;
  J = 1.91075e-4;
  M = 8.50935e-5 - 6.12293e-6 * t + 5.2787e-8 * t2;
  N = -9.9348e-7 + 2.0816e-8 * t + 9.1697e-10 * t2;
  s1p5 = s * sqrt(s);
  pb = p/10;
  K = (E + F*s + G*s1p5) + (H + I*s + J*s1p5) * pb + (M + N*s) * pb * pb;
  d = d0 / (1 - pb/K);
  return d/1000.0;
}

double dens2sal(double d,double t,double p) {
  d=d*1000.0;
  double sal=0,slo=0.0,shi=100.0,err=1.0,dNew,i=0,limit=1000;

  while( (fabs(err) > 0.00001) && (i < limit) ) {
    sal = (shi+slo) / 2.0;
    dNew = sal2dens(sal,t,p)*1000.0;
    err = (dNew-d) / d;
    if (err > 0.0) shi = sal;
    else slo = sal;
    ++i;
  }
  if ( i >= limit ) return 0.0;
  else return sal;
}

double s_Rt(double t,double Rt) {
  double Rt5,t15,dels,sal;
  Rt5 = sqrt( Rt );
  t15 = t - 15;
  dels = t15 / ( 1 + 0.0162 * t15 );
  sal =  ( 14.0941 + dels * -0.0375 ) + Rt5 * ( ( -7.0261 + dels *  0.0636 )
         + Rt5 * ( (  2.7081 + dels *  -0.0144 ) ) );
  sal = ( 0.008 + dels * 0.0005 ) + Rt5 * ( ( -0.1692 + dels * -0.0056 )
        + Rt5 * ( ( 25.3851 + dels * -0.0066 ) + Rt5 * sal ) );
  return sal;
}

double  cond2sal(double c,double t,double p) {
  c=c/10.0;
  double R,rt,Rp,Rt,A,B,C,sal;
  R = c / 4.29140;
  rt = 0.6766097 + t * ( 0.0200564 + t * ( 1.104259e-04 
       + t * ( -6.9698e-07 + t * 1.0031e-09 ) ) );
  A = 0.4215 - 0.003107 * t;
  B = 1 + t * ( 0.03426 + t * 0.0004464 );
  C = p * ( 2.07e-5 + p * ( -6.37e-10 + p * 3.989e-15 ) );
  Rp = 1 + C / ( B + A * R );
  Rt = R / rt / Rp;
  sal = s_Rt(t,Rt);
  return sal;
}

double sal2cond(double s,double t,double p) {
  double R,Rt=0,rt,A,B,C;
  double Rtlo = 0.0,Rthi = 10.0,err = 1.0,sNew,limit = 1000;
  int i = 0;
  
  if( s < 0.01 ) return 0.0;
  
  while( (fabs(err) > 0.00001) && (i < limit) ) {
    Rt = (Rthi+Rtlo) / 2.0;
    sNew = s_Rt(t,Rt);
    err = (sNew-s) / s;
    if (err > 0.0) Rthi = Rt;
    else Rtlo = Rt;
    ++i;
  }
  if ( i >= limit )
    return 0.0;
  
  rt = 0.6766097 + t * ( 0.0200564 + t * ( 1.104259e-04 
       + t * ( -6.9698e-07 + t * 1.0031e-09 ) ) );
  A = 0.4215 - 0.003107 * t;
  B = 1 + t * ( 0.03426 + t * 0.0004464 );
  C = p * ( 2.07e-5 + p * ( -6.37e-10 + p * 3.989e-15 ) );
  
  R = ( sqrt( (A*rt*Rt-B)*(A*rt*Rt-B) + 4*rt*Rt*A*(B+C) ) + (A*rt*Rt-B) ) / (2*A);
  return R * 4.29140 * 10.0;
}

double sal2dens_teos10(double sa,double ct,double p) {
  double v01 =  9.998420897506056e+2, v02 =  2.839940833161907;
  double v03 = -3.147759265588511e-2, v04 =  1.181805545074306e-3;
  double v05 = -6.698001071123802, v06 = -2.986498947203215e-2;
  double v07 =  2.327859407479162e-4, v08 = -3.988822378968490e-2;
  double v09 =  5.095422573880500e-4, v10 = -1.426984671633621e-5;
  double v11 =  1.645039373682922e-7, v12 = -2.233269627352527e-2;
  double v13 = -3.436090079851880e-4, v14 =  3.726050720345733e-6;
  double v15 = -1.806789763745328e-4, v16 =  6.876837219536232e-7;
  double v17 = -3.087032500374211e-7, v18 = -1.988366587925593e-8;
  double v19 = -1.061519070296458e-11, v20 =  1.550932729220080e-10;
  double v21 =  1.0, v22 =  2.775927747785646e-3, v23 = -2.349607444135925e-5;
  double v24 =  1.119513357486743e-6, v25 =  6.743689325042773e-10;
  double v26 = -7.521448093615448e-3, v27 = -2.764306979894411e-5;
  double v28 =  1.262937315098546e-7, v29 =  9.527875081696435e-10;
  double v30 = -1.811147201949891e-11, v31 = -3.303308871386421e-5;
  double v32 =  3.801564588876298e-7, v33 = -7.672876869259043e-9;
  double v34 = -4.634182341116144e-11, v35 =  2.681097235569143e-12;
  double v36 =  5.419326551148740e-6, v37 = -2.742185394906099e-5;
  double v38 = -3.212746477974189e-7, v39 =  3.191413910561627e-9;
  double v40 = -1.931012931541776e-12, v41 = -1.105097577149576e-7;
  double v42 =  6.211426728363857e-10, v43 = -1.119011592875110e-10;
  double v44 = -1.941660213148725e-11, v45 = -1.864826425365600e-14;
  double v46 =  1.119522344879478e-14, v47 = -1.200507748551599e-15;
  double v48 =  6.057902487546866e-17;

  double sqrtsa, v_hat_denominator, v_hat_numerator, gsw_rho;

  sqrtsa = sqrt(sa);

  v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))  
                      + sa*(v05 + ct*(v06 + v07*ct) 
                      + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct)))) 
                      + p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct) 
                      + p*(v17 + ct*(v18 + v19*ct) + v20*sa));
  
  v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) 
                    + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa 
                    + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))  
                    + p*(v37 + ct*(v38 + ct*(v39 + v40*ct))  
                    + sa*(v41 + v42*ct) 
                    + p*(v43 + ct*(v44 + v45*ct + v46*sa) 
                    + p*(v47 + v48*ct)));

  gsw_rho = v_hat_denominator/v_hat_numerator;
  return gsw_rho/1000.0;
}

double dens2sal_teos10(double d,double t,double p) {
  int i=0,limit=1000;
  double sal=0,slo=0.0,shi=100.0,err=1.0,dNew;
  d=d*1000.0;
  while( (fabs(err) > 0.00000001) && (i < limit) ) {
    sal = (shi+slo) / 2.0;
    dNew = sal2dens_teos10(sal,t,p)*1000.0;
    err = (dNew-d) / d;
    if (err > 0.0) shi = sal;
    else slo = sal;
    ++i;
  }
  if ( i >= limit ) return 0.0;
  else return sal;  
}
