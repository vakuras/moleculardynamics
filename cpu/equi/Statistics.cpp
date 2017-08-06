#include "Statistics.h"

static const double rel_error= 1E-12;
static const double Cpi      = 3.1415926535897932384626434;
static const double C2sqrtPi = 1.1283791670955125738961589;
static const double CsqrtPi  = 1.7724538509055160272981675;
static const double Csqrt2   = 1.4142135623730950488016887;
static const double C1sqrt2  = 0.7071067811865475244008444;


double Phi (double z)
{
   return 0.5 * erfc( - C1sqrt2*z );
}

/***********************************************************************

   Inverse of the error function Erf.

   Implementation: Inversion by Newton iteration of erf(x).
      The initial value x0 = 0.
      For |z| <= 0.84 (=erf(1)) at most 4 iterations are necessary.

***********************************************************************/

static double InvErfSmall (const double z)
{
   /* f(x)   = erf(x) - z   */
   /* f'(x)  = c*exp(-x*x)  */
   /* f''(x) = -2 f'(x)     */
   double c = C2sqrtPi;
   double f = -z, f1=c;
   double q = f/f1, x = -q, x0 = 0;

   while (fabs(x-x0) > 1e-12 && fabs(f) > 1e-14 ) {
      /* Newton 2nd order: x <- x - f/f'(1 + f*f''/(2 f'^2)) */
      x0  = x;
      f   = erf(x) - z;
      f1  = c*exp(-x*x);
      q   = f/f1;
      x  -= q*(1-x*q);  /* Newton Step 2nd order */
   }

   return x;
}

/***********************************************************************

 lambert_w2 is that "Dilbert lambda" (inverse x e^x^2) 
 If you don't have W functions, you can probably get away with 

                                     log(log(x)) 
     lambert_w2(x) ~ sqrt(log(x)) - -------------- . 
                                    4 sqrt(log(x)) 

***********************************************************************/

static double lambertW2 (const double z)
{
   /* "Dilbert lambda" (inverse x e^x^2) approximation */
   double logz = log(z);
   double slz  = sqrt(logz);
   return slz - 0.25*log(logz)/slz;
}


/***********************************************************************

   Inverse of the error function Erfc.

   Implementation: Inversion by Newton iteration of erfc(sqrt(log(x))).
      The initial value is computed via lambertW2.
      For z < 0.25 at most 4 iterations are necessary.

***********************************************************************/

static double InvErfcSmall (const double z)
{
   /* f(x)   = erfc(sqrt(log(x))) - z   */
   /* f'(x)  = 1/(c x^2 sqrt(log(x)))   */
   /* f''(x) = -c*x*f'(x)^2*(2*sqrt(log(x))+1/(2*sqrt(log(x))))   */
   double c = CsqrtPi;
   double f = 1, f1i=0;
   double a = lambertW2(1/(c*z));
   double x = erfc(a)/(c * a * z*z), x0 = 0;

   while (fabs(x-x0) > 1e-12*x && fabs(f) > 1e-15*z ) {
      /* Newton 2nd order: x <- x - f/f'(1 + f*f''/(2 f'^2)) */
      double slx;
      x0  = x;
      slx = sqrt(log(x));
      f   = z - erfc(slx);
      f1i = c * x * x * slx;
      x  -= f*f1i*(1 - c*f*x*(slx + 0.25/slx));       /* Newton Step */
   }

   return  sqrt(log(x));
}

/***********************************************************************

   Inverse of the error function Erf.

   Implementation: 
      For small and big values two differnt approximations are used.
      
   Return values:
      If x is +1, InvErf() returns +INFINITY.
      If x is -1, InvErf() returns -INFINITY.
      If x is not in the range [-1, 1],
      InvErf() returns NaN and sets errno to EDOM.

***********************************************************************/

double inverf (const double z)
{
   double az = fabs(z);
   double x = 0;
   if (az>=1)
   {
      double huge = DBL_MAX;
      double inf  = 2*huge;
      /* Infinity of NaN */
      if (z == 1)
	 return  inf; /* +infinity */
      if (z == -1)
	 return -inf; /* -infinity */
      /* out of range, set errno=EDOM, return NaN */
      errno=EDOM;
      return inf-inf;
   }

   if (az < 0.8125)          /* -13/16 < z < 13/16 */
      x =  InvErfSmall(z);
   else if (z > 0)           /* z >=  13/16 */
      x =  InvErfcSmall(1-z);
   else                      /* z <= -13/16 */
      x = -InvErfcSmall(z+1);
   return x;
}

double inverfc (const double z)
{
   double x = 0;
   if ( z<=0 || z>=2 )
   {
      double huge = DBL_MAX;
      double inf  = 2*huge;
      /* Infinity of NaN */
      if (z == 0)
	 return  inf; /* +infinity */
      if (z == 2)
	 return -inf; /* -infinity */
      /* out of range, set errno=EDOM, return NaN */
      errno=EDOM;
      return inf-inf;
   }

   if (z <= 0.1875)          /* z <=   3/16 */
      x =  InvErfcSmall(z);
   else if (z < 1.8125)      /* 3/16 < z < 29/16 */
      x =  InvErfSmall(1-z);
   else                      /* z >=  29/16 */
      x = -InvErfcSmall(2-z);
   return x;
}

/***********************************************************************

   Inverse of the cumulative normal distribution

   Implementation: 

***********************************************************************/

double erf(double x)
//erf(x) = 2/sqrt(pi)*integral(exp(-t^2),t,0,x)
//       = 2/sqrt(pi)*[x - x^3/3 + x^5/5*2! - x^7/7*3! + ...]
//       = 1-erfc(x)
{
    static const double two_sqrtpi=  1.128379167095512574;        // 2/sqrt(pi)
    if (fabs(x) > 2.2) {
        return 1.0 - erfc(x);        //use continued fraction when fabs(x) > 2.2
    }
    double sum= x, term= x, xsqr= x*x;
    int j= 1;
    do {
        term*= xsqr/j;
        sum-= term/(2*j+1);
        ++j;
        term*= xsqr/j;
        sum+= term/(2*j+1);
        ++j;
    } while (fabs(term)/sum > rel_error);
    return two_sqrtpi*sum;
}


double erfc(double x)
//erfc(x) = 2/sqrt(pi)*integral(exp(-t^2),t,x,inf)
//        = exp(-x^2)/sqrt(pi) * [1/x+ (1/2)/x+ (2/2)/x+ (3/2)/x+ (4/2)/x+ ...]
//        = 1-erf(x)
//expression inside [] is a continued fraction so '+' means add to denominator only
{
    static const double one_sqrtpi=  0.564189583547756287;        // 1/sqrt(pi)
    if (fabs(x) < 2.2) {
        return 1.0 - erf(x);        //use series when fabs(x) < 2.2
    }
    if (x<0) {               //continued fraction only valid for x>0
        return 2.0 - erfc(-x);
    }
    double a=1, b=x;                //last two convergent numerators
    double c=x, d=x*x+0.5;          //last two convergent denominators
    double q1, q2= b/d;             //last two convergents (a/c and b/d)
    double n= 1.0, t;
    do {
        t= a*n+b*x;
        a= b;
        b= t;
        t= c*n+d*x;
        c= d;
        d= t;
        n+= 0.5;
        q1= q2;
        q2= b/d;
      } while (fabs(q1-q2)/q2 > rel_error);
    return one_sqrtpi*exp(-x*x)*q2;
}