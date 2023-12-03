#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "TP3.h"

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */


static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */


/***************************************************************************************/
/*                                                                                     */
/* SimPi                     Calcul de PI par la methode de Monte Carlo                */
/*                                                                                     */
/* On tire au hasard des points dans le carré unité [0...1] x[0...1]                   */
/* On compte tous les points qui sont sur le quart de disque unité                     */
/* La proportion de ces points par rapport au nombre total de points tirés             */
/* converge vers PI/4                                                                  */
/*                                                                                     */
/*                                                                                     */
/*  En Entrée : nbr_points le nombre de points pour l'estimation de PI                 */
/*                                                                                     */
/*  En sortie:  une valeur approchée de PI en double précision                         */
/*                                                                                     */
/*                                                                                     */ 
/*                                                                                     */
/*                                                                                     */
/***************************************************************************************/
double simPi(long nbr_points){
  
  long i;
  
  double xr=0;
  
  double yr=0;
  
  double Pi=0;
  
  long inDisk=0;

  for(i=0; i< nbr_points; i++){
    
    xr=genrand_real1();
    
    yr=genrand_real1();

    if((xr * xr + yr * yr) <= 1){
      
      inDisk++;
    }
  }
  
  Pi=((double)inDisk/(double)nbr_points)*4;

  return Pi;
  
  }

/**********************************************************************************************************************************************/
/*                                                                                                                                            */
/* Experience_moyenne             Calcul d'expériences indépendantes et obtention de PI moyen                                                 */
/*                                                                                                                                            */
/* elle utilise une boucle pour réaliser un certain nombre d'expériences avec une fonction appelée "simPi".                                   */
/* Chaque expérience consiste à générer un certain nombre de points et à calculer une approximation de la valeur de pi à partir de ces points.*/
/* Ensuite, la fonction calcule la moyenne des résultats de toutes les expériences et détermine l'erreur absolue                              */
/* et relative par rapport à la valeur réelle de pi (M_PI). Enfin, la fonction affiche les résultats.                                         */ /*                                                                                                                                            */
/*                                                                                                                                            */
/*  En Entrée :Tab_exp[ ] un tableau de double pour stocker les résultats de chaque expérience de simulation de Pi                            */
/*             et nb_exp un entiers qui represente le nombres  d'expériences a realiser                                                       */
/*                                                                                                                                            */
/*                                                                                                                                            */
/*                                                                                                                                            */
/* En sortie:  La fonction ne renvoie rien, Elle affiche simplement les résultats des expériences                                             */
/*                                                                                                                                            */
/*                                                                                                                                            */
/***********************************************************************************************************************************************/

/*question 2*/

void Experience_moyenne(double Tab_exp[], int nb_exp){

// Nombre de points à utiliser dans chaque expérience
long nbPoints[3] = {1000, 1000000,1000000000};


// Boucle pour réaliser toutes les expériences
for(int i=0; i<nb_exp; i++){

    // Appel de la fonction simPi avec le nombre de points spécifié
    Tab_exp[i] = simPi(nbPoints[1]);
}

// Calcul de la moyenne des résultats de toutes les expériences
double pi_means = 0;
  
for(int i=0; i<nb_exp; i++) {
  
    pi_means += Tab_exp[i];
}
  
pi_means /= nb_exp;

// Calcul de l'erreur absolue
double erreur_absolue = fabs(pi_means - M_PI);

// Calcul de l'erreur relative
double erreur_rel = erreur_absolue / M_PI;

// Affichage des résultats
printf("PI moyen = %f\n", pi_means);
  
printf("Erreur absolue = %f\n", erreur_absolue);
  
printf("Erreur relative = %f\n", erreur_rel);

  }

/*question 3*/


/**********************************************************************************************************************************************/
/*                                                                                                                                            */
/* Intervalle _confiance                    Calcul des intervalles de confiance autour de la moyenne simulée                                  */
/*                                                                                                                                            */ 
/*  Notre fonction calcule l'intervalle de confiance à 99% pour un tableau de résultats expérimentaux. Elle commence par calculer la moyenne  */
/*  des résultats de la simulation Pi stocké dans Tab_exp, puis la variance. Ensuite, elle utilise une valeur de t de Student prédéfinie pour */
/*  alpha=0.01 pour calculer la largeur de l'intervalle de confiance. Enfin, elle affiche l'intervalle de confiance calculé.                  */ /*                                                                                                                                            */ /*                                                                                                                                            */
/*                                                                                                                                            */
/*  En Entrée :Tab_exp[ ] un tableau de double qui conient les résultats de chaque expérience de simulation de Pi                             */
/*             et nb_exp un entiers qui represente le nombres  d'expériences realiser                                                         */
/*                                                                                                                                            */
/*                                                                                                                                            */
/*                                                                                                                                            */
/*  En sortie: La fonction ne renvoie rien, elle affiche l'intervalle de confiance sous la forme d'une paire de valeurs numériques [X_bar - R,*/ 
/*    X_bar+R], où X_bar - R est la limite inférieure de l'intervalle et X_bar+R est la limite supérieure de l'intervalle.                    */ /*                                                                                                                                            */ 
/*                                                                                                                                            */
/*                                                                                                                                            */
/*                                                                                                                                            */
/***********************************************************************************************************************************************/

void intervalle_confiance(double Tab_exp[], int nb_exp) {
    // Calcul de la moyenne des résultats de toutes les expériences
    double X_bar = 0;
    for(int i=0; i<nb_exp; i++) {
        X_bar += Tab_exp[i];
    }
    X_bar /= nb_exp;

    // Calcul de la variance
    double s = 0;
    for(int i=0; i<nb_exp; i++) {
        s += pow(Tab_exp[i] - X_bar, 2);
    }
    s = (s / (nb_exp - 1));

    // valeurs attribuées à t pour alpha=0.01 pour n = 10, 20, 30 et 40
  
   // double t = 3.169;  // pour n= 10  
  
    //double t=2.845;  // pour n= 20 
  
    double t =2.750; // pour n= 30  
  
    //double t= 2.704;  // pour n= 40  

    // Calcul de l'intervalle de confiance
    
    double R = t * sqrt(s / nb_exp);
    double IC_lower = X_bar - R;
    double IC_upper = X_bar + R;
        // Affichage des résultats
    printf("Intervalle de confiance : [%f, %f]\n",  IC_lower, IC_upper);
    
}
