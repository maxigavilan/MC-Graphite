///-------------------------------------------------------------------------///
///----------------------Grand Canonical MC for Li-C -----------------------///
///-------------------------------------------------------------------------///

///Developed by: Maximiliano Gavilán
///For its use please contact: maxigavilan@gmail.com

///The code builds a network of intercalation sites, taking as the lowest energy site the one
///that occupies the center of the hexagonal carbon rings and half the distance between two
///graphite sheets. It outputs an intercalation isotherm, the energy of the system, the
///molar partial entropy, the molar partial enthalpy, calculated with fluctuations.
///It is assumed that the network does not change with the SoC.

///LIBRARIES
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <time.h>

///DEFINO CONSTANTES
const int NPUNTOS=50;           ///Isotherm points
const double NPTS=NPUNTOS;      ///"double"
const double EQUIL=10000;        ///Equilibration points step
const int MCS=20000;            ///Processing points step
const double MCSR=MCS;
double mui=-0.15;               ///Initial chemical potential, eV
double muf=-0.04;                ///Final chemical potential, eV
double dmu=(muf-mui)/NPUNTOS;
const int NMUESTRA=1;           ///Number of samples to promediate

///PARAMETROS PARA LA CONSTRUCCION DE LA RED
const int Nlx=12;                            ///Number of sites in X-axis
const int NLY=Nlx;                            ///Number of sites in Y-axis
const int Nly=NLY/2;
const double ladox=Nlx;
const double ladoy=Nly;
const int Nlam=4;                            ///Number of graphite sheets
const double lam=Nlam;
const int Npt=Nlx*NLY*Nlam;                  ///Total number of intercalation sites
const double Nt=Npt;
const double dcx=sqrt((1.42*1.42)-(0.71*0.71));
const double  dcy=2.13;                     //0.86602540378443864676372317075294;//1.0;//sqrt((1.42*1.42)-(0.71*0.71));
const double  dvx=1.0;  		            //Distancia vecino en x
const double  dvz=3.35;
const double  Lx=ladox*2*dcx;  		        ///Box-size in x
const double  Ly=ladoy*2*dcy;  		        ///Box-size in y
const double  Lz=lam*dvz;                   ///Box-size in z
const int  plano=Nlx*Nlam;
const int  Nplanos=2*Nly;

///Parametros
const double  BK=0.00008617385; 		    ///Constante de Boltzmann en eV/K
const double T=296.0;                       ///K
double       kT=BK*T;
const double pi=acos(-1.0);

///Hamiltonian parameters
const double  Eps=BK*296;
const double  Kapp=10*BK*296;
const double  rb=1.42;
const double  rm=4.26;
const double  rn=2.47;
const double  alfa=4.0;
const double  gamm=-1.176*BK*296;
const double  rcortexy=10.0;                    ///cut-off in xy plane
const double  rcortez=3.36;                     ///cut-off in z
const int     vxy=60;                           ///neighbors in the same plane
const int     vz=61*(Nlam-2);                   ///neighbors in different planes

///Out-put files
#define grabavmd_en "vmd.xyz"               ///.xyz file for VMD programm
#define grabared_en "sitios.xyz"            ///Coordinates of the intercalation sites
#define Grabacion_en "datosUVT-3.dat"       ///Data

FILE  *archivo, *archivo1, *archivo2, *archivo3, *archivo7;

///For random number generator
#define ULONGMAX	4294967298.0			// unsigned long max + 1
#define _for(a,b,c) 	for(a=b; a<(c); ++a)		// para simplificar
unsigned int seed[256];					// cosas del generador pseudorandom
unsigned int r;
unsigned char irr;

///Variables
int i,j,k,Kk,TT,M,Ocup[Npt],cont2,Nconst,Npart,ii,jj,J,I,kk,MC,n,n2,numven[Npt][vxy],numven2[Npt][vz];
double l,dx,dy,dz,CORD[Npt][3],k1,k2,k3,k4,k5,k6,cont,R,Ft,s1,ss2,Np,H1[Npt][vxy],H2[Npt][vz];
double Nparticula,HH,tit,HHtot,titprom,HHsitio,NN,HH2,NN2,cum,HN,Ental,Entrop,dxdV;

///Sub-codes declaration
void Grafito();                 ///Graphite intercalation sites generation
void Vecinos();                 ///Neighbors and energies
void Inicializa_generador();    ///Random number generator (initialization)
void vmd();                     ///To generate a VMD movie
void Metropolis();              ///Metropolis algorithm
void escribirDatos();           ///To write the output data

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///	GENERA UN NUMERO REAL PSEUDO-ALEATORIO UNIFORMEMENTE DISTRIBUIDO EN EL INTERVALO [0,1)
/////////////////////////////////////////////////////////////////////////////////////////////////////////

inline double randomm(void)
{
	return (double)(r=seed[irr++]+=seed[r>>24])/ULONGMAX;
}

int main(){
int dd;

double a;

Grafito();  ///Construction of the lattice
Vecinos();  ///Neighbors and energies
Inicializa_generador(); ///random number initialization

///Out-put file titles
(archivo =fopen (Grabacion_en,"a"));
fprintf(archivo,"ChemicalPotential[eV] <SoC> <N> <N>^2 <N^2> <Energy>[eV] <E*N>[eV] dxdV Entalhy[eV] Entropy[eV/K] Samples\n");
fclose(archivo);

printf("****MC begging****\n\n");

for(MC=0;MC<NMUESTRA;MC++){
   for(TT=0;TT<NPUNTOS;TT++){
        tit=NN=NN2=HH2=HN=Ental=Entrop=HHtot=dxdV=0.0;
        ///Erase previous "story" of the lattice
        HH=0.0;
        Nparticula=0;
        for(i=0;i<Npt;i++){Ocup[i]=0;}

        printf("Punto=%d, mu=%f, ",(int)(TT),(float)(mui));

        ///EQUILIBRATION STEP
        for(ii=0;ii<EQUIL;ii++){Metropolis();}

        ///PROCESSING STEP
        for(ii=0;ii<MCS;ii++){Metropolis();
            tit+=3*(Nparticula/Npt);
            NN+=Nparticula;
            NN2+=Nparticula*Nparticula;
            HHtot+=HH;
            HH2+=HH*HH;
            HN+=HH*Nparticula;
        }

        ///Average out-put data
        NN/=MCSR;                   ///Number of particles
        NN2/=MCSR;                  ///Valor medio del num de particulas al cudrado
        tit/=MCSR;                  ///Occupation (State-of-Charge)
        HHtot/=MCSR;                ///Valor medio de la energia
        HH2/=MCSR;                  ///Valor medio de la energia al cuadrado
        HN/=MCSR;                   ///Valor medio del producto de la energia por N
        dxdV=NN2-NN*NN;             ///derivada de la ocupacion respecto al num de particulas
        Ental=(HN-HHtot*NN)/(dxdV); ///Molar enthalpy change
        Entrop=(Ental-mui)/T;       ///Molar entropy change

        printf("x=%f\n",(float)(tit));

        ///Write
        escribirDatos();
        vmd();

        mui+=dmu;
  }///De puntos
}///DE MUESTRAS

}
///--------------------END OF THE CODE----------------------------



///-----------------------SUBPROGRAMAS-----------------------------

void Metropolis(){
int MM,q,qj;
double Es,Es2,dEn,prom,prom2;

    for(MM=0;MM<Npt;MM++){
        Es=Es2=dEn=prom=prom2=0.0;
        q=randomm()*Npt;    ///Elijo un sitio al azar, con random() llamo al generador de numeros aleaotorios

        if(Ocup[q]==0){     ///condicion: si el sitio esta vacio
            for(n=0;n<vxy;n++){qj=numven[q][n];if(Ocup[qj]==1){Es+=H1[q][n];}}
            for(n=0;n<vz;n++){qj=numven2[q][n];if(Ocup[qj]==1){Es+=H2[q][n];}}
            Es+=gamm;

            dEn=Es;
            prom=exp(-((dEn-mui)/kT));
            prom2=randomm();    ///numero aleatorio entre [0,1)
            if(prom2<=prom){ Ocup[q]=1; HH+=dEn; Nparticula++;}
        }
        else{               ///por el contrario: si el sitio esta ocupado
            for(n=0;n<vxy;n++){qj=numven[q][n];if(Ocup[qj]==1){Es2+=H1[q][n];}}
            for(n=0;n<vz;n++){qj=numven2[q][n];if(Ocup[qj]==1){Es2+=H2[q][n];}}
            Es2+=gamm;

            dEn=-Es2;
            prom=exp(-((dEn+mui)/kT));
            prom2=randomm();    ///numero aleatorio entre [0,1)
            if(prom2<=prom){ Ocup[q]=0; HH+=dEn; Nparticula--;}
        }
    }
}

void Vecinos(){
   int n,n2;
   double RR,rr2;
    ///numven[i][n] is a matrix labeling the neighbors j of site i in the same plane
    ///numven2[i][n2] is a matrix labeling the neighbors j of site i in different planes
    ///H1[i][n] is a matrix that stores the in-plane energy of site i with its neighbors j
    ///H2[i][n2] is a matrix that stores the out-plane energy of site i with its neighbors j

    ///numven[i][n] is linked to H1[i][n] by the number n
    ///numven[i][n2] is linked to H2[i][n2] by the number n2

    for(i=0;i<Npt;i++){for(j=0;j<vxy;j++){numven[i][j]=0;H1[i][j]=0.0;}}
    for(i=0;i<Npt;i++){for(j=0;j<vz;j++){numven2[i][j]=0;H2[i][j]=0.0;}}

    for(i=0;i<Npt;i++){
        n=n2=0;
        for (j=0;j<Npt;j++){
                ///CPC
                dx=fabs(CORD[i][0]-CORD[j][0]);if(dx>0.5*Lx){dx=Lx-dx;}
                dy=fabs(CORD[i][1]-CORD[j][1]);if(dy>0.5*Ly){dy=Ly-dy;}
                dz=fabs(CORD[i][2]-CORD[j][2]);if(dz>0.5*Lz){dz=Lz-dz;}

                ///distancia
                RR=sqrt((dx*dx)+(dy*dy)+(dz*dz));
                rr2=sqrt((dx*dx)+(dy*dy));

                ///IN-PLANE (Lennard-Jones)
                if(dz==0.0){if(i!=j){if(rr2<=rcortexy){numven[i][n]=j;
                                                       H1[i][n]=Eps*((pow((rm/RR), 12.0))-2.0*(pow((rm/RR), 6.0)));
                                                        n++;
                                                       }}}
                ///OUT-PLANE (De-Rosa)
                else{if(i!=j){if(rr2<=rcortexy){if(dz<rcortez){numven2[i][n2]=j;
                                                               H2[i][n2]=Kapp*pow(rb/RR, alfa);
                                                               n2++;
                                                               }}}}
        } //del for j
    }//del for i

}///END NEIGHBORS

void escribirDatos(){
if((archivo =fopen (Grabacion_en,"a"))== NULL );
fprintf(archivo,"%4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f\n", (float)(mui), (float)(tit),(float)(NN),(float)(NN*NN),(float)(NN2),
        (float)(HHtot), (float)(HN),(float)(dxdV),(float)(Ental),(float)(Entrop),(float)(NMUESTRA));
fclose(archivo);
}/////

void vmd(){
    (archivo1 =fopen (grabavmd_en,"a"));
   fprintf(archivo1,"%d \n", (int)(Npt));  //\n graba uno por linea
    fprintf(archivo1,"\n");
    for(i=0;i<Npt;i++){
        if(Ocup[i]==1){
               fprintf(archivo1,"O %4.5f %4.5f %4.5f\n",(float)(CORD[i][0]),(float)(CORD[i][1]),(float)(CORD[i][2]));
        }


        else{fprintf(archivo1,"O %4.5f %4.5f -1.0\n",(float)(Lx/2),(float)(Ly/2));}
    }
    fclose(archivo1);
}

void Grafito(){
    /// Coordenadas en x
    for(i=0;i<Nplanos;i++){for(j=0;j<Nlam;j++){for(k=0;k<Nlx;k++){if(i%2==0){CORD[k+j*Nlx+i*plano][0]=k*2*dcx;}else{CORD[k+j*Nlx+i*plano][0]=dcx+k*2*dcx;}}}}
    /// Coordenadas en y
    for(i=0;i<Nplanos;i++){for(j=0;j<Nlam;j++){for(k=0;k<Nlx;k++){if(i%2==0){CORD[k+j*Nlx+i*plano][1]=i*dcy;}else{CORD[k+j*Nlx+i*plano][1]=i*dcy;}}}}
    /// Coordenadas en z
    for(i=0;i<Nplanos;i++){for(j=0;j<Nlam;j++){for(k=0;k<Nlx;k++){if(i%2==0){CORD[k+j*Nlx+i*plano][2]=j*dvz;}else{CORD[k+j*Nlx+i*plano][2]=j*dvz;}}}}

   ///Grabo coordenadas de la red
    if((archivo2=fopen (grabared_en,"a"))== NULL );
        fprintf(archivo2,"%d  \n \n", (int)(Npt));  //\n graba uno por linea
        for(i=0;i<Npt;i++){fprintf(archivo2,"Li %4.5f %4.5f  %4.5f \n",  (float)(CORD[i][0]), (float)(CORD[i][1]), (float)(CORD[i][2]));}
        fclose(archivo2);
}

void Inicializa_generador(void)
{
	int i,j,s;
// inicializo generador pseudorandom: semilla = segundos desde 1970
	srand((unsigned)time(0));
	irr=1;
	_for(i,0,256) seed[i]=rand();
	r=seed[0];
	_for(i,0,70000) r=seed[irr++]+=seed[r>>24];

}


