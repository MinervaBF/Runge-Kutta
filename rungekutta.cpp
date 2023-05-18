#include <cmath>
#include <math.h>
#include <fstream>
#define delta 6.67E-11*5.9736E24/pow(3.844E8, 3)
#define mu 0.07349E24/5.9736E24
#define omega 2.6617E-6
#define disttierraluna 3.844E8
#define g 6.67E-11
#define masaluna 0.07349E24
#define masatierra 5.9736E24

double rpunto(double p);
double phipunto(double p, double r);
double prpunto(double p, double r, double f, double t);
double pphipunto(double r, double f, double t);

int main()
{
    double r, phi, theta, v, pr, pphi, r0, phi0, pr0, pphi0, rprima;
    double posluna[2], posnave[2];
    double k1[4], k2[4], k3[4], k4[4], h, t, hamilton;
    int i, pasos, m;
    FILE *f1, *f2;


    //Abrimos los archivos
    f1=fopen("naveyluna.txt","w");
    f2=fopen("conservada.txt","w");

    //Escribimos el tiempo durante el que se representará la animación:
    pasos=7*24*3600; //Pasamos una semana a segundos.

    //Introducimos las condiciones iniciales reescaladas.
    r=6.37816E6/(3.844E8);

    //La nave parte de la superficie terrestre con velocidad v0 y ángulo theta0 (se cambiarán a voluntad):
    v=(11.2/(3.844E8))*1000; //Reescalamos, asumiendo que se da en km/s (11.2 km/s es la velocidad de escape).
    phi=M_PI/2.; //Ángulo desde el que parte el cohete.
    theta=44*M_PI/180; //Ángulo de despegue.

    //Calculamos los valores iniciales de los momentos pr y pphi:
    pr=v*cos(theta-phi);
    pphi=r*v*sin(theta-phi);

    //Damos una masa aleatoria al cohete para calcular posteriormente el Hamiltoniano.
    m=5;

    //Calculamos las posiciones iniciales de la Luna y del cohete (reescaladas):
    posluna[0]=1.;
    posluna[1]=0.;
    posnave[0]=r*cos(phi);
    posnave[1]=r*sin(phi);

    fprintf(f1, "%e%c\t%e\n", posnave[0], 44, posnave[1]);
    fprintf(f1, "%e%c\t%e\n", posluna[0], 44, posluna[1]);
    fprintf(f1, "%i%c\t%i\n", 0, 44, 0);
    fprintf(f1, "\n");

    h=60.;
    i=0;

    //Comenzamos el ciclo:

    for(t=60; t<=pasos; t=t+h)
    {
        //Calculamos k_1 (empleando las ecuaciones de movimiento tras el reescalado):
        k1[0]=h*rpunto(pr);
        k1[1]=h*phipunto(pphi, r);
        k1[2]=h*prpunto(pphi, r, phi, t);
        k1[3]=h*pphipunto(r, phi, t);

        //Calculamos k_2:
        k2[0]=h*rpunto((pr+k1[3]/2.));
        k2[1]=h*phipunto((pphi+k1[3]/2.), (r+k1[0]/2.));
        k2[2]=h*prpunto((pphi+k1[3]/2.), (r+k1[0]/2.), (phi+k1[1]/2.), (t+h/2.));
        k2[3]=h*pphipunto((r+k1[0]/2.), (phi+k1[1]/2.), (t+h/2.));

        //Calculamos k_3:
        k3[0]=h*rpunto((pr+k2[3]/2.));
        k3[1]=h*phipunto((pphi+k2[3]/2.), (r+k2[0]/2.));
        k3[2]=h*prpunto((pphi+k2[3]/2.), (r+k2[0]/2.), (phi+k2[1]/2.), (t+h/2.));
        k3[3]=h*pphipunto((r+k2[0]/2.), (phi+k2[1]/2.), (t+h/2.));

        //Calculamos k_4:
        k4[0]=h*rpunto((pr+k3[3]));
        k4[1]=h*phipunto((pphi+k3[3]), (r+k3[0]));
        k4[2]=h*prpunto((pphi+k3[3]), (r+k3[0]), (phi+k3[1]), (t+h));
        k4[3]=h*pphipunto((r+k3[0]), (phi+k3[1]), (t+h));

        //Calculamos los nuevos valores de r, phi, pr y pphi:
        r=r+1./6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
        phi=phi+1./6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
        pr=pr+1./6*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
        pphi=pphi+1./6*(k1[3]+2*k2[3]+2*k3[3]+k4[3]);

        //Deshacemos el reescalamiento.
        r0=r*disttierraluna;
        phi0=phi;
        pr0=pr*m*disttierraluna;
        pphi0=pphi*m*pow(disttierraluna,2.);


        //Calculamos la constante de movimiento H'. Hemos dado un valor arbitrario para la masa del cohete, pues solo nos interesa ver que es constante.
        rprima=(pow(r0, 2.)+pow(disttierraluna,2.)-2*r0*disttierraluna*cos(phi0-omega*t), 1./2);
        hamilton=pow(pr0,2.)/(2.*m)+pow(pphi0, 2.)/(2*m*pow(r0, 2.))-g*m*masatierra/r0-g*m*masaluna/rprima-omega*pphi0;
        fprintf(f2, "%e\t%e\n", t, hamilton);

        //Calculamos la posición del cohete y de la Luna.
        posnave[0]=r*cos(phi);
        posnave[1]=r*sin(phi);
        posluna[0]=cos(omega*t);
        posluna[1]=sin(omega*t);

        //Escribimos los datos en un fichero, separados por una coma.
        if(i%10==0)
        {
        fprintf(f1, "%e%c\t%e\n", posnave[0], 44, posnave[1]);
        fprintf(f1, "%e%c\t%e\n", posluna[0], 44, posluna[1]);
        fprintf(f1, "%i%c\t%i\n", 0, 44, 0);
        fprintf(f1, "\n"); //Los datos tienen que estar separados en bloques.
        }
        i++;

        t=t+h;

    }
    fclose(f1);
    fclose(f2);

}

double rpunto(double p)
{
    double r;
    r=p;
    return r;
}

double phipunto(double p, double r)
{
    double phi;
    phi=p/pow(r,2.);
    return phi;
}

double prpunto(double p, double r, double f, double t)
{
    double pr, rprima;
    rprima=pow(1+pow(r,2.)-2*r*cos(f-omega*t),1./2);
    pr=pow(p, 2.)/pow(r, 3.)-delta*(1/pow(r, 2.)+mu/pow(rprima, 3.)*(r-cos(f-omega*t)));
    return pr;
}

double pphipunto(double r, double f, double t)
{
    double pphi, rprima;
    rprima=pow(1+pow(r,2.)-2*r*cos(f-omega*t),1./2);
    pphi=(-1)*delta*mu*r/pow(rprima,3.)*sin(f-omega*t);
    return pphi;
}
