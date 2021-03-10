#include<stdio.h>
#include<string.h>
#include<math.h>
#include<unistd.h>
#include<stdlib.h>
#define BUFSIZE 1000

float typeofwave(float Lw, float h);
float wavelenght(float *Tw, float h, float lat);
float planetarylocalgravity(float lat);
int read_ints();

int velocitiesxd(float *ampa, float *dx, float *dz, float *dt, float *Twa, float Lwa, float h, int in);
int velocitiesxt(float *ampa, float *dx, float *dz, float *dt, float *Twa, float Lwa, float h, int in);


float pi=3.14159;
float conv=1;

int main(int argc, char *argv[])
{
float Lwa2=0, g=9.81, a=0, m=0.02615, mt=0, xt=0,yt=0, op=0, dif=0, amp=0, lat=0;
int condition=0, i=0, t=0, h=0, Lw=0;

//we read the size of the data file to store infromation on the arrays
t=read_ints();
t=t/5;


//Here we create  the dynamic arrays to store the values from the data file
float *Twa = malloc(t * sizeof(*Twa));
if (!Twa) {
    printf("There was a problem with malloc.");
    exit(EXIT_FAILURE);
}
float *ampa = malloc(t * sizeof(*ampa));
if (!ampa) {
    printf("There was a problem with malloc.");
    exit(EXIT_FAILURE);
}

float *dx = malloc(t * sizeof(*dx));
if (!dx) {
    printf("There was a problem with malloc.");
    exit(EXIT_FAILURE);
}

float *dz = malloc(t * sizeof(*dz));
if (!dz) {
    printf("There was a problem with malloc.");
    exit(EXIT_FAILURE);
}

float *dt = malloc(t * sizeof(*dt));
if (!dt) {
    printf("There was a problem with malloc.");
    exit(EXIT_FAILURE);
}

float *Cwa = malloc(t * sizeof(*Cwa));
if (!Cwa) {
    printf("There was a problem with malloc.");
    exit(EXIT_FAILURE);
}

//Here we store the data from the file to each array as it corresponds
FILE *fpa = fopen(argv[1], "r"); /* "r" = open for reading */
char buff[BUFSIZE]; /* a buffer to hold what you read in */

int count=1, ip=0, ia=0, ix=0, iz=0, it=0, condit=1;

while(!feof (fpa))
  {
    if(condit==0)
    {
      count=1;
      condit=1;
    }

    if (count==1)
    {
      fscanf(fpa, "%f", &Twa[ip]);

      ip++;
    }
    else if (count==2)
    {
      fscanf(fpa, "%f", &ampa[ia]);

      ia++;
    }

    else if (count==3)
    {
      fscanf(fpa, "%f", &dx[ix]);

      ix++;
    }

    else if (count==4)
    {
      fscanf(fpa, "%f", &dz[iz]);

      iz++;
    }

    else if (count==5)
    {
      fscanf(fpa, "%f", &dt[it]);

      it++;
      condit=0;
    }
    count++;


  }
fclose(fpa);

ia=0;






//We introduce some paameters bu hand that are defined by the user
printf("Give the wave parameters to calculate if the wave is on deep water, transitional waters or shallow waters\n");
printf("Introduce the average depth at that part of the sea\n");
scanf("%d", &h);
printf("\n" );
printf("Give the latitude of the buoy: \t");
scanf("%f",&lat);
printf("\n" );



while (i<=t)
{
  //The non-dimensional parameters to calculate the wave order are calculated
  xt=(h/(g*pow(Twa[i],2)));
  yt=(ampa[i]/(g*pow(Twa[i],2)));
  mt=(yt-.00005)/(xt-.001369);

//Here we define if our wave is linear or 2nd order
  if (mt<m||((yt<.001)&&(xt>0.37)))
    {
      op=1;
    }
    else
    {
      op=2;
    }

Lwa2=wavelenght(Twa+i,h,lat);
Cwa[i]=typeofwave(Lwa2,h);


if (Cwa[i]==1&&op==1)

{
    //Case for linear deep water waves

    Lwa2=1.56*pow(Twa[i],2);
    velocitiesxd(ampa+i, dx+i, dz+i, dt+i, Twa+i, Lwa2, h, i);



}
else if (Cwa[i]==3&&op==1)
{
  dif=1;
  Lwa2=1.56*pow(Twa[i],2);
  while(dif>0.1)
  {
  Lw=1.56*pow(Twa[i],2)*tanh((2*pi*h)/Lwa2);
  Lwa2=Lw;
  dif=sqrt(pow(pow(Lw,2)-pow(Lwa2,2),2));
  }
  velocitiesxt(ampa+i, dx+i, dz+i, dt+i, Twa+i, Lwa2, h, i);

}

else if(Cwa[i]==2&&op==1)
{
    //Case for linear transitional water waves
}

else if (op==2)
{
  //Case for second order waves

}
i++;
op=0;
}





return 0;
}

float wavelenght(float *Tw, float h, float lat)

{
    float gravitylocal=0, wavelenghtval=0;
    //We calculaye the apprximated gravity to be used on the wavelenght calculation
    gravitylocal=planetarylocalgravity(lat);
    //We calculate the wavelenght of the wave
    wavelenghtval=((gravitylocal*pow(*Tw,2))/(2*pi))* pow( tanh( pow((2*pi)*((sqrt(h/gravitylocal))/(*Tw)),(3/2)) ), (2/3) );
    return wavelenghtval;

}

float typeofwave(float Lw, float h)

{
    float ratio=0, val=0;
    //Here we define if the wave is propagating in deep, transitional or shallow waters
    ratio=Lw/h;

    if (ratio<=2)
    {
        val=1;

    }
    else if (ratio>=20)
    {
        val=2;
    }
    else if (ratio<=20&&ratio>=2)
    {
        val=3;
    }

    return val;
}


float planetarylocalgravity(float lat)

{
    float  g=0;
    g=( 9.780327 * (1+(0.0053024*(pow(sin(lat),2)))-(0.0000058*(pow(sin(lat),2)))) );
    return g;
}


int read_ints (void)
{
   FILE *fp;
  fp=fopen("data.txt", "r");
  float i = 0, count=1;

  fscanf(fp, "%f", &i);
  while(!feof (fp))
    {
      fscanf(fp, "%f", &i);
      count++;
    }
  fclose(fp);

  return count;
}


int velocitiesxd(float *ampa, float *dx, float *dz, float *dt, float *Twa, float Lwa, float h, int in)
{

  FILE * fp;

  char name[FILENAME_MAX];
    float z=0, x=0, t=0, k=0, tpi=0, a0=0, a1=0, a2=0, arg=0, kons=0, kons2=0, ratio=0, trig=0;
    int lenght=0, i=0, i2=0;
      snprintf(name, sizeof(name), "%d.txt", in);
      fp = fopen(name, "w");

   lenght=Lwa/(*dx);
   kons2=pi/2;
   //We initialize the arrays dimensions to store the information on the wave velocities
    float *arrayx = malloc(lenght * sizeof(*arrayx));
    if (!arrayx) {
        printf("There was a problem with malloc.");
        exit(EXIT_FAILURE);
    }
    float *arrayy = malloc(lenght * sizeof(*arrayy));
    if (!arrayy) {
        printf("There was a problem with malloc.");
        exit(EXIT_FAILURE);
    }
    //Constants to be used on the functions to calculate the velocities
    tpi=2*pi;
    a0=((tpi*(*ampa))/(*Twa));
    k=(tpi)/Lwa;
    ratio=Lwa/2;
      //We start to calculate the velocities from t0 to tn=wave period
        while(t<(*Twa))
            {
              a2=(t/(*Twa));
              //We calculate the velocity from the mean water level z to the depth of propaation h
                while(z<h)
                    {
                      kons=a0*exp(-k*z);
                      i=0;
                      //We calculate the velocity from the position x=0 to a xn=waves wavelenght
                        while (x<Lwa)
                            {
                                  if (z>ratio)
                                  {
                                    //If the wave does not reach the bottom then its velocity field is 0
                                    arrayx[i]=0;
                                    arrayy[i]=0;


                                  }
                                else
                                  {
                                    //if the wave field reach the bottom its velocity its calculated
                                    arg=tpi*( -(x/Lwa) +a2 );
                                    arrayx[i]=kons*sin(arg);
                                    arrayy[i]= kons*cos(arg);
                                  }
                                  i++;
                                  x=x+(*dx);
                                }
                                while(i2<=i){
                                  //We store the whole data from the arrays to a file
                                  fprintf(fp, "%.3f,%.3f\t", arrayx[i2], arrayy[i2]);
                                  i2++;
                                }
                                fprintf(fp, "\n");
                        i2=0;
                        x=0;
                        z=z+(*dz);
                    }
                    fprintf(fp, "\n\n");
                    z=0;
                    t=t+(*dt);
            }
              fp=NULL;
              in++;


        return 0;
}


int velocitiesxt(float *ampa, float *dx, float *dz, float *dt, float *Twa, float Lwa, float h, int in)
{

  FILE * fp;

  char name[FILENAME_MAX];
    float z=0, x=0, t=0, k=0, tpi=0, a0=0, a1=0, a2=0, arg=0, kons=0, kons2=0, ratio=0, trig=0;
    int lenght=0, i=0, i2=0;
      snprintf(name, sizeof(name), "%dt.txt", in);
      fp = fopen(name, "w");

   lenght=Lwa/(*dx);
   kons2=pi/2;
   //We initialize the arrays dimensions to store the information on the wave velocities
    float *arrayx = malloc(lenght * sizeof(*arrayx));
    if (!arrayx) {
        printf("There was a problem with malloc.");
        exit(EXIT_FAILURE);
    }
    float *arrayy = malloc(lenght * sizeof(*arrayy));
    if (!arrayy) {
        printf("There was a problem with malloc.");
        exit(EXIT_FAILURE);
    }
    //Constants to be used on the functions to calculate the velocities
    tpi=2*pi;
    a0=((tpi*(*ampa))/(*Twa));
    k=(tpi)/Lwa;
      //We start to calculate the velocities from t0 to tn=wave period
        while(t<(*Twa))
            {
              a2=(t/(*Twa));
              //We calculate the velocity from the mean water level z to the depth of propaation h
                while(z<h)
                    {
                      kons=a0*cosh(k*(z-h))/sinh(k*h);
                      i=0;
                      //We calculate the velocity from the position x=0 to a xn=waves wavelenght
                        while (x<Lwa)
                            {

                                    //if the wave field reach the bottom its velocity its calculated
                                    arg=tpi*( -(x/Lwa) +a2 );
                                    arrayx[i]=kons*sin(arg);
                                    arrayy[i]= kons*cos(arg);

                                  i++;
                                  x=x+(*dx);
                                }
                                while(i2<=i){
                                  //We store the whole data from the arrays to a file
                                  fprintf(fp, "%.3f,%.3f\t", arrayx[i2], arrayy[i2]);
                                  i2++;
                                }
                                fprintf(fp, "\n");
                        i2=0;
                        x=0;
                        z=z+(*dz);
                    }
                    fprintf(fp, "\n\n");
                    z=0;
                    t=t+(*dt);
            }
              fp=NULL;
              in++;


        return 0;



}

