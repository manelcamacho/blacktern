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
int velocities2xd(float *ampa, float *dx, float *dz, float *dt, float *Twa, float Lwa, float h, int in);
int velocities2xt(float *ampa, float *dx, float *dz, float *dt, float *Twa, float Lwa, float h, int in);
int spectral(float *ampa, float *Twa, float Lwa, float h, int in, float lat, float *Cwa);



float pi=3.14159;
float conv=1;
float rho=1023.6;

int main(int argc, char *argv[])
{
float Lwa2=0, g=9.81, a=0, m=0.02615, mt=0, xt=0,yt=0, op=0, dif=0, amp=0, lat=0, xdt=0.084, xst=0.00154, x=0, func1=0, func2=0,func3=0,func4=0, func5=0, func6=0;
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



while (i<t)
{
  //The non-dimensional parameters to calculate the wave order are calculated

  xt=(h/(g*pow(Twa[i],2)));
  yt=((2*ampa[i])/(g*pow(Twa[i],2)));
  x=xt;


  /*Here we set the functions that define each range of te wave theories starting from lineal, second order, and fifth*/
  /*This is the function that divides the Airy space-second stokes space*/
  func1=0.0000322962 - (0.00200355 *(pow(x,(0.25)))) + (0.0111283 *pow(x,0.5)) - (0.0136872*x);

  /*This is the function that divides the Kortewegg de Vries space-second stokes and Fifth-Third stokes space*/
  func2=(0.00124726 - (0.0244293*pow(x,0.25)) + (0.0988321*pow(x,0.5)) - (0.10934*x));

  /*This is the function that divides Fifth-Forth-Third stokes space to breaking waves*/
  func3=(0.00472213  - (0.0966994*pow(x,0.25)) + (0.390791*pow(x,0.5)) - (0.432903*x));

  /*This is the function that divides the Kortewegg de Vries space from second order and the Stokes-Kortewegg de Vries with 3rd and 4th order*/
  func4=(-0.0262843 + (0.2938879*pow(x,0.25)) - (0.92127 *pow(x,0.5)) + (2.77775*x));

  /*This is the function that divides the Third-stream-fifth stokes spaces*/
  func5=(-0.0242947 + (0.114039*pow(x,0.25)) - (0.0975339*pow(x,0.5)) + (0.00202345*x));

  /*This is the function that divides the Third-Forth stokes spaces*/
  func6=(-0.122536 + (0.608611*pow(x,0.25)) - (0.725805 *pow(x,0.5)) + (0.299734*x));

  printf("xt:%f, yt:%f, x%f, y:%f\n", xt, yt, x, func1);



  Lwa2=wavelenght(Twa+i,h,lat);
  Cwa[i]=typeofwave(Lwa2,h);

  if(xt<0.1){

  /*Here we compare our non dimensional values to confirm if we hace a linear wave*/
    if(func1>=yt)
    {


        if (xt>xdt) {
          Lwa2=1.56*pow(Twa[i],2);
          velocitiesxd(ampa+i, dx+i, dz+i, dt+i, Twa+i, Lwa2, h, i);
          Cwa[i]=11;
        }
        else if (xt<=xdt) {

      dif=1;
     Lwa2=1.56*pow(Twa[i],2);
     while(dif>0.01)
      {
      Lw=1.56*pow(Twa[i],2)*tanh((2*pi*h)/Lwa2);
      Lwa2=Lw;
      dif=sqrt(pow(pow(Lw,2)-pow(Lwa2,2),2));
      }
      velocitiesxt(ampa+i, dx+i, dz+i, dt+i, Twa+i, Lwa2, h, i);
      Cwa[i]=12;

        }
        else {
          printf("No wave recognized\n" );
        }

    }
    /*In this part we compare to see if the wave is second order or cnoidal*/


    else if (func1<yt&&func2>yt) {
        /*we compare the slope that divides the cnoidal theory from the second order theory*/
        if (func4>yt) {

          if ( xt>=xdt)
          {
              Lwa2=1.56*pow(Twa[i],2);
              velocities2xd(ampa+i, dx+i, dz+i, dt+i, Twa+i, Lwa2, h, i);
              Cwa[i]=21;
            }
            else if (xt<xdt)
            {
              dif=1;
              Lwa2=1.56*pow(Twa[i],2);
              while(dif>0.01)
              {
              Lw=1.56*pow(Twa[i],2)*tanh((2*pi*h)/Lwa2);
              Lwa2=Lw;
              dif=sqrt(pow(pow(Lw,2)-pow(Lwa2,2),2));
              }
              velocities2xt(ampa+i, dx+i, dz+i, dt+i, Twa+i, Lwa2, h, i);
              Cwa[i]=22;
            }

        }
        else if (func4<=yt) {
          /* code */printf("The wave is a cnoidal wave\n");
          Cwa[i]=6;
        }

    }

    /*Here we search for a 5th, 4th or 3rd order wave*/
    else if (func2<yt&&func3>yt) {
      /*Here we compare for a 5th order wave*/
      if (func4<yt){
        printf("You have a 5th order wave\n" );
        Cwa[i]=5;
      }
        /*Here we compare for a 4th or 3rd order wave*/
      else if (func4>yt) {
        if (func6>=yt){
        printf("You have a third order wave \n" );
        Cwa[i]=3;
        }
        else if (func6<yt) {
          printf("You have a fourth order wave\n" );
          Cwa[i]=4;
        }
      }
    }
      else
      {

        printf("You have shallow waters waves or breaking conditions\n" );
        Cwa[i]=0;
      }

    }

    else if(xt>=0.1)
    {
      if(yt<0.001)
      {
        Lwa2=1.56*pow(Twa[i],2);
        velocitiesxd(ampa+i, dx+i, dz+i, dt+i, Twa+i, Lwa2, h, i);
        Cwa[i]=11;

      }
      else if(0.001<=yt&&yt<0.0073)
      {
        Lwa2=1.56*pow(Twa[i],2);
        velocities2xd(ampa+i, dx+i, dz+i, dt+i, Twa+i, Lwa2, h, i);
        Cwa[i]=21;
      }
      else if(0.0073<=yt&&yt<0.0088)
      {
      /* code */printf("The wave is a 5th order order wave\n");
      Cwa[i]=5;
      }
      else if(0.0088<=yt&&yt<0.0198)
      {
      /* code */printf("The wave is a 3rd order order wave\n");
      Cwa[i]=3;
      }
      else if(0.0198<=yt&&yt<0.0285)
      {
      /* code */printf("The wave is a 4th order order wave\n");
      Cwa[i]=4;
      }
      else{
        /* code */printf("breaking wave or not recgnize pattern\n");
        Cwa[i]=0;
      }

    }

i++;
op=0;
}

i=0;

  while(i<t)
  {
    Lwa2=wavelenght(Twa+i,h,lat);
    if(Lwa2>(h/2))
    {
      dif=1;
      Lwa2=1.56*pow(Twa[i],2);
      while(dif>0.1)
      {
        Lw=1.56*pow(Twa[i],2)*tanh((2*pi*h)/Lwa2);
        Lwa2=Lw;
        dif=sqrt(pow(pow(Lw,2)-pow(Lwa2,2),2));
      }
    }
    spectral(ampa+i, Twa+i, Lwa2, h, i, lat, Cwa+i);
    i++;
  }




return 0;
}

int spectral(float *ampa, float *Twa, float Lwa, float h, int in, float lat, float *Cwa)
{
  FILE * fp;
  char name[FILENAME_MAX];
    float g=0, power=0, rho=0, depthoint=0, mvvx=0, mvvy=0, a0=0, f=0, k=0;
      fp = fopen("spectral.txt", "a");
      if(fp == NULL)
   {
      printf("Error could not open or create file!");
      exit(1);
   }
      g=planetarylocalgravity(lat);
      f=(1/(*Twa));
      power=((rho*g*g)/(64*pi))*4*(*ampa)*(*ampa)*(*Twa);
      depthoint=Lwa/2;
      a0=(*ampa*2*pi)/(*Twa);
      k=(2*pi)/Lwa;
      if((*Cwa)==11)
      {
        mvvx=a0;
        mvvy=a0;
      }
      else if((*Cwa)==12)
      {
        mvvx=a0*(cosh(k*(-h)))/sinh(k*h);
        mvvy=a0*(sinh(k*(-h)))/sinh(k*h);
      }
      else if((*Cwa)==21)
      {
        mvvx=fabs((a0*(-1/4)*(1/sinh(k*h)))*((sinh(k*(-h))) + (3*(*ampa)*k*(pow(1/sinh(k*h),3))*(sinh(k*(-h))))));
        mvvy=fabs((a0*(-1/4)*(1/sinh(k*h)))*((cosh(k*(-h))) + (3*(*ampa)*k*(pow(1/sinh(k*h),3))*(cosh(k*(-h))))));

      }
      else if((*Cwa)==22)
      {
        mvvx=fabs((-a0*cosh(k*(-h))*(1/sinh(k*h))));
        mvvy=fabs((-a0*sinh(k*(-h))*(1/sinh(k*h))));
      }
      else{
        mvvx=0;
        mvvy=0;
      }
      if(in==0)
      {
        fprintf(fp, "H(m),L(m),T(s),f(Hz),P(W),DI(m),Mvx(m/s),Mvy(m/s)\n");
      }
      fprintf(fp, "%.3f,%.3f,%.3f,%.3f,%.3f,%3f,%3f,%3f\n", *ampa, Lwa, *Twa, f, power, depthoint, mvvx, mvvy);
  fp=NULL;

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

//Here we calculate the linear velocities and deep waters
int velocitiesxd(float *ampa, float *dx, float *dz, float *dt, float *Twa, float Lwa, float h, int in)
{

  FILE * fp1, * fp2;

  char namex[FILENAME_MAX];
  char namep[FILENAME_MAX];
    float z=0, x=0, t=0, k=0, tpi=0, a0=0, a1=0, a2=0, arg=0, kons=0, kons2=0, ratio=0, n=0, g=9.8;
    int lenght=0, i=0, i2=0;
    printf("You have a linear deep water wave\n");
      snprintf(namex, sizeof(namex), "%dx1.txt", in);
      snprintf(namep, sizeof(namep), "%dp1.txt", in);
      fp1 = fopen(namex, "w");
      fp2 = fopen(namep, "w");

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
    float *arrayp = malloc(lenght * sizeof(*arrayp));
    if (!arrayp) {
        printf("There was a problem with malloc.");
        exit(EXIT_FAILURE);
    }
    //Constants to be used on the functions to calculate the velocities
    *ampa=*ampa/2;
    tpi=2*pi;
    a0=((tpi*(*ampa))/(*Twa));
    k=(tpi)/Lwa;
    ratio=Lwa/2;


    while(x<Lwa){
      fprintf(fp1, "Vx(%d),Vy(%d),",i,i);
      fprintf(fp2, "P(%d),", i);
      i++;
     x=x+(*dx);
    }
    fprintf(fp1, "\n");
    fprintf(fp2, "\n");
    x=0;
    i=0;



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
                                    n=rho*g*(z+((*ampa)*sin(arg)));
                                    arrayp[i]=n;


                                  }
                                else
                                  {
                                    //if the wave field reach the bottom its velocity its calculated
                                    arg=tpi*( -(x/Lwa) +a2 );
                                    arrayx[i]=kons*sin(arg);
                                    arrayy[i]=kons*cos(arg);
                                    n=rho*g*(z+((*ampa)*sin(arg)));
                                    arrayp[i]=n;
                                  }
                                  i++;
                                  x=x+(*dx);
                                }
                                while(i2<i){
                                  //We store the whole data from the arrays to a file
                                  fprintf(fp1, "%.3f,%.3f,",arrayx[i2],arrayy[i2]);
                                  fprintf(fp2, "%.3f,", arrayp[i2]);
                                  i2++;
                                }
                                fprintf(fp1, "\n");
                                fprintf(fp2, "\n");
                        i2=0;
                        x=0;
                        z=z+(*dz);
                    }
                    fprintf(fp1, "\n\n");
                    fprintf(fp2, "\n\n");
                    z=0;
                    t=t+(*dt);
            }
              fp1=NULL;
              fp2=NULL;
              in++;


        return 0;
}

//Here we calculate the linear velocities at transitional water
int velocitiesxt(float *ampa, float *dx, float *dz, float *dt, float *Twa, float Lwa, float h, int in)
{

  FILE * fp1, * fp2;

  char namex[FILENAME_MAX];
  char namep[FILENAME_MAX];
    float z=0, x=0, t=0, k=0, tpi=0, a0=0, a1=0, a2=0, arg=0, kons=0, kons1=0, kons2=0, ratio=0, n=0, g=9.81;
    int lenght=0, i=0, i2=0;
    printf("You have a linear transitional water wave\n");
      snprintf(namex, sizeof(namex), "%dx1.txt", in);
      snprintf(namep, sizeof(namep), "%dp1.txt", in);
      fp1 = fopen(namex, "w");
      fp2 = fopen(namep, "w");

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
    float *arrayp = malloc(lenght * sizeof(*arrayp));
    if (!arrayp) {
        printf("There was a problem with malloc.");
        exit(EXIT_FAILURE);
    }
    //Constants to be used on the functions to calculate the velocities
    *ampa=*ampa/2;
    tpi=2*pi;
    a0=((tpi*(*ampa))/(*Twa));
    k=(tpi)/Lwa;

    while(x<Lwa){
      fprintf(fp1, "Vx(%d),Vy(%d),",i,i);
      fprintf(fp2, "P(%d),", i);
      i++;
     x=x+(*dx);
    }
    fprintf(fp1, "\n");
    fprintf(fp2, "\n");
    x=0;
    i=0;



      //We start to calculate the velocities from t0 to tn=wave period
        while(t<(*Twa))
            {
              a2=(t/(*Twa));
              //We calculate the velocity from the mean water level z to the depth of propaation h
                while(z<h)
                    {
                      kons=a0*cosh(k*(z-h))/sinh(k*h);
                      kons1=a0*sinh(k*(z-h))/sinh(k*h);
                      i=0;
                      //We calculate the velocity from the position x=0 to a xn=waves wavelenght
                        while (x<Lwa)
                            {

                                    //if the wave field reach the bottom its velocity its calculated
                                    arg=tpi*( -(x/Lwa) +a2 );
                                    arrayx[i]=kons*sin(arg);
                                    arrayy[i]=kons1*cos(arg);
                                    n=rho*g*(z+( (*ampa)*sin(arg) ));
                                    arrayp[i]=n;

                                  i++;
                                  x=x+(*dx);
                                }
                                while(i2<i){
                                  //We store the whole data from the arrays to a file
                                  fprintf(fp1, "%.3f,%.3f,",arrayx[i2],arrayy[i2]);
                                  fprintf(fp2, "%.3f,", arrayp[i2]);
                                  i2++;
                                }
                                fprintf(fp1, "\n");
                                fprintf(fp2, "\n");
                        i2=0;
                        x=0;
                        z=z+(*dz);
                    }
                    fprintf(fp1, "\n\n");
                    fprintf(fp2, "\n\n");
                    z=0;
                    t=t+(*dt);
            }
              fp1=NULL;
              fp2=NULL;
              in++;


        return 0;



}



//Here w calculate the 2nd order velocities at transitional waters
int velocities2xt(float *ampa, float *dx, float *dz, float *dt, float *Twa, float Lwa, float h, int in)
{

  FILE * fp1, * fp2;

  char namex[FILENAME_MAX];
  char namep[FILENAME_MAX];
    float z=0, x=0, t=0, k=0, tpi=0, a0=0, a1=0, a2=0, a3=0, arg=0, kons=0, kons1x=0, kons2x=0, kons1y=0, kons2y=0, ratio=0, hk=0, n=0, g=9.81, nam=0;
    int lenght=0, i=0, i2=0;
    printf("You have a 2nd order transitonal wave\n");
    snprintf(namex, sizeof(namex), "%dx2.txt", in);
    snprintf(namep, sizeof(namep), "%dp2.txt", in);
    fp1 = fopen(namex, "w");
    fp2 = fopen(namep, "w");

   lenght=Lwa/(*dx);
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
    float *arrayp = malloc(lenght * sizeof(*arrayp));
    if (!arrayp) {
        printf("There was a problem with malloc.");
        exit(EXIT_FAILURE);
    }
    //Constants to be used on the functions to calculate the velocities, 2pi, a*w, wave number, depth and wavenumber, cte a at velocities, wavelenght ratio
    *ampa=*ampa/2;
    tpi=2*pi;
    a0=((tpi*(*ampa))/(*Twa));
    k=(tpi)/Lwa;
    hk=h*k;
    a1=(-a0*(1/sinh(hk)))/4;
    ratio=Lwa/2;
    a3=*ampa*3*k;
    nam=tanh(hk);

    while(x<Lwa){
      fprintf(fp1, "Vx(%d),Vy(%d),",i,i);
      fprintf(fp2, "P(%d),", i);
      i++;
     x=x+(*dx);
    }
    fprintf(fp1, "\n");
    fprintf(fp2, "\n");
    x=0;
    i=0;

      //We start to calculate the velocities from t0 to tn=wave period
        while(t<(*Twa))
            {
              a2=(t/(*Twa));
              //We calculate the velocity from the mean water level z to the depth of propaation h
                while(z<h)
                    {
                      kons=k*(-h+z);
                      kons1x=cosh(kons);
                      kons2x=cosh(2*kons);
                      kons1y=sinh(kons);
                      kons2y=sinh(2*kons);
                      i=0;
                      //We calculate the velocity from the position x=0 to a xn=waves wavelenght
                        while (x<Lwa)
                            {
                                  if (z>ratio)
                                  {
                                    //If the wave does not reach the bottom then its velocity field is 0
                                    arrayx[i]=0;
                                    arrayy[i]=0;
                                    n=(*ampa)*( (cos(arg)+(((k*(*ampa))*((3-pow(nam,2))/(4*pow(nam,3))))*cos(2*arg))));
                                    arrayp[i]=g*(z+n)*rho;


                                  }
                                else
                                  {

                                    //if the wave field reach the bottom its velocity its calculated
                                    arg=tpi*( -(x/Lwa) +a2 );
                                    arrayx[i]=a1*((4*cos(arg)*kons1x)+(a3*(1/(sinh(hk)*sinh(hk)*sinh(hk)))*cos((2*arg))*kons2x));
                                    arrayy[i]=a1*((4*sin(arg)*kons1y)+(a3*(1/(sinh(hk)*sinh(hk)*sinh(hk)))*sin((2*arg))*kons2y));
                                    n=(*ampa)*( (cos(arg)+(((k*(*ampa))*((3-pow(nam,2))/(4*pow(nam,3))))*cos( (2*arg) ))));
                                    arrayp[i]=g*(z+n)*rho;

                                  }
                                  i++;
                                  x=x+(*dx);
                                }
                                while(i2<i){
                                  //We store the whole data from the arrays to a file
                                  fprintf(fp1, "%.3f,%.3f,",arrayx[i2],arrayy[i2]);
                                  fprintf(fp2, "%.3f,", arrayp[i2]);
                                  i2++;
                                }
                                fprintf(fp1, "\n");
                                fprintf(fp2, "\n");
                        i2=0;
                        x=0;
                        z=z+(*dz);
                    }
                    fprintf(fp1, "\n\n");
                    fprintf(fp2, "\n\n");
                    z=0;
                    t=t+(*dt);
            }
              fp1=NULL;
              fp2=NULL;
              in++;


        return 0;
}

//Here we calculate the 2nd order velocities at deep waters
int velocities2xd(float *ampa, float *dx, float *dz, float *dt, float *Twa, float Lwa, float h, int in)
{

  FILE * fp1, * fp2;

  char namex[FILENAME_MAX];
  char namep[FILENAME_MAX];
    float z=0, x=0, t=0, k=0, tpi=0, a0=0, a1=0, a2=0, arg=0, kons=0, kons1x=0, kons1y=0, ratio=0, hk=0, nam=0, g=9.81;
    int lenght=0, i=0, i2=0, n=0;
    printf("You have a 2nd order depp water wave\n");
    snprintf(namex, sizeof(namex), "%dx2.txt", in);
    snprintf(namep, sizeof(namep), "%dp2.txt", in);
    fp1 = fopen(namex, "w");
    fp2 = fopen(namep, "w");

   lenght=Lwa/(*dx);
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
    float *arrayp = malloc(lenght * sizeof(*arrayp));
    if (!arrayp) {
        printf("There was a problem with malloc.");
        exit(EXIT_FAILURE);
    }
    //Constants to be used on the functions to calculate the velocities, 2pi, a*w, wave number, depth and wavenumber, cte a at velocities, wavelenght ratio
    *ampa=*ampa/2;
    tpi=2*pi;
    a0=((tpi*(*ampa))/(*Twa));
    k=(tpi)/Lwa;
    hk=h*k;
    a1=-a0*(1/sinh(hk));
    ratio=Lwa/2;
    nam=tanh(hk);

    while(x<Lwa){
      fprintf(fp1, "Vx(%d),Vy(%d),",i,i);
      fprintf(fp2, "P(%d),", i);
      i++;
     x=x+(*dx);
    }
    fprintf(fp1, "\n");
    fprintf(fp2, "\n");
    x=0;
    i=0;

      //We start to calculate the velocities from t0 to tn=wave period
        while(t<(*Twa))
            {
              a2=(t/(*Twa));
              //We calculate the velocity from the mean water level z to the depth of propaation h
                while(z<h)
                    {
                      kons=k*(-h+z);
                      kons1x=cosh(kons);
                      kons1y=sinh(kons);
                      i=0;
                      //We calculate the velocity from the position x=0 to a xn=waves wavelenght
                        while (x<Lwa)
                            {
                                  if (z>ratio)
                                  {
                                    //If the wave does not reach the bottom then its velocity field is 0
                                    arrayx[i]=0;
                                    arrayy[i]=0;
                                    n=(*ampa)*( (cos(arg)+(((k*(*ampa))*((3-pow(nam,2))/(4*pow(nam,3))))*cos( (2*arg) ))));
                                    arrayp[i]=g*(z+n)*rho;


                                  }
                                else
                                  {
                                    //if the wave field reach the bottom its velocity its calculated
                                    arg=tpi*( -(x/Lwa) +a2 );
                                    arrayx[i]=a1*((cos(arg)*kons1x));
                                    arrayy[i]=a1*((sin(arg)*kons1y));
                                    n=(*ampa)*( (cos(arg)+(((k*(*ampa))*((3-pow(nam,2))/(4*pow(nam,3))))*cos( (2*arg) ))));
                                    arrayp[i]=g*(z+n)*rho;

                                  }
                                  i++;
                                  x=x+(*dx);
                                }
                                while(i2<i){
                                  //We store the whole data from the arrays to a file
                                  fprintf(fp1, "%.3f,%.3f,",arrayx[i2],arrayy[i2]);
                                  fprintf(fp2, "%.3f,", arrayp[i2]);
                                  i2++;
                                }
                                fprintf(fp1, "\n");
                                fprintf(fp2, "\n");
                        i2=0;
                        x=0;
                        z=z+(*dz);
                    }
                    fprintf(fp1, "\n\n");
                    fprintf(fp2, "\n\n");
                    z=0;
                    t=t+(*dt);
            }
              fp1=NULL;
              fp2=NULL;
              in++;


        return 0;
}
