#include <iostream>
#include <cstdlib>
#include <math.h>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>

#define EPS 0.1
#define N 2

using namespace std;

// struktura
struct taskas
{
    double x; // pirmas taskas
    double y; // antras taskas
    double f; // funkcija nuo x ir y
};


// Vektoriaus begalines (max) normos funkcijos deklaracija
double Vector_Max_Norm(double v[], int n);

// Greiciausio nusileidimo (angl. Steepest Descent) metodo deklaracija
int  Steepest_Descent(double (*f)(double *), void (*df)(double *, double *),
     int (*stopping_rule)(double*, double, double*, double, double*, int, int),
                          double a[], double *fa, double *dfa, double cutoff,
						double cutoff_scale_factor, double tolerance, int n);

// Generuoja atsitiktini realu skaiciu tarp dLow and dHigh
double GetRandomNumber(double dLow, double dHigh){
    return static_cast<double>(rand())/RAND_MAX*(dHigh-dLow) + dLow;
}

// Apskaiciuoja Six-hump Camel Back funkcijos reiksme taske x
double SixHumpCamelBack(double *x){
    return (4-2.1*x[0]*x[0]+x[0]*x[0]*x[0]*x[0]/3)*x[0]*x[0] + x[0]*x[1] +
    (-4+4*x[1]*x[1])*x[1]*x[1];
}
// Apskaiciuoja Six-hump Camel Back gradiento reiksme taske x
void SixHumpCamelBackGradient(double *x, double *fGrad){
    fGrad[0] = 8*x[0]-8.4*x[0]*x[0]*x[0]+2*x[0]*x[0]*x[0]*x[0]*x[0]+x[1];
    fGrad[1] = x[0]-8*x[1]+16*x[1]*x[1]*x[1];
}

// Algoritmo sustojimo salyga kontroliuojanti funkcija
int StoppingRule(double* a, double fa, double* x, double fx, double* dfa, int
iteration, int n){
	double fEps = abs(fx - fa); // Funkcijos reiksmiu skirtumas
	double xa[n];
	for(int i = 0; i < n; ++i) xa[i] = x[i]-a[i];
	double xEps = Vector_Max_Norm(xa, 2); // Argumento skirtumo norma
	double dfaEps = Vector_Max_Norm(dfa, 2); // Gradiento norma
	if(iteration > 3)
		return -6;
	else
		return 0;
}

// heapsort funkcijos
void swap(taskas *x, taskas *y)
{
    taskas rez;
    rez = *x;
    *x = *y;
    *y = rez;
}

void heapsort(taskas *a,int n)
{
    int i,s,f;
    for(i=1;i< n;++i)
    {
        s=i;
        f=(s-1)/2;
        while( a[f].f< a[s].f)
        {
            swap(&a[f],&a[s]);
            s=f;
            if(s==0)
            break;
            f=(s-1)/2;
        }
    }
    for(i=n-1;i>=1;--i)
    {
        swap(&a[0],&a[i]);
        f=0;
        s=1;
        if(i==1)
        break;
        if(i>2)if(a[2].f>a[1].f)s=2;
        while( a[f].f< a[s].f )
        {
            swap(&a[f],&a[s]);
            f=s;
            s=2*f+1;
            if(i>s+1 )if(a[s+1].f>a[s].f)s=s+1;
            if (s>=i)
            break;
        }
    }
}


// funkcija six hump
double sixhump(double *x)
{
    double f = (4- 2.1 * x[0] * x[0] + (pow(x[0],4))/3) * x[0] * x[0] + x[0] * x[1] + ((-4+4*x[1]*x[1])*(pow(x[1],2)));
    return f;
}

// funkcija rand
double generuoti(double a, double b)
{
return a + (b - a) * rand() / ((double) RAND_MAX);
}

bool didejanciai(const taskas &a, const taskas &b) // funkcijas vektoriaus rusevimui didejanciai
{
    return a.f < b.f;
}

bool mazejanciai(const taskas &a, const taskas &b) // funkcijas vektoriaus rusevimui mazejanciai
{
    return a.f > b.f;
}

// optimizavimas
void optimizavimas(double x, double y, double f)
{
    double x_old = x, y_old=y, x_new=0, y_new=0;
    double eps = 0.01;
    double fr;
    double precision = 0.00001;
    double mas[3];

    while (sqrt((x_new-x_old)*(x_new-x_old)+(y_new-y_old)*(y_new-y_old)) > precision)
    {

        mas[0] = x_old;
        mas[1] = y_old;
        fr = sixhump(&mas[0]);
        x_new = x_old - eps * fr;
        y_new = x_old - eps * fr;

        x_old = x_new;
        y_old = y_new;
    }
    cout << "\nx = " << x_new << " y = " << y_new;
}

int main()
{
    srand(time(0));
    cout << "Monte Carlo realizacijos pradzia" << endl;
    cout << "(RANDOM SEARCH METHOD)" << endl;
    double mas[3]; // laikinoms reiksmes laikyti
    double fspr; // funkcijos reiksme taskuose x ir y
    double max, min, suma = 0, x_min, y_min, x_max, y_max;
    taskas rez; //strukturos kintamasis

    // masyvas su trim maziausiom reiksmem
    double a[6];


    vector <taskas> reiksmes(0);

    mas[0] = generuoti(-1.9, 1.9); // sugeneruoja reiksme
    mas[1] = generuoti(-1.1, 1.1); // sugeneruoja reiksme
    fspr = sixhump(&mas[0]); // apskaiciuoja funkcija taskuose x ir y

    max = fspr;
    min = fspr;

    rez.x = mas[0];
    rez.y = mas[1];
    rez.f = fspr;

    reiksmes.push_back(rez);

    int j=1;
    while (abs(-1.031628453 - fspr) > EPS)
    {
        mas[0] = generuoti(-1.9, 1.9); // sugeneruoja reiksme
        mas[1] = generuoti(-1.1, 1.1); // sugeneruoja reiksme
        fspr = sixhump(&mas[0]); // apskaiciuoja funkcija taskuose x ir y

        max = fspr;
        min = fspr;

        rez.x = mas[0];
        rez.y = mas[1];
        rez.f = fspr;

        reiksmes.push_back(rez);

        suma = suma + fspr;

        if(min > fspr)
        {
            min = fspr;
            x_min = mas[0];
            y_min = mas[1];
        }
        if(max < fspr)
        {
            max = fspr;
            x_max = mas[0];
            y_max = mas[1];
        }
        j++;

    }

    cout<< "====================================="<<endl;
    cout << "\nF-jos minimumas yra "<< min <<endl;
    cout << "taske ("<< x_min <<", "<<y_min<< ")"<<endl;
    cout<< "-------------------------------------"<<endl;
    cout << "\nF-jos maximumas yra "<< max <<endl;
    cout << "taske ("<< x_max <<", "<<y_max<< ")"<<endl;
    cout<< "-------------------------------------"<<endl;
    cout << "\nF-jos vidurkis yra "<< (suma/(j-1)) <<endl;
    cout<< "-------------------------------------"<<endl;
    cout << "\nCiklas prasuktas "<< (j-1)<<" kartu" <<endl;
    cout<< "====================================="<<endl;

    /*for (int i=0; i<3; i++) // iki reiksmes.size()
        {
            cout << "=== " << i+1 << " ===" << endl;
            //cout << reiksmes[i].x << endl;
            //cout << reiksmes[i].y << endl;
            cout << reiksmes[i].f << endl;
            cout << "=============" << endl;
        }
    */

    // rusiuoja per vector su sort

    sort(reiksmes.begin(), reiksmes.end(), didejanciai); // didejanciai

    cout << "\nTrys maziausios reiksmes (su vector sort): " << endl;

    for (int i=0; i<3; i++) // iki reiksmes.size()
        {
            cout << "\n=== " << i+1 << " ===" << endl;
            //cout << reiksmes[i].x << endl;
            //cout << reiksmes[i].y << endl;
            cout << reiksmes[i].f << endl;
        }

    //sort(reiksmes.begin(), reiksmes.end(), mazejanciai); // mazejanciai

    // rusiuoja su heapsort
     heapsort(&reiksmes[0],(j-1));

     cout << "\nTrys maziausios reiksmes (su heapsort): " << endl;

     for (int i=0; i<3; i++)
     {
         cout << "\n=== " << i+1 << " ===" << endl;
         cout << "x = " << reiksmes[i].x << " y = " << reiksmes[i].y << " f = " << reiksmes[i].f << endl;
     }

     // reiksmiu priskyrimas
     a[0] = reiksmes[0].x;
     a[1] = reiksmes[0].y;
     a[2] = reiksmes[1].x;
     a[3] = reiksmes[1].y;
     a[4] = reiksmes[2].x;
     a[5] = reiksmes[2].y;


     //optimizavimas(reiksmes[0].x, reiksmes[0].y, reiksmes[0].f);

     // Steepest Descent prijungiu cia
    cout << "\n===========================================================\n";
    // atspausdiname maziausias reiksmes:
     for(int i=0; i<6; i++)
        cout << a[i] << " -- ";

    double region[] = {-1.9, 1.9, -1.1, 1.1};

    for (int i=0; i<3; i=i+2)
    {

    cout << "\n\n===PRIES===" << endl;
    cout << "x = " << a[i] << " y = " << a[i+1] << " F = " << reiksmes[i].f << endl;




    double fa = SixHumpCamelBack(&a[i]); // Funkcijos reiksme pradiniame taske a
    double dfa[N];
    SixHumpCamelBackGradient(a, dfa); // Funkcijos gradiento reiksme taske a
    double cutoff = 1.0, cutoff_scale_factor = 1.0; // Pap. parametrai
    double tolerance = 0.01;
    int err = Steepest_Descent( SixHumpCamelBack, SixHumpCamelBackGradient, StoppingRule,
    a, &fa, dfa, cutoff, cutoff_scale_factor, tolerance, N);

    switch (err)
    {
        case 0:
        cout << "Success" << endl;
        break;
        case -1:
        cout << "In the line search three points are collinear." << endl;
        break;
        case -2:
        cout << "In the line search the extremum of the parabola through the three points is a maximum." << endl;
        break;
        case -3:
        cout << "Int the line search the initial points failed to satisfy the condition that x1 < x2 < x3 and fx1 > fx2 < fx3." << endl;
        break;
        case -4:
        cout << "Not enough memory." << endl;
        break;
        case -5:
        cout << "The gradient evaluated at the initial point vanishes." << endl;
        case -6:
        cout << "\nExceed maximal number of iterations." << endl;
        break;
    }
    cout << "\nGreiciausio nusileidimo (angl. Steepest Descent) metodu" << endl;
    cout << "surastas sprendinys yra:" << endl;
    cout << "\nxMin = (" << a[0] << ", " << a[1] << ")" << endl;
    cout << "\nf(xMin) = " << fa << endl;

    }


    return 0;
}
