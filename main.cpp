#include <iostream>
#include <cstdlib>
#include <math.h>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>

#define EPS 0.00001

using namespace std;

// struktura
struct taskas
{
    double x; // pirmas taskas
    double y; // antras taskas
    double f; // funkcija nuo x ir y
};


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


int main()
{
    srand(time(0));
    cout << "Monte Carlo realizacijos pradzia" << endl;
    cout << "(RANDOM SEARCH METHOD)" << endl;
    double mas[3]; // laikinoms reiksmes laikyti
    double fspr; // funkcijos reiksme taskuose x ir y
    double max, min, suma = 0, x_min, y_min, x_max, y_max;
    taskas rez; //strukturos kintamasis

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
    cout << "F-jos minimumas yra "<< min <<endl;
    cout << "taske ("<< x_min <<", "<<y_min<< ")"<<endl;
    cout<< "-------------------------------------"<<endl;
    cout << "F-jos maximumas yra "<< max <<endl;
    cout << "taske ("<< x_max <<", "<<y_max<< ")"<<endl;
    cout<< "-------------------------------------"<<endl;
    cout << "F-jos vidurkis yra "<< (suma/(j-1)) <<endl;
    cout<< "-------------------------------------"<<endl;
    cout << "Ciklas prasuktas "<< (j-1)<<" kartu" <<endl;
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

    cout << "Trys maziausios reiksmes (su vector sort): " << endl;

    for (int i=0; i<3; i++) // iki reiksmes.size()
        {
            cout << "=== " << i+1 << " ===" << endl;
            //cout << reiksmes[i].x << endl;
            //cout << reiksmes[i].y << endl;
            cout << reiksmes[i].f << endl;
        }

    //sort(reiksmes.begin(), reiksmes.end(), mazejanciai); // mazejanciai

    // rusiuoja su heapsort
     heapsort(&reiksmes[0],(j-1));

     cout << "Trys maziausios reiksmes (su heapsort): " << endl;

     for (int i=0; i<3; i++)
     {
         cout << "=== " << i+1 << " ===" << endl;
         cout << "x = " << reiksmes[i].x << " y = " << reiksmes[i].y << " f = " << reiksmes[i].f << endl;
     }



    return 0;
}
