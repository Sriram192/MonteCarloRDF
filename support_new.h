#include<iostream>
#include<cstdlib>
#include<fstream>
#include<time.h>
using namespace std;

const int maxColloids=10000;
const int temp_loop = 100;
const int maxbin = 500;
//position coordinates of particles//
double colloid_x[maxColloids], colloid_y[maxColloids], colloid_z[maxColloids];
double colloidn_x[maxColloids], colloidn_y[maxColloids], colloidn_z[maxColloids];

double temp;
void initial_condition();
void writeconf(int ind);
void rdfsample(int t);
double potential(int i,double x,double y);
double randnum();
double periodic(double pos);
double pressurecalculation();
double totalenergy();
double ieq, ipr, mctry, mcacc;
double acceptance = mcacc/mctry;
double finalenergy;
void MC_NVT(),NVT_IN();
int ind;
int temp_count;
double energy;
double hist[maxbin];
int rdf_count;
ofstream of1;
