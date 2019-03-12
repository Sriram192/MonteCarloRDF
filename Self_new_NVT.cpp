//Written by Sriram
//28/09/2018
//Self assembly of honeycomb potential taken from Torquato et. al. 2006
//Also prints rdf whose algorithm has to be verified for the 2D case.
//The temperature schedule for cooling in simulated annealing is linear in this particular code.

#include<iostream>
#include<math.h>
#include "support_new.h"
#include<string>
#include<fstream>
//declared as global variables
int number_of_colloids = 500;
double alpha = 1.45;//1.45 for honeycomb
double box_length = sqrt(number_of_colloids * alpha);
double cutoff = 2.5;
double intemp = 3.0;
double dx=0.1, dy=0.1;// dz=0.1;
double delr = box_length/double(2*maxbin);
int main(int argc, char *argv[])
{

    srand(time(0));
    if(argc > 1) {
        ind = atoi(argv[1]);
    }
    char outfile[100];
    ofstream of;
    sprintf( outfile, "out.%d.txt", ind);//creating a new file
    of.open (outfile, ofstream::out | ofstream::trunc);


cout<<"box length is "<<box_length<<endl;
initial_condition();
//as mentioned in the paper.


//for(int i=0;i<30;i++)
//{
//if(i==0) temp[i] = 0.22;
//temp[i+1] = double(temp[i])*pow(0.95,i);
//temp_count += 1;
//cout<<temp[i]<<"\n";
//}
//cout<<temp_count<<endl;

double finalenergy=0;
for (int i=0; i< number_of_colloids; i++)
{
finalenergy += potential(i,colloid_x[i],colloid_y[i]);

}
cout<<"the total initial energy is-------> "<<finalenergy<<endl;

cout<<"Enter the no of equilibrium cycles"<<endl;
cin>>ieq;
cout<<endl;
cout<<"Enter the no of production cycles"<<endl;
cin>>ipr;
cout<<endl;
mctry=0;mcacc=0;acceptance=0;

//----------------------------Warming----------------------------------//
int warm_cycle = 50000;
for(int i=0;i<warm_cycle;i++) // note the number of iterations for the warming loop.
{
	NVT_IN();
		if(i%100==0)
		{
//		cout<<"progress (warming) : "<<(double(i)/warm_cycle)*100<<"%"<<endl;
		acceptance = mcacc/mctry;
        		if(acceptance>0.5)
        		{
        		dx = dx*1.05;
        		dy = dy*1.05;
//      		dz = dz*1.05;
        		}
        		if(acceptance<0.4)
        		{
        		dx = dx*0.95;
        		dy = dy*0.95;
//      		dz = dz*0.95;
        		}
        		if (dx > box_length/2.) dx = box_length/10.;
        		if (dy > box_length/2.) dy = box_length/10.;
        	cout<<"progress (warming) : "<<(double(i)/warm_cycle)*100<<"%"<<"	"<<" acceptance  "<<acceptance<<endl;
		mctry=0;
		mcacc=0;
		acceptance = 0;

		}
}

char warm[100];

ofstream ofs;
    sprintf( warm, "warm_config.xyz");
    ofs.open (warm, ofstream::out | ofstream::trunc);
    if (ofs.is_open())
    {
        ofs << number_of_colloids<< "\n";
        ofs << "\n";

        for (int i=0; i<number_of_colloids; i++)
        {
            ofs << 'h'<<"  "<<colloid_x[i] << "  " << colloid_y[i]<<"  "<<"0.0"<< " \n ";// << "  "<< colloid_z[i] << "\n";
        }
    }
    else
    {
    cerr << "unable to open file for config output \n";
    }
    ofs.close();
//-------------------------------------MC main---------------------------//

temp = 0.22;
double tol = 1.0e-3;
char logfile[100];
	ofstream of1;
	sprintf(logfile,"log.txt");
	of1.open(logfile, ofstream::out | ofstream::app);
		if(of1.is_open()){
			of1<<"Temperature run"<<"     "<<"Potential energy/particle"<<"     "<<"Pressure"<<"     "<<"area/particle"<<"\n";
			}
int temp_count = 0;			
while (temp > tol)
{


energy = totalenergy();
	mctry=0;mcacc=0;acceptance=0;
	for(int i =0; i<ieq;i++)	
	{
	MC_NVT();

		if(i%100==0)
		{
		acceptance = mcacc/mctry;
			if(acceptance>0.5)
			{
			dx = dx*1.05;
			dy = dy*1.05;
//			dz = dz*1.05;
			}
			if(acceptance<0.4)
			{
			dx = dx*0.95;
			dy = dy*0.95;
//			dz = dz*0.95;
			}

		cout<<"progress(eq)..."<<(i/ieq)*100.0<<"%"<<"       "<<"acceptance ratio "<<acceptance<<"       temp	"<<temp<<endl;
			if (of.is_open())
        	    	{
        	       	of << i << "   "<< acceptance << "   "<<totalenergy() << "    "<< pressurecalculation()<<"	"<<temp<<"   \n";
		    	}
			else{cerr <<"unable to open output file";}
		mctry=0;
		mcacc=0;
		acceptance = 0;
		}	
	}

	int count = 0;
	double inspot=0.0, inspress=0.0, pot2 = 0.0, press2 = 0.0, insalpha = 0.0;
	double avpress = 0.0, avpot = 0.0, avalpha = 0.0, alpha2 = 0.0;
	for(int i=0; i<ipr; i++)
	{
		MC_NVT();

		if(i%100==0)
		{
		acceptance = mcacc/mctry;
		cout<<"progress(pro)..."<<(i/ipr)*100.0<<"%"<<"       "<<"acceptance ratio "<<acceptance<<"	   temp  "<<temp<<endl;
		           
	               of << i+ieq << "   "<< acceptance << "   "<<totalenergy() << "    "<< pressurecalculation()<<"	"<<temp<<"   \n";
	            
	       
		mctry=0;mcacc=0; acceptance=0;
	
		writeconf(1);
		inspot = energy;
		avpot += inspot;
		pot2 += inspot*inspot;
		inspress = pressurecalculation();
		avpress += inspress;
		press2 += inspress*inspress;
		insalpha = number_of_colloids/(box_length*box_length);//*box_length);
		avalpha += insalpha;
		alpha2 += insalpha*insalpha;
		count +=1;
		}
		
		
		
//			
//			rdf_count +=1;
		

	}
	avpot /= double(number_of_colloids *count);
	pot2 /= double(number_of_colloids*number_of_colloids*count);
	avpress /= double(count);
	press2 /= double(count);
	avalpha /= double(count);
	alpha2 /= double(count);
	
	if(of1.is_open()){
			of1<<temp_count<<"       "<< avpot<<"+-"<<sqrt(pot2 - avpot*avpot)<<"      "<<avpress << " +- " << sqrt(press2 - avpress*avpress)<<"      "<< 1.0/avalpha<<"\n";
			}else{cerr<< "unable to open log file";}	

	rdfsample();
//	temp = temp * exp(-0.2);
	temp = temp*0.95;//New temperature schedule tried on 10/11/2018
	temp_count += 1;
	
	

	
//write condition for mcacc and pos updation with delta r.
}
of1.close();
of.close();
return 0;
}

//main ends here;

//-------------------------------MC_NVT------------------------------------//

void MC_NVT()
{
//select a random index for colloid:
mctry+=1;
int index = randnum() * number_of_colloids;
colloidn_x[index] = colloid_x[index] + dx*(randnum() - 0.5);
colloidn_y[index] = colloid_y[index] + dy*(randnum() - 0.5);
//colloidn_z[index] = colloid_z[index] + dz*(randnum() - 0.5);

colloidn_x[index] = periodic(colloidn_x[index]);
colloidn_y[index] = periodic(colloidn_y[index]);
//colloidn_z[index] = periodic(colloidn_z[index]);


double Eold, Enew, dE_T;

Eold = potential(index,colloid_x[index],colloid_y[index]);
Enew = potential(index,colloidn_x[index],colloidn_y[index]);

dE_T= double(-(Enew-Eold)/temp);
	if(randnum() < exp(dE_T))
	{

	colloid_x[index] = colloidn_x[index];
	colloid_y[index] = colloidn_y[index];
//	colloid_z[index] = colloidn_z[index];
	energy += Enew - Eold;
	mcacc += 1;
	}
return;
}

//----------------------------------randum number-----------------------//

double randnum()
{
return double(rand())/(double(RAND_MAX));
}

//--------------------function for initial condition---------------------//

void initial_condition()
{
//dimension of a unit cell for fcc lattice:
/*In an a fcc lattice, each unit cell has 4 particles. Therefore, we first estimate the number of cells for a given number of particles
  which is given by number_cells = (number_of _colloids)/4. Dimension of each cell can be calculated as follows:
Volume of 1 cell * number_cell = volume of box */

cout<<"Enter 1 for FCC or \nEnter 2 for BCC\nor Enter 3 for random initial configuration"<<endl;
int init;
cin>>init;
cout<<"\n";

if(init == 1)//condition for fcc
{
	int index = 0; 
	int baseunit = ceil(pow(number_of_colloids/4.,(1./3.)));
	double cell_length = double(box_length)/baseunit;

//baseunit is the number of cells that can be accomodated in a single axis. Therfore, the length of each cell will be calculated accordingly from the box length. 
	for(int ix=0; ix<baseunit; ix++)//index starts from zero because we consider the coordinates from origin at (0,0,0)
	{
        	for(int iy=0; iy<baseunit; iy++)
        	{
                	for(int iz=0; iz<baseunit; iz++)
                	{
                	colloid_x[index] = ix * cell_length;
                	colloid_y[index] = iy * cell_length;
                	colloid_z[index] = iz * cell_length;
                 
                	colloid_x[index+1] = colloid_x[index]+cell_length/2.0;
                	colloid_y[index+1] = colloid_y[index]+cell_length/2.0;
                	colloid_z[index+1] = colloid_z[index];
                
			colloid_x[index+2] = colloid_x[index]+cell_length/2.0;
                	colloid_y[index+2] = colloid_y[index];
                	colloid_z[index+2] = colloid_z[index]+cell_length/2.0;

                	colloid_x[index+3] = colloid_x[index];
                	colloid_y[index+3] = colloid_y[index]+cell_length/2.0;
                	colloid_z[index+3] = colloid_z[index]+cell_length/2.0;

                	index+=4;
               		}	
        	}
	}
}	

if(init ==2)
{
	int index = 0;
	int baseunit = ceil(pow(number_of_colloids,(1./2.)));
	double cell_length = double(box_length)/baseunit;
	for(int ix=0; ix<baseunit; ix++)
	{
        	for(int iy=0; iy<baseunit; iy++)
        	{
//                	for(int iz=0; iz<baseunit; iz++)
//                	{
                	colloid_x[index] = ix * double(cell_length);
                	colloid_y[index] = iy * double(cell_length);
//               	colloid_z[index] = iz * double(cell_length);
                	index+=1;
//                	}
        	}
	}
}

if(init == 3)
{
double tx, ty, tz, xij, yij, zij, r2, r1;
//random configuration
colloid_x[0] = (randnum() - 0.5)*box_length;
colloid_y[0] = (randnum() - 0.5)*box_length;
//colloid_z[0] = (randnum() - 0.5)*box_length;
for (int i=1;i<number_of_colloids; i++) {
        do {
            tx = (randnum() - 0.5)*box_length;
            ty = (randnum() - 0.5)*box_length;
//          tz = (randnum() - 0.5)*box_length;
            for (int j=0; j<i; j++) {
                xij = colloid_x[j] - tx;
		yij = colloid_y[j] - ty;
//		zij = colloid_z[j] - tz;
		xij = periodic(xij);
		yij = periodic(yij);
//		zij = periodic(zij);

                r2 = (xij*xij + yij*yij);// + zij*zij;
                r1 = sqrt(r2);
cout<<r1<<endl;
                if(r1<1.0) break;
		}
        }while (r1<1.0);

        colloid_x[i] = tx;
        colloid_y[i] = ty;
//      colloid_z[i] = tz;
    }
}
writeconf(0); 

return;
}//initial_condition ends here

//-------------------------Writing configuration into a file------------------------------//

void writeconf(int ind)
{
if(ind==0)
{
    char IntStr[80];
    ofstream of;
    sprintf( IntStr, "initial_config.xyz");
    of.open (IntStr, ofstream::out | ofstream::app);
    if (of.is_open())
    {
        of << number_of_colloids<< "\n";
        of << "\n";

        for (int i=0; i<number_of_colloids; i++)
        {
            of << 'h'<<"  "<<colloid_x[i] << "  " << colloid_y[i]<<"  "<<"0.0"<<" \n";// << "  "<< colloid_z[i] << "\n";
        }
    }

    else
    {
    cerr << "unable to open file for config output \n";
    }
    of.close();
}
if(ind==1)
{
char IntStr[80];
    ofstream of;
    sprintf( IntStr, "running_config.xyz");
    of.open (IntStr, ofstream::out | ofstream::app);
    if (of.is_open())
    {
        of << number_of_colloids<< "\n";
        of << "\n";

        for (int i=0; i<number_of_colloids; i++)
        {
            of << 'h'<<"  "<<colloid_x[i] << "  " << colloid_y[i]<<"  "<<"0.0"<< " \n ";// << "  "<< colloid_z[i] << "\n";
        }
    }
    else
    {
    cerr << "unable to open file for config output \n";
    }
    of.close();
}
    return;
}

//-----------------------------potential energy----------------------------------//

double potential(int i,double x,double y)
{
double xij, yij;
//In actual simulation, get the cutoff radius as an input.
double cutoff = 2.5, cutoff2;
double r2,r,U=0.0,r12,r10;

	for(int j=0; j<number_of_colloids; j++)
	{
		if(j != i){
		xij = double(x-colloid_x[j]);
		yij = double(y-colloid_y[j]);
//		zij = double(z-colloid_z[j]);
//minimum image convention	
		xij = periodic(xij);
		yij = periodic(yij);
//		zij = periodic(zij);

//in actual simulation, check for periodic BC
		r2 = double((xij*xij) + (yij*yij));
		
		cutoff2 = double(cutoff*cutoff);
	
		if(r2 < cutoff2)
		{
		r = sqrt(r2);//defining radial distance from the particle
		r12 =1.0/pow(r2,6);
		r10 = 1.0/pow(r2,5);
		U += (5.0*r12) - (5.89*r10) + 17.9*exp(-2.49*r) - 0.4*exp(-40.0*(pow((r-1.823),2.0)));
		}
		}
	}

return U;
}

//---------------------------------------Periodic boundary conditions------------------------//

double periodic(double pos)
{
double hbox = double(box_length)/2.0;
	
	if(pos<hbox) pos += box_length;
	if(pos>hbox) pos -= box_length;
	
	
return pos;
}

//------------------------------total energy calculation--------------------------------------//

double totalenergy()
{
double pot =0;
	for(int i=0; i<number_of_colloids;i++)
	{
	pot += potential(i,colloid_x[i],colloid_y[i]);
	}
//double r3 = 1.0/(cutoff*cutoff*cutoff);
//double Utail = 8*3.14*density*r3*number_of_colloids*(((1.0/9.0)*(pow(r3,2))) - (1.0/3.0));
return (double(pot/2));// + Utail);
}

//------------------------------pressure calculation-----------------------------------------//

double pressurecalculation()
{
double xij,yij,zij,rij,r2,r12,r10,r;
double vir =0;
double Ptail, Pid;
	for(int i=0;i<number_of_colloids -1;i++)
	{
		for(int j=i+1;j<number_of_colloids;j++)
		{
		xij = colloid_x[i] - colloid_x[j];
		yij = colloid_y[i] - colloid_y[j];
//		zij = colloid_z[i] - colloid_z[j];
		xij = periodic(xij);//minimum image convention
		yij = periodic(yij);
//		zij = periodic(zij);
		r2 = xij*xij + yij*yij;

		if(r2<(cutoff*cutoff))
			{
			r = sqrt(r2);
			r10 = 1.0/(pow(r2,5.0));
			r12 = 1.0/(pow(r2,6.0));
			vir += (60*r12) - (58.9*r10) + (44.571 * r * exp(-2.49*r)) - (32 *(r-1.823) *r*exp(-40*pow((r-1.823),2)));
			}
		}
	}

//double r3 = 1.0/(cutoff*cutoff*cutoff);
//Ptail = (16.0/3.0)*density*density*r3*(((2.0/3.0)*pow(r3,2)) - 1.0);
Pid = alpha*temp;
return (Pid) + (vir/(3*(box_length*box_length)));
}

//----------------------------MC function for warm congfiguration-------------------------------//

void NVT_IN()
{
//select a random index for colloid:
mctry+=1;
int index = randnum() * number_of_colloids;
colloidn_x[index] = colloid_x[index] + dx*(randnum() - 0.5);
colloidn_y[index] = colloid_y[index] + dy*(randnum() - 0.5);
//colloidn_z[index] = colloid_z[index] + dz*(randnum() - 0.5);

colloidn_x[index] = periodic(colloidn_x[index]);
colloidn_y[index] = periodic(colloidn_y[index]);
//colloidn_z[index] = periodic(colloidn_z[index]);


double Eold, Enew, dE_T;

Eold = potential(index,colloid_x[index],colloid_y[index]);
Enew = potential(index,colloidn_x[index],colloidn_y[index]);

dE_T= double(-(Enew-Eold)/intemp);
        if(randnum() < exp(dE_T))
        {

        colloid_x[index] = colloidn_x[index];
        colloid_y[index] = colloidn_y[index];
//      colloid_z[index] = colloidn_z[index];
        mcacc += 1;
        }
return;
}
void rdfsample()
{
char RDFfile[80];//what is this?

ofstream ofrdf;
    sprintf( RDFfile, "rdf%d.dat",temp_count);
    ofrdf.open (RDFfile, ofstream::out | ofstream::app);

double pi = 3.1415926;
double boxinv = 1.0/box_length;//why are we calaculating this value?
double hbox = 0.5*box_length;//half box length
double rij2, r1;//rij2 is rij^square and r1 is the magnitude of the distance between the particles
int nsample = 0, bin;
double rijx,rijy,rijz;

int rdf_count = 0;
	if(rdf_count==0)
	{
		for(int i=0; i<maxbin; i++)
		{
		hist[i]=0;//what is hist?
		}
	}
//this condition just gives "true" for all values of t so why is this expression used here?
nsample += 1;
for(int k=0;k<100;k++)
{	
cout<<"RDF Step "<<k<<"\n";
	for(int m=0;m<10000;m++)
	{
	MC_NVT();
	}
	for (int i = 0; i<number_of_colloids; i++) 
	{
    		for(int j=i+1; j<number_of_colloids; j++) 
    		{
    		rijx = colloid_x[i] - colloid_y[j] ;
       		rijy = colloid_y[i] - colloid_y[j] ;
//       rij.z = colloid_z[i] - colloid_z[j] ;
       
       		rijx = periodic(rijx);
       		rijy = periodic(rijy);
//     rijz = periodic(rijz);
	  
       		rij2 = rijx*rijx + rijy*rijy;// + rij.z*rij.z;
       		r1 = sqrt(rij2);
       /*Instead of calculating square root of "r" outside the if condition
        *we can check whether r1^2 < hbox^2 and then if the condition satisfies
        *we calcuate the square root of "r". this will be computationally 
		*less intensive*/
       		if(r1 <= hbox)
	       		{
        		bin = int(r1/delr);//Address of the particle i.e the bin number to which the particle j belongs wrt the ith partice. 
			hist[bin] += 2.0;//this is because there are pairs. jthe particle belongs to bin and vice versa fo jthe particle
			} 
		}
	}
rdf_count++;
}
	
	double idcon = 3.1415*(1/alpha);
	double rlow, rup, ideal, rbin;
	if(ofrdf.is_open()) {
	for(int i=0; i<maxbin; i++)
	{
	rlow = double(i)*delr;
	rup = rlow + delr;
	ideal = idcon*(rup*rup - rlow*rlow);
	hist[i] = hist[i]/double(number_of_colloids*rdf_count)/ideal;
	rbin = rlow + delr/2.0;
           	 ofrdf << rbin <<"  "<<hist[i]<< "\n";
	}	
	}
	else {cerr << "unable to open file for rdf output \n";}
    	ofrdf.close();

return;
}
//---------------------------------THE END----------------------------------------------------------//





















