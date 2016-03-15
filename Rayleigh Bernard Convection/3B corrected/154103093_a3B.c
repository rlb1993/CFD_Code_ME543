#include<stdio.h>
#include<math.h>
#include<stdlib.h>
int main()
{
	FILE *f, *f1, *f2;
	f=fopen("input.txt","r");
	f1=fopen("CFDAssign3B.plt","w");
	f2=fopen("uvel.txt","w");

	int i,j,M, N, count=0;
	float L=5.0, H=1.0, dx, dy, b, Serror=10, Werror=10, Terror=10, sum=0, Re, x, y, C1, C2, C3, C4, C5, Coeff, ws, ww, wt, cc, Pr, Ra, Pe;
	
	//Reading the values from input file
	fscanf(f,"%d",&M); // M is the no of points along x axis
	fscanf(f,"%d",&N); // N is the no of points along y axis
	fscanf(f,"%f",&Re); // Re is the Reynold Number
	fscanf(f,"%f",&Pr); // Pr is the Prandlt Number
	fscanf(f,"%f",&Ra);// Ra is the Rayleigh Number
	fscanf(f,"%f%f%f",&ws,&ww,&wt); // w is the Relaxation factor
	fscanf(f,"%f",&cc); // cc is the convergence criteria

	Pe=Re*Pr;
	
	dx=L/(M-1);
	
	dy=H/(N-1);
	
	b=dx/dy;

	Coeff=(1.0/2.0)*(1.0/(1+(b*b)));


	float **Sold=(float **)malloc(N*sizeof(float));
	for(i=0;i<N;i++)
		Sold[i]=(float *)malloc(M*sizeof(float));

	float **Snew=(float **)malloc(N*sizeof(float));
	for(i=0;i<N;i++)
		Snew[i]=(float *)malloc(M*sizeof(float));

	float **Wold=(float **)malloc(N*sizeof(float));
	for(i=0;i<N;i++)
		Wold[i]=(float *)malloc(M*sizeof(float));

	float **Wnew=(float **)malloc(N*sizeof(float));
	for(i=0;i<N;i++)
		Wnew[i]=(float *)malloc(M*sizeof(float));

	float **U=(float **)malloc(N*sizeof(float));
	for(i=0;i<N;i++)
		U[i]=(float *)malloc(M*sizeof(float));

	float **V=(float **)malloc(N*sizeof(float));
	for(i=0;i<N;i++)
		V[i]=(float *)malloc(M*sizeof(float));

	float **Told=(float **)malloc(N*sizeof(float));
	for(i=0;i<N;i++)
		Told[i]=(float *)malloc(M*sizeof(float));

	float **Tnew=(float **)malloc(N*sizeof(float));
	for(i=0;i<N;i++)
		Tnew[i]=(float *)malloc(M*sizeof(float));

	//U-V BCs
	for(i=0;i<N;i++)
	{
		U[i][0]=0;
		U[i][M-1]=0;
	}

	for(j=0;j<M;j++)
	{
		U[0][j]=0;
		V[0][j]=0;
		U[N-1][j]=0;
		V[N-1][j]=0;
	}

	//Initialisation of Stream func and Vorticity
	for(j=0;j<M;j++)
	{
		Sold[0][j]=0.01;
		Snew[0][j]=0.01;
	}
	for(i=1;i<N-1;i++)
		for(j=1;j<M-1;j++)
			{
				Sold[i][j]=0;
				Snew[i][j]=0;
			}	
	for(i=1;i<N;i++)
	{
		Snew[i][0] = Snew[i-1][0];
		Sold[i][0] = Sold[i-1][0];
		Snew[i][M-1] = Snew[i-1][M-1];
		Sold[i][M-1] = Sold[i-1][M-1];
	}
	
	for(j=1;j<M-1;j++)
	{
		Snew[N-1][j] = Snew[N-1][j-1];
		Sold[N-1][j] = Sold[N-1][j-1];
	}
	
	for(i=0;i<N;i++)
		for(j=0;j<M;j++)
			{
				Wold[i][j]=0.5;
				Wnew[i][j]=0.5;
			}

	for(i=1;i<N-1;i++)
		for(j=1;j<M-1;j++)
			{
				Told[i][j]=0.5;
				Tnew[i][j]=0.5;
			}
	
	for(j=0;j<M;j++)
	{
		Tnew[0][j]=1;
		Told[0][j]=1;
		Told[N-1][j]=0;
		Tnew[N-1][j]=0;
	}
	
	for(i=0;i<N;i++)
	{
		Told[i][0]=Told[i][1];
		Tnew[i][0]=Tnew[i][1];
		Told[i][M-1]=Told[i][M-2];
		Tnew[i][M-1]=Tnew[i][M-2];
	}

	while(Serror > cc || Werror > cc || Terror > cc ) //Reducing the max error in phi and omega to Convergence Criteria
	{
		//Solving for Stream Function
		for(i=1;i<N-1;i++)
			for(j=1;j<M-1;j++)
				{
					Snew[i][j] = Coeff*((b*b)*Snew[i-1][j] + Snew[i][j-1] + (b*b)*Sold[i+1][j] + Sold[i][j+1] + dx*dx*Wold[i][j]);
					Snew[i][j] = (1-ws)*Sold[i][j] + ws*Snew[i][j];
				} 
		//Calculation of stream function error
		sum=0;
		for(i=1;i<N-1;i++)
		{
			for(j=1;j<M-1;j++)
			{
				sum+=pow((Snew[i][j]-Sold[i][j]),2);
			}
		}
		Serror=sqrt(sum/((M-2)*(N-2)));

		//Updating U-V using stream function
		for(i=1;i<N-1;i++)
			for(j=1;j<M-1;j++)
			{
				V[i][j]=-(Snew[i][j+1]-Snew[i][j-1])/(2.0*dx);
				U[i][j]=(Snew[i+1][j]-Snew[i-1][j])/(2.0*dy);
			}

		for(i=0;i<N;i++)
		{
			V[i][0]=V[i][1];
			V[i][M-1]=V[i][M-2];
		}

		for(i=1;i<N-1;i++)
			for(j=1;j<M-1;j++)
				{
						C1=(b*b)*(1-(V[i][j]*dy*Pe*0.5));
						C2=(1-(U[i][j]*dx*Pe*0.5));
						C3=(b*b)*(1+(V[i][j]*dy*Pe*0.5));
						C4=(1+(U[i][j]*dx*Pe*0.5));
						Tnew[i][j] = Coeff*(C1*Told[i+1][j] + C2*Told[i][j+1] + C3*Tnew[i-1][j] + C4*Tnew[i][j-1]);
						Tnew[i][j] = (1-wt)*Told[i][j] + wt*Tnew[i][j];
				}

		//Calculation of Temperature error
		sum=0;
		for(i=1;i<N-1;i++)
		{
			for(j=1;j<M-1;j++)
			{
				sum+=pow((Tnew[i][j]-Told[i][j]),2);
			}
		}
		Terror=sqrt(sum/((M-2)*(N-2)));

		for(i=0;i<N;i++)
		{
			Tnew[i][0]=Tnew[i][1];
			Tnew[i][M-1]=Tnew[i][M-2];
		}

		//updating vorticity BCs using Stream function
		for(j=0;j<M;j++)
		{
			Wnew[0][j]=-(2.0/(dy*dy))*(Snew[1][j]-Snew[0][j]);
			Wnew[N-1][j]=-(2.0/(dy*dy))*(Snew[N-2][j]-Snew[N-1][j]);
		}

		for(i=0;i<N;i++)
		{
			Wnew[i][0]=Wnew[i][1];
			Wnew[i][M-1]=Wnew[i][M-2];
		}

		//Solving for vorticity
		for(i=1;i<N-1;i++)
			for(j=1;j<M-1;j++)
				{
						C1=(b*b)*(1-(V[i][j]*dy*Re*0.5));
						C2=(1-(U[i][j]*dx*Re*0.5));
						C3=(b*b)*(1+(V[i][j]*dy*Re*0.5));
						C4=(1+(U[i][j]*dx*Re*0.5));
						C5=(Ra/Pe)*dx*0.5;
						Wnew[i][j] = Coeff*(C1*Wold[i+1][j] + C2*Wold[i][j+1] + C3*Wnew[i-1][j] + C4*Wnew[i][j-1] + C5*(Tnew[i][j+1]-Tnew[i][j-1]));
						Wnew[i][j] = (1-ww)*Wold[i][j] + ww*Wnew[i][j];
				}
		//Calculation of vorticity error
		sum=0;
		for(i=1;i<N-1;i++)
		{
			for(j=1;j<M-1;j++)
			{
				sum+=pow((Wnew[i][j]-Wold[i][j]),2);
			}
		}
		Werror=sqrt(sum/((M-2)*(N-2)));
		//updating the value of phi and omega
		for(i=0;i<N;i++)
			for(j=0;j<M;j++)
			{
				Sold[i][j]=Snew[i][j];
				Wold[i][j]=Wnew[i][j];
				Told[i][j]=Tnew[i][j];
			}

		count++;
		printf("\nIteration No: %d Serror: %f Werror: %f Terror: %f",count,Serror,Werror,Terror);
		fprintf(f2,"\n%d\t%f\t%f\t%f",count,Serror,Werror,Terror);
	}
	fprintf(f1,"\nZONE I=%d, J=%d",M,N);
	for(y=0,i=0;i<N;i++,y=y+dy)
	{
		for(x=0,j=0;j<M;j++,x=x+dx)
		{
			fprintf(f1,"\n%f\t%f\t%f\t%f\t%f\t%f\t%f",x,y,U[i][j],V[i][j],Sold[i][j],Wold[i][j],Told[i][j]);
		}
	}
	for(y=0,i=0;i<N;i++,y=y+dy)
	{
		fprintf(f2,"\n%f\t%f\t%f",y,(U[i][40]+U[i][41])*0.5,(U[i][50]+U[i][51])*0.5);
	}		
}