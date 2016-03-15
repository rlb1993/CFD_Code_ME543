/*This program solve the algebraic syatem of equation obtained from 
finite difference formulation of 2D steady heat conduction
by ADI Method.......................................
--->Input taken are: No.of Grid points along X and Y direction
--->Output generated is a file which stores the values of (x,y,T)*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
int main()
{
	clock_t begin,end;
	FILE *f;
	f=fopen("Temp_data_ADI.txt","w");
	int count=0,M, N ;
	float L=1, H=1, delta_x, delta_y, b,e=1,sum=0,x,y;
	printf("\nEnter the number of grid points along x axis:");
	scanf("%d",&M);
	printf("\nEnter the number of grid points along y axis:");
	scanf("%d",&N);
	delta_x=L/(M-1);
	printf("\n%f",delta_x);
	delta_y=H/(N-1);
	printf("\n%f",delta_y);
	b=delta_x/delta_y;
	printf("\n%f",b);
	int i,j,k;
	

	//Temperature is stored as a three dimensional array 
        //where i,j are the coordinate where T is calculated
	//and k keep track of the iteration count. 	
	float ***T=(float ***)malloc(N*sizeof(float));
	for(i=0;i<N;i++)
	T[i]=(float **)malloc(M*sizeof(float));
	for(i=0;i<N;i++)
	for(j=0;j<M;j++)
	T[i][j]=(float *)malloc(50000*sizeof(float));

	float *A=(float *)malloc(M*sizeof(float));
	float *B=(float *)malloc(M*sizeof(float));
	float *C=(float *)malloc(M*sizeof(float));
	float *D=(float *)malloc(M*sizeof(float));
	float *a=(float *)malloc(M*sizeof(float));
	float *c=(float *)malloc(M*sizeof(float));

	//Initialising the value of T at interior nodes to 0
	for(i=1;i<N-1;i++)
	for(j=1;j<M-1;j++)
	T[i][j][0]=0;

	//Application of BCs at the 0th iteration level	
	for(i=0;i<N;i++)
		{
			T[i][0][0]=1;
			T[i][M-1][0]=1;
		}
		for(j=0;j<M;j++)
		{
			T[0][j][0]=0;
			T[N-1][j][0]=1;
		}

	begin=clock();//calculation of execution time starts
	for(k=2;;k=k+2)
	{
	
		for(i=0;i<N;i++)
		{
			T[i][0][k]=1;
			T[i][M-1][k]=1;
			T[i][0][k-1]=1;
			T[i][M-1][k-1]=1;
		}
		for(j=0;j<M;j++)
		{
			T[0][j][k]=0;
			T[N-1][j][k]=1;
			T[0][j][k-1]=0;
			T[N-1][j][k-1]=1;
		}
		//X sweep
		for(i=1;i<N-1;i++)
		{
			A[0]=0;
			A[M-1]=0;
			C[0]=1;
			C[M-1]=1;
			for(j=1;j<M-1;j++)
			{
				B[j]=-1;
				D[j]=-4;
				a[j]=-1;
				c[j]=-(T[i+1][j][k-2]+T[i-1][j][k-1]);
				A[j]=(a[j])/(D[j]-B[j]*A[j-1]);
				C[j]=(B[j]*C[j-1]+c[j])/(D[j]-B[j]*A[j-1]);
			}
			for(j=M-2;j>0;j--)
			{
				T[i][j][k-1]=A[j]*T[i][j+1][k-1]+C[j];
			}
		}
		
		
		//Y-sweep	
		for(j=1;j<M-1;j++)
		{
			A[0]=0;
			A[N-1]=0;
			C[0]=0;
			C[N-1]=1;
			for(i=1;i<N-1;i++)
			{
				B[i]=-1;
				D[i]=-4;
				a[i]=-1;
				c[i]=-(T[i][j+1][k-1]+T[i][j-1][k]);
				A[i]=(a[i])/(D[i]-B[i]*A[i-1]);
				C[i]=(B[i]*C[i-1]+c[i])/(D[i]-B[i]*A[i-1]);
			}
			for(i=N-2;i>0;i--)
			{
				T[i][j][k]=A[i]*T[i+1][j][k]+C[i];
			}
		}
		
		
	count++;
	printf("\nIteration Number:%d\t",count);
		//Calcutaion of error
		sum=0;
		for(i=1;i<N-1;i++)
		{
			for(j=1;j<M-1;j++)
			{
				sum+=(T[i][j][k]-T[i][j][k-2]);
			}
		}
		e=sqrt(((sum*sum)/(M-1)*(N-1)));
		printf("Error:%f\n",e);
		if(e<0.000001)
		break;
}
end=clock();//Calculation of execution time stops
float time_spent=(end-begin)/CLOCKS_PER_SEC;
	printf("\nTime Required is:%f",time_spent);
//File printing of output		
for(y=0,i=N-1;i>-1;i--,y=y+0.01)
		{
			for(x=0,j=0;j<M;j++,x=x+0.01)
			{
				fprintf(f,"\n%f\t%f\t%f",x,y,T[i][j][count-1]);
			}
		fprintf(f,"\n");
		}

}
