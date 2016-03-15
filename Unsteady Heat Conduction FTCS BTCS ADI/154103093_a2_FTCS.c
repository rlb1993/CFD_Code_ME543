#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
int main()
{
	clock_t begin,end;
	FILE *f,*fp;
	f=fopen("Temp_data_FTCS.txt","w");
	fp=fopen("ErrorVsTime_FTCS.txt","w");
	int count=0,M, N ;
	float L=1, H=1, delta_x, delta_y, b,e=1,sum=0,x,y,t=0;
	printf("\nEnter the number of grid points along x axis:");
	scanf("%d",&M);
	printf("\nEnter the number of grid points along y axis:");
	scanf("%d",&N);
	delta_x=L/(M-1);
	printf("\n%f",delta_x);
	delta_y=H/(N-1);
	printf("\n%f",delta_y);
	int i,j,k;
	float ***T=(float ***)malloc(N*sizeof(float));
	for(i=0;i<N;i++)
	T[i]=(float **)malloc(M*sizeof(float));
	for(i=0;i<N;i++)
	for(j=0;j<M;j++)
	T[i][j]=(float *)malloc(50000*sizeof(float));
	//Initialising the value of T at interior nodes to 0
	for(i=1;i<N-1;i++)
	for(j=1;j<M-1;j++)
	T[i][j][0]=10;
	for(i=0;i<N;i++)
		{
			T[i][M-1][0]=50;
		}
		for(j=0;j<M;j++)
		{
			T[0][j][0]=50;
			T[N-1][j][0]=50;
		}
	for(i=0;i<N;i++)
	{
		T[i][0][0]=0.989154*T[i][1][0]+0.010864*10;
	}
		begin=clock();
	for(k=1;;k++)
	{
	
		for(i=0;i<N;i++)
		{
			T[i][M-1][k]=50;
		}
		for(j=0;j<M;j++)
		{
			T[0][j][k]=50;
			T[N-1][j][k]=50;
		}
		for(i=0;i<N;i++)
		{
		T[i][0][k]=0.989154*T[i][1][k-1]+0.010864*10;
		}
		
	
		for(i=1;i<N-1;i++)
		{
			for(j=1;j<M-1;j++)
			{
				T[i][j][k]=0.000064*T[i][j][k-1]+0.249984*(T[i+1][j][k-1]+T[i-1][j][k-1]+T[i][j+1][k-1]+T[i][j-1][k-1]);		
			}
		}
			
		count++;
		fprintf(fp,"\n%f\t",(t+0.62*count));
		
		sum=0;
		for(i=1;i<N-1;i++)
		{
			for(j=1;j<M-1;j++)
			{
				sum+=(T[i][j][k]-T[i][j][k-1]);
			}
		}
		e=sqrt(((sum*sum)/(M-1)*(N-1)));
		fprintf(fp,"%f",e);
		if(e<0.0001)
		break;
	}
	end=clock();
	float time_spent=(end-begin)/CLOCKS_PER_SEC;
	printf("\nTime Required is:%f",time_spent);
	for(y=0,i=N-1;i>-1;i--,y=y+delta_y)
	{
		for(x=0,j=0;j<M;j++,x=x+delta_x)
		{
			fprintf(f,"\n%f\t%f\t%f",x,y,T[i][j][count-1]);
		}
		//fprintf(f,"\n");
	}
}