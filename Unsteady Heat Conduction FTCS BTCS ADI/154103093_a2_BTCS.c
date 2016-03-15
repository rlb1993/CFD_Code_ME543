#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
int main()
{
	clock_t begin,end;
	FILE *f, *fp;
	f=fopen("Temp_data_BTCS.txt","w");
	fp=fopen("error_vs_t_btcs.txt","w");
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
	b=delta_x/delta_y;
	
	
	int i,j,k,l,n,mark=0;
	
	float ***Ttime=(float ***)malloc(N*sizeof(float));
	for(i=0;i<N;i++)
	Ttime[i]=(float **)malloc(M*sizeof(float));
	for(i=0;i<N;i++)
	for(j=0;j<M;j++)
	Ttime[i][j]=(float *)malloc(20000*sizeof(float));

	float ***Tspace=(float ***)malloc(N*sizeof(float));
	for(i=0;i<N;i++)
	Tspace[i]=(float **)malloc(M*sizeof(float));
	for(i=0;i<N;i++)
	for(j=0;j<M;j++)
	Tspace[i][j]=(float *)malloc(20000*sizeof(float));


	//Initialising the value of T at interior nodes to 0
	for(i=1;i<N-1;i++)
		for(j=1;j<M-1;j++)
		{	
			Ttime[i][j][0]=10;
		}
	for(i=0;i<N;i++)
		{
			Ttime[i][M-1][0]=50;
			Tspace[i][M-1][0]=50;
		}
	for(j=0;j<M;j++)
		{
			Ttime[0][j][0]=50;
			Ttime[N-1][j][0]=50;
			Tspace[0][j][0]=50;
			Tspace[N-1][j][0]=50;
		}
	for(i=0;i<N;i++)
		{
			Ttime[i][0][0]=0.989154*Ttime[i][1][0]+0.010864*10;
			Tspace[i][0][0]=0.989154*Tspace[i][1][0]+0.010864*10;
		}

	begin=clock();
	for(n=1;;n++)//time loop begins
	{
		for(k=1;;k++)
		{
			
			for(i=1;i<N-1;i++)
			for(j=1;j<M-1;j++)
			{	
				
				Tspace[i][j][k]=Ttime[i][j][n-1];
			}
			for(i=0;i<N;i++)
			{
				Tspace[i][M-1][k]=50;
			}
			for(j=0;j<M;j++)
			{
				Tspace[0][j][k]=50;
				Tspace[N-1][j][k]=50;
			}	
			for(i=0;i<N;i++)
			{
				Tspace[i][0][k]=0.989154*Tspace[i][1][k-1]+0.010864*10;
			}
		
			for(i=1;i<N-1;i++)
			{
				for(j=1;j<M-1;j++)
				{
					Tspace[i][j][k]=(Ttime[i][j][n-1]+ 0.249984*(Tspace[i+1][j][k-1]+Tspace[i-1][j][k]+Tspace[i][j+1][k-1]+Tspace[i][j-1][k]))/(1+4*0.249984);		
				}
			}
		
			sum=0;
			for(i=1;i<N-1;i++)
			{
				for(j=1;j<M-1;j++)
				{
					sum+=(Tspace[i][j][k]-Tspace[i][j][k-1]);
				}
			}
			e=sqrt(((sum*sum)/(M-1)*(N-1)));
			//fprintf(fp,"\n%f",e);
			if(e<0.000001)
			break;
		}

		for(i=1;i<N-1;i++)
			for(j=1;j<M-1;j++)
				Ttime[i][j][n]=Tspace[i][j][k-1];
		
		count++;
		fprintf(fp,"\n%f\t",(t+0.62*count));
		sum=0;
		for(i=1;i<N-1;i++)
		{
			for(j=1;j<M-1;j++)
			{
				sum+=(Ttime[i][j][n]-Ttime[i][j][n-1]);
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
			fprintf(f,"\n%f\t%f\t%f",x,y,Tspace[i][j][k-1]);
		}
		//fprintf(f,"\n");
	}
}
