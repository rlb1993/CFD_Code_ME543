#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
int main()
{
	clock_t begin,end;
	FILE *f,*fp;
	f=fopen("Temp_data_a2_ADI.txt","w");
	fp=fopen("ErrorVsTime_ADI.txt","w");
	int count=0,M, N ;
	float L=1, H=1, delta_x, delta_y, b,e=1,sum=0,x,y,t=0,as,ks,h,delta_t,Yy,Yx;
	printf("\nEnter the number of grid points along x axis:");
	scanf("%d",&M);
	printf("\nEnter the number of grid points along y axis:");
	scanf("%d",&N);
	printf("\nEnter the value of Alpha_s: ");
	scanf("%f",&as);
	printf("\nEnter the value of K_s: ");
	scanf("%f",&ks);
	printf("enter the value of h: ");
	scanf("%f",&h);
	delta_x=L/(M-1);
	printf("\n%f",delta_x);
	delta_y=H/(N-1);
	printf("\n%f",delta_y);
	b=delta_x/delta_y;
	printf("\n%f",b);
	float p,q;
	p=ks/(ks+h*delta_x);
	q=(h*delta_x)/(ks+h*delta_x);
	delta_t=(1/(2*as))*((delta_x*delta_x*delta_y*delta_y)/((delta_x*delta_x)+(delta_y*delta_y)));
	Yx=(as*delta_t)/(delta_x*delta_x);
	Yy=(as*delta_t)/(delta_y*delta_y);

	int i,j,k;
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
		T[i][0][0]=p*T[i][1][0]+q*10;
	}
begin=clock();
	for(k=2;;k=k+2)
	{
	
		for(i=0;i<N;i++)
		{
			T[i][M-1][k]=50;
			T[i][M-1][k-1]=50;
		}
		for(i=0;i<N;i++)
		{
		T[i][0][k-1]=p*T[i][1][k-2]+q*10;
		}
		
		for(j=0;j<M;j++)
		{
			T[0][j][k]=50;
			T[N-1][j][k]=50;
			T[0][j][k-1]=50;
			T[N-1][j][k-1]=50;
		}
		//X sweep
		for(i=1;i<N-1;i++)
		{
			A[0]=0;
			A[M-1]=0;
			C[0]=T[i][0][k-1];
			C[M-1]=50;
			for(j=1;j<M-1;j++)
			{
				B[j]=-0.5*Yx;
				D[j]=-(1+Yx);
				a[j]=-0.5*Yx;
				c[j]=-((0.5*Yx)*(T[i+1][j][k-2]+T[i-1][j][k-2])+(1-Yx)*T[i][j][k-2]);
				A[j]=(a[j])/(D[j]-B[j]*A[j-1]);
				C[j]=(B[j]*C[j-1]+c[j])/(D[j]-B[j]*A[j-1]);
			}
			for(j=M-2;j>0;j--)
			{
				T[i][j][k-1]=A[j]*T[i][j+1][k-1]+C[j];
			}
		}
		for(i=0;i<N;i++)
		{
		T[i][0][k]=p*T[i][1][k-1]+q*10;
		}
		
		//Y-sweep	
		for(j=1;j<M-1;j++)
		{
			A[0]=0;
			A[N-1]=0;
			C[0]=50;
			C[N-1]=50;
			for(i=1;i<N-1;i++)
			{
				B[i]=-0.5*Yy;
				D[i]=-(1+Yy);
				a[i]=-0.5*Yy;
				c[i]=-((0.5*Yy)*(T[i][j+1][k-1]+T[i][j-1][k-1])+(1-Yy)*T[i][j][k-1]);
				A[i]=(a[i])/(D[i]-B[i]*A[i-1]);
				C[i]=(B[i]*C[i-1]+c[i])/(D[i]-B[i]*A[i-1]);
			}
			for(i=N-2;i>0;i--)
			{
				T[i][j][k]=A[i]*T[i+1][j][k]+C[i];
			}
		}
		
		
	count++;
	fprintf(fp,"\n%f\t",(t+0.62*count));
		
		sum=0;
		for(i=1;i<N-1;i++)
		{
			for(j=1;j<M-1;j++)
			{
				sum+=(T[i][j][k]-T[i][j][k-2]);
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
		fprintf(f,"\n");
		}

}
