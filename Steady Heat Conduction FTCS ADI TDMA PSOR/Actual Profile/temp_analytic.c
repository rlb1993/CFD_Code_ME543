#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int main()
{
FILE *f;
f=fopen("Analytic_temp.txt","w");
int i=0,j=0;
double x,y,sum,product,a,b,n,c;
double T[102][102];
for(x=0,i=0;x<1.01;x=x+0.01,i++)
{
	for(y=0,j=0;y<1.01;y=y+0.01,j++) 
	{
		sum=0;
		for(n=1;n<150;n=n+2)
		{
			a=sinh(n*3.142857*y)/sinh(n*3.142857);
			b=sin(n*3.142857*x);
			c=2/(n*3.142857);
			product=a*b*c;
			sum=sum+product;
		}
		T[i][j]=1-2*sum;
	}
}
for(i=0;i<101;i++)
{
	for(j=100;j>-1;j--)
	{
		fprintf(f,"%lf\t",T[j][i]);
	}
	fprintf(f,"\n");
}
}