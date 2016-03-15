#include<stdio.h>
#include<math.h>
int main()
{
	FILE *fp;
	fp=fopen("test_data.txt","w");
	int i=0,j=0;
	float x,y,n,sum;
	float T[101][101];
	for(x=0;x<1.01;x=x+0.01)
		{
			for(y=0;y<1.01;y=y+0.01)
			{
				sum=0;
			 	for(n=1;n<100;n=n+2)
				sum+=(2/(n*(22/7))*(sinh(n*y*(22/7))/sinh(n*(22/7)))*sin(n*x*(22/7)));
				T[i][j]=sum;
				j++;
			}
			i++;
		}
	for(i=0;i<100;i++)
	{	
		for(j=0;j<100;j++)
		{
			fprintf(fp,"%f ",T[i][j]);
		}
	fprintf(fp,"\n");
	}
}