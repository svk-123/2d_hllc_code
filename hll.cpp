 #include<iostream>
 #include<fstream>
 #include<stdio.h>
 #include<cmath>
 using namespace std;
 /* *******************************************************************************************************************************************
  
									x4,y4________3_______x3,y3
  										|	  			|
  									   4|				|2
  										|_______________|
									x1,y1		1     x2,y2
    
  // mdified Ul & Ur
 ********************************************************************************************************************************************* */

float R=287.0,g=1.4;
int i,j,k;
int i_1=1,i_2=75;
int j_1=1,j_2=30;
int k_1=1,k_2=4;
float xl=3.0,yl=1.0,ni=75.0,nj=30.0;

float L=0.1;
int xb1=26,xb2=50;

 struct cell
 {
	 
	 float x,y;
	 float x1,x2,x3,x4,y1,y2,y3,y4;
	 float dx[5],dy[5],ds[5],A,cx,cy,nx[5],ny[5];
	 float U1,U2,U3,U4,Ul1[5],Ul2[5],Ul3[5],Ul4[5],Ur1[5],Ur2[5],Ur3[5],Ur4[5];
	 float p,r,u,v,a,T,M,rold,R2,z[5],delt;
	 float un[5],ut[5],unl[5],unr[5],utl[5],utr[5],sl[5],sr[5],ss[5],pl[5],rl[5],ul[5],al[5],vl[5],pr[5],rr[5],ur[5],ar[5],vr[5];
	 float F1[5],F2[5],F3[5],F4[5],Fl1[5],Fl2[5],Fl3[5],Fl4[5],Fr1[5],Fr2[5],Fr3[5],Fr4[5],Fs1[5],Fs2[5],Fs3[5],Fs4[5];
	 float f1[5],f2[5],f3[5],f4[5],netflux[5];
	  
 };
 
float nx_b1,ny_b1,nx_b2,ny_b2,nx_b3,ny_b3,nx_b4,ny_b4;
float un_b1,ut_b1,un_b2,ut_b2,un_b3,ut_b3,un_b4,ut_b4; 
 
 void plt ();
 void plt1 ();
 void plt2 ();
 
int main()
{
	cell c[80][35];
	
	char filename[20]="test4.txt";
	char filename1[20]="cpwall.txt";
	char filename2[20]="wallmachno.txt";	
	
	FILE *fp;
	FILE *fp1;
	FILE *fp2;
	
	fp=fopen(filename,"w");
	fp1=fopen(filename1,"w");
	fp2=fopen(filename2,"w");
	
	int w;
	float dt=1.0,rss1,rss2,rss3=1;
	float pres,rres,ures,ares,R1,po,To,r1,r2;
	float cpl[101];
	
	
	float pinf=1.0,rinf=1.4,ainf=1.0,Minf=0.675,uinf=1.0,vinf=0.0;
	
//......................................................................	
pres=1.356; // Inlet condition 
ares=1.04;
ures=0.0;
R1=ures+(2*ares/(g-1.0));
po=1.0;	
//......................................................................
	for (i=i_1;i<=i_2;i++)// 	grid points 1
	{
		for (j=j_1;j<=j_2+1;j++)
			{
			c[i_1][j].x=0;
			c[i+1][j].x=(xl/ni)+c[i][j].x;
			
		}
		}
		
	for (i=i_1;i<=i_2+1;i++)// 	grid points 2
	{
		for (j=j_1;j<=j_2;j++)
			{

			c[i][j_1].y=0;
    		c[i][j+1].y=(yl/nj)+c[i][j].y;
		}
			
		}
		
			
	for (i=xb1;i<=xb2;i++)// 	bump
	{
		for (j=j_1;j<=j_2;j++)
			{
				
				c[i][j_1].y=sqrt((1.3*1.3)-(c[i][j_1].x-1.0-0.5)*(c[i][j_1].x-1.0-0.5))-1.2;
				c[i][j+1].y=(yl-c[i][j_1].y)/(nj)+c[i][j].y;
			
			}
		}
	
	
		
				
	for (i=i_1;i<=i_2+1;i++) // cell vertices......
	{
		for (j=j_1;j<=j_2+1;j++)
			{
				
			c[i][j].x1=c[i][j].x;
			c[i][j].y1=c[i][j].y;
			
			c[i][j].x2=c[i+1][j].x;
			c[i][j].y2=c[i+1][j].y;	
					
			c[i][j].x3=c[i+1][j+1].x;
			c[i][j].y3=c[i+1][j+1].y;
					
			c[i][j].x4=c[i][j+1].x;
			c[i][j].y4=c[i][j+1].y;
					
			}
			
		}
			
	for (i=i_1;i<=i_2+1;i++) // dx,dy,ds & A calculation...
	{
		for (j=j_1;j<=j_2+1;j++)
			{		
		
			c[i][j].dx[1]=c[i][j].x2-c[i][j].x1;
			c[i][j].dy[1]=c[i][j].y2-c[i][j].y1;
			
			c[i][j].dx[2]=c[i][j].x3-c[i][j].x2;
			c[i][j].dy[2]=c[i][j].y3-c[i][j].y2;
			
			c[i][j].dx[3]=c[i][j].x4-c[i][j].x3;
			c[i][j].dy[3]=c[i][j].y4-c[i][j].y3;
			
			c[i][j].dx[4]=c[i][j].x1-c[i][j].x4;
			c[i][j].dy[4]=c[i][j].y1-c[i][j].y4;
				
			c[i][j].ds[1]=sqrt(c[i][j].dx[1]*c[i][j].dx[1]+c[i][j].dy[1]*c[i][j].dy[1]);
			c[i][j].ds[2]=sqrt(c[i][j].dx[2]*c[i][j].dx[2]+c[i][j].dy[2]*c[i][j].dy[2]);
			c[i][j].ds[3]=sqrt(c[i][j].dx[3]*c[i][j].dx[3]+c[i][j].dy[3]*c[i][j].dy[3]);
			c[i][j].ds[4]=sqrt(c[i][j].dx[4]*c[i][j].dx[4]+c[i][j].dy[4]*c[i][j].dy[4]);
			
			c[i][j].cx=0.25*(c[i][j].x1+c[i][j].x2+c[i][j].x3+c[i][j].x4);
			c[i][j].cy=0.25*(c[i][j].y1+c[i][j].y2+c[i][j].y3+c[i][j].y4);
			
			c[i][j].A=0.5*((c[i][j].x2-c[i][j].x4)*(c[i][j].y3-c[i][j].y1)-(c[i][j].x3-c[i][j].x1)*(c[i][j].y2-c[i][j].y4));
			
	
			}
	
	}
	
		for (i=i_1-1;i<=i_2+1;i++) //ghost cells...
		{
			for (j=j_1-1;j<=j_2+1;j++)
			{
					for (k=k_1;k<=k_2;k++)
				{
						c[0][j].dx[k]=c[1][j].dx[k];
						c[0][j].dy[k]=c[1][j].dy[k];
						c[0][j].ds[k]=c[1][j].ds[k];
	
						c[i_2+1][j].dx[k]=c[i_2][j].dx[k];
						c[i_2+1][j].dy[k]=c[i_2][j].dy[k];
						c[i_2+1][j].ds[k]=c[i_2][j].ds[k];
						
						c[i][0].dx[k]=c[i][1].dx[k];
						c[i][0].dy[k]=c[i][1].dy[k];
						c[i][0].ds[k]=c[i][1].ds[k];
						
						c[i][j_2+1].dx[k]=c[i][j_2].dx[k];
						c[i][j_2+1].dy[k]=c[i][j_2].dy[k];
						c[i][j_2+1].ds[k]=c[i][j_2].ds[k];
										
				}
				
			}
		}
		
//......................................................................		
		for (i=i_1-1;i<=i_2+1;i++) //nx,ny calculation...
		{
			for (j=j_1-1;j<=j_2+1;j++)
			{
					for (k=k_1;k<=k_2;k++)
				{
					c[i][j].nx[k]=c[i][j].dy[k]/c[i][j].ds[k];
					c[i][j].ny[k]=-c[i][j].dx[k]/c[i][j].ds[k];

										
				}
				
			}
		}
//......................................................................		
	for (i=i_1-1;i<=i_2+1;i++) // Initial condition...
	{
	for (j=j_1-1;j<=j_2+1;j++)
	{
	c[i][j].p=pinf;
	c[i][j].r=rinf;
	c[i][j].a=ainf;
	c[i][j].v=vinf;
	c[i][j].u=Minf*ainf;
			
	c[i][j].U1=c[i][j].r;
	c[i][j].U2=c[i][j].r*c[i][j].u;
	c[i][j].U3=c[i][j].r*c[i][j].v;
	c[i][j].U4=(c[i][j].p/(g-1))+(0.5*c[i][j].r*(c[i][j].u*c[i][j].u+c[i][j].v*c[i][j].v));
	}
	}
//......................................................................
	for (w=1;rss3>0.000001;w++)//........................//main time loop
	{
		rss2=0;
		
	for (i=i_1;i<=i_2;i++) // BC...
	{
	for (j=j_1;j<=j_2;j++)
	{
	nx_b1=1;// at entry(left)	
	ny_b1=0;
	
	un_b1=c[1][j].u*nx_b1+c[1][j].v*ny_b1;	
	ut_b1=-c[1][j].u*ny_b1+c[1][j].v*nx_b1;	
	
	un_b1=un_b1;
	ut_b1=-ut_b1;
		
	c[0][j].R2=un_b1-(2*c[1][j].a/(g-1.0));
	
	c[0][j].un[2]=0.5*(R1+c[0][j].R2);
	c[0][j].ut[2]=0.0;					
	c[0][j].a=0.1*(R1-c[0][j].R2);
		
	c[0][j].u=c[0][j].un[2]*nx_b1-c[0][j].ut[2]*ny_b1;
	c[0][j].v=c[0][j].un[2]*ny_b1+c[0][j].ut[2]*nx_b1;
	
	c[0][j].M=sqrt(c[0][j].u*c[0][j].u+c[0][j].v*c[0][j].v)/c[0][j].a;
	
	c[0][j].p=pres/pow((1+(0.2*c[0][j].M*c[0][j].M)),3.5);
	c[0][j].r=g*c[0][j].p/(c[0][j].a*c[0][j].a);
	
	
	nx_b2=1.0;//  right
	ny_b2=0;
	
	un_b2=c[i_2][j].u*nx_b2+c[i_2][j].v*ny_b2;	
	ut_b2=-c[i_2][j].u*ny_b2+c[i_2][j].v*nx_b2;
	
	un_b2=-un_b2;
	ut_b2=-ut_b2;
	
	//r1=-un_b2+(2*c[i_2][j].a/(g-1));
	//r2=uinf+(2*ainf/(g-1));
	
	//c[i_2+1][j].un[4]=(r1+r2)/2.0;
	

	c[i_2+1][j].un[4]=un_b2;
    c[i_2+1][j].ut[4]=ut_b2;
    
    //c[i_2+1][j].un[4]=un_b2+5*(c[i_2][j].a-c[i_2+1][j].a);
	//c[i_2+1][j].ut[4]=ut_b2;
		
	c[i_2+1][j].p=po;
	c[i_2+1][j].r=c[i_2][j].r*pow((c[i_2+1][j].p/c[i_2][j].p),0.714);
	c[i_2+1][j].a=sqrt(g*c[i_2+1][j].p/c[i_2+1][j].r);
	
		
	c[i_2+1][j].u=c[i_2+1][j].un[4]*nx_b2-c[i_2+1][j].ut[4]*ny_b2;
	c[i_2+1][j].v=c[i_2+1][j].un[4]*ny_b2+c[i_2+1][j].ut[4]*nx_b2;
	
	
	nx_b3=0;// at top
	ny_b3=1.0;
	
	un_b3=c[i][j_2].u*nx_b3+c[i][j_2].v*ny_b3;	
	ut_b3=-c[i][j_2].u*ny_b3+c[i][j_2].v*nx_b3;
	
	un_b3=un_b3;
	ut_b3=-ut_b3;
	
	c[i][j_2+1].un[1]=un_b3;
	c[i][j_2+1].ut[1]=ut_b3;

	c[i][j_2+1].p=c[i][j_2].p;
	c[i][j_2+1].r=c[i][j_2].r;
	c[i][j_2+1].a=c[i][j_2].a;
	
	c[i][j_2+1].u=c[i][j_2+1].un[1]*nx_b3-c[i][j_2+1].ut[1]*ny_b3;
	c[i][j_2+1].v=c[i][j_2+1].un[1]*ny_b3+c[i][j_2+1].ut[1]*nx_b3;
	
	
	nx_b4=c[i][1].nx[1];// at bottom
	ny_b4=c[i][1].ny[1];
	
	un_b4=c[i][1].u*nx_b4+c[i][1].v*ny_b4;	
	ut_b4=-c[i][1].u*ny_b4+c[i][1].v*nx_b4;
	
	un_b4=un_b4;
	ut_b4=-ut_b4;
	
	c[i][0].un[3]=un_b4;
	c[i][0].ut[3]=ut_b4;
	
	c[i][0].p=c[i][1].p;
	c[i][0].r=c[i][1].r;
	c[i][0].a=c[i][1].a;

	c[i][0].u=c[i][0].un[3]*nx_b4-c[i][0].ut[3]*ny_b4;
	c[i][0].v=c[i][0].un[3]*ny_b4+c[i][0].ut[3]*nx_b4;
	
	}
	}
		printf("%d\t%f\t%f\t%f\n",w,dt,c[0][15].M,rss3);
		
		//printf("%f\t%f\t%f\t%f\n",c[30][j_2+1].un[1],c[0][1].un[2],c[30][0].un[3],c[i_2+1][1].un[4]);
		//printf("%f\t%f\t%f\t%f\n",c[30][j_2+1].ut[1],c[0][1].ut[2],c[30][0].ut[3],c[i_2+1][1].ut[4]);
		
//......................................................................

for (i=i_1;i<=i_2;i++) // Normal & Tangential velocity calculation...
	{
	for (j=j_1;j<=j_2;j++)
	{
	for (k=k_1;k<=k_2;k++)
	{
		c[i][j].un[k]=c[i][j].u*c[i][j].nx[k]+c[i][j].v*c[i][j].ny[k];
		c[i][j].ut[k]=-c[i][j].u*c[i][j].ny[k]+c[i][j].v*c[i][j].nx[k];
	}
	}	
	}

//......................................................................

for (i=i_1;i<=i_2;i++) //	Time step  calculation ...
	{
	for (j=j_1;j<=j_2;j++)
	{
						
		c[i][j].z[1]=(fabs(c[i][j].un[1])+fabs(c[i][j].a))*c[i][j].ds[1];
		c[i][j].z[2]=(fabs(c[i][j].un[2])+fabs(c[i][j].a))*c[i][j].ds[2];
		c[i][j].z[3]=(fabs(c[i][j].un[3])+fabs(c[i][j].a))*c[i][j].ds[3];
		c[i][j].z[4]=(fabs(c[i][j].un[4])+fabs(c[i][j].a))*c[i][j].ds[4];
		
		c[i][j].delt=0.9*2*c[i][j].A/(c[i][j].z[1]+c[i][j].z[2]+c[i][j].z[3]+c[i][j].z[4]);
		
			
		if (dt>c[i][j].delt)
		{
			dt=c[i][j].delt;
		}
		
		}
}
//......................................................................

for (i=i_1;i<=i_2;i++) // left & right state calculation (u,v,a,p,r)...
	{
	for (j=j_1;j<=j_2;j++)
	{
c[i][j].ul[1]=c[i][j].u;
c[i][j].vl[1]=c[i][j].v;
c[i][j].rl[1]=c[i][j].r;
c[i][j].pl[1]=c[i][j].p;
c[i][j].al[1]=c[i][j].a;

c[i][j].ur[1]=c[i][j-1].u;
c[i][j].vr[1]=c[i][j-1].v;
c[i][j].rr[1]=c[i][j-1].r;
c[i][j].pr[1]=c[i][j-1].p;
c[i][j].ar[1]=c[i][j-1].a;

c[i][j].ul[2]=c[i][j].u;
c[i][j].vl[2]=c[i][j].v;
c[i][j].rl[2]=c[i][j].r;
c[i][j].pl[2]=c[i][j].p;
c[i][j].al[2]=c[i][j].a;

c[i][j].ur[2]=c[i+1][j].u;
c[i][j].vr[2]=c[i+1][j].v;
c[i][j].rr[2]=c[i+1][j].r;
c[i][j].pr[2]=c[i+1][j].p;
c[i][j].ar[2]=c[i+1][j].a;

c[i][j].ul[3]=c[i][j].u;
c[i][j].vl[3]=c[i][j].v;
c[i][j].rl[3]=c[i][j].r;
c[i][j].pl[3]=c[i][j].p;
c[i][j].al[3]=c[i][j].a;

c[i][j].ur[3]=c[i][j+1].u;
c[i][j].vr[3]=c[i][j+1].v;
c[i][j].rr[3]=c[i][j+1].r;
c[i][j].pr[3]=c[i][j+1].p;
c[i][j].ar[3]=c[i][j+1].a;

c[i][j].ul[4]=c[i][j].u;
c[i][j].vl[4]=c[i][j].v;
c[i][j].rl[4]=c[i][j].r;
c[i][j].pl[4]=c[i][j].p;
c[i][j].al[4]=c[i][j].a;

c[i][j].ur[4]=c[i-1][j].u;
c[i][j].vr[4]=c[i-1][j].v;
c[i][j].rr[4]=c[i-1][j].r;
c[i][j].pr[4]=c[i-1][j].p;
c[i][j].ar[4]=c[i-1][j].a;
}
}
//......................................................................

for (i=i_1;i<=i_2;i++)//.........left and right states (un & ut)
	{
	for (j=j_1;j<=j_2;j++)
	{
c[i][j].unl[1]=c[i][j].un[1];
c[i][j].unr[1]=-c[i][j-1].un[3];
c[i][j].unl[2]=c[i][j].un[2];
c[i][j].unr[2]=-c[i+1][j].un[4];
c[i][j].unl[3]=c[i][j].un[3];
c[i][j].unr[3]=-c[i][j+1].un[1];
c[i][j].unl[4]=c[i][j].un[4];
c[i][j].unr[4]=-c[i-1][j].un[2];

c[i][j].utl[1]=c[i][j].ut[1];
c[i][j].utr[1]=-c[i][j-1].ut[3];
c[i][j].utl[2]=c[i][j].ut[2];
c[i][j].utr[2]=-c[i+1][j].ut[4];
c[i][j].utl[3]=c[i][j].ut[3];
c[i][j].utr[3]=-c[i][j+1].ut[1];
c[i][j].utl[4]=c[i][j].ut[4];
c[i][j].utr[4]=-c[i-1][j].ut[2];
}
}

//......................................................................
for (i=i_1;i<=i_2;i++) // sl & sr calculation...
	{
	for (j=j_1;j<=j_2;j++)
	{
	for (k=k_1;k<=k_2;k++)
	{
		
c[i][j].sl[k]=c[i][j].unl[k]-c[i][j].al[k];
if((c[i][j].unr[k]-c[i][j].ar[k])<(c[i][j].unl[k]-c[i][j].al[k]))
	{
c[i][j].sl[k]=c[i][j].unr[k]-c[i][j].ar[k];
	}

c[i][j].sr[k]=c[i][j].unl[k]+c[i][j].al[k];
if((c[i][j].unr[k]+c[i][j].ar[k])>(c[i][j].unl[k]+c[i][j].al[k]))
	{
c[i][j].sr[k]=c[i][j].unr[k]+c[i][j].ar[k];
}

}
}
}
//......................................................................

for (i=i_1;i<=i_2;i++) // ... Ul & Ur Calculation
	{
	for (j=j_1;j<=j_2;j++)
	{
	for (k=k_1;k<=k_2;k++)
	{

c[i][j].Ul1[k]=c[i][j].rl[k];
c[i][j].Ul2[k]=c[i][j].rl[k]*c[i][j].unl[k];
c[i][j].Ul3[k]=c[i][j].rl[k]*c[i][j].utl[k];
c[i][j].Ul4[k]=(c[i][j].pl[k]/(g-1))+(0.5*c[i][j].rl[k]*(c[i][j].ul[k]*c[i][j].ul[k]+c[i][j].vl[k]*c[i][j].vl[k]));

c[i][j].Ur1[k]=c[i][j].rr[k];
c[i][j].Ur2[k]=c[i][j].rr[k]*c[i][j].unr[k];
c[i][j].Ur3[k]=c[i][j].rr[k]*c[i][j].utr[k];
c[i][j].Ur4[k]=(c[i][j].pr[k]/(g-1))+(0.5*c[i][j].rr[k]*(c[i][j].ur[k]*c[i][j].ur[k]+c[i][j].vr[k]*c[i][j].vr[k]));  
	}
	}
	}
//......................................................................
for (i=i_1;i<=i_2;i++) // Fl and Fr calculation...
	{
	for (j=j_1;j<=j_2;j++)
	{
	for (k=k_1;k<=k_2;k++)
	{
c[i][j].Fl1[k]=c[i][j].rl[k]*c[i][j].unl[k];
c[i][j].Fl2[k]=c[i][j].rl[k]*c[i][j].unl[k]*c[i][j].unl[k]+c[i][j].pl[k];
c[i][j].Fl3[k]=c[i][j].rl[k]*c[i][j].utl[k]*c[i][j].unl[k];
c[i][j].Fl4[k]=c[i][j].unl[k]*(c[i][j].Ul4[k]+c[i][j].pl[k]);

c[i][j].Fr1[k]=c[i][j].rr[k]*c[i][j].unr[k];
c[i][j].Fr2[k]=c[i][j].rr[k]*c[i][j].unr[k]*c[i][j].unr[k]+c[i][j].pr[k];
c[i][j].Fr3[k]=c[i][j].rr[k]*c[i][j].utr[k]*c[i][j].unr[k];
c[i][j].Fr4[k]=c[i][j].unr[k]*(c[i][j].Ur4[k]+c[i][j].pr[k]);

c[i][j].Fs1[k]=(c[i][j].sr[k]*c[i][j].Fl1[k]-c[i][j].sl[k]*c[i][j].Fr1[k]+c[i][j].sl[k]*c[i][j].sr[k]*(c[i][j].Ur1[k]-c[i][j].Ul1[k]))/(c[i][j].sr[k]-c[i][j].sl[k]);
c[i][j].Fs2[k]=(c[i][j].sr[k]*c[i][j].Fl2[k]-c[i][j].sl[k]*c[i][j].Fr2[k]+c[i][j].sl[k]*c[i][j].sr[k]*(c[i][j].Ur2[k]-c[i][j].Ul2[k]))/(c[i][j].sr[k]-c[i][j].sl[k]);
c[i][j].Fs3[k]=(c[i][j].sr[k]*c[i][j].Fl3[k]-c[i][j].sl[k]*c[i][j].Fr3[k]+c[i][j].sl[k]*c[i][j].sr[k]*(c[i][j].Ur3[k]-c[i][j].Ul3[k]))/(c[i][j].sr[k]-c[i][j].sl[k]);
c[i][j].Fs4[k]=(c[i][j].sr[k]*c[i][j].Fl4[k]-c[i][j].sl[k]*c[i][j].Fr4[k]+c[i][j].sl[k]*c[i][j].sr[k]*(c[i][j].Ur4[k]-c[i][j].Ul4[k]))/(c[i][j].sr[k]-c[i][j].sl[k]);

}
}
}
//......................................................................
for (i=i_1;i<=i_2;i++) // hll condtion
	{
	for (j=j_1;j<=j_2;j++)
	{
	for (k=k_1;k<=k_2;k++)
	{
		
	if(c[i][j].sl[k]>0) 
		{ 
			c[i][j].F1[k]=c[i][j].Fl1[k];
			c[i][j].F2[k]=c[i][j].Fl2[k];
			c[i][j].F3[k]=c[i][j].Fl3[k];
			c[i][j].F4[k]=c[i][j].Fl4[k];
		}
		

		else if (c[i][j].sl[k]<=0 && c[i][j].sr[k]>=0)
		{
			c[i][j].F1[k]=c[i][j].Fs1[k];
			c[i][j].F2[k]=c[i][j].Fs2[k];
			c[i][j].F3[k]=c[i][j].Fs3[k];
			c[i][j].F4[k]=c[i][j].Fs4[k];
			
		}
		else if (c[i][j].sr[k]<0)
		{
			c[i][j].F1[k]=c[i][j].Fr1[k];
			c[i][j].F2[k]=c[i][j].Fr2[k];
			c[i][j].F3[k]=c[i][j].Fr3[k];
			c[i][j].F4[k]=c[i][j].Fr4[k];
		}
}
}
}
//......................................................................
for (i=i_1;i<=i_2;i++) // Inverse tranformation of fluxes...
	{
	for (j=j_1;j<=j_2;j++)
	{
	for (k=k_1;k<=k_2;k++)
	{
			c[i][j].f1[k]=c[i][j].F1[k];		
			c[i][j].f2[k]=c[i][j].F2[k]*c[i][j].nx[k]-c[i][j].F3[k]*c[i][j].ny[k];
			c[i][j].f3[k]=c[i][j].F2[k]*c[i][j].ny[k]+c[i][j].F3[k]*c[i][j].nx[k];
			c[i][j].f4[k]=c[i][j].F4[k];
		}
	}
}
//......................................................................
for (i=i_1;i<=i_2;i++) // netflux calculation...
	{
	for (j=j_1;j<=j_2;j++)
	{
	c[i][j].netflux[1]=c[i][j].f1[1]*c[i][j].ds[1]+c[i][j].f1[2]*c[i][j].ds[2]+c[i][j].f1[3]*c[i][j].ds[3]+c[i][j].f1[4]*c[i][j].ds[4];
	c[i][j].netflux[2]=c[i][j].f2[1]*c[i][j].ds[1]+c[i][j].f2[2]*c[i][j].ds[2]+c[i][j].f2[3]*c[i][j].ds[3]+c[i][j].f2[4]*c[i][j].ds[4];
	c[i][j].netflux[3]=c[i][j].f3[1]*c[i][j].ds[1]+c[i][j].f3[2]*c[i][j].ds[2]+c[i][j].f3[3]*c[i][j].ds[3]+c[i][j].f3[4]*c[i][j].ds[4];
	c[i][j].netflux[4]=c[i][j].f4[1]*c[i][j].ds[1]+c[i][j].f4[2]*c[i][j].ds[2]+c[i][j].f4[3]*c[i][j].ds[3]+c[i][j].f4[4]*c[i][j].ds[4];	
			
	}
}
//......................................................................
for (i=i_1;i<=i_2;i++) //  U update...
	{
	for (j=j_1;j<=j_2;j++)
	{
	
c[i][j].U1=c[i][j].U1-(dt/c[i][j].A)*c[i][j].netflux[1];
c[i][j].U2=c[i][j].U2-(dt/c[i][j].A)*c[i][j].netflux[2];
c[i][j].U3=c[i][j].U3-(dt/c[i][j].A)*c[i][j].netflux[3];
c[i][j].U4=c[i][j].U4-(dt/c[i][j].A)*c[i][j].netflux[4];	

}
}
	
//......................................................................

for (i=i_1;i<=i_2;i++) // p,r,u,v,a update...
	{
	for (j=j_1;j<=j_2;j++)
	{
	c[i][j].rold=c[i][j].r;
	c[i][j].u=c[i][j].U2/c[i][j].U1;
	c[i][j].v=(c[i][j].U3/c[i][j].U1);
	c[i][j].r=c[i][j].U1;
	c[i][j].p=(c[i][j].U4-(0.5*c[i][j].r*(c[i][j].u*c[i][j].u+c[i][j].v*c[i][j].v)))*(g-1);
	c[i][j].a=sqrt(g*fabs(c[i][j].p)/fabs(c[i][j].r));
	
}
}
//......................................................................
	for (i=i_1;i<=i_2;i++) // Residue condition...
	{
	for (j=j_1;j<=j_2;j++)
	{
	rss1=((c[i][j].r-c[i][j].rold)*(c[i][j].r-c[i][j].rold)/c[i][j].rold);
		rss2=rss1+rss2;
	}
	}
		rss3=sqrt(fabs(rss2/(ni*nj)));
		
} //time loop end
//......................................................................
	for (i=i_1;i<=i_2;i++) // Writing file 1...
	{
	for (j=j_1;j<=j_2;j++)
	{
	fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",c[i][j].cx,c[i][j].cy,c[i][j].p,c[i][j].u,c[i][j].v);
	}
	fprintf(fp,"\n");
	}
//......................................................................
	for (i=i_1;i<=i_2;i++) // Writing file for cp...
	{
	for (j=j_1;j<=j_2;j++)
	{
	cpl[i]=(c[i][1].p-pinf)/(0.5*rinf*uinf*uinf);
	fprintf(fp1,"%f\t%f\n",c[i][1].cx,cpl[i]);
	}
	}
//......................................................................
	for (i=i_1;i<=i_2;i++) //Writing file for M...
	{
	for (j=j_1;j<=j_2;j++)
	{
	c[i][1].M=(sqrt(c[i][1].u*c[i][1].u+c[i][1].v*c[i][1].v)/c[i][1].a);
	fprintf(fp2,"%f\t%f\n",c[i][1].cx,c[i][1].M);
	}
	}
//......................................................................
fclose(fp);		
fclose(fp1);
fclose(fp2);

plt();		
plt1();
plt2();

return 0;	
}


void plt ()// contour plot
	{
FILE *plot = popen("gnuplot -persist","w");
fprintf(plot, "set pm3d\n");
fprintf(plot, "unset surface\n");
fprintf(plot, "unset key\n");
fprintf(plot, "set view map\n");
fprintf(plot, "set title 'Mach no'\n");
fprintf(plot, "splot 'test4.txt' using 1:2:5 with lines\n");
fclose(plot);
	}

void plt1 ()// cp plot
	{
FILE *plot = popen("gnuplot -persist","w");
fprintf(plot, "plot 'cpwall.txt' using 1:2 with lines,\'hirsch.txt' using 1:2 with lines\n");
fclose(plot);
	}

void plt2 () // wall Mach no plot
	{
FILE *plot = popen("gnuplot -persist","w");
fprintf(plot, "plot 'wallmachno.txt' using 1:2 with lines,\'Wall_mach0675.csv' using 1:2 with lines title 'shai' ,\'trans.txt' using 1:2:(.01) with circles\n");
fclose(plot);
	}






















