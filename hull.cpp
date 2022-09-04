#include<stdio.h>
#include<string.h>
#include<math.h>
#include<algorithm>
using namespace std;
#define PR 1e-8
#include "myparam.h"


int imin_array(double a[], int L) {
   int i, imin=0;
   double amin=a[0];
   for(i=1;i<L;i++) if(a[i] < amin) { amin=a[i]; imin = i; }
   return imin;

}    

void cross_product_hull(double a[dd],double b[dd], double res[dd]) {
   res[0] = a[1]*b[2] - a[2]*b[1];
   res[1] = a[2]*b[0] - a[0]*b[2];
   res[2] = a[0]*b[1] - a[1]*b[0];
}



void cal_normals(double I[dd], double J[dd], double K[dd], double planeijk[4]) {
   int d;
   double  vij[dd], vik[dd], vtmp[dd], vabs;

   for(d=0;d<dd;d++) vij[d] = J[d] - I[d];
   for(d=0;d<dd;d++) vik[d] = K[d] - I[d];
   cross_product_hull(vij, vik, vtmp);
   vabs = sqrt( vtmp[0]*vtmp[0] + vtmp[1]*vtmp[1] + vtmp[2]*vtmp[2] );
   for(d=0;d<dd;d++) planeijk[d] = vtmp[d]/vabs;
   planeijk[3]  = 0.;
   for(d=0;d<dd;d++)  planeijk[3] += I[d]*planeijk[d];
   planeijk[3]  = -planeijk[3];

}



struct TPoint
{
    double x,y,z;
    TPoint(){}
    TPoint(double _x,double _y,double _z):x(_x),y(_y),z(_z){}
    TPoint operator-(const TPoint p) {return TPoint(x-p.x,y-p.y,z-p.z);}
    TPoint operator*(const TPoint p) {return TPoint(y*p.z-z*p.y,z*p.x-x*p.z,x*p.y-y*p.x);}//叉积
    double operator^(const TPoint p) {return x*p.x+y*p.y+z*p.z;}//点积
};
struct fac//
{
    int a,b,c;//凸包一个面上的三个点的编号
    bool ok;//该面是否是最终凸包中的面
};
struct T3dhull
{
    int n;//初始点数
    TPoint ply[NBmax];//初始点
    int trianglecnt;//凸包上三角形数
    fac tri[10000000];//凸包三角形
    int vis[10000][10000];//点i到点j是属于哪个面
    double dist(TPoint a){return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);}//两点长度
    double area(TPoint a,TPoint b,TPoint c){return dist((b-a)*(c-a));}//三角形面积*2
    double volume(TPoint a,TPoint b,TPoint c,TPoint d){return (b-a)*(c-a)^(d-a);}//四面体有向体积*6
    double ptoplane(TPoint &p,fac &f)//正：点在面同向
    {
        TPoint m=ply[f.b]-ply[f.a],n=ply[f.c]-ply[f.a],t=p-ply[f.a];
        return (m*n)^t; // this is not a distance, because m*n is not normalized
    }
    void deal(int p,int a,int b)
    {
        int f=vis[a][b];//与当前面(cnt)共边(ab)的那个面
        fac add;
        if(tri[f].ok)
        {
            if((ptoplane(ply[p],tri[f]))>PR) dfs(p,f);//如果p点能看到该面f，则继续深度探索f的3条边，以便更新新的凸包面
            else//否则因为p点只看到cnt面，看不到f面，则p点和a、b点组成一个三角形。
            {
                add.a=b,add.b=a,add.c=p,add.ok=1;
                vis[p][b]=vis[a][p]=vis[b][a]=trianglecnt;
                tri[trianglecnt++]=add;
                if( trianglecnt > 9999999 ) { printf("trianglecnt %d > 9999999; too big\n", trianglecnt); exit(0); }
            }
        }
    }
    void dfs(int p,int cnt)//维护凸包，如果点p在凸包外更新凸包
    {
        tri[cnt].ok=0;//当前面需要删除，因为它在更大的凸包里面

//下面把边反过来(先b,后a)，以便在deal()中判断与当前面(cnt)共边(ab)的那个面。即判断与当头面(cnt)相邻的3个面(它们与当前面的共边是反向的，如下图中(1)的法线朝外(即逆时针)的面130和312,它们共边13，但一个方向是13,另一个方向是31)       

        deal(p,tri[cnt].b,tri[cnt].a);
        deal(p,tri[cnt].c,tri[cnt].b);
        deal(p,tri[cnt].a,tri[cnt].c);
    }
    bool same(int s,int e)//判断两个面是否为同一面
    {
        TPoint a=ply[tri[s].a],b=ply[tri[s].b],c=ply[tri[s].c];
        return fabs(volume(a,b,c,ply[tri[e].a]))<PR
            &&fabs(volume(a,b,c,ply[tri[e].b]))<PR
            &&fabs(volume(a,b,c,ply[tri[e].c]))<PR;
    }
    void construct()//构建凸包
    {
        int i,j;
        trianglecnt=0;
        if(n<4) return ;
        bool tmp=true;
        for(i=1;i<n;i++)//前两点不共点
        {
            if((dist(ply[0]-ply[i]))>PR)
            {
                swap(ply[1],ply[i]); tmp=false; break;
            }
        }
        if(tmp) return;
        tmp=true;
        for(i=2;i<n;i++)//前三点不共线
        {
            if((dist((ply[0]-ply[1])*(ply[1]-ply[i])))>PR)
            {
                swap(ply[2],ply[i]); tmp=false; break;
            }
        }
        if(tmp) return ;
        tmp=true;
        for(i=3;i<n;i++)//前四点不共面
        {
            if(fabs((ply[0]-ply[1])*(ply[1]-ply[2])^(ply[0]-ply[i]))>PR)
            {
                swap(ply[3],ply[i]); tmp=false; break;
            }
        }
        if(tmp) return ;


        // printf("OK1\n");

        fac add;
        for(i=0;i<4;i++)//构建初始四面体(4个点为ply[0],ply[1],ply[2],ply[3])
        {
            add.a=(i+1)%4,add.b=(i+2)%4,add.c=(i+3)%4,add.ok=1;
            if((ptoplane(ply[i],add))>0) swap(add.b,add.c);//保证逆时针，即法向量朝外，这样新点才可看到。
            vis[add.a][add.b]=vis[add.b][add.c]=vis[add.c][add.a]=trianglecnt;//逆向的有向边保存
            tri[trianglecnt++]=add;
        }

        // printf("OK2\n");

        for(i=4;i<n;i++)//构建更新凸包
        {
            for(j=0;j<trianglecnt;j++)//对每个点判断是否在当前3维凸包内或外(i表示当前点,j表示当前面)
            {
                 // printf("OK3 i %5d n %5d j %5d trianglecnt %5d\n", i, n, j, trianglecnt);
                if(tri[j].ok&&(ptoplane(ply[i],tri[j]))>PR)//对当前凸包面进行判断，看是否点能否看到这个面
                {
                    dfs(i,j); break;//点能看到当前面，更新凸包的面(递归，可能不止更新一个面)。当前点更新完成后break跳出循环

                }
            }
        }

        // printf("OK4\n");
        int cnt=trianglecnt;//这些面中有一些tri[i].ok=0，它们属于开始建立但后来因为在更大凸包内故需删除的，所以下面几行代码的作用是只保存最外层的凸包
        trianglecnt=0;
        for(i=0;i<cnt;i++)
        {
            if(tri[i].ok) {
                tri[trianglecnt++]=tri[i];
               // printf(" %5d  %5d   %5d  %5d %5d %8.3f %8.3f %8.3f\n",i,trianglecnt, tri[i].a, tri[i].b, tri[i].c, ply[tri[i].a].x, ply[tri[i].a].y,ply[tri[i].a].z);
            }
        }
    }
    double area()//表面积
    {
        double ret=0;
        for(int i=0;i<trianglecnt;i++)
            ret+=area(ply[tri[i].a],ply[tri[i].b],ply[tri[i].c]);
        return ret/2;
    }
    double volume()//体积
    {
        TPoint p(0,0,0);
        double ret=0;
        for(int i=0;i<trianglecnt;i++)
            ret+=volume(p,ply[tri[i].a],ply[tri[i].b],ply[tri[i].c]);
        return fabs(ret/6);
    }
    int facetri() {return trianglecnt;}//表面三角形数
    int facepolygon()//表面多边形数
    {
        int ans=0,i,j,k;
        for(i=0;i<trianglecnt;i++)
        {
            for(j=0,k=1;j<i;j++)
            {
                if(same(i,j)) {k=0;break;}
            }
            ans+=k;
        }
        return ans;
    }

}hull;



extern "C" int hull_ends(double x[NBmax][dd],  int L, double point1[dd], double point2[dd] ) {
   FILE   *fp_tmp;
   int    d, i, j, imin1, imin2, Nplanes;
   double xplane[10000][3][dd], planes[10000][4], dis_1end_to_plane[10000],  mindis1, dis_2end_to_plane[10000],  mindis2;
   double vext1[dd], vext2[dd], xcen[dd];
   double dis_n2n, dis_tmpA, dis_tmpB, frescale1[dd], frescale2[dd], frescale_max1, frescale_max2, xdev_max1[dd], xdev_max2[dd], xdev_tmp;
   int    flag_close_type = 1;

   // printf("OK4\n");   
  
  // L      = 100;
  // fp_tmp = fopen("pos.txt","r");
   hull.n =  L;   

   if(L>10000) {printf("L %10d; > 10000; so it is hard to declare vis[L][L]\n", L );  exit(0); }


  // for(i=0;i<L;i++) fscanf(fp_tmp,"%lf%lf%lf",&x[i][0], &x[i][1], &x[i][2] );
   for(d=0;d<dd;d++) xcen[d] = 0;
   for(i=0;i<L;i++) for(d=0;d<dd;d++)  xcen[d] += x[i][d];
   for(d=0;d<dd;d++) xcen[d] /= L;
   
   for(i=0;i<L;i++) { hull.ply[i].x = x[i][0]; hull.ply[i].y = x[i][1];  hull.ply[i].z = x[i][2];  }

   // printf("OK5\n"); 

   hull.construct();
  
   // printf("OK6\n");

   Nplanes = hull.trianglecnt;

   // printf("Nplanes is %d\n",Nplanes);




   for(i=0;i<Nplanes;i++) {
      xplane[i][0][0] = hull.ply[hull.tri[i].a].x;    xplane[i][0][1] = hull.ply[hull.tri[i].a].y;     xplane[i][0][2] = hull.ply[hull.tri[i].a].z; // #1 point of a plane
      xplane[i][1][0] = hull.ply[hull.tri[i].b].x;    xplane[i][1][1] = hull.ply[hull.tri[i].b].y;     xplane[i][1][2] = hull.ply[hull.tri[i].b].z; // #2 point of a plane
      xplane[i][2][0] = hull.ply[hull.tri[i].c].x;    xplane[i][2][1] = hull.ply[hull.tri[i].c].y;     xplane[i][2][2] = hull.ply[hull.tri[i].c].z; // #3 point of a plane

      cal_normals( xplane[i][0],  xplane[i][1], xplane[i][2], planes[i]);
      


      dis_1end_to_plane[i]  =  fabs(  x[0][0]*planes[i][0] +   x[0][1]*planes[i][1] +   x[0][2]*planes[i][2] + planes[i][3]); 
      dis_2end_to_plane[i]  =  fabs(x[L-1][0]*planes[i][0] + x[L-1][1]*planes[i][1] + x[L-1][2]*planes[i][2] + planes[i][3]);


      // printf("i %5d %10.3f %10.3f\n",i, dis_1end_to_plane[i], dis_2end_to_plane[i]);
   }


    dis_n2n = 0;
    for(d=0;d<dd;d++)  dis_n2n +=  (x[L-1][d] - x[0][d]) * (x[L-1][d] - x[0][d]);
    dis_n2n = sqrt(dis_n2n);

    imin1  = imin_array( dis_1end_to_plane, Nplanes );     mindis1  =  dis_1end_to_plane[imin1];  for(d=0;d<dd;d++) vext1[d] = planes[imin1][d];
    imin2  = imin_array( dis_2end_to_plane, Nplanes );     mindis2  =  dis_2end_to_plane[imin2];  for(d=0;d<dd;d++) vext2[d] = planes[imin2][d];

    // it is shorter to close two ends by a line
    if( dis_n2n < (mindis1 + mindis2) )  { 
       flag_close_type  = 0; 
      // printf("close ends by a line\n") ; 
       return 0; 
    }
   
    if(mindis1 > 1e-6) {
       for(d=0;d<dd;d++) point1[d] = x[0][d] + mindis1 * vext1[d]; 
       dis_tmpA = point1[0]*planes[imin1][0] + point1[1]*planes[imin1][1] + point1[2]*planes[imin1][2] + planes[imin1][3];       

       for(d=0;d<dd;d++) point1[d] = x[0][d] - mindis1 * vext1[d];
       dis_tmpB = point1[0]*planes[imin1][0] + point1[1]*planes[imin1][1] + point1[2]*planes[imin1][2] + planes[imin1][3];

       if(fabs(dis_tmpA) > 1e-4 &&  fabs(dis_tmpB) > 1e-4   )  {printf("disA %10.5f or disB %10.5f should be zero\n",dis_tmpA, dis_tmpB); return 0; /*exit(0);*/ }

       if( fabs(dis_tmpB) < fabs(dis_tmpA)  )  { for(d=0;d<dd;d++) vext1[d] = -vext1[d]; }
       else {  for(d=0;d<dd;d++)  point1[d] = x[0][d] + mindis1 * vext1[d];  }
    }

    if(mindis2 > 1e-6) {
       for(d=0;d<dd;d++) point2[d] = x[L-1][d] + mindis2 * vext2[d];                 
       dis_tmpA = point2[0]*planes[imin2][0] + point2[1]*planes[imin2][1] + point2[2]*planes[imin2][2] + planes[imin2][3];                     

       for(d=0;d<dd;d++) point2[d] = x[L-1][d] - mindis2 * vext2[d];
       dis_tmpB = point2[0]*planes[imin2][0] + point2[1]*planes[imin2][1] + point2[2]*planes[imin2][2] + planes[imin2][3];

      if(fabs(dis_tmpA) > 1e-4 &&  fabs(dis_tmpB) > 1e-4   )  {printf("disA %10.5f or disB %10.5f should be zero\n",dis_tmpA, dis_tmpB); return 0; /*exit(0);*/ }

       if( fabs(dis_tmpB) < fabs(dis_tmpA)  )  { for(d=0;d<dd;d++) vext2[d] = -vext2[d]; }
       else {for(d=0;d<dd;d++)  point2[d] = x[L-1][d] + mindis2 * vext2[d];  }
    }


    if(mindis1 < 1e-6)  { for(d=0;d<dd;d++) vext1[d] = x[0][d]   - xcen[d]; }
    if(mindis2 < 1e-6)  { for(d=0;d<dd;d++) vext2[d] = x[L-1][d] - xcen[d]; }
           

    // measure max distance to 1st and last ends
    for(d=0;d<dd;d++)  { xdev_max1[d] = 10.0;  xdev_max2[d] = 10.0; }
    for(d=0;d<dd;d++)  {
	    for(i=0;i<L;i++) {
               xdev_tmp = fabs(x[i][d] - x[0][d]);
               if( xdev_tmp > xdev_max1[d])   xdev_max1[d] = xdev_tmp;
	    }
    } 

    for(d=0;d<dd;d++)  {
	    for(i=0;i<L;i++) {
                xdev_tmp = fabs(x[i][d] - x[L-1][d]);
                if( xdev_tmp > xdev_max2[d])   xdev_max2[d] = xdev_tmp;
	    }
    }


    for(d=0;d<dd;d++) frescale1[d] = 2*xdev_max1[d]/(fabs(vext1[d]));
    for(d=0;d<dd;d++) frescale2[d] = 2*xdev_max2[d]/(fabs(vext2[d]));

    frescale_max1 = frescale1[0]; 
    if( frescale1[1] > frescale_max1 ) frescale_max1 = frescale1[1];  
    if( frescale1[2] > frescale_max1 ) frescale_max1 = frescale1[2];

    frescale_max2 = frescale2[0];
    if( frescale2[1] > frescale_max2 ) frescale_max2 = frescale2[1];
    if( frescale2[2] > frescale_max2 ) frescale_max2 = frescale2[2];

    for(d=0;d<dd;d++)  point1[d]  = x[0][d]   + frescale_max1 * vext1[d];
    for(d=0;d<dd;d++)  point2[d]  = x[L-1][d] + frescale_max2 * vext2[d];

  
    // for(i=0;i<L;i++) printf("%10.3f %10.3f %10.3f\n", x[i][0], x[i][1], x[i][2] );
    // printf("%10.3f %10.3f %10.3f\n", point1[0], point1[1],point1[2] );
    // printf("%10.3f %10.3f %10.3f\n", point2[0], point2[1],point2[2] );

    flag_close_type = 1;



    return flag_close_type;
}


