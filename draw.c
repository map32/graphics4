#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"

/*======== void add_sphere() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double r
	    double step  
  Returns: 

  adds all the points for a sphere with center 
  (cx, cy) and radius r.

  should call generate_sphere to create the
  necessary points

  jdyrlandweaver
  ====================*/
void add_sphere( struct matrix * points, 
		 double cx, double cy, double r, 
		 double step ) {
  struct matrix * pts;
  pts = new_matrix(4,100);
  generate_sphere(pts,cx,cy,r,step);
  int i=0;
  for(i=0;i<pts->lastcol-1;i+=2){
    add_edge(points,pts->m[0][i],pts->m[1][i],pts->m[2][i],
	     pts->m[0][i+1],pts->m[1][i+1],pts->m[2][i+1]);
  }
  free_matrix(pts);
}

/*======== void generate_sphere() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double r
	    double step  
  Returns: 

  Generates all the points along the surface of a 
  sphere with center (cx, cy) and radius r

  Adds these points to the matrix parameter

  03/22/12 11:30:26
  jdyrlandweaver
  ====================*/
void generate_sphere( struct matrix * points, 
		      double cx, double cy, double r, 
		      double step ) {
  double a=0;
  double b=0;
  double x,y,z;
  x = r*cos(M_PI*2*b)+cx;
  y = r*sin(M_PI*2*b)*cos(M_PI*a)+cy;
  z = r*sin(M_PI*2*b)*sin(M_PI*a);
  for(a=step;a<1+step;a+=step){
    for(b=step;b<1+step;b+=step){
      add_point(points,x,y,z);
      x = r*cos(M_PI*2*b)+cx;
      y = r*sin(M_PI*2*b)*cos(M_PI*a)+cy;
      z = r*sin(M_PI*2*b)*sin(M_PI*a);
      add_point(points,x,y,z);
    }
  }
  a=0;
  b=0;
  x = r*cos(M_PI*2*b)+cx;
  y = r*sin(M_PI*2*b)*cos(M_PI*a)+cy;
  z = r*sin(M_PI*2*b)*sin(M_PI*a);
  for(b=0;b<1+step;b+=step){
    for(a=0;a<1+step;a+=step){
      add_point(points,x,y,z);
      x = r*cos(M_PI*2*b)+cx;
      y = r*sin(M_PI*2*b)*cos(M_PI*a)+cy;
      z = r*sin(M_PI*2*b)*sin(M_PI*a);
      add_point(points,x,y,z);
    }
  }
}    

/*======== void add_torus() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double r1
	    double r2
	    double step  
  Returns: 

  adds all the points required to make a torus
  with center (cx, cy) and radii r1 and r2.

  should call generate_torus to create the
  necessary points

  03/22/12 13:34:03
  jdyrlandweaver
  ====================*/
void add_torus( struct matrix * points, 
		double cx, double cy, double r1, double r2, 
		double step ) {
  struct matrix * pts;
  pts = new_matrix(4,100);
  generate_torus(pts,cx,cy,r1,r2,step);
  int i=0;
  for(i=0;i<pts->lastcol-1;i+=2){
    add_edge(points,pts->m[0][i],pts->m[1][i],pts->m[2][i],
	     pts->m[0][i+1],pts->m[1][i+1],pts->m[2][i+1]);
  }
  free_matrix(pts);
}

/*======== void generate_torus() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double r
	    double step  
  Returns: 

  Generates all the points along the surface of a 
  tarus with center (cx, cy) and radii r1 and r2

  Adds these points to the matrix parameter

  jdyrlandweaver
  ====================*/
void generate_torus( struct matrix * points, 
		     double cx, double cy, double r1, double r2, 
		     double step ) {
  double a=0;
  double b=0;
  double x,y,z;
  x = r1*cos(M_PI*2*b)+cx;
  y = cos(M_PI*2*a)*(r1*sin(M_PI*2*b)+r2)+cy;
  z = sin(M_PI*2*a)*(r1*sin(M_PI*2*b)+r2);
  for(a=step;a<1;a+=step){
    for(b=step;b<1;b+=step){
      add_point(points,x,y,z);
      x = r1*cos(M_PI*2*b)+cx;
      y = cos(M_PI*2*a)*(r1*sin(M_PI*2*b)+r2)+cy;
      z = sin(M_PI*2*a)*(r1*sin(M_PI*2*b)+r2);
      add_point(points,x,y,z);
    }
  }
  a=0;
  b=0;
  x = r1*cos(M_PI*2*b)+cx;
  y = cos(M_PI*2*a)*(r1*sin(M_PI*2*b)+r2)+cy;
  z = sin(M_PI*2*a)*(r1*sin(M_PI*2*b)+r2);
  for(b=step;b<1;b+=step){
    for(a=step;a<1;a+=step){
      add_point(points,x,y,z);
      x = r1*cos(M_PI*2*b)+cx;
      y = cos(M_PI*2*a)*(r1*sin(M_PI*2*b)+r2)+cy;
      z = sin(M_PI*2*a)*(r1*sin(M_PI*2*b)+r2);
      add_point(points,x,y,z);
    }
  }
}

/*======== void add_box() ==========
  Inputs:   struct matrix * points
            double x
	    double y
	    double z
	    double width
	    double height
	    double depth
  Returns: 

  add the points for a rectagular prism whose 
  upper-left corner is (x, y, z) with width, 
  height and depth dimensions.

  jdyrlandweaver
  ====================*/
void add_box( struct matrix * points,
	      double x, double y, double z,
	      double width, double height, double depth ) {
  add_edge(points,x,y,z,x+width,y,z);
  add_edge(points,x,y,z,x,y-height,z);
  add_edge(points,x,y,z,x,y,z-depth);
  add_edge(points,x+width,y-height,z,x,y-height,z);
  add_edge(points,x+width,y-height,z,x+width,y,z);
  add_edge(points,x+width,y-height,z,x+width,y-height,z-depth);
  add_edge(points,x+width,y,z-depth,x,y,z-depth);
  add_edge(points,x+width,y,z-depth,x+width,y-height,z-depth);
  add_edge(points,x+width,y,z-depth,x+width,y,z);
  add_edge(points,x,y-height,z-depth,x+width,y-height,z-depth);
  add_edge(points,x,y-height,z-depth,x,y,z-depth);
  add_edge(points,x,y-height,z-depth,x,y-height,z);
}

/*======== void add_circle() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double y
	    double step  
  Returns: 


  03/16/12 19:53:52
  jdyrlandweaver
  ====================*/
void add_circle( struct matrix * points, 
		 double cx, double cy, 
		 double r, double step ) {
  double x , y, t;
  double i;
  int x0 = cx;
  int y0 = r+cy;
  for(i=step;i<=1.;i+=step){
    t = 2*M_PI*i;
    x = r*sin(t)+cx;
    y = r*cos(t)+cy;
    add_edge(points,x0,y0,0,x,y,0);
    x0 = x;
    y0 = y;
  }
  add_edge(points,x0,y0,0,cx,r+cy,0);
}

/*======== void add_curve() ==========
Inputs:   struct matrix *points
         double x0
         double y0
         double x1
         double y1
         double x2
         double y2
         double x3
         double y3
         double step
         int type  
Returns: 

Adds the curve bounded by the 4 points passsed as parameters
of type specified in type (see matrix.h for curve type constants)
to the matrix points

03/16/12 15:24:25
jdyrlandweaver
====================*/
void add_curve( struct matrix *points, 
		double x0, double y0, 
		double x1, double y1, 
		double x2, double y2, 
		double x3, double y3, 
		double step, int type ) {
  struct matrix *xcfs, *ycfs;
  double x,y,xs,ys;
  xs = x0;
  ys = y0;
  if(type==0){ //hermy
    double s1 = (y1-y0)/(x1-x0);
    double s2 = (y2-y3)/(x2-x3);
    if(x1==x0){
      if(y1>y0){
	s1 = y1-y0;
      } else {
	s1 = y0-y1;
      }
    }
    if(x3==x2){
      if(y3>y2){
	s2 = y2-y3;
      } else {
	s2 = y3-y2;
      }
    }
    //printf("%lf %lf %lf %lf %lf %lf %lf %lf all of thse\n",y1,y0,x1,x0,y3,y2,x3,x2);
    xcfs = generate_curve_coefs(x0,x2,s1,s2,0);
    ycfs = generate_curve_coefs(y0,y2,s1,s2,0);
  } else { //bezzy
    xcfs = generate_curve_coefs(x0,x1,x2,x3,1);
    ycfs = generate_curve_coefs(y0,y1,y2,y3,1);
  }
  double i;
  for(i=step;i<1;i+=step){
      x = ((xcfs->m[0][0]*i+xcfs->m[1][0])*i+xcfs->m[2][0])*i+xcfs->m[3][0];
      y = ((ycfs->m[0][0]*i+ycfs->m[1][0])*i+ycfs->m[2][0])*i+ycfs->m[3][0];
      //print_matrix(xcfs);
      //print_matrix(ycfs);
      add_edge(points,xs,ys,0,x,y,0);
      xs = x;
      ys = y;
    }
  //add_edge(points,xs,ys,0,x3,y3,0);
    //printf("\n");
}

/*======== void add_point() ==========
Inputs:   struct matrix * points
         int x
         int y
         int z 
Returns: 
adds point (x, y, z) to points and increment points.lastcol
if points is full, should call grow on points
====================*/
void add_point( struct matrix * points, double x, double y, double z) {
  
  if ( points->lastcol == points->cols )
    grow_matrix( points, points->lastcol + 100 );

  points->m[0][points->lastcol] = x;
  points->m[1][points->lastcol] = y;
  points->m[2][points->lastcol] = z;
  points->m[3][points->lastcol] = 1;

  points->lastcol++;
}

/*======== void add_edge() ==========
Inputs:   struct matrix * points
          int x0, int y0, int z0, int x1, int y1, int z1
Returns: 
add the line connecting (x0, y0, z0) to (x1, y1, z1) to points
should use add_point
====================*/
void add_edge( struct matrix * points, 
	       double x0, double y0, double z0, 
	       double x1, double y1, double z1) {
  add_point( points, x0, y0, z0 );
  add_point( points, x1, y1, z1 );
}

/*======== void draw_lines() ==========
Inputs:   struct matrix * points
         screen s
         color c 
Returns: 
Go through points 2 at a time and call draw_line to add that line
to the screen
====================*/
void draw_lines( struct matrix * points, screen s, color c) {

  int i;
 
  if ( points->lastcol < 2 ) {
    
    printf("Need at least 2 points to draw a line!\n");
    return;
  }

  for ( i = 0; i < points->lastcol - 1; i+=2 ) {
    draw_line( points->m[0][i], points->m[1][i], 
	       points->m[0][i+1], points->m[1][i+1], s, c);
  } 	       
}

void draw_line(int x0, int y0, int x1, int y1, screen s, color c) {
 
  int x, y, d, dx, dy;

  x = x0;
  y = y0;
  
  //swap points so we're always draing left to right
  if ( x0 > x1 ) {
    x = x1;
    y = y1;
    x1 = x0;
    y1 = y0;
  }

  //need to know dx and dy for this version
  dx = (x1 - x) * 2;
  dy = (y1 - y) * 2;

  //positive slope: Octants 1, 2 (5 and 6)
  if ( dy > 0 ) {

    //slope < 1: Octant 1 (5)
    if ( dx > dy ) {
      d = dy - ( dx / 2 );
  
      while ( x <= x1 ) {
	plot(s, c, x, y);

	if ( d < 0 ) {
	  x = x + 1;
	  d = d + dy;
	}
	else {
	  x = x + 1;
	  y = y + 1;
	  d = d + dy - dx;
	}
      }
    }

    //slope > 1: Octant 2 (6)
    else {
      d = ( dy / 2 ) - dx;
      while ( y <= y1 ) {

	plot(s, c, x, y );
	if ( d > 0 ) {
	  y = y + 1;
	  d = d - dx;
	}
	else {
	  y = y + 1;
	  x = x + 1;
	  d = d + dy - dx;
	}
      }
    }
  }

  //negative slope: Octants 7, 8 (3 and 4)
  else { 

    //slope > -1: Octant 8 (4)
    if ( dx > abs(dy) ) {

      d = dy + ( dx / 2 );
  
      while ( x <= x1 ) {

	plot(s, c, x, y);

	if ( d > 0 ) {
	  x = x + 1;
	  d = d + dy;
	}
	else {
	  x = x + 1;
	  y = y - 1;
	  d = d + dy + dx;
	}
      }
    }

    //slope < -1: Octant 7 (3)
    else {

      d =  (dy / 2) + dx;

      while ( y >= y1 ) {
	
	plot(s, c, x, y );
	if ( d < 0 ) {
	  y = y - 1;
	  d = d + dx;
	}
	else {
	  y = y - 1;
	  x = x + 1;
	  d = d + dy + dx;
	}
      }
    }
  }
}

