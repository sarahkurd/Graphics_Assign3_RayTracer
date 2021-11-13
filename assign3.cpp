/*
CSCI 420
Assignment 3 Raytracer

Name: Sarah Kurdoghlian
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "pic.h"

// For Linux
// #include <GL/gl.h>
// #include <GL/glu.h>
// #include <GL/glut.h>
// For Mac
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0
#define PI 3.14159265
#define FOV_RADIANS (fov*PI)/180

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

float aspect_ratio = 640/480;

/**
 *  1. Generate the ray going through pixel (x, y)
 *  2. Calculate the values for at^2 + bt + c = 0 
 *  3. Return the intersection point
 **/
float ray_sphere_intersection(unsigned int x, unsigned int y, Sphere sphere) {
  float ray_origin[3] = {0, 0, 0};

  float x_dir = -aspect_ratio * tan(FOV_RADIANS/2) + x;
  float y_dir = tan(FOV_RADIANS/2) - y;
  float z_dir = -1;

  float a = (x_dir * x_dir) + (y_dir * y_dir) + (z_dir * z_dir);
  float b = 2 * ((x_dir * (ray_origin[0] - sphere.position[0])) + (y_dir * (ray_origin[1] - sphere.position[1])) + (z_dir * (ray_origin[2] - sphere.position[2])));
  float c= ((ray_origin[0] - sphere.position[0]) * (ray_origin[0] - sphere.position[0])) +
           ((ray_origin[1] - sphere.position[1]) * (ray_origin[1] - sphere.position[1])) +
           ((ray_origin[2] - sphere.position[2]) * (ray_origin[2] - sphere.position[2])) -
           (sphere.radius * sphere.radius);

  float discriminant = (b * b) - (4 * a * c);
  if (discriminant < 0 ) {
      return -1.0;
  } else {
    return (-b - sqrt(discriminant)) / (2 * a);
  }
}

/**
 * Iterate over objects and calculate their intersection point with the ray
 * Return the closest intersection point
 * 
 * Params: The ray going through pixel (x, y)
 **/
void calc_closest_intersection(unsigned int x, unsigned int y) {
  double closest_intersection[2];
  if (num_spheres != 0) {
    Sphere current_sphere;
    for(int i = 0; i < num_spheres; i++) {
      current_sphere = spheres[i];
    }
  }

  if (num_triangles != 0) {
    Triangle current_triangle;
    for(int i = 0; i < num_triangles; i++) {
      current_triangle = triangles[i];
    }
  }
}

//MODIFY THIS FUNCTION
void draw_scene()
{
  unsigned int x,y;
  // for each roay
  for(x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++)
    {
      // for each object calculate closest interection
      plot_pixel(x,y,x%256,y%256,(x+y)%256);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); 
  fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

void parse_check(char *expected,char *found)
{
  if(strcasecmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);

      if(strcasecmp(type,"triangle")==0) {

        printf("found triangle\n");
        int j;

        for(j=0;j < 3;j++)
          {
            parse_doubles(file,"pos:",t.v[j].position);
            parse_doubles(file,"nor:",t.v[j].normal);
            parse_doubles(file,"dif:",t.v[j].color_diffuse);
            parse_doubles(file,"spe:",t.v[j].color_specular);
            parse_shi(file,&t.v[j].shininess);
          }

        if(num_triangles == MAX_TRIANGLES)
          {
            printf("too many triangles, you should increase MAX_TRIANGLES!\n");
            exit(0);
          }
	        triangles[num_triangles++] = t;
	}
      else if(strcasecmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(strcasecmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void reshape(int w, int h) 
{
  glViewport(0, 0, w, h);
  aspect_ratio = w/h;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  // setup projection
  gluPerspective(fov, w/h, 0.01, 100);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void display()
{
  draw_scene();
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
