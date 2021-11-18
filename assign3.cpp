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

Sphere current_sphere;

float aspect_ratio = WIDTH/HEIGHT;
float ray_origin[3] = {0, 0, 0};
float x_dir = 0;
float y_dir = 0;
float z_dir = 0;

float current_intersection_point[3] = {0, 0, 0};
float closest_intersection[3] = {10000, 10000, 10000};
bool does_intersect = false;

float normal_vector[3] = {0, 0, 0};

float current_color_at_point[3] = {0, 0, 0};

float min_val(float a, float b) {
  return a < b ? a : b;
}

float magnitide(float x, float y, float z) {
  return sqrt((x * x) + (y * y) + (z * z));
}

/**
 *  1. Generate the ray going through pixel (x, y)
 *  2. Calculate the values for at^2 + bt + c = 0 
 *  3. Return the intersection point
 **/
float ray_sphere_intersection(int x, int y, Sphere sphere) {
  printf("\npixel(%i, %i)\n", x, y);

  // calculate the length of each pixel using the WIDTH and HEIGHT of the window
  float x_pixel_length = (aspect_ratio * tan(FOV_RADIANS/2) - (-aspect_ratio * tan(FOV_RADIANS/2))) / WIDTH;
  float y_pixel_length = (tan(FOV_RADIANS/2) - (-tan(FOV_RADIANS/2))) / HEIGHT;

  // calculate direction vector 
  x_dir = (-aspect_ratio * tan(FOV_RADIANS/2)) + ((x + 0.5) * x_pixel_length);
  y_dir = tan(FOV_RADIANS/2) - ((y + 0.5) * y_pixel_length);
  z_dir = -1;

  // normalize the direction
  float mag = magnitide(x_dir, y_dir, z_dir);
  x_dir = x_dir/mag;
  y_dir = y_dir/mag;
  z_dir = z_dir/mag;

  printf("xdir: %f    , ydir: %f     , zdir: %f    \n", x_dir, y_dir, z_dir);

  float a = (x_dir * x_dir) + (y_dir * y_dir) + (z_dir * z_dir);
  float b = 2 * ((x_dir * (ray_origin[0] - sphere.position[0])) + (y_dir * (ray_origin[1] - sphere.position[1])) + (z_dir * (ray_origin[2] - sphere.position[2])));
  float c= ((ray_origin[0] - sphere.position[0]) * (ray_origin[0] - sphere.position[0])) +
           ((ray_origin[1] - sphere.position[1]) * (ray_origin[1] - sphere.position[1])) +
           ((ray_origin[2] - sphere.position[2]) * (ray_origin[2] - sphere.position[2])) -
           (sphere.radius * sphere.radius);
  printf("a: %f    , b: %f     , c: %f  \n", a, b, c);

  float discriminant = (b * b) - (4 * a * c);
  printf("discriminant: %f\n", discriminant);
  // negative discriminant means ray missed sphere, abort rest of calculations
  if (discriminant < 0) {
      return -1.0;
  } else {
    does_intersect = true;
    float t0 = (-b - sqrt(discriminant)) / (2 * a); // nearest intersection point
    //float t1 = (-b + sqrt(discriminant)) / (2 * a);
    printf("t0: %f\n", t0);
    return t0;
  }
}

float get_intersect_coordinates(float min_root) {
  current_intersection_point[0] = ray_origin[0] + (min_root * x_dir);
  current_intersection_point[1] = ray_origin[1] + (min_root * y_dir);
  current_intersection_point[2] = ray_origin[2] + (min_root * z_dir);
  printf("current intersection point: (%f, %f, %f)\n", current_intersection_point[0], current_intersection_point[1], current_intersection_point[2]);
}

float distance_to_cop(float x, float y, float z) {
  return sqrt((x * x) + (y * y) + (z * z));
}

/**
 * Iterate over objects and calculate their intersection point with the ray
 * Return the closest intersection point
 * 
 * Params: The ray going through pixel (x, y)
 * Return: The closest intersection point (x, y z) out of all the objects
 **/
void iterate_over_objects(int x, int y) {
  float closest_distance = distance_to_cop(closest_intersection[0], closest_intersection[1], closest_intersection[2]);
  if (num_spheres != 0) {
    for(int i = 0; i < num_spheres; i++) {
      current_sphere = spheres[i];
      float min_root = ray_sphere_intersection(x, y, current_sphere);

      // if ray intersects the sphere at 1 or 2 locations
      if (min_root != -1) {
        get_intersect_coordinates(min_root);

        // continue updating closest intersection point
        float dis = distance_to_cop(current_intersection_point[0], current_intersection_point[1], current_intersection_point[2]);
        printf("distance to current point: %f , vs. closest distance: %f\n", dis, closest_distance);
        if (dis < closest_distance) {
          closest_distance = dis;
          closest_intersection[0] = current_intersection_point[0];
          closest_intersection[1] = current_intersection_point[1];
          closest_intersection[2] = current_intersection_point[2];
        }
      }
    }
  }

  if (num_triangles != 0) {
    Triangle current_triangle;
    for(int i = 0; i < num_triangles; i++) {
      current_triangle = triangles[i];
    }
  }

  printf("closest intersection point: (%f, %f, %f)\n", closest_intersection[0], closest_intersection[1], closest_intersection[2]);
}

/**
 * Use the closest intersection point and calculate sphere surface normal
 * */
void sphere_normal() {
  normal_vector[0] = (closest_intersection[0] - current_sphere.position[0])/current_sphere.radius;
  normal_vector[1] = (closest_intersection[1] - current_sphere.position[1])/current_sphere.radius;
  normal_vector[2] = (closest_intersection[2] - current_sphere.position[2])/current_sphere.radius;
  printf("Sphere normal vector at point (%f, %f, %f): (%f, %f, %f)\n", closest_intersection[0], closest_intersection[1], closest_intersection[2], normal_vector[0], normal_vector[1], normal_vector[2]);
  printf("Sphere normal vector magnitude: %f\n", magnitide(normal_vector[0], normal_vector[1], normal_vector[2]));
}

/**
 * Compute if shadow ray hits the given sphere
 * Return a bool true if ray is in shadow, else return false
 **/
bool shadow_ray_sphere_intersection(Sphere sphere, float shadow_ray[]) {
  // calculate direction vector 
  float dir_x = shadow_ray[0];
  float dir_y = shadow_ray[1];
  float dir_z = shadow_ray[2];

  // origin point of the shadow ray is the current intersection point on the sphere
  float origin_x = current_intersection_point[0];
  float origin_y = current_intersection_point[1];
  float origin_z = current_intersection_point[2];

  // normalize the direction
  float mag = magnitide(dir_x, dir_y, dir_z);
  dir_x = dir_x/mag;
  dir_y = dir_y/mag;
  dir_z = dir_z/mag;

  printf("SHADOW RAY dir_x: %f    , dir_y: %f     , dir_z: %f    \n", dir_x, dir_y, dir_z);

  float a = (dir_x * dir_x) + (dir_y * dir_y) + (dir_z * dir_z);
  float b = 2 * ((dir_x * (origin_x - sphere.position[0])) + (dir_y * (origin_y - sphere.position[1])) + (dir_z * (origin_z - sphere.position[2])));
  float c= ((origin_x - sphere.position[0]) * (origin_x - sphere.position[0])) +
           ((origin_y - sphere.position[1]) * (origin_y - sphere.position[1])) +
           ((origin_z - origin_z) * (origin_z - origin_z)) -
           (sphere.radius * sphere.radius);
  printf("a: %f    , b: %f     , c: %f  \n", a, b, c);

  float discriminant = (b * b) - (4 * a * c);
  printf("discriminant: %f\n", discriminant);

  // negative discriminant means ray missed sphere, abort rest of calculations
  return discriminant < 0 ? true : false;
}

/**
 * Given a shadow ray, do ray-sphere and ray-plane intersection
 * for all the objects in the scene
 * 
 * return: a bool reprenting if the shadow ray is in shadow
 **/
bool shadow_ray_intersection(float shadow_ray[]) {
  bool isInShadow = false;

  if (num_spheres != 0) {
    for(int i = 0; i < num_spheres; i++) {
      Sphere s = spheres[i];
      isInShadow = shadow_ray_sphere_intersection(s, shadow_ray);
      if (isInShadow) {
        return true;
      }
    }
  }

  if (num_triangles != 0) {
      // TODO()
  }

  // if shadow ray never intersected with an object in the scene
  return false;
}


/**
 * Shadow_ray: unit vector to light
 * normal_vector: surface normal
 * q: the distance to light source (shadow ray's magnitude)
 **/
void phong_illumination(float shadow_ray[], float q_distance_to_light_source, Light current_light) {
  float diffuse;
  float specular;
  float k_diffuse = 0.5;
  float k_specular = 0.75;
  float shininess = 2.0;

  // unit vector to camera
  float v_to_camera[3] = {0, 0, 0};
  v_to_camera[0] = 0 - current_intersection_point[0];
  v_to_camera[1] = 0 - current_intersection_point[1];
  v_to_camera[2] = 0 - current_intersection_point[2];

  float mag = magnitide(v_to_camera[0], v_to_camera[1], v_to_camera[2]);
  v_to_camera[0] = v_to_camera[0]/mag;
  v_to_camera[1] = v_to_camera[1]/mag;
  v_to_camera[2] = v_to_camera[2]/mag;

  float reflected_vector[3] = {0, 0, 0};

  float light_dot_normal = (shadow_ray[0] * normal_vector[0]) + (shadow_ray[1] * normal_vector[1]) + (shadow_ray[2] * normal_vector[2]);
  // clamp value to zero if necessary
  if (light_dot_normal < 0) {
    light_dot_normal = 0;
  }
  float illumination_diffuse = k_diffuse * light_dot_normal;

  reflected_vector[0] = (2 * light_dot_normal + normal_vector[0]) - shadow_ray[0];
  reflected_vector[1] = (2 * light_dot_normal + normal_vector[1]) - shadow_ray[1];
  reflected_vector[2] = (2 * light_dot_normal + normal_vector[2]) - shadow_ray[2];
  printf("REFLECTED VECTOR: (%f, %f, %f)\n", reflected_vector[0], reflected_vector[1], reflected_vector[2]);
  printf("REFLECTED VECTOR MAGNITUDE: %f\n", magnitide(reflected_vector[0], reflected_vector[1], reflected_vector[2]));

  float reflected_dot_viewer = (v_to_camera[0] * reflected_vector[0]) + (v_to_camera[1] * reflected_vector[1]) + (v_to_camera[2] * reflected_vector[2]);
  // clamp value to zero if necessary
  if (reflected_dot_viewer < 0) {
    reflected_dot_viewer = 0;
  }
  float illumination_specular = k_specular * pow(reflected_dot_viewer, shininess);

  // Add to the RGB color channels separately
  current_color_at_point[0] += current_light.color[0] * (illumination_diffuse + illumination_specular);
  current_color_at_point[1] += current_light.color[1] * (illumination_diffuse + illumination_specular);
  current_color_at_point[2] += current_light.color[2] * (illumination_diffuse + illumination_specular);
}

/**
 * For each light source, fire a shadow ray
 * Determine if light hits the surface
 **/
void fire_shadow_rays() {
  for (int i = 0; i < num_lights; i++) {
    Light current_light = lights[i];
    float shadow_ray[3] = {0, 0, 0};
    shadow_ray[0] = current_light.position[0] - closest_intersection[0];
    shadow_ray[1] = current_light.position[1] - closest_intersection[1];
    shadow_ray[2] = current_light.position[2] - closest_intersection[2];

    bool isInShadow = shadow_ray_intersection(shadow_ray);
    if (isInShadow) {
      current_color_at_point[0] += 0;
      current_color_at_point[1] += 0;
      current_color_at_point[2] += 0;
    } else {
      // 6) for each unblocked shadow ray, evaluate local Phong model for that light
      
      // normalize the shadow ray
      float mag = magnitide(shadow_ray[0], shadow_ray[1], shadow_ray[2]);
      shadow_ray[0] = shadow_ray[0]/mag;
      shadow_ray[1] = shadow_ray[1]/mag;
      shadow_ray[2] = shadow_ray[2]/mag;

      phong_illumination(shadow_ray, mag, current_light);
    }
    // The final color of the point is the sum of the contributions from all lights, plus the ambient color
    current_color_at_point[0] += ambient_light[0];
    current_color_at_point[1] += ambient_light[1];
    current_color_at_point[2] += ambient_light[2];

    // clamp color to 1.0 if necessary
    if (current_color_at_point[0] > 1) {
      current_color_at_point[0] = 1;
    } 
    if (current_color_at_point[1] > 1) {
      current_color_at_point[1] = 1;
    } 
    if (current_color_at_point[2] > 1) {
      current_color_at_point[2] = 1;
    } 
  }
}

/**
 *  Ray Casting Algorithm broken down step by step
 **/ 
void draw_scene()
{
  unsigned int x,y;
  // 1) for each ray
  for(x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++)
    {
      // 2) for each object, calculate closest interection 
      // 3) and set the global closest intersection point
      iterate_over_objects(x, y);
      if (does_intersect) {
        // 4) for closest intersection point p, calculate surface normal
        sphere_normal();

        // 5) for each light source, fire a shadow ray
        fire_shadow_rays();

        // 6) use global color value to color the pixel
        plot_pixel(x, y, current_color_at_point[0], current_color_at_point[1], current_color_at_point[2]);

        does_intersect = false;
      }
      // reset the global color variable to all zeros
      // reset the global closest intersection value for the next pixel we will be evaluating
      current_color_at_point[0] = 0;
      current_color_at_point[1] = 0;
      current_color_at_point[2] = 0;

      closest_intersection[0] = 10000;
      closest_intersection[1] = 10000;
      closest_intersection[2] = 10000;
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

  glClearColor(255, 255, 255,0);
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
