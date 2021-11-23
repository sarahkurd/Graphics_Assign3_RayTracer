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
#define WIDTH 320
#define HEIGHT 240

//the field of view of the camera
#define fov 60.0
#define PI 3.14159265
#define FOV_RADIANS fov / 2 * PI / 180

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
Triangle current_triangle;

float aspect_ratio = WIDTH/HEIGHT;
float ray_origin[3] = {0, 0, 0};
float x_dir = 0;
float y_dir = 0;
float z_dir = 0;

double current_intersection_point[3] = {0, 0, 0};
float closest_intersection[3] = {10000, 10000, 10000};
bool does_intersect = false;

float normal_vector[3] = {0, 0, 0};

float current_color_at_point[3] = {0, 0, 0};

float shadow_ray[3] = {0, 0, 0};

float normal_of_plane[3] = {0, 0, 0};

double triangle_area_cross_product[3] = {0, 0, 0};

float triangle_intersection_normal[3] = {0, 0, 0};
float interpolated_diffuse[3] = {0, 0, 0};
float interpolated_specular[3] = {0, 0, 0};

float current_alpha = 0.0;
float current_beta = 0.0;
float current_gamma = 0.0;

typedef enum { SPHERE, TRIANGLE } OBJECT;
OBJECT object_type;

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
float ray_sphere_intersection(int x, int y) {
  printf("\npixel(%i, %i)\n", x, y);

  // calculate the length of each pixel using the WIDTH and HEIGHT of the window
  float x_pixel_length = (aspect_ratio * tan(FOV_RADIANS) - (-aspect_ratio * tan(FOV_RADIANS))) / WIDTH;
  float y_pixel_length = (tan(FOV_RADIANS) - (-tan(FOV_RADIANS))) / HEIGHT;

  // calculate direction vector 
  x_dir = (-aspect_ratio * tan(FOV_RADIANS)) + ((x + 0.5) * x_pixel_length);
  y_dir = tan(FOV_RADIANS) - ((y - 0.5) * y_pixel_length);
  z_dir = -1;

  // normalize the direction
  float mag = magnitide(x_dir, y_dir, z_dir);
  x_dir = x_dir/mag;
  y_dir = y_dir/mag;
  z_dir = z_dir/mag;

  printf("xdir: %f    , ydir: %f     , zdir: %f    \n", x_dir, y_dir, z_dir);

  float a = (x_dir * x_dir) + (y_dir * y_dir) + (z_dir * z_dir);
  float b = 2 * ((x_dir * (ray_origin[0] - current_sphere.position[0])) + (y_dir * (ray_origin[1] - current_sphere.position[1])) + (z_dir * (ray_origin[2] - current_sphere.position[2])));
  float c= ((ray_origin[0] - current_sphere.position[0]) * (ray_origin[0] - current_sphere.position[0])) +
           ((ray_origin[1] - current_sphere.position[1]) * (ray_origin[1] - current_sphere.position[1])) +
           ((ray_origin[2] - current_sphere.position[2]) * (ray_origin[2] - current_sphere.position[2])) -
           (current_sphere.radius * current_sphere.radius);
  printf("a: %f    , b: %f     , c: %f  \n", a, b, c);

  float discriminant = (b * b) - (4 * a * c);
  printf("discriminant: %f\n", discriminant);
  // negative discriminant means ray missed sphere, abort rest of calculations
  if (discriminant < 0) {
      return -1.0;
  } else {
    float t0 = (-b - sqrt(discriminant)) / (2 * a); // nearest intersection point
    float t1 = (-b + sqrt(discriminant)) / (2 * a); // nearest intersection point
    printf("t0: %f\n", t0);
    if (t0 > 0 && t1 > 0 && discriminant > 0) {
      does_intersect = true;
      return min_val(t0, t1);
    } else {
      return -1;
    }
  }
}

/**
 * Helper function to calculate the normal vector of a plane
 * given 3 vertices on the plane
 **/
void plane_normal(Vertex v1, Vertex v2, Vertex v3) {
  float ray1[3] = {0, 0, 0};
  float ray2[3] = {0, 0, 0};

  // v1 - v2
  ray1[0] = v1.position[0] - v2.position[0];
  ray1[1] = v1.position[1] - v2.position[1];
  ray1[2] = v1.position[2] - v2.position[2];

  // v3 - v2
  ray2[0] = v3.position[0] - v2.position[0];
  ray2[1] = v3.position[1] - v2.position[1];
  ray2[2] = v3.position[2] - v2.position[2];

  // cross product = (v1 - v2) x (v3 - v2)
  normal_of_plane[0] = (ray1[1] * ray2[2]) - (ray1[2] * ray2[1]);
  normal_of_plane[1] = (ray1[2] * ray2[0]) - (ray1[0] * ray2[2]);
  normal_of_plane[2] = (ray1[0] * ray2[1]) - (ray1[1] * ray2[0]);

  float mag = magnitide(normal_of_plane[0], normal_of_plane[1], normal_of_plane[2]);
  normal_of_plane[0] = normal_of_plane[0]/mag;
  normal_of_plane[1] = normal_of_plane[1]/mag;
  normal_of_plane[2] = normal_of_plane[2]/mag;

  printf("Plane Normal: (%f, %f, %f)\n", normal_of_plane[0], normal_of_plane[1], normal_of_plane[2]);
}

float ray_plane_intersection(int x, int y, Triangle triangle) {
  // calculate the length of each pixel using the WIDTH and HEIGHT of the window
  float x_pixel_length = (aspect_ratio * tan(FOV_RADIANS) - (-aspect_ratio * tan(FOV_RADIANS))) / WIDTH;
  float y_pixel_length = (tan(FOV_RADIANS) - (-tan(FOV_RADIANS))) / HEIGHT;

  // calculate direction vector 
  x_dir = (-aspect_ratio * tan(FOV_RADIANS)) + ((x + 0.5) * x_pixel_length);
  y_dir = tan(FOV_RADIANS) - ((y - 0.5) * y_pixel_length);
  z_dir = -1;

  // normalize the direction
  float mag = magnitide(x_dir, y_dir, z_dir);
  x_dir = x_dir/mag;
  y_dir = y_dir/mag;
  z_dir = z_dir/mag;
  printf("DIRECTION VECTOR TO PLANE: %f, %f, %f\n", x_dir, y_dir, z_dir);

  // Retrieve normal vector of the plane using cross product of 2 rays that lie on the plane
  plane_normal(triangle.v[0], triangle.v[1], triangle.v[2]);
  printf("Normal of Plane: (%f, %f, %f)\n", normal_of_plane[0], normal_of_plane[1], normal_of_plane[2]);

  float numerator = 0;
  float denominator = 0;
  double point_on_plane[3] = {triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]};
  printf("POINT ON PLANE: %f, %f, %f \n", point_on_plane[0], point_on_plane[1], point_on_plane[2]);
  numerator = (point_on_plane[0] * normal_of_plane[0]) + (point_on_plane[1] * normal_of_plane[1]) + (point_on_plane[2] * normal_of_plane[2]);
  denominator = (x_dir * normal_of_plane[0]) + (y_dir * normal_of_plane[1]) + (z_dir * normal_of_plane[2]);
  printf("numerator: %f    denominator: %f\n", numerator, denominator);

  if (denominator == 0) {
    // line and plane do not intersect
    printf("Triangle intersection: Denominator == 0\n");
    return -1;
  } else {
    // return the intersection point
    printf("Triangle intersection value: %f\n", numerator / denominator);
    float t = numerator / denominator;
    if (t <= 0) {
      return -1;
    } else {
      return t;
    }
  }
}

/**
 * Given a root,
 * plug root into ray equation to derive
 * the current intersection point
 **/
float get_intersect_coordinates(float min_root) {
  current_intersection_point[0] = ray_origin[0] + (min_root * x_dir);
  current_intersection_point[1] = ray_origin[1] + (min_root * y_dir);
  current_intersection_point[2] = ray_origin[2] + (min_root * z_dir);
  printf("current intersection point: (%f, %f, %f)\n", current_intersection_point[0], current_intersection_point[1], current_intersection_point[2]);
}

/**
 * Helper function
 * Compute cross product given two rays
 * in the form of two arrays of size 3
 **/
void cross_product(double ray1[3], double ray2[3]) {
  triangle_area_cross_product[0] = (ray1[1] * ray2[2]) - (ray1[2] * ray2[1]);
  triangle_area_cross_product[1] = (ray1[2] * ray2[0]) - (ray1[0] * ray2[2]);
  triangle_area_cross_product[2] = (ray1[0] * ray2[1]) - (ray1[1] * ray2[0]);
  printf("Cross product: (%f, %f, %f)\n", triangle_area_cross_product[0], triangle_area_cross_product[1], triangle_area_cross_product[2]);
}

/**
 * Compute the area of a triangle in 3D
 * using cross product
 **/
float triangle_area(double v1[], double v2[], double v3[]) {
  double ray1[3];
  double ray2[3];

  // v1 - v3
  ray1[0] = v1[0] - v3[0];
  ray1[1] = v1[1] - v3[1];
  ray1[2] = v1[2] - v3[2];

  // v2 - v3
  ray2[0] = v2[0] - v3[0];
  ray2[1] = v2[1] - v3[1];
  ray2[2] = v2[2] - v3[2];

  cross_product(ray1, ray2);
  printf("Cross product: (%f, %f, %f)\n", triangle_area_cross_product[0], triangle_area_cross_product[1], triangle_area_cross_product[2]);

  float mag = magnitide(triangle_area_cross_product[0], triangle_area_cross_product[1], triangle_area_cross_product[2]);
  printf("MAGNITUDE of Cross product: %f\n", mag);

  return 0.5 * mag;
}

float two_d_triangle_area(float x1, float y1, float x2, float y2, float x3, float y3) {
  return abs((x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2))/2.0);
}

/**
 * Use the current intersection point and the current triangle
 * and barycentric coordinates to determine if point is inside
 * the triangle 
 **/
bool is_intersection_inside_triangle(Triangle triangle) {
  float total_area = triangle_area(triangle.v[0].position, triangle.v[1].position, triangle.v[2].position);

  float alpha = triangle_area(current_intersection_point, triangle.v[1].position, triangle.v[2].position) / total_area;

  float beta = triangle_area(triangle.v[0].position, current_intersection_point, triangle.v[2].position) / total_area;

  float gamma = triangle_area(triangle.v[0].position, triangle.v[1].position, current_intersection_point) / total_area;

  printf("ALPHA: %f, BETA:, %f, GAMMA: %f\n", alpha, beta, gamma);

  float sum = alpha + beta + gamma;
  printf("SUM OF AREAS: %f\n", sum);
  if (sum <= 1.001 && sum >= 0.999
          && alpha >= 0 && alpha <= 1
          && beta >= 0 && beta <= 1
          && gamma >= 0 && gamma <= 1) {
    current_alpha = alpha;
    current_beta = beta;
    current_gamma = gamma;
    return true;
  } else {
    return false;
  }
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
      float min_root = ray_sphere_intersection(x, y);

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
          object_type = SPHERE;
        }
      }
    }
  }

  if (num_triangles != 0) {
    for(int i = 0; i < num_triangles; i++) {
      current_triangle = triangles[i];
      float intersection = ray_plane_intersection(x, y, current_triangle);

      if (intersection != -1) {
        get_intersect_coordinates(intersection);

        // use barycentric coordinates to determine if
        // current intersection is inside triangle
        bool isInTriangle = is_intersection_inside_triangle(current_triangle);

        if (isInTriangle) {
          does_intersect = true;
          // continue updating closest intersection point
          float dis = distance_to_cop(current_intersection_point[0], current_intersection_point[1], current_intersection_point[2]);
          printf("distance to current point: %f , vs. closest distance: %f\n", dis, closest_distance);
          if (dis < closest_distance) {
            closest_distance = dis;
            closest_intersection[0] = current_intersection_point[0];
            closest_intersection[1] = current_intersection_point[1];
            closest_intersection[2] = current_intersection_point[2];
            object_type = TRIANGLE;
          }
        }
      }
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
bool shadow_ray_sphere_intersection(Sphere sphere) {
  // calculate direction vector 
  float dir_x = shadow_ray[0];
  float dir_y = shadow_ray[1];
  float dir_z = shadow_ray[2];

  // origin point of the shadow ray is the current intersection point on the sphere
  float origin_x = closest_intersection[0];
  float origin_y = closest_intersection[1];
  float origin_z = closest_intersection[2];

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
           ((origin_z - sphere.position[2]) * (origin_z - sphere.position[2])) -
           (sphere.radius * sphere.radius);
  printf("a: %f    , b: %f     , c: %f  \n", a, b, c);

  float discriminant = (b * b) - (4 * a * c);
  printf("discriminant: %f\n", discriminant);

  // negative discriminant means ray missed sphere, abort rest of calculations
  return discriminant < 0 ? false : true;
}

bool shadow_ray_triangle_intersection(Triangle t) {
  // calculate direction vector 
  float dir_x = shadow_ray[0];
  float dir_y = shadow_ray[1];
  float dir_z = shadow_ray[2];

  // origin point of the shadow ray is the current intersection point on the plane
  float origin_x = closest_intersection[0];
  float origin_y = closest_intersection[1];
  float origin_z = closest_intersection[2];

  // normalize the direction
  float mag = magnitide(dir_x, dir_y, dir_z);
  dir_x = dir_x/mag;
  dir_y = dir_y/mag;
  dir_z = dir_z/mag;

  // Retrieve normal vector of the plane using cross product of 2 rays that lie on the plane
  plane_normal(t.v[0], t.v[1], t.v[2]);
  printf("Normal of Plane: (%f, %f, %f)\n", normal_of_plane[0], normal_of_plane[1], normal_of_plane[2]);

  float numerator = 0;
  float denominator = 0;
  double point_on_plane[3] = {t.v[1].position[0], t.v[1].position[1], t.v[1].position[2]};
  printf("POINT ON PLANE: %f, %f, %f \n", point_on_plane[0], point_on_plane[1], point_on_plane[2]);
  numerator = (point_on_plane[0] * normal_of_plane[0]) + (point_on_plane[1] * normal_of_plane[1]) + (point_on_plane[2] * normal_of_plane[2]);
  denominator = (x_dir * normal_of_plane[0]) + (y_dir * normal_of_plane[1]) + (z_dir * normal_of_plane[2]);
  printf("numerator: %f    denominator: %f\n", numerator, denominator);

  if (denominator == 0) {
    // line and plane do not intersect
    printf("Triangle intersection: Denominator == 0\n");
    return -1;
  } else {
    // return the intersection point
    printf("Triangle intersection value: %f\n", numerator / denominator);
    float t = numerator / denominator;
    if (t <= 0) {
      return -1;
    } else {
      return t;
    }
  }
}

bool sphere_equality(Sphere s1, Sphere s2) {
  printf("s1.position[0] == s2.position[0]: %d\n", s1.position[0] == s2.position[0]);
  printf("s1.position[1] == s2.position[1]: %d\n", s1.position[1] == s2.position[1]);
  printf("s1.position[2] == s2.position[2]: %d\n", s1.position[2] == s2.position[2]);
  printf("s1.color_diffuse[0] == s2.color_diffuse[0]: %d\n", s1.color_diffuse[0] == s2.color_diffuse[0]);
  printf("s1.color_diffuse[1] == s2.color_diffuse[1]: %d\n", s1.color_diffuse[1] == s2.color_diffuse[1]);
  printf("s1.color_diffuse[2] == s2.color_diffuse[2]: %d\n", s1.color_diffuse[2] == s2.color_diffuse[2]);

  printf("s1.color_specular[0] == s2.color_specular[0]: %d\n", s1.color_specular[0] == s2.color_specular[0]);
  printf("s1.color_specular[1] == s2.color_specular[1]: %d\n", s1.color_specular[1] == s2.color_specular[1]);
  printf("s1.color_specular[2] == s2.color_specular[2]: %d\n", s1.color_specular[2] == s2.color_specular[2]);

  printf("s1.radius == s2.radius: %d\n", s1.radius == s2.radius);
  printf("s1.shininess == s2.shininess: %d\n", s1.shininess == s2.shininess);


  return s1.position[0] == s2.position[0] && s1.position[1] == s2.position[1] && s1.position[2] == s2.position[2]
      && s1.radius == s2.radius
      && s1.shininess == s2.shininess
      && s1.color_diffuse[0] == s2.color_diffuse[0] && s1.color_diffuse[1] == s2.color_diffuse[1] && s1.color_diffuse[2] == s2.color_diffuse[2]
      && s1.color_specular[0] == s2.color_specular[0] && s1.color_specular[1] == s2.color_specular[1] && s1.color_specular[2] == s2.color_specular[2];
}

/**
 * Given a shadow ray, do ray-sphere and ray-plane intersection
 * for all the objects in the scene
 * 
 * return: a bool reprenting if the shadow ray is in shadow
 **/
bool shadow_ray_intersection() {
  bool isInShadow = false;

  if (num_spheres != 0) {
    for(int i = 0; i < num_spheres; i++) {
      Sphere s = spheres[i];

      // do not intersect the shadow ray with the current sphere you are on
      // or else you will get all the sphere is in shadow
      printf("Sphere equality: %d\n", sphere_equality(s, current_sphere));
      if (object_type == SPHERE) {
        if (!(sphere_equality(s, current_sphere))) {
          isInShadow = shadow_ray_sphere_intersection(s);
          if (isInShadow) {
            return true;
          }
        }
      } else {
        isInShadow = shadow_ray_sphere_intersection(s);
          if (isInShadow) {
            return true;
          }
      }
    }
  }

  if (num_triangles != 0) {
      for (int i = 0; i < num_triangles; i++) {
        Triangle t = triangles[i];

        float intersection = shadow_ray_triangle_intersection(t);

        if (intersection != -1) {
          get_intersect_coordinates(intersection);

          // use barycentric coordinates to determine if
          // current intersection is inside triangle
          bool isInTriangle = is_intersection_inside_triangle(current_triangle);
          // this means ray intersects with a triangle, thus point is in shadow
          if (isInTriangle) {
            return true;
          }
        }
      }
  }

  // if shadow ray never intersected with an object in the scene
  return false;
}

/**
 * Shadow_ray: unit vector to light
 * normal_vector: surface normal
 * q: the distance to light source (shadow ray's magnitude)
 **/
void phong_illumination(Light current_light) {
  // unit vector to camera
  float v_to_camera[3] = {0, 0, 0};
  v_to_camera[0] = 0 - closest_intersection[0];
  v_to_camera[1] = 0 - closest_intersection[1];
  v_to_camera[2] = 0 - closest_intersection[2];

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
  float illumination_diffuse_r = current_sphere.color_diffuse[0] * light_dot_normal;
  float illumination_diffuse_g = current_sphere.color_diffuse[1] * light_dot_normal;
  float illumination_diffuse_b = current_sphere.color_diffuse[2] * light_dot_normal;


  reflected_vector[0] = (2 * light_dot_normal * normal_vector[0]) - shadow_ray[0];
  reflected_vector[1] = (2 * light_dot_normal * normal_vector[1]) - shadow_ray[1];
  reflected_vector[2] = (2 * light_dot_normal * normal_vector[2]) - shadow_ray[2];
  printf("REFLECTED VECTOR: (%f, %f, %f)\n", reflected_vector[0], reflected_vector[1], reflected_vector[2]);
  printf("REFLECTED VECTOR MAGNITUDE: %f\n", magnitide(reflected_vector[0], reflected_vector[1], reflected_vector[2]));

  float reflected_dot_viewer = (v_to_camera[0] * reflected_vector[0]) + (v_to_camera[1] * reflected_vector[1]) + (v_to_camera[2] * reflected_vector[2]);
  // clamp value to zero if necessary
  if (reflected_dot_viewer < 0) {
    reflected_dot_viewer = 0;
  }
  float illumination_specular_r = current_sphere.color_specular[0] * pow(reflected_dot_viewer, current_sphere.shininess);
  float illumination_specular_g = current_sphere.color_specular[1] * pow(reflected_dot_viewer, current_sphere.shininess);
  float illumination_specular_b = current_sphere.color_specular[2] * pow(reflected_dot_viewer, current_sphere.shininess);

  // Add to the RGB color channels separately
  current_color_at_point[0] += current_light.color[0] * (illumination_diffuse_r + illumination_specular_r);
  current_color_at_point[1] += current_light.color[1] * (illumination_diffuse_g + illumination_specular_g);
  current_color_at_point[2] += current_light.color[2] * (illumination_diffuse_b + illumination_specular_b);
}

/**
 * Interpolate the x,y,z coordinates of the normals given
 * at the triangle vertices, and then normalize the length
 **/ 
void normal_at_triangle_intersection() {
  triangle_intersection_normal[0] = (current_beta * current_triangle.v[0].normal[0]) + (current_alpha * current_triangle.v[1].normal[0] + (current_gamma * current_triangle.v[2].normal[0]));
  triangle_intersection_normal[1] = (current_beta * current_triangle.v[0].normal[1]) + (current_alpha * current_triangle.v[1].normal[1] + (current_gamma * current_triangle.v[2].normal[1]));
  triangle_intersection_normal[2] = (current_beta * current_triangle.v[0].normal[2]) + (current_alpha * current_triangle.v[1].normal[2] + (current_gamma * current_triangle.v[2].normal[2]));
  
  float mag = magnitide(triangle_intersection_normal[0], triangle_intersection_normal[1], triangle_intersection_normal[2]);
  triangle_intersection_normal[0] = triangle_intersection_normal[0] / mag;
  triangle_intersection_normal[1] = triangle_intersection_normal[1] / mag;
  triangle_intersection_normal[2] = triangle_intersection_normal[2] / mag;

  printf("Normal at triangle intersection: %f, %f, %f\n", triangle_intersection_normal[0], triangle_intersection_normal[1], triangle_intersection_normal[2]); 
}

void interpolate_diffuse() {
  interpolated_diffuse[0] = (current_alpha * current_triangle.v[0].color_diffuse[0]) + (current_beta * current_triangle.v[1].color_diffuse[0] + (current_gamma * current_triangle.v[2].color_diffuse[0]));
  interpolated_diffuse[1] = (current_alpha * current_triangle.v[0].color_diffuse[1]) + (current_beta * current_triangle.v[1].color_diffuse[1] + (current_gamma * current_triangle.v[2].color_diffuse[1]));
  interpolated_diffuse[2] = (current_alpha * current_triangle.v[0].color_diffuse[2]) + (current_beta * current_triangle.v[1].color_diffuse[2] + (current_gamma * current_triangle.v[2].color_diffuse[2]));

  float mag = magnitide(interpolated_diffuse[0], interpolated_diffuse[1], interpolated_diffuse[2]);
  interpolated_diffuse[0] = interpolated_diffuse[0] / mag;
  interpolated_diffuse[1] = interpolated_diffuse[1] / mag;
  interpolated_diffuse[2] = interpolated_diffuse[2] / mag;

  printf("DIFFUSE at triangle intersection: %f, %f, %f\n", interpolated_diffuse[0], interpolated_diffuse[1], interpolated_diffuse[2]); 
}

void interpolate_specular() {
  interpolated_specular[0] = (current_alpha * current_triangle.v[0].color_specular[0]) + (current_beta * current_triangle.v[1].color_specular[0] + (current_gamma * current_triangle.v[2].color_specular[0]));
  interpolated_specular[1] = (current_alpha * current_triangle.v[0].color_specular[1]) + (current_beta * current_triangle.v[1].color_specular[1] + (current_gamma * current_triangle.v[2].color_specular[1]));
  interpolated_specular[2] = (current_alpha * current_triangle.v[0].color_specular[2]) + (current_beta * current_triangle.v[1].color_specular[2] + (current_gamma * current_triangle.v[2].color_specular[2]));

  float mag = magnitide(interpolated_specular[0], interpolated_specular[1], interpolated_specular[2]);
  interpolated_specular[0] = interpolated_specular[0] / mag;
  interpolated_specular[1] = interpolated_specular[1] / mag;
  interpolated_specular[2] = interpolated_specular[2] / mag;

  printf("SPECULAR at triangle intersection: %f, %f, %f\n", interpolated_specular[0], interpolated_specular[1], interpolated_specular[2]); 
}

void phong_illumination_triangle(Light current_light) {
  // unit vector to camera
  float v_to_camera[3] = {0, 0, 0};
  v_to_camera[0] = 0 - closest_intersection[0];
  v_to_camera[1] = 0 - closest_intersection[1];
  v_to_camera[2] = 0 - closest_intersection[2];

  float mag = magnitide(v_to_camera[0], v_to_camera[1], v_to_camera[2]);
  v_to_camera[0] = v_to_camera[0]/mag;
  v_to_camera[1] = v_to_camera[1]/mag;
  v_to_camera[2] = v_to_camera[2]/mag;

  float reflected_vector[3] = {0, 0, 0};

  float light_dot_normal = (shadow_ray[0] * triangle_intersection_normal[0]) + (shadow_ray[1] * triangle_intersection_normal[1]) + (shadow_ray[2] * triangle_intersection_normal[2]);
  // clamp value to zero if necessary
  if (light_dot_normal < 0) {
    light_dot_normal = 0;
  }

  interpolate_diffuse();
  float illumination_diffuse_r = interpolated_diffuse[0] * light_dot_normal;
  float illumination_diffuse_g = interpolated_diffuse[1] * light_dot_normal;
  float illumination_diffuse_b = interpolated_diffuse[2] * light_dot_normal;

  reflected_vector[0] = (2 * light_dot_normal * triangle_intersection_normal[0]) - shadow_ray[0];
  reflected_vector[1] = (2 * light_dot_normal * triangle_intersection_normal[1]) - shadow_ray[1];
  reflected_vector[2] = (2 * light_dot_normal * triangle_intersection_normal[2]) - shadow_ray[2];
  printf("REFLECTED VECTOR: (%f, %f, %f)\n", reflected_vector[0], reflected_vector[1], reflected_vector[2]);
  printf("REFLECTED VECTOR MAGNITUDE: %f\n", magnitide(reflected_vector[0], reflected_vector[1], reflected_vector[2]));

  float reflected_dot_viewer = (v_to_camera[0] * reflected_vector[0]) + (v_to_camera[1] * reflected_vector[1]) + (v_to_camera[2] * reflected_vector[2]);
  // clamp value to zero if necessary
  if (reflected_dot_viewer < 0) {
    reflected_dot_viewer = 0;
  }

  interpolate_specular();
  float illumination_specular_r = interpolated_specular[0] * pow(reflected_dot_viewer, current_sphere.shininess);
  float illumination_specular_g = interpolated_specular[1] * pow(reflected_dot_viewer, current_sphere.shininess);
  float illumination_specular_b = interpolated_specular[2] * pow(reflected_dot_viewer, current_sphere.shininess);

  // Add to the RGB color channels separately
  current_color_at_point[0] += current_light.color[0] * (illumination_diffuse_r + illumination_specular_r);
  current_color_at_point[1] += current_light.color[1] * (illumination_diffuse_g + illumination_specular_g);
  current_color_at_point[2] += current_light.color[2] * (illumination_diffuse_b + illumination_specular_b);
}


/**
 * For each light source, fire a shadow ray
 * Determine if light hits the surface
 **/
void fire_shadow_rays() {
  for (int i = 0; i < num_lights; i++) {
    Light current_light = lights[i];
    shadow_ray[0] = current_light.position[0] - closest_intersection[0];
    shadow_ray[1] = current_light.position[1] - closest_intersection[1];
    shadow_ray[2] = current_light.position[2] - closest_intersection[2];

    bool isInShadow = shadow_ray_intersection();
    if (isInShadow) {
      current_color_at_point[0] += 0;
      current_color_at_point[1] += 0;
      current_color_at_point[2] += 0;
    } 
    else {
      // 6) for each unblocked shadow ray, evaluate local Phong model for that light
      
      // normalize the shadow ray
      float mag = magnitide(shadow_ray[0], shadow_ray[1], shadow_ray[2]);
      shadow_ray[0] = shadow_ray[0]/mag;
      shadow_ray[1] = shadow_ray[1]/mag;
      shadow_ray[2] = shadow_ray[2]/mag;

      if (object_type == SPHERE) {
        phong_illumination(current_light); 
      } else {
        // calculate normal at trianlge intersection point
        normal_at_triangle_intersection();
        phong_illumination_triangle(current_light);
      }
    }
  }
  // The final color of the point is the sum of the contributions from all lights, plus the ambient color
  // only add ambient color once
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
  // final color after going through all the lights
  printf("FINAL COLOR: r: %f   g: %f    b: %f\n", current_color_at_point[0], current_color_at_point[1], current_color_at_point[2]);
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
        if (object_type == SPHERE) {
          sphere_normal();

          // 5) for each light source, fire a shadow ray
          fire_shadow_rays();

          // 6) use global color value to color the pixel
          plot_pixel(x, y, current_color_at_point[0] * 255, current_color_at_point[1] * 255, current_color_at_point[2] * 255);
        } 
        else {
          fire_shadow_rays();
          
          plot_pixel(x, y, current_color_at_point[0] * 255, current_color_at_point[1] * 255, current_color_at_point[2] * 255);
        }
      } else {
        plot_pixel(x, y , 255, 255, 255);
      }
      // reset the global color variable to all zeros
      // reset the global closest intersection value for the next pixel we will be evaluating
      current_color_at_point[0] = 0;
      current_color_at_point[1] = 0;
      current_color_at_point[2] = 0;

      closest_intersection[0] = 10000;
      closest_intersection[1] = 10000;
      closest_intersection[2] = 10000;

      current_intersection_point[0] = 0;
      current_intersection_point[1] = 0;
      current_intersection_point[2] = 0;

      does_intersect = false;
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
  glVertex2i(x, y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT - y-1][x][0]=r;
  buffer[HEIGHT - y-1][x][1]=g;
  buffer[HEIGHT - y-1][x][2]=b;
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

  in = pic_alloc(WIDTH, HEIGHT, 3, NULL);
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

  glClearColor(255, 255, 255, 0);
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
