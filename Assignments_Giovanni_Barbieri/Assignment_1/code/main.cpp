/**
@file main.cpp
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"
#include "Image.h"
#include "Material.h"

using namespace std;

/**
 Class representing a single ray.
 */
class Ray
{
public:
	glm::vec3 origin;	 ///< Origin of the ray
	glm::vec3 direction; ///< Direction of the ray
						 /**
						  Contructor of the ray
						  @param origin Origin of the ray
						  @param direction Direction of the ray
						  */
	Ray(glm::vec3 origin, glm::vec3 direction) : origin(origin), direction(direction)
	{
	}
};

class Object;

/**
 Structure representing the even of hitting an object
 */
struct Hit
{
	bool hit;				///< Boolean indicating whether there was or there was no intersection with an object
	glm::vec3 normal;		///< Normal vector of the intersected object at the intersection point
	glm::vec3 intersection; ///< Point of Intersection
	float distance;			///< Distance from the origin of the ray to the intersection point
	Object *object;			///< A pointer to the intersected object
};

/**
 General class for the object
 */
class Object
{

public:
	glm::vec3 color;   ///< Color of the object
	Material material; ///< Structure describing the material of the object
					   /** A function computing an intersection, which returns the structure Hit */
	virtual Hit intersect(Ray ray) = 0;

	/** Function that returns the material struct of the object*/
	Material getMaterial()
	{
		return material;
	}
	/** Function that set the material
	 @param material A structure describing the material of the object
	*/
	void setMaterial(Material material)
	{
		this->material = material;
	}
};

/**
 Implementation of the class Object for sphere shape.
 */
class Sphere : public Object
{
private:
	float radius;	  ///< Radius of the sphere
	glm::vec3 center; ///< Center of the sphere

public:
	/**
	 The constructor of the sphere
	 @param radius Radius of the sphere
	 @param center Center of the sphere
	 @param color Color of the sphere
	 */
	Sphere(float radius, glm::vec3 center, glm::vec3 color) : radius(radius), center(center)
	{
		this->color = color;
	}
	Sphere(float radius, glm::vec3 center, Material material) : radius(radius), center(center)
	{
		this->material = material;
	}
	/** Implementation of the intersection function*/
	Hit intersect(Ray ray)
	{

		Hit hit;
		hit.hit = false;

		/* ------------------ Exercise 1 -------------------

		Place for your code: ray-sphere intersection. Remember to set all the fields of the hit structure:

		 hit.intersection =
		 hit.normal =
		 hit.distance =
		 hit.object = this;

		------------------------------------------------- */
		glm::vec3 o = ray.origin;
		glm::vec3 d = ray.direction;
		glm::vec3 C = center;

		glm::vec3 c = C - o; // vector from center of sphere to origin of ray

		// notes from pdf on ray-sphere intersection:
		// (slightly modded in case origin is not 0,0,0)

		//  a = <c,d>
		//  D^2 = ||c||^2 - a^2

		float a = glm::dot(c, d);	// // distance along ray direction to closest approach of sphere center
		float c2 = glm::dot(c, c);	// this is ||c||^2, represents the squared length of vector C from origin to center of sphere
		float D = sqrt(c2 - a * a); // the squared distance from the center of the sphere to the closest point on the ray

		// D < r -> 2 intersections
		// D = r -> 1 intersection
		// D > r -> no intersection

		// now we can determine the number of intersections
		if (D > radius)
		{
			// no intersection
			hit.hit = false;
			return hit;
		}
		else if (fabs(D - radius) < 1e-6f) // this comparison was suggested by chatGPT to avoid precision issues with ==
		{
			// 1 intersection
			float t = a; // distance from origin to intersection point along ray direction
			if (t < 0)
			{
				hit.hit = false;
				return hit; // intersection is behind the ray origin
			}
			hit.hit = true;
			hit.intersection = o + t * d; // intersection point
			hit.distance = glm::length(hit.intersection - o);
			hit.normal = glm::normalize(hit.intersection - C); // normal at intersection point
			hit.object = this;
		}
		else
		{
			// 2 intersections
			// b = sqrt(radius^2 - D^2)
			float b = sqrt(radius * radius - D * D);
			float t1 = a - b; // distance to first intersection point
			float t2 = a + b; // distance to second intersection point
			float t;
			if (t1 >= 0)
			{
				t = t1; // take the first intersection if it's in front of the ray origin
			}
			else if (t2 >= 0)
			{
				t = t2; // otherwise take the second intersection if it's in front of the ray origin
			}
			else
			{
				hit.hit = false;
				return hit; // both intersections are behind the ray origin
			}
			hit.hit = true;
			hit.intersection = o + t * d; // intersection point
			hit.distance = glm::length(hit.intersection - o);
			hit.normal = glm::normalize(hit.intersection - C); // normal at intersection point
			hit.object = this;
		}
		return hit;
	}
};

/**
 Light class
 */
class Light
{
public:
	glm::vec3 position; ///< Position of the light source
	glm::vec3 color;	///< Color/intentisty of the light source
	Light(glm::vec3 position) : position(position)
	{
		color = glm::vec3(1.0);
	}
	Light(glm::vec3 position, glm::vec3 color) : position(position), color(color)
	{
	}
};

vector<Light *> lights; ///< A list of lights in the scene
glm::vec3 ambient_light(0.1, 0.1, 0.1);
vector<Object *> objects; ///< A list of all objects in the scene

/** Function for computing color of an object according to the Phong Model
 @param point A point belonging to the object for which the color is computer
 @param normal A normal vector the the point
 @param view_direction A normalized direction from the point to the viewer/camera
 @param material A material structure representing the material of the object
*/
glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal, glm::vec3 view_direction, Material material)
{

	glm::vec3 color(0.0);

	/* ------------------Excercise 3--------------------

	 Phong model.
	 Your code should implement a loop over all the lightsourses in the array lights and agredate the contribution of each of them to the final color of the object.
	 Outside of the loop add also the ambient component from ambient_light.

	 ------------------------------------------------- */

	// Loop over each light source
	for (Light *light : lights)
	{
		// we calculate the vectors L, V, R becuase they are needed for the diffuse and specular components
		glm::vec3 L = glm::normalize(light->position - point); // light direction
		glm::vec3 V = glm::normalize(view_direction);		   // view direction
		glm::vec3 R = glm::reflect(-L, normal);				   // reflection direction (equivalent to r from lesson)

		// Diffuse component
		// from lesson Id = pd * dot(n,l) * I
		// need to clamp the dot product to avoid negative values otheriwse I get more shading than expected
		glm::vec3 diffuse = material.diffuse * glm::max(glm::dot(normal, L), 0.0f) * light->color;
		color += diffuse;

		// Specular component
		// from lesson Is = ps * cos(alpha)^k * I = ps * dot(r,v)^k * I

		// with reference to lecture 03 slide 14:
		float cos_alpha = glm::max(glm::dot(R, V), 0.0f);

		glm::vec3 specular = material.specular * powf(cos_alpha, material.shininess) * light->color;
		color += specular;
	}

	// From lesson -> ambient term = pa * Ia
	// we add it directly to the final color
	color += ambient_light * material.ambient;

	// The final color has to be clamped so the values do not go beyond 0 and 1.
	color = glm::clamp(color, glm::vec3(0.0), glm::vec3(1.0));
	return color;
}

/**
 Functions that computes a color along the ray
 @param ray Ray that should be traced through the scene
 @return Color at the intersection point
 */
glm::vec3 trace_ray(Ray ray)
{

	Hit closest_hit;

	closest_hit.hit = false;
	closest_hit.distance = INFINITY;

	for (int k = 0; k < objects.size(); k++)
	{
		Hit hit = objects[k]->intersect(ray);
		if (hit.hit == true && hit.distance < closest_hit.distance)
			closest_hit = hit;
	}

	glm::vec3 color(0.0);

	if (closest_hit.hit)
	{
		/* ------------------Excercise 3--------------------

		 Use the second line when you implement PhongModel function - Exercise 3

		 ------------------------------------------------- */
		color = closest_hit.object->color;
		color = PhongModel(closest_hit.intersection, closest_hit.normal, glm::normalize(-ray.direction), closest_hit.object->getMaterial());
	}
	else
	{
		color = glm::vec3(0.0, 0.0, 0.0);
	}
	return color;
}
/**
 Function defining the scene
 */
void sceneDefinition()
{

	/* ------------------Excercise 2--------------------

	Place for your code: additional sphere

	------------------------------------------------- */

	/* ------------------Excercise 3--------------------

	 Add here all the objects to the scene. Remember to add them using the constructor for the sphere with material structure.
	 You will also need to define the materials.
	 Example of adding one sphere:

	 Material red_specular;
	 red_specular.diffuse = glm::vec3(1.0f, 0.3f, 0.3f);
	 red_specular.ambient = glm::vec3(0.01f, 0.03f, 0.03f);
	 red_specular.specular = glm::vec3(0.5);
	 red_specular.shininess = 10.0;

	 objects.push_back(new Sphere(0.5, glm::vec3(-1,-2.5,6), red_specular));

	 Remember also about adding some lights. For example a white light of intensity 0.4 and position in (0,26,5):

	 lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(0.4)));

	------------------------------------------------- */

	Material red_specular;
	red_specular.diffuse = glm::vec3(1.0f, 0.3f, 0.3f);
	red_specular.ambient = glm::vec3(0.01f, 0.03f, 0.03f);
	red_specular.specular = glm::vec3(0.5);
	red_specular.shininess = 10.0;

	objects.push_back(new Sphere(0.5, glm::vec3(-1, -2.5, 6), red_specular));

	Material green_specular;
	green_specular.diffuse = glm::vec3(0.7, 0.9, 0.7);
	green_specular.ambient = glm::vec3(0.07, 0.09, 0.07);
	green_specular.specular = glm::vec3(0.0);
	green_specular.shininess = 0.0;

	objects.push_back(new Sphere(1.0, glm::vec3(2, -2, 6), green_specular));

	Material blue_specular;
	blue_specular.diffuse = glm::vec3(0.7, 0.7, 1.0);
	blue_specular.ambient = glm::vec3(0.07, 0.07, 0.1);
	blue_specular.specular = glm::vec3(0.6);
	blue_specular.shininess = 100.0;

	objects.push_back(new Sphere(1.0, glm::vec3(1, -2, 8), blue_specular));

	lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(0.4)));
	lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.4)));
	lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.4)));
}

int main(int argc, const char *argv[])
{
	clock_t t = clock(); // variable for keeping the time of the rendering

	int width = 1024; // width of the image
	int height = 768; // height of the image
	float fov = 90;	  // field of view

	sceneDefinition(); // Let's define a scene

	Image image(width, height); // Create an image where we will store the result

	/* ------------------ Exercise 1 -------------------

	Place for your code: Loop over pixels to form and traverse the rays through the scene

	------------------------------------------------- */

	/* ------------------ Exercise 1 -------------------

			Place for your code: ray definition for pixel (i,j), ray traversal

			 ------------------------------------------------- */

	// Definition of the ray
	// glm::vec3 origin(0, 0, 0);
	// glm::vec3 direction(?, ?, ?);               // fill in the correct values
	// direction = glm::normalize(direction);

	// Ray ray(origin, direction);  // ray traversal

	// image.setPixel(i, j, trace_ray(ray));

	// prof notes
	//  (X,Y) is the top left corner of the image plane
	//  X = -width/2 * s
	//  Y = height/2 * s
	//  s = 2*tan(fov/2)/width -> size of a pixel

	float s = (2 * tan((fov / 2) * M_PI / 180)) / width; // chatgpt helped me with this line with the conversion to radians
	float X = -width / 2.0f * s;
	float Y = height / 2.0f * s;

	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
		{
			glm::vec3 origin(0, 0, 0);

			// prof notes -> dx = X + i*s + 0.5s   dy = Y - j*s - 0.5s   dz = 1
			float dx = X + i * s + 0.5f * s;
			float dy = Y - j * s - 0.5f * s;
			float dz = 1.0;

			glm::vec3 direction = glm::normalize(glm::vec3(dx, dy, dz));
			Ray ray(origin, direction);			  // creating ray with defiend origin and direction
			image.setPixel(i, j, trace_ray(ray)); // tracing the ray and setting the pixel color
		}

	t = clock() - t;
	cout << "It took " << ((float)t) / CLOCKS_PER_SEC << " seconds to render the image." << endl;
	cout << "I could render at " << (float)CLOCKS_PER_SEC / ((float)t) << " frames per second." << endl;

	// Writing the final results of the rendering
	if (argc == 2)
	{
		image.writeImage(argv[1]);
	}
	else
	{
		image.writeImage("./result.ppm");
	}

	return 0;
}
