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
	#include "Colors.h" // color names with corresponding RGB values

	using namespace std;

	/**
	 Class representing a single ray.
	*/
	class Ray{
	public:
		glm::vec3 origin; ///< Origin of the ray
		glm::vec3 direction; ///< Direction of the ray
		/**
		 Contructor of the ray
		@param origin Origin of the ray
		@param direction Direction of the ray
		*/
		Ray(glm::vec3 origin, glm::vec3 direction) : origin(origin), direction(direction){
		}
	};


	class Object;

	/**
	 Structure representing the even of hitting an object
	*/
	struct Hit{
		bool hit; ///< Boolean indicating whether there was or there was no intersection with an object
		glm::vec3 normal; ///< Normal vector of the intersected object at the intersection point
		glm::vec3 intersection; ///< Point of Intersection
		float distance; ///< Distance from the origin of the ray to the intersection point
		Object *object; ///< A pointer to the intersected object
	};

	/**
	 General class for the object
	*/
	class Object{
		
	protected:
		glm::mat4 transformationMatrix; ///< Matrix representing the transformation from the local to the global coordinate system
		glm::mat4 inverseTransformationMatrix; ///< Matrix representing the transformation from the global to the local coordinate system
		glm::mat4 normalMatrix; ///< Matrix for transforming normal vectors from the local to the global coordinate system
		
	public:
		glm::vec3 color; ///< Color of the object
		Material material; ///< Structure describing the material of the object
		/** A function computing an intersection, which returns the structure Hit */
		virtual Hit intersect(Ray ray) = 0;

		/** Function that returns the material struct of the object*/
		Material getMaterial(){
			return material;
		}
		/** Function that set the material
		 @param material A structure describing the material of the object
		*/
		void setMaterial(Material material){
			this->material = material;
		}
		/** Functions for setting up all the transformation matrices
		@param matrix The matrix representing the transformation of the object in the global coordinates */
		void setTransformation(glm::mat4 matrix){
				
			transformationMatrix = matrix;
				
			/* ----- Exercise 2 ---------
			Set the two remaining matrices
				
			inverseTransformationMatrix =
			normalMatrix =
				
			*/

			inverseTransformationMatrix = glm::inverse(transformationMatrix); // inverse
			normalMatrix = glm::transpose(inverseTransformationMatrix); // transpose of the inverse

		}
	};

	/**
	 Implementation of the class Object for sphere shape.
	*/
	class Sphere : public Object{
	private:
		float radius; ///< Radius of the sphere
		glm::vec3 center; ///< Center of the sphere

	public:
		/**
		 The constructor of the sphere
		@param radius Radius of the sphere
		@param center Center of the sphere
		@param color Color of the sphere
		*/
		Sphere(float radius, glm::vec3 center, glm::vec3 color) : radius(radius), center(center){
			this->color = color;
		}
		Sphere(float radius, glm::vec3 center, Material material) : radius(radius), center(center){
			this->material = material;
		}
		/** Implementation of the intersection function*/
		Hit intersect(Ray ray){

			glm::vec3 c = center - ray.origin;

			float cdotc = glm::dot(c,c);
			float cdotd = glm::dot(c, ray.direction);

			Hit hit;

			float D = 0;
			if (cdotc > cdotd*cdotd){
				D =  sqrt(cdotc - cdotd*cdotd);
			}
			if(D<=radius){
				hit.hit = true;
				float t1 = cdotd - sqrt(radius*radius - D*D);
				float t2 = cdotd + sqrt(radius*radius - D*D);

				float t = t1;
				if(t<0) t = t2;
				if(t<0){
					hit.hit = false;
					return hit;
				}

				hit.intersection = ray.origin + t * ray.direction;
				hit.normal = glm::normalize(hit.intersection - center);
				hit.distance = glm::distance(ray.origin, hit.intersection);
				hit.object = this;
			}
			else{
				hit.hit = false;
			}
			return hit;
		}
	};

	class Plane : public Object{

	private:
		glm::vec3 normal;
		glm::vec3 point;

	public:
		Plane(glm::vec3 point, glm::vec3 normal) : point(point), normal(normal){
		}
		Plane(glm::vec3 point, glm::vec3 normal, Material material) : point(point), normal(normal){
			this->material = material;
		}
		Hit intersect(Ray ray){
			
			Hit hit; 		 // create a Hit object
			hit.hit = false; // set it as false by default
			
			/*
			
			Excercise 1 - Plane-ray intersection
			
			>>>Idea behind our solution<<<
			A plane is defined by a point p0 and a normal vector n
			A point p is on the plane, if the angle between n and (p - p0) is 90 degrees --> cosine (angle) equals zero --> dot product equals zero
			<n, p - p0> = 0   (1)

			The ray is defined as gamma(t) = o + td, so the intersection point p can be expressed as o +td, with o and d known, but t unknown.
			The equation from (1) can therefore be rewritten as:
			<n, gamma(t) - p0> = 0
			<n, o + td - p0> = 0
			<n, o -p0> + t*<n, d>  (t moved outside the dot product, as it is a scalar)
			t = - <n, o - p0> / <n, d>    (2)

			Because the plane is infinite, the ray always hits unless the ray direction is parallel to the plane.
			The ray is parallel to the plane if the dot product between the ray direction and the plane's normal equals zero.
			We calculate this dot product already in (2), so if it equals zero, we should return the hit as false.
			Also, as the dot product is the denominator of a fraction in (2), it would lead to a division by zero.

			*/

			glm::vec3 o_p0 = ray.origin - point; 	  // prepare o - p0
			float numerator = glm::dot(normal, o_p0); // <n, o - p0>
			float denominator = glm::dot(normal, ray.direction); // <n, d>

			// Check if ray direction is parallel to the plane --> no hit
			if (fabs(denominator) < 1e-6f) { 
				// Note: Instead of exactly zero, we instead use a small number close to zero because of floating point rounding errors
				return hit; // return the hit; remember that hit.hit is false by default
			}

			float t = - numerator / denominator; 	// t = - <n, o - p0> / <n, d>

			// Check if intersection is behind the origin (camera)
			if (t < 0) {
				return hit;
			}

			// Populate the hit structure
			hit.hit = true;
			hit.intersection = ray.origin + t * ray.direction;
			hit.normal = normal; // normal of the plane
			hit.distance = glm::distance(ray.origin, hit.intersection);
			hit.object = this;
			
			return hit;
		}
	};

	class Cone : public Object{
	private:
		Plane *plane;
		
	public:
		Cone(Material material){
			this->material = material;
			plane = new Plane(glm::vec3(0,1,0), glm::vec3(0.0,1,0));
		}
		Hit intersect(Ray ray){
			
			Hit hit;
			hit.hit = false; // set hit as false by default
			
		
			/*  ---- Exercise 2 -----
			
			Implement the ray-cone intersection. Before intersecting the ray with the cone,
			make sure that you transform the ray into the local coordinate system.
			Remember about normalizing all the directions after transformations.
			
			*/
		
			/* If the intersection is found, you have to set all the critical fields in the Hit strucutre
			Remember that the final information about intersection point, normal vector and distance have to be given
			in the global coordinate system.
			
			hit.hit = true;
			hit.object = this;
			hit.intersection =
			hit.normal =
			hit.distance =
			
			*/

			// ------------------- Overview -------------------
			// The main three components of this code are:
			// (a) Transformations (at the beginning and end)
			// (b) Plane-ray intersection for the closed cone's disk
			// (c) Cone-ray intersection
			// ------------------------------------------------

			// ------------- (a) Transformations --------------
			// Strategy for transformations: 
			// 1) Apply the inverse transformation on the ray (origin and direction) --> multiply by inverseTransformationMatrix
			// 2) Compute intersections as usual
			// 3) Apply transformation on the results: intersection point (transformationMatrix) and normal vector (normalMatrix)
			
			// OLD: Remove this after the transformations are implemented
			// glm::vec3 o = ray.origin;
			// glm::vec3 d = ray.direction;

			// Transform the ray's origin and direction:
			// Add the homogenous coordinate (vec3 --> vec4) and then apply the inverseTransformationMatrix
			// Note: Going back to vec3 after the transformation, because it is not needed anymore and operations like normalize expect vec3
			// Direction will need to be normalized, as transformation may impact magnitude
			// Note: Make sure to use local_ray and not ray for the next steps!

			glm::vec3 o = glm::vec3(inverseTransformationMatrix * glm::vec4(ray.origin, 1.0f)); // point: homogenous coordinate = 1
			glm::vec3 d = glm::normalize(glm::vec3(inverseTransformationMatrix * glm::vec4(ray.direction, 0.0f))); // direction: homogenous coordinate = 0
			Ray local_ray(o, d); // Build new ray with the local origin and the local direction

			// Preparing indices for easier readability. e.g. vector[x] takes the vector's first value, which is x
			int x = 0;
			int y = 1;
			int z = 2;

			const float eps = 1e-6f;

			// -- (b) Plane-ray intersection for closed cone --
			// Plane-ray intersection for closed cone
			// The plane covers the open cone at its end at y=1
			// We intersect the plane and then check if the intersection is within the radius
			// If it is, we return the hit, else it is a miss
			Hit plane_hit = plane->intersect(local_ray); // Use the plane class's intersect function
			if (plane_hit.hit == true) {
				glm::vec3 plane_intersection = plane_hit.intersection;
				glm::vec3 plane_normal = plane_hit.normal;
				float distance_to_plane_hit = glm::distance(plane_intersection, o);
				float distance_to_center = plane_intersection[x] * plane_intersection[x] + plane_intersection[z] * plane_intersection[z]; // Disk's center is at (0, 1, 0)
				if (distance_to_center <= 1.0f) { // intersection is on the disk --> valid intersection
					// Transform from local to global coordinates, before populating the hit struct
					glm::vec3 plane_intersection_global = glm::vec3(transformationMatrix * glm::vec4(plane_intersection, 1.0f));
					glm::vec3 plane_normal_global = glm::normalize(glm::vec3(normalMatrix * glm::vec4(plane_normal, 0.0f)));

					plane_hit.hit = true;
					plane_hit.intersection = plane_intersection_global;
					plane_hit.normal = plane_normal_global;
					plane_hit.distance = glm::distance(ray.origin, plane_hit.intersection); // careful, use global ray's origin as we are back to global coordinates
					plane_hit.object = this;
				}
				else {
					plane_hit.hit = false;
				}
			}

			// --------- (c) Cone-ray intersection) -----------

			// ---------------------------------------
			// Compute the quadratic formula to find t
			// ---------------------------------------
			// Implicit form of the cone's surface:
			// The cone starts at the origin (0,0,0) and is centered around the y-axis, with an opening angle of 90 degrees
			// radius r = sqrt(x^2 + z^2)
			// y / r = tan(45 degrees) = 1 --> y = r --> y = sqrt(x^2 + z^2) --> y^2 = x^2 + z^2 --> x^2 + z^2 - y^2 = 0
			// Insert ray: gamma(t) = o + t*d     o and d have x, y and z components (ox, oy, oz), (dx, dy, dz)
			// (ox + t*dx)^2 + (oz + t*dz)^2 - (oy + t*dy)^2 = 0
			// expand and order this term to get the quadratic term (in t):
			// t^2 * (dx^2 + dz^2 - dy^2) + t * 2*(ox*dx + oz*dz - oy*dy) + (ox^2 + oz^2 -oy^2)
			float A = d[x] * d[x] + d[z] * d[z] - d[y] * d[y]; 			// A = (dx^2 + dz^2 - dy^2)
			float B = 2.0f * (o[x] * d[x] + o[z] * d[z] - o[y] * d[y]); // B = 2*(ox*dx + oz*dz - oy*dy)
			float C = o[x] * o[x] + o[z] * o[z] - o[y] * o[y]; 			// C = (ox^2 + oz^2 -oy^2)
			float D = B * B - 4 * A * C;	// Discriminant D = B^2 - 4*A*C

			// The discriminant gives information, if the ray hits the open cone's surface
			// When D is negative, there is no intersection
			// When D equals zero (numerically: very close to zero), there is one intersection
			// When D is positive, there are two intersections
			
			// Negative D: No intersection with the open surface
			if (D < 0) {
				// if plane disk was hit, use that
				if (plane_hit.hit == true) {
					return plane_hit;
				}
				else {
					return hit;
				}
			}

			float t;

			// D equals zero: One intersection
			if (fabs(D) < eps) {
				t = - B / (2 * A);
			}
			// D is positive: Two intersections
			else {
				float t1 = ((-B) - sqrt(D)) / (2 * A);
				float t2 = ((-B) + sqrt(D)) / (2 * A);
				// t1 is closer to the viewer, as sqrt(D) is subtracted
				t = t1;
				// Check if t1 is negative --> behind the viewer --> no hit --> switch to t2
				if (t < 0) {
					t = t2;
				}
				// Check if t2 is negative --> behind the viewer --> no hit --> both don't hit
				if (t < 0) {
					if (plane_hit.hit == true) {
						// if plane disk was hit, use that
						return plane_hit;
					}
					else {
						return hit;
					}
				}
				// if both t1 and t2 are valid (>= 0), we need to check if their intersection is within the cone's y limit
				// For example, it is possible that t1's intersection is on the infinite cone above the y limit, but t2's intersection is within the y limit
				if (t1 >= 0 && t2 >= 0) {
					glm::vec3 intersection1 = o + t1 * d;
					glm::vec3 intersection2 = o + t2 * d;
					if (intersection1[y] >= 0 && intersection1[y] <= 1) { // intersection 1 is within the cone's y limit
						t = t1;
					}
					else if (intersection2[y] >= 0 && intersection2[y] <= 1) { // intersection 2 is within the cone's y limit
						t = t2;
					}
					else {
						// both intersections are outside the cone's y limits --> no hit
						if (plane_hit.hit == true) {
							// if plane disk was hit, use that
							return plane_hit;
						}
						else {
							return hit;
						}
					}
				}

			}
			// Now, the correct t is selected

			glm::vec3 intersection = o + t * d; // Compute intersection
			// Check if the intersection is within the cone's y limit 0 <= y <= 1 --> else, no hit
			if (intersection[y] < 0.0f - eps || intersection[y]> 1.0f + eps) {
				if (plane_hit.hit == true) {
					// if plane disk was hit, use that
					return plane_hit;
				}
				else {
					return hit;
				}
			}

			// Find normal vector: 
			// The implicit form of the cone's surface (x^2 + z^2 - y^2 = 0) dictates, that we are on the cone's surface if this function equals zero
			// The function's gradient (2x, -2y, 2z) shows the direction, in which the function's value increases the most
			// This means, that the gradient shows the (not yet normalized) direction of the normal vector
			// Below, we normalize the gradient. Also, we need to make sure the sign of the normal is correct.
			glm::vec3 gradient(2 * intersection[x], -2 * intersection[y], 2 * intersection[z]); // compute gradient
			gradient = glm::normalize(gradient); // normalize gradient

			if (glm::dot(gradient, d) > 0) {
				gradient = -gradient; // Make sure normal and ray direction point against each other
			}

			// Main results are ready: intersection and normal
			// Transform from local to global coordinates, before populating the hit struct
			glm::vec3 intersection_global = glm::vec3(transformationMatrix * glm::vec4(intersection, 1.0f));
			glm::vec3 normal_global = glm::normalize(glm::vec3(normalMatrix * glm::vec4(gradient, 0.0f)));

			// Populate hit struct
			hit.hit = true;
			hit.intersection = intersection_global;
			hit.normal = normal_global;
			hit.distance = glm::distance(ray.origin, hit.intersection); // careful, use global ray's origin as we are back to global coordinates
			hit.object = this;

			// if the plane disk was hit and that intersection is closer to the viewer, choose that intersection
			if (plane_hit.hit == true && plane_hit.distance < hit.distance) {
				return plane_hit;
			}
			// else (cone intersection was closer), choose the cone intersection
			return hit;
		}
	};

	/**
	 Light class
	*/
	class Light{
	public:
		glm::vec3 position; ///< Position of the light source
		glm::vec3 color; ///< Color/intentisty of the light source
		Light(glm::vec3 position): position(position){
			color = glm::vec3(1.0);
		}
		Light(glm::vec3 position, glm::vec3 color): position(position), color(color){
		}
	};

	vector<Light *> lights; ///< A list of lights in the scene
	glm::vec3 ambient_light(0.05,0.05,0.05);
	vector<Object *> objects; ///< A list of all objects in the scene


	/** Function for computing color of an object according to the Phong Model
	 @param point A point belonging to the object for which the color is computer
	@param normal A normal vector the the point
	@param view_direction A normalized direction from the point to the viewer/camera
	@param material A material structure representing the material of the object
	*/
	glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal, glm::vec3 view_direction, Material material){

		glm::vec3 color(0.0);
		for(int light_num = 0; light_num < lights.size(); light_num++){

			glm::vec3 light_direction = glm::normalize(lights[light_num]->position - point);
			glm::vec3 reflected_direction = glm::reflect(-light_direction, normal);

			float NdotL = glm::clamp(glm::dot(normal, light_direction), 0.0f, 1.0f);
			float VdotR = glm::clamp(glm::dot(view_direction, reflected_direction), 0.0f, 1.0f);

			glm::vec3 diffuse_color = material.diffuse;
			glm::vec3 diffuse = diffuse_color * glm::vec3(NdotL);
			glm::vec3 specular = material.specular * glm::vec3(pow(VdotR, material.shininess));
			
			/*  ---- Exercise 3-----
			
			Include light attenuation due to the distance to the light source.
			
			*/

			// Inverse-square law: intensity of light is inversely proportional to the square of the distance from the source
			// intensity_object = intensity_emitted  / distance^2
			// The attenuation factor introduces additional terms to 'fine-tune' this:
			// attenuation = 1 / (a + b * r + c * r^2)
			// The three terms influence the attenuation quadratically, linearly and constantly

			float a = 0.5f;
			float b = 0.09f;
			float c = 0.002f;
			float distance = glm::length(lights[light_num]->position - point);
			float attenuation = 1.0f / (a + b * distance + c * distance * distance);
			
			// old
			// color += lights[light_num]->color * (diffuse + specular);

			// updated
			color += lights[light_num]->color * (diffuse + specular) * attenuation; // multiplied attenuation factor
		}
		color += ambient_light * material.ambient;
		color = glm::clamp(color, glm::vec3(0.0), glm::vec3(1.0));
		return color;
	}

	/**
	 Functions that computes a color along the ray
	@param ray Ray that should be traced through the scene
	@return Color at the intersection point
	*/
	glm::vec3 trace_ray(Ray ray){

		Hit closest_hit;

		closest_hit.hit = false;
		closest_hit.distance = INFINITY;

		for(int k = 0; k<objects.size(); k++){
			Hit hit = objects[k]->intersect(ray);
			if(hit.hit == true && hit.distance < closest_hit.distance)
				closest_hit = hit;
		}

		glm::vec3 color(0.0);

		if(closest_hit.hit){
			color = PhongModel(closest_hit.intersection, closest_hit.normal, glm::normalize(-ray.direction), closest_hit.object->getMaterial());
		}else{
			color = glm::vec3(0.0, 0.0, 0.0);
		}
		return color;
	}
	/**
	 Function defining the scene
	*/
	void sceneDefinition (){

		/*  ---- All Exercises -----
		
		Modify the scene definition according to the exercises
		
		*/

		Material green_diffuse;
		green_diffuse.ambient = glm::vec3(Colors::Green);
		green_diffuse.diffuse = glm::vec3(Colors::Green);

		Material red_specular;
		red_specular.ambient = glm::vec3(Colors::Red);
		red_specular.diffuse = glm::vec3(Colors::Red);
		red_specular.specular = glm::vec3(0.5);
		red_specular.shininess = 10.0;

		Material blue_specular;
		blue_specular.ambient = glm::vec3(Colors::Blue);
		blue_specular.diffuse = glm::vec3(Colors::Blue);
		blue_specular.specular = glm::vec3(0.6);
		blue_specular.shininess = 100.0;
		
		
		objects.push_back(new Sphere(1.0, glm::vec3(1,-2,8), blue_specular));
		objects.push_back(new Sphere(0.5, glm::vec3(-1,-2.5,6), red_specular));
		// objects.push_back(new Sphere(1.0, glm::vec3(2,-2,6), green_diffuse));
			
		lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(0.8)));
		lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.8)));
		lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.8)));

		// --------------------------------------------------------
		// Exercise 1: Add six planes to the scene to build a "box"
		// --------------------------------------------------------

		// Materials for the planes
		// See Colors.h for a list of predefined colors | Declaration of generative AI: The list of colors was generated by ChatGPT
		// Defined a new function "make_box_material" in Material.h
		Material darkgreen_wall = make_box_material(Colors::DarkGreen);
		Material teal_wall = make_box_material(Colors::Teal);
		Material purple_wall = make_box_material(Colors::Purple);
		Material lightblue_wall = make_box_material(Colors::LightBlue);
		Material grey_wall = make_box_material(Colors::LightGray);
		Material brown_wall = make_box_material(Colors::Brown);

		// Construct plains using a point and direction
		// Using the 8 corner points, we can easily construct normal directions for the planes
		// Example: The "roof" of the box can be described by a point in the top-right corner (as the point) and 
		// the normalized vector between this point and its counterpart on the bottom right.
		glm::vec3 behind_top_right(15.0f, 27.0f, 30.0f);
		glm::vec3 front_top_right(15.0f, 27.0f, -0.01f);
		glm::vec3 behind_bottom_right(15.0f, -3.0f, 30.0f);
		glm::vec3 behind_top_left(-15.0f, 27.0f, 30.0f);
		glm::vec3 behind_bottom_left(-15.0f, -3.0f, 30.0f);

		objects.push_back(new Plane(behind_top_right, glm::normalize(front_top_right - behind_top_right), purple_wall)); //front
		objects.push_back(new Plane(behind_top_right, glm::normalize(behind_bottom_right - behind_top_right), brown_wall)); // up
		objects.push_back(new Plane(behind_top_right, glm::normalize(behind_top_left - behind_top_right), teal_wall)); // right
		objects.push_back(new Plane(behind_bottom_left, glm::normalize(behind_top_left - behind_bottom_left), grey_wall)); // bottom
		objects.push_back(new Plane(behind_bottom_left, glm::normalize(behind_bottom_right - behind_bottom_left), lightblue_wall)); // left
		objects.push_back(new Plane(front_top_right, glm::normalize(behind_top_right - front_top_right), darkgreen_wall)); // behind

		// --------------------------------------------------------
		// Exercise 2: Add cones to the scene
		// --------------------------------------------------------

		// Yellow cone
		// Material - requirements: yellow color, highly specular
		Material yellow_highly_specular;
		yellow_highly_specular.ambient = glm::vec3(Colors::Yellow);
		yellow_highly_specular.diffuse = glm::vec3(Colors::Yellow);
		yellow_highly_specular.specular = glm::vec3(1.2);
		yellow_highly_specular.shininess = 150.0;

		Cone* cone_yellow = new Cone(yellow_highly_specular); // make new cone
		// Build matrices for transformations
		glm::mat4 S_yellow = glm::scale(glm::vec3(3.0f, 12.0f, 3.0f)); // scale y by 12 (height 1 --> 12), x & z by 3 (disk radius 1 --> 3)
		glm::mat4 R_yellow = glm::rotate(glm::radians(180.0f), glm::vec3(1,0,0)); // rotate 180 degrees around x-axis
		glm::mat4 T_yellow = glm::translate(glm::vec3(5.0f, 9.0f, 14.0f)); 		 // translate (move)
		glm::mat4 M_yellow = T_yellow * R_yellow * S_yellow; // combine transformations; first scale, then rotate, then translate
		cone_yellow->setTransformation(M_yellow); // 
		objects.push_back(cone_yellow);

		Cone* cone_green = new Cone(green_diffuse);
		glm::mat4 S_green = glm::scale(glm::vec3(1.0f, 3.0f, 1.0f)); // scale y by 3 (height 1 --> 3), x & z by 1 (disk radius 1 --> 1)
		glm::mat4 R_green = glm::rotate(glm::radians(71.5f), glm::vec3(0,0,1)); // rotate around z-axis by atan(3/1) = approx. 71.5 degrees
		glm::mat4 T_green = glm::translate(glm::vec3(6.0f, -3.0f, 7.0f)); 		// translate (move)
		glm::mat4 M_green = T_green * R_green * S_green;
		cone_green->setTransformation(M_green);
		objects.push_back(cone_green);

	}
	glm::vec3 toneMapping(glm::vec3 intensity){

		/*  ---- Exercise 3-----
		
		Implement a tonemapping strategy and gamma correction for a correct display.
		
		*/

		// Tone mapping (simple)
		// Techniques discussed in the lecture: simple scaling, logarithm function, power function, exponential scaling
		// Simple scaling: scaled_intensity = alpha * intensity   Choose alpha, e.g. alpha = 1 / I_max
		// Exponential: scaled_intensity = 1 - exp(-a * intensity)

		// Another (apparently) widely used simple tone mapping is the Reinhard tone mapping (Source: https://en.wikipedia.org/wiki/Tone_mapping)
		// scaled_intensity = intensity / (intensity + 1)

		// Exponential tone mapping
		// float a = 1.0f;
		// intensity = 1.0f - glm::exp(-a * intensity);

		// Reinhard tone mapping
		intensity = intensity / (intensity + glm::vec3(1.0f));

		
		// Gamma correction: display_intensity = intensity^(1 / gamma)
		// Choose gamma value: 1.8 <= gamma <= 2.4
		float gamma = 2.2f;
		intensity = glm::pow(intensity, glm::vec3(1.0f / gamma)); // gamma correction
		
		return intensity;
	}
	int main(int argc, const char * argv[]) {

		clock_t t = clock(); // variable for keeping the time of the rendering

		int width = 1024; //width of the image
		int height = 768; // height of the image
		float fov = 90; // field of view

		sceneDefinition(); // Let's define a scene

		Image image(width,height); // Create an image where we will store the result

		float s = 2*tan(0.5*fov/180*M_PI)/width;
		float X = -s * width / 2;
		float Y = s * height / 2;

		for(int i = 0; i < width ; i++)
			for(int j = 0; j < height ; j++){

				float dx = X + i*s + s/2;
				float dy = Y - j*s - s/2;
				float dz = 1;

				glm::vec3 origin(0, 0, 0);
				glm::vec3 direction(dx, dy, dz);
				direction = glm::normalize(direction);

				Ray ray(origin, direction);

				image.setPixel(i, j, glm::clamp(toneMapping(trace_ray(ray)), glm::vec3(0.0), glm::vec3(1.0)));

			}

		t = clock() - t;
		cout<<"It took " << ((float)t)/CLOCKS_PER_SEC<< " seconds to render the image."<< endl;
		cout<<"I could render at "<< (float)CLOCKS_PER_SEC/((float)t) << " frames per second."<<endl;

		// Writing the final results of the rendering
		if (argc == 2){
			image.writeImage(argv[1]);
		}else{
			image.writeImage("./result.ppm");
		}

		
		return 0;
	}
