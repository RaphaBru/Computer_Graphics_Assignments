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
#include <sstream>
#include <array>

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

protected:
	glm::mat4 transformationMatrix;		   ///< Matrix representing the transformation from the local to the global coordinate system
	glm::mat4 inverseTransformationMatrix; ///< Matrix representing the transformation from the global to the local coordinate system
	glm::mat4 normalMatrix;				   ///< Matrix for transforming normal vectors from the local to the global coordinate system

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
	/** Functions for setting up all the transformation matrices
	@param matrix The matrix representing the transformation of the object in the global coordinates */
	void setTransformation(glm::mat4 matrix)
	{

		transformationMatrix = matrix;

		inverseTransformationMatrix = glm::inverse(matrix);
		normalMatrix = glm::transpose(inverseTransformationMatrix);
	}
};
class Plane : public Object
{

private:
	glm::vec3 normal;
	glm::vec3 point;

public:
	Plane(glm::vec3 point, glm::vec3 normal) : point(point), normal(normal)
	{
	}
	Plane(glm::vec3 point, glm::vec3 normal, Material material) : point(point), normal(normal)
	{
		this->material = material;
	}
	Hit intersect(Ray ray)
	{

		Hit hit;
		hit.hit = false;

		float DdotN = glm::dot(ray.direction, normal);
		if (DdotN < 0)
		{

			float PdotN = glm::dot(point - ray.origin, normal);
			float t = PdotN / DdotN;

			if (t > 0)
			{
				hit.hit = true;
				hit.normal = normal;
				hit.distance = t;
				hit.object = this;
				hit.intersection = t * ray.direction + ray.origin;
			}
		}

		return hit;
	}
};
// Assignment 3: Class for ray triangle intersection
class Triangle : public Object
{
private:
	std::array<glm::vec3, 3> vertices;
	std::array<glm::vec3, 3> vertex_normals;
	bool has_vertex_normals = false;

public:
	// Flat shading
	Triangle(std::array<glm::vec3, 3> v, Material material) : vertices(v)
	{
		this->material = material;
		// triangles live in mesh-local; leave Object's transforms as identity
		setTransformation(glm::mat4(1.0f));
	}

	// Smooth shading (per-vertex normals)
	Triangle(std::array<glm::vec3, 3> v, std::array<glm::vec3, 3> vn, Material material)
		: vertices(v), vertex_normals(vn), has_vertex_normals(true)
	{
		this->material = material;
		setTransformation(glm::mat4(1.0f));
	}

	Hit intersect(Ray ray) override
	{
		// ray and triangle are in the SAME (mesh-local) space.
		// Ray-plane intersection: n = (b-a) x (c-a)
		const glm::vec3 &a = vertices[0];
		const glm::vec3 &b = vertices[1];
		const glm::vec3 &c = vertices[2];

		glm::vec3 ab = b - a;
		glm::vec3 ac = c - a;
		glm::vec3 n = glm::cross(ab, ac);
		float nlen2 = glm::dot(n, n);
		Hit out{};
		out.hit = false;

		// Degenerate?
		if (nlen2 <= 0.0f)
			return out;

		float denom = glm::dot(n, ray.direction);
		if (fabs(denom) < 1e-8f)
			return out; // parallel to plane

		float t = glm::dot(n, a - ray.origin) / denom;
		if (t <= 0.0f)
			return out; // behind or at origin

		glm::vec3 p = ray.origin + t * ray.direction; // intersection point on plane

		// Barycentric test via subtriangle areas (using cross products)
		glm::vec3 n1 = glm::cross(b - p, c - p);
		glm::vec3 n2 = glm::cross(c - p, a - p);
		glm::vec3 n3 = glm::cross(a - p, b - p);

		float w1 = glm::dot(n, n1) / nlen2;
		float w2 = glm::dot(n, n2) / nlen2;
		float w3 = glm::dot(n, n3) / nlen2;

		// Inside if all >= 0 (no epsilon for simplicity; add small -1e-6f if needed)
		if (w1 < 0.0f || w2 < 0.0f || w3 < 0.0f)
			return out;

		glm::vec3 normal;
		if (has_vertex_normals)
		{
			normal = glm::normalize(w1 * vertex_normals[0] + w2 * vertex_normals[1] + w3 * vertex_normals[2]);
		}
		else
		{
			normal = glm::normalize(n);
		}

		// Fill hit in LOCAL coords; distance is t along THIS ray
		out.hit = true;
		out.object = this;	  // not used by Mesh (Mesh will overwrite), but fine
		out.intersection = p; // local
		out.normal = normal;  // local
		out.distance = t;	  // IMPORTANT: t in local space

		return out;
	}
};

/**
 Class representing a mesh object made up of triangles
 */
class Mesh : public Object
{
public:
	std::vector<Triangle *> triangles;			   // pointers to all triangles of the mesh
	std::vector<glm::vec3> vertices;			   // vertices from v lines
	std::vector<glm::vec3> normals;				   // normals from vn lines (if provided)
	std::vector<std::array<int, 3>> faces;		   // triangles from f lines
	std::vector<std::array<int, 3>> normalIndices; // indices of the three vertex normals of a triangle
	bool hasNormals = false;					   // true if vertex normals are provided, else false

	Mesh(std::string filename, Material m)
	{
		this->material = m;
		loadOBJ(filename, m);
	}

	void loadOBJ(const std::string &filepath, Material m)
	{
		std::ifstream in(filepath);

		// Check if the file is open, if not print an error message and return
		if (!in.is_open())
		{
			std::cerr << "Error: could not open OBJ file " << filepath << std::endl;
			return;
		}

		// initialize temporary storage for positions and normals
		std::vector<glm::vec3> posnTmp;
		std::vector<glm::vec3> nrmTmp;
		std::string row;

		// helpers to parse integers and split strings
		auto parseInt = [](const std::string &s) -> int
		{
			int v = 0;
			std::istringstream is(s);
			is >> v;
			return v;
		};

		auto splitBySlash = [](const std::string &token) -> std::vector<std::string>
		{
			std::vector<std::string> parts;
			std::string cur;
			for (char ch : token)
			{
				if (ch == '/')
				{
					parts.push_back(cur);
					cur.clear();
				}
				else
					cur.push_back(ch);
			}
			parts.push_back(cur);
			return parts;
		};

		// while we read lines from the file
		while (std::getline(in, row))
		{
			if (row.size() < 2)
				continue;

			// here we check for a space after the first character to distinguish between "v " and "vn "
			if (row[0] == 'v' && row[1] == ' ')
			{
				// vertex position
				std::istringstream ss(row.substr(2)); // we need to skip the "v " part of the line
				glm::vec3 p;
				ss >> p.x >> p.y >> p.z; // we read the three components from the rest of ss stream
										 // and store them in a glm::vec3
										 // finally we add the position to the positions vector
				posnTmp.push_back(p);
			}
			else if (row.rfind("vn ", 0) == 0)
			{
				// vertex normal
				std::istringstream ss(row.substr(3)); // skip "vn "
				glm::vec3 n;
				ss >> n.x >> n.y >> n.z; // read normal components
										 // normalize and store
				nrmTmp.push_back(glm::normalize(n));
			}
			else if (row.rfind("f ", 0) == 0)
			{
				std::istringstream ss(row.substr(2));
				std::string a, b, c;
				ss >> a >> b >> c;
				if (a.empty() || b.empty() || c.empty())
					continue;

				auto parseTriplet = [&](const std::string &tok, int &vi, int &vti, int &vni)
				{
					auto parts = splitBySlash(tok);
					// OBJ indices are 1-based; we’ll store 0-based later
					vi = parts.size() >= 1 && !parts[0].empty() ? parseInt(parts[0]) : 0;
					vti = parts.size() >= 2 && !parts[1].empty() ? parseInt(parts[1]) : 0;
					vni = parts.size() >= 3 && !parts[2].empty() ? parseInt(parts[2]) : 0;
				};

				int av = 0, bv = 0, cv = 0, at = 0, bt = 0, ct = 0, an = 0, bn = 0, cn = 0;

				// Handle v/vt/vn, v//vn, or v
				if (a.find('/') != std::string::npos ||
					b.find('/') != std::string::npos ||
					c.find('/') != std::string::npos)
				{
					parseTriplet(a, av, at, an);
					parseTriplet(b, bv, bt, bn);
					parseTriplet(c, cv, ct, cn);

					if (av > 0 && bv > 0 && cv > 0)
					{
						faces.push_back({av - 1, bv - 1, cv - 1});
						if (an > 0 && bn > 0 && cn > 0)
						{
							normalIndices.push_back({an - 1, bn - 1, cn - 1});
							hasNormals = true;
						}
					}
				}
				else
				{
					// positions only: "f v v v"
					int iv1 = parseInt(a), iv2 = parseInt(b), iv3 = parseInt(c);
					if (iv1 > 0 && iv2 > 0 && iv3 > 0)
					{
						faces.push_back({iv1 - 1, iv2 - 1, iv3 - 1});
					}
				}
			}
		}
		// move parsed data into member storage
		vertices = std::move(posnTmp);
		normals = std::move(nrmTmp);
		// create Triangle objects now (identity transform; you can set a real one later)
		buildTriangles(m, glm::mat4(1.0f));
		// we must close the file when we are done
		in.close();
	}

	Hit intersect(Ray ray) override
	{
		// 1) World → Mesh-local
		glm::vec3 oL = glm::vec3(inverseTransformationMatrix * glm::vec4(ray.origin, 1.0f));
		glm::vec3 dL = glm::normalize(glm::vec3(inverseTransformationMatrix * glm::vec4(ray.direction, 0.0f)));
		Ray rLocal(oL, dL);

		Hit best{};
		best.hit = false;
		best.distance = INFINITY;

		// 2) Test all triangles (they expect mesh-local rays)
		for (auto *tri : triangles)
		{
			Hit h = tri->intersect(rLocal);
			if (h.hit && h.distance < best.distance)
				best = h;
		}

		if (!best.hit)
			return best;

		// 3) Mesh-local → World (for the chosen hit only)
		glm::vec3 pW = glm::vec3(transformationMatrix * glm::vec4(best.intersection, 1.0f));
		glm::vec3 nW = glm::normalize(glm::vec3(normalMatrix * glm::vec4(best.normal, 0.0f)));

		best.intersection = pW;
		best.normal = nW;
		best.distance = glm::length(pW - ray.origin); // true world-space distance
		best.object = this;							  // ensure Phong uses Mesh material

		return best;
	}

	// builds triangle objects from the parsed vertex and face data
	void buildTriangles(Material material, const glm::mat4 &transformation_matrix)
	{
		// clear previous triangles (if any) and pre-allocate space for all faces
		triangles.clear();
		triangles.reserve(faces.size());

		// loop over all faces
		for (size_t i = 0; i < faces.size(); ++i)
		{
			// get the vertex indices for this face (assumed 0-based)
			const auto &face = faces[i];

			// gather the three vertex positions for this triangle
			std::array<glm::vec3, 3> triangle_vertices = {
				vertices[face[0]],
				vertices[face[1]],
				vertices[face[2]]};

			Triangle *new_triangle = nullptr; // initialize Triangle outside if/else

			// if vertex normals are provided --> smooth shading
			bool faceHasNormals = (i < normalIndices.size());
			if (faceHasNormals)
			{
				const auto &normal_index = normalIndices[i];
				std::array<glm::vec3, 3> triangle_vertex_normals = {
					normals[normal_index[0]],
					normals[normal_index[1]],
					normals[normal_index[2]]};
				new_triangle = new Triangle(triangle_vertices, triangle_vertex_normals, material);
			}

			// else --> flat shading
			else
			{
				new_triangle = new Triangle(triangle_vertices, material);
			}

			// apply the same transformation matrix to all triangles of this mesh
			new_triangle->setTransformation(transformation_matrix);

			// store pointer to triangle
			triangles.push_back(new_triangle);
		}
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

		glm::vec3 c = center - ray.origin;

		float cdotc = glm::dot(c, c);
		float cdotd = glm::dot(c, ray.direction);

		Hit hit;

		float D = 0;
		if (cdotc > cdotd * cdotd)
		{
			D = sqrt(cdotc - cdotd * cdotd);
		}
		if (D <= radius)
		{
			hit.hit = true;
			float t1 = cdotd - sqrt(radius * radius - D * D);
			float t2 = cdotd + sqrt(radius * radius - D * D);

			float t = t1;
			if (t < 0)
				t = t2;
			if (t < 0)
			{
				hit.hit = false;
				return hit;
			}

			hit.intersection = ray.origin + t * ray.direction;
			hit.normal = glm::normalize(hit.intersection - center);
			hit.distance = glm::distance(ray.origin, hit.intersection);
			hit.object = this;
		}
		else
		{
			hit.hit = false;
		}
		return hit;
	}
};

class Cone : public Object
{
private:
	Plane *plane;

public:
	Cone(Material material)
	{
		this->material = material;
		plane = new Plane(glm::vec3(0, 1, 0), glm::vec3(0.0, 1, 0));
	}
	Hit intersect(Ray ray)
	{

		Hit hit;
		hit.hit = false;

		glm::vec3 d = inverseTransformationMatrix * glm::vec4(ray.direction, 0.0); // implicit cast to vec3
		glm::vec3 o = inverseTransformationMatrix * glm::vec4(ray.origin, 1.0);	   // implicit cast to vec3
		d = glm::normalize(d);

		float a = d.x * d.x + d.z * d.z - d.y * d.y;
		float b = 2 * (d.x * o.x + d.z * o.z - d.y * o.y);
		float c = o.x * o.x + o.z * o.z - o.y * o.y;

		float delta = b * b - 4 * a * c;

		if (delta < 0)
		{
			return hit;
		}

		float t1 = (-b - sqrt(delta)) / (2 * a);
		float t2 = (-b + sqrt(delta)) / (2 * a);

		float t = t1;
		hit.intersection = o + t * d;
		if (t < 0 || hit.intersection.y > 1 || hit.intersection.y < 0)
		{
			t = t2;
			hit.intersection = o + t * d;
			if (t < 0 || hit.intersection.y > 1 || hit.intersection.y < 0)
			{
				return hit;
			}
		};

		hit.normal = glm::vec3(hit.intersection.x, -hit.intersection.y, hit.intersection.z);
		hit.normal = glm::normalize(hit.normal);

		Ray new_ray(o, d);
		Hit hit_plane = plane->intersect(new_ray);
		if (hit_plane.hit && hit_plane.distance < t && length(hit_plane.intersection - glm::vec3(0, 1, 0)) <= 1.0)
		{
			hit.intersection = hit_plane.intersection;
			hit.normal = hit_plane.normal;
		}

		hit.hit = true;
		hit.object = this;
		hit.intersection = transformationMatrix * glm::vec4(hit.intersection, 1.0); // implicit cast to vec3
		hit.normal = (normalMatrix * glm::vec4(hit.normal, 0.0));					// implicit cast to vec3
		hit.normal = glm::normalize(hit.normal);
		hit.distance = glm::length(hit.intersection - ray.origin);

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
// glm::vec3 ambient_light(0.1,0.1,0.1);
//  new ambient light
glm::vec3 ambient_light(0.001, 0.001, 0.001);
vector<Object *> objects; ///< A list of all objects in the scene

static const float EPSILON = 1e-4f; // offset to avoid self-intersection acne
static const int MAX_DEPTH = 1;		// recursion limit for reflections

/** Function for computing color of an object according to the Phong Model
 @param point A point belonging to the object for which the color is computer
 @param normal A normal vector the the point
 @param view_direction A normalized direction from the point to the viewer/camera
 @param material A material structure representing the material of the object
*/
glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal, glm::vec3 view_direction, Material material)
{
	glm::vec3 color(0.0f);

	// for each light source
	for (int light_num = 0; light_num < (int)lights.size(); ++light_num)
	{
		// Direction to light
		glm::vec3 Lvec = lights[light_num]->position - point; // vector from point to light
		float Ldist = glm::length(Lvec);					  // distance to light
		if (Ldist <= 0.0f)									  // if distance is zero or negative, skip this light
			continue;
		glm::vec3 L = Lvec / Ldist; // L is the normalized direction to light

		// Shadow test: cast a ray towards the light
		bool inShadow = false;
		{
			// Create shadow ray with origin slightly offset along normal to avoid self-intersection
			Ray shadowRay(point + normal * EPSILON, L);

			// for each object, check if it intersects the shadow ray before it reaches the light
			// if so, the point is in shadow
			for (int k = 0; k < (int)objects.size(); ++k)
			{
				Hit h = objects[k]->intersect(shadowRay);
				if (h.hit && h.distance > 0.0f && h.distance < Ldist - EPSILON)
				{
					inShadow = true;
					break;
				}
			}
		}
		// if the point is in shadow, skip this light’s contribution
		if (inShadow)
			continue;

		// Phong terms
		glm::vec3 lightDir = L;
		glm::vec3 reflectDir = glm::reflect(-lightDir, normal); // direction of perfect reflection

		float NdotL = glm::clamp(glm::dot(normal, lightDir), 0.0f, 1.0f);			// cosine of angle between normal and light direction
		float VdotR = glm::clamp(glm::dot(view_direction, reflectDir), 0.0f, 1.0f); // cosine of angle between view direction and reflection direction

		glm::vec3 diffuse = material.diffuse * NdotL;
		glm::vec3 specular = material.specular * powf(VdotR, material.shininess);

		// simple distance attenuation like you had (1/r^2)
		float r = glm::max(Ldist, 0.1f);
		color += lights[light_num]->color * (diffuse + specular) / (r * r);
	}

	// finally we add ambient term
	color += ambient_light * material.ambient;
	return glm::clamp(color, glm::vec3(0.0f), glm::vec3(1.0f)); // we return the color clamped between 0 and 1
}

/**
 Functions that computes a color along the ray
 @param ray Ray that should be traced through the scene
 @return Color at the intersection point
 */
glm::vec3 trace_ray(const Ray &ray, int depth = 0)
{
	Hit closest_hit;
	closest_hit.hit = false;
	closest_hit.distance = INFINITY;

	// Find nearest hit
	for (int k = 0; k < (int)objects.size(); ++k)
	{
		Hit hit = objects[k]->intersect(ray);
		if (hit.hit && hit.distance < closest_hit.distance)
		{
			closest_hit = hit;
		}
	}

	// Miss therefore in the background
	if (!closest_hit.hit)
	{
		return glm::vec3(0.0f);
	}

	// Local shading
	const Material mat = closest_hit.object->getMaterial();
	glm::vec3 viewDir = glm::normalize(-ray.direction);
	glm::vec3 local = PhongModel(
		closest_hit.intersection,
		closest_hit.normal,
		viewDir,
		mat);

	// Stop if max depth or non-reflective
	float kr = glm::clamp(mat.reflectivity, 0.0f, 1.0f);
	if (kr <= 0.0f || depth >= MAX_DEPTH)
	{
		return local;
	}

	// Perfect mirror reflection
	glm::vec3 Rdir = glm::normalize(glm::reflect(ray.direction, closest_hit.normal));
	// Create reflected ray with origin slightly offset along normal to avoid self-intersection
	Ray reflRay(closest_hit.intersection + closest_hit.normal * EPSILON, Rdir);
	// we trace the reflected ray recursively
	glm::vec3 reflCol = trace_ray(reflRay, depth + 1);

	// Mix reflection with local color because of reflectivity
	return (1.0f - kr) * local + kr * reflCol;
}
/**
 Function defining the scene
 */
void sceneDefinition()
{

	Material green_diffuse;
	green_diffuse.ambient = glm::vec3(0.03f, 0.1f, 0.03f);
	green_diffuse.diffuse = glm::vec3(0.3f, 1.0f, 0.3f);

	Material red_specular;
	red_specular.diffuse = glm::vec3(1.0f, 0.2f, 0.2f);
	red_specular.ambient = glm::vec3(0.01f, 0.02f, 0.02f);
	red_specular.specular = glm::vec3(0.5);
	red_specular.shininess = 10.0;

	Material blue_specular;
	blue_specular.ambient = glm::vec3(0.7f, 0.7f, 1.0f);
	blue_specular.diffuse = glm::vec3(0.7f, 0.7f, 1.0f);
	blue_specular.specular = glm::vec3(0.6);
	blue_specular.shininess = 50.0;
	blue_specular.reflectivity = 0.75f;

	// spheres
	objects.push_back(new Sphere(1.0, glm::vec3(1, -2, 8), blue_specular));
	objects.push_back(new Sphere(0.5, glm::vec3(-1, -2.5, 6), red_specular));

	// meshes
	// auto *armadillo = new Mesh("meshes/armadillo_with_normals.obj", green_diffuse);

	// glm::mat4 S = glm::scale(glm::vec3(1.5f));
	// glm::mat4 R = glm::mat4(1.0f);
	// glm::mat4 T = glm::translate(glm::vec3(-5.0f, -3.0f, 10.0f));
	// armadillo->setTransformation(T * R * S);

	// objects.push_back(armadillo);

	// auto *bunny = new Mesh("meshes/bunny_with_normals.obj", green_diffuse);

	// S = glm::scale(glm::vec3(1.5f));
	// R = glm::mat4(1.0f);
	// T = glm::translate(glm::vec3(0.0f, -3.0f, 10.0f));
	// bunny->setTransformation(T * R * S);

	// objects.push_back(bunny);

	// auto *lucy = new Mesh("meshes/lucy_with_normals.obj", green_diffuse);

	// S = glm::scale(glm::vec3(1.5f));
	// R = glm::mat4(1.0f);
	// T = glm::translate(glm::vec3(5.0f, -3.0f, 10.0f));
	// lucy->setTransformation(T * R * S);

	// objects.push_back(lucy);

	// lights
	lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(1.0, 1.0, 1.0)));
	lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.1)));
	lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.4)));

	Material red_diffuse;
	red_diffuse.ambient = glm::vec3(0.09f, 0.06f, 0.06f);
	red_diffuse.diffuse = glm::vec3(0.9f, 0.6f, 0.6f);

	Material blue_diffuse;
	blue_diffuse.ambient = glm::vec3(0.06f, 0.06f, 0.09f);
	blue_diffuse.diffuse = glm::vec3(0.6f, 0.6f, 0.9f);
	objects.push_back(new Plane(glm::vec3(0, -3, 0), glm::vec3(0.0, 1, 0)));
	objects.push_back(new Plane(glm::vec3(0, 1, 30), glm::vec3(0.0, 0.0, -1.0), green_diffuse));
	objects.push_back(new Plane(glm::vec3(-15, 1, 0), glm::vec3(1.0, 0.0, 0.0), red_diffuse));
	objects.push_back(new Plane(glm::vec3(15, 1, 0), glm::vec3(-1.0, 0.0, 0.0), blue_diffuse));
	objects.push_back(new Plane(glm::vec3(0, 27, 0), glm::vec3(0.0, -1, 0)));
	objects.push_back(new Plane(glm::vec3(0, 1, -0.01), glm::vec3(0.0, 0.0, 1.0), green_diffuse));

	// Cones
	Material yellow_specular;
	yellow_specular.ambient = glm::vec3(0.1f, 0.10f, 0.0f);
	yellow_specular.diffuse = glm::vec3(0.4f, 0.4f, 0.0f);
	yellow_specular.specular = glm::vec3(1.0);
	yellow_specular.shininess = 0.0;
	yellow_specular.reflectivity = 0.9f;

	Cone *cone = new Cone(yellow_specular);
	glm::mat4 translationMatrix = glm::translate(glm::vec3(5, 9, 14));
	glm::mat4 scalingMatrix = glm::scale(glm::vec3(3.0f, 12.0f, 3.0f));
	glm::mat4 rotationMatrix = glm::rotate(glm::radians(180.0f), glm::vec3(1, 0, 0));
	cone->setTransformation(translationMatrix * scalingMatrix * rotationMatrix);
	objects.push_back(cone);

	Cone *cone2 = new Cone(green_diffuse);
	translationMatrix = glm::translate(glm::vec3(6, -3, 7));
	scalingMatrix = glm::scale(glm::vec3(1.0f, 3.0f, 1.0f));
	rotationMatrix = glm::rotate(glm::atan(3.0f), glm::vec3(0, 0, 1));
	cone2->setTransformation(translationMatrix * rotationMatrix * scalingMatrix);
	objects.push_back(cone2);
}
glm::vec3 toneMapping(glm::vec3 intensity)
{
	float gamma = 1.0 / 2.0;
	float alpha = 12.0f;
	return glm::clamp(alpha * glm::pow(intensity, glm::vec3(gamma)), glm::vec3(0.0), glm::vec3(1.0));
}
int main(int argc, const char *argv[])
{

	clock_t t = clock(); // variable for keeping the time of the rendering

	int width = 1024; // width of the image, actual width is 1024 pixels, we put 200 for testing
	int height = 768; // height of the image, actual height is 768 pixels, we put 150 for testing
	float fov = 90;	  // field of view

	sceneDefinition(); // Let's define a scene

	Image image(width, height); // Create an image where we will store the result
	vector<glm::vec3> image_values(width * height);

	float s = 2 * tan(0.5 * fov / 180 * M_PI) / width;
	float X = -s * width / 2;
	float Y = s * height / 2;

	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
		{

			float dx = X + i * s + s / 2;
			float dy = Y - j * s - s / 2;
			float dz = 1;

			glm::vec3 origin(0, 0, 0);
			glm::vec3 direction(dx, dy, dz);
			direction = glm::normalize(direction);

			Ray ray(origin, direction);
			image.setPixel(i, j, toneMapping(trace_ray(ray)));
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
// Hey Rafa here is the documentation for the added features:
// ================================================================
// Secondary Effects Added (Shadows + Reflections)
// ================================================================
//
// Overall:
// --------
// Added support for two core ray-tracing effects:
//   1) Shadows -> points only receive light if visible to the light
//   2) Reflections -> reflective materials spawn secondary reflection rays
//
// Introduced two global constants:
//   EPSILON   -> small offset to avoid self-intersection (“acne”)
//   MAX_DEPTH -> recursion limit for reflections
//
// Material was extended with:
//   float reflectivity;  // 0 = no reflection, 1 = perfect mirror
//
// =================================================================
// PhongModel() — Shadow Rays
// =================================================================
//
// For every light, we now:
//   - Compute a shadow ray from the surface point toward the light
//   - Offset the origin by normal * EPSILON to avoid self-shadowing
//   - Check whether any object intersects that shadow ray before the light
//     • If an object blocks the light → skip diffuse+specular for that light
//   - Otherwise compute:
//       diffuse = kd * max(dot(N, L), 0)
//       specular = ks * pow(dot(V, R), shininess)
//   - Still apply ambient term as before (ambient light is not shadowed)
//
// Result: Proper hard shadows cast by all opaque objects.
//
// =================================================================
// trace_ray() — Recursive Perfect Reflections
// =================================================================
//
// trace_ray now takes an additional parameter: recursion depth.
// The logic is:
//   1) Find closest object hit (unchanged)
//   2) Compute local lighting using PhongModel (unchanged)
//   3) If material.reflectivity == 0, return local shading
//   4) Otherwise:
//        - Compute reflection direction using glm::reflect()
//        - Spawn a new ray from (point + normal * EPSILON)
//        - Recursively call trace_ray(reflRay, depth + 1)
//        - Mix local+reflection via:
//            color = (1 - kr) * local + kr * reflection
//   5) Prevent infinite recursion via MAX_DEPTH
//
// Result: Objects with reflectivity > 0 act as mirrors.
// ================================================================
