//
//  Material.h
//  Raytracer
//
//  Created by Piotr Didyk on 14.07.21.
//

#ifndef Material_h
#define Material_h

#include "glm/glm.hpp"

/**
 Structure describing a material of an object
 */
struct Material{
    glm::vec3 ambient = glm::vec3(0.0);
    glm::vec3 diffuse = glm::vec3(1.0);
    glm::vec3 specular = glm::vec3(0.0);
    float shininess = 0.0;
};

Material make_box_material(glm::vec3 color) {
    Material wall_material;
	wall_material.ambient = color;
	wall_material.diffuse = color;
	wall_material.specular = glm::vec3(0.2f);
	wall_material.shininess = 10.0;
    return wall_material;
}

#endif /* Material_h */
