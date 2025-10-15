// Colors.h
// 
// Predefined color names and their corresponding RGB values
// Declaration of generative AI: used ChatGPT to propagate the list


#ifndef COLORS_H
#define COLORS_H

#include "glm/glm.hpp"

namespace Colors {

    // ===== Basic Colors =====
    const glm::vec3 White     = glm::vec3(1.0f, 1.0f, 1.0f);
    const glm::vec3 Black     = glm::vec3(0.0f, 0.0f, 0.0f);
    const glm::vec3 Gray      = glm::vec3(0.5f, 0.5f, 0.5f);
    const glm::vec3 LightGray = glm::vec3(0.8f, 0.8f, 0.8f);
    const glm::vec3 DarkGray  = glm::vec3(0.3f, 0.3f, 0.3f);

    // ===== Reds =====
    const glm::vec3 Red       = glm::vec3(1.0f, 0.0f, 0.0f);
    const glm::vec3 DarkRed   = glm::vec3(0.5f, 0.0f, 0.0f);
    const glm::vec3 LightCoral= glm::vec3(0.94f, 0.5f, 0.5f);
    const glm::vec3 Salmon    = glm::vec3(0.98f, 0.5f, 0.45f);
    const glm::vec3 Crimson   = glm::vec3(0.86f, 0.08f, 0.24f);

    // ===== Greens =====
    const glm::vec3 Green     = glm::vec3(0.0f, 1.0f, 0.0f);
    const glm::vec3 DarkGreen = glm::vec3(0.0f, 0.5f, 0.0f);
    const glm::vec3 LightGreen= glm::vec3(0.5f, 1.0f, 0.5f);
    const glm::vec3 Lime      = glm::vec3(0.7f, 1.0f, 0.0f);
    const glm::vec3 Olive     = glm::vec3(0.5f, 0.5f, 0.0f);
    const glm::vec3 SeaGreen  = glm::vec3(0.18f, 0.55f, 0.34f);

    // ===== Blues =====
    const glm::vec3 Blue      = glm::vec3(0.0f, 0.0f, 1.0f);
    const glm::vec3 DarkBlue  = glm::vec3(0.0f, 0.0f, 0.5f);
    const glm::vec3 LightBlue = glm::vec3(0.4f, 0.7f, 1.0f);
    const glm::vec3 SkyBlue   = glm::vec3(0.53f, 0.81f, 0.92f);
    const glm::vec3 Teal      = glm::vec3(0.0f, 0.5f, 0.5f);
    const glm::vec3 Cyan      = glm::vec3(0.0f, 1.0f, 1.0f);

    // ===== Purples / Magentas =====
    const glm::vec3 Purple    = glm::vec3(0.5f, 0.0f, 0.5f);
    const glm::vec3 Violet    = glm::vec3(0.93f, 0.51f, 0.93f);
    const glm::vec3 Magenta   = glm::vec3(1.0f, 0.0f, 1.0f);
    const glm::vec3 Orchid    = glm::vec3(0.85f, 0.44f, 0.84f);
    const glm::vec3 Indigo    = glm::vec3(0.29f, 0.0f, 0.51f);

    // ===== Yellows / Oranges =====
    const glm::vec3 Yellow    = glm::vec3(1.0f, 1.0f, 0.0f);
    const glm::vec3 Gold      = glm::vec3(1.0f, 0.84f, 0.0f);
    const glm::vec3 Orange    = glm::vec3(1.0f, 0.5f, 0.0f);
    const glm::vec3 DarkOrange= glm::vec3(0.8f, 0.4f, 0.0f);
    const glm::vec3 Khaki     = glm::vec3(0.76f, 0.69f, 0.57f);
    const glm::vec3 Brown       = glm::vec3(0.65f, 0.16f, 0.16f);

    // ===== Pastels & Neutrals =====
    const glm::vec3 Beige     = glm::vec3(0.96f, 0.96f, 0.86f);
    const glm::vec3 Tan       = glm::vec3(0.82f, 0.71f, 0.55f);
    const glm::vec3 LightPink = glm::vec3(1.0f, 0.71f, 0.76f);
    const glm::vec3 Peach     = glm::vec3(1.0f, 0.85f, 0.73f);
    const glm::vec3 Mint      = glm::vec3(0.6f, 1.0f, 0.8f);

} // namespace Colors

#endif // COLORS_H