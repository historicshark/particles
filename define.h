#pragma once

#include "Vector2D.h"

// #define DIM 2
//

// // #define DRAG
// // #define BUOYANCY
// #define WALL

// #ifdef DRAG
// #   define DRAG_(...) __VA_ARGS__
// #else
// #   define DRAG_(...)
// #endif

// #ifdef BUOYANCY
// #   define BUOYANCY_(...) __VA_ARGS__
// #else
// #   define BUOYANCY_(...)
// #endif

// #ifdef WALL
// #   define WALL_(...) __VA_ARGS__
// #else
// #   define WALL_(...)
// #endif

// inline auto sum(auto... ts)
// {
    // return (... + ts);
// }

// inline auto sum() {return 0.0;}

double dot(Vector v1, Vector v2)
{
    auto v = v1 * v2;
    return v[0] + v[1];
}
