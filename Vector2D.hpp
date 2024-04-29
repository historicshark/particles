#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <numeric>
#include <random>

class Vector
{
    std::array<double, 2> vector;
    
    public:
    Vector(double x, double y)
    : vector{x, y}
    {};
    
    Vector()
    : vector{0,0}
    {};
    
    Vector& operator=(const Vector& other)
    {
        if(this == &other)
        {
            return *this;
        }
        else
        {
            vector = other.vector;
            return *this;
        }
    }
    
    Vector operator+(Vector other) { return {vector[0] + other.vector[0], vector[1] + other.vector[1]}; };
    Vector operator-(Vector other) { return {vector[0] - other.vector[0], vector[1] - other.vector[1]}; };
    Vector operator/(Vector other) { return {vector[0] / other.vector[0], vector[1] / other.vector[1]}; };
    Vector operator*(Vector other) { return {vector[0] * other.vector[0], vector[1] * other.vector[1]}; };
    
    Vector operator+=(Vector other) { return {vector[0] += other.vector[0], vector[1] += other.vector[1]}; };
    
    friend Vector operator+(double c, Vector v) { return {v.vector[0] + c, v.vector[1] + c}; };
    friend Vector operator+(Vector v, double c) { return {v.vector[0] + c, v.vector[1] + c}; };
    friend Vector operator-(double c, Vector v) { return {c - v.vector[0], c - v.vector[1]}; };
    friend Vector operator-(Vector v, double c) { return {v.vector[0] - c, v.vector[1] - c}; };
    friend Vector operator/(double c, Vector v) { return {c / v.vector[0], c / v.vector[1]}; };
    friend Vector operator/(Vector v, double c) { return {v.vector[0] / c, v.vector[1] / c}; };
    friend Vector operator*(double c, Vector v) { return {v.vector[0] * c, v.vector[1] * c}; };
    friend Vector operator*(Vector v, double c) { return {v.vector[0] * c, v.vector[1] * c}; };
    
    friend Vector operator^(Vector v, double c) { return {std::pow(v.vector[0], c), std::pow(v.vector[1], c)}; };
    
    double operator[](size_t i) { return vector[i]; };
    
    auto norm_squared() { return std::accumulate(vector.begin(), vector.end(), 0.0, [](auto a, auto x){return std::move(a) + pow(x,2);}); };
    
    auto norm() { return std::sqrt(norm_squared()); };
    
    double dot(Vector& other) { return vector[0] * other.vector[0] + vector[1] * other.vector[1]; };
    
    auto string()
    { 
        std::ostringstream os;
        os << vector[0] << "," << vector[1];
        return os.str();
    };
    
    void print() { std::cout << "(" << string() << ")" << std::endl; }
};