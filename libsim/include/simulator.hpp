//
//  simulator.hpp
//  simulator_galaxie
//
//  Created by Nolwen Dolléans on 08/01/2026.
//

#ifndef simulator_hpp
#define simulator_hpp
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <memory>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <chrono>
#include <thread>
#include <numeric>
#include <mutex>

#define eps3 1e-6
#define N 100
#define NTHREADS 10


#define real double

namespace Simulator {

const double G = 5e-3; // constante de Gauss (UA, jours)
const double prsc = 3e13;       // km

class Array{
	
public:
	std::array<real, 3> vec3;
	~Array() = default;
	Array() {vec3 = {0.0, 0.0, 0.0};}
	Array(real x, real y, real z) : vec3{x, y, z}{}
	
	void print(std::ofstream& file) const;
	
	real& operator[](const std::size_t i);
	const real& operator[](const std::size_t i) const;
	
	Array(const Array&) = default;
	Array& operator=(const Array&) = default;
	
	Array operator+(const Array& other) const;
	Array operator+=(const Array& other);
	Array operator-(const Array& other) const;
	Array operator-=(const Array& other);
	Array operator*(const real alpha) const;
	Array operator*=(const real alpha);
	
	real norm2() const;
	
	auto begin() const;
	auto begin();
	auto end() const;
	auto end();
};
Array operator*(const real alpha, const Array& arr);

Array cross(const Array& a, const Array& b);

real dot(const Array& a, const Array& b);


struct Body{
	Array positions;
	Array velocities;
	real mass;
	Array force;
	std::size_t id;
	void print() const;
};


Array gravity(const Body& s1, const Body& s2);
void compute_force(std::vector<Body>& bodies);
void compute_force_parallel(std::vector<Body>& bodies);
std::vector<real> compute_energy(const std::vector<Body>& bodies);

struct Integrator{
	virtual void step(std::vector<Body>& bodies, const real dt) = 0;
	virtual void step_parallel(std::vector<Body>& bodies, const real dt) = 0;
	virtual ~Integrator() = default;
};

struct Euler : Integrator{
	void step(std::vector<Body>& bodies, const real dt) override;
	void step_parallel(std::vector<Body>& bodies, const real dt) override;
};

struct Leapfrog : Integrator{
	void step(std::vector<Body>& bodies, const real dt) override;
	void step_parallel(std::vector<Body>& bodies, const real dt) override;
};

struct Barnes_Hut : Integrator{
	void step(std::vector<Body>& bodies, const real dt) override;
	void step_parallel(std::vector<Body>& bodies, const real dt) override;
};


enum class IntChoice{Euler, Leapfrog, Barnes_Hut};

class galaxy{
	std::unique_ptr<Integrator> integrator = nullptr;
	std::size_t Nbodies;
	IntChoice intchoice;
public:
	std::vector<Body> bodies;
	galaxy() = default;
	
	void initialize(const std::size_t i, const IntChoice& choice);
	
	void compute(const std::string& path);
	void compute_parallel(const std::string& path);
	~galaxy() {
		Integrator * raw = integrator.release();
		delete raw;
	}
};

class Chrono{
	std::chrono::time_point<std::chrono::system_clock> time{};
	std::chrono::milliseconds rt{};
	bool running = false;
public:
	Chrono() = default;
	~Chrono() = default;
	/**
	 * @brief Store the actual time and start the chrono
	 */
	void start();
	/**
	 * @brief Store the runtime the chrono
	 */
	void stop();
	/**
	 * @brief Retrun the runtime computed
	 * @return runtime stored
	 */
	std::chrono::milliseconds runtime() const;
	
	void print() const;
};

}/* end simulator */

#endif /* simulator_hpp */

