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
#include <barrier>
#include "omp.h"
#include <arm_neon.h>
#include <cassert>
#include <iomanip>

#define eps2 1e-3
#define MAXDEPTH 13
#define MAXBODY_LEAF 8
#define N 1000
#define NTHREADS 10


#define real double




namespace Simulator {

const real G = 5e-3; // constante de Gauss (parsecs, Myrs)


class Array{
	std::array<real, 3> vec3;
	
public:
	~Array() = default;
	Array() {vec3 = {0.0, 0.0, 0.0};}
	Array(real x, real y, real z) : vec3{x, y, z}{}
	
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
	Array force;
	real mass;
	std::size_t id;
	void print() const;
};

static inline void gravity(const Body& s1, const Body& s2, Array * force){
	Array r = s1.positions - s2.positions;
	real r2 = r.norm2();
	real inv_r = 1/std::sqrt(r2+eps2);
	real inv_r2 = inv_r * inv_r;
	real inv_r3 = inv_r2 * inv_r;
	*force = (-Simulator::G * s1.mass * s2.mass * inv_r3) * r;
}
void compute_csts();
void init_mass(size_t nBodies, real &m);
void init_position(Simulator::Array& p0);
void init_velocity(real x, real y, real z, Simulator::Array& v0);

void compute_force(std::vector<Body>& bodies);
void compute_force_parallel(std::vector<Body>& bodies);
void compute_energy(const std::vector<Body>& bodies, std::vector<real> * energy);
void compute_energy_parallel(real * x, real * y, real * z, real * vx, real * vy, real * vz, real * m, size_t n, std::vector<real> * energy);

struct Integrator{
	virtual void step(std::vector<Body>& bodies, const real dt) = 0;
	virtual void step_parallel(std::vector<Body>& bodies, const real dt) = 0;
	virtual void step_parallel(real* x, real* y, real* z, real* vx, real* vy, real* vz, real* fx, real* fy, real* fz, real* m, const real dt, size_t n) = 0;
	virtual ~Integrator() = default;
};

struct Euler : Integrator{
	void step(std::vector<Body>& bodies, const real dt) override;
	void step_parallel(std::vector<Body>& bodies, const real dt) override;
	void step_parallel(real* x, real* y, real* z, real* vx, real* vy, real* vz, real* fx, real* fy, real* fz, real* m, const real dt, size_t n) override;
};

struct Leapfrog : Integrator{
	void step(std::vector<Body>& bodies, const real dt) override;
	void step_parallel(std::vector<Body>& bodies, const real dt) override;
	void step_parallel(real* x, real* y, real* z, real* vx, real* vy, real* vz, real* fx, real* fy, real* fz, real* m, const real dt, size_t n) override;
};


//##################################################################

struct OctreeNode {
	Array center;          // centre du cube
	Array com;			   // centre de masse
	real half_size;      // demi-taille
	real mass;           // masse totale

	struct OctreeNode* children[8];
	Body* bodies[MAXBODY_LEAF];            // NULL ou 1 particule
	
	int body_count;
};


struct Barnes_Hut : Integrator{
	void step(std::vector<Body>& bodies, const real dt) override;
	void step_parallel(std::vector<Body>& bodies, const real dt) override;
	void step_parallel(real* x, real* y, real* z, real* vx, real* vy, real* vz, real* fx, real* fy, real* fz, real* m, const real dt, size_t n) override;
};

int get_octant(OctreeNode* node, Body* b);

OctreeNode* create_child(OctreeNode* parent, int oct);

OctreeNode* create_node(const Array& c_in, real half_size);

void insert_body(OctreeNode* node, Body* b);

void compute_mass(OctreeNode* node);

void compute_force_BH(OctreeNode* node, Body* b, real theta);

OctreeNode* build_octree(std::vector<Body>& bodies);

void free_tree(OctreeNode* node);




struct OctreeNode_t {
	Array center;          // centre du cube
	Array com;			   // centre de masse
	real half_size;      // demi-taille
	real mass;           // masse totale

	struct OctreeNode_t* children[8];
	size_t index[MAXBODY_LEAF];
	
	size_t body_count;
};

int get_octant_t(OctreeNode_t* node, real x, real y, real z);

OctreeNode_t* create_child_t(OctreeNode_t* parent, int oct);

OctreeNode_t* create_node_t(const Array& c_in, real half_size);

void insert_body_t(OctreeNode_t* node, size_t id_body, real *x, real *y, real *z);

void compute_mass_t(OctreeNode_t* node, real * x, real * y, real * z, real * m);

void compute_force_BH_t(OctreeNode_t* node, real * x, real * y, real * z, real * fx, real * fy, real * fz, real * m, size_t id_body, real theta);

OctreeNode_t* build_octree_t(real * x, real * y, real * z, real * m, size_t n);

void free_tree_t(OctreeNode_t* node);

//##################################################################

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

