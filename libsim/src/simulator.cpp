//
//  simulator.cpp
//  simulator_galaxie
//
//  Created by Nolwen Dolléans on 08/01/2026.
//

#include "simulator.hpp"

struct args_t{
	std::vector<std::size_t> bodies_id;
	std::vector<Simulator::Body> * bodies;
	
};

double invSqrt(double number) {
	union {
		double f;
		uint64_t i;
	} conv;

	double x2 = number * 0.5;
	const double threehalfs = 1.5;

	conv.f = number;
	conv.i = 0x5fe6ec85e7de30daULL - (conv.i >> 1); // constante magique pour double
	conv.f = conv.f * (threehalfs - x2 * conv.f * conv.f); // 1ère itération
	// conv.f = conv.f * (threehalfs - x2 * conv.f * conv.f); // optionnel 2ème itération pour plus de précision

	return conv.f;
}

/* ---------- ARRAY ---------*/
void Simulator::Array::print(std::ofstream& file) const{
	/*
	if (!file) {
		std::cerr << "Error: File doesn't open.\n";
		exit(EXIT_FAILURE);
	}
	std::for_each(vec3.begin(), vec3.end(), [&file](auto& value){file << (real)value <<",";});
	file << std::endl;*/
}

real& Simulator::Array::operator[](const std::size_t i){
	return vec3[i];
}

const real& Simulator::Array::operator[](const std::size_t i) const{
	return vec3[i];
}

Simulator::Array Simulator::Array::operator+(const Array& other) const{
	return {
		vec3[0] + other[0],
		vec3[1] + other[1],
		vec3[2] + other[2]
	};
}

Simulator::Array Simulator::Array::operator+=(const Array& other){
	vec3[0] += other[0];
	vec3[1] += other[1];
	vec3[2] += other[2];
	return *this;
}

Simulator::Array Simulator::Array::operator-(const Array& other) const{
	return {
		vec3[0] - other[0],
		vec3[1] - other[1],
		vec3[2] - other[2]
	};
}

Simulator::Array Simulator::Array::operator-=(const Array& other){
	vec3[0] -= other[0];
	vec3[1] -= other[1];
	vec3[2] -= other[2];
	return *this;
}

Simulator::Array Simulator::Array::operator*(const real alpha) const{
	return {
		vec3[0] * alpha,
		vec3[1] * alpha,
		vec3[2] * alpha
	};
}

Simulator::Array Simulator::Array::operator*=(const real alpha){
	vec3[0] *= alpha;
	vec3[1] *= alpha;
	vec3[2] *= alpha;
	return *this;
}

Simulator::Array Simulator::operator*(const real alpha, const Array& other){
	return other * alpha;
}

real Simulator::Array::norm2() const{
	return dot(*this, *this);
}

auto Simulator::Array::begin() const{
	return vec3.begin();
}

auto Simulator::Array::begin(){
	return vec3.begin();
}

auto Simulator::Array::end() const{
	return vec3.end();
}

auto Simulator::Array::end(){
	return vec3.end();
}

Simulator::Array Simulator::cross(const Array& a, const Array& b){
	return {
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0]
	};
}

real Simulator::dot(const Array& a, const Array& b){
	real res = 0.0;
	for (int i = 0; i<3; ++i) res += a[i]*b[i];
	return res;
}


/* ------------- SIMULATOR -------------*/
void Simulator::Body::print() const{
	std::cout << "Position: ";
}

real random_real(real min, real max) {
	return min + (max - min) * (rand() / (real)RAND_MAX);
}

Simulator::Array vitesseInitiale3D(double x, double y, double v_max = 220.0, double r0 = 3.0) {
	// Calcul de la distance au centre dans le plan XY
	double r = std::sqrt(x*x + y*y);
	double inv_r = invSqrt(x*x + y*y);
	
	// Vitesse orbitale selon le profil de rotation
	double v = v_max * r / (r + r0);
	
	// Dispersion aléatoire ±10%
	double dispersion = 0.1 * v;
	double v_finale = v + ((rand() / double(RAND_MAX)) * 2 - 1) * dispersion;

	// Direction tangente à l'orbite dans le plan XY
	double vx = -v_finale * (y * inv_r);
	double vy =  v_finale * (x * inv_r);

	// Vitesse verticale aléatoire ±10 km/s
	double vz = ((rand() / double(RAND_MAX)) * 2 - 1) * 10.0;

	return {vx, vy, vz};
}

void Simulator::galaxy::initialize(const std::size_t i, const IntChoice& choice){
	Nbodies = i;
	bodies = std::vector<Body>(i);
	srand((unsigned)time(0));
	real R_max = 1.66e7;
	
	bodies[0].mass = 1e6;
	bodies[0].positions = {0.0, 0.0, 0.0};
	bodies[0].velocities = {0.0, 0.0, 0.0};
	bodies[0].id = 0;
	
	for (std::size_t j = 1; j < i; ++j) {
		real r = R_max * sqrt(random_real(0.0, 1.0));
		real theta = random_real(0.0, 2*M_PI);
		real z = random_real(-0.33e6, 0.33e6);

		bodies[j].mass = random_real(0.1f, 1.0f);
		bodies[j].positions = {r * cos(theta), r * sin(theta), z};
		bodies[j].velocities = vitesseInitiale3D(r * cos(theta), r * sin(theta));
		bodies[j].id = j;

	}
	switch (choice) {
		case IntChoice::Euler:
			integrator = std::make_unique<Euler>();
			break;
		case IntChoice::Leapfrog:
			integrator = std::make_unique<Leapfrog>();
			break;
		case IntChoice::Barnes_Hut:
			integrator = std::make_unique<Barnes_Hut>();
			break;
	}
}

void Simulator::galaxy::compute(const std::string& path){
	std::cout << "--- Simulator started ---" << std::endl;
	Chrono timer;
	timer.start();
	std::ofstream file_(path+"stars_data.csv", std::ios::trunc);
	file_.close();
	std::ofstream file2(path+"stars_data.csv");
	file2 << "temps,id,x,y,z,vx,vy,vz,mass;\n";
	file2.close();
	
	std::ofstream file3(path+"energy.csv", std::ios::trunc);
	file3.close();
	std::ofstream file4(path+"energy.csv");
	file4 << "temps,T,V,E;\n";
	file4.close();
	
	real t = 0;
	const real end = 10.0f;
	const real dt = end/(real)N;
	std::ofstream file(path+"stars_data.csv", std::ios::app);
	std::ofstream energy_file(path+"energy.csv", std::ios::app);
	
	while (t<end) {
		for (std::size_t i = 0; i<Nbodies; ++i) {
			   file << t << ","
					<< i << ","
					<< bodies[i].positions[0] 	<< ","
					<< bodies[i].positions[1] 	<< ","
					<< bodies[i].positions[2] 	<< ","
					<< bodies[i].velocities[0]	<< ","
					<< bodies[i].velocities[1]	<< ","
					<< bodies[i].velocities[2] 	<< ","
					<< bodies[i].mass 			<< ";\n";
		}
		
		integrator->step(bodies, dt);
		
		std::vector<real> energy = compute_energy(bodies);
		energy_file << t << ","
					<< energy[0] 	<< ","
					<< energy[1] 	<< ","
					<< energy[2] 	<< ";\n";
		t += dt;
	}
	file.close();
	energy_file.close();
	timer.stop();
	timer.print();
}

void Simulator::galaxy::compute_parallel(const std::string& path){
	std::cout << "--- Simulator started ---" << std::endl;
	Chrono timer;
	timer.start();
	std::ofstream file_(path+"stars_data.csv", std::ios::trunc);
	file_.close();
	std::ofstream file2(path+"stars_data.csv");
	file2 << "temps,id,x,y,z,vx,vy,vz,mass;\n";
	file2.close();
	
	std::ofstream file3(path+"energy.csv", std::ios::trunc);
	file3.close();
	std::ofstream file4(path+"energy.csv");
	file4 << "temps,T,V,E;\n";
	file4.close();
	
	real t = 0;
	const real end = 10.0f;
	const real dt = end/(real)N;
	std::ofstream file(path+"stars_data.csv", std::ios::app);
	std::ofstream energy_file(path+"energy.csv", std::ios::app);
	
	while (t<end) {
		for (std::size_t i = 0; i<Nbodies; ++i) {
			   file << t << ","
					<< i << ","
					<< bodies[i].positions[0] 	<< ","
					<< bodies[i].positions[1] 	<< ","
					<< bodies[i].positions[2] 	<< ","
					<< bodies[i].velocities[0]	<< ","
					<< bodies[i].velocities[1]	<< ","
					<< bodies[i].velocities[2] 	<< ","
					<< bodies[i].mass 			<< ";\n";
		}
		
		integrator->step_parallel(bodies, dt);
		
		std::vector<real> energy = compute_energy(bodies);
		energy_file << t << ","
					<< energy[0] 	<< ","
					<< energy[1] 	<< ","
					<< energy[2] 	<< ";\n";
		t += dt;
	}
	file.close();
	energy_file.close();
	timer.stop();
	timer.print();
}

void Simulator::Euler::step(std::vector<Body>& bodies, const real dt){
	std::size_t n = bodies.size();
	compute_force(bodies);
	for (std::size_t j = 0; j<n; ++j) {
		// positions
		bodies[j].positions += bodies[j].velocities * dt;
		// velocities
		bodies[j].velocities +=1.0/bodies[j].mass * bodies[j].force * dt;
	}
}

void Simulator::Euler::step_parallel(std::vector<Body>& bodies, const real dt){
	std::size_t n = bodies.size();
	compute_force(bodies);
	for (std::size_t j = 0; j<n; ++j) {
		// positions
		bodies[j].positions += bodies[j].velocities * dt;
		// velocities
		bodies[j].velocities +=1.0/bodies[j].mass * bodies[j].force * dt;
	}
}

void Simulator::Leapfrog::step(std::vector<Body>& bodies, const real dt){
	std::size_t n = bodies.size();
	compute_force(bodies);
	for (std::size_t j = 0; j<n; ++j) {
		// velocities
		bodies[j].velocities +=1.0/bodies[j].mass * bodies[j].force * (dt/2.0);
		// positions
		bodies[j].positions += bodies[j].velocities * dt;
		// velocities
	}
	compute_force(bodies);
	for (std::size_t j = 0; j<n; ++j) {
		bodies[j].velocities +=1.0/bodies[j].mass * bodies[j].force * (dt/2.0);
	}
}

void Simulator::Leapfrog::step_parallel(std::vector<Body>& bodies, const real dt){
	std::size_t n = bodies.size();
	compute_force_parallel(bodies);
	for (std::size_t j = 0; j<n; ++j) {
		// velocities
		bodies[j].velocities +=1.0/bodies[j].mass * bodies[j].force * (dt/2.0);
		// positions
		bodies[j].positions += bodies[j].velocities * dt;
		// velocities
	}
	compute_force_parallel(bodies);
	for (std::size_t j = 0; j<n; ++j) {
		bodies[j].velocities +=1.0/bodies[j].mass * bodies[j].force * (dt/2.0);
	}
}

void Simulator::Barnes_Hut::step(std::vector<Body>& bodies, const real dt){
	return;
}

void Simulator::Barnes_Hut::step_parallel(std::vector<Body>& bodies, const real dt){
	return;
}

Simulator::Array Simulator::gravity(const Body& s1, const Body& s2){
	Array r = s1.positions - s2.positions;  // vecteur séparant s1 et s2
	real r2 = r.norm2();                     // distance au carré
	real soft2 = 1e-3;                       // ajuster selon l’échelle
	real inv_r = invSqrt(r2+soft2);        // 1/sqrt(r^2 + soft^2)
	real inv_r2 = inv_r * inv_r;             // 1/(r^2 + soft^2)
	real inv_r3 = inv_r2 * inv_r;            // 1/(r^2 + soft^2)^(3/2)
	return -Simulator::G * s1.mass * s2.mass * inv_r3 * r; // force sur s1
}


void Simulator::compute_force(std::vector<Body>& bodies){
	std::size_t n = bodies.size();
	for (std::size_t i = 0; i<n; ++i) {
		bodies[i].force = {0,0,0};
	}
	for (std::size_t i = 0; i<n; ++i) {
		for (std::size_t j = i+1; j<n; ++j) {
			Array grav = gravity(bodies[i], bodies[j]);
			bodies[i].force += grav;
			bodies[j].force -= grav;
		}
	}
}

void Simulator::compute_force_parallel(std::vector<Body>& bodies){
	
	
	std::size_t n = bodies.size();
	for(auto& body : bodies){
		body.force = Array{0,0,0};
	}
	
	std::size_t leftovers = n%NTHREADS;
	std::size_t tasks_per_thread = n/NTHREADS;
	std::vector<args_t> tasks(NTHREADS);
	
	for (std::size_t i = 0; i<NTHREADS-1; ++i) {
		tasks[i].bodies_id.resize(tasks_per_thread);
		for (std::size_t j = 0; j<tasks_per_thread; ++j) {
			tasks[i].bodies_id[j] = j+i*tasks_per_thread;
		}
		tasks[i].bodies = &bodies;
	}
	tasks[NTHREADS-1].bodies_id.resize(tasks_per_thread+leftovers);
	tasks[NTHREADS-1].bodies = &bodies;
	for (std::size_t j = 0; j<tasks_per_thread+leftovers; ++j) {
		tasks[NTHREADS-1].bodies_id[j] = j+(NTHREADS-1)*tasks_per_thread;
	}
	
	
	std::thread threads[NTHREADS];
	for (int i = 0; i<NTHREADS; ++i) {
		threads[i] = std::thread([](auto* args){
			args_t input = *(args_t *) args;
			
			
			std::vector<Array> forces(input.bodies_id.size());
			
			std::size_t idx=0;
			for (std::size_t i : input.bodies_id){
				Array force(0,0,0);
				for (std::size_t j = 0; j<i; ++j) {
					force += gravity((*input.bodies)[i], (*input.bodies)[j]);
				}
				for (std::size_t j = i+1; j<(*input.bodies).size(); ++j) {
					force += gravity((*input.bodies)[i], (*input.bodies)[j]);
				}
				forces[idx] = force;
				++idx;
			}
			idx = 0;
			for (std::size_t i : input.bodies_id){
				(*input.bodies)[i].force = forces[idx];
				++idx;
			}
			
		}, &tasks[i]);
	}
	for (int i = 0; i<NTHREADS; ++i) threads[i].join();
	
}

std::vector<real> Simulator::compute_energy(const std::vector<Body>& bodies){
	std::size_t n = bodies.size();
	real T = 0.0;
	real V = 0.0;
	for (std::size_t i = 0; i<n; ++i) {
		const Body& current_body = bodies[i];
		T += 0.5*current_body.mass*current_body.velocities.norm2();
	}
	for (std::size_t i = 0; i<n; ++i) {
		for (std::size_t j = i+1; j<n; ++j) {
			Array r = bodies[i].positions-bodies[j].positions;
			real r2 = r.norm2() + eps3;
			V -= Simulator::G * bodies[i].mass * bodies[j].mass * 1/std::sqrt(r2);
		}
	}
	return {T,V,T+V};
}

void Simulator::Chrono::start(){
	time = std::chrono::system_clock::now();
	running = true;
}

void Simulator::Chrono::stop(){
	if (running) {
		const auto end = std::chrono::system_clock::now();
		rt = std::chrono::duration_cast<std::chrono::milliseconds>(end - time);
		running = false;
	}
}

std::chrono::milliseconds Simulator::Chrono::runtime() const{
	return rt;
}

void Simulator::Chrono::print() const{
	std::cout << "--- Simulator Runtime: " << (real)rt.count()/1000.0f << "s ---" << std::endl;
}
