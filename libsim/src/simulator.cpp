//
//  simulator.cpp
//  simulator_galaxie
//
//  Created by Nolwen Dolléans on 08/01/2026.
//

#include "simulator.hpp"


const real M_disc 	= 5e10; 	//Solar mass
const real M_halo 	= 1e12; 	//Solar mass

const real R_d 		= 3e3;		//parcecs
const real R_max 	= 1.5e4;	//parcecs
const real R_z 		= 3e2;		//parcecs

const real R_s 		= 2e4;		//parcecs
const real c 		= 10.0;
real rho_0 			= 0.0;

const real V_c		= 225;		//percecs/Myr

real dt				= 0;		//Myr
real end			= 0;		//Myr


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

static inline real gaussian(real sigma){
	real u1 = drand48();
	real u2 = drand48();
	return sigma * sqrt(-2.0*log(u1)) * cos(2*M_PI*u2);
}

void Simulator::compute_csts(){
	real T_orbite = (R_max - R_d)/2 * 2 * M_PI / V_c;
	
	end = 10*T_orbite;
	dt = end / N;
	
	
	real a = log(c+1) - c/(c+1);
	rho_0 = M_halo/(4*M_PI*R_s*R_s*R_s*a);
}

void Simulator::init_mass(size_t nBodies, real &m){
	real sigma = 0.05;
	m = M_disc/nBodies * std::max(0.0, 1.0 + sigma*gaussian(1.0));
}


void Simulator::init_position(Simulator::Array& p0){
	real u1 = drand48();
	real u2 = drand48();
	real u3 = drand48();
	real theta = 2*M_PI*u2;
	
	real r_circ = R_d * sqrt(-2.0 * log(1.0 - u1));
	r_circ = std::min(r_circ, R_max);
	
	real x = r_circ*cos(theta);
	real y = r_circ*sin(theta);
	real z = gaussian(R_z);
	
	p0 = {x, y, z};
	
}

static inline real mass_disc(real r){
	real r_Rd = r/R_d;
	return M_disc * (1 - exp(-r_Rd)*(1+r_Rd));
}

static inline real mass_halo(real r){
	real x = r / R_s;
	return 4.0 * M_PI * rho_0 * R_s*R_s*R_s * ( log(1.0 + x) - x/(1.0 + x) );
}

void Simulator::init_velocity(real x, real y, real z, Simulator::Array& v0){
	real r = sqrt(x*x+y*y);
	real inv_r = 1.0 / (r + eps2);
	
	real M_discr = mass_disc(r);
	real M_halor = mass_halo(r);
	real M_total = M_discr + M_halor;
	
	real v_circ = sqrt(Simulator::G * M_total / r);
	real sigma_r = 0.25 * v_circ;
	real sigma_phi = sigma_r / sqrt(2.0);
	real v_phi = v_circ + gaussian(sigma_phi);
	real v_r = gaussian(sigma_r);
	
	v0[0] = (x*v_r-y*v_phi)*inv_r;
	v0[1] = (y*v_r+x*v_phi)*inv_r;
	v0[2] = gaussian(0.5*sigma_r);
	v0 = {0.0,0.0,0.0};
	
}


void Simulator::galaxy::initialize(std::size_t i, const IntChoice& choice){
	compute_csts();
	
	Nbodies = i;
	bodies.resize(i);

	bodies[0].mass = 4e6;
	bodies[0].positions = {0.0, 0.0, 0.0};
	bodies[0].velocities = {0.0, 0.0, 0.0};
	bodies[0].id = 0;
	
	for (std::size_t j = 1; j < i; ++j) {
		init_mass(i, bodies[j].mass);
		init_position(bodies[j].positions);
		init_velocity(bodies[j].positions[0], bodies[j].positions[1], bodies[j].positions[2], bodies[j].velocities);
		bodies[0].id = j;
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
	std::ofstream file(path+"stars_data.csv", std::ios::app);
	std::ofstream energy_file(path+"energy.csv", std::ios::app);
	
	while (t<end) {
		std::cout << "--- Compute Simulator at " << t << " ---\r";
		std::cout.flush();
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
		
		std::vector<real> energy;
		compute_energy(bodies, &energy);
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
	std::ofstream file_(path+"stars_data", std::ios::trunc);
	file_.close();
	std::ofstream file2(path+"stars_data");
	file2 << Nbodies << "\n";
	file2.close();
	
	std::ofstream file3(path+"energy", std::ios::trunc);
	file3.close();
	std::ofstream file4(path+"energy");
	file4 << "temps,T,V,E;\n";
	file4.close();
	
	real t = 0;
	std::ofstream file(path+"stars_data", std::ios::app);
	std::ofstream energy_file(path+"energy", std::ios::app);
	
	
	// --- Allocation ---
	real* x  = (real*) malloc(sizeof(real) * Nbodies);
	real* y  = (real*) malloc(sizeof(real) * Nbodies);
	real* z  = (real*) malloc(sizeof(real) * Nbodies);

	real* vx = (real*) malloc(sizeof(real) * Nbodies);
	real* vy = (real*) malloc(sizeof(real) * Nbodies);
	real* vz = (real*) malloc(sizeof(real) * Nbodies);

	real* fx = (real*) malloc(sizeof(real) * Nbodies);
	real* fy = (real*) malloc(sizeof(real) * Nbodies);
	real* fz = (real*) malloc(sizeof(real) * Nbodies);

	real* m  = (real*) malloc(sizeof(real) * Nbodies);
	
	for (size_t i = 0; i < Nbodies; ++i) {
		x[i]  = bodies[i].positions[0];
		y[i]  = bodies[i].positions[1];
		z[i]  = bodies[i].positions[2];

		vx[i] = bodies[i].velocities[0];
		vy[i] = bodies[i].velocities[1];
		vz[i] = bodies[i].velocities[2];

		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;

		m[i]  = bodies[i].mass;
	}
	
	while (t<end) {
		std::cout << "--- Compute Simulator at " << std::fixed << std::setprecision(1) << t/end*100.0f << "% ---\r";
		std::cout.flush();
		for (std::size_t i = 0; i<Nbodies; ++i) {
			file << x[i] << "," << y[i] << "," << z[i] << ",";
		}
		file << "\n";
		
		integrator->step_parallel(x, y, z, vx, vy, vz, fx, fy, fz, m, dt, Nbodies);
		/*
		std::vector<real> energy;
		compute_energy_parallel(x, y, z, vx, vy, vz, m, Nbodies, &energy);*/
		
		
		t += dt;
	}
	
	for (size_t i = 0; i < Nbodies; ++i) {
		bodies[i].positions[0]  =  x[i];
		bodies[i].positions[1]  =  y[i];
		bodies[i].positions[2]  =  z[i];
		
		bodies[i].velocities[0] = vx[i];
		bodies[i].velocities[1] = vy[i];
		bodies[i].velocities[2] = vz[i];
		
		bodies[i].force[0]	    = fx[i];
		bodies[i].force[1]	    = fy[i];
		bodies[i].force[2]	    = fz[i];


		bodies[i].mass = m[i];
	}
	
	file.close();
	energy_file.close();
	timer.stop();
	timer.print();
	
	// --- Free ---
	free(x);  free(y);  free(z);
	free(vx); free(vy); free(vz);
	free(fx); free(fy); free(fz);
	free(m);
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

void Simulator::Euler::step_parallel(real* x, real* y, real* z, real* vx, real* vy, real* vz, real* fx, real* fy, real* fz, real* m, const real dt, size_t n){
	return;
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

void Simulator::Leapfrog::step_parallel(real* x, real* y, real* z, real* vx, real* vy, real* vz, real* fx, real* fy, real* fz, real* m, const real dt, size_t n){
	return;
}

//##################################################################

int Simulator::get_octant(OctreeNode* node, Body* b) {
	int oct = 0;
	if (b->positions[0] >= node->center[0]) oct |= 1;
	if (b->positions[1] >= node->center[1]) oct |= 2;
	if (b->positions[2] >= node->center[2]) oct |= 4;
	return oct;
}

Simulator::OctreeNode* Simulator::create_child(OctreeNode* parent, int oct) {
	OctreeNode* c = new OctreeNode{};
	real hs = parent->half_size / 2.0;

	c->center[0] = parent->center[0] + hs * ((oct & 1) ? 1 : -1);
	c->center[1] = parent->center[1] + hs * ((oct & 2) ? 1 : -1);
	c->center[2] = parent->center[2] + hs * ((oct & 4) ? 1 : -1);

	c->half_size = hs;
	c->mass = 0.0;
	c->body_count = 0;
	
	
	c->com = {0.0, 0.0, 0.0};
	for(int i = 0; i<MAXBODY_LEAF; ++i) c->bodies[i] = NULL;
	for (int i = 0; i < 8; i++) c->children[i] = NULL;

	return c;
}

Simulator::OctreeNode* Simulator::create_node(const Array& c_in, real half_size){
	OctreeNode* c = (OctreeNode*) new OctreeNode{};
	
	c->center = c_in;
	c->half_size = half_size;
	c->mass = 0.0;
	
	c->body_count = 0;
	for(int i = 0; i<MAXBODY_LEAF; ++i) c->bodies[i] = NULL;
	
	c->com = {0.0, 0.0, 0.0};
	
	for(int i = 0; i<8; ++i) c->children[i] = NULL;
	
	
	return c;
}

void Simulator::insert_body(OctreeNode* root, Body* b) {
	OctreeNode* node = root;
	int depth = 0;

	while (true) {
		// Si c'est une feuille
		if (node->children[0] == nullptr) {
			// Si on peut stocker le corps dans cette feuille
			if (node->body_count < MAXBODY_LEAF) {
				node->bodies[node->body_count++] = b;
				return;
			}
			// Si on a atteint la profondeur max, ne pas dépasser le buffer
			if (depth >= MAXDEPTH) {
				// Éviter tout dépassement mémoire : on remplace le dernier si plein
				if (node->body_count < MAXBODY_LEAF) {
					node->bodies[node->body_count++] = b;
				} else {
					node->bodies[MAXBODY_LEAF - 1] = b;
				}
				return;
			}

			// Subdivision du nœud courant : initialiser les enfants à nullptr
			for (int i = 0; i < 8; ++i) node->children[i] = nullptr;

			// Sauvegarder les corps existants
			Body* existing[MAXBODY_LEAF];
			int count = node->body_count;
			for (int i = 0; i < count; ++i) existing[i] = node->bodies[i];
			node->body_count = 0; // la feuille devient interne

			// Réinsérer proprement tous les corps existants via insert_body
			for (int i = 0; i < count; ++i) {
				insert_body(node, existing[i]);
			}
			// Puis réinsérer le nouveau corps b (la boucle while va continuer)
		}

		// Descendre vers le bon enfant pour le nouveau corps
		int oct = get_octant(node, b);
		if (!node->children[oct])
			node->children[oct] = create_child(node, oct);
		node = node->children[oct];
		depth++;
	}
}

void Simulator::compute_mass(Simulator::OctreeNode* node) {
	if (!node) return;
	
	node->mass = 0.0;
	node->com = {0.0, 0.0, 0.0};
	if (node->children[0] == NULL && node->body_count > 0) {
		for (int i = 0; i<node->body_count; ++i) {
			node->mass += node->bodies[i]->mass;
			node->com += node->bodies[i]->positions * node->bodies[i]->mass;
		}
		node->com *= 1.0/node->mass;
		return;
	}
	

	for (int i = 0; i < 8; i++) {
		if (node->children[i]) {
			compute_mass(node->children[i]);
			node->mass += node->children[i]->mass;
			node->com += node->children[i]->mass * node->children[i]->com;
		}
	}

	if (node->mass > 0) {
		node->com *= 1.0/node->mass;
	}
}

void Simulator::compute_force_BH(OctreeNode* node, Body* b, real theta) {
	if (!node || node->mass == 0.0) return;

	// Distance au centre de masse
	Array d = node->com - b->positions;
	real dist2 = d.norm2() + eps2;
	real dist = sqrt(dist2);

	// Cas feuille : calcul exact
	if (node->children[0] == NULL) {
		for (int i = 0; i < node->body_count; ++i) {
			Body* other = node->bodies[i];
			if (other == b) continue;

			Array r = other->positions - b->positions;
			real r2 = r.norm2() + eps2;
			real inv_r = 1.0 / sqrt(r2);

			real F = G * b->mass * other->mass * inv_r * inv_r;
			b->force += F * r * inv_r;
		}
		return;
	}

	// Critère Barnes–Hut
	if ((node->half_size * 2) / dist < theta) {
		// Approximation multipolaire monopôle
		real inv_dist = 1.0 / dist;
		real F = G * b->mass * node->mass * inv_dist * inv_dist;
		b->force += F * d * inv_dist;
		return;
	}

	// Sinon : descente récursive
	for (int i = 0; i < 8; ++i) {
		if (node->children[i])
			compute_force_BH(node->children[i], b, theta);
	}
}

Simulator::OctreeNode* Simulator::build_octree(std::vector<Body>& bodies) {

	// 1. Bounding box
	Array min = bodies[0].positions;
	Array max = bodies[0].positions;

	for (const Body& b : bodies) {
		min = {std::min(min[0], b.positions[0]),std::min(min[1], b.positions[1]), std::min(min[2], b.positions[2])};
		max = {std::max(max[0], b.positions[0]),std::max(max[1], b.positions[1]), std::max(max[2], b.positions[2])};
	}

	// 2. Cube englobant
	Array center = 0.5 * (min + max);
	real size = std::max({ max[0] - min[0],
							 max[1] - min[1],
							 max[2] - min[2] });

	// 3. Racine
	OctreeNode* root = create_node(center, size * 0.5);

	// 4. Insertion
	for (Body& b : bodies) {
		insert_body(root, &b);
	}

	// 5. Masse & centre de masse
	compute_mass(root);
	
	return root;
}

//##################################################################

int Simulator::get_octant_t(OctreeNode_t* node, real x, real y, real z) {
	int oct = 0;
	if (x >= node->center[0]) oct |= 1;
	if (y >= node->center[1]) oct |= 2;
	if (z >= node->center[2]) oct |= 4;
	return oct;
}

Simulator::OctreeNode_t* Simulator::create_child_t(OctreeNode_t* parent, int oct) {
	OctreeNode_t* c = new OctreeNode_t{};
	real hs = parent->half_size / 2.0;

	c->center[0] = parent->center[0] + hs * ((oct & 1) ? 1 : -1);
	c->center[1] = parent->center[1] + hs * ((oct & 2) ? 1 : -1);
	c->center[2] = parent->center[2] + hs * ((oct & 4) ? 1 : -1);

	c->half_size = hs;
	c->mass = 0.0;
	c->body_count = 0;
	
	
	c->com = {0.0, 0.0, 0.0};
	for(int i = 0; i<MAXBODY_LEAF; ++i) c->index[i] = 0;
	for (int i = 0; i < 8; i++) c->children[i] = NULL;

	return c;
}

Simulator::OctreeNode_t* Simulator::create_node_t(const Array& c_in, real half_size){
	OctreeNode_t* c = new OctreeNode_t{};
	
	c->center = c_in;
	c->half_size = half_size;
	c->mass = 0.0;
	
	c->body_count = 0;
	for(int i = 0; i<MAXBODY_LEAF; ++i) c->index[i] = 0;
	
	c->com = {0.0, 0.0, 0.0};
	
	for(int i = 0; i<8; ++i) c->children[i] = nullptr;
	
	
	return c;
}

void Simulator::insert_body_t(OctreeNode_t* root, size_t id_body, real* x, real* y, real* z) {
	OctreeNode_t* node = root;
	int depth = 0;
	constexpr real eps_jitter = 1e-12; // pour éviter superpositions exactes

	// Ajouter un petit jitter si le corps est exactement sur le centre
	real px = x[id_body];
	real py = y[id_body];
	real pz = z[id_body];

	if (px == node->center[0]) px += eps_jitter;
	if (py == node->center[1]) py += eps_jitter;
	if (pz == node->center[2]) pz += eps_jitter;

	while (true) {
		assert(node != nullptr); // sécurité

		// Feuille
		if (node->children[0] == nullptr) {
			if (node->body_count < MAXBODY_LEAF || depth >= MAXDEPTH) {
				// Insérer le corps ici
				if (node->body_count < MAXBODY_LEAF) {
					node->index[node->body_count++] = id_body;
				} else {
					// feuille saturée, garder le corps ici malgré tout
					// Evite dépassement mémoire
					// Optionnel : stocker un compteur de corps supplémentaires
				}
				return;
			}

			// Subdivision stricte : redistribuer les corps existants
			for (int i = 0; i < node->body_count; ++i) {
				size_t old = node->index[i];
				real ox = x[old];
				real oy = y[old];
				real oz = z[old];

				// Eviter corps exactement sur le centre
				if (ox == node->center[0]) ox += eps_jitter;
				if (oy == node->center[1]) oy += eps_jitter;
				if (oz == node->center[2]) oz += eps_jitter;

				int oct = get_octant_t(node, ox, oy, oz);
				assert(oct >= 0 && oct < 8);

				if (!node->children[oct])
					node->children[oct] = create_child_t(node, oct);

				OctreeNode_t* child = node->children[oct];
				if (child->body_count < MAXBODY_LEAF)
					child->index[child->body_count++] = old;
				// sinon : corps coincés, restent dans la feuille actuelle
			}
			node->body_count = 0;
		}

		// Descendre vers le bon octant pour le nouveau corps
		int oct = get_octant_t(node, px, py, pz);
		assert(oct >= 0 && oct < 8);

		if (!node->children[oct])
			node->children[oct] = create_child_t(node, oct);

		node = node->children[oct];
		depth++;

		// Sécurité : ne jamais dépasser MAXDEPTH
		if (depth > MAXDEPTH) {
			// Stocker le corps dans cette feuille même si saturée
			if (node->body_count < MAXBODY_LEAF)
				node->index[node->body_count++] = id_body;
			return;
		}
	}
}
void Simulator::compute_mass_t(Simulator::OctreeNode_t* node, real * x, real * y, real * z, real * m) {
	if (!node) return;
	
	node->mass = 0.0;
	node->com = {0.0, 0.0, 0.0};
	if (node->children[0] == NULL && node->body_count > 0) {
		for (int i = 0; i<node->body_count; ++i) {
			node->mass += m[node->index[i]];
			node->com += {x[node->index[i]]*m[node->index[i]], y[node->index[i]]*m[node->index[i]], z[node->index[i]]*m[node->index[i]]};
		}
		node->com *= 1.0/node->mass;
		return;
	}
	

	for (int i = 0; i < 8; i++) {
		if (node->children[i]) {
			compute_mass_t(node->children[i],x, y, z, m);
			node->mass += node->children[i]->mass;
			node->com += node->children[i]->mass * node->children[i]->com;
		}
	}

	if (node->mass > 0) {
		node->com *= 1.0/node->mass;
	}
}

void Simulator::compute_force_BH_t(OctreeNode_t* node, real * x, real * y, real * z, real * fx, real * fy, real * fz, real * m, size_t id_body, real theta) {
	if (!node || node->mass == 0.0) return;
	
	real px = x[id_body];
	real py = y[id_body];
	real pz = z[id_body];
	
	// Distance au centre de masse
	real dx = node->com[0] - px;
	real dy = node->com[1] - py;
	real dz = node->com[2] - pz;
	real dist2 = dx*dx + dy*dy + dz*dz + eps2;
	real dist = sqrt(dist2);

	// Cas feuille : calcul exact
	if (node->children[0] == NULL) {
		for (int i = 0; i < node->body_count; ++i) {
			size_t other = node->index[i];
			if (other == id_body) continue;
			real rx = x[other] - px;
			real ry = y[other] - py;
			real rz = z[other] - pz;
			real r2 = rx*rx+ry*ry+rz*rz + eps2;
			real inv_r = 1.0 / sqrt(r2);
			
			real F = G * m[id_body] * m[other] * inv_r * inv_r;
			fx[id_body] += F * rx * inv_r;
			fy[id_body] += F * ry * inv_r;
			fz[id_body] += F * rz * inv_r;
		}
		return;
	}

	// Critère Barnes–Hut
	if ((node->half_size * 2) / dist < theta) {
		// Approximation multipolaire monopôle
		real inv_dist = 1.0 / dist;
		real F = G * m[id_body] * node->mass * inv_dist * inv_dist;
		fx[id_body] += F * dx * inv_dist;
		fy[id_body] += F * dy * inv_dist;
		fz[id_body] += F * dz * inv_dist;
		return;
	}

	// Sinon : descente récursive
	for (int i = 0; i < 8; ++i) {
		if (node->children[i])
			compute_force_BH_t(node->children[i], x, y, z, fx, fy, fz, m, id_body, theta);
	}
}

Simulator::OctreeNode_t* Simulator::build_octree_t(real * x, real * y, real * z, real * m, size_t n) {
	
	// 1. Bounding box
	real min[3] = {x[0], y[0], z[0]};
	real max[3] = {x[0], y[0], z[0]};

	float64x2_t minx = vdupq_n_f64(1e308);
	float64x2_t miny = vdupq_n_f64(1e308);
	float64x2_t minz = vdupq_n_f64(1e308);
	float64x2_t maxx = vdupq_n_f64(-1e308);
	float64x2_t maxy = vdupq_n_f64(-1e308);
	float64x2_t maxz = vdupq_n_f64(-1e308);
	for (size_t b = 0; b<n-2; b+=2) {
		float64x2_t vx = vld1q_f64(x + b);
		float64x2_t vy = vld1q_f64(y + b);
		float64x2_t vz = vld1q_f64(z + b);
		minx = vminq_f64(minx, vx);
		miny = vminq_f64(miny, vy);
		minz = vminq_f64(minz, vz);
		
		maxx = vmaxq_f64(maxx, vx);
		maxy = vmaxq_f64(maxy, vy);
		maxz = vmaxq_f64(maxz, vz);
	}
	min[0] = fmin(vgetq_lane_f64(minx, 0), vgetq_lane_f64(minx, 1));
	min[1] = fmin(vgetq_lane_f64(miny, 0), vgetq_lane_f64(miny, 1));
	min[2] = fmin(vgetq_lane_f64(minz, 0), vgetq_lane_f64(minz, 1));
	
	max[0] = fmax(vgetq_lane_f64(maxx, 0), vgetq_lane_f64(maxx, 1));
	max[1] = fmax(vgetq_lane_f64(maxy, 0), vgetq_lane_f64(maxy, 1));
	max[2] = fmax(vgetq_lane_f64(maxz, 0), vgetq_lane_f64(maxz, 1));
	for (size_t b = (n - 2); b<n; ++b) {
		min[0] = (real) std::min(min[0], x[b]);
		min[1] = (real) std::min(min[1], y[b]);
		min[2] = (real) std::min(min[2], z[b]);
		max[0] = (real) std::max(max[0], x[b]);
		max[1] = (real) std::max(max[1], y[b]);
		max[2] = (real) std::max(max[2], z[b]);
	}
	
	// 2. Cube englobant
	Simulator::Array center = {0.5f * (min[0] + max[0]), 0.5f * (min[1] + max[1]), 0.5f * (min[2] + max[2])};
	real size = std::max({ max[0] - min[0], max[1] - min[1], max[2] - min[2] });

	// 3. Racine
	OctreeNode_t* root = create_node_t(center, size * 0.5);

	// 4. Insertion
	for (size_t b = 0; b<n; ++b) {
		insert_body_t(root, b, x, y, z);
	}

	// 5. Masse & centre de masse
	compute_mass_t(root, x, y, z, m);
	
	return root;
}

void Simulator::free_tree(OctreeNode* node) {
	if (!node) return;
	for (int i = 0; i < 8; ++i)
		free_tree(node->children[i]);
	delete node;
}

void Simulator::free_tree_t(OctreeNode_t* node) {
	if (!node) return;

	for (int i = 0; i < 8; ++i) {
		if (node->children[i]) {
			free_tree_t(node->children[i]);
			node->children[i] = nullptr;
		}
	}

	delete node;
}

void Simulator::Barnes_Hut::step(std::vector<Body>& bodies, real dt) {
	real thta = 0.75;
	// 1. Forces
	OctreeNode* tree = build_octree(bodies);

	for (Body& b : bodies){
		b.force = {0.0, 0.0, 0.0};
		compute_force_BH(tree, &b, thta);
	}

	// 2. Kick (dt/2)
	for (Body& b : bodies){
		b.velocities += (b.force * (1.0/b.mass)) * (dt * 0.5);
		b.positions += b.velocities * dt;
	}


	free_tree(tree);

	// 4. Forces au temps t + dt
	tree = build_octree(bodies);

	for (Body& b : bodies){
		b.force = {0.0, 0.0, 0.0};
		compute_force_BH(tree, &b, thta);
	}
	// 5. Kick (dt/2)
	for (Body& b : bodies)
		b.velocities += (b.force * (1.0/b.mass)) * (dt * 0.5);

	free_tree(tree);
}

struct ForceTask {
	std::vector<Simulator::Body>* bodies;
	Simulator::OctreeNode* tree;
	real theta;
	std::atomic<int>* next;
	int n;
};

struct LeapfrogTask {
	std::vector<Simulator::Body>* bodies;
	size_t next;
	real dt;
	int n;
	Simulator::OctreeNode * tree;
};

struct LeapfrogTask_BH {
	Simulator::Body body;
	real dt;
	Simulator::OctreeNode tree;
};

void Simulator::Barnes_Hut::step_parallel(std::vector<Body>& bodies, const real dt){
	return;
}



real computeHaloAcceleration(real * x, real * y, real * z, size_t id_body)
{
	real r = std::sqrt(x[id_body]*x[id_body] + y[id_body]*y[id_body] + z[id_body]*z[id_body]);
	if (r < eps2) return 0.0;

	real M_r = M_halo * (std::log(1 + r/R_s) - r/(R_s + r)) / (std::log(1 + c) - c/(1 + c));

	return Simulator::G * M_r / (r*r*r);
}


void Simulator::Barnes_Hut::step_parallel(real* x, real* y, real* z, real* vx, real* vy, real* vz, real* fx, real* fy, real* fz, real* m, const real dt, size_t n)
{
    real thta = 0.75;
    OctreeNode_t* tree = nullptr;

    // Parallel region over the whole step
    #pragma omp parallel shared(tree)
    {
		// Build tree once
        #pragma omp single
        {
            if (tree) { free_tree_t(tree); tree = nullptr; }
            tree = build_octree_t(x, y, z, m, n);
        }
        #pragma omp barrier

        // 1. Forces at t
        #pragma omp for schedule(static)
        for (size_t i = 0; i < n; ++i) {
            fx[i] = 0.0;
            fy[i] = 0.0;
            fz[i] = 0.0;
            compute_force_BH_t(tree, x, y, z, fx, fy, fz, m, i, thta);
			real factor_DM = computeHaloAcceleration(x, y, z, i);
			fx[i] += -factor_DM * x[i] * m[i];
			fy[i] += -factor_DM * y[i] * m[i];
			fz[i] += -factor_DM * z[i] * m[i];
        }

        // Free tree once
        #pragma omp single
        {
            free_tree_t(tree);
            tree = nullptr;
        }
        #pragma omp barrier

        // 2. Kick + Drift
        const size_t vecEnd = (n >= 2) ? (n - (n % 2)) : 0;
        float64x2_t ones = vdupq_n_f64(1.0);
        float64x2_t udt  = vdupq_n_f64(dt * 0.5);
        float64x2_t Udt  = vdupq_n_f64(dt);

        #pragma omp for schedule(static)
        for (size_t i = 0; i < vecEnd; i += 2) {
            float64x2_t ux = vld1q_f64(x + i);
            float64x2_t uy = vld1q_f64(y + i);
            float64x2_t uz = vld1q_f64(z + i);

            float64x2_t fux = vld1q_f64(fx + i);
            float64x2_t fuy = vld1q_f64(fy + i);
            float64x2_t fuz = vld1q_f64(fz + i);

            float64x2_t um  = vld1q_f64(m + i);

            float64x2_t vux = vld1q_f64(vx + i);
            float64x2_t vuy = vld1q_f64(vy + i);
            float64x2_t vuz = vld1q_f64(vz + i);

            float64x2_t inv_um = vdivq_f64(ones, um);

            vux = vaddq_f64(vux, vmulq_f64(fux, vmulq_f64(inv_um, udt)));
            vuy = vaddq_f64(vuy, vmulq_f64(fuy, vmulq_f64(inv_um, udt)));
            vuz = vaddq_f64(vuz, vmulq_f64(fuz, vmulq_f64(inv_um, udt)));

            ux = vaddq_f64(ux, vmulq_f64(vux, Udt));
            uy = vaddq_f64(uy, vmulq_f64(vuy, Udt));
            uz = vaddq_f64(uz, vmulq_f64(vuz, Udt));

            vst1q_f64(x + i, ux);
            vst1q_f64(y + i, uy);
            vst1q_f64(z + i, uz);

            vst1q_f64(vx + i, vux);
            vst1q_f64(vy + i, vuy);
            vst1q_f64(vz + i, vuz);
        }

        #pragma omp for schedule(static)
        for (size_t i = vecEnd; i < n; ++i) {
            const real inv_m = 1.0 / m[i];
            vx[i] += fx[i] * inv_m * (dt * 0.5);
            vy[i] += fy[i] * inv_m * (dt * 0.5);
            vz[i] += fz[i] * inv_m * (dt * 0.5);
            x[i]  += vx[i] * dt;
            y[i]  += vy[i] * dt;
            z[i]  += vz[i] * dt;
        }

        // 3. Forces at t + dt
        #pragma omp single
        {
            if (tree) { free_tree_t(tree); tree = nullptr; }
            tree = build_octree_t(x, y, z, m, n);
        }
        #pragma omp barrier

        #pragma omp for schedule(static)
        for (size_t i = 0; i < n; ++i) {
            fx[i] = 0.0;
            fy[i] = 0.0;
            fz[i] = 0.0;
            compute_force_BH_t(tree, x, y, z, fx, fy, fz, m, i, thta);
			real factor_DM = computeHaloAcceleration(x, y, z, i);
			fx[i] += -factor_DM * x[i] * m[i];
			fy[i] += -factor_DM * y[i] * m[i];
			fz[i] += -factor_DM * z[i] * m[i];
        }

        // 4. Kick final
        #pragma omp for schedule(static)
        for (size_t j = 0; j < vecEnd; j += 2) {
            float64x2_t fux = vld1q_f64(fx + j);
            float64x2_t fuy = vld1q_f64(fy + j);
            float64x2_t fuz = vld1q_f64(fz + j);

            float64x2_t um  = vld1q_f64(m + j);

            float64x2_t vux = vld1q_f64(vx + j);
            float64x2_t vuy = vld1q_f64(vy + j);
            float64x2_t vuz = vld1q_f64(vz + j);

            float64x2_t inv_um = vdivq_f64(ones, um);

            vux = vaddq_f64(vux, vmulq_f64(fux, vmulq_f64(inv_um, udt)));
            vuy = vaddq_f64(vuy, vmulq_f64(fuy, vmulq_f64(inv_um, udt)));
            vuz = vaddq_f64(vuz, vmulq_f64(fuz, vmulq_f64(inv_um, udt)));

            vst1q_f64(vx + j, vux);
            vst1q_f64(vy + j, vuy);
            vst1q_f64(vz + j, vuz);
        }

        #pragma omp for schedule(static)
        for (size_t j = vecEnd; j < n; ++j) {
            const real inv_m = 1.0 / m[j];
            vx[j] += fx[j] * inv_m * (dt * 0.5);
            vy[j] += fy[j] * inv_m * (dt * 0.5);
            vz[j] += fz[j] * inv_m * (dt * 0.5);
        }

        // Free the tree once before leaving parallel region
        #pragma omp single
        {
            if (tree) { free_tree_t(tree); tree = nullptr; }
        }
    }
}

//##################################################################

void Simulator::compute_force(std::vector<Body>& bodies){
	std::size_t n = bodies.size();
	for (std::size_t i = 0; i<n; ++i) {
		bodies[i].force = {0,0,0};
	}
	for (std::size_t i = 0; i<n; ++i) {
		for (std::size_t j = i+1; j<n; ++j) {
			Array grav;
			gravity(bodies[i], bodies[j], &grav);
			bodies[i].force += grav;
			bodies[j].force -= grav;
		}
	}
}

std::vector<std::pair<size_t, size_t>> balanced_row_partition(size_t n, size_t nthreads){
	std::vector<std::pair<size_t, size_t>> result;

	real total = n * (n - 1) / 2.0;
	real chunk = total / nthreads;

	auto prefix = [&](size_t i) {
		return i * (n - 1) - i * (i - 1) / 2.0;
	};

	size_t current = 0;

	for (size_t t = 0; t < nthreads; ++t) {
		size_t start = current;

		real target = (t + 1) * chunk;

		while (current + 1 < n && prefix(current + 1) < target)
			++current;

		size_t end = (t == nthreads - 1) ? n - 1 : current;

		result.emplace_back(start, end);
		++current;
	}

	return result;
}

void Simulator::compute_force_parallel(std::vector<Body>& bodies)
{
	
	return;
}

void Simulator::compute_energy(const std::vector<Body>& bodies, std::vector<real> * energy){
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
			real r2 = r.norm2() + eps2;
			V -= Simulator::G * bodies[i].mass * bodies[j].mass * 1/std::sqrt(r2);
		}
	}
	*energy = {T,V,T+V};
}

void Simulator::compute_energy_parallel(real * x, real * y, real * z, real * vx, real * vy, real * vz, real * m, size_t n, std::vector<real> * energy){
	real T = 0.0;
	real V = 0.0;
#pragma omp parallel reduction(+:T)
{
	float64x2_t uT = vdupq_n_f64(0.0);

#pragma omp for schedule(static)
	for (size_t i = 0; i < n-2; i += 2) {

		float64x2_t vux = vld1q_f64(vx + i);
		float64x2_t vuy = vld1q_f64(vy + i);
		float64x2_t vuz = vld1q_f64(vz + i);
		float64x2_t um  = vld1q_f64(m  + i);

		float64x2_t uv2 = vaddq_f64(vmulq_f64(vux, vux),vaddq_f64(vmulq_f64(vuy, vuy),vmulq_f64(vuz, vuz)));

		float64x2_t half = vdupq_n_f64(0.5);

		uT = vfmaq_f64(uT, um, vmulq_f64(uv2, half));
	}

	T += vgetq_lane_f64(uT, 0) + vgetq_lane_f64(uT, 1);
}
#pragma omp parallel for reduction(+:T) schedule(static)
	for (std::size_t i = n-2; i<n; ++i) {
		real r2 = vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
		T += 0.5*m[i]*r2;
	}
	
#pragma omp parallel for schedule(static) reduction(+:V)
	for (std::size_t i = 0; i<n; ++i) {
		for (std::size_t j = i+1; j<n; ++j) {
			real dx = x[i] - x[j];
			real dy = y[i] - y[j];
			real dz = z[i] - z[j];

			real r2 = dx*dx + dy*dy + dz*dz + eps2;
			real inv_r2 = 1/std::sqrt(r2);
			V -= Simulator::G * m[i] * m[j] * inv_r2;
		}
	}
	*energy = {T,V,T+V};
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

