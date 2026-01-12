//
//  main.cpp
//  simulator_galaxie
//
//  Created by Nolwen Dolléans on 08/01/2026.
//

#include <iostream>
#include "simulator.hpp"

int main() {
	Simulator::galaxy galaxy;
	
	galaxy.initialize(5e3, Simulator::IntChoice::Leapfrog);
	
	std::string path = "results/Leapfrog_";
	galaxy.compute(path);
	galaxy.compute_parallel(path);
	
	return EXIT_SUCCESS;
}
