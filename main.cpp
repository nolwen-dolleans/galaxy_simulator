//
//  main.cpp
//  simulator_galaxie
//
//  Created by Nolwen Dolléans on 08/01/2026.
//

#include <iostream>
#include "simulator.hpp"

int main(int argc, char* argv[]) {
	Simulator::galaxy galaxy;
	
	std::size_t size = 1e3;
	
	if (argc == 2) size = atoi(argv[1]);
	
	/*
	galaxy.initialize(size, Simulator::IntChoice::Euler);
	std::string path = "results/Euler_";
	galaxy.compute_parallel(path);
	
	galaxy.initialize(size, Simulator::IntChoice::Leapfrog);
	path = "results/Leapfrog_";
	galaxy.compute_parallel(path);*/
	
	//galaxy.initialize(size, Simulator::IntChoice::Barnes_Hut);
	std::string path = "results/Barnes_Hut_";
	//galaxy.compute(path);
	
	galaxy.initialize(size, Simulator::IntChoice::Barnes_Hut);
	std::string path2 = "results/Barnes_Hut_parallel";
	galaxy.compute_parallel(path2);
	
	return EXIT_SUCCESS;
}
