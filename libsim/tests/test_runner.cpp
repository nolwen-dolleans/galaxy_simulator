//
//  test_runner.cpp
//  simulator_galaxie
//
//  Created by Nolwen Dolléans on 08/01/2026.
//

#include <stdio.h>
#include <gtest/gtest.h>
#include "simulator.hpp"


TEST(ArrayTests, ArrayTestsConstructor){
	Simulator::Array vec;
	Simulator::Array vec2 = Simulator::Array();
	for (int i = 0; i<3; ++i) EXPECT_EQ(vec[i], 0.0);
	for (int i = 0; i<3; ++i) EXPECT_EQ(vec2[i], 0.0);
}

TEST(ArrayTests, ArrayTestsOpertators){
	Simulator::Array vec;
	Simulator::Array vec2(1.0,2.0,3.0);
	Simulator::Array vec3 = vec + vec2;
	Simulator::Array vec4 = vec3*2.0;
	Simulator::Array vec5 = 2.0*vec3;
	
	for (int i = 0; i<3; ++i) EXPECT_EQ(vec2[i], (double)(i+1));
	for (int i = 0; i<3; ++i) EXPECT_EQ(vec3[i], (double)(i+1));
	for (int i = 0; i<3; ++i) EXPECT_EQ(vec4[i], (double)2*(i+1));
	for (int i = 0; i<3; ++i) EXPECT_EQ(vec5[i], (double)2*(i+1));
}

TEST(ArrayTests, ArrayTestsFunctions){
	Simulator::Array vec(1.0,0.0,0.0);
	Simulator::Array vec2(0.0,1.0,0.0);
	
	Simulator::Array cross(0.0,0.0,1.0);
	Simulator::Array res = Simulator::cross(vec, vec2);
	for (int i = 0; i<3; ++i) EXPECT_EQ(cross[i], res[i]);
	
	double dot = Simulator::dot(vec, vec2);
	EXPECT_EQ(dot, 0.0);
	
	double norm = vec.norm2();
	double norm2 = vec2.norm2();
	
	EXPECT_EQ(norm, norm2);
}
