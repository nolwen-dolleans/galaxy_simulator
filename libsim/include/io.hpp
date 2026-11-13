//
//  simulator.hpp
//  simulator_galaxie
//
//  Created by Nolwen Dolléans on 08/01/2026.
//

#ifndef io_hpp
#define io_hpp

#include <iostream>
#include <type_traits>
#include <fstream>
#include <string>
#include <cstdint>

namespace my_ios{

constexpr std::size_t NPY_HEADER_TOTAL = 128;

template <typename T> constexpr const char* npy_descr();
template <> constexpr const char* npy_descr<float>()  { return "<f4"; }
template <> constexpr const char* npy_descr<double>() { return "<f8"; }

template <typename T>
void npy_write_header(std::ofstream &f, std::size_t steps, std::size_t n){
	std::string header = "{'descr': '" + std::string(npy_descr<T>()) + 
						 "', 'fortran_order': False, 'shape': (" +
						 std::to_string(steps) + ", 3, " + std::to_string(n) + "), }";

	const std::size_t header_len = NPY_HEADER_TOTAL - 10; // "-10" laisse un espace pour la signature
	header.resize(header_len-1, ' ');
	header += '\n';
	
	f.seekp(0);
    f.write("\x93NUMPY\x01\x00", 8);
    const uint16_t hl = static_cast<uint16_t>(header_len);   // little-endian sur Mac
    f.write(reinterpret_cast<const char*>(&hl), 2);
    f.write(header.data(), header.size());

}

template <typename T>
void npy_append_step(std::ofstream& f, const T* x, const T* y, const T* z, std::size_t n){
    f.write(reinterpret_cast<const char*>(x), n * sizeof(T));
    f.write(reinterpret_cast<const char*>(y), n * sizeof(T));
    f.write(reinterpret_cast<const char*>(z), n * sizeof(T));
}




} /* io */


#endif /* io_hpp */

