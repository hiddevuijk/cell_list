#ifndef GUARD_VERLET_LIST_H
#define GUARD_VERLET_LIST_H

#include <vector>
#include "vec3.h"

void get_verlet_list(
  double Lx, double Ly, double Lz,
  double r_verlet,
  const std::vector<Vec3>& positions,
  std::vector<std::vector<unsigned int> >& verlet_list,
  std::vector<unsigned int>& number_of_neighbors,
  bool double_bond)
{

  unsigned int n_particles = positions.size();
  verlet_list = std::vector<std::vector<unsigned int> >(
                    n_particles,
                std::vector<unsigned int>(n_particles) );
  number_of_neighbors = std::vector<unsigned int>(n_particles);

  double dist_sq;
  for (unsigned int i = 0; i < n_particles; ++i) {
  for (unsigned int j = i + 1; j < n_particles; ++j) {

    dist_sq = dist2(positions[i], positions[j], Lx, Ly, Lz);
    if (dist_sq < r_verlet * r_verlet) {

      verlet_list[i][ number_of_neighbors[i] ] = j;
      ++number_of_neighbors[i];

      if (double_bond == true) {
        verlet_list[j][ number_of_neighbors[j] ] = i;
        ++number_of_neighbors[j];
      }

    }
  }} 

}

#endif //GUARD_VERLET_LIST_H
