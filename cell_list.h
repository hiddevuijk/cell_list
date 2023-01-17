#ifndef GUARD_CELL_LIST_H
#define GUARD_CELL_LIST_H

#include <list>

#include "vec3.h"

namespace neighbor_list {

double distance_squared(Vec3 r1, Vec3 r2,
  double xmin, double Lx,
  double ymin, double Ly,
  double zmin, double Lz)
{
  r1 -= r2;
  if (Lx > 0) r1.x -= Lx * round( (r1.x - xmin) / Lx);
  if (Ly > 0) r1.y -= Ly * round( (r1.y - ymin) / Ly);
  if (Lz > 0) r1.z -= Lz * round( (r1.z - zmin) / Lz);
  return r1.LengthSquared();
}

class VerletList
{

  VerletList(double rv, double xmin, double Lx,
            double ymin, double Ly, double zmin, double Lz)
    : rv_(rv), xmin_(xmin), xmax_(xmin + Lx), Lx_(Lx),
      ymin_(ymin), ymax_(ymin + Ly), Ly_(Ly),
      zmin_(zmin), zmax_(zmin + Lz), Lz_(Lz) {};


  
 private:
  double rv_;
  double xmin_, xmax_;
  double ymin_, ymax_;
  double zmin_, zmax_;
  double Lx_, Ly_, Lz_;

  
  

};


class CellList
{

};

//void generate_verlet_list(
//  const std::vector<Vec3>& positions,
//  double Lx, double Ly, double Lz,
//  double rv,
//  std::vector<std::vector<unsigned int> >& verlet_list,
//  std::vector<unsigned int> number_of_neighbors)
//{
//  unsigned int Nx =  std::floor(Lx / rv);
//  unsigned int Ny =  std::floor(Ly / rv);
//  unsigned int Nz =  std::floor(Lz / rv);
//  unsigned int N = Nx * Ny * Nz;
//  double dx = Lx / Nx;
//  double dy = Ly / Ny;
//  double dz = Lz / Nz;
//
//  std::vector<std::vector<int> > boxes(N);
//
//  unsigned int xi, yi, zi, bi;
//  for (unsigned int i = 0; i < positions.size(); ++i) {
//    xi = std::floor(positions[i].x / dx); 
//    yi = std::floor(positions[i].y / dy); 
//    zi = std::floor(positions[i].z / dz); 
//    bi = zi * Nx * Ny + yi * Nx + xi;
//    boxes[bi].push_back(i);
//  }
//
//};
//


}; // namespace: neighbor_list

#endif //GUARD_CELL_LIST_H
