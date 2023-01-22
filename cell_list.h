#ifndef GUARD_CELL_LIST_H
#define GUARD_CELL_LIST_H


#include <vector>
#include <list>
#include "vec3.h"

struct Cell
{
  std::list<unsigned int> particles;
  std::list<*Cell> neighboring_cells;
};

std::vector<std::list<unsigned int> > get_cell_list(
    double Lx, double Ly, double Lz,
    double r_verlet,
    const std::vector<Vec3>& positions)
{

  unsigned int n_particles = positions.size();

  unsigned int n_cells_x = std::floor(Lx / r_verlet);
  unsigned int n_cells_y = std::floor(Ly / r_verlet);
  unsigned int n_cells_z = std::floor(Lz / r_verlet);
  unsigned int n_cells = n_cells_x * n_cells_y * n_cells_z;

  double delta_x = Lx / n_cells_x;
  double delta_y = Ly / n_cells_y;
  double delta_z = Lz / n_cells_z;
 
  std::vector<Cell> cell_list;

  unsigned int xi, yi, zi, cell_index_i, cell_index_j;
  // assign neighboring cells to cells in cell_list
  for (unsigned int xj = 0; xj < n_cells_x; ++xj) {
  for (unsigned int yj = 0; yj < n_cells_y; ++yj) {
  for (unsigned int zj = 0; zj < n_cells_z; ++zj) {
    cell_index_j = zj*n_cells_x*n_cells_y + yj*n_cells_x + xj;

    xi = (xj + 1) % n_cells_x;  
    yi = yj;
    zi = zj;
    cell_index_i = zi*n_cells_x*n_cells_y + yi*n_cells_x + xi;
    if (cell_index_i > cell_index_j) {
      cell_list[cell_index_i].push_back(&cell_list[cell_index_j]);
    }

    xi = (xj - 1) % n_cells_x;  
    yi = yj;
    zi = zj;
    cell_index_i = zi*n_cells_x*n_cells_y + yi*n_cells_x + xi;
    if (cell_index_i > cell_index_j) {
      cell_list[cell_index_i].push_back(&cell_list[cell_index_j]);
    }
    xi = (xj - 1) % n_cells_x;  
    yi = yj;
    zi = zj;
    cell_index_i = zi*n_cells_x*n_cells_y + yi*n_cells_x + xi;
    if (cell_index_i > cell_index_j) {
      cell_list[cell_index_i].push_back(&cell_list[cell_index_j]);
    }
  
  }}}

  // assign cells to particles 
  // xi = xi_th cell in the x direction
  for (unsigned int pi = 0; pi < n_particles; ++pi) {
    xi = std::floor( positions[pi].x / delta_x);
    yi = std::floor( positions[pi].y / delta_y);
    zi = std::floor( positions[pi].z / delta_z);
    cell_index_i = zi*n_cells_x*n_cells_y + yi*n_cells_x + xi;

    // add particle number to particle list of cell 
    cell_list[cell_index_i].particles.push_back(pi);;
  }

  return cell_list;
}


std::vector<unsigned int> get_cell_neighbor_list(
                          unsigned int n_cells_x,
                          unsigned int n_cells_y,
                          unsigned int n_cells_z)
{
  unsigned int n_cells = n_cells_x * n_cells_y * n_cells_z;

  std::vector<std::vector<unsigned int> >
      cell_neighbor_list(n_cells, std::vector<unsigned int>(27));

  for (unsigned int xi = 0; xi < n_cells_x; ++xi) {
  for (unsigned int yi = 0; yi < n_cells_y; ++yi) {
  for (unsigned int zi = 0; zi < n_cells_z; ++zi) {
      ci =  zi * n_cells_y * n_cells_x
          + yi * n_cells_x + xi;
      cell_neighbor_list[ci][0] = ci;

      ci =  zi * n_cells_y * n_cells_x
          + yi * n_cells_x + xi;
  }}}

  return cell_neighbor_list;
}

//std::vector<std::list<unsigned int> > get_neighbor_list(
//    double Lx, double Ly, double Lz,
//    double r_verlet,
//    const std::vector<Vec3>& positions)
//{
//  unsigned int n_particles = positions.size();
//
//  std::vector<std::list<unsigned int> > cell_list =
//      get_cell_list(Lx, Ly, Lz, r_verlet, positions);
//
//  unsigned int n_cells = cell_list.size();
//
//  std::vector<std::list<unsigned int> >
//              neighbor_list(n_particles);
//
//
//  for (unsigned ci = 0; ci < n_cells; ++ci) { 
//    
//  }
//
//}

void get_neighbor_list2(
      double system_size_x,
      double system_size_y,
      double system_size_z,
      double verlet_list_radius, 
      const std::vector<Vec3>& positions,
      std::vector<std::vector<unsigned int> > neighbor_list,
      std::vector<unsigned int> n_neighbor_list)
{

}



#endif  // GUARD_CELL_LIST_H
