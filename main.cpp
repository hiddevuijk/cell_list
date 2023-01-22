

#include "config_file.h"
#include "vec3.h"
#include "systemEDBD.h"
#include "initialize_positions.h"

#include "pair_correlation.h"
#include "density.h"
#include "read_positions.h"

#include "cell_list.h"

#include <iostream>
#include <vector>
#include <string>


using namespace std;



int main()
{

  Config params("input.txt");

  double system_size_x =
		params.get_parameter<double>("system_size_x");
  double system_size_y =
		params.get_parameter<double>("system_size_y");
  double system_size_z =
		params.get_parameter<double>("system_size_z");


  double verlet_list_radius =
		params.get_parameter<double>("verlet_list_radius");


  std::vector<Vec3> positions = read_positions("positions.dat");

  std::vector<std::vector<unsigned int> > neighbor_list; 
  std::vector<unsigned int> n_neighbor_list; 

  get_neighbor_list2(system_size_x, system_size_y, system_size_z,
            verlet_list_radius,
            positions, neighbor_list, n_neighbor_list);

  
  return 0;
}
