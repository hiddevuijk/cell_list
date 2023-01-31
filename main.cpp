

#include "config_file.h"
#include "vec3.h"
#include "systemEDBD.h"
#include "initialize_positions.h"

#include "pair_correlation.h"
#include "density.h"
#include "read_positions.h"

#include "neighbor_list.h"


#include <iostream>
#include <fstream>
#include <vector>
#include <string>


using namespace std;



int main()
{


  Config params("input.txt");

  double Lx =
		params.get_parameter<double>("system_size_x");
  double Ly =
		params.get_parameter<double>("system_size_y");
  double Lz =
		params.get_parameter<double>("system_size_z");


  double r_verlet = 1.1;

  bool pbc_x = false;
  bool pbc_y = pbc_x;
  bool pbc_z = pbc_x;

  bool double_bond = false;
  
  vector<Vec3> positions = read_positions("positions.dat");
  
  for (unsigned int pi = 0; pi < positions.size(); ++pi) {
    Vec3 pos = positions[pi];
    if (!pbc_x) {
      if (pos.x < 0) pos.x = -1 * pos.x;
      pos.x -= Lx * floor(pos.x / Lx);
    }
    if (!pbc_y) {
      if (pos.y < 0) pos.y = -1 * pos.y;
      pos.y -= Ly * floor(pos.y / Ly);
    }
    if (!pbc_z) {
      if (pos.z < 0) pos.z = -1 * pos.z;
      pos.z -= Lz * floor(pos.z / Lz);
    }
    positions[pi] = pos;
  }

  vector<list<unsigned int> > neighbor_list;
  for (unsigned int i = 0; i < 5000; ++i) {
    //neighbor_list = get_neighbor_list(Lx,Ly,Lz,pbc_x,pbc_y,pbc_z,double_bond,r_verlet,positions);
    neighbor_list = get_verlet_list(Lx,Ly,Lz,pbc_x,pbc_y,pbc_z,double_bond,r_verlet,positions);
  }

  for (unsigned int i=0 ; i< positions.size(); ++i) {
    cout << i << ": \t"; 
    for( list<unsigned int>::iterator it =
                            neighbor_list[i].begin();
        it != neighbor_list[i].end(); ++it) {
      cout << *it << " ";  
    }
    cout << "\n";
  }

  
  return 0;
}
