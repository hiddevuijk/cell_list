

#include "config_file.h"
#include "vec3.h"
#include "systemEDBD.h"
#include "initialize_positions.h"

#include "pair_correlation.h"
#include "density.h"
#include "read_positions.h"

#include "cell_list.h"
#include "verlet_list.h"


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


  double r_verlet = 1.2;
  bool double_bond = false;
  
  vector<Vec3> positions = read_positions("positions.dat");

  vector<list<unsigned int> > neighbor_list;
  for (unsigned int i = 0; i < 100; ++i) {
  neighbor_list = get_neighbor_list(Lx,Ly,Lz,r_verlet,positions, double_bond);
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

  vector<vector<unsigned int> > verlet_list;
  vector<unsigned int> n_neighbors;
  for (unsigned int i = 0; i < 100; ++i) {
    get_verlet_list(Lx,Ly,Lz,r_verlet,positions,verlet_list,n_neighbors, double_bond);
  }


  for (unsigned int i=0; i<positions.size(); ++i) {
    cout << i << ":\t";
    for(unsigned int ni=0; ni < n_neighbors[i]; ++ni) {
      cout << verlet_list[i][ni] << " ";
    }
    cout << "\n";
  }
  
  return 0;
}
