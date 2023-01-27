

#include "config_file.h"
#include "vec3.h"
#include "systemEDBD.h"
#include "initialize_positions.h"

#include "pair_correlation.h"
#include "density.h"
#include "read_positions.h"

#include "cell_list.h"

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


  double r_verlet = 1.01;
  
  vector<Vec3> positions;
  positions.push_back( Vec3(0,0,0));
  positions.push_back( Vec3(1.002, 0.0, 0.0));

  vector<list<unsigned int> > neighbor_list =  
      get_neighbor_list(Lx,Ly,Lz,r_verlet,positions);

  ofstream pout;
  pout.open("positions.dat");
  for (unsigned i = 0; i < positions.size(); ++i) {
    pout << positions[i]. x << "\t" 
         << positions[i]. y << "\t" 
         << positions[i]. z << "\n" ;
  }
  pout.close();


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
