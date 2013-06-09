#include <string>
#include <iostream>
#include "openMVG/multiview/solver_essential_kernel.hpp"

int main(int argc, char *argv[])
{
  typedef openMVG::essential::kernel::EightPointKernel Kernel;	
	
  std::vector<openMVG::Mat3> Es, Rs;  // Essential, Rotation matrix.
  std::vector<openMVG::Vec3> ts;      // Translation matrix.
  openMVG::Mat3 K; K << 2969.23348018,0,1632,0,2956.21052631,1224,0,0,1;
	
  openMVG::Mat x0(2,8);
  x0 << 1463.93,2465.00,870.66,1753.69,1466.39,1530.39,2083.59,2131.50,
1168.22,501.67,377.56,1090.73,1094.44,347.25,268.57,1403.58;
  openMVG::Mat x1(2,8);
  x1 << 1527.114,2499.633,1343.952,2029.253,1851.069,1987.251,2512.056,2068.323,
1170.306,307.911,489.691,1080.936,1093.521,316.804,79.326,1428.752;
  
  std::cout<< "x0: " << std::endl << x0 << std::endl;	
  std::cout<< "x1: " << std::endl << x1 << std::endl;	
	
  Kernel kernel(x0, x1, K, K);
  std::vector<size_t> samples;
  for (size_t k = 0; k < Kernel::MINIMUM_SAMPLES; ++k) {
    samples.push_back(k);
  }
  kernel.Fit(samples, &Es);
	
  // Recover rotation and translation from E.
  Rs.resize(Es.size());
  ts.resize(Es.size());
  for (int s = 0; s < Es.size(); ++s)
  {
    std::cout << "E: " << std::endl << Es[s] << std::endl;
	
	for (int j = 0; j < 4 ; j++)
	{
	  openMVG::MotionFromEssential(Es[s],&Rs,&ts);
	  std::cout << "R" << j << ": " << std::endl << Rs[s] << std::endl;
	  std::cout << "t" << j << ": " << std::endl << ts[s] << std::endl;
    }
  }
	
  return EXIT_SUCCESS;
}
