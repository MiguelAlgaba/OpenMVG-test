#include <string>
#include <iostream>
#include "openMVG/multiview/solver_resection_kernel.hpp"

int main(int argc, char *argv[])
{
  openMVG::Mat2X x = openMVG::Mat2X(2,14); x << 320.00,366.87,273.13,366.27,273.73,350.45,289.55,324.51,315.49,348.53,291.47,357.17,320.00,282.83,
                                               240.00,232.85,232.85,263.15,263.15,281.34,281.34,282.51,282.51,304.50,304.50,349.97,359.05,349.97;
  openMVG::Mat3X X = openMVG::Mat3X(3,14); X << 0.00000,15.27000,-15.27000,14.81000,-14.81000,9.17000,-9.17000,1.37000,-1.37000,8.74000,-8.74000,11.71000,0.00000,-11.71000,
                                               0.00000,-2.33000,-2.33000,7.41000,7.41000,12.45000,12.45000,12.90000,12.90000,19.76000,19.76000,34.64000,35.80000,34.64000,
                                               0.00000,7.54000,7.54000,3.89000,3.89000,-8.15000,-8.15000,-6.68000,-6.68000,-4.85000,-4.85000,0.66000,-8.44000,0.66000;
  std::cout<< "x: " << std::endl << x << std::endl;
  std::cout<< "X: " << std::endl << X << std::endl;
  openMVG::resection::kernel::PoseResectionKernel kernel(x, X);

  size_t samples_[14]={0,1,2,3,4,5,6,7,8,9,10,11,12,13};
  std::vector<size_t> samples(samples_,samples_+14);
  std::vector<openMVG::Mat34> Ps;
  kernel.Fit(samples, &Ps);

  for(size_t i=0; i<Ps.size(); i++)
  {
    std::cout << "Ps: " << std::endl << Ps[i] << std::endl;
  }

  return EXIT_SUCCESS;
}
