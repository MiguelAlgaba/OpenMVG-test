#include <string>
#include <iostream>
#include "openMVG/multiview/solver_resection_kernel.hpp"

#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/multiview/solver_resection_p3p.hpp"
#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"

using namespace openMVG;

#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

using namespace openMVG::robust;


int main(int argc, char *argv[])
{
  openMVG::Mat x = openMVG::Mat(2,14); x << 320.00,366.87,273.13,366.27,273.73,350.45,289.55,324.51,315.49,348.53,291.47,357.17,320.00,282.83,
                                               240.00,232.85,232.85,263.15,263.15,281.34,281.34,282.51,282.51,304.50,304.50,349.97,359.05,349.97;
  openMVG::Mat X = openMVG::Mat(3,14); X << 0.00000,15.27000,-15.27000,14.81000,-14.81000,9.17000,-9.17000,1.37000,-1.37000,8.74000,-8.74000,11.71000,0.00000,-11.71000,                                               0.00000,-2.33000,-2.33000,7.41000,7.41000,12.45000,12.45000,12.90000,12.90000,19.76000,19.76000,34.64000,35.80000,34.64000,                                               0.00000,7.54000,7.54000,3.89000,3.89000,-8.15000,-8.15000,-6.68000,-6.68000,-4.85000,-4.85000,0.66000,-8.44000,0.66000;

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

  {
    typedef openMVG::resection::kernel::SixPointResectionSolver SolverType;

    typedef ACKernelAdaptorResection<
      SolverType, SolverType, UnnormalizerResection, Mat34>
      KernelType;

    std::vector<size_t> pvec_inliers;
    Mat34 P;
    KernelType kernel(x, 640, 480, X);
    double dPrecision = std::numeric_limits<double>::infinity();
    // Robustly estimation of the Projection matrix and it's precision
    std::pair<double,double> ACRansacOut = ACRANSAC(kernel, pvec_inliers,
      1024, &P, dPrecision, true);

    std::cout << "AC-RANSAC have found :" << pvec_inliers.size()
      << " inliers from " << x.cols() << " putatives data." << std::endl
      << " with the upper bound of the precision being : " << ACRansacOut.first << " pixels.\n";

    //-- Extract K, R, t
    Mat3 K, R;
    Vec3 t;
    KRt_From_P(P, &K, &R, &t);
    std::cout << "\nFound K :\n" << K << std::endl
      << "Found R : \n" << R << "\nFound t:" << t.transpose() << std::endl;

    // Display residuals errors :
    Vec residuals(X.cols());
    for (size_t i = 0; i < X.cols(); ++i)
    {
      residuals(i) = (x.col(i) - Project(P, Vec3(X.col(i)))).norm();
    }
    std::cout << "Approximative L2 residual (pixels): " << residuals.norm() << std::endl;

}

  return EXIT_SUCCESS;
}
