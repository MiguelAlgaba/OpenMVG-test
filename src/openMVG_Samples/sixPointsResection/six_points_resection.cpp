#include <string>
#include <iostream>
#include <fstream>

#include <Eigen/QR>
#include "openMVG/multiview/solver_resection_kernel.hpp"

int read3DPoints(const std::string & points3DFileName,
                           openMVG::Mat3X & points3D)
{
  std::vector<openMVG::Vec3> points3D_aux;

  // Read the coordinates of the projections of the 3D points
  std::string line;
  std::ifstream points3DFile (points3DFileName.c_str());
  if (points3DFile.is_open())
  {
    while ( points3DFile.good())
    {
      std::getline (points3DFile,line);
      std::istringstream buffer(line); size_t idx; buffer >> idx; // TODO: Do not ignore the 3D point index
      float x,y,z; buffer >> x; buffer >> y; buffer >> z; buffer.clear();
      if(points3DFile.eof()){break;}

      points3D_aux.push_back(openMVG::Vec3(x,y,z));
    }
    points3DFile.close();
  }
  else
  {
    std::cerr<<"Could not open file: " << points3DFileName << std::endl;
    return EXIT_FAILURE;
  }

  points3D.resize(3,points3D_aux.size());
  for(size_t i=0;i<points3D_aux.size();i++)
  {
    points3D.block<3,1>(0,i) = points3D_aux[i];
  }

  return 0;
}

int read2DPoints(const std::string & projections2DFileName,
                           openMVG::Mat2X & points2D)
{
    std::vector<openMVG::Vec2> points2D_aux;

  // Read the coordinates of the projections of the 3D points
  std::string line;
  std::ifstream points2DFile (projections2DFileName.c_str());
  if (points2DFile.is_open())
  {
    while ( points2DFile.good())
    {
      std::getline (points2DFile,line);
      std::istringstream buffer(line); size_t idx; buffer >> idx; // TODO: Do not ignore the 3D point index
      float x,y; buffer >> x; buffer >> y; buffer.clear();
      if(points2DFile.eof()){break;}

      points2D_aux.push_back(openMVG::Vec2(x,y));
    }
    points2DFile.close();
  }
  else
  {
    std::cerr<<"Could not open file: " << projections2DFileName << std::endl;
    return EXIT_FAILURE;
  }

  points2D.resize(2,points2D_aux.size());
  for(size_t i=0;i<points2D_aux.size();i++)
  {
    points2D.block<2,1>(0,i) = points2D_aux[i];
  }

  return 0;
}

void RQDecomposition(const openMVG::Mat34 & P,
                                     openMVG::Mat3 & R,
                                     openMVG::Mat3 & Q)
{
  openMVG::Mat4 P_aux;
  P_aux.block<3,4>(0,0) = P.block<3,4>(0,0);
  P_aux.block<1,4>(3,0) << 0,0,0,1;
  openMVG::Mat4 ReverseRows;

  ReverseRows << 0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0;

  Eigen::HouseholderQR<openMVG::Mat4> qr( (ReverseRows * P_aux).transpose() );
  openMVG::Mat4 Q_aux = qr.householderQ();
  openMVG::Mat4 R_aux = qr.matrixQR().triangularView<Eigen::Upper>();
  R_aux = ReverseRows * R_aux.transpose() * ReverseRows;
  Q_aux = ReverseRows * Q_aux.transpose();

  // Adjust R and Q so that the diagonal elements of R are +ve
  for(size_t j=0;j<4;j++)
  {
    if(R_aux(j,j) < 0)
    {
      R_aux.block<4,1>(0,j) = -R_aux.block<4,1>(0,j);
      Q_aux.block<1,4>(j,0) = -Q_aux.block<1,4>(j,0);
      }
  }

  R.block<3,3>(0,0) = R_aux.block<3,3>(0,0);
  Q.block<3,3>(0,0) = Q_aux.block<3,3>(0,0);
}

int main(int argc, char *argv[])
{
  if(argc!=3)
  {
    std::cout<<"Usage: "<< argv[0] << " <points2D> <points3D>" << std::endl;
    return EXIT_FAILURE;
  }
  openMVG::Mat2X x;
  read2DPoints(argv[1],x);
  openMVG::Mat3X X;
  read3DPoints(argv[2],X);

  std::cout<< "x: " << std::endl << x << std::endl;
  std::cout<< "X: " << std::endl << X << std::endl;
  openMVG::resection::kernel::PoseResectionKernel kernel(x, X);

  std::vector<size_t> samples(x.cols(),x.cols());
  for(size_t i=0;i<x.cols();i++)
  {
    samples[i]=i;
  }

  std::vector<openMVG::Mat34> Ps;
  kernel.Fit(samples, &Ps);

  for(size_t i=0; i<Ps.size(); i++)
  {
    std::cout << "Ps: " << std::endl << Ps[i] << std::endl;

    openMVG::Mat3 R,Q;
    RQDecomposition(Ps[i],R,Q);

    openMVG::Mat3 Rot = Q;
    std::cout << "Rot: " << std::endl << Rot << std::endl;
  }

  return EXIT_SUCCESS;
}
