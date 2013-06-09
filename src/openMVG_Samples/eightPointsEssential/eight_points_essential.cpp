#include <string>
#include <iostream>
#include <fstream>
#include "openMVG/multiview/solver_essential_kernel.hpp"

int read2DPoints(const std::string & projections2DFileName,
                           openMVG::Mat & points2D)
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

int readCalibrationMatrix(const std::string & calibrationMatrixFileName,
                          openMVG::Mat3 & K)
{
  std::string line;
  std::ifstream calibrationMatrixFile (calibrationMatrixFileName.c_str());
  if (calibrationMatrixFile.is_open())
  {
	size_t i=0;
	while(calibrationMatrixFile.good() && i<9)
    {
      std::getline (calibrationMatrixFile,line);	  
      if(calibrationMatrixFile.eof()){break;}
      std::istringstream buffer(line); 
	  while(buffer.good() && i/3<3)
	  {
	    float k;
	    buffer >> k;
		K(i/3,i%3)=k;
	    i++;
	  }
    }
    calibrationMatrixFile.close();
	if(i!=9)
	{
	  std::cerr<<"Invalid format in file: " << calibrationMatrixFileName << std::endl;
	  return EXIT_FAILURE;
	}
  }
  else
  {
    std::cerr<<"Could not open file: " << calibrationMatrixFileName << std::endl;
    return EXIT_FAILURE;
  }

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc!=4)
  {
    std::cout<<"Usage: "<< argv[0] << " <points2D_0> <points2D_1> <K>" << std::endl;
    return EXIT_FAILURE;
  }
  openMVG::Mat x0,x1;
  openMVG::Mat3 K;
  if(read2DPoints(argv[1],x0)){return EXIT_FAILURE;}
  if(read2DPoints(argv[2],x1)){return EXIT_FAILURE;}
  if(readCalibrationMatrix(argv[3],K)){return EXIT_FAILURE;}
	
  typedef openMVG::essential::kernel::EightPointKernel Kernel;	
	
  std::vector<openMVG::Mat3> Es, Rs;  // Essential, Rotation matrix.
  std::vector<openMVG::Vec3> ts;      // Translation matrix.	
	 
  std::cout<< "x0: " << std::endl << x0 << std::endl;	
  std::cout<< "x1: " << std::endl << x1 << std::endl;
  std::cout<< "K: " << std::endl << K << std::endl;	
	
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
