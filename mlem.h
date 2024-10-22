#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include "image.h"
#include "acqphy.h"
#include "bgeot_ftool.h"
 

class MLEM
{
 private :
  // iterations
 // int data_type;
  int nb_iter;
  int first_iter;

  double alpha_TV;
  int Niter_TV;
  // image dimensions
  int nb_vox[3];
  int size_filter[3];
  int size_zeropad[3];
  Point voxel_length;
  Point corner;
  Image image;

 public:

 MLEM(): nb_iter(0){}

  void set_iterations(int iter, int first)
  {
    nb_iter = iter;
    first_iter = first;
  }

  void set_volume(int nbvox[], const Point& voxellength, const Point& corn)
  {
    nb_vox[0]=nbvox[0];
    nb_vox[1]=nbvox[1];
    nb_vox[2]=nbvox[2];
    voxel_length = voxellength;
    corner = corn;
  }


  void set_filter(int nbvox[])
  {
    size_filter[0]=nbvox[0];
    size_filter[1]=nbvox[1];
    size_filter[2]=nbvox[2];
   // std::cout <<"test filter size"<< size_filter[2]<<". \n" ;
  }

  void set_zeropad(int nbvox[])
  {
    size_zeropad[0]=nbvox[0];
    size_zeropad[1]=nbvox[1];
    size_zeropad[2]=nbvox[2];
  }



  void set_TV(double  alpha, int Niter)
  {
    Niter_TV=Niter;
    alpha_TV=alpha;
  }


  void one_iter_without_selection(Image &,  const std::vector<int>&,  const std::vector<int>&, const std::vector<int>&, const std::vector<double>& ,Image &);

  void one_iter_without_selection_pid(Image &,  std::vector<int>&,  const std::vector<int>&, const std::vector<int>&, const std::vector<double>& ,Image &);

  bool run(const char *results_file, const std::vector<int>& ,const std::vector<int>& , const std::vector<int>& ,const std::vector<double>& , Image &);

  void one_iter_without_selection_psf(Image &,  const std::vector<int>&,  const std::vector<int>&, const std::vector<int>&, const std::vector<double>& ,Image &,Image & filter);

  bool run_psf(const char *results_file, const std::vector<int>& ,const std::vector<int>& , const std::vector<int>& ,const std::vector<double>& , Image &,Image & filter);

  void projection(Image &, const std::vector<int>&, const std::vector<int>&, std::vector<double>&);
  void forwardprojection(const std::vector<double>&, std::vector<double>,const std::vector<int>&, const std::vector<int>&, std::vector<double>&);

};
