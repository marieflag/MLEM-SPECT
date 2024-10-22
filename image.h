#pragma once

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>


class Point
{

 public:

  double x,y,z;

  Point(double xx=0, double yy=0, double zz=0);
  Point(Point const& P);


  void setPoint(double xx, double yy, double zz);

  void display(std::ostream& flux) const;

  void operator=(Point const& P);
  Point& operator+=(const Point& P);
  Point& operator-=(const Point& P);
  Point& operator*=(const int* C);
  Point& operator*=(const double C);
  Point& operator/=(const int* C);
  Point& operator/=(const double C);


  double get_max() {return std::max(x,std::max(y,z)) ;}
  Point get_local_coord( const Point &O, const Point &Ox, const Point &Oy, const Point &Oz) const;
  Point get_absolute_coord(const Point &O, const Point &Ox, const Point &Oy, const Point &Oz) const;
  Point get_spherical_coord(const Point &Ox, const Point &Oy, const Point &Oz) const;
  double norm2() const  {return sqrt(x*x+y*y+z*z);}


  bool inVolume(const Point& corner, const Point& volume_dim) const;

  bool inDetector(const Point& dim_detector, const Point& center_detector) const ;
};

std::ostream &operator<<( std::ostream& flux, Point const& P);
Point operator+(Point const &P1, Point const &P2);
Point operator-(Point const &P1, Point const &P2);
double operator*(const Point& P1, const Point& P2);  // dot product
Point operator*(const int* C, Point const& P);  // point by point product
Point operator*(double C, Point const& P);  // point by point product
Point operator/(Point const& P, const int* C);  // point by point division
Point operator/(Point const& P, double C);  // point by point division





class Image
{

 public:

  int DimInVoxels[3]; // image dimensions in voxels
  Point DimInCm; // image dimensions in cm
  Point VoxelSize; // voxel dimensions in cm
  Point Corner; // coordinates of the corner bottom left, front
  int NbVoxels;
  std::vector<double> Value;



  Image() : NbVoxels(0) {}
  Image(const int* image_size, const Point& voxel_size, const Point& corner);


  bool index_1Dto3D(int index_voxel, int& i, int& j, int& k) const;
  int index_3Dto1D(int i,int j, int k) const;
  bool coord2index_3D(Point const& P, int& i, int& j, int& k) const;
  int coord2index_1D(Point const& P) const;
	
  void initialize(double value)
  {
    Value.assign( NbVoxels, value);
  }

  void threshold(void);
  void filter(void);
  void multiply(const std::vector<double> mult, const std::vector<double> divisor);
  bool setIntensity(const std::vector<int> index_voxel, const std::vector<double> value);
  double maximum() const;
  double maxabs() const;
  double minimum() const;

  int readFile(const char *image_file);
  const int writeFile(const char *image_file);
  bool readFromFile(std::ifstream&  image_file);  // returns eof
  const void writeToFile(std::ofstream& image_file);
  void arraypad(Image &im_pad, const Image &im);
 // void calcul_filter(Image &filter, int a, int b, int c, int number, const std::vector<double>& A_psf, const std::vector<double>& B, const std::vector<double>& D,const std::vector<double>& E, const std::vector<double>& sigmax_left, const std::vector<double>& sigmax_right, const std::vector<double>& sigmay_left, const std::vector<double>& sigmay_right, const std::vector<double>& sigmaz_left, const std::vector<double>& sigmaz_right);
  void convolution(Image &im_conv , const Image &im_pad, Image &filter);
  void correlation(Image &im_corr , const Image &im_pad, Image &filter);
  void TV_regularization(Image &im, Image &im_TV, const Image &sens,Image &VolumeTV, Image &VolumeTV_pad, double alpha_TV, int Niter_TV, double d);
  void TV_regularization_ca(double log_likelihood, Image &im, Image &im_TV, const Image &sens, Image &VolumeTV, Image &divp, Image &lastdivp, Image &TV_im, Image &p3, Image &p2 ,Image &p1, double alpha_TV, int Niter_TV);
  void TV_regularization_vm(double log_likelihood, Image &im,const Image &sens, Image &VolumeTV, Image &divp, Image &lastdivp, Image &TV_im, Image &p3, Image &p2 ,Image &p1,Image &z3, Image &z2 ,Image &z1, double alpha_TV, int Niter_TV, double EPSILON);
  void FISTA(Image &im_fista, Image &im_0,Image &im,Image &t1,Image &t2);
};


double operator*(const Image& I, const Image& J);  //dot product of J by I


