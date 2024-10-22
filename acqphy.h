#pragma once

#include "bgeot_ftool.h"
#include "gmm_std.h"
#include "image.h"



typedef enum { Gate} data_format;

class Head
{

 private:
  Point pin_dim;
  Point pin_centre;
  Point det_dim;
  Point det_centre;
  Point det_corner;
  Point det_vox_size;
  Point det_vox;
  int nb_cameras ;
  double psize;



 public:

 Head() : nb_cameras(0) {}

  // initialise camera parameters 

   void set_det_dim(double length, double width, double thick) { det_dim.setPoint(length, width, thick);}
  void set_psize(double i){psize = i;}
  void set_det_corner() { det_corner = det_centre-det_dim/2;}
  void set_det_vox(int length, int width, int thick) { det_vox.setPoint(length, width, thick);}
  void set_nb_cameras(int i){nb_cameras = i;}
  void set_det_vox_size() { det_vox_size.setPoint(det_dim.x/det_vox.x, det_dim.y/det_vox.y,det_dim.z/det_vox.z);}

  void add_det(double centre) { det_centre.setPoint( 0,0, centre); }



  const int get_nb_cameras() {return nb_cameras ; }
  const double get_psize(){return psize;}


 
  const Point& get_det_dim() {return det_dim; }
  const Point& get_pin_dim() { return pin_dim;}
  const Point& get_det_vox() {return det_vox;}
  const Point& get_det_centre() const {return det_centre;}
  const Point& get_pin_centre() {return pin_centre;}
  const double get_det_pos() {return det_centre.z;}
  const Point& get_det_corner() {return det_corner;}



  const Point& get_det_vox_size(){ return det_vox_size;}


};





class ParamReco
{
 protected:

  std::string run_name;
  std::string data_file;//projection *bin
  std::string simu_file;//system matrix *bin

  data_format df;
  std::string results_path; 
  std::string events_file;

  Point corner;
  int nb_vox[3];
  int size_filter[3];
  Point voxel_length;


 public:

  ParamReco() {}


  std::string get_run_name() { return run_name;}
  std::string get_data_file() { return data_file;}
  std::string get_simu_file() { return simu_file;}

  data_format get_df() { return df;}
  std::string get_results_path() { return results_path;}
  std::string get_events_file() { return events_file;}


  // the volume
  const Point & get_corner() { return corner ;}
  void get_nb_vox(int v[]) { v[0]=nb_vox[0]; v[1]=nb_vox[1]; v[2]=nb_vox[2]; }
  void get_nb_filter(int v[]) { v[0]=size_filter[0]; v[1]=size_filter[1]; v[2]=size_filter[2]; }
  void get_nb_zeropad(int v[]){v[0]=size_filter[0]+nb_vox[0]-1; v[1]=size_filter[1]+nb_vox[1]-1; v[2]=size_filter[2]+nb_vox[2]-1;}
  const Point & get_voxel_length() { return voxel_length;}


  // read parameters
  int readParamInput(int argc, char * argv[], bgeot::md_param & PARAM);
  int readParamOutput(bgeot::md_param & PARAM);
  int readDataType(bgeot::md_param & PARAM);
  int readParamCamera(bgeot::md_param & PARAM, Head & head);
  int readParamVolume(bgeot::md_param & PARAM);

};


class ParamMLEM : public ParamReco
{
 protected:

  bool psf_flag;//analytic calculation of system matrix is ON or OFF
  bool psf_deconv_flag;//Rl deconv is ON or OFF
  int iter_first;
  int iter_nb;
  double alpha_TV;
  int Niter_TV;
  std::string sens_file;


  int data_type;
 public:
  ParamMLEM() {}
  int get_data_type(){return data_type;}
  bool get_psf_flag() {return psf_flag;}
  bool get_psf_deconv_flag() {return psf_deconv_flag;}
  int get_iter_first() {return iter_first ;}
  int get_iter_nb() {return iter_nb;}
  int readParamMLEM(bgeot::md_param & PARAM, Head & head);
  int readParam(int argc, char * argv[], Head & head);
  int readParamPSF(bgeot::md_param & PARAM);
  int readParamPSF_deconv(bgeot::md_param & PARAM);
  
  double get_alpha_TV(){return alpha_TV;}
  int get_Niter_TV(){return Niter_TV;}


};
