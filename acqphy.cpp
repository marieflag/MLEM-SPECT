
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <vector>

#include "acqphy.h"



int ParamReco::readParamInput(int argc, char * argv[], bgeot::md_param & PARAM)
{
  char *extension = strrchr( argv[1], '.');
  char * pwd = strrchr( argv[1], '/');
  run_name = argv[1];
  
  if (extension ==NULL || strncmp(extension, ".m",2)!=0)
    {
      std::cout<<"extension error, expected .m file. \n";
      return -1;
    }
  if (pwd!=NULL)
    {
      run_name.erase(0,pwd-argv[1]+1);
      run_name=run_name.substr(0,extension-pwd-1);
    }
  else
    {
      run_name=run_name.substr(0,extension-argv[1]);
    }
  std::cout<< "\n~~~~~~ * Run name and data *~~~~~~~\n\n";
  std::cout <<"Run name: " <<run_name<< '\n';

  data_file = PARAM.string_value("DATA_FILE","Data file ");
  std::cout << "Data from " << data_file << '\n';

  simu_file = PARAM.string_value("SIMU_FILE", "Simulation geometry file");
  std::cout << "Simulated system matrix is being read from: " <<simu_file<<"\n";

  return 0;
}

int ParamReco::readParamOutput(bgeot::md_param & PARAM)
{
 
  std::cout<< "\n~~~~~~ * Results *~~~~~~~\n\n";
  results_path = PARAM.string_value("RESULTS_DIR","directory with the results ");
  std::cout<< "The results will be put in "<<results_path <<'\n';

  struct stat buf;
  int ret;
  ret = stat (results_path.c_str(), &buf);
  if (ret == -1 && errno == ENOENT) 
    {
      ret = mkdir(results_path.c_str(), 0755);
      if (ret == -1) 
	{
	  std::cout<< "Could't create directory "<< results_path << '\n';
	  return -1;
	}
    } 
  events_file = results_path+run_name;

  return 0;
}

int ParamReco::readDataType(bgeot::md_param & PARAM)
{
  int err;
  std::string string_content = PARAM.string_value("DATA_TYPE","Data type ");
  if (string_content.compare("Gate") ==0)
    df=Gate;
  return 0;
}



int  ParamReco::readParamCamera(bgeot::md_param & PARAM, Head & head)
{
  std::cout<<"\n\n~~~~~* Configuration of the  DC SPECT geometry*~~~~~\n";
  std::vector<bgeot::md_param::param_value> tmp ;

  double psize= PARAM.real_value("PINHOLE_SIZE", " pinhole size ");

  head.set_psize(psize);


  tmp= PARAM.array_value("DET_SIZE", " dimensions of the detector ");
  head.set_det_dim(tmp[0].real(),tmp[1].real(),tmp[2].real());
  tmp= PARAM.array_value("DET_VOXELS", " nb detector units ");
  head.set_det_vox(int(tmp[0].real()),int(tmp[1].real()),int(tmp[2].real()));
  head.set_det_corner();

  head.set_det_vox_size();
  int nb_cameras = PARAM.int_value("NB_CAMERAS", "Number of heads in the system");
  head.set_nb_cameras(nb_cameras);
  char parameter_name[40] ;


  std::cout<<"Size of the pinhole: "<< head.get_pin_dim() << std::endl;

  std::cout<<"Size of the detector: "<< head.get_det_dim()<< std::endl;


  return 0;
}


int ParamReco::readParamVolume(bgeot::md_param & PARAM)
{
  std::cout<< "\n~~~~~~ * VOLUME *~~~~~~~\n";
  std::vector<bgeot::md_param::param_value> tmp ;
  tmp = PARAM.array_value("VOLUME_DIMENSIONS", "dimensions of the volume");
  Point volume_length(tmp[0].real(),tmp[1].real(),tmp[2].real());
  std::cout << "Volume dimensions: " << volume_length << std::endl;
  tmp= PARAM.array_value("VOXELS", "nb voxels of the volume");
  std::cout << "Voxels in the volume: [ " ;
  for (int i=0; i<3; i++) 
    {
      nb_vox[i]=int(tmp[i].real());
      std::cout<<nb_vox[i] << " ";
    }
  std::cout << "]\n";
  int total_nb_voxels=nb_vox[0]*nb_vox[1]*nb_vox[2]; 
  voxel_length = volume_length/nb_vox ;
  std::cout << "Voxel length: "<< voxel_length<< std::endl;

  tmp= PARAM.array_value("VOLUME_CENTRE", "centre of the volume");
  Point frame_center(tmp[0].real(),tmp[1].real(),tmp[2].real());
  std::cout << "Volume centre: " << frame_center <<std::endl;
  corner = frame_center-volume_length/2.0;
  std::cout << "Corner: " << corner << std::endl;


    tmp= PARAM.array_value("VOXELS_FILTER", "nb voxels of the filter");
     for (int i=0; i<3; i++)
       {
         size_filter[i]=int(tmp[i].real());
        // std::cout<<size_filter[i] << " ";
       }

  return 0;
}


int ParamMLEM::readParamPSF(bgeot::md_param & PARAM)
{
	  std::string string_content = PARAM.string_value("PSF","PSF is ON or OFF ");
	  psf_flag = string_content.compare("OFF");
	  if (psf_flag)
	    std::cout << "Analytical calculation of system matrix is ON"<< std::endl;
	  else
	    std::cout << "Analytical calculation of system matrix is OFF"<< std::endl;
      return 0;
}

int ParamMLEM::readParamPSF_deconv(bgeot::md_param & PARAM)
{
    std::string string_content = PARAM.string_value("PSF_DECONV","RL deconvolution is ON or OFF ");
    psf_deconv_flag = string_content.compare("OFF");
    if (psf_deconv_flag)
        std::cout << "RL deconvolution is ON"<< std::endl;
    else
        std::cout << "RL deconvolution is OFF"<< std::endl;
    return 0;
}
int  ParamMLEM::readParamMLEM(bgeot::md_param & PARAM, Head & head)
{
  std::cout<< "\n~~~~~~ * ITERATIONS *~~~~~~~\n";
  iter_first = PARAM.int_value("FIRST_ITERATION", "first iteration");
  if (iter_first<0) 
    {
      std::cout << "FIRST_ITERATION should not be < 0; automatically set to 0\n";
      iter_first=0;
    }
  std::cout << "First iteration: " << iter_first << "\n";
  iter_nb = PARAM.int_value("ITERATIONS", "number of iterations");
  std::cout << "Number of iterations: " << iter_nb << "\n";

  alpha_TV = PARAM.real_value("ALPHA_TV", "TV regularization");
  Niter_TV = PARAM.int_value("TV_ITERATION", "number of TV iteration");


  return 0;
}



int  ParamMLEM::readParam(int argc, char * argv[], Head & head)
{
  int err;
  bgeot::md_param PARAM;
  PARAM.read_command_line(argc, argv);
  err = readParamInput(argc, argv, PARAM);
  if (err) return -1;
  err = readParamOutput(PARAM);
  if (err) return -1;
  err = readDataType(PARAM);
  err = readParamVolume(PARAM);
  err = readParamCamera(PARAM, head);  //to be done before readParamMLEM
  err = readParamMLEM(PARAM, head);
  err = readParamPSF(PARAM);
  err = readParamPSF_deconv(PARAM);

  return 0;
}


