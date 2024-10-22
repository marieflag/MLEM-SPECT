#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>


#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <unistd.h>

#include "image.h"
#include "acqphy.h"
#include "mlem.h"

#include <omp.h>

int main(int argc, char *argv[])
{



	 //recon parameters------------------------------------------------------
	  ParamMLEM list_param;
	  if (argc<2)
	    {
	      std::cout<<"missing argument, needs parameters .m file.\n";
	      return -1;
	    }




	  Head head;
	  int err =  list_param.readParam(argc, argv, head);

	  if (err) return -2;



//open & read binary file simulated from UAMC

	std::string geo_file = list_param.get_simu_file().c_str();
    std::ostringstream os;
    os << geo_file;
    std::string geofilename = os.str();
    std::ifstream inputFile(geofilename, std::ios::binary);

    std::cout << "\n\n~~~~~ * Reading simulated system matrix * ~~~~~\n\n";

	    if (!inputFile) {
	        std::cerr << "Unable to open file for reading" << std::endl;
	        return 1;
	    }

	    std::vector<uint16_t> z;
	    uint16_t value;
	    while (inputFile.read(reinterpret_cast<char*>(&value), sizeof(uint16_t))) {
	        z.push_back(value);
	    }
	    inputFile.close();

	    // Reshape the vector into a 2D vector
	    size_t numRows = z.size() / 4;
	    std::vector<std::vector<int32_t>> blah(numRows, std::vector<int32_t>(3));
           // std::vector<std::vector<uint8_t>> blah(numRows, std::vector<uint8_t>(3));	   
 for (size_t i = 0; i < numRows; ++i) {
	        blah[i][0] = z[i * 4 + 1] * 65536 + z[i * 4];
	        blah[i][1] = z[i * 4 + 2];
	        blah[i][2] = z[i * 4 + 3];
	    }



	    std::cout << "number of non-zeros in system matrix: "<< numRows<<std::endl;

	    std::string filename1 ="PID_int32.bin";
	    std::string filename2 ="VID_int32.bin";
	    std::string filename3 ="counts_int32.bin";

	    std::ofstream vidFile(filename2, std::ios::binary);
	        if (!vidFile) {
	            std::cerr << "Unable to open file for writing" << std::endl;
	            return 1;
	        }
	        for (const auto& row : blah) {
	            int32_t val = static_cast<int32_t>(row[0]);
	            vidFile.write(reinterpret_cast<char*>(&val), sizeof(int32_t));
                 /*   uint8_t val = static_cast<uint8_t>(row[0]);
                    vidFile.write(reinterpret_cast<char*>(&val), sizeof(uint8_t));*/
	        }
	        vidFile.close();

	        std::ofstream pidFile(filename1, std::ios::binary);
	        if (!pidFile) {
	            std::cerr << "Unable to open file for writing" << std::endl;
	            return 1;
	        }
	        for (const auto& row : blah) {
	            int32_t val = static_cast<int32_t>(row[1]);
	            pidFile.write(reinterpret_cast<char*>(&val), sizeof(int32_t));
                   /* uint8_t val = static_cast<uint8_t>(row[1]);
                    pidFile.write(reinterpret_cast<char*>(&val), sizeof(uint8_t));*/	       
 }
	        pidFile.close();

	        std::ofstream countsFile(filename3, std::ios::binary);
	        if (!countsFile) {
	            std::cerr << "Unable to open file for writing" << std::endl;
	            return 1;
	        }
	        for (const auto& row : blah) {
	            int32_t val = static_cast<int32_t>(row[2]);
	            countsFile.write(reinterpret_cast<char*>(&val), sizeof(int32_t));
//                    uint8_t val = static_cast<uint8_t>(row[2]);
//                    countsFile.write(reinterpret_cast<char*>(&val), sizeof(uint8_t));
	        }
	        countsFile.close();
z.clear();
z.shrink_to_fit();
blah.clear();
blah.shrink_to_fit();


  std::cout << "\n\n~~~~~ * counts, index of voxel and pixel are saved in: * ~~~~~\n\n";
  std::cout << filename1<<std::endl;
  std::cout << filename2<<std::endl;
  std::cout << filename3<<std::endl;

  // Prepare MLEM ------------------------------------------------------

  clock_t t;
  t = clock();


  double PI = 3.14;
  int nb_vox[3];
  int size_filter[3];
  int size_zeropad[3];
  list_param.get_nb_vox(nb_vox);
  list_param.get_nb_filter(size_filter);
  list_param.get_nb_zeropad(size_zeropad);

  Point detsize = head.get_det_vox();
  //int param_size = numRows;
  size_t param_size = numRows;

  /*std::vector<int> PID_15 (param_size*sizeof(int),1);
  std::vector<int> VID_15 (param_size*sizeof(int),1);


    std::vector<int> counts_15 (param_size*sizeof(int),1);*/
    std::vector<int> PID_15 (param_size,1);
     std::vector<int> VID_15 (param_size,1);


       std::vector<int> counts_15 (param_size,1);
    std::vector<double> projd (head.get_nb_cameras()*detsize.x*detsize.y*sizeof(double),1);
    

    std::cout << "\n\n~~~~~ * reading counts, index of voxels and pixels from processed files * ~~~~~\n\n";
    

    std::ifstream F1;
    F1.open(filename1.c_str(), std::ios::in|std::ios::binary);
    F1.read((char *)(&(PID_15[0])), param_size*sizeof(int));


    std::ifstream F2;
    F2.open(filename2.c_str(), std::ios::in|std::ios::binary);
    F2.read((char *)(&(VID_15[0])), param_size*sizeof(int));


    std::ifstream F3;
    F3.open(filename3.c_str(), std::ios::in|std::ios::binary);
    F3.read((char *)(&(counts_15[0])), param_size*sizeof(int));

    std::cout << "\n\n~~~~~ * reading projection * ~~~~~\n\n";

    std::string data_file = list_param.get_data_file();
    std::ostringstream oss;
    oss << data_file;
    std::string filename4 = oss.str();
    std::ifstream F4;
    F4.open(filename4.c_str(), std::ios::in|std::ios::binary);
    F4.read((char *)(&(projd[0])), head.get_nb_cameras()*detsize.x*detsize.y*sizeof(double));
    

  int sample_id=0;
  char name_results_file[500];



  std::cout << "\n\n~~~~~ * Initialise Recon * ~~~~~\n\n";


  MLEM mlem;
  mlem.set_iterations(list_param.get_iter_nb(),list_param.get_iter_first());
  mlem.set_volume(nb_vox,list_param.get_voxel_length(),list_param.get_corner());
  mlem.set_TV(list_param.get_alpha_TV(), list_param.get_Niter_TV());

  mlem.set_zeropad(size_zeropad);
  Image sens(nb_vox,list_param.get_voxel_length(),list_param.get_corner());

std::cout <<"\n\n~~~~~ * translating FOV coordinates from UAMC * ~~~~~~~~\n\n";

for (int j=0;j<PID_15.size();j++){

        int z = int(VID_15[j]/(sens.DimInVoxels[0]*sens.DimInVoxels[1]));
        int x = int((VID_15[j]-z*sens.DimInVoxels[0]*sens.DimInVoxels[1])/sens.DimInVoxels[1]);
        int y = sens.DimInVoxels[1]-int(VID_15[j] - x*sens.DimInVoxels[1]-z*sens.DimInVoxels[0]*sens.DimInVoxels[1])-1;
        VID_15[j] =int((sens.DimInVoxels[0]*sens.DimInVoxels[1])*z + sens.DimInVoxels[0]*y + x);
    }
 
    //prepare psf for RL-deconvolution---------------------------------------------
    
    Image filter(size_filter,list_param.get_voxel_length(), list_param.get_corner());





    //run MLEM---------------------------------------------
     std::cout << "\n\n~~~~~ * Start MLEM * ~~~~~\n\n";

     sprintf(name_results_file, "%s.sample%d",list_param.get_events_file().c_str(),sample_id);


      if(list_param.get_psf_deconv_flag()){
    	  std::cout << "\n\n~~~~~ * RL-deconvolution is called in MLEM * ~~~~~\n\n";
            std::string filename ="image_point.bin";
            std::ifstream F;
            F.open(filename.c_str(), std::ios::in|std::ios::binary);
            if (F.is_open()){
                std::cout <<"opening filter file "<< filename<<". \n" ;

                F.read((char *)(&(filter.Value[0])), size_filter[0]*size_filter[1]*size_filter[2]*sizeof(double));

            }

            if (!F){
                std::cout <<"error opening filter file"<< ". \n" ;}
           if (!mlem.run_psf(name_results_file, PID_15, VID_15, counts_15, projd, sens,filter)) return 1;}

      else{
    	  if (!mlem.run(name_results_file, PID_15, VID_15, counts_15, projd, sens));
      return 1;
      }

  return 0;
}


