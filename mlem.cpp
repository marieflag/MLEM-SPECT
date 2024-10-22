
#include <math.h>
#include <stdlib.h>
#include "mlem.h"
#include <algorithm>
#include <numeric>

#include <omp.h>

#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = omp_orig)



void MLEM::projection(Image & lambda, const std::vector<int>& VID_15, const std::vector<int>& counts_15, std::vector<double>& tmp)
{
	double start_time = omp_get_wtime();


	#pragma omp parallel for
	for (int j = 0; j < VID_15.size(); j++) {
		int vid = VID_15[j];
	    tmp[j] = lambda.Value[vid] * counts_15[j];

	}


	double end_time = omp_get_wtime();
	std::cout << "projection time " << end_time-start_time <<" seconds...\n\n";

}

void MLEM::forwardprojection(const std::vector<double>& ProjD, std::vector<double> Proj,const std::vector<int>& PID_15, const std::vector<int>& counts_15, std::vector<double>& tmp)
{
	double start_time = omp_get_wtime();


        #pragma omp parallel for
	    for (int j=0;j<PID_15.size();j++){
	    	int pid = PID_15[j];
	  	tmp[j] =  counts_15[j]*ProjD[pid]/Proj[pid];
	      }
	    double end_time = omp_get_wtime();
    	std::cout << "forward projection time " << end_time-start_time <<" seconds...\n\n";

}

void MLEM::one_iter_without_selection(Image & lambda, const std::vector<int>& PID_15, const std::vector<int>& VID_15, const std::vector<int>& counts_15, const std::vector<double>& ProjD, Image & sens)

// Iterate with already selected good events
{
	  Image VolumeTV(nb_vox,voxel_length, corner);

	  Image p1(nb_vox,voxel_length, corner);
	  Image p2(nb_vox,voxel_length, corner);
	  Image p3(nb_vox,voxel_length, corner);
	  Image z1(nb_vox,voxel_length, corner);
	  Image z2(nb_vox,voxel_length, corner);
	  Image z3(nb_vox,voxel_length, corner);

	  Image TV_im(nb_vox,voxel_length, corner);
	  Image divp(nb_vox,voxel_length, corner);
	  Image lastdivp(nb_vox,voxel_length, corner);
	  Image ForwardProj(nb_vox,voxel_length, corner);

	  std::vector<double> & Val = ForwardProj.Value;

    std::vector<double> Proj (ProjD.size(),0);

    std::vector<double> tmp (VID_15.size(),1);

    Point el, pinholeloc;
    Point OjV1, OjV1_init= lambda.Corner;


    projection(lambda,VID_15,counts_15, tmp);

   #pragma omp parallel for reduction(vec_double_plus:Proj)
   for (int j=0;j<PID_15.size();j++){
   	int pid = PID_15[j];
    Proj[pid] += tmp[j];
   }

   forwardprojection(ProjD, Proj,PID_15, counts_15, tmp);


   #pragma omp parallel for reduction(vec_double_plus:Val)
   for (int j=0;j<PID_15.size();j++){
    int vid = VID_15[j];
   	Val[vid] += tmp[j];
   }


    #pragma omp parallel for
    for (int k=0;k<lambda.NbVoxels;k++){
    	lambda.Value[k] = ForwardProj.Value[k]*lambda.Value[k]/sens.Value[k];
    }


/*   double EPSILON= 0.001;
    double alpha_TV_last, alpha_TV, log_likelihood;
    int Niter_TV;
    for (int k=0; k<lambda.NbVoxels; k++){
        		      lambda.Value[k]=lambda.Value[k]+EPSILON;}
                 alpha_TV_last = alpha_TV;

                 image.TV_regularization_vm(log_likelihood,lambda,sens,VolumeTV,divp,lastdivp,TV_im,p3,p2,p1,z3,z2,z1,alpha_TV_last, Niter_TV, EPSILON);
                 for (int k=0; k<lambda.NbVoxels; k++){
                 lambda.Value[k]=lambda.Value[k]-EPSILON;}*/


}


void MLEM::one_iter_without_selection_psf(Image & lambda, const std::vector<int>& PID_15, const std::vector<int>& VID_15, const std::vector<int>& counts_15, const std::vector<double>& ProjD, Image & sens,Image & filter)
// Iterate with already selected good events
{


    Image ForwardProj(nb_vox,voxel_length, corner);

    std::vector<double> Proj (ProjD.size(), 0);
    Image lambda_conv(nb_vox,voxel_length, corner);
    Image ForwardProj_conv(nb_vox,voxel_length, corner);


    Image lambda_pad(size_zeropad,voxel_length, corner);
    Image ForwardProj_pad(size_zeropad,voxel_length, corner);

  image.arraypad(lambda_pad,lambda);
  image.convolution(lambda_conv,lambda_pad,filter);


   for (int j=0;j<PID_15.size();j++){
        int pid = PID_15[j];
        Proj[pid] += lambda.Value[VID_15[j]]*counts_15[j];

    }


    for (int j=0;j<PID_15.size();j++){
                int pid = PID_15[j];
	ForwardProj.Value[VID_15[j]] +=  counts_15[j]*ProjD[pid]/Proj[pid];
    }

    image.arraypad(ForwardProj_pad,ForwardProj);

    image.correlation(ForwardProj_conv,ForwardProj_pad,filter);


    for (int k=0;k<ForwardProj.NbVoxels;k++){
    	lambda.Value[k] = ForwardProj_conv.Value[k]*lambda.Value[k]/sens.Value[k];
    }



}


bool MLEM::run( const char *results_file, const std::vector<int>& PID_15, const std::vector<int>& VID_15, const std::vector<int>& counts_15, const std::vector<double>& projd, Image &sens)
// run all iterations, without hodoscope
{
  Image lambda(nb_vox,voxel_length, corner);
  lambda.initialize(1);

  char results_file_iter[500];
  char results_file_sens[500];

  clock_t t = clock();




  for (int j=0;j<PID_15.size();j++){

	  int z = int(VID_15[j]/(lambda.DimInVoxels[1]*lambda.DimInVoxels[0]));
	          int x = int((VID_15[j]-z*lambda.DimInVoxels[1]*lambda.DimInVoxels[0])/lambda.DimInVoxels[1]);
	          int y = lambda.DimInVoxels[1]-int(VID_15[j] - x*lambda.DimInVoxels[1]-z*lambda.DimInVoxels[0]*lambda.DimInVoxels[1])-1;

        //  double c = counts_15[j];
      sens.Value[sens.index_3Dto1D(x,y,z)] += counts_15[j];
  }

  t = clock()-t;
  std::cout<< "Time spent in sensitivity matrix calculation: "<< ((float)t)/CLOCKS_PER_SEC<<" seconds \n\n";

  sprintf(results_file_sens, "%s.sens.bin",results_file);
  if (sens.writeFile(results_file_sens)) return false;

	  for (int iter = first_iter; iter <= first_iter+nb_iter; iter++)
    {

      std::cout << "Iteration " << iter <<" ...\n\n";
      clock_t t = clock();
      

      one_iter_without_selection(lambda, PID_15, VID_15, counts_15, projd, sens);
      sprintf(results_file_iter, "%s.iter%d.bin",results_file, iter);

    	 if (lambda.writeFile(results_file_iter)) return false;
    	            t = clock()-t;
    	           // std::cout<< "Time spent in iteration: "<< ((float)t)/CLOCKS_PER_SEC<<" seconds \n\n";

    }


  return true;
}


bool MLEM::run_psf( const char *results_file, const std::vector<int>& PID_15, const std::vector<int>& VID_15, const std::vector<int>& counts_15, const std::vector<double>& projd, Image &sens,Image & filter)
// run all iterations, without hodoscope
{
  Image lambda(nb_vox,voxel_length, corner);
  lambda.initialize(1);

  char results_file_iter[500];
  char results_file_sens[500];


  sprintf(results_file_sens, "%s.sens.bin",results_file);
  if (sens.writeFile(results_file_sens)) return false;
  sprintf(results_file_iter, "%s.iter%d.bin",results_file, 0);
  lambda.writeFile(results_file_iter);

	  for (int iter = 1; iter <= first_iter+nb_iter; iter++)
    {

      std::cout << "Iteration " << iter <<" ...\n\n";
      clock_t t = clock();


      one_iter_without_selection_psf(lambda, PID_15, VID_15, counts_15, projd, sens,filter);

      sprintf(results_file_iter, "%s.iter%d.bin",results_file, iter);

    	 if (lambda.writeFile(results_file_iter)) return false;
    	            t = clock()-t;
    	            std::cout<< "Time spent in iteration: "<< ((float)t)/CLOCKS_PER_SEC<<" seconds \n\n";

    }


  return true;
}
