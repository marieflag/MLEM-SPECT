#include <iostream>
#include <fstream>
#include <math.h>
#include "image.h"



Point::Point(double xx, double yy, double zz) : x(xx), y(yy), z(zz)
{
}

Point::Point(Point const& P) 
{
  x = P.x;
  y = P.y;
  z = P.z;
}


void Point::setPoint(double xx, double yy, double zz)
{
  x = xx;
  y = yy;
  z = zz;
}

void Point::operator=(Point const& P)
{
  x = P.x;
  y = P.y;
  z = P.z;
}

Point& Point::operator+=(const Point& P)
{
  x += P.x;
  y += P.y;
  z += P.z;
  return *this;
}

Point& Point::operator-=(const Point& P)
{
  x -= P.x;
  y -= P.y;
  z -= P.z;
  return *this;
}


Point& Point::operator*=(const int* C)
{
  x *= C[0];
  y *= C[1];
  z *= C[2];
  return *this;
}

Point& Point::operator*=(const double C)
{
  x *= C;
  y *= C;
  z *= C;
  return *this;
}

Point& Point::operator/=(const int* C)
{
  x /= C[0];
  y /= C[1];
  z /= C[2];
  return *this;
}

Point& Point::operator/=(const double C)
{
  x /= C;
  y /= C;
  z /= C;
  return *this;
}

void Point::display(std::ostream& flux) const
{
  flux << x<< ' '<< y << ' '<< z;
}

Point Point::get_local_coord(const Point &O, const Point &Ox, const Point &Oy, const Point &Oz) const
  {
    Point translated((*this)-O);
    Point P_loc(translated*Ox, translated*Oy, translated*Oz);
    return P_loc;
  }

Point Point::get_absolute_coord(const Point &O, const Point &Ox, const Point &Oy, const Point &Oz) const
{
  Point P_abs = O+x*Ox+y*Oy+z*Oz;
  return P_abs;
}

Point Point::get_spherical_coord(const Point &Ox, const Point &Oy, const Point &Oz) const
{
  Point P;
  Point P_loc((*this)*Ox, (*this)*Oy, (*this)*Oz);
  P.x = sqrt(x*x+y*y+z*z);
  P.y = acos(P_loc.z/P.x);
  P.z = atan2(P_loc.y,P_loc.x);
  return P;
}

bool Point::inVolume(const Point& corner, const Point& volume_dim) const
{
  bool x_ok=(x>=corner.x) && (x<=corner.x+volume_dim.x);
  bool y_ok=(y>=corner.y) && (y<=corner.x+volume_dim.y);
  bool z_ok=(z>=corner.z) && (z<=corner.z+volume_dim.z);
  return (x_ok && y_ok && z_ok);
}

bool Point::inDetector(const Point& dim_detector, const Point& center_detector) const
{
  bool x_ok=fabs(x-center_detector.x)<=dim_detector.x/2.0+0.0001;
  bool y_ok=fabs(y-center_detector.y)<=dim_detector.y/2.0+0.0001;
  bool z_ok=fabs(z-center_detector.z)<=dim_detector.z/2.0+0.0001;
  return (x_ok && y_ok && z_ok);
}


std::ostream &operator<<(std::ostream&flux, const Point & P)
{
  P.display(flux);
  return flux;
}

Point operator+(const Point &P1, const Point &P2)
{
  Point copy(P1);
  copy += P2;
  return copy;
}

Point operator-(const Point &P1, const Point &P2)
{
  Point copy(P1);
  copy -= P2;
  return copy;
}

double operator*(const Point& P1, const Point& P2)
{
  return P1.x*P2.x+P1.y*P2.y+P1.z*P2.z;
}

Point operator*(const int* C, const Point & P)
{
  Point copy(P);
  copy *= C;
  return copy;
}

Point operator*(double C, const Point & P)
{
  Point copy(P);
  copy *= C;
  return copy;
}

Point operator/(const Point & P, const int* C)
{
  Point copy(P);
  copy /= C;
  return copy;
}

Point operator/(const Point & P, double C)
{
  Point copy(P);
  copy /= C;
  return copy;
}





Image::Image(const int* image_size, const Point& voxel_size, const Point& corner)
{
  for (int i=0;i<3; i++)
    {
      DimInVoxels[i]=image_size[i];
    }
  VoxelSize = voxel_size;
  Corner = corner;
  DimInCm =  DimInVoxels*VoxelSize;
  NbVoxels=DimInVoxels[0]*DimInVoxels[1]*DimInVoxels[2];
  Value.assign( NbVoxels,0.0);
}


bool Image::index_1Dto3D(int index_voxel, int& i, int& j, int& k) const
{
  if ((index_voxel<0) || (index_voxel>=NbVoxels)) return false;
  k=index_voxel / (DimInVoxels[0]*DimInVoxels[1]);
  j=index_voxel % (DimInVoxels[0]*DimInVoxels[1]);
  i= j % DimInVoxels[0];
  j = j/DimInVoxels[0];
  return true;
}

int Image::index_3Dto1D(int i,int j, int k) const
{
  if ((i<0)||(i>=DimInVoxels[0])||(j<0)||(j>=DimInVoxels[1])||(k<0)||(k>=DimInVoxels[2])) return -1;
  return (i+j*DimInVoxels[0]+k*DimInVoxels[0]*DimInVoxels[1]);
}

bool Image::coord2index_3D(Point const& P, int& i, int& j, int& k) const
{
  bool inside = P.inVolume(Corner, DimInCm);
  //((x>=Corner[0]) && (x<=Corner[0]+DimInCm[0]) && (y>=Corner[1]) && (y<=Corner[1]+DimInCm[1]) && (z>=Corner[2]) && (z<=Corner[2]+DimInCm[2]));
  i=int(std::min(int(floor((P.x-Corner.x)/VoxelSize.x)), DimInVoxels[0]));
  j=int(std::min(int(floor((P.y-Corner.y)/VoxelSize.y)), DimInVoxels[1]));
  k=int(std::min(int(floor((P.z-Corner.z)/VoxelSize.z)), DimInVoxels[2]));
  // std::cout<< i<< ' ' << j<< ' ' << k << ' ' << '\n';
  return(inside);
}

int Image::coord2index_1D(Point const& P) const
{
  int i,j,k;
  bool inside=Image::coord2index_3D(P,i,j,k);
  if (inside) 
    return Image::index_3Dto1D(i,j,k);
  else
    return -1;
}


void Image::threshold(void)
{
  double thr=maxabs()/1000;
  for (int i=0; i<NbVoxels; i++)
    if (fabs(Value[i])<thr) Value[i]=0;
}


void Image::filter(void)
{
  int i,k,l=0, step;
  Image copy( DimInVoxels, VoxelSize, Corner);
  for(k=0; k<DimInVoxels[1]*DimInVoxels[2]; k++)
    {
      copy.Value[l] = Value[l]/2+Value[l+1]/4;
      l++;
      for ( ; l<(k+1)*DimInVoxels[0]-1 ; l++)
	copy.Value[l] = Value[l-1]/4+Value[l]/2+Value[l+1]/4;
      copy.Value[l] = Value[l-1]/4+Value[l]/2;
      l++;
    }

  step= DimInVoxels[0];
  //Value.assign (copy.Value.begin(),copy.Value.end());
  for (k=0; k<DimInVoxels[2]; k++)
    for(i=0; i<DimInVoxels[0]; i++)
      {
	l=k*DimInVoxels[0]*DimInVoxels[1]+i;
	Value[l] = copy.Value[l]/2+copy.Value[l+step]/4;
	l +=step;
	for ( ; l<k*DimInVoxels[0]*DimInVoxels[1]+DimInVoxels[0]*(DimInVoxels[1]-1)+i ; l=l+step)
	  Value[l] = copy.Value[l-step]/4+copy.Value[l]/2+copy.Value[l+step]/4;
	Value[l] = copy.Value[l-step]/4+copy.Value[l]/2;
      }
  
  step=DimInVoxels[0]*DimInVoxels[1];
  for(i=0; i<DimInVoxels[0]*DimInVoxels[1]; i++)
    {
      l=i;
      copy.Value[l]=Value[l]/2+Value[l+step]/4;
      l +=step;
      for ( ; l<i+(DimInVoxels[2]-1)*DimInVoxels[0]*DimInVoxels[1]; l=l+step)
	copy.Value[l]=Value[l-step]/4+Value[l]/2+Value[l+step]/4;
      copy.Value[l]=Value[l-step]/4+Value[l]/2;
    }
  Value.assign (copy.Value.begin(),copy.Value.end());
}


void Image::multiply(const std::vector<double> mult, const std::vector<double> divisor)
{
  for(int i=0; i<NbVoxels; i++)
    Value[i] *= mult[i]/divisor[i];
}

bool Image::setIntensity(const std::vector<int> index_voxel, const std::vector<double> value)
{
  if (index_voxel.size()!=value.size()) return false;
  for(int i=0; i<index_voxel.size(); i++)
    Value[index_voxel[i]] = value[i];
  return true;
}

double Image::maximum() const
{
  double max = 0.0;
  for (int i = 0; i<NbVoxels;i++)
    if (Value[i] > max) max = Value[i];
  return max;
}


double Image::minimum() const
{
  double min = Value[0];
  for (int i = 1; i<NbVoxels;i++)
    if (Value[i] <= min) min = Value[i];
  return min;
}


double Image::maxabs() const
{
  double max = 0.0;
  for (int i = 0; i<NbVoxels;i++)
    if (fabs(Value[i]) > max) max = fabs(Value[i]);
  return max;
}

int Image::readFile(const char *image_file)
{
  std::ifstream fp;
  fp.open(image_file, std::ios::in|std::ios::binary);
  if (!fp.is_open())
    {
      std::cout<<"Unknown file "<<image_file<< '\n';
      return 1;
    }
  fp.read((char *)(&(Value[0])), NbVoxels*sizeof(double));
  if (!fp) 
    {
      std::cout << "Error in " << image_file << ": only "<< fp.gcount()<< " items could be read.\n" ;
      return 1;
    }
  fp.close();
  return 0;
}

const int Image::writeFile(const char *image_file)
{
  std::ofstream fp;
  fp.open(image_file, std::ios::out|std::ios::binary);
  if (!fp.is_open())
    {
      std::cout<<"Unknown file "<<image_file<< '\n';
      return 1;
    }
  fp.write((char*)(&(Value[0])), NbVoxels*sizeof(double));
  fp.close();
  return 0;
}


bool Image::readFromFile(std::ifstream& fp)
{
  fp.read((char *)(&(Value[0])), NbVoxels*sizeof(double));
  return fp.eof();
  /*
  if (!fp) 
    {
      std::cout << "Error in " << image_file << ": only "<< fp.gcount()<< " items could be read.\n" ;
      return 1;
    }
  */
}

const void Image::writeToFile(std::ofstream& fp)
{
  fp.write((char*)(&(Value[0])), NbVoxels*sizeof(double));
}





double operator*(const Image& I, const Image & J)
// dot product of two images
{
  double sum = 0;
  if (!(J.DimInVoxels[0]==I.DimInVoxels[0] && J.DimInVoxels[1]==I.DimInVoxels[1] && J.DimInVoxels[2]==I.DimInVoxels[2]))
    {
      std::cout << "Error in operator* for images: to calculate dot product the two vectors should have the same length. Return 0"<< std::endl;
      return 0;
    }
  for (int i = 0; i<I.NbVoxels; i++)
    sum += J.Value[i]*I.Value[i];
  return sum;
}

void Image::arraypad(Image &im_pad, const Image &im)

{

	im_pad.initialize(0);

	div_t divresult;
	std::vector<int> infb(3,0);
	infb[0] = div(im_pad.DimInVoxels[0]-im.DimInVoxels[0],2).quot;
	infb[1] = div(im_pad.DimInVoxels[1]-im.DimInVoxels[1],2).quot;
	infb[2] = div(im_pad.DimInVoxels[2]-im.DimInVoxels[2],2).quot;
	//std::cout <<"correlation function being called during the iteration"<<im_pad.DimInVoxels[2]<< '\n';
	for(int z = infb[2]; z<im_pad.DimInVoxels[2]-infb[2];z++){
			for(int y = infb[1]; y<im_pad.DimInVoxels[1]-infb[1];y++){
				for(int x = infb[0]; x<im_pad.DimInVoxels[0]-infb[0];x++){
					im_pad.Value[im_pad.index_3Dto1D(x,y,z)] = im.Value[im.index_3Dto1D(x-infb[0],y-infb[1],z-infb[2])];}}}


	for(int z = 0; z<infb[2];z++){
				for(int y = infb[1]; y<im_pad.DimInVoxels[1]-infb[1];y++){
					for(int x = infb[0]; x<im_pad.DimInVoxels[0]-infb[0];x++){
						im_pad.Value[im_pad.index_3Dto1D(x,y,z)] = im.Value[im.index_3Dto1D(x-infb[0],y-infb[1],infb[2]-z-1)];}}}

	for(int z = im_pad.DimInVoxels[2]-infb[2]; z<im_pad.DimInVoxels[2];z++){
					for(int y = infb[1]; y<im_pad.DimInVoxels[1]-infb[1];y++){
						for(int x = infb[0]; x<im_pad.DimInVoxels[0]-infb[0];x++){
							im_pad.Value[im_pad.index_3Dto1D(x,y,z)] = im.Value[im.index_3Dto1D(x-infb[0],y-infb[1],im.DimInVoxels[2]-(im_pad.DimInVoxels[2]-z))];}}}

}

void Image::convolution(Image &im_conv, const Image &im_pad, Image &filter)
{

	clock_t time;
    time = clock();

    im_conv.initialize(0);
	//std::cout <<"convolution function being called during the iteration"<< '\n';


	div_t divresult;
	std::vector<int> infb(3,0), supb(3,0);
	divresult = div (filter.DimInVoxels[0],2);
	if (divresult.rem)
		infb[0] = divresult.quot;
	else
		infb[0] = divresult.quot-1;
	supb[0] = divresult.quot;

	divresult = div (filter.DimInVoxels[1],2);
	if (divresult.rem)
		infb[1] = divresult.quot;
	else
		infb[1] = divresult.quot-1;
	supb[1] = divresult.quot;

	divresult = div (filter.DimInVoxels[2],2);
	if (divresult.rem)
		infb[2] = divresult.quot;
	else
		infb[2] = divresult.quot-1;
	supb[2] = divresult.quot;


	for(int z = infb[2]; z<im_pad.DimInVoxels[2]-infb[2];z++){
			for(int y = infb[1]; y<im_pad.DimInVoxels[1]-infb[1];y++){
				for(int x = infb[0]; x<im_pad.DimInVoxels[0]-infb[0];x++){
					int number = im_conv.index_3Dto1D(x-infb[0],y-infb[1],z-infb[2]);
	//				//std::cout << "corr maximum "<< number<< std::endl;
	//				//std::cout << "corr maximum "<<sigmax[number]<< std::endl;
	//				calcul_filter(filter,x-infb[0],y-infb[1],z-infb[2],number,A_psf,B,D,E, sigmax_left,sigmax_right, sigmay_left,sigmay_right, sigmaz_left,sigmaz_right);
	//				//parameter filter = file of parameter filter(x-infb[0],y-infb[0],z-infb[2]);
	//				//kp =aa.*exp(-0.5*(((yy-bb)./sigmax_x).^2+((xx-dd)./sigmay_y).^2+((zz-ee)./sigmaz_z).^2));
					for(int k=z-infb[2]; k<z+supb[2]+1; k++){
						int ind_z = -k+z+supb[2];
						for(int j=y-infb[1]; j<y+supb[1]+1; j++){
							int ind_y = -j+y+supb[1];
							for(int i=x-infb[0]; i<x+supb[0]+1; i++){
								im_conv.Value[im_conv.index_3Dto1D(x-infb[0],y-infb[1],z-infb[2])] += im_pad.Value[im_pad.index_3Dto1D(i,j,k)]*filter.Value[filter.index_3Dto1D(-i+x+supb[0],ind_y,ind_z)];
                                //std::cout << "corr maximum "<< filter.Value[filter.index_3Dto1D(-i+x+supb[0],ind_y,ind_z)]<< std::endl;
							}}}

       

			}
	}
	}

	time = clock() - time;
	printf ("convolution took %ld clicks (%f seconds).\n",time,((float)time)/CLOCKS_PER_SEC);
}

void Image::correlation(Image &im_corr, const Image &im_pad, Image &filter)
{
	clock_t time;
        time = clock();
	im_corr.initialize(0);


	//std::cout <<"correlation function being called during the iteration"<< '\n';


	div_t divresult;
	std::vector<int> infb(3,0), supb(3,0);
	divresult = div (filter.DimInVoxels[0],2);
	if (divresult.rem)
		infb[0] = divresult.quot;
	else
		infb[0] = divresult.quot-1;
	supb[0] = divresult.quot;

	divresult = div (filter.DimInVoxels[1],2);
	if (divresult.rem)
		infb[1] = divresult.quot;
	else
		infb[1] = divresult.quot-1;
	supb[1] = divresult.quot;

	divresult = div (filter.DimInVoxels[2],2);
	if (divresult.rem)
		infb[2] = divresult.quot;
	else
		infb[2] = divresult.quot-1;
	supb[2] = divresult.quot;



	for(int z = infb[2]; z<im_pad.DimInVoxels[2]-infb[2];z++){
		for(int y = infb[1]; y<im_pad.DimInVoxels[1]-infb[1];y++){
			for(int x = infb[0]; x<im_pad.DimInVoxels[0]-infb[0];x++){
				int number = im_corr.index_3Dto1D(x-infb[0],y-infb[1],z-infb[2]);

			//	calcul_filter(filter,x-infb[0],y-infb[1],z-infb[2],number,A_psf,B,D,E, sigmax_left,sigmax_right, sigmay_left,sigmay_right, sigmaz_left,sigmaz_right);
			//	}
				for(int k=z-infb[2]; k<z+supb[2]+1; k++){
					int ind_z = k-z+infb[2];
					for(int j=y-infb[1]; j<y+supb[1]+1; j++){
						int ind_y = j-y+infb[1];
						for(int i=x-infb[0]; i<x+supb[0]+1; i++){
							im_corr.Value[im_corr.index_3Dto1D(x-infb[0],y-infb[1],z-infb[2])] += im_pad.Value[im_pad.index_3Dto1D(i,j,k)]*filter.Value[filter.index_3Dto1D(i-x+infb[0],ind_y,ind_z)];
							
									}}}
                




			}
	}
	}

	time = clock() - time;
    printf ("Correlation took %ld clicks (%f seconds).\n",time,((float)time)/CLOCKS_PER_SEC);
}






/*void Image::calcul_filter(Image &filter, int a, int b, int c, int number, const std::vector<double>& A_psf, const std::vector<double>& B, const std::vector<double>& D,const std::vector<double>& E, const std::vector<double>& sigmax_left, const std::vector<double>& sigmax_right, const std::vector<double>& sigmay_left, const std::vector<double>& sigmay_right, const std::vector<double>& sigmaz_left, const std::vector<double>& sigmaz_right)
{

	clock_t time;
    time = clock();
	//std::cout << "volume"<<filter.DimInVoxels[0]<<std::endl;
	double sum=0.0;
	double somme;
	div_t divresult;
	std::vector<int> infb(3,0), supb(3,0);
	divresult = div (filter.DimInVoxels[0],2);
	if (divresult.rem)
		infb[0] = divresult.quot;
	else
		infb[0] = divresult.quot-1;
	supb[0] = divresult.quot;

	divresult = div (filter.DimInVoxels[1],2);
	if (divresult.rem)
		infb[1] = divresult.quot;
	else
		infb[1] = divresult.quot-1;
	supb[1] = divresult.quot;

	divresult = div (filter.DimInVoxels[2],2);
	if (divresult.rem)
		infb[2] = divresult.quot;
	else
		infb[2] = divresult.quot-1;
	supb[2] = divresult.quot;
	//std::cout << "corr maximum "<<sigmax[number]<< std::endl;
	//std::cout << "corr maximum "<< B[number]<< std::endl;

    double sigmax = sigmax_left[number];
    double sigmay = sigmay_left[number];
    double sigmaz = sigmaz_left[number];

    double rotx = sigmax_right[number];
    double rotz = sigmaz_right[number];

    double cos_z = cos(rotz);
    double sin_z = sin(rotz);
    double cos_x = cos(rotx);
    double sin_x = sin(rotx);

    
	for(int l=0; l<filter.DimInVoxels[2]; l++){
						for(int m=0; m<filter.DimInVoxels[1]; m++){
							for(int n=0; n<filter.DimInVoxels[0]; n++){

                                
                                
                                //                                if ((n+a-supb[0]+1)>B[number])
                                //                                sigmax = sigmax_right[number];
                                //                                if ((m+b-supb[1]+1)>D[number])
                                //                                sigmay = sigmay_right[number];
                                //                                if ((l+c-supb[2]+1)>E[number])
                                //                                sigmaz = sigmaz_right[number];
								//std::cout << "corr maximum "<<sigmaz<< std::endl;
								 //filter.Value[filter.index_3Dto1D(n,m,l)] = exp(-0.5*(pow((n-supb[0])/2,2)+pow((m-supb[1])/2,2)+pow((l-supb[2])/2,2)));
                                
                                
                                
                                double xdata = cos_z*n - sin_z*m;
                                double ydata = cos_x*sin_z*n + cos_x*cos_z*m - sin_x*l;
                                double zdata = sin_x*sin_z*n + sin_x*cos_z*m + cos_x*l;
                                
                                double x0rot = cos_z*B[number]- sin_z*D[number];
                                double y0rot = cos_x*sin_z*B[number]+cos_x*cos_z*D[number]- sin_x*E[number];
                                double z0rot = sin_x*sin_z*B[number]+sin_x*cos_z*D[number]+ cos_x*E[number];
                                
                             //   filter.Value[filter.index_3Dto1D(n,m,l)]=A_psf[number]*exp(-0.5*(pow((n-4)/sigmax,2)+pow((m-4)/sigmay,2)+pow((l-4)/sigmaz,2)));
                                								//filter.Value[filter.index_3Dto1D(n,m,l)]=1;
                                 filter.Value[filter.index_3Dto1D(n,m,l)]=A_psf[number]*exp(-0.5*(pow((xdata-x0rot)/sigmax,2)+pow((ydata-y0rot)/sigmay,2)+pow((zdata-z0rot)/sigmaz,2)));
								//filter.Value[filter.index_3Dto1D(n,m,l)]=1;
								//filter.Value[filter.index_3Dto1D(n,m,l)]=exp(-0.5*(pow((n-supb[0]+a-B[number]+1)/sigmax,2)+pow((m-supb[1]+b-D[number]+1)/sigmay,2)+pow((l-supb[2]+c-E[number]+1)/sigmaz,2)));
								sum  += filter.Value[filter.index_3Dto1D(n,m,l)];


							}}}
//	filter.Value[filter.index_3Dto1D(10,10,10)] = 1;
//	std::cout << "sum "<< sum<< std::endl;
	for(int l=0; l<filter.DimInVoxels[2]; l++){
						for(int m=0; m<filter.DimInVoxels[1]; m++){
							for(int n=0; n<filter.DimInVoxels[0]; n++){
									filter.Value[filter.index_3Dto1D(n,m,l)] = filter.Value[filter.index_3Dto1D(n,m,l)]/sum;
									//somme += filter.Value[filter.index_3Dto1D(i,j,k)];
								}}}

	time = clock() - time;
    

}

*/

/*

void Image::calcul_filter(Image &filter, int a, int b, int c, int number, const std::vector<double>& A_psf, const std::vector<double>& B, const std::vector<double>& D,const std::vector<double>& E, const std::vector<double>& sigmax_left, const std::vector<double>& sigmax_right, const std::vector<double>& sigmay_left, const std::vector<double>& sigmay_right, const std::vector<double>& sigmaz_left, const std::vector<double>& sigmaz_right)
{
    
    clock_t time;
    time = clock();
    //std::cout << "volume"<<filter.DimInVoxels[0]<<std::endl;
    double sum=0.0;
    double somme;
    div_t divresult;
    std::vector<int> infb(3,0), supb(3,0);
    divresult = div (filter.DimInVoxels[0],2);
    if (divresult.rem)
        infb[0] = divresult.quot;
    else
        infb[0] = divresult.quot-1;
    supb[0] = divresult.quot;
    
    divresult = div (filter.DimInVoxels[1],2);
    if (divresult.rem)
        infb[1] = divresult.quot;
    else
        infb[1] = divresult.quot-1;
    supb[1] = divresult.quot;
    
    divresult = div (filter.DimInVoxels[2],2);
    if (divresult.rem)
        infb[2] = divresult.quot;
    else
        infb[2] = divresult.quot-1;
    supb[2] = divresult.quot;
    //std::cout << "corr maximum "<<sigmax[number]<< std::endl;
    //std::cout << "corr maximum "<< B[number]<< std::endl;
    
    double sigmax = sigmax_left[number];
    double sigmay = sigmay_left[number];
    double sigmaz = sigmaz_left[number];
    
    
    
    
    for(int l=0; l<filter.DimInVoxels[2]; l++){
        for(int m=0; m<filter.DimInVoxels[1]; m++){
            for(int n=0; n<filter.DimInVoxels[0]; n++){
                if (n>B[number])
                    sigmax = sigmax_right[number];
                if (m>D[number])
                    sigmay = sigmay_right[number];
                if (l>E[number])
                    sigmaz = sigmaz_right[number];
                
                
                //                                if ((n+a-supb[0]+1)>B[number])
                //                                sigmax = sigmax_right[number];
                //                                if ((m+b-supb[1]+1)>D[number])
                //                                sigmay = sigmay_right[number];
                //                                if ((l+c-supb[2]+1)>E[number])
                //                                sigmaz = sigmaz_right[number];
                //std::cout << "corr maximum "<<sigmaz<< std::endl;
                //filter.Value[filter.index_3Dto1D(n,m,l)] = exp(-0.5*(pow((n-supb[0])/2,2)+pow((m-supb[1])/2,2)+pow((l-supb[2])/2,2)));
                filter.Value[filter.index_3Dto1D(n,m,l)]=A_psf[number]*exp(-0.5*(pow((n-B[number])/sigmax,2)+pow((m-D[number])/sigmay,2)+pow((l-E[number])/sigmaz,2)));
                //filter.Value[filter.index_3Dto1D(n,m,l)]=1;
                //filter.Value[filter.index_3Dto1D(n,m,l)]=exp(-0.5*(pow((n-supb[0]+a-B[number]+1)/sigmax,2)+pow((m-supb[1]+b-D[number]+1)/sigmay,2)+pow((l-supb[2]+c-E[number]+1)/sigmaz,2)));
                sum  += filter.Value[filter.index_3Dto1D(n,m,l)];

            }}}
    //    filter.Value[filter.index_3Dto1D(10,10,10)] = 1;
    //    std::cout << "sum "<< sum<< std::endl;
    for(int l=0; l<filter.DimInVoxels[2]; l++){
        for(int m=0; m<filter.DimInVoxels[1]; m++){
            for(int n=0; n<filter.DimInVoxels[0]; n++){
                filter.Value[filter.index_3Dto1D(n,m,l)] = filter.Value[filter.index_3Dto1D(n,m,l)]/sum;
                //somme += filter.Value[filter.index_3Dto1D(i,j,k)];
            }}}
    
    time = clock() - time;
    
    
 
 }*/


void Image::TV_regularization(Image &im, Image &im_TV, const Image &sens, Image &VolumeTV, Image &VolumeTV_pad, double alpha_TV, int Niter_TV, double d)
{


	double A, B;
	A = 1;
	B = 1;
	double EPSILON = 1e-8;
    double sum_im;
    double t[im.NbVoxels];
    std::cout<< "alpha TV: "<<alpha_TV <<" \n";

    double max;

    for (int k=0; k<im.NbVoxels; k++){
    	 if((im_TV.Value[k]/sens.Value[k])>max)
    	           max =im_TV.Value[k]/sens.Value[k];
    }
    for (int k=0; k<im.NbVoxels; k++){
    	t[k] = im_TV.Value[k]/sens.Value[k]/max;
    	if(t[k]<EPSILON)
    	           t[k] =EPSILON;
    }
	for (int ii = 0; ii < Niter_TV; ii++)
	{


    	VolumeTV.initialize(0);
    	VolumeTV_pad.initialize(0);
		double v1, v2, v3, v4;//compute gradient of ||f||_tv
		double sum = 0;
        double TV_im=0;

		for (int k = 0; k < im.DimInVoxels[2]; k++)
			for (int j = 0; j < im.DimInVoxels[1]; j++)
				for (int i = 0; i < im.DimInVoxels[0]; i++)
				{
					VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i+1,j+1,k+1)] = im.Value[im.index_3Dto1D(i,j,k)];

				}

				for (int k = 0; k < im.DimInVoxels[2]; k++)
					for (int j = 0; j < im.DimInVoxels[1]; j++)
						for (int i = 0; i < im.DimInVoxels[0]; i++)
						{
												double a = VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i,j+1,k+1)];
												double b = VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i+1,j,k+1)];
												double c = VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i+1,j+1,k)];
												double d = VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i+1,j+1,k+1)];
											    double e = VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i+2,j+1,k+1)];
											    double f = VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i+1,j+2,k+1)];
											    double g = VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i+1,j+1,k+2)];
							v1 = (A*(d-a)+B*(d-b)+d-c)/sqrt(A*pow(d-a,2) + B*pow(d-b,2) + pow(d-c,2) + EPSILON);

							v2 = (A*(e-d))/sqrt(A*pow(e-d,2) + B*pow(e-VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i+2,j,k+1)],2) + pow(e-VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i+2,j+1,k)],2) + EPSILON);

							v3 = (B*(f-d))/sqrt(A*pow(f-VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i,j+2,k+1)],2) + B*pow(f-d,2) + pow(f-VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i+1,j+2,k)],2) + EPSILON);

							v4 = (g-d)/sqrt(A*pow(g-VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i,j+1,k+2)],2) + B*pow(g-VolumeTV_pad.Value[VolumeTV_pad.index_3Dto1D(i+1,j,k+2)],2) + pow(g-d,2) + EPSILON);

							VolumeTV.Value[VolumeTV.index_3Dto1D(i,j,k)] = v1 - v2 - v3 - v4;

							//TV_im += sqrt(pow(d-c,2)+pow(d-a,2)+pow(d-b,2));



						}




		for (int k=0; k<im.NbVoxels; k++){
			sum += pow(VolumeTV.Value[k],2);}

		for (int k=0; k<im.NbVoxels; k++){
			     im.Value[k] = im.Value[k]- alpha_TV*d*VolumeTV.Value[k]/sqrt(sum)*t[k];

			     if(im.Value[k]<0)
           		 im.Value[k]=EPSILON;
		}

	}
}


void Image::TV_regularization_ca(double log_likelihood, Image &im, Image &im_TV, const Image &sens, Image &VolumeTV, Image &divp, Image &lastdivp, Image &TV_im, Image &p3, Image &p2 ,Image &p1, double alpha_TV, int Niter_TV)
{



	double EPSILON = 1e-8;
    double dt=0.125;
    double Tol=0.1;
    //double TV=1;
    double sup=0;
    double last_sup=0;
    double im_0[im.NbVoxels];
    double t[im.NbVoxels];

    double TV=0;
    double last_TV=0;
    double im_0_sum=0;
    double alpha_ideal =0;

    double max=0;
    double max_sens = sens.maximum();

    for (int k=0; k<im.NbVoxels; k++){
         im_0[k] =im.Value[k];
         im_0_sum +=im.Value[k];
    }
    im_0_sum = im_0_sum/im.NbVoxels;
    for (int k=0; k<im.NbVoxels; k++){
    	 if((im_TV.Value[k]/sens.Value[k])>max)
    	           max =im_TV.Value[k]/sens.Value[k];
    }
	//for (int ii = 0; ii < Niter_TV; ii++){
    for (int k=0; k<im.NbVoxels; k++){
    	t[k] = im_TV.Value[k]/sens.Value[k]/max;
        //t[k] =1;
    	if(t[k]<EPSILON)
    	           t[k] =EPSILON;
    }

                last_TV = TV;
		last_sup = sup;
		std::cout<< "alpha_TV: "<< alpha_TV <<" \n";
		double norm=1;
		TV_im.initialize(0);
		VolumeTV.initialize(0);
		divp.initialize(0);
		lastdivp.initialize(0);
		p1.initialize(0);
		p2.initialize(0);
		p3.initialize(0);

		while(norm>Tol){
    	 for (int k=0; k<im.NbVoxels; k++){
    	 			        lastdivp.Value[k] =divp.Value[k];}


		for (int k=0; k<im.NbVoxels; k++){
			        VolumeTV.Value[k] =divp.Value[k]-im_0[k]/(alpha_TV*(t[k]));

		}
				for (int k = 0; k < im.DimInVoxels[2]-1; k++)
					for (int j = 0; j < im.DimInVoxels[1]-1; j++)
						for (int i = 0; i < im.DimInVoxels[0]-1; i++)
						{
							                    double d = VolumeTV.Value[VolumeTV.index_3Dto1D(i,j,k)];///(alpha_TV*(im.Value[im.index_3Dto1D(i,j,k)]/sens.Value[sens.index_3Dto1D(i,j,k)]/max));
											    double e = VolumeTV.Value[VolumeTV.index_3Dto1D(i+1,j,k)];///(alpha_TV*(im.Value[im.index_3Dto1D(i+1,j,k)]/sens.Value[sens.index_3Dto1D(i+1,j,k)]/max));
											    double f = VolumeTV.Value[VolumeTV.index_3Dto1D(i,j+1,k)];///(alpha_TV*(im.Value[im.index_3Dto1D(i,j+1,k)]/sens.Value[sens.index_3Dto1D(i,j+1,k)]/max));
											    double g = VolumeTV.Value[VolumeTV.index_3Dto1D(i,j,k+1)];///(alpha_TV*(im.Value[im.index_3Dto1D(i,j,k+1)]/sens.Value[sens.index_3Dto1D(i,j,k+1)]/max));
											    	TV_im.Value[TV_im.index_3Dto1D(i,j,k)]= 1+dt*sqrt(pow(d-e,2)+pow(d-f,2)+pow(d-g,2));
							p1.Value[p1.index_3Dto1D(i,j,k)] = (p1.Value[p1.index_3Dto1D(i,j,k)]+dt*(e-d))/TV_im.Value[TV_im.index_3Dto1D(i,j,k)];
							p2.Value[p2.index_3Dto1D(i,j,k)] = (p2.Value[p2.index_3Dto1D(i,j,k)]+dt*(f-d))/TV_im.Value[TV_im.index_3Dto1D(i,j,k)];
							p3.Value[p3.index_3Dto1D(i,j,k)] = (p3.Value[p3.index_3Dto1D(i,j,k)]+dt*(g-d))/TV_im.Value[TV_im.index_3Dto1D(i,j,k)];
						}

				for (int k = 1; k < im.DimInVoxels[2]; k++)
								for (int j = 1; j < im.DimInVoxels[1]; j++)
									for (int i = 1; i < im.DimInVoxels[0]; i++)
									{
										double a = p1.Value[p1.index_3Dto1D(i-1,j,k)];
										double b = p2.Value[p2.index_3Dto1D(i,j-1,k)];
									    double c = p3.Value[p3.index_3Dto1D(i,j,k-1)];

										divp.Value[divp.index_3Dto1D(i,j,k)] = p1.Value[p1.index_3Dto1D(i,j,k)]-a + p2.Value[p2.index_3Dto1D(i,j,k)]-b + p3.Value[p3.index_3Dto1D(i,j,k)]-c;


									}

				double sum=0;
				for (int k=0; k<im.NbVoxels; k++){
						        sum += pow((lastdivp.Value[k]-divp.Value[k]),2);}
				norm = sqrt(sum);


		}
    double im_sum=0;
	for (int k=0; k<im.NbVoxels; k++){
		 im.Value[k] = im_0[k]-divp.Value[k]*(alpha_TV*(t[k]));

	     im_sum +=im.Value[k];}
	im_sum = im_sum/im.NbVoxels;
	 sup=0;
	 TV=0;
    double norm_divp=0;
	for (int k = 1; k < im.DimInVoxels[2]; k++)
		for (int j = 1; j < im.DimInVoxels[1]; j++)
			for (int i = 1; i < im.DimInVoxels[0]; i++)
			{
									double h = im.Value[im.index_3Dto1D(i,j,k)];
								    double l = im.Value[im.index_3Dto1D(i-1,j,k)];
								    double m = im.Value[im.index_3Dto1D(i,j-1,k)];
								    double n = im.Value[im.index_3Dto1D(i,j,k-1)];
				TV += sqrt(pow(h-l,2)+pow(h-m,2)+pow(h-n,2));
			}
	std::cout<< "TV: "<< TV/im_sum <<" \n";
	for (int k=0; k<im.NbVoxels; k++){
		         sup += pow((im_0[k]/im_0_sum- im.Value[k]/im_sum),2);}
	sup = sup/(2*TV/im_sum);

	for (int k = 1; k < divp.DimInVoxels[2]; k++)
		for (int j = 1; j < divp.DimInVoxels[1]; j++)
			for (int i = 1; i < divp.DimInVoxels[0]; i++)
			{
									double h = divp.Value[divp.index_3Dto1D(i,j,k)];
								    double l = divp.Value[divp.index_3Dto1D(i-1,j,k)];
								    double m = divp.Value[divp.index_3Dto1D(i,j-1,k)];
								    double n = divp.Value[divp.index_3Dto1D(i,j,k-1)];
				norm_divp += pow(h-l,2)+pow(h-m,2)+pow(h-n,2);
			}
	norm_divp = sqrt(norm_divp);
    alpha_ideal = log_likelihood/norm_divp;
    std::cout<< "proportion: "<< alpha_ideal <<" \n";


}

void Image::TV_regularization_vm(double log_likelihood, Image &im, const Image &sens, Image &VolumeTV, Image &divp, Image &lastdivp, Image &TV_im, Image &p3, Image &p2 ,Image &p1,Image &z3, Image &z2 ,Image &z1, double alpha_TV, int Niter_TV, double EPSILON)
{



	//double EPSILON = 1e-8;
    double l_iter = 0;
    double Tol=0.02;


    //double sens_norm[im.NbVoxels];

    double alpha_ideal =0;

    double max=0;



    double max_sens = sens.maximum();
    double min_sens = sens.minimum();




    for (int k=0; k<im.NbVoxels; k++){
         if((im.Value[k]*sens.Value[k])>max)
           max = im.Value[k]*sens.Value[k];
    }


		double norm=1;
		TV_im.initialize(1);
		VolumeTV.initialize(0);
		divp.initialize(1);
		lastdivp.initialize(0);
		p1.initialize(1);
		p2.initialize(1);
		p3.initialize(1);
		z1.initialize(0);
		z2.initialize(0);
		z3.initialize(0);
		//sens.initialize(1);
		double dt=pow(min_sens-6*alpha_TV, 2.0)/(12*alpha_TV*max);
        double dtt=dt;
		while((norm>Tol*dtt)&&(l_iter<2000)){
		//while(l_iter<20000){
			l_iter = l_iter+1;
			for (int k=0; k<im.NbVoxels; k++){
				    	 			        lastdivp.Value[k] =divp.Value[k];}
		 for (int k=0; k<im.NbVoxels; k++){
			        VolumeTV.Value[k] = im.Value[k]*sens.Value[k]/(sens.Value[k]+alpha_TV*divp.Value[k]);}
		 //std::cout<< "min VolumeTV: "<<  VolumeTV.maximum() <<" \n";
         dtt=dt;
		 bool flag_tv=0;
			if(VolumeTV.minimum()<=0){
				dt=dt/2;
				 for (int k=0; k<im.NbVoxels; k++){
					        VolumeTV.Value[k] = im.Value[k]*sens.Value[k]/(sens.Value[k]+alpha_TV*lastdivp.Value[k]);}
                flag_tv=1;
                //std::cout<< "VolumeTV_flag: "<< flag_tv <<" \n";
			}
			else
				{
				//if (flag_tv==1)
				dt=dt*2;
				for (int k=0; k<im.NbVoxels; k++){
									        VolumeTV.Value[k] = im.Value[k]*sens.Value[k]/(sens.Value[k]+alpha_TV*lastdivp.Value[k]);}
			}
			//else
			//	dt=dt;
			//std::cout<< "VolumeTV: "<< VolumeTV.minimum() <<" \n";
			//std::cout<< "dt: "<< dt <<" \n";


		 for (int k = 0; k < im.DimInVoxels[2]; k++)
							for (int j = 0; j < im.DimInVoxels[1]; j++)
								for (int i = 0; i < im.DimInVoxels[0]-1; i++)
								{
									z1.Value[z1.index_3Dto1D(i,j,k)] = VolumeTV.Value[VolumeTV.index_3Dto1D(i+1,j,k)]-VolumeTV.Value[VolumeTV.index_3Dto1D(i,j,k)];

								}

for (int k = 0; k < im.DimInVoxels[2]; k++)
							for (int j = 0; j < im.DimInVoxels[1]-1; j++)
								for (int i = 0; i < im.DimInVoxels[0]; i++)
								{
									z2.Value[z2.index_3Dto1D(i,j,k)] = VolumeTV.Value[VolumeTV.index_3Dto1D(i,j+1,k)]-VolumeTV.Value[VolumeTV.index_3Dto1D(i,j,k)];

								}
								for (int k = 0; k < im.DimInVoxels[2]-1; k++)
															for (int j = 0; j < im.DimInVoxels[1]; j++)
																for (int i = 0; i < im.DimInVoxels[0]; i++)
																{
																	z3.Value[z3.index_3Dto1D(i,j,k)] = VolumeTV.Value[VolumeTV.index_3Dto1D(i,j,k+1)]-VolumeTV.Value[VolumeTV.index_3Dto1D(i,j,k)];

																}



		 for (int k = 0; k < im.DimInVoxels[2]; k++)
					for (int j = 0; j < im.DimInVoxels[1]; j++)
						for (int i = 0; i < im.DimInVoxels[0]; i++)
						{

														TV_im.Value[TV_im.index_3Dto1D(i,j,k)]= 1+dt*sqrt(pow(z1.Value[z1.index_3Dto1D(i,j,k)],2.0)+pow(z2.Value[z2.index_3Dto1D(i,j,k)],2.0)+pow(z3.Value[z3.index_3Dto1D(i,j,k)],2.0));
														p1.Value[p1.index_3Dto1D(i,j,k)] = (p1.Value[p1.index_3Dto1D(i,j,k)]-dt*(z1.Value[z1.index_3Dto1D(i,j,k)]))/TV_im.Value[TV_im.index_3Dto1D(i,j,k)];
														p2.Value[p2.index_3Dto1D(i,j,k)] = (p2.Value[p2.index_3Dto1D(i,j,k)]-dt*(z2.Value[z2.index_3Dto1D(i,j,k)]))/TV_im.Value[TV_im.index_3Dto1D(i,j,k)];
														p3.Value[p3.index_3Dto1D(i,j,k)] = (p3.Value[p3.index_3Dto1D(i,j,k)]-dt*(z3.Value[z3.index_3Dto1D(i,j,k)]))/TV_im.Value[TV_im.index_3Dto1D(i,j,k)];


						}
			//	divp.initialize(0);
				z1.initialize(0);
				z2.initialize(0);
				z3.initialize(0);



				for (int k = 0; k < im.DimInVoxels[2]; k++)
												for (int j = 0; j < im.DimInVoxels[1]; j++)
													for (int i = 1; i < im.DimInVoxels[0]; i++)
													{
														double a = p1.Value[p1.index_3Dto1D(i-1,j,k)];
														z1.Value[z1.index_3Dto1D(i,j,k)] = p1.Value[p1.index_3Dto1D(i,j,k)]-a; }

				for (int k = 0; k < im.DimInVoxels[2]; k++)
												for (int j = 1; j < im.DimInVoxels[1]; j++)
													for (int i = 0; i < im.DimInVoxels[0]; i++)
													{
														double b = p2.Value[p2.index_3Dto1D(i,j-1,k)];
														z2.Value[z2.index_3Dto1D(i,j,k)] =  p2.Value[p2.index_3Dto1D(i,j,k)]-b;
													}

				for (int k = 1; k < im.DimInVoxels[2]; k++)
												for (int j = 0; j < im.DimInVoxels[1]; j++)
													for (int i = 0; i < im.DimInVoxels[0]; i++)
													{
													    double c = p3.Value[p3.index_3Dto1D(i,j,k-1)];
														z3.Value[z3.index_3Dto1D(i,j,k)] = p3.Value[p3.index_3Dto1D(i,j,k)]-c;
													}



				for (int k = 0; k < im.DimInVoxels[2]; k++)
												for (int j = 0; j < im.DimInVoxels[1]; j++)
													for (int i = 0; i < im.DimInVoxels[0]; i++)
													{
														divp.Value[divp.index_3Dto1D(i,j,k)] = z3.Value[z3.index_3Dto1D(i,j,k)]+z2.Value[z2.index_3Dto1D(i,j,k)]+z1.Value[z1.index_3Dto1D(i,j,k)];
													}


				double sum=0;
				for (int k=0; k<im.NbVoxels; k++){
						        sum += pow((lastdivp.Value[k]-divp.Value[k]),2);}
				norm = sqrt(sum);
				//std::cout<< "norm: "<< norm <<" \n";
		}
		//std::cout<< "iteration: "<< l_iter <<" \n";
		for (int k=0; k<im.NbVoxels; k++){
		 im.Value[k] = im.Value[k]*sens.Value[k]/(sens.Value[k]+alpha_TV*divp.Value[k]);
		        if(im.Value[k]<=EPSILON)
	       		im.Value[k]=EPSILON;
		 }

}

void Image::FISTA(Image &im_fista, Image &im_0,Image &im,Image &t1,Image &t2)
{
	for (int k=0; k<t2.NbVoxels; k++){
							        t2.Value[k]=(1+sqrt(4*pow(t1.Value[k],2)+1))/2;}
	for (int k=0; k<im.NbVoxels; k++){
		im_fista.Value[k]=im.Value[k]+(t1.Value[k]-1)/t2.Value[k]*(im.Value[k]-im_0.Value[k]);
		 if(im_fista.Value[k]<=0)
			       		im_fista.Value[k]=1e-32;
	}
	for (int k=0; k<t2.NbVoxels; k++){
								        t1.Value[k]=t2.Value[k];}
	std::cout<< "maximum: "<< im.maximum() <<" \n";
	std::cout<< "maximum t2: "<< t2.maximum() <<" \n";
}
