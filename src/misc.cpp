/*
 * Miscellaneous functions
 */

#include <cmath>
#include <stdlib.h>
#include <numeric>
#include <iomanip>
#include <mpi.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_eigen.h>

#include "misc.hpp"

using namespace std;

/*
 * Get the date in string format YYMMDD
 */
string getDate ()
{
  char s[10];
  time_t t = time(0);
  strftime(s, 10, "%Y%m%d", localtime(&t));
  string outTime(s,s+8);
  return outTime;
}

/*
 * Trim a string from spaces
 */
string trim (const string &orig)
{
  string::size_type p1=orig.find_first_not_of(" \t");
  if (p1==string::npos) return "";
  string::size_type p2=orig.find_last_not_of(" \t");
  return orig.substr(p1,p2-p1+1);
}

vector<double> angToCart(double theta, double phi)
{
  vector<double> temp(3,0.0);
  
  temp[0] = sin(theta)*cos(phi);
  temp[1] = sin(theta)*sin(phi);
  temp[2] = cos(theta);
  return temp;
}

vector<double> cartToAng(vector<double> cart)
{
  vector<double> output(2,0.0);
  // Theta
  output[0] = atan2(sqrt(cart[0]*cart[0]+cart[1]*cart[1]),cart[2]);
  // Phi
  output[1] = atan2(cart[1],cart[0]);
  if (output[1]<0.) 
    output[1] += TWOPI;
  
  return output;
}

vector<string> 
getDetectorIds(int frequency)
{
  vector<string> output;
  switch (frequency)
    {
    case 30:
      output.push_back("LFI27M-00");
      output.push_back("LFI27S-10");
      output.push_back("LFI28M-00");
      output.push_back("LFI28S-10");
      break;
    case 44:
      output.push_back("LFI24M-00");
      output.push_back("LFI24S-10");
      output.push_back("LFI25M-00");
      output.push_back("LFI25S-10");
      output.push_back("LFI26M-00");
      output.push_back("LFI26S-10");
      break;
    case 70:
      output.push_back("LFI18M-00");
      output.push_back("LFI18S-10");
      output.push_back("LFI19M-00");
      output.push_back("LFI19S-10");
      output.push_back("LFI20M-00");
      output.push_back("LFI20S-10");
      output.push_back("LFI21M-00");
      output.push_back("LFI21S-10");
      output.push_back("LFI22M-00");
      output.push_back("LFI22S-10");
      output.push_back("LFI23M-00");
      output.push_back("LFI23S-10");
      break;
    default:
      cout << "Invalid frequency" << endl;
    }
  return output;
}

string 
getDetectorId(int diode)
{
  switch(diode)
    {
    case 0:
      return "M-00";
      break;
    case 1:
      return "M-01";
      break;
    case 2:
      return "S-10";
      break;
    case 3:
      return "S-11";
      break;
    default:
      cout << "Invalid diode" << endl;
    }
  return string("");
}

vector<string> 
getAllDetectorIds()
{
  vector<string> output;
  output.push_back("LFI18M-00");
  output.push_back("LFI18S-10");
  output.push_back("LFI19M-00");
  output.push_back("LFI19S-10");
  output.push_back("LFI20M-00");
  output.push_back("LFI20S-10");
  output.push_back("LFI21M-00");
  output.push_back("LFI21S-10");
  output.push_back("LFI22M-00");
  output.push_back("LFI22S-10");
  output.push_back("LFI23M-00");
  output.push_back("LFI23S-10");
  output.push_back("LFI24M-00");
  output.push_back("LFI24S-10");
  output.push_back("LFI25M-00");
  output.push_back("LFI25S-10");
  output.push_back("LFI26M-00");
  output.push_back("LFI26S-10");
  output.push_back("LFI27M-00");
  output.push_back("LFI27S-10");
  output.push_back("LFI28M-00");
  output.push_back("LFI28S-10");
  return output;
}

vector<string> 
getAllDiodes()
{
  vector<string> output;
  output.push_back("LFI18M-00");
  output.push_back("LFI18M-01");
  output.push_back("LFI18S-10");
  output.push_back("LFI18S-11");
  output.push_back("LFI19M-00");
  output.push_back("LFI19M-01");
  output.push_back("LFI19S-10");
  output.push_back("LFI19S-11");
  output.push_back("LFI20M-00");
  output.push_back("LFI20M-01");
  output.push_back("LFI20S-10");
  output.push_back("LFI20S-11");
  output.push_back("LFI21M-00");
  output.push_back("LFI21M-01");
  output.push_back("LFI21S-10");
  output.push_back("LFI21S-11");
  output.push_back("LFI22M-00");
  output.push_back("LFI22M-01");
  output.push_back("LFI22S-10");
  output.push_back("LFI22S-11");
  output.push_back("LFI23M-00");
  output.push_back("LFI23M-01");
  output.push_back("LFI23S-10");
  output.push_back("LFI23S-11");
  output.push_back("LFI24M-00");
  output.push_back("LFI24M-01");
  output.push_back("LFI24S-10");
  output.push_back("LFI24S-11");
  output.push_back("LFI25M-00");
  output.push_back("LFI25M-01");
  output.push_back("LFI25S-10");
  output.push_back("LFI25S-11");
  output.push_back("LFI26M-00");
  output.push_back("LFI26M-01");
  output.push_back("LFI26S-10");
  output.push_back("LFI26S-11");
  output.push_back("LFI27M-00");
  output.push_back("LFI27M-01");
  output.push_back("LFI27S-10");
  output.push_back("LFI27S-11");
  output.push_back("LFI28M-00");
  output.push_back("LFI28M-01");
  output.push_back("LFI28S-10");
  output.push_back("LFI28S-11");
  return output;
}

int getDetectorIdasInt(string detectorId)
{
  int horn = stringToData<int>(detectorId.substr(3,2));
  int rad = stringToData<int>(detectorId.substr(7,1));
  int detector = stringToData<int>(detectorId.substr(8,1));

  return horn+rad+detector;
}

string intToString(int64 x, unsigned int width)
{
  ostringstream strstrm;
  (x>=0) ? strstrm << setw(width) << setfill('0') << x
    : strstrm << "-" << setw(width-1) << setfill('0') << -x;
  string res = strstrm.str();

  return trim(res);
}

string getExchangeName (const string prefix, const string od, 
			const string suffix)
{
  time_t rawtime;
  time ( &rawtime );
  tm *res = localtime(&rawtime);
  string date_time_zero = intToString(res->tm_year+1900,4)
    +intToString(res->tm_mon+1,2)
    +intToString(res->tm_mday,2);
  
  string filename = prefix
    +"-"
    +od
    +suffix
    +date_time_zero
    +".fits";

  return filename;
}
     
void flagConversion (vector<int> inputFlag, vector<int> & outFlag, vector<int> & commonFlag, bool isHk)
{
  for (unsigned int i=0; i<inputFlag.size(); ++i) 
    {
      // Science quality bit0 flags
      if ( (inputFlag[i] & (1 << 14)) > 0 ) // 100 0000 0000 0000
	outFlag[i] |= (1<<0); 
      else if ( isHk && ((inputFlag[i] & (1 << 2)) > 0) )  // 100
	outFlag[i] |= (1<<0); 
      
      
      // Planet Cross bit2 flags
      if ( (inputFlag[i] & (1 << 18)) > 0 ){ // 100 0000 0000 0000 0000
      	outFlag[i] |= (1<<2);  
	//	cout << "set planet flag" << endl;
      }

      // Moving Objects bit3 flags
      if ( ((inputFlag[i] & (1 << 19)) > 0 ) || 
	   ((inputFlag[i] & (1 << 20)) > 0 ) )
	outFlag[i] |= (1<<3); 

      // Gap bit4 flags
      if ( (inputFlag[i] & (1 << 16)) > 0 )
	outFlag[i] |= (1<<4); 



      
      // -- common flags --//

      // Stable pointing bit0
      if ( (inputFlag[i] & (1<<4) ) > 0 )
	commonFlag[i] |= (1<<0);           // set bit0 IF bit4 = 1 0000
      
      // Time quality bit1
      if ( ( (inputFlag[i] & (1<<5) ) > 0 )||
	   ( (inputFlag[i] & (1<<6) ) > 0 ) ) // check bit5 or bit6 => 110 0000    
	commonFlag[i] |= (1<<1);     // set bit 1 IF bit5 OR bit6 = 1    

      // Special observation bit2
      if ( ( (inputFlag[i] & (1<<22) ) > 0 ) ) // check bit22 => 10 0000 0000 0000 0000 0000    
	commonFlag[i] |= (1<<2);     // set bit 1 IF bit22 = 1      
    }
}

void mpiError (int rankMPI, int start, int stop)
{
  int error = 0;
  if (start>=stop)
    {
      cout << rankMPI << ": too many processes for the requested time range" << endl;
      error=1;
    }

  MPI::COMM_WORLD.Allreduce(&error, &error, 1, MPI::INT, MPI::SUM);
  MPI::COMM_WORLD.Barrier();
}

double computeMean(vector<double> & input)
{
  double sum = 0.0;
  double correction = 0.0;
  
  for(size_t index = 0; index < input.size(); ++index)
    {
      double cur_element_corrected = input[index] - correction;
      double new_sum = sum + cur_element_corrected;
      correction = (new_sum - sum) - cur_element_corrected;
      sum = new_sum;
    }
  
  return sum/static_cast<double>(input.size());
}

double computeVariance(vector<double> & input)
{
  double mean = computeMean(input);

  double sum = 0.0;  
  for(size_t index = 0; index < input.size(); ++index)
    {
      sum += pow((input[index]-mean),2);
    }
  
  return sqrt(sum/static_cast<double>(input.size()));
}

void invert2_eig (vector<double> & cc)
{
  // invert the matrix
  gsl_matrix * components = gsl_matrix_alloc(2,2);
  gsl_matrix * eigenvectors = gsl_matrix_alloc(2,2);
  gsl_vector * eigenvalues = gsl_vector_alloc(2);
  gsl_eigen_symmv_workspace * work = gsl_eigen_symmv_alloc(2);
  
  // set the matrix elements
  gsl_matrix_set(components, 0, 0, cc[0]);
  gsl_matrix_set(components, 1, 0, cc[1]);
  gsl_matrix_set(components, 0, 1, cc[1]);
  gsl_matrix_set(components, 1, 1, cc[2]);
  
  gsl_eigen_symmv (components, eigenvalues, eigenvectors, work);
  
  gsl_eigen_symmv_free(work);
  gsl_matrix_free(components);

  // Get eigenvalues and eigenvectors  
  double d0 = gsl_vector_get(eigenvalues,0);
  double d1 = gsl_vector_get(eigenvalues,1);
  double xd0=d0>1.e-30 ? 1./d0 : 0;
  double xd1=d1>1.e-30 ? 1./d1 : 0;

  double v00 = gsl_matrix_get (eigenvectors, 0, 0);
  double v01 = gsl_matrix_get (eigenvectors, 0, 1);
  double v10 = gsl_matrix_get (eigenvectors, 1, 0);
  double v11 = gsl_matrix_get (eigenvectors, 1, 1);

  cc[0] = xd0*v00*v00 +xd1*v01*v01;
  cc[1] = xd0*v00*v10 +xd1*v01*v11;
  cc[2] = xd0*v10*v10 +xd1*v11*v11;
}

/**
 * Sort pixel numbers and combine equal values
 **/
vector<int> sortAndCount(vector<int> & pixels)
{
  sort (pixels.begin(), pixels.end());

  vector<int> compressedPixels;
  int pix0=-1;
  for (size_t i=0; i<pixels.size(); ++i)
    if (pixels[i]>pix0)
      {
	compressedPixels.push_back(pixels[i]);
	pix0  = pixels[i];
      }

  return compressedPixels;
}

void twiddle (unsigned int &v)
{
  for (int i=0; i<9; ++i)
    {
      v ^= v<<13;
      v ^= v>>17;
      v ^= v<<5;
    }
}
    
/*! Returns uniformly distributed random integer numbers from the
  interval [0;0xFFFFFFFF]. */
unsigned int intRandUniform(unsigned int x, unsigned int y,
			    unsigned int z, unsigned int w)
{
  unsigned int t = x^(x<<11);
  x = y;
  y = z;
  z = w;
  
  return w=(w^(w>>19))^(t^(t>>8));
}

/* Get the Temperature termometer given the horn */
std::string getTFem (int horn)
{
  switch(horn)
    {
    case 18:
      return std::string("LM302332");
    case 19:
      return std::string("LM302332");
    case 20:
      return std::string("LM302332");
    case 21:
      return std::string("LM206332");
    case 22:
      return std::string("LM206332");
    case 23:
      return std::string("LM206332");
    case 24:
      return std::string("LM203332");
    case 25:
      return std::string("LM204332");
    case 26:
      return std::string("LM306332");
    case 27:
      return std::string("LM202332");
    case 28:
      return std::string("LM203332");
    default:
      std::cout << "Invalid horn" << std::endl;
    }
  return std::string("");
}
