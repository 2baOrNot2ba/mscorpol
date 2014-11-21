#include <iostream>
#include <cmath>
#include <cstdlib>
#include "ElementResponse.h"

using namespace std;

int main(int argc, char* argv[]){
  double freq;
  double az;
  double el;
  complex<double> response[2][2];

  if (argc != 4){
    cerr << "Wrong number of arguments" << endl
         << "Usage: ElementResponseMain freq az el" << endl;
    exit(1);
  }
  freq=atof(argv[1]);
  az=atof(argv[2]);
  el=atof(argv[3]);

  //Updated by Andrea Horneffer 2013-02-15
  if ((freq>=100e6) && (freq<=300e6)) {
    LOFAR::element_response_hba(freq, az, el,(response));
    cout << response[0][0] << ';' << response[0][1] << endl 
	 << response[1][0] << ';' << response[1][1] << endl;
  } else if ((freq>=10e6) && (freq<=100e6)) {
    LOFAR::element_response_lba(freq, az, el,(response));
    cout << response[0][0] << ';' << response[0][1] << endl 
	 << response[1][0] << ';' << response[1][1] << endl;
  } else {
    cerr << "Unsupported frequency value: " << freq << endl
	 << "frequency needs to be within 10e6 and 300 e6 Hz" << endl;
  }
}
