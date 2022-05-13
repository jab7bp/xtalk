#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
using namespace std::chrono;

void vector_test(){
	int cnt = 1;
	vector<double> v(1,0.0);
	vector<vector<double> > v2(10, v);
	for(int i = 0; i < 10; i++){

	
		for(int j = 0; j < 7; j++){

			v2[i].push_back(0);
		}
		cnt++;
	}
	// v2[12].push_back(99);
	// v2[3][10] = 99;

for (int i = 0; i < v2.size(); i++)       // loops through each row of vy
   {  for (int j = 0; j < v2[i].size(); j++) // loops through each element of each row 
          cout << " " << v2[i][j];           // prints the jth element of the ith row
      cout << endl;
   }



	
	
}