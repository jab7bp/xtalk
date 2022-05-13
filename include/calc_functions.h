#ifndef calc_functions
#define calc_functions

double calc_mean(double *array){
	int n = (sizeof(array)/sizeof(array[0]));
	double sum = 0.0, mean;

	for(int i = 0; i < n; i++){
		sum += array[i];
	}
	mean = sum/n;
	return mean;
}

double calc_StdDev(double *array){
	int n = sizeof(array)/sizeof(array[0]);
	cout << "n: " << n << endl;

	return 0;
}

#endif