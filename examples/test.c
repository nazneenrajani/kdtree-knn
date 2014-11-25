/*! gcc -Wall -g -o test test.c libkdtree.a */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>
#include "kdtree.h"

unsigned int get_msec(void)
{
	static struct timeval timeval, first_timeval;

	gettimeofday(&timeval, 0);

	if(first_timeval.tv_sec == 0) {
		first_timeval = timeval;
		return 0;
	}
	return (timeval.tv_sec - first_timeval.tv_sec) * 1000 + (timeval.tv_usec - first_timeval.tv_usec) / 1000;
}


int main(int argc, char **argv)
{
	int i, j, vcount = 10;
	struct kdres *set;
    struct kdtree *kd;
	unsigned int msec, start;

	if(argc > 1 && isdigit(argv[1][0])) {
		vcount = atoi(argv[1]);
	}
	fflush(stdout);
    int n = 1965;                 // # data points
    int D = 560;                  // data dimension
    char fname[] = "/Users/nrajani/Downloads/knn/data/frey.dat"; // data file
    double **X = NULL;
    X= (double*)malloc(sizeof(double)*n);
    for (i=0; i<n; i++) {
        X[i] = NULL;
        X[i] = (double*)malloc(sizeof(double)*D);
    }
    //new_2d_array(X, n, D);
    if (read_X_from_file(X, n, D, fname)==0)
    {
        printf("error open file.\n");
        exit(1);
    }
	kd = kd_create(D);
    double pos[D];
	start = get_msec();
	for(i=0; i<n; i++) {
        for (j=0; j<D; j++) {
            pos[j] = X[i][j];
        }
        assert(kd_insert(kd, pos, i+1) == 0);
	}
	msec = get_msec() - start;
	printf("%.3f sec for creating tree \n", (float)msec / 1000.0);
	start = get_msec();
    for (i=0; i<n; i++) {
        set = kd_nearest_n(kd, X[i],12);
    }
	msec = get_msec() - start;
	printf("%.5f sec for querying\n", (float)msec / 1000.0);
    for (i =0; i<kd_res_size(set); i++) {
        printf("%i \n",kd_res_next(set));
    }
	kd_res_free(set);

	kd_free(kd);
	return 0;
}
int read_X_from_file(double **X, int n, int D, char *filename)
{
    FILE *fp = NULL;
    if (!(fp = fopen(filename, "rb")))
    return 0;
    int i;
    int num_in, num_total;
    for (i = 0; i < n; i++)
    {
        num_total = 0;
        while (num_total < D)
        {
            num_in = fread(X[i]+num_total, sizeof(double), D, fp);
            num_total += num_in;
        }
    }
    
    fclose(fp);
    
    //printf("done\n");
    
    return 1;
}




