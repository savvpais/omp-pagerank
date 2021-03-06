#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>

int ThreadNumber;

double calcError(double *a, double *b, int N);
void pageRankGaussSeidel(int N,int **adjMat,int *degreeArray,int *inboundArray,double *pageRankVec);
void test (int N, double *pageRankVec, char *testfile);

int main(int argc , char **argv){
	if (argc < 3 || argc > 4) {
        printf("Error using %s! Usage is: ./%s <dataset>.txt <ThreadNumber> <testfile>.txt (testfile is optional) \n", argv[0], argv[0]);
        exit(1);
    }
    
    ThreadNumber = atoi(argv[2]);
    
    int i, j, N;     
    char *filename,line[512],s[2] = " ", *token;
    filename = argv[1];

    FILE *file = fopen(filename, "r");
    
    if(file == NULL){
        printf("Error opening file %s\n",filename);   
        exit(1);             
    }
    
    while (fgets(line, 512, file) != NULL){
        if (line[0] == '#'){
            token = strtok(line,s);
            token = strtok(NULL,s);
            if ( !strcmp(token,"Nodes:") ){
                token = strtok(NULL,s);
                
                N = atoi(token);
                
            }
        } 
    }
    
    //Debug
    //~ printf("Number of nodes: %d \n",N);
    
    fseek(file, 0, SEEK_SET);
    //Check if N is correct
	int from, to;
    while (fgets(line, 512, file) != NULL){
        token = strtok(line,s);
       
        if (strcmp(token, "#")){
            sscanf(line, "%d %d \n", &from, &to);
            if (from>N)
                N = from;
            if (to>N)
                N = to;
        } 
       
    }
    
    int *degreeArray;
    degreeArray = malloc(sizeof(int) * N);
    if (degreeArray == NULL){
        printf("Could not allocate memory for the degreeArray \n");
        exit(-1);
    
    }
    
    int *inboundArray;
    inboundArray = malloc(sizeof(int) * N);
    if (inboundArray == NULL){
        printf("Could not allocate memory for the inboundArray \n");
        exit(-1);
    }
    
    for (i = 0; i < N; i++)
        inboundArray[i] = 0;
    
    

    fseek(file, 0, SEEK_SET);
	//Count outbound and inbound links for each page
    while (fgets(line, 512, file) != NULL){
        
        token = strtok(line,s);
        if (strcmp(token, "#")){
            sscanf(line, "%d %d \n", &from, &to);         
            degreeArray[from]++;    
            inboundArray[to]++;     
        }
    }
   
  
	int **adjMat;
    adjMat = malloc(sizeof(int*) * N);
    if (adjMat == NULL){
        printf("Could not allocate memory for the adjMat array\n");
        exit(-1);
    }

    for (i = 0; i <= N;i++)
        adjMat[i] = malloc(sizeof(int) * inboundArray[i]);    
    
    int *counter;
    counter = malloc(sizeof(int) * N);
    if (counter == NULL){
        printf("Could not allocate memory for the counter array \n");
        exit(-1);
    }
    
    for (i = 0; i < N; i++)
        counter[i] = 0;
    
    fseek(file, 0, SEEK_SET);
    //Get the adjecency matrix
    while (fgets(line, 512, file) != NULL){
        
        token = strtok(line, s);
        if (strcmp(token,"#")){
            sscanf(line, "%d %d \n", &from, &to);
            adjMat[to][counter[to]] = from;
            counter[to]++;
        }

    }
    
    double *pageRankVec;
    pageRankVec = malloc(sizeof(double) * N);
    if (pageRankVec == NULL){
        printf("Could not allocate memory for the pageRank vector \n");
        exit(-1);
    }
        
    struct timeval startwtime, endwtime;
    gettimeofday (&startwtime, NULL);
    pageRankGaussSeidel(N, adjMat, degreeArray, inboundArray, pageRankVec);
    gettimeofday (&endwtime, NULL);  
    
    
    double time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

    //Debug
    //~ for (i=0;i<10;i++)
        //~ printf("%.20f \n",pageRankVec[i]);
   
    
    printf("Serial PageRank time: %f \n", time);
    
    
    free(inboundArray);
    free(degreeArray);
    free(pageRankVec);
    free(counter);
    
    
    for (i = 0; i < N; i++)
        free(adjMat[i]);
    
    free(adjMat);
    
   
    return 0;
    
}

double calcError(double *a, double *b, int N){
	double maximum = -500;
	int i;
	double absDiff;
	
	omp_set_num_threads(ThreadNumber);
	#pragma omp parallel private(i,absDiff) shared(maximum,a,b)
	{ 
		#pragma omp for
		for (i = 0; i < N; ++i){
		absDiff = fabs(a[i] - b[i]);

		if (maximum < absDiff)
			maximum = absDiff;
		}
	}
	return maximum;
}

void pageRankGaussSeidel(int N, int **adjMat, int *degreeArray, int *inboundArray, double *pageRankVec){
    int i, j;
    double error=1, danglingSum, iterationSum;
    
    //Damping factor
    double alpha = 0.85;
    double delta = 1 - alpha;
    
    //PageRank vector of previous iteration is saved in z
    double *z;
    z = malloc(sizeof(double) * N);
    if (z == NULL){
        printf("Could not allocate memory for the z array");
        exit(-1);
    }
    
    //PageRank vector initiallization
    for (i = 0; i < N; i++)
        pageRankVec[i] = 1.0 / N;
    
    int count = 0;
    
    omp_set_num_threads(ThreadNumber);
    while (error > 1e-6){
        
        danglingSum=0;
        
        //Compute dangling node sum
		#pragma omp parallel private(j) shared(degreeArray,N)
        {
			#pragma omp for reduction(+: danglingSum)
			for(j = 0; j < N; j++){
				if(degreeArray[j] == 0){
					danglingSum = danglingSum + pageRankVec[j] / N;
				}
			}
		}
		
		#pragma omp parallel shared(N, inboundArray, pageRankVec, adjMat, degreeArray, alpha, delta) private(i, j, iterationSum)
		{
			#pragma omp for
			for(i = 0; i < N; i++){
				z[i] = pageRankVec[i];
				iterationSum = 0;
				for (j=0;j<inboundArray[i];j++) 
					iterationSum = iterationSum  + pageRankVec[adjMat[i][j]] / degreeArray[adjMat[i][j]];
				
				pageRankVec[i] = alpha * (danglingSum + iterationSum) + delta / (double)N;
				
			}
		}
        count++;
        error = calcError(pageRankVec, z, N);
    }
	//~ printf("PageRank iterations:%d \n",count);
}

void test (int N, double *pageRankVec, char *testfile){
	FILE *file = fopen(testfile, "r");
	if(file==NULL){
		printf("Could not open %s\n",testfile);
		exit(1);
	}
	
	double *testVec;
	testVec = malloc(sizeof(double) * N);
	if (testVec == NULL){
        printf("Could not allocate memory for the testVec array \n");
        exit(-1);
    }
    
    int i, k;
    for (i = 0; i < N; i++)
		k = fscanf(file,"%lf",&testVec[i]);

	int failed = 0;
	for(i = 0; i < N; i++){
		if (fabs(pageRankVec[i]-(double)testVec[i])/fabs((double)testVec[i]) > 0.1)
			failed++;
	}
	
	if((double)(N-failed) / (double)N * 100 < 95.0)
		printf("Test failed!\n");
	else
		printf("Test passed!\n");
}
