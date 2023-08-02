#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <pthread.h>


double *pr;
mwIndex *ir;
mwIndex *jc;
double *Br;
double *Bi;
double *Cr;
double *Ci;
size_t m,n;      /* matrix dimensions */

struct ThreadData{
    unsigned int colStart, colEnd;
};

void* multiplyReal(void *td){
    mwIndex c_idx;
    mwIndex sparseStart;
    mwIndex sparseEnd;
    mwIndex sparseIdx;
    mwIndex b_idx;
    double sparseVal;
    struct ThreadData *threadDat = td;
    
    // Loop through columns (which can be though of as rows of A since we have A-transpose) in a thread-safe way
    for (c_idx=threadDat->colStart; c_idx<=threadDat->colEnd; ++c_idx){
        // Get starting/ending indexes into sparse ir/pr matrices for this column 
        sparseStart = jc[c_idx];
        sparseEnd = jc[c_idx+1];
//         printf("Calculating output[%d] from jc[%d:5d]\n",c_idx,sparseStart,sparseEnd);
        
        // Initialize to zero
        Cr[c_idx] = 0;
        
        // Loop through nonzero entries in sparse matrix and accumulate output matrix
        for (sparseIdx=sparseStart; sparseIdx<sparseEnd; ++sparseIdx){
            b_idx = ir[sparseIdx];
            sparseVal = pr[sparseIdx];
            Cr[c_idx] += Br[b_idx]*sparseVal;
//             printf("\tAdded Br[ir[%d]]*pr[%d]=Br[%d]*%f=%f*%f=%f to give Cr[%d]=%f\n",sparseIdx,sparseIdx,b_idx,sparseVal,Br[b_idx],sparseVal,Br[b_idx]*sparseVal,c_idx,Cr[c_idx]);
        }
    }
    return NULL;
}

void* multiplyComplex(void *td){
    mwIndex c_idx;
    mwIndex sparseStart;
    mwIndex sparseEnd;
    mwIndex sparseIdx;
    mwIndex b_idx;
    double sparseVal;
    struct ThreadData *threadDat = td;
    
    // Loop through columns (which can be though of as rows of A since we have A-transpose) in a thread-safe way
    for (c_idx=threadDat->colStart; c_idx<=threadDat->colEnd; ++c_idx){
        // Get starting/ending indexes into sparse ir/pr matrices for this column 
        sparseStart = jc[c_idx];
        sparseEnd = jc[c_idx+1];
        
        // Initialize to zero
        Cr[c_idx] = 0;
        Ci[c_idx] = 0;
        
        // Loop through nonzero entries in sparse matrix and accumulate output matrix
        for (sparseIdx=sparseStart; sparseIdx<sparseEnd; ++sparseIdx){
            b_idx = ir[sparseIdx];
            sparseVal = pr[sparseIdx];
            Cr[c_idx] += Br[b_idx]*sparseVal;
            Ci[c_idx] += Bi[b_idx]*sparseVal;
        }
    }
    return NULL;
}

// This function is used to solve problems of the form: C = A*B.
// The inputs to this function are:
//     1) a sparse real matrix, A-transpose (NOT A)
//     2) a (optionally complex) vector (or matrix, which will be treated as a vector) B 
//     3) a (optionally complex) vector (or matrix, which will be treated as a vector) C
// Note that since A is stored in the compressed column storage format, 
// it is only threadsafe to compute C = (b*A')', thus this function 
// requires A-transpose (similar to just storing in row compressed storage format)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    unsigned int *nthreads;
    
    struct ThreadData *threadData;
    pthread_t *thread;
    unsigned int nPerThread, iThread;;
    bool complex; 
        
    pr = (double*)mxGetData(prhs[0]);/* sparse matrix A-transpose is assumed to be all real */ 
    ir = mxGetIr(prhs[0]);
    jc = mxGetJc(prhs[0]);
    m = mxGetM(prhs[0]); //number of elements in B since we passed in A-transpose
    n = mxGetN(prhs[0]); //number of elements in C since we passed in A-transpose
        
    Br = (double*)mxGetData(prhs[1]); /* second full matrix - can be real or complex, but only complex for now*/
    Bi = (double*)mxGetImagData(prhs[1]); /* second full matrix - can be real or complex, but only complex for now*/
    complex = mxIsComplex(prhs[1]);
    
    Cr = (double*)mxGetData(prhs[2]); /* third full matrix - can be real or complex, but only complex for now*/
    Ci = (double*)mxGetImagData(prhs[2]); /* third full matrix - can be real or complex, but only complex for now*/
    
    if(mxIsComplex(prhs[2]) != complex){
        mexErrMsgIdAndTxt("MATLAB:odearguments:InconsistentDataType","B and C must both either be real or complex.");
    }
    
    /* INPUT 4 - NTHREADS */
    nthreads = (unsigned int*)mxGetData(prhs[3]);
    threadData = calloc(*nthreads,sizeof(struct ThreadData));
    thread = calloc(*nthreads,sizeof(pthread_t));
    
    nPerThread = floor(((double)n)/(*nthreads));
//     printf("CREATING %d THREADS with %d PTS TO GRID EACH\n",*nthreads,nPerThread);
    for(iThread=0; iThread<*nthreads; iThread++){
        threadData[iThread].colStart = iThread*nPerThread;
        threadData[iThread].colEnd = (iThread+1)*nPerThread-1;
    }
    threadData[*nthreads-1].colEnd = n-1;
    
    for(iThread=0; iThread<*nthreads; iThread++){
//         printf("   STARTING THREADS %d\n",iThread);
        if(complex){
            pthread_create(&thread[iThread], NULL, multiplyComplex, &threadData[iThread]);
        }else{
            pthread_create(&thread[iThread], NULL, multiplyReal, &threadData[iThread]);
        }
    }
    
    /* Wait for threads to finish */
    for(iThread=0; iThread<*nthreads; iThread++){
        pthread_join(thread[iThread],NULL);
//         printf("   THREADS %d FINISHED\n",iThread);
    }
    
    free(threadData);
    free(thread);
}