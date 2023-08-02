#include "mex.h"
#include <math.h>
#include <tgmath.h>
#include <pthread.h>
#include "matrix.h"

unsigned int min(unsigned int a, unsigned int b) {
    return ((a < b) ? a : b);
}

unsigned int max(unsigned int a, unsigned int b) {
    return ((a > b) ? a : b);
}

double *coords_nonCart;
double *half_kernel_extent;
double *kernel_sigma;
double *temp_output_kernel_weights;
double *inv_kern_mult;
mwIndex *temp_output_cart_idx;
mwIndex *output_nNeighbors;
const mwSize *dims_nonCart;
unsigned int *dims_cart;
unsigned int *idx_convert;
mwSize nPts_nonCart;
mwSize nDims;
unsigned int maxNonSparsePerSample;
unsigned int maxNonSparseTotal;

struct ThreadData{
    unsigned int startSample, endSample;
    unsigned int *nonSparseVals;
};

void* go(void *td){
    struct ThreadData *threadDat = td;
    unsigned int *kernel_bounds;
    unsigned int *coord_neighbor;
    double *cur_coord_noncart;
    unsigned int idx_nonCart, idx_cart, nNeighbors, iNeighbor, dim;
    double kernel_weight;
    unsigned int sparseMatIdx;
    unsigned int dim_thread;
    unsigned int *nonSparseThreadCount = threadDat->nonSparseVals;
    double accumVar;
    unsigned int maxIdx;
    
    // Allocate memory
    kernel_bounds = calloc(nDims*2,sizeof(unsigned int));
    coord_neighbor = calloc(nDims,sizeof(unsigned int));
    cur_coord_noncart = calloc(nDims,sizeof(double));
    
    // Loop though each sample point
    *nonSparseThreadCount = 0;
    maxIdx = 0;
    for(idx_nonCart=threadDat->startSample; idx_nonCart<=threadDat->endSample; idx_nonCart++){
//         printf("idx_nonCart = %u\n",idx_nonCart);
        // Calculate local cartesian neighborhood defining kernel extent
        nNeighbors = 1;
        for(dim_thread=0; dim_thread<nDims; dim_thread++){
            cur_coord_noncart[dim_thread] = (double)coords_nonCart[idx_nonCart*nDims+dim_thread];
            
            // Lower bounds
            kernel_bounds[2*dim_thread] = max(ceil(cur_coord_noncart[dim_thread]-half_kernel_extent[dim_thread]),0);
            if((kernel_bounds[2*dim_thread]+half_kernel_extent[dim_thread]) == cur_coord_noncart[dim_thread]){
                // Eliminates points that are exactly extent away
                kernel_bounds[2*dim_thread] = kernel_bounds[2*dim_thread]+1;
            }
            
            // Upper Bounds
            kernel_bounds[2*dim_thread+1] = min(floor(cur_coord_noncart[dim_thread]+half_kernel_extent[dim_thread]),dims_cart[dim_thread]);
            if((kernel_bounds[2*dim_thread+1]-half_kernel_extent[dim_thread]) == cur_coord_noncart[dim_thread]){
                // Eliminates points that are exactly extent away
                kernel_bounds[2*dim_thread+1] = kernel_bounds[2*dim_thread+1]-1;
            }
            
            coord_neighbor[dim_thread] = kernel_bounds[2*dim_thread]; // Start at lower bounds
            nNeighbors = nNeighbors*max(kernel_bounds[2*dim_thread+1]-kernel_bounds[2*dim_thread]+1,1);
        }
        
        if(nNeighbors > maxNonSparsePerSample){
            printf("    TOO MANY NONSPARSE (%d>%d)!\n",nNeighbors,maxNonSparsePerSample);
        }
        
        // Store number of cartesian voxels contributing to this non-Cartesian sample
        output_nNeighbors[idx_nonCart] = nNeighbors;
        
        // Loop through neighborhood and grid
        for(iNeighbor=0; iNeighbor<nNeighbors; iNeighbor++){
            //Grid neighbor
            kernel_weight = 0;
            idx_cart = 0;
            for(dim_thread=0; dim_thread<nDims; dim_thread++){
                // Only calculate exponential once per point
                accumVar = cur_coord_noncart[dim_thread]-coord_neighbor[dim_thread];
                accumVar = accumVar*accumVar*inv_kern_mult[dim_thread];
                kernel_weight = kernel_weight+accumVar;
                // Calculate index for cartesian voxel (neighbor)
                idx_cart = idx_cart + coord_neighbor[dim_thread]*idx_convert[dim_thread];
            }
            kernel_weight = exp(kernel_weight);
            
            // Store kernel values and weights
            sparseMatIdx = idx_nonCart*maxNonSparsePerSample+iNeighbor;         
            temp_output_kernel_weights[sparseMatIdx] = kernel_weight;
            temp_output_cart_idx[sparseMatIdx] = idx_cart;
            
//             printf("Neighbor %d -> vox[x=%d y=%d z=%d = %d], weight=%f, non_cart=%d\n", iNeighbor, 
//                     coord_neighbor[0], coord_neighbor[1], coord_neighbor[2],
//                     idx_cart, kernel_weight),idx_nonCart;
            
            //Update to next neighbor
            for(dim_thread=0; dim_thread<nDims; dim_thread++){
                if(coord_neighbor[dim_thread]==kernel_bounds[2*dim_thread+1]){
                    coord_neighbor[dim_thread] = kernel_bounds[2*dim_thread]; //Reset dimension to lower bound
                }else{
                    coord_neighbor[dim_thread]++;
                    break;
                }
            }
        }
        *nonSparseThreadCount = *nonSparseThreadCount + nNeighbors;
        
        if(*nonSparseThreadCount > maxNonSparseTotal){
            printf("    MATRIX TOO SMALL (%d>%d)!\n",*nonSparseThreadCount,maxNonSparseTotal);
        }
    }
    
    // Free up memory
    free(kernel_bounds);
    free(coord_neighbor);
    free(cur_coord_noncart);
    
    return NULL;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    unsigned int dim;
    unsigned int iThread;
    unsigned int *nthreads;
    unsigned int nonSparseTotalCount;
    mwSignedIndex tempSize[2];
    unsigned int numPtsPerThread;
    double *kernel_extent;
    size_t nPts_cart;
    unsigned int sparseIdx = 0;
    double *output_kernel_weights;
    mwIndex *output_cart_idx;
    unsigned int iSample, iNeighbor;
    mwIndex *output_nNeighbors_copy;
            
    struct ThreadData *threadData;
    pthread_t *thread;
    
    /* INPUT 0 - NONCARTESIAN COORDINATES */
    dims_nonCart = mxGetDimensions(prhs[0]);
    nPts_nonCart = dims_nonCart[1];
    nDims = dims_nonCart[0];
    coords_nonCart = (double*) mxGetData(prhs[0]);
    
    /* INPUT 1 - CARTESIAN OUTPUT MATRIX SIZE */
    dims_cart = (unsigned int*)mxGetData(prhs[1]);
            
    /* INPUT 2 - KERNEL SIGMA */
    kernel_sigma = (double*)mxGetData(prhs[2]);
    
    /* INPUT 3 - KERNEL EXTENT */
    kernel_extent = (double*)mxGetData(prhs[3]);
    
    /* INPUT 4 - NTHREADS */
    nthreads = (unsigned int*)mxGetData(prhs[4]);
    
    // Allocate memory
    idx_convert = calloc(nDims,sizeof(unsigned int));
    inv_kern_mult = calloc(nDims, sizeof(double));
    threadData = calloc(*nthreads,sizeof(struct ThreadData));
    thread = calloc(*nthreads,sizeof(pthread_t));
    half_kernel_extent = calloc(nDims, sizeof(double));
    
    // Create index conversion
    idx_convert[0] = 1;
    maxNonSparsePerSample = ceil(kernel_extent[0]);
    inv_kern_mult[0] = -0.5/(kernel_sigma[0]*kernel_sigma[0]);
    half_kernel_extent[0] = 0.5*kernel_extent[0];
    for(dim=1; dim<nDims; dim++){
        idx_convert[dim] = idx_convert[dim-1]*dims_cart[dim-1];
        maxNonSparsePerSample = maxNonSparsePerSample*(ceil(kernel_extent[dim])); // maybe add one? - This should be less sparse!
        inv_kern_mult[dim] = -0.5/(kernel_sigma[dim]*kernel_sigma[dim]);
        half_kernel_extent[dim] = 0.5*kernel_extent[dim];
    }
    nPts_cart = idx_convert[nDims-1]*dims_cart[nDims-1];
    
    maxNonSparseTotal = maxNonSparsePerSample*nPts_nonCart;
   
    //allocate temporary matrices
    temp_output_kernel_weights = calloc(maxNonSparseTotal, sizeof(double));
    temp_output_cart_idx = calloc(maxNonSparseTotal, sizeof(mwIndex));
    output_nNeighbors = calloc(nPts_nonCart, sizeof(mwIndex));
    nonSparseTotalCount = 0;
    
    numPtsPerThread = floor(((double)nPts_nonCart)/(*nthreads));
//     printf("CREATING %d THREADS with %d PTS TO GRID EACH\n",*nthreads,numPtsPerThread);
    for(iThread=0; iThread<*nthreads; iThread++){
        threadData[iThread].startSample = iThread*numPtsPerThread;
        threadData[iThread].endSample = (iThread+1)*numPtsPerThread-1;
        threadData[iThread].nonSparseVals = calloc(1,sizeof(unsigned int));
    }
    threadData[*nthreads-1].endSample = nPts_nonCart-1;
    
    for(iThread=0; iThread<*nthreads; iThread++){
//         printf("   STARTING THREADS %d\n",iThread);
        pthread_create(&thread[iThread], NULL, go, &threadData[iThread]);
    }
    
    /* Wait for threads to finish */
    for(iThread=0; iThread<*nthreads; iThread++){
        pthread_join(thread[iThread],NULL);
        nonSparseTotalCount = nonSparseTotalCount+(*threadData[iThread].nonSparseVals);
        
        //Free up memory
        free(threadData[iThread].nonSparseVals);
//         printf("   THREADS %d FINISHED\n",iThread);
    }
    
    /* Free up memory */
//     printf("FREEING MEMORY\n");
    free(idx_convert);
    free(threadData);
    free(thread);
    free(inv_kern_mult);
    
    /* Calculate actual sparce matrix */
//     printf("Creating sparse matrix of size %dx%d with %d nonsparse values (%f sparse)\n",nPts_cart,nPts_nonCart,nonSparseTotalCount,((double)nonSparseTotalCount)/((double)nPts_cart*nPts_nonCart));
    plhs[0] = mxCreateSparse(nPts_cart, nPts_nonCart, nonSparseTotalCount, mxREAL);
    output_kernel_weights  = (double*)mxGetData(plhs[0]);
    output_cart_idx        = mxGetIr(plhs[0]);
    output_nNeighbors_copy = mxGetJc(plhs[0]);
    
    output_nNeighbors_copy[0] = 0;
    for(iSample=0; iSample<nPts_nonCart;iSample++){
        output_nNeighbors_copy[iSample+1] = output_nNeighbors_copy[iSample]+output_nNeighbors[iSample];
//         printf("output_nNeighbors_copy[%d]=%d\n",iSample,output_nNeighbors_copy[iSample]);
        //THis could be done in main loop
    }
    
    /* Copy everything into the sparse array */
//     printf("Copying to sparse matrix...\n");
    sparseIdx = 0;
    for(iSample=0; iSample<nPts_nonCart; iSample++){
        /* Loop through each neighbor, skipping spaces */
        for(iNeighbor=0; iNeighbor<output_nNeighbors[iSample]; iNeighbor++){
            output_kernel_weights[sparseIdx] = temp_output_kernel_weights[iSample*maxNonSparsePerSample + iNeighbor];
            output_cart_idx[sparseIdx] = temp_output_cart_idx[iSample*maxNonSparsePerSample + iNeighbor];
            sparseIdx++;
            
//             *output_kernel_weights = temp_output_kernel_weights[iSample*maxNonSparsePerSample + iNeighbor];
//             output_kernel_weights++; // increment to next index
        }
    }    
    
    // Free up memory
    free(temp_output_kernel_weights);
    free(temp_output_cart_idx);
    free(output_nNeighbors);
}