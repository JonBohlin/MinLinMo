class VariableSelector{
    private:
        // Global variables
        float ybar=0.0f, sumydiffsq=0.0f;
        float cor_threshold, R2_threshold, resCor_threshold;
        const int startIndex = 0;
        size_t allColumns = universal::matrix[0].size();
        size_t allRows = universal::matrix.size();
        std::vector<float> ydiff, ydiffsq;
        std::mutex mtx;

        // Element in the priority queue
        struct corIndexElement{
            float cor;
            int index;
        };

        struct corIndexElement corListElement;
        std::vector<corIndexElement> corList, sortedCorList;

        // Sort priority queue descendingly
        std::vector<corIndexElement> sort_correlation_list(std::vector<corIndexElement> list){
            std::sort(list.begin(), list.end(), [this](const corIndexElement& e1, const corIndexElement& e2){
                return (e1.cor > e2.cor);});
            return list;
        }

        // Apply submatrices for multi-thread handling
    public:
    // Constructor
    VariableSelector(float cor, float R2, float resCor){
        cor_threshold = cor;
        R2_threshold = R2;
        resCor_threshold = resCor;
        std::cout<<"\nCalculating correlation matrix"<<std::endl;
        calculate_response_stats(); // Pre-calculate all variables needed from outcome to do correlations
        // Time Pearson correlations
        auto s2 = std::chrono::high_resolution_clock::now();
        compute_correlations( universal::matrix );
        auto e2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> d2 = e2 - s2;
        std::cout<<"Time taken (sec/millisec):"<<d2.count()<<std::endl;
        std::cout<<"\nNumber of predictors correlating with outcome:"<<sortedCorList.size()<<std::endl;
        // Time model fitting
        auto s3 = std::chrono::high_resolution_clock::now();
        fit_model( universal::matrix );
        auto e3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> d3 = e3 - s3;
        std::cout<<"Time taken (sec/millisec):"<<d3.count()<<std::endl;
    }

    private:
    // Pre-calculate variables to be used in Pearson correlation phase
    void calculate_response_stats( void ){
        for(int i = 0; i < allRows; i++)
            ybar += universal::outcome[i];
        ybar/=allRows;

        for(int i=0; i < allRows; i++){
            ydiff.push_back( universal::outcome[i] - ybar );
            ydiffsq.push_back( ydiff[i] * ydiff[i] );
            sumydiffsq += ( ydiff[i] * ydiff[i] );
        }        
    }

        // SIMD part of Pearson correlations
    void compute_avx_correlation( Submatrix &predMatrix, const int offset ){
        int currCounter, prevCounter=0, roundCount;
        int predMatrixIndex;
        int rowSize = predMatrix.numRows();
        int colSize = predMatrix.numCols();
        corIndexElement CI;

        const int numLanes = 8; // Each lane is a 32 bit floating point
        int avxCycles = colSize / numLanes;
        const int avxRestCycle = colSize % numLanes;

        // Add one to counter if rest not zero
        if(avxRestCycle != 0)
            avxCycles++;

        // All SIMD variables start with underscore
        __m256 _x, _y, _xbar, _xdiv, _xdiff, _sumyx, _sumyxsq, _rootsumyxsq, _rho;
        __m256 _xdiffsq, _ydiff, _ydiffsq, _sumydiffsq;

        _xdiv = _mm256_set1_ps( (float) rowSize );
        _x = _mm256_set1_ps( 0.0f );
        _xbar = _mm256_set1_ps( 0.0f );
        _sumyx = _mm256_set1_ps( 0.0f );
        _sumyxsq = _mm256_set1_ps( 0.0f );
        _sumydiffsq = _mm256_set1_ps( sumydiffsq );

        for(int k=0; k < avxCycles; k++){
            // Need to calculate _xbar to proceed.
            for(int i=0; i < rowSize; i++){
                for(int j=0; j < numLanes; j++){
                    predMatrixIndex = k * numLanes + j;
                    if( predMatrixIndex >= colSize ){
                        predMatrixIndex--;
                        break;
                    }
                    _x[j] = predMatrix(i, predMatrixIndex);
                }
                _xbar = _mm256_add_ps( _xbar, _x );
            }

            _xbar = _mm256_div_ps( _xbar, _xdiv );

            // Calculate Pearson correlation between outcome and predictors number of lanes at a time
            for(int i=0; i < rowSize; i++){
                for(int j=0; j < numLanes; j++){
                    predMatrixIndex = k * numLanes + j;
                    if( predMatrixIndex >= colSize ){
                        predMatrixIndex--;
                        break;
                    }
                    _x[j] = predMatrix(i, predMatrixIndex);
                }
                _ydiff = _mm256_set1_ps( ydiff[i] );
                _xdiff = _mm256_sub_ps( _x, _xbar );
                _xdiffsq = _mm256_mul_ps( _xdiff, _xdiff );
                _sumyx = _mm256_fmadd_ps( _ydiff, _xdiff, _sumyx);  // a = Sum_i (y_i - ybar)(x_i - xbar)
                _sumyxsq = _mm256_fmadd_ps(_sumydiffsq, _xdiffsq, _sumyxsq);    // b = Sum_i (y_i - ybar)^2 Sum_i (x_i - xbar)^2
            }

            _rootsumyxsq = _mm256_sqrt_ps( _sumyxsq );  // rho = (a / sqrt(b)) = Pearson correlation
            _rho = _mm256_div_ps( _sumyx, _rootsumyxsq );
            _x = _mm256_set1_ps( 0.0f );
            _xbar = _mm256_set1_ps( 0.0f );
            _sumyx = _mm256_set1_ps( 0.0f );
            _sumyxsq = _mm256_set1_ps( 0.0f );
            _sumydiffsq = _mm256_set1_ps( sumydiffsq );

            currCounter = predMatrixIndex + 1;
            roundCount = currCounter - prevCounter;
            
            for(int j=0; j < roundCount; j++){
                if( fabs(_rho[j]) > cor_threshold ){
                    CI.cor = fabs(_rho[j]);
                    CI.index = offset + (k*numLanes + j + 1);
                    std::unique_lock<std::mutex> lock(mtx); // Need this to avoid race conditions
                    corList.push_back( CI );
                }
            }
            prevCounter = currCounter;
        }
    }

    // Compute Pearson correlation between outcome and predictors
    void compute_correlations( universal::FloatMatrix &matrix ){
        std::vector<std::thread> threads;
        int offset=0;
        // Matrix divided up into section corresponding to number threads available
        int numIntervals = std::thread::hardware_concurrency();
        int subMatrixSize = (int) (allColumns / numIntervals);
        std::vector<int> offsets;
        int rest = (allColumns % numIntervals);
        std::vector<Submatrix> submatrices;
        
        for(int i=0; i < numIntervals; i++){
            submatrices.push_back( Submatrix(matrix, 0, i*subMatrixSize, matrix.size(), (i+1) * subMatrixSize) );
            offsets.push_back( (i * subMatrixSize) - 1 );
        }

        submatrices.push_back( Submatrix(matrix, 0, (numIntervals * subMatrixSize), matrix.size(), numIntervals * subMatrixSize + rest) );
        offsets.push_back( numIntervals * subMatrixSize - 1 );

        // One thread for each submatrix
        for( int i=0; i < submatrices.size(); i++ )
            threads.push_back(std::thread(&VariableSelector::compute_avx_correlation, this, std::ref(submatrices[i]), offsets[i]));

        for(auto& thr : threads)
            thr.join();

        // Priority queue
        sortedCorList = sort_correlation_list( corList );
        
        // If nothing to do quit
        if( sortedCorList.size() == 0 ){
            std::cout<<"No predictors found to correlate with response"<<std::endl;
            exit( 0 );
        }     
    }

    // Calculate adjusted (with respect to degrees of freedom (DoF) coefficient of variation
    double calculate_r2(const gsl_vector *res, const gsl_vector *y, const double DoF){
        double SS_res = 0.0, SS_tot=0.0, R2=0.0, adjR2=0.0;
        // Calculate y mean
        for(int i = 0; i < y->size; i++)
            SS_res += ( res->data[i] * res->data[i] );
        
        for(int i = 0; i < y->size; i++)
            SS_tot += ((y->data[i] - ybar) * (y->data[i] - ybar));

        R2 = ( 1 - SS_res/SS_tot );
        adjR2 = (1.0 - (1.0 - R2)*((res->size - 1.0)/(res->size - (DoF + 1.0))));
        return adjR2;
    }

    // Calculate Pearson correlation between residuals and predictors
    double pearson_correlation(const gsl_vector *res, const gsl_vector *x){
        double resbar=0.0, xbar=0.0, numerator=0.0, denumeratorRes=0.0, denumeratorX=0.0;
        // Calculate the averages of res and x

        for(int i = 0; i < res->size; i++){
            resbar+= res->data[i];
            xbar += x->data[i];
        }
        resbar /= res->size;
        xbar /= x->size;

        for(int i = 0; i < res->size; i++){
            numerator += (res->data[i] - resbar) * (x->data[i] - xbar);
            denumeratorRes += (res->data[i] - resbar) * (res->data[i] - resbar);
            denumeratorX += (x->data[i] - xbar) * (x->data[i] - xbar);
        }
        return numerator / sqrt(denumeratorRes * denumeratorX );
    }

    // Main model fitting function
    void fit_model(universal::FloatMatrix predMatrix){
        std::cout<<"\nComputing model..."<<std::endl;
        std::vector<float> slope;
        const double cor = 0.1;
        int numColsToInclude = 2, destIndex = 0;
        int tau_size, status, predIndex=0;
        double adjR2 = 0.0, prev_adjR2 = 0.0, last_adjR2 = 0.0, resCor = 0.0;
        std::vector<corIndexElement> significantPredictors;
        gsl_matrix *qrmatrix, *tempMatrix;
        gsl_vector *tau, *res, *resid, *beta, *betahat, *predVec, *yy, *tempVec;
        corIndexElement element;
        bool leaveUnchainged = false;
        tempVec = gsl_vector_alloc( allRows );
        yy = gsl_vector_alloc( allRows );
        res = gsl_vector_alloc( allRows );
        resid = gsl_vector_alloc( allRows );
        predVec = gsl_vector_alloc( allRows );
        
        for(int i=0;i < allRows; i++){
            gsl_vector_set(predVec, i, 1.0); // Intercept
            gsl_vector_set(yy, i, (double) universal::outcome[i]); // The response in gsl vector format; must convert to double
        }

        /*  Loop to perform all regressions using the following pseudo-code:
            
            1. Have a matrix ready with the first predictor, which is the intercept
            2. Add a vector to that matrix and the first predictor
            3. Perform linear regression using QR decomposition with y as response
            4. If added predictor is significant add to matrix; if not, overwrite with next predictor
        */

        // The first QR matrix consists of the intercept and the first predictor
        qrmatrix = gsl_matrix_alloc( allRows, numColsToInclude);
        tempMatrix = gsl_matrix_alloc( allRows, numColsToInclude);
        gsl_matrix_set_col(qrmatrix, destIndex, predVec); // Intercept

        // destIndex points to next predictor, predIndex to highest correlating element in priority queue
        destIndex++;
        element = sortedCorList[ predIndex ]; // Get most correlating predictor from priority queue (try "pop")

        // Copy element from matrix into gsl type vector/matrix
        for(int i=0;i < allRows; i++)
            predVec->data[i] = (double) universal::matrix[i][element.index]; // Get column index to predictor from priority queue

        gsl_matrix_set_col(qrmatrix, destIndex, predVec); // First, highest correlating predictor
        
        // tau_size is the smalles dimension of either matrix rows or columns
        if( qrmatrix->size1 > qrmatrix->size2 )
            tau_size = (int) qrmatrix->size2;
            else tau_size = (int) qrmatrix->size1;
        
        tau = gsl_vector_alloc( tau_size );
        // regression coefficients, always nubmer of columns/predictors
        beta = gsl_vector_alloc( tau_size );
        betahat = gsl_vector_alloc( tau_size );

        gsl_matrix_memcpy(tempMatrix, qrmatrix); // Need a copy predictor matrix before qrmatrix is overwritten after QR decomposition

        status = gsl_linalg_QR_decomp( qrmatrix, tau); // Perform new QR decomposition
        status = gsl_linalg_QR_lssolve(qrmatrix, tau, yy, beta, res); // Use previous QR decomposition in least sqaure estimates

        // Results from initial regression model
        adjR2 = calculate_r2(res, yy, destIndex); // Calculate adj. R2 - destIndex represents Degrees of Freedon (DoF)
        significantPredictors.push_back( element ); // Might as well take note of the most correlating predictor, if model is too "weak" it won't be reported anyway.
        last_adjR2 = adjR2;

        gsl_vector_memcpy( resid, res );
        gsl_vector_memcpy( betahat, beta );

        predIndex++; // Get next predictor, this will only be correlated against the residuals of the previous model

        if(adjR2<R2_threshold){
            std::cout<<"No predictor correlates sufficiently with outcome."<<std::endl;
            return;
        }
        else{
            do{
                // "Pop" element from priority queue
                element = sortedCorList[ predIndex ];
                for(int i=0;i < allRows; i++)
                    predVec->data[i] = (double) predMatrix[i][element.index]; // Get column index to predictor from priority queue
                resCor = pearson_correlation( resid, predVec ); // Correlate against residuals so we won't have to do a full QR decomposition
            
                if( fabs(resCor) >= resCor_threshold ){
                    // Need to copy elements from tempMatrix into new largermatrix
                    // Get rid of the old vectors and matrix
                    gsl_vector_free( tau );
                    gsl_vector_free( beta );
                    gsl_matrix_free( qrmatrix );
                    numColsToInclude++;
                    destIndex++; // Predictor correlates with residuals, increase model matrix
                    
                    // Allocate memory for new vectors and matrix
                    if( numColsToInclude > res->size)
                        tau_size = (int) res->size;
                    else tau_size = numColsToInclude;

                    tau = gsl_vector_alloc( tau_size );
                    beta = gsl_vector_alloc( tau_size );
                    qrmatrix = gsl_matrix_alloc( allRows, numColsToInclude );

                    for(int j=0; j < numColsToInclude - 1; j++){
                        gsl_matrix_get_col( tempVec, tempMatrix, j); 
                        gsl_matrix_set_col( qrmatrix, j, tempVec );
                    }

                    gsl_matrix_set_col( qrmatrix, numColsToInclude - 1, predVec );

                    // Can now free up old tempMatrix before a new larger one is allocated
                    if(!leaveUnchainged){
                        gsl_matrix_free( tempMatrix );
                        tempMatrix = gsl_matrix_alloc( allRows, numColsToInclude );
                    }

                    gsl_matrix_memcpy(tempMatrix, qrmatrix); // qrmatrix will hold Q and R matrices after QR decomposition, tempMatrix all predictors + intercept
                    
                    status = gsl_linalg_QR_decomp(qrmatrix, tau); // Perform QR decomposition
                    status = gsl_linalg_QR_lssolve(qrmatrix, tau, yy, beta, res); // Use previous QR decomposition in least sqaure estimates
                    
                    adjR2 = calculate_r2(res, yy, destIndex);
                    //std::cout<<"R2:"<<(adjR2 - last_adjR2)<<"\tIndex:"<<element.index<<std::endl;

                    if( adjR2 - last_adjR2 >= R2_threshold ){
                        // we're have to overwrite the last predictor in the model
                        gsl_vector_memcpy( resid, res );
                        gsl_vector_free( betahat );
                        betahat = gsl_vector_alloc( tau_size );
                        gsl_vector_memcpy( betahat, beta );
                        significantPredictors.push_back( element );
                        last_adjR2 = adjR2;
                        leaveUnchainged = false;
                    }
                    else {
                        // Indexes to model matrix stay the same, just over write the last column
                        numColsToInclude--;
                        destIndex--;
                        leaveUnchainged = true; // Don't spend time recreating tempMatrix
                    }
                }
                predIndex++;
            } while( predIndex < sortedCorList.size());

            std::cout<<betahat->data[0]<<' '<<"intercept"<<std::endl;
            int ind=1;
            for(auto &hd : significantPredictors){
                std::cout<<betahat->data[ind]<<' '<<universal::matrixHeaders[hd.index]<<std::endl;
                ind++;  
            }
            std::cout<<"\nFinal model R2:"<<last_adjR2<<std::endl; //destIndex = Degrees of freedom
            std::cout<<"\nNumber of predictors:"<<significantPredictors.size()<<std::endl;

            if( universal::output_results){
                std::ofstream outp( universal::output_file );
                std::cout<<"Outputting results to file:"<<universal::output_file<<std::endl;
                outp<<"est id"<<std::endl;
                outp<<betahat->data[0]<<' '<<"intercept"<<std::endl;
                int ind=1;
                for(auto &hd : significantPredictors){
                    outp<<betahat->data[ind]<<' '<<universal::matrixHeaders[hd.index]<<std::endl;
                    ind++;  
                }
                outp.close();
            }

            // Free up
            gsl_vector_free( tau );
            gsl_vector_free( yy );
            gsl_vector_free( res );
            gsl_vector_free( resid );
            gsl_vector_free( predVec );
            gsl_vector_free( tempVec );
            gsl_vector_free( beta );
            gsl_vector_free( betahat );
            gsl_matrix_free( qrmatrix );
            gsl_matrix_free( tempMatrix );
        }
    }
};
