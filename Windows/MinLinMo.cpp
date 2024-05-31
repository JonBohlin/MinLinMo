#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <arm_neon.h>
#include <thread>
#include <chrono>
#include <mutex>
#include <atomic>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>
#include "MinLinMo.hpp"
#include "VariableSelector.hpp"

// Main MinLinMo class

int main(int argv, char* args[]){
    // Default correlation values
    
    float cor=0.1f, R2=0.01f, resCor=0.1f;
    std::string matrixFile, outcomeFile, outputFile;
    
    std::cout<<"============="<<std::endl;
    std::cout<<"MinLinMo v1.1"<<std::endl;
    std::cout<<"============="<<std::endl;

    if( argv <= 1 ){
        std::cout<<"\nAdditional arguments required. Use -h for help."<<std::endl;

        exit( 0 );
    }

    // Activate CMD user interface
    CMDUI ui(argv, args);

    outcomeFile = ui.y_file;
    matrixFile = ui.X_file;
    cor = ui.r1;
    R2 = ui.R2;
    resCor = ui.r2;

    if( ui.output_results ){
        universal::output_file = ui.output_file;
        universal::output_results = true;
    }
        
    std::cout<<"Outcome vector filename entered:"<<outcomeFile<<std::endl;
    std::cout<<"Matrix filename entered:"<<matrixFile<<std::endl;
    std::cout<<"Correlation entered:"<<cor<<std::endl;
    std::cout<<"Delta Rsq entered:"<<R2<<std::endl;
    std::cout<<"Residual correlation entered:"<<resCor<<std::endl;
    std::cout<<"Results output file entered:"<<universal::output_file<<std::endl;

    // Load data with DataLoader class
    DataLoader dl(outcomeFile, matrixFile);
    // Perform variable selection and linear modeling with the VariableSelector class.
    VariableSelector vs(cor, R2, resCor);
    
    return 0;
}
