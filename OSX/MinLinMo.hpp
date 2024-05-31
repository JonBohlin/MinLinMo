// Global variables
namespace universal{
    typedef std::vector<float> FloatRow;
    typedef std::vector<FloatRow> FloatMatrix;
    static FloatMatrix matrix;
    static std::vector<float> outcome;
    static std::vector<std::string> matrixHeaders;
    static std::string outcomeHeader;
    static std::string output_file;
    static bool output_results = false;
}

// Matrix class that allows for submatrix interrogations
class Submatrix{
    universal::FloatMatrix &matrix;
    size_t startRow, startCol, endRow, endCol;

    public:
    // Submatrix constructor initialising variables directly
    Submatrix(universal::FloatMatrix &matrix, size_t startRow, size_t startCol, size_t endRow, size_t endCol)
        :matrix(matrix), startRow(startRow), startCol(startCol), endRow(endRow), endCol(endCol) {
    }

    // Allow for "natural" matrix access
    float& operator()(size_t row, size_t col) {
        if (row + startRow >= endRow || col + startCol >= endCol) {
            throw std::out_of_range("Index out of range for submatrix.");
        }
        float &temp = matrix[startRow + row][startCol + col];
        return temp;
    }

    size_t numRows() const {
        return endRow - startRow;
    }
    
    size_t numCols() const {
        return endCol - startCol;
    }

    // Allow for streaming of output
    friend std::ostream& operator<<(std::ostream& os, Submatrix& sm) {
        for (size_t i = 0; i < sm.endRow - sm.startRow; ++i) {
            for (size_t j = 0; j < sm.endCol - sm.startCol; ++j) {
                os << sm(i, j) << " ";
            }
            os << "\n";
        }
        return os;
    }
};

// Class for loading data: outcome and matrix
class DataLoader{
    public:

    // Constructor
    DataLoader(std::string ofile, std::string mfile){
        text_to_vector( ofile );
        if( &universal::outcome!=nullptr)
            std::cout<<"\nLoaded outcome:\nNumber of rows:"<<universal::outcome.size()<<std::endl;
            else{
                std::cout<<"Couldn't load outcome."<<std::endl;
                exit( 1 );
            }

        csv_to_float_matrix( mfile );

        if( &universal::matrix!=nullptr)
            std::cout<<"Loaded matrix data:\nNumber of rows:"<<universal::matrix.size()<<"\tNumber of columns:"<<universal::matrix[0].size()<<std::endl;
            else{
                std::cout<<"Matrix not defined"<<std::endl;
                exit(1);
            }
        
        if( universal::outcome.size() != universal::matrix.size() ){
            std::cout<<"Number of outcome rows:"<<universal::outcome.size()<<std::endl;
            std::cout<<"Number of matrix rows:"<<universal::matrix.size()<<std::endl;
            std::cout<<"...don't match"<<std::endl;
            exit( 1 );
        }
    }

    private:

    // Load predictor matrix
    void csv_to_float_matrix(const std::string& filename){
        std::ifstream file(filename);
        //FloatMatrix matrix;
        std::string word;

        if (file) {
            std::string line;
            // Read the column headers
            if (std::getline(file, line)) {
                std::stringstream ss(line);
                std::string colHeader;
                while (std::getline(ss, colHeader, ','))
                    universal::matrixHeaders.push_back(colHeader);
            }
            // Read the data rows
            while (std::getline(file, line)) {
                universal::FloatRow row;
                std::stringstream ss(line);
                std::string value;

                while (std::getline(ss, value, ',')) {
                    row.push_back(std::stof(value));
                }

                universal::matrix.push_back(row);
            }
        } else {
            std::cerr << "Failed to open file:" << filename << std::endl;
            exit(1);
        }
    }

    // Get outcome vector
    void text_to_vector(const std::string& filename ){
        float num;

        std::ifstream inputFile(filename);

        if (!inputFile.is_open()) {
            std::cerr << "Can't fetch outcome from file:" << filename << std::endl;
            exit(1);
        }

        std::getline(inputFile, universal::outcomeHeader);

        while (inputFile >> num)
            universal::outcome.push_back(num);

        inputFile.close();
    }
};

class CMDUI{
// Variables for internal usage
private:
    int numArguments = 0, command;
    float testValue;
    std::string enteredCommand, testArgument;

// Accepted commands, sort them according to the arguments they take (e.g. 3 floats, 3 strings, and nothing)
    std::vector<std::string> acceptedCommands {"-R2", "-r1", "-r2", "-O", "-X", "-y", "-h"};
    enum CMD {c_R2, c_r1, c_r2, c_O, c_X, c_y, c_h}; // enum of all commands

// Keep track whether a command line argument has been entered multiple times     
    std::vector<bool> enteredCommands { false, false, false, false, false, false, false };

// Numerical values representing the different commands
    int numericArguments=3, stringArgument=3, noArgument=4;

    int getCommand(std::string token ){
        int j;
        for( j=CMD::c_R2; j<=CMD::c_h;++j)
            if( token.compare( acceptedCommands[ j ] ) == 0 )
                return j;
        return j+1;
    }

    void parseInput(int argv, char* args[]){
        int i = 1;
        while(i < argv){
            enteredCommand = args[ i ];
// Are the provided command line arguments in the command dictionary?
            command = getCommand( enteredCommand );
            if( command > CMD::c_h ){
// If no, exit.                
                std::cout<<"Command not found"<<std::endl;
                exit( 1 );
            } else {
// Command found, accept first entry only
                    if( enteredCommands[ command ] ){
                        std::cout<<"Argument entered multiple times"<<std::endl;
                        exit( 4 );
                        break;
                    } else enteredCommands[ command ] = true;

                    if(command < CMD::c_h){
// These commands take one argument
                    if( i+1 < argv ){
// These commands take a float argument
                        if(command < CMD::c_O ){
                            try{
                                testValue = std::stof( args[i + 1] );
                                switch( command ){
                                    case CMD::c_R2: R2=testValue; break;
                                    case CMD::c_r1: r1=testValue; break;
                                    case CMD::c_r2: r2=testValue; break;
                                }
                            } catch( const std::invalid_argument& ){
                                std::cerr<<"Invalid float argument"<<std::endl;
                                exit( 2 );
                            } catch( const std::logic_error& ){
                                std::cerr<<"Argument not accepted"<<std::endl;
                                exit( 3 );
                            }
                            
                        } else if( command >= CMD::c_O && command <= CMD::c_y ){
// These commands takes a string argument
                            try{
                                testArgument = args[ i + 1 ];
                                
                                switch( command ){
                                    case CMD::c_O: output_results = true; output_file = testArgument; break;
                                    case CMD::c_X: X_file = testArgument; break;
                                    case CMD::c_y: y_file = testArgument; break;
                                }

                            } catch( std::logic_error& ){
                                std::cerr<<"Argument not accepted"<<std::endl;
                                exit( 3 );
                            }
                            
                        }
                    } else {
                        std::cout<<"Error entering argument"<<std::endl;
                        exit( 1 );
                    }
                    i++; // Increase entered commands by 1
                } else if( command >= CMD::c_h ){ 
// These commands don't take any arguments
                        switch( command ){
                            case CMD::c_h: 
                            std::cout<<"Created by Jon Bohlin, June 2024"<<std::endl;

                            std::cout<<"\nList of MinLinMo arguments:\n";
                            std::cout<<"---------------------------\n"<<std::endl;

                            std::cout<<"-h\tOverview of possible command line arguments"<<std::endl;
                            std::cout<<"-X\t<filename> filename to prediction matrix in .csv format - required"<<std::endl;
                            std::cout<<"-y\t<filename> filename to outcome vector in .csv format - required"<<std::endl;
                            std::cout<<"-O\t<filename> filename for tab-separated prediction results - optional"<<std::endl;
                            std::cout<<"-R2\t<float> (0.01) required model improvement for variable selection - optional"<<std::endl;
                            std::cout<<"-r1\t<float> (0.1) Pearson correlation threshold for admission to priority queue - optional"<<std::endl;
                            std::cout<<"-r2\t<float> (0.1) residual Pearson correlation threshold for model inclusion - optional"<<std::endl;
                            exit( 0 );
                            break;
                        }
                    }
            }
            i++;
        }
    }

public:
    // Variables that can be accessed by user
    
    // Default values if no command line argument provided
    float R2 = 0.01f, r1 = 0.1f, r2 = 0.1f;
    bool output_results = false;
    std::string output_file, X_file, y_file;

    CMDUI(int argv, char* args[]){
        std::cout<<std::endl;
        parseInput(argv, args);

    }
};
