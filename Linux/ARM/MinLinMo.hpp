// Global variables
namespace universal{
    typedef std::vector<float> FloatRow;
    typedef std::vector<FloatRow> FloatMatrix;
    static FloatMatrix matrix;
    static std::vector<float> outcome;
    static std::vector<std::string> matrixHeaders;
    static std::string outcomeHeader;
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
            std::cout<<"Loaded outcome:\nNumber of rows:"<<universal::outcome.size()<<std::endl;
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
            std::cerr << "Failed to open file: " << filename << std::endl;
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
