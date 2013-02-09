#include <string>
#include <iostream>


typedef std::string String;

template <typename T>
class Matrix {

private:
    unsigned int rows_;
    unsigned int cols_;

public:

    Matrix(unsigned int rows, unsigned int cols): rows_(rows), cols_(cols) {
        T *matrix_ = new T[rows*cols];
    };

    ~Matrix() {
        delete matrix_;
    }

    unsigned int rows() const {return rows_}
    unsigned int cols() const {return cols_}
    T &operator(unsigned int i, unsigned int j) {
        return &matrix_[i * rows_ + j];
};


template<typename T, class iterator>
T min(const iterator &begin, const iterator &end) {
    T minimum = std::nueric_limits<T>::max();
    for(iterator it = begin; it != end; ++it) {
        minimum = std::min(it, minimum);
}

class SubstitutionScore {
    dobule operator(char a, char b) {
        if( a = '-' or b = '-') return 0; // no cost for a gap
        if(a = b) return 1; //match
        return -1; // mismatch
    }
};


template<typename T>
String global_aligment(const String &A, const String &B) {
    Matrix M<T>(A.size()+1, B.size()+1);
    std::vector< std::pair<unsigned int,unsigned int> > traceback;
    // Initialization conditions
    for(unsigned int i = 0; i <= A.size(); ++i) M(i,0) = score(A[i],'-');
    for(unsigned int j = 0; j <= B.size(); ++j) M(0,j) = score('-',B[j]);
    for(unsigned int i = 1; i <= A.size(); ++i) {
        for(unsigned int j = 1; j <= B.size(); ++j) {
            std::pair p;
            M(i,j) = score(i,j) + M(i-1,j-1);
            std::pair p = std::make_pair(i-1,j-1);

            T val = score(i-1,j) + M(i-1,j);
            if(val > M(i,j) {
                M(i,j) = val;
                p = std::make_pair(i-1,j);
            }

            val = score(i,j-1) + M(i,j-1);
            if(val > M(i,j) {
                M(i,j) = val;
                p = std::make_pair(i,j-1);
            }
            traceback.pushback(p);
        }
    }

    // recover the best alignment from A
    Sttring max_aligned = "";
    for(unsigned int k = traceback.size() - 1; k >= 0; --k) {

    }
}












int main(int argc, char **argv) {
    using namespace std;

    if(argc != 2) {
        cout << "Parameters" << endl;
        cout << "[1] - Input file" << endl;
        exit(0);
    }

    std::ifstream f_in;
    f_in.open(argv[1]);
    if( ! f_in.good() ) {
        throw ios_base::failure("Error opening the file ");
    }

    String line = "";
    while(std::getline(f_in, line)) {
        if(line.empty()) continue;
        Strings words = split(line, ' ');
        TokenNode *node = build_tree(words);
        std::cout << evaluate_expression(node) << std::endl;
        node->clear();
    }
    exit(0);
}
