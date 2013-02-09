#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>

typedef std::string String;
typedef std::vector< String > Strings;


Strings split(const String &s, const char delimiter) {
  std::istringstream ss(s);
  String item;
  Strings words;
  while(std::getline(ss, item, delimiter)) {
    if(! item.empty()) words.push_back(item);
  }
  return words;
}


template <typename T>
class Matrix {

private:
    unsigned int rows_;
    unsigned int cols_;
    T *matrix_;
    String name_;

public:

    Matrix(unsigned int rows, unsigned int cols): rows_(rows), cols_(cols) {
        matrix_ = new T[rows*cols];
    };

    ~Matrix() {
        delete matrix_;
    }

    void set_name(String name) {name_ =name;}



    unsigned int rows() const {return rows_;}
    unsigned int cols() const {return cols_;}

    T &operator()(unsigned int i, unsigned int j) {
        return matrix_[i * rows_ + j];
    }

    T &operator()(const std::pair<unsigned int, unsigned int> &p) {
        return matrix_[p.first * rows_ + p.second];
    }
};



template <typename T>
class SubstitutionScore {
public:

    T insertion_score()  const { return 1;}
    T match_score() const { return 100;}
    T mismatch_score() { return -5000;}
};

String global_aligment(const String &A, const String &B) {

    SubstitutionScore<int> score;
    Matrix<int> M(A.size()+1, B.size()+1);
    typedef std::pair<int,int> Index;
    Matrix< Index > traceback(A.size()+1, B.size()+1);
    // Initialization conditions
    M(0,0) = 0;
    for(unsigned int i = 1; i <= A.size(); ++i) {
        M(i,0) = M(i-1,0) + score.insertion_score();
    }
    for(unsigned int j = 1; j <= B.size(); ++j) {
        M(0,j) = M(0,j-1) + score.insertion_score();
    }
    for(unsigned int i = 1; i <= A.size(); ++i) {
        for(unsigned int j = 1; j <= B.size(); ++j) {
            assert((i-1)>=0);
            assert((j-1)>=0);
            if(A[i-1] == B[j-1] ) {
                Index ind = std::make_pair(i-1,j-1);
                M(i,j) = score.match_score() + M(ind);
                traceback(i,j) = ind;
            } else {
                Index ind1 = std::make_pair(i-1,j);
                int val1 = score.insertion_score() + M(ind1);
                Index ind2 = std::make_pair(i,j-1);
                int val2 = score.insertion_score() + M(ind2);
                if (val1 > val2) {
                    M(i,j) = val1;
                    traceback(i,j) = ind1;
                } else {
                    traceback(i,j) = ind2;
                    M(i,j) = val2;
                }
            }
        }
    }

    // recover sequence from the alignment
    String longest_subsequence = "";
    Index p = std::make_pair(A.size(), B.size());
    while( p.first != 0 && p.second != 0) {
        assert((p.first-1)>=0);
        assert((p.second-1)>=0);
        char c = A[p.first - 1];
        char d = B[p.second - 1];
        Index a =  traceback(p);
        if( A[p.first - 1 ] == B[p.second - 1]) {
            longest_subsequence.append(1,c);
        }
        p = traceback(p);
    }
    std::reverse(longest_subsequence.begin(), longest_subsequence.end());
    return longest_subsequence;
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
        Strings words = split(line, ';');
        String longest_subsequence;
        longest_subsequence = global_aligment(words[0],words[1]);
        cout << longest_subsequence << endl;
    }
    exit(0);
}
