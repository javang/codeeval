
#include <stdexcept>
#include <functional>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>


typedef std::string String;
typedef std::vector<String> Strings;
typedef std::vector<int> Ints;


Strings split(const String &s, const char delimiter) {
  std::istringstream ss(s);
  String item;
  Strings words;
  while(std::getline(ss, item, delimiter)) {
    if(! item.empty()) words.push_back(item);
  }
  return words;
}

namespace dst {




// Node for a binary tree
template <typename T>
class TreeNode {
    
private:
    T value_;
    TreeNode<T> *right_;
    TreeNode<T> *left_;
    TreeNode<T> *parent_;
    
    void swap_pointer(TreeNode<T> *t, TreeNode<T> *p) {
        TreeNode<T> *temp = t;
        t = p;
        p = temp;
    }

public:

    TreeNode(T value, TreeNode *right=NULL, TreeNode *left=NULL):
         value_(value), right_(right), left_(left), parent_(NULL) {}

    void delete_left() {
        delete right_;    
        right_ = NULL;
    }
       
    void delete_right() {
        delete left_;
        left_ = NULL;
    }

    void set_right(const T &value) {
        if(right_ != NULL) delete_right();
        right_ = new TreeNode(value);
        right_->set_parent(this);
    }

    void set_right(TreeNode<T> *t) {
        if(right_ != NULL) delete_right();
        right_ = t;
        right_->set_parent(this);
    }

    void set_parent(TreeNode <T> *t) {
        parent_ = t;
    }
    void set_left(const T &value) {
        if(left_ != NULL) delete_left();
        left_ = new TreeNode(value);
        left_->set_parent(this);
    }

    void set_left(TreeNode<T> *t) {
        if(left_ != NULL) delete_left();
        left_ = t;
        left_->set_parent(this);
    }

    void swap(TreeNode<T> *t) {
        std::swap(value_,t->value);
        swap_pointer(right_, t->right_);
        swap_pointer(left_, t->left_);
    }
    
    TreeNode<T> *left() {
        return left_;
    }

    TreeNode<T> *parent() {
        return parent_;
    }

    TreeNode<T> *right() {
        return right_;
    }

    T get_value() const { return value_; }
    
    bool is_leaf() const {
        if( left_ == NULL && right_ == NULL) return true;
        return false;
    }

    void add_child(TreeNode<T> *t) {
        if (left_ == NULL) {
            set_left(t);
        } else if (right_ == NULL) {
            set_right(t);
        } else {
            throw std::length_error("It is not possible to add more children");
        }
    }

    //! true if the node has two children already
    bool is_full() const {
        if (left_==NULL ||right_==NULL ) return false;
        return true;
    }

    void clear() {
        if(left_ != NULL) {
            left_->clear(); 
            delete_left();
            left_ = NULL;
        }
        if(right_ != NULL) {
            right_->clear(); 
            delete_right();
            right_ = NULL;
        }
    }
};
}

double plus(double x, double y) {return x+y;}
double minus(double x, double y) {return x-y;}
double multiply(double x, double y) {return x*y;}
double divide(double x, double y) {return x/y;}

struct Token {
private:
    double number_;
    bool is_number_;
    std::pointer_to_binary_function<double, double, double> binary_operation_;     

public:
    Token(const String &x) {
        is_number_ = false;
        if( x == "*") {
            binary_operation_ = std::pointer_to_binary_function<double, double, double>(multiply);
        } else if ( x == "/") {
            binary_operation_ = std::pointer_to_binary_function<double, double, double>(divide);
        } else if( x == "-") {
            binary_operation_ = std::pointer_to_binary_function<double, double, double>(minus);
        } else if( x == "+") {
            binary_operation_ = std::pointer_to_binary_function<double, double, double>(plus);
        } else {
            number_ = std::atof(x.c_str());        
            is_number_ = true;
        }
    }

    double evaluate( double x, double y) {
        return binary_operation_(x,y);
    }

    double evaluate() {
        return number_;
    }
    
    bool is_number() const {
        return is_number_;
    } 
};

typedef dst::TreeNode<Token> TokenNode; 

// Converts a set of strings forming a operation in prefix notation into an evaluation tree
TokenNode *build_tree(const Strings &words) {
    if(words.size() == 0) {
        throw std::length_error("There are no tokens to process");    
    }
    TokenNode *current = NULL;
    TokenNode *root = NULL;
    for (unsigned int i = 0; i < words.size(); ++i) {
        Token token(words[i]);
        TokenNode *node = new TokenNode(token);
        if(root == NULL) {
            root = node;
            current = root;
        } else {
            while(current->is_full()) { 
                current = current->parent();
                if (current == NULL) { // root
                    throw std::invalid_argument("Wrong expression");
                }
            }
            current->add_child(node);                    
            if( ! token.is_number()) {
                current = node;
            }
        }
    }
    return root;
}

// Evaluate the tree formed by the expression
double evaluate_expression(TokenNode *node) {
    if(node->is_leaf()) return node->get_value().evaluate();
    return node->get_value().evaluate(evaluate_expression(node->left()),
                                     evaluate_expression(node->right()) );
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
