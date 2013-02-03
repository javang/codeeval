#include <iostream>
#include <stdexcept>
#include <vector>
#include <sstream>
#include <fstream>
#include <ios>
#include <string>
#include <cstdlib>

typedef std::string String;
typedef std::vector<String> Strings;

Strings split(const String &s, const char delimiter) {
  std::istringstream ss(s);
  String item;
  Strings words;
  while(std::getline(ss, item, delimiter)) {
    if(! item.empty()) words.push_back(item);
  }
  return words;
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
        int n_words = words.size();
        int mth_index = atoi(words[n_words-1].c_str());
        if(mth_index > (n_words-1)) continue;
        cout << words[n_words-mth_index-1] << endl;
    }
    exit(0);

}