#include <string>

using namespace std;

class LinkInfo_entry {

 public:

  LinkInfo_entry() {};
  LinkInfo_entry(string seq1, string seq2, int direction);

  string get_seq1();
  string get_seq2();
  int get_direction(); 
 
 private:
  string _seq1;
  string _seq2;
  int _direction;
  
};

