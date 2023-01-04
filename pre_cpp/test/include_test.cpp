#include <iostream>
using namespace std;
int test(int i){
    return i+1;
}

class Test{
public:
    int testtest;
    int print();
    Test(int);
};

Test::Test(int tt=1){
    testtest = tt;
}

int Test::print(){
    cout << testtest << endl;
}