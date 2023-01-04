#include <iostream>
#include "include_test.cpp"
using namespace std;
int main(){
    cout << test(1) << endl;
    Test t(2);
    t.print();
    return 0;
}