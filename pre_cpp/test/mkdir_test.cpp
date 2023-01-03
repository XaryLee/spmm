#include <iostream>
#include <direct.h>
#include <string>
using namespace std;

int main(){
    string folderpath = "mkdir_test";
    string command = "mkdir " + folderpath;
    system(command.c_str());
    cout << "Success." << endl;
    return 0;
}