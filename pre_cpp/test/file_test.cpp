#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main(int argc, char* argv[]){
    // char* path = argv[1];
    char path[10] = "hello.txt";
    ifstream fin;
    fin.open(path, ios::in);
    if(!fin.is_open())
        cout << "Failed.";
    string buf;
    while(getline(fin, buf))
        cout << buf << endl;
    fin.close();

    ofstream fout;
    fout.open(path);
    if(!fout.is_open())
        cout << "Failed.";
    fout << "Genshin Impact" << endl << "Hello world";
    fout.close();

    return 0;
}