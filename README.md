# Serial
## Version
The latest version of Serial is in PreProcessing folder.
## Make
Run the following commands to compile the code.
```
cd ./PreProcessing
make
```
## Run
Write all the names of matrices to test into matrix.txt.

e.g.
```
web-Stanford
web-Google
sx-askubuntu
```
Then put the matrices(.mtx) into the corresponding folder in ./PreProcessing/mat/mtx/<matrix name>

e.g.
```
mat/mtx/web-Stanford
```
Finally, run the following command to run the code.

```
./serial
```

Additionally, if needing to compute large matrices, the following command needs to be run to prevent stack overflow.
```
ulimit -s 65536
```
