macOS Build with gfortran https://gcc.gnu.org/
Tested on macOS Monterey 12.5

git clone https://github.com/noaa-ngs/HTDP
cd HTDP
gfortran -std=legacy -Wall -Wtabs *.f -o htdp
./htdp     