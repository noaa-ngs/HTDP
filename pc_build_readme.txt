Install MYS2 (https://www.msys2.org/).

Install both 64- and 32-bit (https://www.msys2.org/wiki/MSYS2-installation/).  

Install following package for 32-bit:  mingw-w64-i686-toolchain

Details on package management at https://www.msys2.org/wiki/Using-packages/.


Assume here that msys2 installed at C:/msys2/ 
and mingw32 installed at C:/msys2/mingw2/bin/


Once installed with packages, launch msys2.exe
(at C:/msys2/, or create shortcut).


Change path to location of source code (this path just an example):
[user computer] MSYS ~
$ cd /c/temp/htdp

Set system path to compiler (must go to where msys2 and mingw32 installed):
[user computer] MSYS /c/temp/htdp
$ export PATH=$PATH:/c/msys2/mingw32/bin

Run compiler with following command:
[user computer] MSYS /c/temp/htdp
$ gfortran -DNGS_PC_ENV -static -m32 -Wall -Wtabs *.f -o htdp