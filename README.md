The program runs a MC simulation following the work of Perassi and Leiva Electrochemistry Communications 65 (2016) 48–52. Several output files are generated, according to the user's needs. A .dat file with the results, a .xyz file to display the coordinates of the lithium ions that have been intercalated in the lattice sites and another .xyz file that prints the coordinates of all the intercalation sites.
There are notes commenting on the simulation parameters, the outputs and the different parts of the code. I have left a couple of notes in Spanish, accidentally.

The code is set to generate a square lattice in the x,y planes and 4 graphite sheets in the z axis. By default it generates a 12x12x4 lattice sites.
As I have always run the code in Windows, I have not yet tried to do it in Linux. But if you use Linux, you could try using:

g++ -O3 -o "name of the executable" MC-LiC-UVT.cpp

-O3 is a flag that optimizes some things in the code when compiling, I recommend testing how it works with -O3 and without -O3. It may be faster with the flag but I don't know if it causes some error in the code. It would be a matter of testing.
I also recommend testing how the random number generator works in Linux. If necessary you could use another one.

Also, in case you want to use Windows, you can use CodeBlocks (https://www.codeblocks.org/) which is an easy-to-use program to compile and run C++
