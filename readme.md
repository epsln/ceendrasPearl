# Ceendra's Pearl
Implementation of Indra's Pearl depth first search in C. Yes, C. I followed the classic algorithm, using no particular data structures. As such, the code is a tad messy, but should be fairly understandable. 

## Compilation
Compile using `make`. This will also create an `out/` folder that will contain the images.

## Running
Run the `./ceendra` executable. This will create about 53 images of 1000x1000 resolution using random initial parameters. 

## Modifying the parameters
The general parameters are located in the `main.c` file. 

### The ta and tb seeds
ta and tb are the initial complex numbers used by the 'Grandma's formula' to create 2 generators and their inverse. The program will then compute the limit set of those generators. 
