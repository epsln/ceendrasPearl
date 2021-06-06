# Ceendra's Pearl
Implementation of Indra's Pearl depth first search in C. I followed the classic algorithm. As such, the code is a tad messy, but should be fairly understandable. I wrote an article underlying the basic principles behind it on my blog, which is available [here](https://epsln.github.io/blog/indraspearl_pt1/)

## Compilation

Compile using `make`. This will also create an `out/` folder that will contain the images.

## Running
Run the `./ceendra` executable. 

## Parameters: How to run the thing !
The general parameters are located in the `main.c` file. 

### Hyperparameters:
* ANTIALPOW: Must be a power of 2, determine the power of antialiasing. 4 is fine
* WIDTH, HEIGHT: Resolution of the image.
* BOUNDS: Bounds of the window into the complex plane. 1 is fine and will create a [-1;1] window.
* RANDBOUNDS: Bounds of the random complex number generator.
* EPSI: Espilon value used for the bailout function. The lower this value is, the more precise the image is, but will also take much longer to compute. 0.001 is fine.
* MAXWORD: The other bailout condition. Maximum lenght of a word before bailing out. Long word allows for more precise image, but increase compute time.
* LINE: Draw a line between endpoints. 1 for line drawing, 0 for just points.
* BITWISE: Experimental. Use a bitwise array for the image, allowing for a really good speedup.
* DEBUG: Will throw a lot of info on the current level, the current word, matrix, etc.


### ta and tb traces 
ta and tb are the initial complex traces used by the 'Grandma's formula' to create 2 generators and their inverse. The program will then compute the limit set of those generators. 

### Making animations  
One can also create sequences of images by varying the ta/tb parameters from a starting value to an end point. One way to do that is to use the easing functions I have implemented. 

### Farey and other Sequences
Another cool way is to use farey sequences to compute the mu parameter and get a group that is on the border between chaos and stability. The Farey Sequence is simply all the fractions between 0/1 and 1/1 with a specified maximum denominator. Using the makeFareySeq you can fill an array of rationals with the Farey Sequence, and then use a newton root finder to find the associated mu parameter and get a cool looking fraction out of it. You can also use the Fibonacci serie to get an approximation for irrational groups.

### Special Words
I have implemented a special word algorithm that finds the fixed point of some special sequences of words. Those points are the fixed point cyclic permutations of the special word. Given a special word in integer representation, you can get the associated fixed point array that will allow you to make better and faster plot.

### Fast
It's in C, so it's fast, despite my attempt at slowing the thing down. Moreover, I have (jankily) implemented multithreading. For now, the program uses 4 threads by default, but I should implement a way to find how many threads are available, and work with that.

## Documentation
As I have said the code is a bit messy, but I'll try to make some documentation explaining what's going on in details, and especially some guide in how to modify the program to make it do what you want it to do.

