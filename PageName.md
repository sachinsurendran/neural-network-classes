#Experiments with Genetic Algorithm

# Introduction #

Goal:
> Build a self evolving learning engine, which should be flexible to learn any problem presented to it. Only input the engine should expect is a number which grades its performance.

Genetics:
> Behind any great achievements of mankind, be it the linux kernel or ability to launch robots in Mars, Credit goes to the designer of Man, ie evolution.
> If evolution could design something as complex as man, who achieves great feats. There is a lot we could learn from evolution. Hence most of the attempts to build this engine borrows from nature and copies evolution at every possible juncture.

Neural networks:
> Neural networks (NN) is a simple representation of neurons and can solve quite a few problems. Key aspects about neural networks which makes it a good fit for our problem is:

- Very distributed, cut a chunk of a NN and it still works in many ways.
> A C program if cut in few random places, just fails to work.
> This feature makes NN suited for operation like mutation, crossover...

- Very Parallel
> Extreme parallelism, makes us improvise the scale of the system and run it in a grid of linux machines

Code for NN:
> there is a well written library for NN done by a statistical institute which I have used for the engine,
The code base is heavily modified by us for fitting into our system, the library can be found here:

Flood: An Open Source Neural Networks C++ Library

Problem:

The current problem I have chosen is Tennix. Tennix is not a very simple problem neither a very complicated one.Hence it suited well for initial experiments with Genetics and to understand and modify it better.
> We will use it as an example to refine our system before we move on to a more interesting problem.

Code Infrastructure:

Tennix: http://code.google.com/p/tennix/source/browse/#svn%2Ftrunk

Neural Network & Genetics Engine: http://code.google.com/p/neural-network-classes/source/browse/#svn%2Ftrunk