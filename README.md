COCOS
=====

COncurrent system with COnstraints for protein Structure prediction (COCOS).
This is a Constraint Solver for Protein Structure Prediction (PSP) based on local search strategies.
We present an approach to the Protein Structure Prediction (PSP) problem that relies on a Multi-Agents System (MAS) perspective, where concurrent agents explore the folding of different parts of a protein. The strength of the approach lies in the agents’ ability to apply different types of knowledge, expressed in the form of declarative constraints, to prune the search space of folding alternatives. We demonstrate the suitability of a General-Purpose Graphical Processing Unit (GPGPU) approach to implement such MAS infrastructure, with significant performance improvements over the sequential implementation and other methods.

Software:
1) The CPU version of the MAS tool (C-COCOS).
2) The GPU version of the MAS tool (G-COCOS).

How to install C/G-COCOS:
In order to compile this tool, you need a standard Unix toolchain including a bash-compatible shell, and GNU make. On MacOS X, you need to install the XCode package. XCode is available from the Mac OS install DVD, or from the Apple Mac Dev Center. For Windows you can use the Cygwin environment that provides all necessary tools.

Don't forget that GMAS is based on CUDA. Therefore you need the CUDA framework installed in your computer!

Download and unzip MAS_CPU.zip (MAS_GPU.zip). The folder MAS_CPU (MAS_GPU) will contain all the files needed to run the tool. Then, do the following steps to set your compiler (if needed):

1) cd jnetsrc/src
2) edit CC and CFLAGS lines in the Makefile to reflect your compiler and optimiser
3) cd..
4) edit CC and CFLAGS lines in the Makefile to reflect your compiler and optimiser

Then you just need to run the script CompileAndRun.sh:
$./CompileAndRun.sh

Running examples:

After compiling CMAS(GMAS) , some examples can be run directly from the command line by invoking the executable file. The executable file is called cocos (COncurrent system with COnstraints for protein Structure prediction).

Use:

$ ./cocos -h

to print a help message.

In the folder proteins/ there are some examples of input files.

In particular, you may want to try:

$ ./cocos -i proteins/1ZDD.fa -a -v

that takes just a FASTA file 1ZDD.fa as input.

Otherwise, you can try:

$ ./cocos -i proteins/1ZDD.in.cocos -v

that uses 1ZDD.in.cocos as input file containing the descriptions of the agents to create.

Tests and benchmarks:

We tested the system on several proteins of different type and length. For each target we provide the input file and the best protein predicted by tool saved as a pdb file. These results can be downloaded from here (148.3 MB).

Thanks for reading this page. For any general question about the system feel free to contact write to campe8@nmsu.edu.

