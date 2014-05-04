simutil
=======

Scripts and small utilitary programs to assist in simulations (using DL_POLY, GROMACS or Towhee). The utilitaries are based on a specific program's output or input files.
Every tool is designed to work with as little as possible of input. Idealy the programs run without the need to provide any input, the tool runs taking what it need from the folder in which it is called.

## Folder structure
The folder names are self explanatory. None the less it might be useful to lay out the reasoning behind it.
* `bin´ - In the bin folder you will find compiled programs.
* `source´ - In the source folder live the .f90 files with the respective code for each tool. 
* `test´ - Finaly, in the test folder exists a folder with the necessary files to run the compiled code as a test, as well as the supposed output for that test.