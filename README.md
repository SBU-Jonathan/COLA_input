# COLA_input
Input files used to configure and run our COLA simulations. These include transfer functions, lua files, transferinfo files, files with cosmological parameters, and files containing information about the background. The latter 2 sets of files are not necessary but we store them here. 

There is a directory for each model: 'wCDM', 'LCDM'. Within each model there are directories for each set of simulations, or 'project'. Each project is named by the number of cosmologies, then an underscore, then a 1 or 2 to indicate COLA-precision (default or high respectively). 

For example LCDM/800_1 indicates 800 paired-and-fixed LCDM simulations in default precision.
