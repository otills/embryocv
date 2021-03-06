# ReadMe: Instructions on setting up EmbryoCV on your machine.

EmbryoCV has been developed using the Anaconda package manager and Spyder - an editor that can run with Python.

Setting up EmbryoCV requires several steps:

> Download and install Anaconda with Python 2.7 https://www.anaconda.com/download/#macos
> Download the EmbryoCV source code from the Github repository www.embryocv.org. Unpack in a suitable location.
> Run the setup.py script in Terminal (Mac) or Command prompt (PC). Run like so (source is required to active the virtual environment):

source setup.py

This will install OpenCV and its dependencies within a Python virtual environment called opencv. Use this virtual environment when working with EmbryoCV.

> Once installed, to run EmbryoCV you run the following in your command line editor (Terminal or Command prompt):
# Launch the virtual environment
> source activate opencv  
# Open Spyder a nice editor for working with EmmbryoCV.
> spyder
      
This should open Spyder. On the first time of running you need to add the folder into which you downloaded the EmbryoCV source code. Click Python (top left hand corner), PYTHONPATH manager and then 'Add path'. Find the folder containing the source code and add this. You can check this has worked by typing import EmbryoCV and if you do not receive an error message saying it wasnÕt found then you have successfully installed EmbryoCV and its dependencies (you may receive some depreciation warnings). If you encounter trouble give some description on the EmbryoPhenomics Google Group and we will help.
      

		

