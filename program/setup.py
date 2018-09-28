conda create -n -y opencvtest numpy scipy scikit-learn matplotlib python=2.7
source activate opencvtest
conda install -c -y menpo opencv3 
conda install -c groakat -y pathos=0.2a1
conda install -y pyqtgraph qtpy qt