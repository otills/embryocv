conda create -n opencvtest2 -y numpy scipy scikit-learn matplotlib python=2.7

source activate opencvtest2

conda install -c menpo -y opencv3 
conda install -c groakat -y pathos=0.2a1
conda install -y pyqtgraph qtpy qt