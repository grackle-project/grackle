# install conda and python dependencies
wget --quiet https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash ./Miniconda2-latest-Linux-x86_64.sh -b -p ./grackle-conda -f
export PATH=$PWD/grackle-conda/bin:$PATH
conda install -q -y mercurial cython h5py matplotlib sympy numpy pytest flake8 yt

# install OS dependencies
sudo apt-get update
sudo apt-get install csh libhdf5-serial-dev gfortran libtool

cd $BITBUCKET_CLONE_DIR
hg up tip

echo $PATH

./configure
cd src/clib
make machine-linux-gnu
make
mkdir -p $HOME/local
make install
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
cd ../python
python setup.py develop
make test
