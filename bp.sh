# install conda and python dependencies
wget --quiet https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash ./Miniconda2-latest-Linux-x86_64.sh -b -p ./grackle-conda -f
export PATH=$PWD/grackle-conda/bin:$PATH
conda install -q -y mercurial cython h5py matplotlib sympy numpy pytest flake8 yt

# install OS dependencies
sudo apt-get update
sudo apt-get install -y csh libhdf5-serial-dev gfortran libtool

# download test dataset
wget --quiet http://yt-project.org/data/IsolatedGalaxy.tar.gz
tar xzf IsolatedGalaxy.tar.gz
export YT_DATA_DIR=$PWD

echo "backend : Agg" > $HOME/matplotlibrc
export MATPLOTLIBRC=$HOME

cd $BITBUCKET_CLONE_DIR
hg up tip

./configure
cd src/clib
make machine-linux-gnu
make
mkdir -p $HOME/local
make install
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
cd ../python
python setup.py develop
cd ../
py.test python/tests
