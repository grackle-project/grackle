# install conda and python dependencies
wget --quiet https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash ./Miniconda2-latest-Linux-x86_64.sh -b -p ./grackle-conda -f
export PATH=$PWD/grackle-conda/bin:$PATH
conda install -q -y mercurial cython h5py matplotlib sympy numpy pytest flake8 yt

# install OS dependencies
sudo apt-get update
sudo apt-get install csh libhdf5-serial-dev gfortran

cd $BITBUCKET_CLONE_DIR
hg up tip

csh configure
cd src/clib
make machine-linux-gnu
make
cd ..
make build-python
cd python
make test
