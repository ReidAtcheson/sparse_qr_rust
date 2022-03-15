
this_script=$(realpath $0)
path=$(dirname $this_script)
cd $path
git clone https://github.com/KarypisLab/GKlib.git
cd $path/GKlib
make config cc=gcc prefix=$path
make install
cd ..
git clone https://github.com/KarypisLab/METIS.git
cd $path/METIS
make config cc=gcc i64=1 prefix=$path glib_path=$path
make install
