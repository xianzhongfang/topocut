#! /bin/bash

cd `dirname $0`
cur_dir=`pwd`

## input directory
in_dir=${cur_dir}/example

## output directory
tmpdir=./tmpdir
build_dir=${tmpdir}/build
## the executable file
exe=${build_dir}/bin/topocut

if [ ! -f ${exe} ] ; then
    #sudo apt install libeigen3-dev libboost-program-options-dev libgmp-dev libcgal-dev
    mkdir -p ${build_dir}
    cd ${build_dir}
    cmake ${cur_dir} -DCMAKE_BUILD_TYPE=Release
    make -j5
    cd -
fi

### set mesh name
#mesh=eight
mesh=cylinder

## set output dir
out_dir=${tmpdir}/output
out_mesh_dir=${out_dir}/${mesh}
mkdir -p ${out_mesh_dir}

## output files
obj=${in_dir}/${mesh}.obj
plane=${out_mesh_dir}/${mesh}.plane
cut_mesh=${out_mesh_dir}/${mesh}.cm.vtk



#### Generate cut planes
## axis-aligned cuts
${exe} plane -i ${obj} -o ${plane} -s bound_box_xyz --nx 5 --ny 5 --nz 5
## or use random cuts
#${exe} plane -i ${obj} -o ${plane} -s random -P 10

## Generate cut-mesh by using cut-planes
${exe} cut -i ${obj} -p ${plane} -o ${cut_mesh} --is-triangulate 0

