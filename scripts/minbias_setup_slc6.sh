#!/bin/sh

source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6-gcc46-opt/setup.sh

export HEPMC=/afs/cern.ch/sw/lcg/external/HepMC/2.06.08/x86_64-slc6-gcc46-opt
export PYTHIA8DATA=/afs/cern.ch/cms/slc6_amd64_gcc472/external/pythia8/185/xmldoc
export LHAPDF=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.9/x86_64-slc6-gcc46-opt
export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/HepMC/2.06.08/x86_64-slc6-gcc46-opt/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/afs/cern.ch/cms/slc6_amd64_gcc472/external/pythia8/185/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.9/x86_64-slc6-gcc46-opt/lib:$LD_LIBRARY_PATH
