LatinoTreesDelphes
==================

Latino tree dumper from Delphes

* Minimum bias samples generation with pythia8:
```   
    - source scripts/minbias_setup_slc6
    - gcc -I$HEPMC/include/ -L$HEPMC/lib -I$PYTHIA8DATA/../include -L$PYTHIA8DATA/../lib -I$LHAPDF/include -L$LHAPDF/lib
          -lLHAPDF -lHepMC -lpythia8tohepmc -lpythia8 -o genMinBias_14TeV genMinBias_14TeV.cpp
    - scripts/Pythia8Jobs.py -j 10 -n 10000 -q 1nh -f minBias_test -d /afs/cern.ch/users/s/someone/...
```
* Pythia8Jobs.py: starts #j jobs of #n events on queue #q, produces a file named #f.hepmc and copies it to #d 
* genMinBias_14TeV needs 2 arguments:
```
    - ./genMinBias_14TeV number_of_events output_file
``` 
