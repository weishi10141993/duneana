# Run FD CAFMaker with PRISM GEC code from dunegpvm

## Setup using duneana develop branch

```
mkdir FDCAFIntegration
cd FDCAFIntegration

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

# The dunesw version below is what ups says in the develop branch when I cloned it Mar 7, 2023
# it changes as the develop branches update
# After integration work finish, need to push and resolve conflicts
setup dunesw v09_68_00d00 -q e20:debug  

mrb newDev
source /dune/app/users/weishi/FDCAFIntegration/localProducts_larsoft_v09_68_00_debug_e20/setup

# Get FD maker code
cd /dune/app/users/weishi/FDCAFIntegration/srcs

# On local laptop:
#   git clone https://github.com/weishi10141993/duneana.git -b develop
# the version was v09_68_00d00 when local copy is cloned, same version as in dunesw
# sync the changes in CAFMaker
cd ${LAPTOP_WORK_DIR}/duneana
rsync -e ssh -avSz  ./* weishi@dunegpvm13.fnal.gov:/dune/app/users/weishi/FDCAFIntegration/srcs/duneana

# Build GEC code inside CAFMaker
cd duneana/duneana/CAFMaker
./build.sh

# Build/rebuild FD CAF maker code
cd /dune/app/users/weishi/FDCAFIntegration/srcs
mrb uc
cd ${MRB_BUILDDIR}       
mrb z
mrbsetenv
mrb b
```

Run CAF Maker interactively:
```
cd /dune/app/users/weishi/FDCAFIntegration/srcs/duneana/duneana/CAFMaker
lar -c select_ana_dune10kt_nu.fcl -n 10 /pnfs/dune/tape_backed/dunepro/mcc11/protodune/mc/full-reconstructed/07/51/31/11/nu_dune10kt_1x2x6_13009312_0_20181104T221530_gen_g4_detsim_reco.root
```

Relogin:
```
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunesw v09_68_00d00 -q e20:debug  
source /dune/app/users/weishi/FDCAFIntegration/localProducts_larsoft_v09_68_00_debug_e20/setup
mrbsetenv
```

Recompile:
```
cd ${MRB_BUILDDIR}                       
mrb z                               
mrbsetenv                            
mrb b
```
