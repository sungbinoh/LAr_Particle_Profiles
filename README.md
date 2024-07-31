# LAr Particle Profiles

## Purpose
 To estimate particle profiles inside LAr

## Copy and setup
```
git clone git@github.com:sungbinoh/LAr_Particle_Profiles.git
source setup.sh
```

## Compile
```
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build . -- install
```

## Run
```
root -l -b -q run_profile.C
```