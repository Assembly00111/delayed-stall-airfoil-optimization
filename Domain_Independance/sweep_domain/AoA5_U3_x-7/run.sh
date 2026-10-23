#!/bin/bash

#mesh generation
foamCleanTutorials
foamCleanCase

blockMesh | tee log.blockMesh

snappyHexMesh -overwrite | tee log.snappyHexMesh

extrudeMesh | tee log.extrudeMesh

createPatch -overwrite | tee log.createPatch

#mesh checking
checkMesh -latestTime | tee log.checkMesh

#running case file
pimpleFoam | tee log.pimpleFoam

cd ..
