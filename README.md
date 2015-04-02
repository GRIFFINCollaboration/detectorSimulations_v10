detectorSimulations
===================

The detectorSimulations package contains the Geant4 simulations for GRIFFIN, TIGRESS, and all of their auxilary detectors.  Please note that in order to respect our non-disclosure agreements, all source containing third party IP have been omitted from this repo, and can be obtained from your colleagues directly at the lab.


#Setup

###Requirements
detectorSimulations is confirmed to run on geant4.9.6.p01 and ROOT 5.34.11.

###Getting the code

To setup the simulation package on a computer with GEANT4 already present, just copy the code to your machine:

    git clone https://github.com/GRIFFINCollaboration/detectorSimulations.git
    
Then you'll need to get the files containing our NDA-protected parameters from one of your colleagues at the lab; place these in the /src directory, and everything should compile and run as expected. 

###Building

Before building the simulation itself, some extra stuff needs to be compiled: in `detectorSimulations/dataRootClass`, do `make`.

The rest of the build process is pretty standard for a geant simulation; in a build directory (ie any clean new directory that isn't the source directory), do 

```
cmake -DGeant4_DIR=/path/to/geant4.9.6.p01-install/lib64/Geant4-9.6.1 ~/path/to/detectorSimulations
make clean
make
```

Keep in mind that cmake does not regenerate all the files it uses every time it runs!  So if something changes and this build process suddenly fails, try deleting the build directory and starting over.

###Setup FAQ

- Yes, you need both the secret suppressed files AND their unsuppressed equivalents, not just either / or.


#Usage

###Particle Emission

####General
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/gun/energy double unit ``` | Set energy of particle | 1000 keV |
| ``` /DetSys/gun/particle string ``` | Set particle type (e-, e+, gamma, proton, alpha) | gamma |
| ``` /DetSys/gun/ion Z A E ``` | Set ion type (excitation energy E in keV) | |
| ``` /DetSys/gun/direction x y z ``` | Set momentum direction | |
| ``` /DetSys/gun/position x y z  unit``` | Set particle position | 0.0 0.0 0.0 mm |
| ``` /DetSys/gun/radius r unit``` | Set source radius | 0.0 mm |
| ``` /DetSys/gun/energyrange min max step ``` | Set energy (keV) of particle, loops from min to max | |

####Decay Schemes
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/gun/betaPlusEmission filename ``` | Simulate beta plus decay with energy distribution input file | |
| ``` /DetSys/gun/betaMinusEmission filename ``` | Simulate beta negative decay with energy distribution input file | |
| ``` /DetSys/gun/polarization double``` | Set Polarization of Nuclei (before radioactiveBetaDecay) | |
| ``` /DetSys/gun/radioactiveBetaDecay directory ``` | Simulate complete beta negative decay with simulation directory | |
| ``` /DetSys/gun/emitBetaParticle 0/1 ``` | Emit Beta Particle? True/False | |
| ``` /DetSys/gun/includeXRayInputFileKShell 0/1 ``` | Emit X-rays from K-shell vacancies using input file? True/False | |
| ``` /DetSys/gun/includeXRayInputFileLShell 0/1 ``` | Emit X-rays from L-shell vacancies using input file? True/False | |
| ``` /DetSys/gun/includeXRayInputFileMShell 0/1 ``` | Emit X-rays from M-shell vacancies using input file? True/False | |
| ``` /DetSys/gun/radioactiveDecayHalflife double ``` | Half-life of radioactive isotope simulation (seconds) | |
| ``` /DetSys/gun/numberOfRadioactiveNuclei int ``` | Set the number of radioactive nuclei | |
| ``` /DetSys/gun/radioactiveSourceDecay filename ``` | Simulate source decay with a file containing the decay data | |

####Kinematics
Whenever an electron or positron is emitted, you can choose to simulate kinematic effects (i.e. shifting of energy due to relativistic speeds of ions).  You must specify the ion with ``` /DetSys/gun/ion Z A E ``` first.

| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/gun/simulateKinematics 0/1 ``` | Choose to simulate kinematic effects | False |
| ``` /DetSys/gun/kinematicsBetaValue double ```* | Specify beta value of incident ion | 0 |
| ``` /DetSys/gun/kinematicsIonEnergy double unit ```* | Specify energy of heavy ion | 0 MeV |

*Only one of these should be specified

###Detector Specific

####GRIFFIN
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/det/addGriffinForward int ``` | Add Detection System GriffinForward |  |
| ``` /DetSys/det/addGriffinForwardDetector int ``` | Add GriffinForward Detector |  |
| ``` /DetSys/det/addGriffinBack int ``` | Add Detection System GriffinBack |  |
| ``` /DetSys/det/addGriffinBackDetector int ``` | Add GriffinBack Detector |  |
| ``` /DetSys/det/UseTIGRESSPositions ``` | Use TIGRESS detector positions rather than GRIFFIN | False |
| ------ | ---------------- | ------ |
| ``` /DetSys/det/addGriffinCustomDetector 0 ``` | Adds a detector using the paramaters specified |  |
| ``` /DetSys/det/SetCustomShieldsPresent 0/1 ``` | Selects whether or not the detector suppressors are included | True |
| ``` /DetSys/det/SetCustomRadialDistance double unit ``` | Selects the radial distance for the detector from the origin |  |
| ``` /DetSys/det/SetCustomExtensionSuppressorLocation 0/1 ``` | Selects a position for the extension suppressors. Either forward (0) or back (1) | Forward (0) |
| ``` /DetSys/det/SetCustomDeadLayer det_num pos_num null ``` | Sets the dead layer and position for the detector placed in the next call to addGriffinCustom |  |
| ``` /DetSys/det/includeGriffinHevimet 0/1 ``` | Includes the Hevimet for a Griffin detector | False |
| ``` /DetSys/det/addGriffinCustom int ``` | Adds a detection system using the paramaters specified |  |

####SPICE
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/app/addSpiceTargetChamber ``` | Add SPICE target chamber |  |
| ``` /DetSys/Spice/setResolution double double ``` |Set resolution of SPICE Si(Li)  |  |
| ``` /DetSys/det/addSpice int``` | Add Si(Li) detector |  |
| ``` /DetSys/det/addS3 ``` | Add SPICE S3 detector |  |

####8PI
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/det/add8pi int ``` | Add Detection System 8pi |  |
| ``` /DetSys/det/add8piDetector int ``` | Add 8pi Detector |  |
| ``` /DetSys/app/add8piVacuumChamber ``` | Add 8pi vacuum chamber |  |
| ``` /DetSys/app/add8piVacuumChamberAuxMatShell double unit ``` | Add AuxMat shell around 8pi vacuum chamber with specified thickness |  |

####Other
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/det/addGammaTracking int ``` | Add Detection System GammaTracking |  |
| ``` /DetSys/det/addSodiumIodide int ``` | Add Detection System SodiumIodide |  |
| ``` /DetSys/det/addLanthanumBromide int ``` | Add Detection System LanthanumBromide |  |
| ``` /DetSys/det/addSceptar int ``` | Add Detection System Sceptar |  |
| ``` /DetSys/det/addPaces int ``` | Add Detection System Paces |  |

###Detector General

``` /DetSys/det/update ``` must be called before using ``` beamOn ``` if the geometrical values
have been altered.

####World Volume (i.e. experimental hall)
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/world/material string ``` | Select material for the world | G4_AIR |
| ``` /DetSys/world/dimensions x y z unit ``` | Set world dimensions | 10m x 10m x 10m |
| ``` /DetSys/world/vis 0/1``` | Set world visibility (depreciated) | False |

####Generic Target
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/app/genericTarget string ``` | Create a target with specified material |  |
| ``` /DetSys/app/genericTargetDimensions x y z unit ``` | Set target dimensions |  |
| ``` /DetSys/app/genericTargetPosition x y z unit ``` | Set target position |  |
Requires all three commands to build

####Field Box
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/app/fieldBoxMaterial string ``` | Create a field box with specified material |  |
| ``` /DetSys/app/fieldBoxDimensions x y z unit ``` | Set field box dimensions |  |
| ``` /DetSys/app/fieldBoxPosition x y z unit ``` | Set field box position |  |
| ``` /DetSys/app/fieldBoxMagneticField x y z unit ``` | Set field box magnetic field |  |
Requires all four commands to build

####Box
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/det/boxMat string ``` | Set box material | G4_WATER |
| ``` /DetSys/det/boxThickness double unit ``` | Set box thickness | 0.0 mm |
| ``` /DetSys/det/boxInnerDimensions x y z unit ``` | Set box inner dimensions | 0.0 0.0 0.0 mm |
| ``` /DetSys/det/boxColour r g b``` | Set box colour | 0.0 0.0 1.0 |
| ``` /DetSys/det/addBox ``` | Build/add box (if thickness is not 0) |  |

####Grid
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/det/gridMat string ``` | Set grid material | G4_WATER |
| ``` /DetSys/det/gridSize double unit ``` | Set grid size | 0.0 mm |
| ``` /DetSys/det/gridDimensions x y z unit ``` | Set grid dimensions | 0.0 0.0 0.0 mm |
| ``` /DetSys/det/gridColour r g b ``` | Set grid colour | 1.0 0.0 0.0 |
| ``` /DetSys/det/addGrid ``` | Build/add grid (if grid size is not 0) |  |

####Magnetic Fields
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/world/magneticField x y z unit``` | Set world magnetic field (depreciated) | 0, 0, 0 |
| ``` /DetSys/world/tabMagneticField filename ``` | Set tabulated magnetic field* | Disabled | 

*Used for SPICE: these files are 50Mb each and are different for each lens so until a better solution comes along they will be kept on the network at TRIUMF (email Mohamad or Lee for precise location)

| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ```  ``` |  |  |
| ```  ``` |  |  |
| ```  ``` |  |  |
| ```  ``` |  |  |
