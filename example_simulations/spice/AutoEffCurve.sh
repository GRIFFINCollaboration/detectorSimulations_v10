g++ SpiceMiniSort.C `root-config --cflags --libs` -o sMini

### Because of the complex shape of the efficiency curves
### An uneven distribution of points is best, add more points in complex areas
### Low Energy Lens
#energypoints="50 100 125 150 175 200 225 250 275 300 325 350 375 400 425 450 500 550 650 700 800 900 1000 1200 1400 1600"
### Medium Energy Lens
energypoints="100 200 250 300 350 375 400 425 450 475 500 525 550 575 600 650 700 750 800 850 900 1000 1100 1200 1400 1600 2000"

rm -rf SPICE_*.root

### Change based on time/precision
### Set as low as 10,000 for a very rough fast answer and a maximum of 10,000,000 
NumberOfPoints=10000

echo $NumberOfPoints > EffPoints.txt

for E in $energypoints;
do

	### Look at helpfile and spice.mac for full details

	### Dont change ###
	echo "/run/initialize" > spiceauto.mac 
	echo "/process/em/fluo false" >> spiceauto.mac 
	echo "/process/em/auger false" >> spiceauto.mac 
	echo "/process/em/pixe false" >> spiceauto.mac 
	echo "/DetSys/world/material Vacuum" >> spiceauto.mac 
	echo "/DetSys/gun/particle e-" >> spiceauto.mac 
	echo "/DetSys/gun/efficiencyEnergy $E keV" >> spiceauto.mac 
	echo "/DetSys/gun/BeamSpot 0.5 mm" >> spiceauto.mac 
	### Dont change ###

	### Select the desired options and the corresponding field file field
	echo "/DetSys/app/addSpiceTargetChamber MedBunnyCol" >> spiceauto.mac 
	echo "/DetSys/world/TabMagneticField MagneticField3D.MED.TABLE -1 -45" >> spiceauto.mac 
	#echo "/DetSys/world/TabMagneticField MagneticField3D.LOW.TABLE -1 -45" >> spiceauto.mac 
	echo "/DetSys/det/addSpiceDetector" >> spiceauto.mac 

	### Set target details, thickness in mg/cm2
	### Position should match options (bunny nominally Z=-8.)
	echo "/DetSys/gun/position 0.0 0.0 -8. mm " >> spiceauto.mac 
	echo "/DetSys/app/LayeredTargetAddLayer Gold 1.8" >> spiceauto.mac 
	echo "/DetSys/gun/TargetLayer 0" >> spiceauto.mac 

	echo "/run/beamOn $NumberOfPoints" >> spiceauto.mac 

	#running the sim with macro created here
	./Griffinv10 spiceauto.mac 

	counts=$(./sMini $E 0.05)
    echo "$E $counts" >> EffPoints.txt
	mv g4out.root SPICE_"$E".root
	
done

root -b -l -q MakeEfficiency.c
