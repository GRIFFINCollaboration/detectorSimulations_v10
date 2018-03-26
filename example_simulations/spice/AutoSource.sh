energypoints=(481.6935 553.8372 565.8473 975.651 1047.795 1059.805 1682.224 1754.367)
relint=(21709 6243 1568 100000 25989 6215 336 48)
Multiplyer=1
## Change Multiplier to any integer to increase statistics based on time/precision 1-1000

rm -rf SourceTmp_*.root
# 
for ((i=0;i<${#energypoints[@]};++i));
do
	### Look at helpfile and spice.mac for full details

	### Dont change ###
	echo "/run/initialize" > spiceauto.mac 
	echo "/process/em/fluo false" >> spiceauto.mac 
	echo "/process/em/auger false" >> spiceauto.mac 
	echo "/process/em/pixe false" >> spiceauto.mac 
	echo "/DetSys/world/material Vacuum" >> spiceauto.mac 
	echo "/DetSys/gun/particle e-" >> spiceauto.mac 
	echo "/DetSys/gun/efficiencyEnergy ${energypoints[i]} keV" >> spiceauto.mac 
	echo "/DetSys/gun/BeamSpot 0.5 mm" >> spiceauto.mac 
	### Dont change ###

	echo "/DetSys/gun/position 0.0 0.0 0.0 mm " >> spiceauto.mac 

	echo "/DetSys/app/addSpiceTargetChamber Med" >> spiceauto.mac 
	echo "/DetSys/det/addSpiceDetector" >> spiceauto.mac 
	echo "/DetSys/world/TabMagneticField MagneticField3D.MED.TABLE -1 -45" >> spiceauto.mac 

	### Standard disk source, thickness in mg/cm2
	echo "/DetSys/app/LayeredTargetAddLayer G4_Al 200" >> spiceauto.mac 
	echo "/DetSys/app/LayeredTargetAddLayer Bismuth 0.05" >> spiceauto.mac 
	echo "/DetSys/app/LayeredTargetAddLayer Acrylic 0.2" >> spiceauto.mac 
	echo "/Detsys/gun/TargetLayer 1" >> spiceauto.mac 

	echo "/run/beamOn ${relint[i]*Multiplyer}" >> spiceauto.mac 

	#running the sim with macro created here
	./Griffinv10 spiceauto.mac 

	mv g4out.root SourceTmp_"$i".root

done

rm -rf SPICESource.root
hadd SPICESource.root SourceTmp_*.root

