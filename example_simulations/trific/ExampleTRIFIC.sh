cat part_trif.mac > tmp.mac
echo "/DetSys/det/addTrificDetector 98.0"  >> tmp.mac
echo "/gun/ion 19 41 0 0.0"  >> tmp.mac
echo "/DetSys/gun/efficiencyEnergy 250 MeV"  >> tmp.mac
echo "/run/beamOn 20000"  >> tmp.mac
./Griffinv10 tmp.mac
root -l -q 'LatestTrificSort.C(false,"TRIFIC_A.root")'

cat part_trif.mac > tmp.mac
echo "/DetSys/det/addTrificDetector 98.0"  >> tmp.mac
echo "/gun/ion 20 41 0 0.0"  >> tmp.mac
echo "/DetSys/gun/efficiencyEnergy 250 MeV"  >> tmp.mac
echo "/run/beamOn 20000"  >> tmp.mac
./Griffinv10 tmp.mac
root -l -q 'LatestTrificSort.C(false,"TRIFIC_B.root")'


cat part_trif.mac > tmp.mac
echo "/DetSys/det/TrificFlatWindow true"  >> tmp.mac
echo "/DetSys/det/addTrificDetector 90.0"  >> tmp.mac
echo "/gun/ion 19 41 0 0.0"  >> tmp.mac
echo "/DetSys/gun/efficiencyEnergy 250 MeV"  >> tmp.mac
echo "/run/beamOn 20000"  >> tmp.mac
./Griffinv10 tmp.mac
root -l -q 'LatestTrificSort.C(true,"TRIFIC_C.root")'

cat part_trif.mac > tmp.mac
echo "/DetSys/det/TrificFlatWindow true"  >> tmp.mac
echo "/DetSys/det/addTrificDetector 90.0"  >> tmp.mac
echo "/gun/ion 20 41 0 0.0"  >> tmp.mac
echo "/DetSys/gun/efficiencyEnergy 250 MeV"  >> tmp.mac
echo "/run/beamOn 20000"  >> tmp.mac
./Griffinv10 tmp.mac
root -l -q 'LatestTrificSort.C(true,"TRIFIC_D.root")'

hadd TRIFIC_AB.root TRIFIC_A.root TRIFIC_B.root
hadd TRIFIC_CD.root TRIFIC_C.root TRIFIC_D.root

# /DetSys/det/TrificAluWindow true
# /DetSys/det/TrificWindowThickness 6 um
