WORKDIR=/path/to/detectorSimulations_v10_build/

#TARGET="G4_Au" #if want to declare a target material

# declare -i number=0
# # 
# for k in $(seq 37 43);
# do
# number=$(( k*25 ))
# 
# cat FRAGS/frag1.txt > spiceauto.mac 
# #echo "/DetSys/app/SpiceTargetThickness" "0.9" "mg/cm2" >> spiceauto.mac
##if want to change the thickness across a number of runs 
# echo "/DetSys/gun/efficiencyEnergy" $number "keV"  >>  spiceauto.mac
## changing energies
# cat FRAGS/effFragLow.txt >> spiceauto.mac
# # 
# ./Griffinv10 spiceauto.mac #running the sim with macro created here
# #mv spiceauto.mac spiceauto$number.mac #if want to ave macro used etc
# mv g4out.root /path/to/new/"$number"_root.root #saving ROOT output
# mv kin.dat /path/to/new/"$number"_kin.dat #saving kinematic output
# done
