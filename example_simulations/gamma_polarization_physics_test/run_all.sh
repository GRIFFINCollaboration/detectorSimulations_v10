/home/local1/Builds/v10_GRIFFIN/Griffinv10 run.mac
mv g4*wrl /mnt/hgfs/swap/

rm -fr histo

mkdir histo
/home/local1/Programs/AsciiFromRoot/AsciiFromRoot -rf g4out.root -hn histo/griffin_crystal_unsup_edep_det0
/home/local1/Programs/AsciiFromRoot/AsciiFromRoot -rf g4out.root -hn histo/griffin_crystal_unsup_edep_det4
/home/local1/Programs/AsciiFromRoot/AsciiFromRoot -rf g4out.root -hn histo/griffin_crystal_unsup_edep_det5
/home/local1/Programs/AsciiFromRoot/AsciiFromRoot -rf g4out.root -hn histo/griffin_crystal_unsup_edep_det6
/home/local1/Programs/AsciiFromRoot/AsciiFromRoot -rf g4out.root -hn histo/griffin_crystal_unsup_edep_det12


