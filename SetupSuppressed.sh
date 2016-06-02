#!/bin/bash

# check if the user provided a directory where the suppressed files are installed, otherwise use standard 'suppressed'
if [ -z "$1" ]; then
	supp_dir="suppressed"
else
	supp_dir=$1
fi

# if the directory with the suppressed files doesn't exist, clone it
if [ ! -d $supp_dir ]; then
	if git clone ssh://git@gitlab.com/GRIFFINCollaboration/suppressed.git $supp_dir; then
		echo "Successfully cloned the suppressed files!"
	else
		echo "Failed to clone the suppressed files. Please check that your gitlab account has the rights to access this project and that your ssh-keys are properly set up!"
		exit;
	fi
fi

# to set the link properly, we need to convert the potentially relative path to an absolute path
supp_dir=$(cd $supp_dir; pwd)

# check that the two suppressed files exist (in case someone created the directory for the suppressed files but didn't clone)
if [[ ! -f $supp_dir/DetectionSystemGriffinSuppressed.cc || ! -f $supp_dir/DetectorConstructionSuppressed.cc ]]; then
	if git clone ssh://git@gitlab.com/GRIFFINCollaboration/suppressed.git $supp_dir; then
		echo "Successfully cloned the suppressed files!"
	else
		echo "Failed to clone the suppressed files. Please check that your gitlab account has the rights to access this project, that your ssh-keys are properly set up, and that $supp_dir is empty!"
		exit;
	fi 
fi

# at this point the directory with the suppressed files should have been set up correctly, so we just need to set the links
ln -sf $supp_dir/DetectionSystemGriffinSuppressed.cc src/DetectionSystemGriffinSuppressed.cc
ln -sf $supp_dir/DetectorConstructionSuppressed.cc src/DetectorConstructionSuppressed.cc

echo "The suppressed files are set up, you can now compile the simulation!"

exit;
