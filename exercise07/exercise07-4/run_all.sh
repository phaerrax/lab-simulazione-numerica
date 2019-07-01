#!/usr/bin/sh

phases=("solid" "liquid" "gaseous")
exe="./exercise07-4"

cleanup()
{
	echo "Cleaning up old frame files..."
	for phase in ${phases[@]}
	do
		# Cleanup...
		frames_dir=$phase/frames
		if [ -d $frames_dir ]
		then
			rm -f $frames_dir/*
		else
			mkdir $frames_dir
		fi
	done
}

run()
{
	for phase in ${phases[@]}
	do
		${exe}_$1 ${phase}
		echo "Simulation of $phase phase complete."
	done
}

# The programs write on and read from some config.*
# files, so it is not safe to launch all of them
# simultaneously.
if [ $# -eq 0 ]
then
	cleanup
	run verlet
	run metropolis
elif [ $# -eq 1 ]
then
	if [ $1 == "verlet" ]
	then
		cleanup
		run verlet
	elif [ $1 == "metropolis" ]
	then
		cleanup
		run metropolis
	else
		echo "Invalid input."
	fi
else
	echo "Invalid input."
fi

