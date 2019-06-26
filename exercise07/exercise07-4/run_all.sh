#!/usr/bin/sh

phases=("solid" "liquid" "gaseous")
exe="./exercise07-4"

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

for phase in ${phases[@]}
do
	# The programs write on and read from some config.*
	# files, so it is not safe to launch all of them
	# simultaneously.
	${exe}_verlet ${phase}
	${exe}_metropolis ${phase}
	echo "Simulation of $phase phase complete."
done
