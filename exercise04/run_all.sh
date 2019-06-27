#!/bin/sh

phases=("solid" "liquid" "gaseous")
exe="./exercise04"

for phase in ${phases[@]}
do
	# Cleanup...
	echo "Cleaning up old frame files..."
	frames_dir=$phase/frames
	if [ -d $frames_dir ]
	then
		rm -f $frames_dir/*
	else
		mkdir $frames_dir
	fi

	# ...and run.
	${exe} ${element} ${phase}
done
