#!/bin/sh

start_type=("close" "far")
exe="./exercise05"

for type in ${start_type[@]}
do
	${exe}_raw ${type}
	${exe}_raw_equilibrated ${type}
	${exe} ${type}
done
