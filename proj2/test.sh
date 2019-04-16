#!/bin/bash

FILENAME="vuv"

#pocet procesorů bud zadam nebo 16
if [ $# -lt 1 ];then 
    p=16;
else
	#kontrola validity vstupu
	if [[ $1 =~ ^[1-9][0-9]*$ ]]; then
		p=$1;
	else
		echo "Error: Invalid argument! (Unsigned integer value higher than zero required)"
		exit 1
	fi;
fi;

#preklad cpp
mpic++ --prefix /usr/local/share/OpenMPI -o ${FILENAME} ${FILENAME}.cpp

#spusteni
mpirun --prefix /usr/local/share/OpenMPI -np $p ${FILENAME}

#uklid
rm -f ${FILENAME} numbers

