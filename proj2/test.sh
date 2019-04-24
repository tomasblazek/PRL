#!/bin/bash

FILENAME="vuv"


if [ $# -eq 1 ];then 
	input=$1;
	size=${#input}
	if [ $size -eq 1 ];then
		let p="1"
	else	
		let p="2*$size-2"
	fi;
else
    echo "Error: Invalid count of arguments! (One argument requiered of type string)"
	exit 1
fi;


#preklad cpp
mpic++ --prefix /usr/local/share/OpenMPI -o ${FILENAME} ${FILENAME}.cpp

#spusteni
mpirun --prefix /usr/local/share/OpenMPI -np $p ${FILENAME} $input

#uklid
rm -f ${FILENAME}

