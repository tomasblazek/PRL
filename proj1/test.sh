#!/bin/bash

FILENAME="bks"

#pocet cisel bud zadam nebo 16
if [ $# -lt 1 ];then 
    numbers=16;
else
	#kontrola validity vstupu
	if [[ $1 =~ ^[1-9][0-9]*$ ]]; then
		numbers=$1;
	else
		echo "Error: Invalid argument! (Unsigned integer value higher than zero required)"
		exit 1
	fi;
fi;


#preprocesing pro algoritmus - vypocet procesoru
if [ "$numbers" -gt 0 ];then
	logNumber=`echo "l($numbers)/l(2)" | bc -l`
	logNumber=`python3 -c "import math; print(math.ceil($logNumber))"`

	i=0
	pow2=1
	while [[ "$logNumber" -gt "$pow2" ]]; do
		let i=$i+1 
		pow2=`echo "2^$i" | bc`
	done	
	m=`echo "2^$i" | bc`
	let p="$m*2-1"
fi;


#preklad cpp
mpic++ --prefix /usr/local/share/OpenMPI -o ${FILENAME} ${FILENAME}.cpp


#vyrobeni souboru s random cisly
dd if=/dev/urandom bs=1 count=$numbers of=numbers 2>/dev/null

#spusteni
mpirun --prefix /usr/local/share/OpenMPI -np $p ${FILENAME}

#uklid
rm -f ${FILENAME} numbers

