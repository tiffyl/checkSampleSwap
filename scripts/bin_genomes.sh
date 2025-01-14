#! /bin/bash

window=$1
chr=$2
start=$3
regionend=$4

end=$(($start + $window - 1))

while [ $start -lt $regionend ]
do
    echo -e "$chr:$start-$end"
    start=$(( $start + $window ))
    end=$(( $end + $window ))
done
