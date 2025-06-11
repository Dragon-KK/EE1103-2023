#!/bin/bash 
# EE23B135 Kaushik G Iyer
# 14/11/2023
# Basically creates a file `ee23b135_hem.txt` which contains the values (i, hem(i), hem(i)/hem(i+1))
# Then uses the previously created file to plot it as required (plot is defined in ee23b135_hemgp.txt)
# NOTE: We expect our ee23b135_hem.c file to be compiled into ee23b135_hem.out!

touch ee23b135_hem.txt
rm ee23b135_hem.txt
for i in $(seq 1 1 40);
do
    iplus1=$(echo "$i + 1" | bc)
    ival=$(./ee23b135_hem.out $i)
    iplus1val=$(./ee23b135_hem.out $iplus1)
    echo $i $ival $(echo "scale=10; $ival / $iplus1val" | bc) >> ee23b135_hem.txt
done;

gnuplot -p "ee23b135_hemgp.txt"