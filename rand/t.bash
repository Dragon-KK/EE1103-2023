#!bin/bash
for i in $(awk '{print $1}' "data.txt");
do
    ping $i -i 0.02 -c 100 | grep "rtt min" | awk -F= '{print $2}' | awk -F/ '{print $2}'
done;