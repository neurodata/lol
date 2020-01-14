#!/bin/bash -e
# run this in the background with nohup ./memlog.sh > mem.txt &
#

echo "time$(free -m | grep total | sed -E 's/^    (.*)/\1/g')"
while true; do
    echo "$(free -m | grep Mem)"
    sleep 1
done
