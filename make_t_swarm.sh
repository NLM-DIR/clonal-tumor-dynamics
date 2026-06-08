# Don't go over
cmin=0.02
cmax=1.0
cinc=0.02

# Don't go over
amin=0.1
amax=2.0
ainc=0.1

# Increase max by 1
fmin=0.1
fmax=2.1
finc=0.1

# Increase max by 1
bmin=0.01
bmax=0.11
binc=0.01

config="config_files/config_t.json"
save_path="save_path/"
rm swarm_t.swarm

echo "#SWARM --gb-per-process 16" >> "swarm_t.swarm"
echo "#SWARM --gres lscratch:10" >> "swarm_t.swarm"
echo "#SWARM --time 24:00:00" >> "swarm_t.swarm"
echo "#SWARM --logdir output/" >> "swarm_t.swarm"

for c in $(seq $cmin $cinc $cmax); do
    for a in $(seq $amin $ainc $amax); do
        echo python runner_t.py -c $config --cmin $c --cinc $c --cmax $c --amin $a --ainc $a --amax $a --fmin $fmin --finc $finc --fmax $fmax --bmin $bmin --binc $binc --bmax $bmax --lmin 0.01 --linc 0.01 --lmax 0.01 -s $save_path  >> swarm_t.swarm
    done
done