cmin=0.1
cmax=1.0
cinc=0.1

amin=0
amax=2.0
ainc=0.2

fmin=0
fmax=2.0
finc=0.2

config="config_files/config_t.json"
bests_file="bests/t_1.csv"
save_path="save_path/"
rm swarm_t_3d.swarm

echo "#SWARM --gb-per-process 16" >> "swarm_t_3d.swarm"
echo "#SWARM --gres lscratch:10" >> "swarm_t_3d.swarm"
echo "#SWARM --time 24:00:00" >> "swarm_t_3d.swarm"
echo "#SWARM --logdir output/" >> "swarm_t_3d.swarm"

for c in $(seq $cmin $cinc $cmax); do
    echo python runner_t_3d.py -c $config -b $bests_file --cmin $c --cinc $c --cmax $c --amin $amin --ainc $ainc --amax $amax --fmin $fmin --finc $finc --fmax $fmax -s $save_path  >> swarm_t_3d.swarm
done