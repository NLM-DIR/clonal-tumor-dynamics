m_min=-1
m_max=0
m_inc=0.005

k_min=-1
k_max=1
k_inc=0.005

config="config_files/config_km.json"
save_path="save_path/"
rm swarm_km.swarm

echo "#SWARM --gb-per-process 16" >> "swarm_km.swarm"
echo "#SWARM --gres lscratch:10" >> "swarm_km.swarm"
echo "#SWARM --time 24:00:00" >> "swarm_km.swarm"
echo "#SWARM --logdir output/" >> "swarm_km.swarm"

for m in $(seq $m_min $m_inc $m_max); do
        echo python runner_km.py -c $config --mmin $m --minc $m_inc --mmax $m --kmin $k_min --kinc $k_inc --kmax $k_max -s $save_path  >> swarm_km.swarm
done