#!/bin/bash

for dir in velset_mini velset_2units_new three_unit_new two_pair_new; do 

cd /data/home/shubhangi/sm_frmwrk
cd $dir
pwd

for data in human hprd flickr; do
edgef=$data
mkdir -p output/$data
echo -ne "$dir:\t$data\n"
./subgraph ~/benchmark_dataset/$data/${data}_labels ~/benchmark_dataset/$data/${edgef}_edges <( grep -v v3 ../qfiles/qfile_${data}) > output/$data/output_all &
~/plugins/memory_monitor.sh $! > output/$data/${data}_mem &
# echo "done"
done

data=small_ppi
mkdir -p output/$data
echo -ne "$dir:\t$data - "
edgef=ppi_dataset/experiments/small_ppi/ip_v12k_edges_1
./subgraph ppi_dataset/experiments/small_ppi/ip_v12k_label_1 $edgef ../qfiles/qfile_${data} > output/$data/output_all &
~/plugins/memory_monitor.sh $! > output/$data/${data}_mem &

done
