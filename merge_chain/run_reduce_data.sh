#!/bin/bash

origin_dir=$PWD
chain_id=$1
stop_id=$2

#source /ccs/home/mojia/.bashrc
cd /gpfs/alpine/phy171/proj-shared/mojia/MaCh3/MaCh3_CorrFDSys
#source setup.sh

#chain_dir=/gpfs/alpine/proj-shared/phy171/Asimov_fit_chains
chain_dir=/gpfs/alpine/proj-shared/phy171/Data_fit_chains

while [ $chain_id -lt $stop_id ]; do
  #./AtmJointFit_Bin/reduceDataSet ${chain_dir}/MaCh3_MCMC_chain_${chain_id}.root ${chain_dir}/MaCh3_MCMC_chain_${chain_id}_reduced.root &> ${chain_dir}/reduce_chain_${chain_id}.log
  ./AtmJointFit_Bin/reduceDataSet ${chain_dir}/MaCh3_MCMC_chain_${chain_id}_Iter_0_7.root ${chain_dir}/MaCh3_MCMC_chain_${chain_id}_Iter_0_7_reduced.root &> ${chain_dir}/reduce_chain_${chain_id}.log
 
  let chain_id+=1
done

cd $origin_dir

