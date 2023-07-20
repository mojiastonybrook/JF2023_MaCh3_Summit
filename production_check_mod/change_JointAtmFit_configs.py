import os
import re
from glob import glob

def change_configs(chain_list,itr_start):
    for chain in chain_list:
        print(chain, " : ")
        for i in range(itr_start,11):
            config = chain+"/config/AtmConfig_Iter_"+str(i)+".cfg"
            sed_command="sed -i 's|NSTEPS.*|NSTEPS = 1000000 |' "+config
            os.system(sed_command)
            sed_command="sed -i 's|STEPSCCALECORRFD.*|STEPSCCALECORRFD = 0.008 |' "+config
            os.system(sed_command)
            print("changed config for iteration: ", i)


if __name__=="__main__":
    output_dir="/gpfs/alpine/scratch/mojia/phy171/MaCh3_results/JointAtmFit/job_jointfit_May17"
    
    chains=glob(output_dir+"/chain_*")
    chains_job_0 = [c for c in chains if re.search(r'chain_(9[0-5]|[1-8][0-9]|[0-9])$',c)]
    chains_job_1 = [c for c in chains if re.search(r'chain_(19[0-1]|1[0-8][0-9]|9[6-9])$',c)]
    chains_job_2 = [c for c in chains if re.search(r'chain_(28[0-7]|2[0-7][0-9]|19[2-9])$',c)]
    chains_job_3 = [c for c in chains if re.search(r'chain_(38[0-3]|3[0-7][0-9]|29[0-9]|28[8-9])$',c)]

    change_configs(chains_job_3,2)
