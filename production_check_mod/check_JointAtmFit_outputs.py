import os
import re
from glob import glob
import numpy as np

def take_input():
    job_id=int(input("which job?[0,1,2,3]: "))
    itr_id=int(input("which iteration?[0-11]: "))

    return {"job":job_id,
            "iter":itr_id
            }

if __name__=="__main__":
    output_dir="/gpfs/alpine/scratch/mojia/phy171/MaCh3_results/JointAtmFit/job_jointfit_Jun13"
    output_num = 1
    
    usr_input=take_input()
    itr = usr_input["iter"]
    job_id = usr_input["job"]

    chains=glob(output_dir+"/chain_*")
    chains_job_0 = [c for c in chains if re.search(r'chain_(9[0-5]|[1-8][0-9]|[0-9])$',c)]
    chains_job_1 = [c for c in chains if re.search(r'chain_(19[0-1]|1[0-8][0-9]|9[6-9])$',c)]
    chains_job_2 = [c for c in chains if re.search(r'chain_(28[0-7]|2[0-7][0-9]|19[2-9])$',c)]
    chains_job_3 = [c for c in chains if re.search(r'chain_(38[0-3]|3[0-7][0-9]|29[0-9]|28[8-9])$',c)]
    jobs=[chains_job_0,chains_job_1,chains_job_2,chains_job_3]

    print("chains' number: ", len(chains))
    print("job_0 chains: ", len(chains_job_0))
    print("job_1 chains: ", len(chains_job_1))
    print("job_2 chains: ", len(chains_job_2))
    print("job_3 chains: ", len(chains_job_3))
   
    fail_chains=[]
    partial_chains=[]
    c_chains=[]
    ac_steps=[]
    for chain in jobs[job_id]:
        mcs = glob(chain+"/output/MaCh3-Atmospherics-MCMC_Iter_"+str(itr)+".root")
        if len(mcs) < output_num:
            print("no output root file: ", chain, "!")
            fail_chains.append(chain)
            continue
        else:
            mc_log = chain+"/output/MaCh3-Atmospherics-MCMC_Iter_"+str(itr)+".root.log"
            with open(mc_log,'r') as f:
                last_two_lines = f.readlines()[-2:]
                last_line = last_two_lines[1]
                if last_line.split()[1] != "steps":
                    partial_chains.append(chain)
                else:
                    c_chains.append(chain)
                    ac_steps.append(int(last_line.split()[0]))
                    stl_line = last_two_lines[0] 
                    print(stl_line.split()[7].rstrip('s').strip('('))

    print("Failed chains: ")
    for fc in fail_chains: print(fc)

    print("complete chains: ")
    for cc in c_chains: print(cc)

    print("failed chains: ", len(fail_chains))
    print("partially chains: ", len(partial_chains))
    print("completed chains: ", len(c_chains))

    if len(ac_steps)!=0: print(np.array(ac_steps).mean())
