#source .bashrc to setup ROOT environment before using

import os
import argparse

def merge_chain_iters(chain_id, itr_start, itr_stop):
    iters = []
    #rootfile_dir="/gpfs/alpine/scratch/mojia/phy171/MaCh3_results/JointAtmFit/job_jointfit_Jun13"
    rootfile_dir="/gpfs/alpine/scratch/mojia/phy171/MaCh3_results/JointAtmFit/job_data_jointfit_Jun28"
    
    #output_dir="/gpfs/alpine/proj-shared/phy171/Asimov_fit_chains"
    output_dir="/gpfs/alpine/proj-shared/phy171/Data_fit_chains"

    if(os.path.exists(output_dir)==False):
        os.system("mkdir "+output_dir)

    for i in range(itr_start, itr_stop+1):
        output = rootfile_dir+"/chain_"+str(chain_id)+"/output/MaCh3-Atmospherics-MCMC_Iter_"+str(i)+".root"
        if(os.path.exists(output)==False):
            raise Exception("no valid input root file!")
            
        else:
            print(output)
            iters.append(output)

    merged_chain = output_dir+"/MaCh3_MCMC_chain_"+str(chain_id)+"_Iter_"+str(itr_start)+"_"+str(itr_stop)+".root"
    #merged_chain = output_dir+"/MaCh3_MCMC_chain_"+str(chain_id)+".root" 
    haddComd = "hadd -f "+merged_chain+" "
    for itr in iters:
        haddComd += itr+" "

    print(haddComd)
    os.system(haddComd) 


if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--start',type = int, help="the first chain_id")
    parser.add_argument('-n','--num',type = int, help="the number of chains to be merged")

    args = parser.parse_args()

    print("chain_{0} to chain_{1}".format(args.start, args.num+args.start-1))
    start = args.start
    stop = args.start + args.num

    for i in range(start,stop):
        # Asimov chains:
        #if i<288:
        #    merge_chain_iters(i,0,9)
        #else:
        #    merge_chain_iters(i,0,8)
        
        # data fit chains:
        merge_chain_iters(i,0,7)

    print("DONE.")

