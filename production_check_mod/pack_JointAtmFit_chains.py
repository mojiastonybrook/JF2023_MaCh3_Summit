import os
from glob import glob

def copy_chain_output(chain_id, itr_end):
    output_dir="/gpfs/alpine/scratch/mojia/phy171/MaCh3_results/JointAtmFit/job_jointfit_May17"
    # make directory
    chain_dir = "chain_"+str(chain_id)
    if(os.path.exists("./"+chain_dir)==False):
        os.system("mkdir "+chain_dir)
        os.system("mkdir "+chain_dir+"/output")
    # copy config
    cp_command = "cp -r "+output_dir+"/chain_"+str(chain_id)+"/config ./"+chain_dir+"/"
    os.system(cp_command)
    # copy outout
    for i in range(0,itr_end+1):
        root_file=glob(output_dir+"/"+chain_dir+"/output/MaCh3-Atmospherics-MCMC_Iter_"+str(i)+".root")[0]
        cp_command = "cp "+root_file+" ./"+chain_dir+"/output/"
        os.system(cp_command)

if __name__=="__main__":

    pack_dir="/gpfs/alpine/phy171/scratch/mojia/MaCh3_results/jointfit_packed"
     
    os.chdir(pack_dir)
    #chains in job_0
    #for i in range(0,96):
        #copy_chain_output(i,7)
    copy_chain_output(7,7)
    copy_chain_output(8,7)
    copy_chain_output(9,7)
    #pack chains
    print("packing chains...")
    pack_name="chain_7_9"

    tar_command="tar -czvf jointAtmFit_"+pack_name+".tar.gz chain_*"
    os.system(tar_command)
    print("packing done, delete files")
    del_command="rm -rf chain_*"
    os.system(del_command)

    print("DONE!")
        
