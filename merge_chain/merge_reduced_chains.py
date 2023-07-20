#source .bashrc to setup ROOT environment before using

import os

if __name__=="__main__":

    #file_dir = "/gpfs/alpine/proj-shared/phy171/Asimov_fit_chains"
    file_dir = "/gpfs/alpine/proj-shared/phy171/Data_fit_chains"
 
    chains=[]

    for i in range(384):
        #chain=file_dir+"/MaCh3_MCMC_chain_"+str(i)+"_reduced.root"
        chain=file_dir+"/MaCh3_MCMC_chain_"+str(i)+"_Iter_0_7_reduced.root"

        if os.path.exists(chain)==False:
            raise Exception("No valid input root file!")
        else:
            chains.append(chain)
            print(chain)
    
    merged_chain = file_dir+"/MaCh3_MCMC_Iter_0_7_reduced.root"
    haddComd = "hadd -f "+merged_chain+" "
    for ch in chains:
        haddComd += ch+" "

    print(haddComd)

    os.system(haddComd)

    print("DONE.")

