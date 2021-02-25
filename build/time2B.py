import subprocess
from subprocess import Popen, PIPE
import numpy as np
import re


if __name__ == "__main__":


    ns=np.logspace(np.log10(40),np.log10(30000),num=10)
    nSamples=1000 * 100**2 / ns**2
    print( "{} {} {} {} {}".format("N","method","timer","time" ,"nSamples") )

    
    for method in ["direct","acc_index","acc_copy"]:
        for n,nSample in zip(ns,nSamples):
            nSample=max(int(nSample),1)

            with Popen(['./timers', method ,str(int(n)) , str(int(nSample))], stdout=PIPE, universal_newlines=True) as process:
                for line in process.stdout:
                    
                    m=re.match("Timer (.+): (.+)",line)
                    if(m):
                        print( "{} {} {} {} {}".format(int(n),method,m[1],float(m[2])/nSample * 1000 ,nSample) )

        




