import subprocess
import os
from time import sleep
subprocess.Popen(['make','clean'])
sleep(2)
subprocess.Popen(['make'])
sleep(2)
# running the code
data = {
    "profile":[],
    "val":[],
}
for e in range(10):
    for P in [1,2]: #number of nodes
        for ppn in [1,2,4]:
            print("trying for P = " + str(P) + " ppn = " + str(ppn)) 
            subprocess.Popen([os.path.expanduser('~/UGP/allocator/src/allocator.out'),str(P*ppn),str(ppn)],stdout=subprocess.DEVNULL)
            sleep(2)
            try:
                temp = subprocess.check_output(['mpiexec', '-np', str(P*ppn),"-hostfile=hosts",'./code', "tdata.csv"],timeout=100) 
            except:
                try:
                    temp = subprocess.check_output(['mpiexec', '-np', str(P*ppn),"-hostfile=hostsimproved",'./code', "tdata.csv" ],timeout=100)
                except:
                    # print(temp)
                    print("bad")
                    exit()
                    break
            # temp = temp.decode("utf-8") 
            # temp = temp.split('\n')
            # print('ok')
            time = subprocess.check_output(["tail","-n",'1',"output.txt"],timeout=10)
            subprocess.check_output(['cp', 'output.txt', "outputp" + str(P) + "ppn" + str(ppn) + "epoch" + str(e) + ".txt" ],timeout=10)
            
            pro = "p = " + str(P) + ", ppn = " + str(ppn)
            print( str(pro) + "  = "  + str(float(time)))
            data["profile"].append(pro)
            data["val"].append(float(time))


# for collective in ['Bcast','Gather','Reduce','Alltoallv']:
#     subprocess.Popen(['python3','plot.py',collective])
# print("plots complete")
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt

# import matplotlib
# matplotlib.use('Agg')
plt.figure(figsize=(13,5))
df = pd.DataFrame.from_dict(data)
y = sns.boxplot(x = "profile", y = "val" ,data= df)
# y.set(yscale="log")
plt.savefig('plot.png')





