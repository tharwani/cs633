import subprocess
import os
from time import sleep
subprocess.Popen(['make','clean'])
sleep(2)
subprocess.Popen(['make'])
sleep(2)
f1 = open("output.Bcast","w")
f2 = open("output.Gather","w")
f3 = open("output.Reduce","w")
f4 = open("output.Alltoallv","w")
# running the code
for e in range(10):
    for P in [4,16]: #number of nodes
        for ppn in [1,8]: 
            subprocess.Popen([os.path.expanduser('~/UGP/allocator/src/allocator.out'),str(P*ppn),str(ppn)],stdout=subprocess.DEVNULL)
            sleep(2)
            for D in [16,256,2048]:
                try:
                    temp = subprocess.check_output(['mpiexec', '-np', str(P*ppn),"-hostfile","hosts",'./src', str(D) ],timeout=200) 
                except:
                    try:
                        temp = subprocess.check_output(['mpiexec', '-np', str(P*ppn),"-hostfile","hostsimproved",'./src', str(D) ],timeout=200)
                    except:
                        break
                temp = temp.decode("utf-8") 
                temp = temp.split('\n')
                print('ok')
                print(str(e)+','+str(P)+','+str(ppn)+","+str(D)+' '+temp[0],file=f1)
                print(str(e)+','+str(P)+','+str(ppn)+","+str(D)+' '+temp[1],file=f2)
                print(str(e)+','+str(P)+','+str(ppn)+","+str(D)+' '+temp[2],file=f3)
                print(str(e)+','+str(P)+','+str(ppn)+","+str(D)+' '+temp[3],file=f4)
f1.close()
f2.close()
f3.close()
f4.close()

for collective in ['Bcast','Gather','Reduce','Alltoallv']:
    subprocess.Popen(['python3','plot.py',collective])
print("plots complete")
