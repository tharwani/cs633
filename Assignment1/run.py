
import subprocess
import os
from time import sleep
subprocess.Popen(['make','clean'])
sleep(2)
subprocess.Popen(['make'])
sleep(2)
output = {}
data = {}
# due to error after make call ,try is used. For further explaintion, refer to Readme.pdf Issues faced section
try:
    f = open("output.data","w")
    # running the code
    for e in [1,2,3,4,5]:
        for p in [4,6,7,8]:
            subprocess.Popen([os.path.expanduser('~/UGP/allocator/src/allocator.out'),str(p*p),'8'],stdout=subprocess.DEVNULL)
            for n in [16,32,64,128,256,512,1024]:
                try:
                    temp = subprocess.check_output(['mpiexec', '-np', str(p*p),"-hostfile","hosts",'./a.out', str(p), str(n), str(50) ],timeout=25) 
                except:
                    temp = subprocess.check_output(['mpiexec', '-np', str(p*p),"-hostfile","hostsimproved",'./a.out', str(p), str(n), str(50) ])
                temp = temp.decode("utf-8") 
                temp = list(map(float,temp.split('\n')[0].split(',')))
                output[str(e)+":"+str(p)+":"+str(n)] = temp
                print('ok')
                print("iteration no : "+str(e)+",no of processes : "+str(p)+"^2,size : "+str(n)+"^2"+"\ntime1 :"+str(temp[0])+"\ntime2 :"+str(temp[1])+"\ntime3 :"+str(temp[2]),file=f)
    f.close()
except:
    sleep(5)
    f = open("output.data","w")
    # running the code
    for e in [1,2,3,4,5]:
        for p in [4,6,7,8]:
            subprocess.Popen([os.path.expanduser('~/UGP/allocator/src/allocator.out'),str(p*p),'8'],stdout=subprocess.DEVNULL)
            for n in [16,32,64,128,256,512,1024]:
                try:
                    temp = subprocess.check_output(['mpiexec', '-np', str(p*p),"-hostfile","hosts",'./a.out', str(p), str(n), str(50) ],timeout=25) 
                except:
                    temp = subprocess.check_output(['mpiexec', '-np', str(p*p),"-hostfile","hostsimproved",'./a.out', str(p), str(n), str(50) ])
                temp = temp.decode("utf-8") 
                temp = list(map(float,temp.split('\n')[0].split(',')))
                output[str(e)+":"+str(p)+":"+str(n)] = temp
                print('ok')
                print("iteration no : "+str(e)+",no of processes : "+str(p)+"^2,size : "+str(n)+"^2"+"\ntime1 :"+str(temp[0])+"\ntime2 :"+str(temp[1])+"\ntime3 :"+str(temp[2]),file=f)
    f.close()

for p in [4,6,7,8]:
    for n in [16,32,64,128,256,512,1024]:
        temp1 = []
        temp2 = []
        temp3 = []
        for e in [1,2,3,4,5]:
            temp1.append(output[str(e)+":"+str(p)+":"+str(n)][0]) 
            temp2.append(output[str(e)+":"+str(p)+":"+str(n)][1]) 
            temp3.append(output[str(e)+":"+str(p)+":"+str(n)][2]) 
        data[str(p)+":"+str(n)] = [temp1,temp2,temp3]


# graph plotting part
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
fig, ax = plt.subplots(figsize=(30, 10))
bps = []
def box_plot(data, labels, edge_color, fill_color, positions, color_label):
    bp = ax.boxplot(data, positions = positions ,labels = labels, widths = 0.075, patch_artist=True)    
    bps.append(bp)
    ax.set_yscale('log')
    ax.set_ylabel('time')
    ax.set_xlabel('array size per process')
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bp[element], color=edge_color)

    for patch in bp['boxes']:
        patch.set(facecolor=fill_color)  
color = [['red', 'tan'],['blue', 'cyan'],['darkgreen','lightgreen']]
pos = [[1, 1.6, 2.2, 2.8, 3.4, 4.0, 4.6],[1.2, 1.8, 2.4, 3.0, 3.6, 4.2, 4.8],[1.4, 2.0, 2.6, 3.2, 3.8, 4.4, 5.0]]
for p in [4,6,7,8]:
    for method in range(3):
        temp = []
        label = []
        for n in [16,32,64,128,256,512,1024]:
            temp.append(data[str(p)+":"+str(n)][method])
            label.append(str(n)+"*"+str(n))
        box_plot(temp, label,color[method][0], color[method][1],pos[method], method+1)
    ax.legend([bp["boxes"][0] for bp in bps] , ["method 1","method 2","method 3" ], loc='upper right')
    fig.savefig('plot'+str(p*p)+'.png', dpi=fig.dpi)
    fig, ax = plt.subplots(figsize=(30, 10))

