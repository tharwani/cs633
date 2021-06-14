#!/usr/bin/env python
# coding: utf-8

# In[47]:


import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import sys

sns.set()


# In[48]:


demo_input_format = pd.DataFrame.from_dict({
    "D": [],
    "P": [],
    "ppn": [],
    "mode": [],  # 1 --> optimized, 0 --> standard
    "time": [],
})


# In[49]:

collective = sys.argv[1]
f=open("output."+collective,'r')
lines=f.readlines()
data={}
for line in lines:
    k=line.split()
    data[k[0]]=list(map(float,k[1].split(',')))

for execution in range(10):
    for P in [4, 16]:
        for ppn in [1, 8]:
            for D in [16, 256,2048]:
                # Change with the actual data
                a1 = data.get(str(execution)+','+str(P)+','+str(ppn)+","+str(D),None)
                if a1 is not None:
                    demo_input_format = demo_input_format.append({
                        "D": D, "P": P, "ppn": ppn, "mode": 1, "time": data[str(execution)+','+str(P)+','+str(ppn)+","+str(D)][1]
                    }, ignore_index=True)
                    demo_input_format = demo_input_format.append({
                        "D": D, "P": P, "ppn": ppn, "mode": 0, "time": data[str(execution)+','+str(P)+','+str(ppn)+","+str(D)][0]
                    }, ignore_index=True)

demo_input_format["(P, ppn)"] = list(map(lambda x, y: ("(" + x + ", " + y + ")"), map(str, demo_input_format["P"]), map(str, demo_input_format["ppn"])))

# print(demo_input_format)

# In[50]:


y= sns.catplot(x="(P, ppn)", y="time", data=demo_input_format, kind="bar", col="D", hue="mode",capsize=.2)
y.set(yscale='log')
# plt.show()
plt.savefig('plot_'+collective+'.jpg')

# In[ ]:




