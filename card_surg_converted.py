import re
import numpy as np
import os
import matplotlib.pyplot as plt
import csv
pad=30
size=(4245,pad+1)


##rows = p.split('\n')
##for row in rows:
##    ll.append(re.findall('[\d+-\.e]+', row))
##np.array(ll, dtype=float)
lt = []
lo = []
lp = []
lta = np.empty(size)
lta[:] = np.nan
lto = np.empty(size)
lto[:] = np.nan
ltp = np.empty(size)
ltp[:] = np.nan
fups=np.zeros(32)
loc = os.path.expanduser('~/Documents/dlnmRevData/deidentified_Mar8.csv')
with open(loc) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            lt=re.findall('[\d+-\.e]+', row[2])
            if len(lt)<=pad:
                lta[line_count-1,0:len(lt)]=np.array(lt, dtype=float)
                fups[len(lt)]+=int(row[1])
            else:
                print("t"+str(line_count))
                lta[line_count-1,:-1]=np.array(lt[0:pad], dtype=float)
                fups[pad+1]+=int(row[1])
            lo=re.findall('[\d+-\.e]+', row[3])
            if len(lo)<=pad:
                lto[line_count-1,0:len(lo)]=np.array(lo, dtype=float)
            else:
                print("o"+str(line_count))
                lto[line_count-1,:-1]=np.array(lo[0:pad], dtype=float)
            lp=re.findall('[\d+-\.e]+', row[4])
            if len(lp)<=pad:
                ltp[line_count-1,0:len(lp)]=np.array(lp, dtype=float)
            else:
                print("p"+str(line_count))
                ltp[line_count-1,:-1]=np.array(lp[0:pad], dtype=float)
            lta[line_count-1,pad]=int(row[0])
            lto[line_count-1,pad]=int(row[0])
            ltp[line_count-1,pad]=int(row[0])
            line_count += 1
    print(f'Processed {line_count} lines.')
headera=""
for i in range(30):
    headera=headera+'t'+ str(i)+','
headera=headera+'ID'
np.savetxt('~/Documents/dlnmRevData/Mar8_envs_tem.csv',lta,delimiter = ',',header=headera)
np.savetxt('~/Documents/dlnmRevData/Mar8_envs_o3.csv',lto,delimiter = ',',header=headera)
np.savetxt('~/Documents/dlnmRevData/Mar8_envs_pm25.csv',ltp,delimiter = ',',header=headera)
print(fups)
plt.plot(range(0,32),fups,'bo')
plt.title("# of post-surgical days for mortality/readmission pts")
plt.show()
