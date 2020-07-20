import os


output=[]

for i in os.listdir():
    if(i[-3:]=='hng'):
        for j in open(i):
            temp=j.split()
            if(temp[1][0]=='H'):
                output.append(j)


fh=open('All_Hinges.txt','w',encoding='UTF-8')
for i in output:
    fh.write(i)
fh.flush()
fh.close()