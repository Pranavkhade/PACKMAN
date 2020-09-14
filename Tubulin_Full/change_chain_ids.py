import string


output=[]

for numi,i in enumerate( open('chain_merged.pdb').read().split('TER') ):
    print(numi)
    if(i.find(' A ')!=-1):
        output.append(i.replace(' A ', ' '+string.printable[numi]+' ')+'TER\n')

fh=open('output.pdb','w')
for i in output:
    fh.write(i)
fh.flush()
fh.close()


#B to A chain
'''
output=[]
ter_flag=False
for numi,i in enumerate( open('cluster1_5sheet_pymol.pdb').readlines()):
    #440
    if(i.find(' A ')!=-1):
        output.append(i)

    if(i.find(' B ')!=-1):
        temp=i.replace(' B ',' A ')
        
        temp=temp.replace(temp[22:27], ' '+str(int(temp[22:26])+440)+' ')
        output.append(temp)
    
    if(i.find('GLU B 427')!=-1):
        ter_flag=True
        continue
    
    if(ter_flag):
        output.append(i)
        ter_flag=False
'''

fh=open('chain_merged.pdb','w')
for i in output:
    fh.write(i.strip())
fh.flush()
fh.close()