
output=[]

for i in open('output.pdb'):
    try:
        if(i.strip().split()[2]=='CA'):
            output.append(i)
    except:
        output.append(i)


fh=open('output2.pdb','w')
for i in output:
    fh.write(i)
fh.flush()
fh.close()