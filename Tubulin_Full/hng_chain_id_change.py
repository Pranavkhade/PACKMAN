import string

output=[]

hng_format=open('output.hng').read()

for i in string.printable:
    output.append(hng_format.replace('_A','_'+i))

fh=open('final_output.hng','w')
for i in output:
    fh.write(i)
fh.flush()
fh.close()