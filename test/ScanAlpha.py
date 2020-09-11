'''
Notes:
    * Do something about trash.txt
'''

from packman.apps import predict_hinge
from packman import molecule

import numpy
import argparse


def IO():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',  help='Name of the PDB file')
    parser.add_argument('--chain',     help='Chain ID')
    
    parser.add_argument('--begin',      default=0,
                        help='Starting Alpha Value (Recommended: 0.0)')
    parser.add_argument('--end',        default=10,
                        help='Ending Alpha Value (Recommended: 10.0)')
    parser.add_argument('--step_size',  default=0.5,
                        help='Progression amount from start alpha value to end alpha value')
    parser.add_argument('--outfilename', default='output.hng',
                        help='Name of the output file.')
    args = parser.parse_args()
    return args


def hinge_scanner(atoms,filename,alpha_start, alpha_stop, step_size):
    """
    """
    chain=atoms[0].get_parent().get_parent().get_id()

    for i in numpy.arange(alpha_start, alpha_stop, step_size):
        i = numpy.around(i, decimals=1)
        try:
            predict_hinge(atoms, Alpha=i, outputfile=open('trash.txt', 'w'))
        except:
            continue

    hinges = atoms[0].get_parent().get_parent().get_hinges()

    # Remove the insignificant hinges having p value less than 0.05

    pre_overlap = []
    for i in range(0, len(hinges)):
        for j in range(i+1, len(hinges)):
            h1 = [k.get_id() for k in hinges[i].get_elements()]
            h2 = [k.get_id() for k in hinges[j].get_elements()]
            if(hinges[i].get_pvalue() <= 0.05 and hinges[j].get_pvalue() <= 0.05 and len(set(h1).intersection(h2)) > 0):
                hinge_union = list(set(h1).union(h2))

                flag = True
                for k in range(0, len(pre_overlap)):
                    if(len(set(pre_overlap[k]).intersection(hinge_union)) > 0):
                        pre_overlap[k] = list(
                            set(pre_overlap[k]).union(hinge_union))
                        flag = False
                        break

                if(flag):
                    pre_overlap.append(hinge_union)

    overlap = []
    for i in range(0, len(pre_overlap)):

        flag = True
        for j in range(0, len(overlap)):
            if(len(set(overlap[j]).intersection(pre_overlap[i])) > 0):
                overlap[j] = list(set(overlap[j]).union(pre_overlap[i]))
                flag = False

        if(flag):
            overlap.append(pre_overlap[i])

    # Sort overlap
    overlap = sorted(overlap, key=lambda x: x[0])
    overlap = [sorted(i) for i in overlap]

    output=[]
    for i in range(0, len(overlap)):

        if(i == 0):
            output.append(filename+'_'+chain+'\tD'+str(i+1)+'\t1:'+str(overlap[i][0]-1) + '\n'
                              + filename+'_'+chain+'\tH' +
                              str(i+1)+'\t' +
                              str(overlap[i][0])+':'+str(overlap[i][-1])+'\n'
                              + filename+'_'+chain+'\tD'+str(i+2)+'\t'+str(overlap[i][-1]+1)+':'+str(overlap[i+1][0]-1)+'\n')
        else:
            try:
                output.append(filename+'_'+chain+'\tH'+str(i+1)+'\t'+str(overlap[i][0]) + ':'+str(overlap[i][-1]) + '\n'
                                  + filename+'_'+chain+'\tD'+str(i+2)+'\t'+str(overlap[i][-1]+1)+':'+str(overlap[i+1][0]-1)+'\n')
            except:
                output.append(filename+'_'+chain+'\tH'+str(i+1)+'\t'+str(overlap[i][0]) + ':'+str(overlap[i][-1]) + '\n'
                                  + filename+'_'+chain+'\tD'+str(i+2)+'\t'+str(overlap[i][-1]+1)+':'+'Inf\n')

    return output


def main():
    ARGS = IO()
    filename = ARGS.filename
    chain = ARGS.chain
    alpha_start = float(ARGS.begin)
    alpha_stop = float(ARGS.end)
    step_size = float(ARGS.step_size)

    if(ARGS.outfilename=='output.hng'):
        output_file = open(filename+'.hng', 'w', encoding='UTF-8')
    else:
        output_file = open(ARGS.outfilename, 'w', encoding='UTF-8')

    mol = molecule.load_structure(filename)


    if(chain==None):
        for i in [i.get_id() for i in mol[0].get_chains()]:
            backbone = [k for j in mol[0][i].get_backbone() for k in j if k is not None]
            output=hinge_scanner(backbone,filename,alpha_start,alpha_stop,step_size)
            for j in output:output_file.write(j)
        
    else:
        backbone = [j for i in mol[0][chain].get_backbone() for j in i if j is not None]
        output=hinge_scanner(backbone,filename,alpha_start,alpha_stop,step_size)
        for i in output:output_file.write(i)
        
    output_file.flush()
    output_file.close()
    
    return True


if(__name__ == '__main__'):
    main()
