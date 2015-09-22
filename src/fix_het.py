import sys
def run(in_fn,out_fn):
    inf=open(in_fn,'r')
    outf=open(out_fn,'w')
    for line in inf.readlines():
        if line[0:4]=='ATOM':
            new_line=line[0:17]
            if line[17:20]=='MSE':
                new_line+='MET'
            elif line[17:20]=='HSD':
                 new_line+='HIS'
            else:
                new_line+=line[17:20]
            new_line+=line[20:]
            outf.write(new_line)
        elif line[0:6]=='HETATM':
            if line[17:20]=='HSD':
                new_line='ATOM  '+line[6:17]+'HIS'+line[20:]
                outf.write(new_line)
            else:
                outf.write(line)
        elif line.split()[0:5]==['REMARK','6','MODELLER','BLK','RESIDUE']:
            continue
        else:
            outf.write(line)
    inf.close()
    outf.close()
if __name__=="__main__":
    if len(sys.argv)!=3:
        print "Fixes weird residues. Currently only MSE->MET and HSD->HIS"
        print "USAGE <in_fn> <out_fn>"
    else:
        run(*sys.argv[1:])
