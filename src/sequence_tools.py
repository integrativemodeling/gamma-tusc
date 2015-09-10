from modeller import *
import yaml

def get_seq_from_aln(aln,seqnum):
    seq=''
    for p in aln.positions:
        if p.get_residue(aln[seqnum])==None:
            seq+='-'
        else:
            seq+=p.get_residue(aln[seqnum]).code
    return seq

def get_seqs_from_pir(aln_fn):
    """Manually retrieve exact sequences from alignment file.
    Returns dictionary, keys are template codes
    """
    inf = open(aln_fn)
    seqs = {}
    cur_seq = ''
    start_seq = False
    for l in inf:
        if l[0]=='>':
            if cur_seq!='':
                seqs[seqname] = cur_seq.strip('*')
                cur_seq = ''
            seqname = l.strip('\n')[4:]
            start_seq = True
            continue
        if start_seq:
            start_seq = False
            continue
        else:
            cur_seq+=l.strip('\n')
    if cur_seq!='':
        seqs[seqname] = cur_seq.strip('*')
        cur_seq = ''
    return seqs

def get_seqs_from_fasta(aln_fn):
    inf = open(aln_fn)
    seqs = {}
    cur_seq = ''
    start_seq = False
    for l in inf:
        if l[0]=='>':
            if cur_seq!='':
                seqs[seqname] = cur_seq
                cur_seq = ''
            seqname = l.strip('\n')[1:]
            continue
        else:
            cur_seq+=l.strip('\n')
    if cur_seq!='':
        seqs[seqname] = cur_seq
        cur_seq = ''
    return seqs


def write_align(aln,prefix):
    aln.write(prefix+'.pir',alignment_format="PIR")
    aln.write(prefix+'.pap',alignment_format="PAP")

def pir_to_pap(env,prefix):
    aln = alignment(env,file=prefix+'.pir')
    write_align(aln,prefix)

def append_second_model_options_to_first(opts1,opts2):
    for sse_type in ['helices','sheets','loops']:
        opts1['sses'][sse_type]+=opts2['sses'][sse_type]
    opts1['insertions']+=opts2['insertions']

class SequenceShifter:
    """A class to shift sequence numbers based on an alignment-to-orig sequence.
    Also stores list of insertions (starting point in final sequence and length)"""
    def __init__(self,aln_to_orig):
        self.orig_to_final_dict = {}
        seq0 = get_seq_from_aln(aln_to_orig,0)
        seq1 = get_seq_from_aln(aln_to_orig,1)
        n_final = 1
        self.insertions = []
        cur_ins = [1,0]
        on_insert = False

        # make sequence map and store insertions
        for n,s1 in enumerate(seq1):
            if s1!='-':
                self.orig_to_final_dict[n+1] = n_final
                n_final+=1
                if on_insert:
                    if cur_ins[0]!=0:
                        self.insertions.append(cur_ins)
                    on_insert = False
            else:
                if on_insert:
                    cur_ins[1]+=1
                else:
                    cur_ins = [n_final-1,1]
                    on_insert = True

        # don't add C-term insert
        #if on_insert:
        #    self.insertions.append(cur_ins)
        self.final_len = n_final-1

    def orig_to_final(self,resnum):
        if resnum not in self.orig_to_final_dict:
            print('WARNING: %i not in final sequence'%resnum)
            return -1
        else:
            return self.orig_to_final_dict[resnum]

    def get_insertions(self):
        """Format: [position (in final seq), length]"""
        return self.insertions

    def get_final_length(self):
        return self.final_len

    def shift_model_options(self,opts,offset=0):
        """Return a new ModelOptions with adjusted sequences
        AND ADDED INSERTIONS!
        Doesn't copy symmetries"""
        ret = ModelOptions()
        rsse = ret.get_sses()
        rins = ret.get_insertions()
        for h in opts.get_sses()['helices']:
            rsse['helices'].append([self.orig_to_final(h[0])+offset,
                                    self.orig_to_final(h[1])+offset])
        for h in opts.get_sses()['loops']:
            rsse['loops'].append([self.orig_to_final(h[0])+offset,
                                  self.orig_to_final(h[1])+offset])
        for sheet in opts.get_sses()['sheets']:
            tmp = []
            for h in sheet:
                tmp.append([self.orig_to_final(h[0])+offset,
                            self.orig_to_final(h[1])+offset])
            rsse['sheets'].append(tmp)
        ret.insertions = self.get_insertions()
        return ret



class ModelOptions:
    """A class to store some options for modelling.
    Allowed info: sses, insertions, symmetries
    Allowed sses: helices (lists), sheets (list of lists), loops (lists)
    Allowed insertions: [ [ins position,length] ... ]
    Allowed symmetries: [ [range1,range2], [range3,range4] ] where each range is a symmetry copy
    """
    def __init__(self,fn=''):
        self.sses = {'helices':[],
                     'sheets':[],
                     'loops':[]}
        self.insertions = []
        self.symmetries = []
        if fn!='':
            self.read(fn)

    def read(self,fn):
        cats = set(['sses','insertions','symmetries'])
        sse_cats = set(['helices','sheets','loops'])
        inf = open(fn,'r')
        self._data = yaml.load(inf)
        if not set(self._data) <= cats:
            raise Exception("ERROR: Weird categories in options file")
        if 'sses' in self._data:
            if not set(self._data['sses']) <= sse_cats:
                raise Exception("ERROR: Weird categories in SSEs of options file")
            for cat in self.sses:
                if cat in self._data['sses'] and self._data['sses'][cat] is not None:
                    self.sses[cat]+=self._data['sses'][cat]

    def get_sses(self):
        return self.sses

    def get_data(self):
        return {'sses':self.sses,
                'insertions':self.insertions,
                'symmetries':self.symmetries}

    def get_all_sse_residues(self):
        """Return all SSE residues in a flat set"""
        all_res = set()
        for segment in self.sses['helices']+self.sses['loops']:
            all_res|=set(range(segment[0],segment[1]+1))
        for sheet in self.sses['sheets']:
            for segment in sheet:
                all_res|=set(range(segment[0],segment[1]+1))
        return all_res

    def get_insertions(self):
        return self.insertions

    def get_symmetry_breaks(self):
        return self.symmetry_breaks

    def append(self,mo,offset=0):
        """Append another ModelOptions (sses and insertions only)
        With optional offset."""
        ret = ModelOptions()
        for h in mo.get_sses()['helices']:
            self.sses['helices'].append([h[0]+offset,
                                    h[1]+offset])
        for h in mo.get_sses()['loops']:
            self.sses['loops'].append([h[0]+offset,
                                  h[1]+offset])
        for sheet in mo.get_sses()['sheets']:
            tmp = []
            for h in sheet:
                tmp.append([h[0]+offset,
                            h[1]+offset])
            self.sses['sheets'].append(tmp)
        for ins in mo.get_insertions():
            self.insertions.append([ins[0]+offset,ins[1]])

    def write(self,out_fn):
        stream = file(out_fn,'w')
        yaml.dump(self.get_data(),stream)

class InsertionRemover:
    """Store a template and a target, remove insertions above certain size.
    Can modify what's removed with: limit_to_pdb(), enable_sses()
    Use get_alignment() to get final alignment
    Use get_current_target_alignment() for an alignment of final target to original
    """
    class _Position:
        def __init__(self,ntemp,temp,ntarg,targ):
            self.ntemp_orig = ntemp
            self.ntemp = ntemp
            self.temp = temp
            self.ntarg_orig = ntarg
            self.ntarg = ntarg
            self.targ = targ
            self.do_not_delete = False
            self.to_delete = False
        def __repr__(self):
            return '%s%i %s%i'%(self.temp,self.ntemp_orig,self.targ,self.ntarg_orig)

    def __init__(self,aln,max_ins_len): #,nterm_buffer=0):
        """Any insertions above max_ins_len are removed"""
        template_pos = get_seq_from_aln(aln,0)
        target_pos = get_seq_from_aln(aln,1)
        self.template_orig_seq = ''.join([t.code for t in aln[0].residues])
        self.target_orig_seq = ''.join([t.code for t in aln[1].residues])
        self.orig_aln = aln
        self.positions = []
        self.max_ins_len = max_ins_len
        self.pdb_aln = None
        #self.nterm_buffer = nterm_buffer
        ntemp= 0
        ntarg = 0
        for temp,targ in zip(template_pos,target_pos):
            if temp!='-':
                ntemp+=1
            if targ!='-':
                ntarg+=1
            self.positions.append(self._Position(ntemp,temp,ntarg,targ))

    def limit_to_pdb(self,pdb_aln):
        """Replace template residues with gaps if they are not in the PDB alignment.
        Expecting first sequence in alignment to be the full sequence, second to be PDB"""
        self.pdb_aln = pdb_aln
        complete_pos = get_seq_from_aln(pdb_aln,0)
        pdb_pos = get_seq_from_aln(pdb_aln,1)
        complete_seq = ''.join([t.code for t in pdb_aln[0].residues])
        if complete_seq!=self.template_orig_seq:
            raise Exception("ERROR: The full sequence in pdb_aln does not match aln")
        bad_template_indexes = set()
        ntemp = 1
        for p,o in zip(pdb_pos,complete_pos):
            if p=='-':
                bad_template_indexes.add(ntemp)
            ntemp+=1
        self._set_to_insertions(bad_template_indexes)
        self._update_indexes()

    def enable_sses(self,model_options):
        """For all SSEs, mark as non-insertion.
        Must provide model_options object, read with ModelOptions.read()"""
        self._set_do_not_delete(model_options.get_all_sse_residues())
        self._update_indexes()

    def get_insertions(self):
        """Return column numbers of insertions"""
        return self._update_to_delete()

    def get_alignment(self):
        """Create alignment from current positions after deleting insertions.
        This should ignore insertions at the beginning and end.
        """
        insertions = self._update_to_delete()
        temp_seq = ''
        targ_seq = ''
        on_delete = False
        for npos,pos in enumerate(self.positions):
            if pos.temp=='-' and pos.targ=='-':
                continue
            if pos.to_delete:
                if temp_seq!='':
                    on_delete = True
            else:
                if on_delete:
                    targ_seq += '/'
                temp_seq += pos.temp
                targ_seq += pos.targ
                on_delete = False

        aln = alignment(self.orig_aln.env)
        aln.append_sequence(temp_seq)
        aln.append_sequence(targ_seq)
        aln[0].code = self.orig_aln[0].code
        aln[0].name = self.orig_aln[0].name
        aln[1].code = self.orig_aln[1].code
        aln[1].name = self.orig_aln[1].name
        if self.pdb_aln:
            aln[0].range = self.pdb_aln[1].range
            aln[0].atom_file = self.pdb_aln[1].atom_file
            aln[0].source = self.pdb_aln[1].source
            aln[0].prottyp = self.pdb_aln[1].prottyp
            aln[0].resolution = self.pdb_aln[1].resolution
            aln[0].rfactor = self.pdb_aln[1].rfactor
        return aln

    def get_current_target_alignment(self):
        """Create alignment of current target to the original target sequence.
        This is useful if you need to map to the original indexes."""
        insertions = self._update_to_delete()
        orig_target = []
        new_target = []
        for pos in self.positions:
            if pos.targ!='-':
                orig_target.append(pos.targ)
                if not pos.to_delete:
                    new_target.append(pos.targ)
                else:
                    new_target.append('-')
        new_aln = alignment(self.orig_aln.env)
        new_aln.append_sequence(''.join(orig_target))
        new_aln.append_sequence(''.join(new_target))
        new_aln[0].code = "ORIG_TARGET"
        new_aln[1].code = "FINAL_TARGET"
        return new_aln

    def _update_to_delete(self):
        """Based on current status and max_ins_len, update what should be deleted.
        Also returns list of grouped insertion positions"""
        cur_ins = []
        insertions = []
        for npos,pos in enumerate(self.positions):
            pos.to_delete = False
            if pos.temp=='-' and pos.targ=='-':
                continue
            if pos.temp=='-' and not pos.do_not_delete:
                cur_ins.append(npos)
            else:
                if len(cur_ins)>self.max_ins_len:
                    insertions.append(cur_ins)
                    for nipos in cur_ins:
                        self.positions[nipos].to_delete = True
                cur_ins = []
        if len(cur_ins)>self.max_ins_len:
            insertions.append(cur_ins)
            for nipos in cur_ins:
                self.positions[nipos].to_delete = True
        return insertions

    def _update_indexes(self):
        """Based on current sequences, update indexes"""
        ntemp = 0
        ntarg = 0
        for pos in self.positions:
            if pos.temp!='-':
                ntemp+=1
            if pos.targ!='-':
                ntarg+=1
            pos.ntemp = ntemp
            pos.ntarg = ntarg

    def _set_to_insertions(self,template_indexes):
        """Based on raw template indexes, set positions to insertion."""
        for pos in self.positions:
            if pos.ntemp_orig in template_indexes:
                pos.temp = '-'

    def _set_do_not_delete(self,target_indexes):
        """Based on raw target indexes, set to DO NOT REMOVE"""
        for pos in self.positions:
            if pos.targ!='-' and pos.ntarg_orig in target_indexes:
                pos.do_not_delete = True

    def __repr__(self):
        ret = ''
        for pos in self.positions:
            ret+=pos.__repr__()+'\n'
        return ret
