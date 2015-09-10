import sys
sys.path.append('../src/')
from sequence_tools import *
from modeller import *
import unittest
import tempfile,os

def get_temp_aln(aln):
    tfname = os.path.join(tempfile._get_default_tempdir(),
                          next(tempfile._get_candidate_names()))
    aln.write(tfname,alignment_format='PIR')
    seqs = get_seqs_from_pir(tfname)
    os.unlink(tfname)
    return seqs

class TestInsertionRemover(unittest.TestCase):
    def test_basic_remover(self):
        template = '-ABC---DEF--'
        target   = 'YABCYYYDEFYY'
        env = environ()
        aln = alignment(env)
        aln.append_sequence(template)
        aln.append_sequence(target)
        aln[0].code = 'TEMPLATE'
        aln[1].code = 'TARGET'

        iremove = InsertionRemover(aln,2)
        fix_aln = iremove.get_alignment()
        seqs = get_temp_aln(fix_aln)

        self.assertEqual(seqs['TEMPLATE'],'-ABC-DEF--')
        self.assertEqual(seqs['TARGET'],  'YABC/DEFYY')

        temp_aln = iremove.get_current_target_alignment()
        self.assertEqual(get_seq_from_aln(temp_aln,0),'YABCYYYDEFYY')
        self.assertEqual(get_seq_from_aln(temp_aln,1),'YABC---DEFYY')

    def test_pdb_adjust(self):
        template  = 'ABC--DEF--'
        target    = 'ABCDEFGHKL'
        env = environ()
        aln = alignment(env)
        aln.append_sequence(template)
        aln.append_sequence(target)
        aln[0].code='TEMP'
        aln[1].code='TARG'

        temp_orig = 'ABCDEF'
        temp_pdb  = 'ABC-EF'
        pdb_aln = alignment(env)
        pdb_aln.append_sequence(temp_orig)
        pdb_aln.append_sequence(temp_pdb)

        iremove = InsertionRemover(aln,2)
        iremove.limit_to_pdb(pdb_aln)
        fix_aln = iremove.get_alignment()
        seqs = get_temp_aln(fix_aln)
        self.assertEqual(seqs['TEMP'],'ABC-EF--')
        self.assertEqual(seqs['TARG'],'ABC/GHKL')

        temp_aln = iremove.get_current_target_alignment()
        self.assertEqual(get_seq_from_aln(temp_aln,0),'ABCDEFGHKL')
        self.assertEqual(get_seq_from_aln(temp_aln,1),'ABC---GHKL')
    def test_pdb_adjust_deletion(self):
        """Test adjusting when you already have a partial deletion"""
        template  = 'EEFFIKQGPSSGNVSAQPEEDEEDLGIGGLTGKQLRELQDLRLIEEENMLAPSLK'
        target    = 'GEFFIAENTDTNG---------------------------------TDDDFIYHI'
        env = environ()
        aln = alignment(env)
        aln.append_sequence(template)
        aln.append_sequence(target)
        aln[0].code='TEMP'
        aln[1].code='TARG'

        temp_orig = 'EEFFIKQGPSSGNVSAQPEEDEEDLGIGGLTGKQLRELQDLRLIEEENMLAPSLK'
        temp_pdb  = 'EEFFIKQG--------------------------------------------SLK'
        pdb_aln = alignment(env)
        pdb_aln.append_sequence(temp_orig)
        pdb_aln.append_sequence(temp_pdb)

        # first just test the PDB part
        ir1 = InsertionRemover(aln,1000)
        ir1.limit_to_pdb(pdb_aln)
        fix_aln = ir1.get_alignment()
        seqs = get_temp_aln(fix_aln)
        self.assertEqual(seqs['TEMP'],'EEFFIKQG-----------SLK')
        self.assertEqual(seqs['TARG'],'GEFFIAENTDTNGTDDDFIYHI')

        # now remove the insert
        ir2 = InsertionRemover(aln,6)
        ir2.limit_to_pdb(pdb_aln)
        fix_aln = ir2.get_alignment()
        seqs = get_temp_aln(fix_aln)
        self.assertEqual(seqs['TEMP'],'EEFFIKQG-SLK')
        self.assertEqual(seqs['TARG'],'GEFFIAEN/YHI')

        #temp_aln = iremove.get_current_target_alignment()
        #self.assertEqual(get_seq_from_aln(temp_aln,0),'ABCDEFGHKL')
        #self.assertEqual(get_seq_from_aln(temp_aln,1),'ABC---GHKL')

    def test_real_problem(self):
        env = environ()
        aln = alignment(env,file='input/test_aln.pir')
        pdb_aln = alignment(env,file='input/test_pdblim.pir')

        iremove = InsertionRemover(aln,max_ins_len=5)
        iremove.limit_to_pdb(pdb_aln)
        fix_aln = iremove.get_alignment()
        seqs = get_temp_aln(fix_aln)
        self.assertEqual(seqs['GCP4_HUMAN'],'MIHELLLALSGYPGSIFTWNKR-SGLQVSQD'
                         'FPFLHPSETSVLNRLCRLGTDYIRFTEFIEQYTGHGGLHG')
        self.assertEqual(seqs['GCP2_YEAST'],'VVKDLLNVLIGLEGTYIRYFND/IEFKIAKK'
                         'M---DPSFKTFSRRIVRYGKQYMILTRAYEKWSDT--SFG')


    def test_real_sses(self):
        env = environ()
        aln = alignment(env,file='input/test_aln.pir')
        pdb_aln = alignment(env,file='input/test_pdblim.pir')
        mo = ModelOptions('input/options.txt')

        sses = mo.get_sses()
        self.assertEqual(sses['helices'],[[56,67]])
        self.assertEqual(sses['sheets'],[[[86,89]]])
        self.assertEqual(sses['loops'],[[55,55]])
        self.assertEqual(mo.get_all_sse_residues(),set(range(55,68)+range(86,90)))

        iremove = InsertionRemover(aln,max_ins_len=3)
        iremove.limit_to_pdb(pdb_aln)
        iremove.enable_sses(mo)
        fix_aln = iremove.get_alignment()

        seqs = get_temp_aln(fix_aln)
        self.assertEqual(seqs['GCP4_HUMAN'],'----MIHELLLALSGYPGSIFTWNKR-----'
                         'SGLQVSQDFPFLHPSETSVLNRLCRLGTDYIRFTEFIEQYTGHGGLHG')
        self.assertEqual(seqs['GCP2_YEAST'],'QEALVVKDLLNVLIGLEGTYIRYFND/PETP'
                         'IEFKIAKKM---DPSFKTFSRRIVRYGKQYMILTRAYEKWSDT--SFG')

    def test_shift_sequence(self):
        """Test shifting residue numbers and ModelOptions"""
        template = '-ABC---DEF--'
        target   = 'YABCYYYDEFYY'
        env = environ()
        aln = alignment(env)
        aln.append_sequence(template)
        aln.append_sequence(target)
        aln[0].code = 'TEMPLATE'
        aln[1].code = 'TARGET'

        iremove = InsertionRemover(aln,2)
        aln_to_orig = iremove.get_current_target_alignment()
        #self.assertEqual(get_seq_from_aln(temp_aln,0),'YABCYYYDEFYY')
        #self.assertEqual(get_seq_from_aln(temp_aln,1),'YABC---DEFYY')
        ss = SequenceShifter(aln_to_orig)
        self.assertEqual(ss.orig_to_final(4),4)
        self.assertEqual(ss.orig_to_final(5),-1) # should get warning
        self.assertEqual(ss.orig_to_final(8),5)
        self.assertEqual(ss.get_insertions(),[[4,3]])

        mo = ModelOptions()
        mo.get_sses()['helices'].append([8,10])
        mo.add_inserts_from_shifter(ss)
        nmo = ss.shift_model_options(mo,offset=10)
        self.assertEqual(nmo.get_sses()['helices'][0],[15,17])
        self.assertEqual(nmo.get_insertions(),[[14,3]])

if __name__=="__main__":
    unittest.main()
