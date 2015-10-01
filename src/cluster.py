import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.tools
from mpi4py import MPI
import sys,os,shutil

def splitall(path):
    allparts = []
    while 1:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
number_of_processes = comm.size
mdl = IMP.Model()

file_prefix = sys.argv[1]
num_files = int(sys.argv[2])
fname = sys.argv[3]
num_clusters = int(sys.argv[4])
out_dir = sys.argv[5]
display_plot = True


try:
    os.mkdir(out_dir)
except:
    pass

file_list = [os.path.join(file_prefix,str(n),fname) for n in range(1,num_files+1) if os.path.isfile(os.path.join(file_prefix,str(n),fname))]
my_file_list = IMP.pmi.tools.chunk_list_into_segments(file_list,number_of_processes)[rank]
my_coords = []
print 'rank',rank,'len files',len(my_file_list)
for n,fn in enumerate(my_file_list):
    mh = IMP.atom.read_pdb(fn,mdl,IMP.atom.CAlphaPDBSelector())
    coords = [list(IMP.core.XYZ(p).get_coordinates()) for p in IMP.core.get_leaves(mh)]
    my_coords.append({'gtusc':coords})
    IMP.atom.destroy(mh)
    del mh
    print 'rank',rank,'%.2f'%(float(n+1)/len(my_file_list)*100),'% done'

print 'rank',rank,'scatter 1'
gathered = IMP.pmi.tools.scatter_and_gather(zip(my_file_list,my_coords))
#print 'rank',rank,'scatter 2'
#gather_names = IMP.pmi.tools.scatter_and_gather(my_file_list)
print 'rank',rank,'cluster fill'
clusters = IMP.pmi.analysis.Clustering()
for name,coords in gathered:
    clusters.fill(name,coords)

print "Global calculating distance matrix"
clusters.dist_matrix()
if rank==0:
    clusters.do_cluster(num_clusters)
    if display_plot:
        clusters.plot_matrix(figurename=out_dir+"/clusters.png")
    clusters.save_distance_matrix_file(out_dir+"/clusters.dat")
    for n, cl in enumerate(clusters.get_cluster_labels()):
        print "cluster %s " % str(n)
        print "cluster label %s " % str(cl)
        print clusters.get_cluster_label_names(cl)

        dircluster = out_dir + "/cluster." + str(n) + "/"
        try:
            os.mkdir(dircluster)
        except:
            pass
        rmsd_dict = {"AVERAGE_RMSD":
                     str(clusters.get_cluster_label_average_rmsd(cl))}
        clusstat = open(dircluster + "stat.out", "w")
        clusstat.write(str(rmsd_dict)+'\n')
        for k, structure_name in enumerate(clusters.get_cluster_label_names(cl)):
            allparts = splitall(structure_name)
            print allparts[-2],structure_name
            out_name = os.path.join(dircluster,allparts[-2]+'.pdb')
            shutil.copy(structure_name,out_name)
        clusstat.close()
