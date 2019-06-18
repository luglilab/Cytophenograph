import flowkit as fk
import glob

fnames = glob.glob('/mnt/hpcserver1_datadisk2_spuccio/SP008_Phenograph_BMT/phenograph/test_fcs/*.fcs')


sample = fk.Sample(fnames[0], subsample_count=None)

for i in range(len(sample.pns_labels)):
    print(sample.pns_labels[i])

