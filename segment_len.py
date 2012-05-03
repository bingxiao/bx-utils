#!/usr/bin/env python
''' basic statistics to paired-end mapping '''

import numpy as np
import pysam
import sys
import os

try:
    input_file = sys.argv[1]
    prefix = os.path.split( input_file )[-1]
    if not os.path.exists( input_file + '.bai' ):
        print >>sys.stderr, '.bai file missing, you have to index the bam file first:\n# samtools index %s'%input_file
        sys.exit()
    #samfile = pysam.Samfile( 'a.bam', 'rb' )
    samfile = pysam.Samfile( input_file, 'rb' )
except:
    print 'usage: %s <bam_file.bam>'%sys.argv[0]
    exit(0)

def get_seg_len( read1, read2 ):
    return max( read1.aend, read2.aend ) - min( read1.pos, read2.pos )

read1_dic = {}
read2_dic = {}

segment_length_lst = []
with open( prefix+'.seglen.all', 'w' ) as out:
    for n, read in enumerate( samfile.fetch() ): 
        # if read.is_read1 and
        if read.is_proper_pair and not read.is_qcfail:
            if read.is_read1:
                if read.qname in read2_dic: # good
                    seg_len = get_seg_len( read, read2_dic[ read.qname ] )
                    segment_length_lst.append( seg_len )
                    print >>out, '%s\t%d'%( read.qname, seg_len )
                    del read2_dic[ read.qname ]
                else:
                    read1_dic[ read.qname ] = read
            else:
                assert read.is_read2
                if read.qname in read1_dic:
                    seg_len = get_seg_len( read1_dic[ read.qname ], read )
                    segment_length_lst.append( seg_len )
                    print >>out, '%s\t%d'%( read.qname, seg_len )
                    del read1_dic[ read.qname ]
                else:
                    read2_dic[ read.qname ] = read
        if n % 100000 == 0:
            print >>sys.stderr, 'Parsed %d reads.\r'%n,

samfile.close()

segment_length_lst = np.array( segment_length_lst )
print '-'*32
print 'segment_length_lst array shape', segment_length_lst.shape
print 'mean:', segment_length_lst.mean()        ## these are polluted by gaps. infer it by distributions in R.
print 'std:', segment_length_lst.std()
print 'min:', segment_length_lst.min()
print 'max:', segment_length_lst.max()
print '-'*32


with open( out.name ) as handler:
    a = [  line.strip().split('\t') for line in handler ]

#from collections import Counter
cnt = dict()
for i in zip( *a )[0]:
    if i not in cnt: cnt[ i ] = 0
    cnt[ i ] += 1

with open( prefix + '.seglen.uniq', 'w' ) as out1:
    for b in a:
        if cnt[ b[0] ] > 1: continue
        print >>out1, '\t'.join( b )


