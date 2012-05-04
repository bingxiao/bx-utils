#!/usr/bin/env python
''' basic statistics to paired-end mapping '''

import numpy as np
import pysam
import sys
import os

def usage():
    print >>sys.stdout, 'Usage: %s <bam_file.bam> <output_directory>'%sys.argv[0]

try:
    input_file = sys.argv[1]
    prefix = os.path.basename( input_file )
    if not os.path.exists( input_file + '.bai' ):
        print >>sys.stderr, '.bai file missing, you have to index the bam file first:\n# samtools index %s'%input_file
        sys.exit()
    #samfile = pysam.Samfile( 'a.bam', 'rb' )
    samfile = pysam.Samfile( input_file, 'rb' )
except:
    usage()
    exit(0)

out_dir = sys.argv[2]
if not os.path.exists( out_dir ):
    os.mkdir( out_dir )

def get_seg_len( read1, read2 ):
    return max( read1.aend, read2.aend ) - min( read1.pos, read2.pos )

read1_dic = {}
read2_dic = {}

segment_length_lst = []
with open( os.path.join( out_dir, prefix+'.seglen.all'), 'w' ) as out:
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
            print >>sys.stderr, '\rParsing %d reads.'%n,

print >>sys.stderr, "\nOutput file:", out.name, "writen."
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
print >>sys.stderr, "Preparing unique file..."

cnt = dict()
for i in zip( *a )[0]:
    if i not in cnt: 
        cnt[ i ] = 1
    else:
        cnt[ i ] += 1

print >>sys.stderr, "Output file:", out.name, "writing.."
with open( os.path.join( out_dir, prefix+'.seglen.uniq' ), 'w' ) as out1:
    for b in a:
        if cnt[ b[0] ] > 1: continue
        print >>out1, '\t'.join( b )
print >>sys.stderr, "Done."


