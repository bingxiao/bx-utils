#!/usr/bin/env python

"""
BX Note: 
    this script was modified/bug-fixed based on GALAXY: galaxy-central/tools/filters/ucsc_gene_bed_to_exon_bed.py
    http://bitbucket.org/galaxy/galaxy-central/ 

Read a table dump in the UCSC gene table format and print a tab separated
list of intervals corresponding to requested features of each gene.

usage: ucsc_gene_table_to_intervals.py [options]

options:
  -h, --help                  show this help message and exit
  -rREGION, --region=REGION
                              Limit to region: one of coding, utr3, utr5, codon, intron, transcribed [default]
  -e, --exons                 Only print intervals overlapping an exon
  -i, --input=inputfile       input file
  -o, --output=outputfile     output file
"""

import optparse, string, sys

assert sys.version_info[:2] >= ( 2, 4 )

def main():

    # Parse command line    
    parser = optparse.OptionParser( usage="%prog [options] " )
    parser.add_option( "-r", "--region", dest="region", default="transcribed",
                       help="Limit to region: one of coding, utr3, utr5, transcribed [default]" )
    parser.add_option( "-e", "--exons",  action="store_true", dest="exons", default=True,
                       help="Only print intervals overlapping an exon" )
#    parser.add_option( "-s", "--strand",  action="store_true", dest="strand",
#                       help="Print strand after interval" )
    parser.add_option( "-i", "--input",  dest="input",  default=None,
                       help="Input file, or - for pipe in" )
    parser.add_option( "-o", "--output", dest="output", default=None,
                       help="Output file, or - for pipe out" )
    options, args = parser.parse_args()
    assert options.region in ( 'coding', 'utr3', 'utr5', 'transcribed', 'intron', 'codon' ), "Invalid region argument"

    print >>sys.stderr, sys.argv[1]+":", options
    assert options.output != options.input, "Plz correct the input&output file"

    try:
        if options.output in [ "-", ""]:
            out_file = sys.stdout
        else:
            out_file = open (options.output, "w")
    except:
        print >> sys.stderr, "Bad output file."
        sys.exit(0)
    
    try:
        if options.input == "-":
            in_file = sys.stdin
        else:
            in_file = open (options.input)
    except:
        print >> sys.stderr, "Bad input file."
        sys.exit(0)
    
    print >>sys.stderr, "Region:", options.region+";"
    """print "Only overlap with Exons:",
    if options.exons:
        print "Yes"
    else:
        print "No"
    """
    
    # Read table and handle each gene
    for line in in_file:
        #try:
            if line[0:1] == "#":
                continue
            # Parse fields from gene tabls
            fields = line.split( '\t' )
            chrom     = fields[0]
            tx_start  = int( fields[1] )
            tx_end    = int( fields[2] )
            name      = fields[3]
            strand    = fields[5].replace(" ","_")
            cds_start = int( fields[6] )
            cds_end   = int( fields[7] )

            # Determine the subset of the transcribed region we are interested in
            if options.region == 'utr3':
                if strand == '-': region_start, region_end = tx_start, cds_start
                else: region_start, region_end = cds_end, tx_end 
            elif options.region == 'utr5':
                if strand == '-': region_start, region_end = cds_end, tx_end
                else: region_start, region_end = tx_start, cds_start
            elif options.region == 'coding' or options.region == 'codon':
                region_start, region_end = cds_start, cds_end
            else:
                region_start, region_end = tx_start, tx_end

            # If only interested in exons, print the portion of each exon overlapping
            # the region of interest, otherwise print the span of the region
        # options.exons is always TRUE
            if options.exons:
                exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
                exon_starts = map((lambda x: x + tx_start ), exon_starts)
                exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
                exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);

        #for Intron regions:
            if options.region == 'intron':
                i=0
                while i < len(exon_starts)-1:
                    intron_starts = exon_ends[i]
                    intron_ends = exon_starts[i+1]
                    if strand: print_tab_sep(out_file, chrom, intron_starts, intron_ends, name, "0", strand )
                    else: print_tab_sep(out_file, chrom, intron_starts, intron_ends )
                    i+=1
        #for non-intron regions:
            else:
                for start, end in zip( exon_starts, exon_ends ):
                    start = max( start, region_start )
                    end = min( end, region_end )
                    if start < end:
                        if options.region == 'codon':
                            start += (3 - ((start-region_start)%3))%3
                            c_start = start 
                            while c_start+3 <= end:
                                if strand:
                                    print_tab_sep(out_file, chrom, c_start, c_start+3, name, "0", strand )
                                else:
                                    print_tab_sep(out_file, chrom, c_start, c_start+3)
                                c_start += 3
                        else:
                            if strand:
                                print_tab_sep(out_file, chrom, start, end, name, "0", strand )
                            else: 
                                print_tab_sep(out_file, chrom, start, end )
                    """
                    else:
                        if options.region == 'codon':
                            c_start = start
                            c_end = end
                            if c_start > c_end:
                                t = c_start
                                c_start = c_end
                                c_end = t
                            while c_start+3 <= c_end:
                                if strand:
                                    print_tab_sep(out_file, chrom, c_start, c_start+3, name, "0", strand )
                                else:
                                    print_tab_sep(out_file, chrom, c_start, c_start+3)
                                c_start += 3
                        else:
                            if strand:
                                print_tab_sep(out_file, chrom, region_start, region_end, name, "0", strand )
                            else: 
                                print_tab_sep(out_file, chrom, region_start, region_end )
                    """
#         except:
#             print >>sys.stderr, 
#             continue

def print_tab_sep(out_file, *args ):
    """Print items in `l` to stdout separated by tabs"""
    print >>out_file, string.join( [ str( f ) for f in args ], '\t' )

if __name__ == "__main__": main()
