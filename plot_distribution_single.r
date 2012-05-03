
input.sd <- 'th_sd.bam.seglen' 
x <- read.delim( input.sd, header=F, row.names=1 )
y <- x[, 2]
y <- y[ y<500 ]
print( sd(y) )
#pdf( paste( input.file, '.pdf', sep='' ) )
pdf( 'segment_length_distri_merged.pdf' )
a <- hist( y, breaks=200, xlab='len/bp', main='Dist of Segment Lengths' )
#text( 400, 1e+5, sprintf( 'Peak at: %s bp', a$mids[which.max( a$counts )] ) )
dev.off()

sprintf( 'Peak at: %s bp', a$mids[which.max( a$counts )] )

