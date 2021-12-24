###
### Read a fasta file and its blast-output (versus Silva132)
### and report somes statistics
###
### Mestivier Dens - December 2021
###

import sys

###
### command line
###

if len( sys.argv )!=3:
    print( "Syntax: %s file.fa outblast.csv" % sys.argv[0])
    sys.exit(1)

###
### read fasta
### Construct a dico of seqname
###

d = {}

