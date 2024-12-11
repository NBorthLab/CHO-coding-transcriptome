#!/bin/awk -f
#
# Linearize FASTA file into lines of headers and sequences separated by a tab.
{
    if ((NR>1)&&($0~/^>/)) {
        printf("\n%s", $0); 
    } else if (NR==1) {
        printf("%s", $0); 
    } else {
        printf("\t%s", $0);
    }
}

