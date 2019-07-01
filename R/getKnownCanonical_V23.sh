# /usr/bin/bash
# Author: Komal S. Rathi
# Date: 05/17/2019
# Function: Get known canonical transcripts for Gencode v23

mysql -h genome-mysql.soe.ucsc.edu -u genome -Ne "select g.name, a.geneId, g.txEnd-g.txStart from wgEncodeGencodeBasicV23 g, wgEncodeGencodeAttrsV23 a where g.name = a.transcriptId" hg38 | sort -rnk 3 | awk '{if (!found[$2]) print ; found[$2] = 1}' | awk '{print $2}' > knownCanonicalV23.txt