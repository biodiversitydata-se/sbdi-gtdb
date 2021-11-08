all_gtdb-sbdi-sativa: gtdb-sbdi-sativa.assignTaxonomy.fna.gz gtdb-sbdi-sativa.addSpecies.fna.gz

gtdb-sbdi-sativa.assignTaxonomy.fna.gz: $(wildcard *.hmmer.rfmask.sativafilt5pt.fna.gz)
	gunzip -c $^ | sed 's/Reversed: *//' | sed '/^>/s/>\([^ ]\+\) \([^[]\+\) \[.*/>\2(\1\)/' | sed '/^>/s/;s__.*//' | sed 's/[a-z]__//g' | sed 's/ /_/g' | sed '/^>/s/\(Archaea\)\|\(Bacteria\)/&;&/' | gzip -c > $@

gtdb-sbdi-sativa.addSpecies.fna.gz: $(wildcard *.hmmer.rfmask.sativafilt5pt.fna.gz)
	gunzip -c $^ | sed 's/Reversed: *//' | sed '/^>/s/>\([^ ]\+\) .*;s__\([^[]\+\) \[.*/>\1 \2/' | gzip -c > $@

gtdb-sbdi-sativa.r06rs202.fna.gz: $(wildcard *.hmmer.rfmask.sativafilt5pt.fna.gz)
	cat $^ > $@
