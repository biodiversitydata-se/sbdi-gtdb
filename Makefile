GTDB_RELEASE = 207

ASSIGN_TAXONOMY = gunzip -c $^ | sed 's/Reversed: *//' | sed '/^>/s/>[^ ]\+ \([^[]\+\) \[.*/>\1/' | sed '/^>/s/ \[.*//' | sed 's/[a-z]__//g' | sed '/^>/s/\(Archaea\)\|\(Bacteria\)/&;&/' | gzip -c > $@
ADD_SPECIES     = gunzip -c $^ | sed 's/Reversed: *//' | sed '/^>/s/>\([^ ]\+\) .*;s__\([^[]\+\) \[.*/>\1 \2/' | gzip -c > $@

all_gtdb-sbdi.r207: all_gtdb-sbdi.207.assignTaxonomy all_gtdb-sbdi.207.addSpecies

all_gtdb-sbdi.207.assignTaxonomy: gtdb-sbdi-sativa.r207.1genome.assignTaxonomy.fna.gz gtdb-sbdi-sativa.r207.5genomes.assignTaxonomy.fna.gz gtdb-sbdi-sativa.r207.20genomes.assignTaxonomy.fna.gz

all_gtdb-sbdi.207.addSpecies: gtdb-sbdi-sativa.r207.1genome.addSpecies.fna.gz gtdb-sbdi-sativa.r207.5genomes.addSpecies.fna.gz gtdb-sbdi-sativa.r207.20genomes.addSpecies.fna.gz

rackham_download:
	for n in 1 5 20; do \
	    for t in archaea bacteria; do \
	    	scp rackham:/proj/snic2020-6-126/projects/sbdi/sativa/gtdb_16S/release$(GTDB_RELEASE)/$${t}_ssu_all.hmmer.rfmask.sativafilt$${n}pt.fna.gz $${t}_ssu_all_r$(GTDB_RELEASE).hmmer.rfmask.sativafilt$${n}pt.fna.gz; \
	    done; \
	done
	touch $@

gtdb-sbdi-sativa.r07rs207.1genome.assignTaxonomy.fna.gz: $(wildcard *.hmmer.rfmask.sativafilt1pt.fna.gz)
	$(ASSIGN_TAXONOMY)

gtdb-sbdi-sativa.r07rs207.5genomes.assignTaxonomy.fna.gz: $(wildcard *.hmmer.rfmask.sativafilt5pt.fna.gz)
	$(ASSIGN_TAXONOMY)

gtdb-sbdi-sativa.r07rs207.20genomes.assignTaxonomy.fna.gz: $(wildcard *.hmmer.rfmask.sativafilt20pt.fna.gz)
	$(ASSIGN_TAXONOMY)

gtdb-sbdi-sativa.r07rs207.1genome.addSpecies.fna.gz: $(wildcard *.hmmer.rfmask.sativafilt1pt.fna.gz)
	$(ADD_SPECIES)

gtdb-sbdi-sativa.r07rs207.5genomes.addSpecies.fna.gz: $(wildcard *.hmmer.rfmask.sativafilt5pt.fna.gz)
	$(ADD_SPECIES)

gtdb-sbdi-sativa.r07rs207.20genomes.addSpecies.fna.gz: $(wildcard *.hmmer.rfmask.sativafilt20pt.fna.gz)
	$(ADD_SPECIES)

MANIFEST.txt: README.txt $(wildcard gtdb-sbdi-sativa.r07rs207.*.fna.gz)
	for f in $^; do \
	    echo "$$f ($$(LC_ALL=C ls -lhL $$f | sed 's/.* \([0-9.]\+\)\([KM]\) .*/\1 \2iB/'))"; \
	done > $@
