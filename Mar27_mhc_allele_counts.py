from epitopes import iedb
mhc_df = iedb.load_mhc(human = True, mhc_class = 1)
alleles = mhc_df['MHC Allele Name']
counts = dict(alleles.value_counts())

for gene in ('A', 'B', 'C'):
    for major in xrange(1,100):
        allele_group = "HLA-%s%d" % (gene, major)
        if allele_group not in counts:
            continue
        group_count = counts[allele_group]
        print "Allele group", allele_group, group_count

        for minor in xrange(1,100):
            allele = 'HLA-%s*%02d:%02d' % (gene, major, minor)
            if allele not in counts:
                continue
            allele_count =  counts[allele]
            print "- Allele", allele, allele_count
            group_count += allele_count
        print "Total for group:", group_count
        print