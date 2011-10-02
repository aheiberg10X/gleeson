#TODOi: finish
#return filtered SNPList for each patient
def filterCoverage( SNPList ) :
    #want to run through each call of an snp.
    #if it meets certain criteria, add the patient responsible for 
    #that call the the SNPs 'shares' or 'unsure' list.
    #what this code would be doing here I don't get

    hm_cov_thresh = 8
    ht_cov_thresh = 15

    for snp in SNPList :
        for k,call in enumerate( snp.calls ) :
            if (call.genotype == 2 and call.depth >= hm_cov_thresh) or \
               (INHERITANCE_MODE == 't' and call.genotype == 3 \
                                        and call.depth >= ht_cov_thresh) :
                snp.shares.append( pat )
            elif call.depth < hm_cov_thresh or \
                 (call.genotype == 3 and call.depth < ht_cov_thresh ) :
                snp.unsure.append( pat )

    for i,pat in enumerate(patients) :
        patSNPs = []
        for j,snp in enumerate(SNPList) :
            gt = snp.calls[i].genotype
            na, hm_ref, hm_mut, ht = gt == 0, gt == 1, gt == 2, gt == 3

            depth = snp.calls[i].depth
            low_hm_cov = depth < hm_cov_thresh
            low_ht_cov = depth < ht_cov_thresh
            high_hm_cov, high_ht_cov = not low_hm_cov, not low_ht_cov

            if INHERITANCE_MODE == 'm' :
                low_conditions = hm_mut or na or (hm_ref and low_hm_cov) \
                                 or (ht and low_ht_cov)
                high_conditons = hm_mut and cov_gte

                if (low_conditions and ACCEPT_LOW_COVERAGE) or \
                   (high_conditions and not ACCEPT_LOW_COVERAGE) :
                    patSNPs.append( snp )

            elif INHERITANCE_MODE == 't' :
                low_conditions = hm_mut or na or (hm_mut and high_hm_cov) or ht
                high_conditions = (hm_mut and high_hm_cov) or \
                                  (ht and high_ht_cov)

                if (low_conditions and ACCEPT_LOW_COVERAGE) or \
                   (high_condtions and not ACCEPT_LOW_COVERAGE) :
                    patSNPs.append( snp )

        print "Patient: %s has %d SNPs" % (pat,len(patSNP))
        filterPlink( patSNP, i, plink_level )
            #print patients??! Why here?


