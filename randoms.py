def bgs_mask_randoms(random):    
    
    

    # Apply custom imaging mask around bright stars etc.   
    bitnamelist = ["BRIGHT", "CLUSTER"] 
    bits = get_imaging_maskbits(bitnamelist) 
    retain_random = np.ones(len(random['MASKBITS']), dtype=bool) #got rid of .data as didnt work below 

    for bit, ttype in zip(bits, bitnamelist): 
        # Keep random if bit not set for bits corresponding to BRIGHT and CLUSTER. 
        retain_random &= ((random['MASKBITS'] & 2**bit) == 0)  

        #print to show amount of randoms cut 
        print(ttype, bit, np.mean(retain_random))

    #other cuts
    NOBS_mask = ((random['NOBS_G'] > 0) | (random['NOBS_R'] > 0) | (random['NOBS_Z'] > 0))
    
    #final mask with all masks incorporated 
    retain_random = retain_random & NOBS_mask
 
    #again print to show updated cuts 
    print('NOBS', np.mean(retain_random))

    return random[retain_random]
