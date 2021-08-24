#!/usr/bin/env python


from collections import defaultdict


infile = "~/Raw_allele_data_grouped_pub.txt"
outfile = "~/Raw_allele_data_grouped_pub.gen"


openfile = open(infile, "r")
genfile = open(outfile, "w")


popdict = defaultdict(list)


def format_alleles(alleles):
    #print(alleles)
    addedallele = [(f'000{allele}') if ":" not in allele else allele for allele in alleles]
    NArplaceallele = [allele.replace("NA","000") if "NA" in allele else allele for allele in addedallele]
    replaceallele = [allele.replace(":","") if ":" in allele else allele for allele in NArplaceallele]
    #print(replaceallele)
    return replaceallele
    
         

genfile.write(f'Title line: Diploid Humulus\n')
for line in openfile:
    line = line.strip()
    #print(line)
    if line.startswith("Sample_name"):
        header = line
        #print(header)
        headerparts = header.split("\t")
        print(headerparts)
        locusnames = headerparts[3:]
        print(locusnames)
        formatlocusnames = [(f'{locus},') for i,locus in zip(range(0,len(locusnames)-1), locusnames) ]
        formatlocusnames.append(locusnames[-1])
        print(formatlocusnames)
        joinedlocusnames = " ".join(formatlocusnames)
        print(joinedlocusnames)
        genfile.write(f'{joinedlocusnames}\n')
        genfile.write(f'pop\n')
    else:
        lineparts = line.split("\t")
        #print(lineparts)
        ploidy = int(lineparts[2])
        pop = lineparts[1]
        alleledat = lineparts[3:]
        #print(alleledat)

        if ploidy <= 2:
            #print(lineparts)
            popdict[pop].append(alleledat)
            
#print(popdict)

for k,v in popdict.items():
   # print(k)
    for alleles in v:
        alleles = format_alleles(alleles)
        joinedalleles = " ".join(alleles)
        #print(alleles)
        #print(joinedalleles)
        genfile.write(f'All , {joinedalleles}\n')
                      
        
                
genfile.close()    





