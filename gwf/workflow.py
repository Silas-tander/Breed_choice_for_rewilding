"""
------------------------------------------------------------------------------------------------------------------------
This workflow handles raw files (paired end fastqs or bams), maps reads, calles variants and genotypes variants
------------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------------
Author: Moisès Coll Macià and Juraj Bergman
Date: 08/05/2024

Do not use this pipeline without permision of the authors
------------------------------------------------------------------------------------------------------------------------
"""


#A. Importing
import os
import numpy as np
import subprocess
from gwf import Workflow, AnonymousTarget
from templates import *
import pandas as pd

gwf = Workflow()

#B. Code
## B.1. Initialize output directory
subprocess.run(["mkdir", "-p", "mcg/done"])
## B.2. Parse metadata and regions file
metadata, grouplist, groupindlist, groupref = parse_metadata()

## B.3. Map and Call
for group in grouplist:
    call_done = []
    regions, regions_list, chromosomes, chromosomes_list, batches, batches_list = parse_regions(group)
    for ind in groupindlist[group]:
        ## B.3.1. Processing raw data
        ## B.3.1.1. In the individual has bam files as raw inputs, extract the fastq files
        if len(metadata[group][ind]["bam"]):
            jobid = f"bam2fastq_{group}_{ind.replace("-", "_")}"
            gwf.target_from_template(jobid, 
                                     bam2fastq(group = group,
                                               ind   = ind, 
                                               bams  = metadata[group][ind]["bam"], 
                                               done  = f"mcg/done/{jobid}"))
            metadata[group][ind]["fastq"]["all"].append(f"mcg/{group}/samples/{ind}/fastq/b2f_{ind}_R1.fastq.gz")
            metadata[group][ind]["fastq"]["1"].append(f"mcg/{group}/samples/{ind}/fastq/b2f_{ind}_R1.fastq.gz")
            metadata[group][ind]["fastq"]["all"].append(f"mcg/{group}/samples/{ind}/fastq/b2f_{ind}_R2.fastq.gz")
            metadata[group][ind]["fastq"]["2"].append(f"mcg/{group}/samples/{ind}/fastq/b2f_{ind}_R2.fastq.gz")
        ## B.3.1.2. In the individual has multiple fastqs, concatenate them all together
        if len(metadata[group][ind]["fastq"]["all"]) > 2:
            prev_done = []
            if len(metadata[group][ind]["bam"]):
                prev_done.append(f"mcg/done/bam2fastq_{group}_{ind.replace("-", "_")}")
            jobid = f"concatfastqs_{group}_{ind.replace("-", "_")}"
            gwf.target_from_template(jobid, 
                                     concatfastqs(group     = group,
                                                  ind       = ind, 
                                                  fastqs_1  = metadata[group][ind]["fastq"]["1"],
                                                  fastqs_2  = metadata[group][ind]["fastq"]["2"],
                                                  prev_done = metadata[group][ind]["fastq"]["1"] + metadata[group][ind]["fastq"]["2"] + prev_done, 
                                                  done      = f"mcg/done/{jobid}"))
        ## B.3.1.3. Keep track of job dependencies and right fastq paths
        if len(metadata[group][ind]["fastq"]["all"]) > 2:
            prev_done = [f"mcg/done/concatfastqs_{group}_{ind.replace("-", "_")}"]
            fastq1 = f"mcg/{group}/samples/{ind}/fastq/{ind}_R1.fastq.gz"
            fastq2 = f"mcg/{group}/samples/{ind}/fastq/{ind}_R2.fastq.gz"
        elif len(metadata[group][ind]["fastq"]["all"]) == 2:
            prev_done = []
            if len(metadata[group][ind]["bam"]):
                prev_done = [f"mcg/done/bam2fastq_{group}_{ind.replace("-", "_")}"]
            fastq1 = metadata[group][ind]["fastq"]["1"][0]
            fastq2 = metadata[group][ind]["fastq"]["2"][0]
        else:
            raise Exception(f"Individual with no parired end fastqs : ind {ind} which belongs to group {group} has raw files {metadata[group][ind]["fastq"]["all"]}, which are not paired end")
        ## B.3.2. Make uBAMs from fastqs
        jobid = f"makeuBAM_{group}_{ind.replace("-", "_")}"
        gwf.target_from_template(jobid, 
                                 makeuBAM(group     = group,
                                          ind       = ind,
                                          fastq1    = fastq1, 
                                          fastq2    = fastq2, 
                                          prev_done = prev_done,
                                          done      = f"mcg/done/{jobid}"))
        prev_jobid = jobid
        ## B.3.3. split uBAMs
        jobid = f"splituBAM_{group}_{ind.replace("-", "_")}"
        gwf.target_from_template(jobid, 
                                splituBAM(group     = group,
                                          ind       = ind,
                                          prev_done = [f"mcg/done/{prev_jobid}"],
                                          done      = f"mcg/done/{jobid}"))
        splituBAM_jobid = jobid
        
        ## B.3.4. Once they are split...
        if os.path.isfile(f"mcg/done/{splituBAM_jobid}"):
            with open(f"mcg/{group}/samples/{ind}/bam/{ind}_nsplitubams.txt") as f:
                for l in f:
                    nbams = int(l.strip())
            for ishard in range(nbams):
                shard = shardstr(ishard)
                ## B.3.4.1. Mark the adapters in each bam
                jobid = f"markadapt_{group}_{ind.replace("-", "_")}_{shard}"
                gwf.target_from_template(jobid, 
                                        markadapt(group     = group,
                                                  ind       = ind,
                                                  shard     = shard,
                                                  prev_done = [f"mcg/done/{splituBAM_jobid}"],
                                                  done      = f"mcg/done/{jobid}"))
                prev_jobid = jobid

                ## B.3.4.2. Map each bam
                jobid = f"mapBAM_{group}_{ind.replace("-", "_")}_{shard}"
                gwf.target_from_template(jobid, 
                                         mapBAM(group     = group,
                                                ind       = ind,
                                                shard     = shard,
                                                ref       = metadata[group][ind]["ref_genome"],
                                                prev_done = [f"mcg/done/{prev_jobid}"],
                                                done      = f"mcg/done/{jobid}"))
                

            ## B.3.4.3. Map each bam
            jobid = f"mergeBAMs_{group}_{ind.replace("-", "_")}"
            gwf.target_from_template(jobid,
                                     mergeBAMs(group     = group,
                                               ind       = ind,
                                               bams      = [f"mcg/{group}/samples/{ind}/bam/split_uBAM/shard_{shardstr(ishard)}_markadapt_mapped.bam" for ishard in range(nbams)],
                                               prev_done = [f"mcg/done/mapBAM_{group}_{ind.replace("-", "_")}_{shardstr(ishard)}" for ishard in range(nbams)],
                                               done      = f"mcg/done/{jobid}"))
            prev_jobid = jobid
            ## B.3.4.4. Mark Duplicated in mapped bams
            jobid = f"markduplicates_{group}_{ind.replace("-", "_")}"
            gwf.target_from_template(jobid, 
                                     markduplicates(group     = group,
                                                    ind       = ind, 
                                                    prev_done = [f"mcg/done/{prev_jobid}"],
                                                    done      = f"mcg/done/{jobid}"))                                        
            prev_jobid = jobid
            ## B.3.4.5. Sort by coordinates
            jobid = f"coordsort_{group}_{ind.replace("-", "_")}"
            gwf.target_from_template(jobid, 
                                    coordsort(group     = group,
                                              ind       = ind, 
                                              prev_done = [f"mcg/done/{prev_jobid}"],
                                              done      = f"mcg/done/{jobid}"))
            coordsort_jobid = jobid
            ## B.3.4.6. Get coverage statistics
            jobid = f"cov_{group}_{ind.replace("-", "_")}"
            gwf.target_from_template(jobid,
                                     cov(group       = group,
                                         ind         = ind,
                                         regions     = regions_list,
                                         chromosomes = [regions[region]["chrom"]                     for region in regions_list],
                                         ploidies    = [regions[region][metadata[group][ind]["sex"]] for region in regions_list],
                                         starts      = [regions[region]["start"]+1                   for region in regions_list],
                                         ends        = [regions[region]["end"]                       for region in regions_list],
                                         prev_done   = [f"mcg/done/{coordsort_jobid}"],
                                         done        = f"mcg/done/{jobid}"))

            for batch in batches_list:
                jobid = f"call_{group}_{ind.replace("-", "_")}_{batch}"
                gwf.target_from_template(jobid, 
                                         call_batch(group       = group,
                                                    ind         = ind,
                                                    regions     = batches[batch][metadata[group][ind]["sex"]]["region"],
                                                    chromosomes = batches[batch][metadata[group][ind]["sex"]]["chrom"],
                                                    ploidies    = batches[batch][metadata[group][ind]["sex"]]["ploidy"],
                                                    starts      = [s+1 for s in batches[batch][metadata[group][ind]["sex"]]["start"]],
                                                    ends        = batches[batch][metadata[group][ind]["sex"]]["end"],
                                                    ref         =  metadata[group][ind]["ref_genome"],
                                                    prev_done   = [f"mcg/done/{coordsort_jobid}"],
                                                    done        = f"mcg/done/{jobid}"))
                call_done.append(f"mcg/done/{jobid}")

    ## B.4. Genotype
    ## B.4.1. Check if all individuals have been mapped and called
    if len(call_done) and not np.array([not os.path.isfile(x) for x in call_done]).sum():
        for region in regions_list:
            ## B.4.2. For every region create the folders where the outputs will be stored
            subprocess.run(f'''mkdir -p mcg/{group}/GenomicsDB/{region}/''', shell=True, text=True)
            subprocess.run(f'''mkdir -p mcg/{group}/gVCF/{region}/''',       shell=True, text=True)

            ## B.4.3. Check the individuals that should be on the genotyping according to the metadata file
            inds = [ind for ind in groupindlist[group] if regions[region][metadata[group][ind]["sex"]]]

            ## B.4.4. If there is already a cohort.sample_map file
            if os.path.isfile(f"mcg/{group}/GenomicsDB/{region}/cohort.sample_map"):
                ## B.4.4.1. Compare the individuals from the metadata file with the ones in the cohort.sample_map
                ind_in_GenomicsDB = get_ind_from_cohort(group, region)
                groupindinGenomicsDB = np.array([ind in ind_in_GenomicsDB   for ind in inds])
                GenomicsDBingroupind = np.array([ind in inds for ind in ind_in_GenomicsDB])
                ## B.4.4.2. If there are samples in the cohort.sample_map that are not in the metadata, rise an error
                if np.sum(~GenomicsDBingroupind):
                    print(f"{" ".join([ind_in_GenomicsDB[i] for i in range(len(ind_in_GenomicsDB)) if GenomicsDBingroupind[i] == 0])}")
                    raise Exception(f"GenomicsDB cohort file has individuals that are not present in the metadata file. You should check why is so and consider deleting GenomicsDB manually so that they are computed again from scratch. The individuals in GenomicsDB not in the metadata file are printed above.")
                ## B.4.4.3. If there are samples in the metadata that are not in cohort.sample_map, create a file to add those samples in the already created GenomicsDB
                elif np.sum(~groupindinGenomicsDB):
                    pd.DataFrame({"ind"  : np.array(inds)[~groupindinGenomicsDB],
                                  "path" : [f"mcg/{group}/samples/{ind}/gvcf/{ind}_{region}.gvcf.gz" for ind in np.array(inds)[~groupindinGenomicsDB]]}).to_csv(f"mcg/{group}/GenomicsDB/{region}/cohort_update.sample_map", sep='\t', header=False, index=False)
            ## B.4.5. If not, create it
            else:
                pd.DataFrame({"ind"  : inds,
                              "path" : [f"mcg/{group}/samples/{ind}/gvcf/{ind}_{region}.gvcf.gz" for ind in inds]}).to_csv(f"mcg/{group}/GenomicsDB/{region}/cohort.sample_map", sep='\t', header=False, index=False)
            
            ## B.4.6. Define vriables for the region being processed
            window = 10_000_000
            chrom  = regions[region]["chrom"]

            ## B.4.6. Divide each region in windows of length "window" to genotype
            GenomicsDB_region_done    = []
            GenotypeGVCFs_region_done = []
            GenotypeGVCFs_region_gvcf = []
            for start in range((regions[region]["start"]//window)*window, regions[region]["end"], window):
                end = start+window

                s = int(start//1e6)
                e = int(end//1e6)

                if start < regions[region]["start"]:
                    start = regions[region]["start"]
                
                if end > regions[region]["end"]:
                    end = regions[region]["end"]

                GenomicsDB_done = []
                jobid = f'GenomicsDB_{group}_{region}_{chrom}_{s}_{e}' 
                gwf.target_from_template(jobid, 
                                         GenomicsDB(group     = group,
                                                    region    = region,
                                                    chrom     = chrom,
                                                    start     = start+1,
                                                    s         = s,
                                                    end       = end,
                                                    e         = e,
                                                    prev_done = call_done,
                                                    done      = f"mcg/done/{jobid}"))

                GenomicsDB_region_done.append(f"mcg/done/{jobid}")
                GenomicsDB_done.append(f"mcg/done/{jobid}")

                if os.path.isfile(f"mcg/{group}/GenomicsDB/{region}/cohort_update.sample_map"):
                    jobid = f'GenomicsDBupdate_{group}_{region}_{chrom}_{s}_{e}' 
                    gwf.target_from_template(jobid, 
                                            GenomicsDBupdate(group     = group,
                                                             region    = region,
                                                             chrom     = chrom,
                                                             start     = start+1,
                                                             s         = s,
                                                             end       = end,
                                                             e         = e,
                                                             prev_done = GenomicsDB_region_done[0],
                                                             done      = f"mcg/done/{jobid}"))
                    GenomicsDB_region_done.append(f"mcg/done/{jobid}")
                    GenomicsDB_done.append(f"mcg/done/{jobid}")
                
                jobid = f'GenotypeGVCFs_{group}_{region}_{chrom}_{s}_{e}' 
                gwf.target_from_template(jobid, 
                                         GenotypeGVCFs(group     = group,
                                                       region    = region,
                                                       chrom     = chrom,
                                                       start     = start+1,
                                                       s         = s,
                                                       end       = end,
                                                       e         = e,
                                                       ref       = groupref[group],
                                                       prev_done = GenomicsDB_done,
                                                       done      = f"mcg/done/{jobid}"))
                GenotypeGVCFs_region_done.append(f"mcg/done/{jobid}")
                GenotypeGVCFs_region_gvcf.append(f"mcg/{group}/gVCF/{region}/{region}_{chrom}_{s}_{e}.gvcf.gz")
            
            if os.path.isfile(f"mcg/{group}/GenomicsDB/{region}/cohort_update.sample_map"):
                jobid = f'cohortupdate_{group}_{region}' 
                gwf.target_from_template(jobid, 
                                         cohortupdate(group     = group,
                                                      region    = region,
                                                      prev_done = GenomicsDB_region_done,
                                                      done      = f"mcg/done/{jobid}"))
            jobid = f'bcftoolsconcat_region_{group}_{region}'
            gwf.target_from_template(jobid, 
                                     bcftoolsconcat(vcfs      = " ".join(GenotypeGVCFs_region_gvcf),
                                                    vcf       = f"mcg/{group}/gVCF/{region}/{region}.gvcf.gz",
                                                    prev_done = GenotypeGVCFs_region_done,
                                                    done      = f"mcg/done/{jobid}"))
        for chrom in chromosomes_list:
            jobid = f'bcftoolsconcat_chromosome_{group}_{chrom}'
            gwf.target_from_template(jobid, 
                                     bcftoolsconcat(vcfs      = " ".join([f"mcg/{group}/gVCF/{region}/{region}.gvcf.gz" for region in chromosomes[chrom]]),
                                                    vcf       = f"mcg/{group}/gVCF/{chrom}.gvcf.gz",
                                                    prev_done = [f"mcg/done/bcftoolsconcat_region_{group}_{region}" for region in chromosomes[chrom]],
                                                    done      = f"mcg/done/{jobid}"))
