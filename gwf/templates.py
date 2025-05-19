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
from gwf import Workflow, AnonymousTarget
gwf = Workflow()

#B. Functions
def parse_metadata():
    metadata     = {}
    grouplist    = []
    groupindlist = {}
    groupref     = {}
    with open("metadata.txt") as f:
        for i, l in enumerate(f):
            if i:
                ind, sex, group, ref_genome, raw_file = l.strip().split("\t")
                if sex not in ["M", "F"]:
                    sex = "F"
                if group not in grouplist:
                    metadata[group]     = {}
                    groupindlist[group] = []
                    grouplist.append(group)
                if ind not in groupindlist[group]:
                    groupindlist[group].append(ind)
                    metadata[group][ind] = {}
                    metadata[group][ind]["sex"]          = sex
                    metadata[group][ind]["group"]        = group
                    metadata[group][ind]["ref_genome"]   = ref_genome
                    metadata[group][ind]["bam"]          = []
                    metadata[group][ind]["fastq"]        = {}
                    metadata[group][ind]["fastq"]["1"]   = []
                    metadata[group][ind]["fastq"]["2"]   = []
                    metadata[group][ind]["fastq"]["all"] = []
                if ".bam" in raw_file:
                    metadata[group][ind]["bam"].append(raw_file)
                elif ".fastq" in raw_file:
                    metadata[group][ind]["fastq"]["all"].append(raw_file)
                    if   "1.fastq" in raw_file:
                        metadata[group][ind]["fastq"]["1"].append(raw_file)
                    elif "2.fastq" in raw_file:
                        metadata[group][ind]["fastq"]["2"].append(raw_file)
                    else:
                        raise Exception(f"Fastq file not recognized : ind {ind} which belongs to group {group} has raw file {raw_file}, which does not have the pattern '1.fastq' nor '2.fastq'")
                else:
                    raise Exception(f"Raw file not recognized : ind {ind} which belongs to group {group} has raw file {raw_file}")
                if group not in groupref:
                    groupref[group] = ref_genome
                else:
                    if groupref[group] != ref_genome:
                        raise Exception(f"Different reference genome : ind {ind} which belongs to group {group} has a different reference genome {ref_genome} than the others with reference {groupref[group]}")

    return metadata, grouplist, groupindlist, groupref

def parse_regions(group):
    regions           = {}
    regions_list      = []
    chromosomes       = {}
    chromosomes_list  = []
    batches      = {}
    batches_list           = []
    with open(f"regions_{group}.txt") as f:
        for i, l in enumerate(f):
            if i:
                region, chrom, start, end, batch, female_ploidy, male_ploidy = l.strip().split("\t")
                female_ploidy = int(female_ploidy)
                male_ploidy   = int(male_ploidy)
                start         = int(start)
                end           = int(end)


                regions_list.append(region)
                regions[region] = {}
                regions[region]["chrom"] = chrom
                regions[region]["start"] = start
                regions[region]["end"]   = end
                regions[region]["F"]     = female_ploidy
                regions[region]["M"]     = male_ploidy

                if chrom not in chromosomes:
                    chromosomes[chrom] = []
                chromosomes[chrom].append(region)
                
                if chrom not in chromosomes_list:
                    chromosomes_list.append(chrom)

                if batch not in batches:
                    batches[batch] = {}
                    batches_list.append(batch)
                    for sex, ploidy in zip(["F", "M"], [female_ploidy, male_ploidy]):
                        batches[batch][sex] = {}
                        batches[batch][sex]["region"] = []
                        batches[batch][sex]["chrom"]  = []
                        batches[batch][sex]["start"]  = []
                        batches[batch][sex]["end"]    = []
                        batches[batch][sex]["ploidy"] = []
                for sex, ploidy in zip(["F", "M"], [female_ploidy, male_ploidy]):
                    if ploidy:
                        batches[batch][sex]["region"].append(region)
                        batches[batch][sex]["chrom"].append(chrom)
                        batches[batch][sex]["start"].append(start)
                        batches[batch][sex]["end"].append(end)
                        batches[batch][sex]["ploidy"].append(ploidy)

                
    return regions, regions_list, chromosomes, chromosomes_list, batches, batches_list

def shardstr(ishard):
    return (4-len(str(ishard+1)))*"0"+str(ishard+1)

def get_ind_from_cohort(group, region):
    ind_in_GenomicsDB = []
    with open(f"mcg/{group}/GenomicsDB/{region}/cohort.sample_map") as f:
        for l in f:
            ind, _ = l.strip().split()
            ind_in_GenomicsDB.append(ind)
    return ind_in_GenomicsDB

#C. Templates
job_header = '''
    echo "JOBID:" $PBS_JOBID
    echo "NODE :" $HOSTNAME
    echo "USER :" $USER
    source ~/.bashrc
    conda activate gwf
    echo "CONDA:" $CONDA_DEFAULT_ENV
'''

default_options = {"cores" : 1, 'memory': "8g", 'walltime': "1-00:00:00", 'account': "Coregonus"}

def all_done(jobid, inputs, output):
    '''
    A dummy job to just flag that the same job for different inputs have been done
    '''
    outputs = [output]
    options = default_options.copy()
    options["walltime"] = "00:10:00"
    spec    = job_header+f'''	
    touch {output}'''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def bam2fastq(group, ind, bams, done):
    '''
    Convert bams to fastqs
    '''
    inputs  = bams
    outputs = [done]
    options = default_options.copy()
    spec    = job_header+f'''
    out_dir=mcg/{group}/samples/{ind}/fastq

    mkdir -p ${{out_dir}}

    samtools cat {" ".join(bams)} \
        | samtools sort  -n                        \
                         -T /scratch/$SLURM_JOB_ID \
                         --threads 1               \
        | samtools fastq -1 ${{out_dir}}/b2f_{ind}_R1.fastq.gz \
                         -2 ${{out_dir}}/b2f_{ind}_R2.fastq.gz \
                         -s ${{out_dir}}/b2f_{ind}_s.fastq.gz  
    touch {done}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def concatfastqs(group, ind, fastqs_1, fastqs_2, prev_done, done):
    '''
    Concatenate fastqs
    '''
    inputs    = prev_done
    outputs   = [done]
    options   = default_options.copy()
    spec      = job_header+f'''
    out_dir=mcg/{group}/samples/{ind}/fastq

    mkdir -p ${{out_dir}}

    cat {" ".join(fastqs_1)} > ${{out_dir}}/{ind}_R1.fastq.gz
    cat {" ".join(fastqs_2)} > ${{out_dir}}/{ind}_R2.fastq.gz

    touch {done}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def makeuBAM(group, ind, fastq1, fastq2, prev_done, done):
    '''
    Make uBAM file
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec    = job_header+f"""
    out_dir=mcg/{group}/samples/{ind}/bam

    mkdir -p ${{out_dir}}

    picard FastqToSam --FASTQ  {fastq1}                         \
                      --FASTQ2 {fastq2}                         \
                      --OUTPUT mcg/{group}/samples/{ind}/bam/u{ind}.bam \
                      --SAMPLE_NAME {ind}                       \
                      --TMP_DIR /scratch/$SLURM_JOB_ID
    

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def splituBAM(group, ind, prev_done, done):
    '''
    Split uBAMs
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec    = job_header+f"""
    splitubamdir=mcg/{group}/samples/{ind}/bam/split_uBAM
    rm    -fr ${{splitubamdir}}
    mkdir -p  ${{splitubamdir}} 

    picard SplitSamByNumberOfReads -I mcg/{group}/samples/{ind}/bam/u{ind}.bam \
                                   -O ${{splitubamdir}}                \
                                   --CREATE_INDEX true                 \
                                   -N_READS 48000000

    ls ${{splitubamdir}} | wc -l > mcg/{group}/samples/{ind}/bam/{ind}_nsplitubams.txt

    rm -f mcg/{group}/samples/{ind}/fastq/*
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def markadapt(group, ind, shard, prev_done, done):
    '''
    Mark adaptor sequences in uBAM file
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    picard MarkIlluminaAdapters -I mcg/{group}/samples/{ind}/bam/split_uBAM/shard_{shard}.bam           \
                                -O mcg/{group}/samples/{ind}/bam/split_uBAM/shard_{shard}_markadapt.bam \
                                -M mcg/{group}/samples/{ind}/bam/split_uBAM/shard_{shard}_markadapt.txt \
                                --TMP_DIR /scratch/$SLURM_JOB_ID

    rm -f mcg/{group}/samples/{ind}/bam/u{ind}.bam
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mapBAM(group, ind, shard, ref, prev_done, done):
    '''
    Mapping bam file
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    picard SamToFastq -I           mcg/{group}/samples/{ind}/bam/split_uBAM/shard_{shard}_markadapt.bam       \
                      --FASTQ      /dev/stdout                                                        \
                      --INTERLEAVE true                                                               \
                      --TMP_DIR    /scratch/$SLURM_JOB_ID                                             \
        | bwa mem -M                    \
                  -t {options["cores"]} \
                  -p {ref}              \
                  /dev/stdin            \
        | picard MergeBamAlignment --ALIGNED_BAM                  /dev/stdin                                                          \
                                   --UNMAPPED_BAM                 mcg/{group}/samples/{ind}/bam/split_uBAM/shard_{shard}_markadapt.bam        \
                                   --OUTPUT                       mcg/{group}/samples/{ind}/bam/split_uBAM/shard_{shard}_markadapt_mapped.bam \
                                   -R                             {ref}                                                               \
                                   --CREATE_INDEX                 true                                                                \
                                   --INCLUDE_SECONDARY_ALIGNMENTS false                                                               \
                                   --TMP_DIR                      /scratch/$SLURM_JOB_ID

    rm -f mcg/{group}/samples/{ind}/bam/split_uBAM/shard_{shard}.bam
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mergeBAMs(group, ind, bams, prev_done, done):
    '''
    Merge bam files
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    picard MergeSamFiles {" ".join([f"-I {bam} " for bam in bams])}                  \
                         -O mcg/{group}/samples/{ind}/bam/{ind}_markadapt_mapped_merged.bam  \
                         --SORT_ORDER   queryname                                    \
                         --TMP_DIR      /scratch/$SLURM_JOB_ID

    rm -f mcg/{group}/samples/{ind}/bam/split_uBAM/shard_*_markadapt.bam
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def markduplicates(group, ind, prev_done, done):
    '''
    Mark and remove duplicates
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    picard MarkDuplicates -I mcg/{group}/samples/{ind}/bam/{ind}_markadapt_mapped_merged.bam                \
                          -M mcg/{group}/samples/{ind}/bam/{ind}_markadapt_mapped_merged_markduplicates.txt \
                          -O mcg/{group}/samples/{ind}/bam/{ind}_markadapt_mapped_merged_markduplicates.bam \
                          --REMOVE_DUPLICATES true                                                  \
                          --CREATE_INDEX true                                                       \
                          --TMP_DIR /scratch/$SLURM_JOB_ID

    rm -fr mcg/{group}/samples/{ind}/bam/split_uBAM
    rm -f mcg/{group}/samples/{ind}/bam/{ind}_markadapt_mapped_merged_markduplicates.txt
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def coordsort(group, ind, prev_done, done):
    '''
    Sort BAM by coordinates
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    picard SortSam -I mcg/{group}/samples/{ind}/bam/{ind}_markadapt_mapped_merged_markduplicates.bam           \
                   -O mcg/{group}/samples/{ind}/bam/{ind}_markadapt_mapped_merged_markduplicates_coordsort.bam \
                   -SO coordinate                                                                      \
                   --CREATE_INDEX true                                                                 \
                   --TMP_DIR /scratch/$SLURM_JOB_ID

    rm -f mcg/{group}/samples/{ind}/bam/{ind}_markadapt_mapped_merged.bam
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def cov(group, ind, regions, chromosomes, ploidies, starts, ends, prev_done, done):
    '''
    Get average coverage across chromosomes at covered sites, for files already coordinate-sorted
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    covdir=mcg/{group}/samples/{ind}/cov
    rm    -fr ${{covdir}}
    mkdir -p  ${{covdir}} 

    regions=({" ".join(regions)})
    chromosomes=({" ".join(chromosomes)})
    ploidies=({" ".join([str(p) for p in ploidies])})
    starts=({" ".join([str(s) for s in starts])})
    ends=({" ".join([str(e) for e in ends])})
    length=${{#regions[@]}}

    # Iterate over both arrays simultaneously
    for ((i = 0; i < length; i++)); do
        region=${{regions[i]}}
        chrom=${{chromosomes[i]}}
        ploidy=${{ploidies[i]}}
        start=${{starts[i]}}
        end=${{ends[i]}}


        if [ ${{ploidy}} -gt 0 ];
        then
            echo ${{region}} ${{chrom}}
            date

            samtools depth -r ${{chrom}}:${{start}}-${{end}} mcg/{group}/samples/{ind}/bam/{ind}_markadapt_mapped_merged_markduplicates_coordsort.bam \
                | awk '{{sum += $3}} END {{if(sum == 0 || NR == 0){{cov=0}}else{{cov=sum/NR}};print "'${{region}}'\t'${{chrom}}'\t'${{start}}'\t'${{end}}'\t"NR"\t"sum"\t"cov}}' >> ${{covdir}}/{ind}.cov

        fi

    done

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def call_batch(group, ind, regions, chromosomes, ploidies, starts, ends, ref, prev_done, done):
    '''
    Call a batch of chromosomes per individual
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    dir=mcg/{group}/samples/{ind}
    mkdir -p  ${{dir}}/gvcf

    regions=({" ".join(regions)})
    chromosomes=({" ".join(chromosomes)})
    ploidies=({" ".join([str(p) for p in ploidies])})
    starts=({" ".join([str(s) for s in starts])})
    ends=({" ".join([str(e) for e in ends])})
    length=${{#regions[@]}}

    # Iterate over both arrays simultaneously
    for ((i = 0; i < length; i++)); do
        region=${{regions[i]}}
        chrom=${{chromosomes[i]}}
        ploidy=${{ploidies[i]}}
        start=${{starts[i]}}
        end=${{ends[i]}}

        echo "   1. Running" {ind} ${{region}}...
        echo "#############################"
        date

        if [ ! -f {done}_${{region}} ] && [ ${{ploidy}} -gt 0 ];
        then
            echo "   2. Getting the mapped reads mapped to the chromosome..."
            date

            samtools view -b ${{dir}}/bam/{ind}_markadapt_mapped_merged_markduplicates_coordsort.bam ${{chrom}}:${{start}}-${{end}} > /scratch/$SLURM_JOB_ID/{ind}_${{region}}.bam

            echo "   3. Add or Replace read groups..."
            date
            gatk AddOrReplaceReadGroups -I             /scratch/$SLURM_JOB_ID/{ind}_${{region}}.bam           \
                                        -O             /scratch/$SLURM_JOB_ID/{ind}_${{region}}_readgroup.bam \
                                        -LB            lib1                                                   \
                                        -PL            ILLUMINA                                               \
                                        -PU            unit1                                                  \
                                        -SM            {ind}                                                  \
                                        --CREATE_INDEX true                                                   \
                                        --TMP_DIR      /scratch/$SLURM_JOB_ID

            echo "   4. Calling..."
            date
            gatk HaplotypeCaller -R                        {ref}                                                  \
                                 -I                        /scratch/$SLURM_JOB_ID/{ind}_${{region}}_readgroup.bam \
                                 -L                        ${{chrom}}:${{start}}-${{end}}                          \
                                 -ploidy                   ${{ploidy}}                                            \
                                 --native-pair-hmm-threads {options["cores"]}                                     \
                                 -ERC                      BP_RESOLUTION                                          \
                                 -O                        ${{dir}}/gvcf/{ind}_${{region}}.gvcf.gz 

            echo "   5. Indexing..."
            date
            gatk IndexFeatureFile -I ${{dir}}/gvcf/{ind}_${{region}}.gvcf.gz

            echo "   6. DONE!"
            date

            touch {done}_${{region}}

        else
            if [ -f {done}_${{region}} ];
            then
                echo "   "2. {done}_${{region}}" is created. Skiping this chromosome..."
            else
                echo "   "2. ploidy is ${{ploidy}}". Skipping this chromosome..."
            fi
        fi

    done

    rm -f mcg/{group}/samples/{ind}/bam/{ind}_markadapt_mapped_merged_markduplicates.bam

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def GenomicsDB(group, region, chrom, start, s, end, e, prev_done, done):
    '''
    Create a GenomicsDB
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    gatk GenomicsDBImport --sample-name-map                         mcg/{group}/GenomicsDB/{region}/cohort.sample_map \
                          --genomicsdb-workspace-path               /scratch/$SLURM_JOB_ID/{region}_{chrom}_{s}_{e}   \
                          --tmp-dir                                 /scratch/$SLURM_JOB_ID      \
                          --genomicsdb-shared-posixfs-optimizations true                        \
                          --genomicsdb-vcf-buffer-size              4194304                     \
                          --intervals                               {chrom}:{start}-{end}       \
                          --reader-threads                          {options["cores"]}

    
    chmod -R 777 /scratch/$SLURM_JOB_ID/{region}_{chrom}_{s}_{e}

    mv /scratch/$SLURM_JOB_ID/{region}_{chrom}_{s}_{e} mcg/{group}/GenomicsDB/{region}/{region}_{chrom}_{s}_{e}
    
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def GenomicsDBupdate(group, region, chrom, start, s, end, e, prev_done, done):
    '''
    Genotype
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    cp -r mcg/{group}/GenomicsDB/{region}/{region}_{chrom}_{s}_{e} /scratch/$SLURM_JOB_ID/{region}_{chrom}_{s}_{e}
    
    gatk GenomicsDBImport --sample-name-map                         mcg/{group}/GenomicsDB/{region}/cohort_update.sample_map \
                          --genomicsdb-update-workspace-path        /scratch/$SLURM_JOB_ID/{region}_{chrom}_{s}_{e}          \
                          --tmp-dir                                 /scratch/$SLURM_JOB_ID      \
                          --genomicsdb-shared-posixfs-optimizations true                        \
                          --genomicsdb-vcf-buffer-size              4194304                     \
                          --intervals                               {chrom}:{start}-{end}       \
                          --reader-threads                          {options["cores"]}

    chmod -R 777 /scratch/$SLURM_JOB_ID/{region}_{chrom}_{s}_{e}
    mv mcg/{group}/GenomicsDB/{region}/{region}_{chrom}_{s}_{e} mcg/{group}/GenomicsDB/{region}/{region}_{chrom}_{s}_{e}_tmp
    mv /scratch/$SLURM_JOB_ID/{region}_{chrom}_{s}_{e} mcg/{group}/GenomicsDB/{region}/{region}_{chrom}_{s}_{e}
    rm -r mcg/{group}/GenomicsDB/{region}/{region}_{chrom}_{s}_{e}_tmp
    
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def GenotypeGVCFs(group, region, chrom, start, s, end, e, ref, prev_done, done):
    '''
    Genotype
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    cp -r mcg/{group}/GenomicsDB/{region}/{region}_{chrom}_{s}_{e} /scratch/$SLURM_JOB_ID/{region}_{chrom}_{s}_{e}
    
    gatk GenotypeGVCFs --include-non-variant-sites                                                                \
                       -R                                 {ref}                                                   \
                       --variant                          gendb:///scratch/$SLURM_JOB_ID/{region}_{chrom}_{s}_{e} \
                       -O                                 /scratch/$SLURM_JOB_ID/{region}_{chrom}_{s}_{e}.gvcf.gz \
                       -L                                 {chrom}:{start}-{end}                                   \
                       --max-genotype-count               3000                                                    \
                       --genomicsdb-max-alternate-alleles 3001                                                    \
                       --max-alternate-alleles            200

    chmod -R 777 /scratch/$SLURM_JOB_ID/{region}_{chrom}_{s}_{e}.gvcf.gz
    mv /scratch/$SLURM_JOB_ID/{region}_{chrom}_{s}_{e}.gvcf.gz mcg/{group}/gVCF/{region}/{region}_{chrom}_{s}_{e}.gvcf.gz
    
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def cohortupdate(group, region, prev_done, done):
    '''
    update cohort files
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""

    cat mcg/{group}/GenomicsDB/{region}/cohort.sample_map mcg/{group}/GenomicsDB/{region}/cohort_update.sample_map > mcg/{group}/GenomicsDB/{region}/cohort.sample_map.tmp
    mv mcg/{group}/GenomicsDB/{region}/cohort.sample_map.tmp mcg/{group}/GenomicsDB/{region}/cohort.sample_map
    
    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def bcftoolsconcat(vcfs, vcf, prev_done, done):
    '''
    concatenate vcfs
    '''
    inputs  = prev_done
    outputs = [done]
    options = default_options.copy()
    spec = job_header+f"""
    bcftools concat {vcfs} -O z -o {vcf}

    touch {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
