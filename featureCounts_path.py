#!/usr/bin/env python3

import argparse
import subprocess
from datetime import datetime
from pathlib import Path
from sys import exit

def generate_slurm_script(job_name, queue, nodes, tasks, memory, files, outname, gtf, paired, strandness):
    if not Path(gtf).is_file():
        exit('GTF file does not appear to exist. Check path')

    if outname == "getcurrenttime":
        outname = str("featCounts_" + datetime.now().strftime("%Y-%m-%d_%H.%M.%S") + ".csv")

    if paired:
        paired = " -p -C"
    else:
        paired = ""

    slurm_script = \
        f"""#!/bin/bash -x
#SBATCH -p {queue}
#SBATCH -N {nodes}
#SBATCH -n {tasks}
#SBATCH -J {job_name}
#SBATCH --mem={memory}
#SBATCH -o jobs/{job_name}-%A_%a.out
#SBATCH -e jobs/{job_name}-%A_%a.err

module load Subread/2.0.0-IGB-gcc-8.2.0 

filelist=$(ls -1 {files} |awk 'ORS=" "{{print}}')

"""
    slurm_script = slurm_script + f'featureCounts -T 4{paired} -t exon -g gene_id -s {strandness} -a {gtf} ' \
                                  f'-o {outname} $filelist\n\n'

    return slurm_script


def submit_job(slurm_script, job_name):
    if Path("./src").exists():
        srcpath = "./src/"
    else:
        srcpath = "./"

    # Write the Slurm script to a file
    with open(f"{srcpath}{job_name}.sh", "w") as file:
        file.write(slurm_script)

    # Submit the job using sbatch
    try:
        subprocess.run(['sbatch', f'{srcpath}{job_name}.sh'], check=True)
        print('sbatch ' + f'{srcpath}{job_name}.sh')
        print("Job submitted successfully!")
    except subprocess.CalledProcessError as e:
        print("Error submitting job:", e)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Submit a featureCounts job to Slurm',
                                     usage='python featureCounts.py Files GTF [Options]')
    parser.add_argument('files', type=str, help='Wildcard pattern that lists the sorted.bam files, put in quotes'
                                                'example: "results/aligned_hisat/sample_*/*.sorted.bam" (required)')
    parser.add_argument('gtf', type=str, help='Path to the genome GTF annotation the reads were aligned to (required)')

    fcounts = parser.add_argument_group("FeatureCounts options")
    fcounts.add_argument('-o', '--outname', type=str, default="getcurrenttime",
                         help='The name for the output counts table (default:featCounts_[DateTime].csv)', metavar='')
    fcounts.add_argument('-S', '--single-end', dest='paired', action='store_false',
                         help='Toggle single-end sequencing (default:paired-end)')
    fcounts.add_argument('-r', '--strandness', type=str, default="2",
                         help='Strandness of sequenced reads, can be 0 (unstranded), 1 (ISF/FR/F ligation method) or 2 '
                              '(RF/R/ISR dUTP method) (Default:2)', metavar='')
    slurm = parser.add_argument_group("SLURM options")
    slurm.add_argument('-j', '--job_name', type=str, default='featureCounts',
                       help='Name of the job (default: featureCounts)', metavar='')
    slurm.add_argument('-q', '--queue', type=str, default='lowmem',
                       help='Cluster queue to submit the job (default:lowmem)', metavar='')
    slurm.add_argument('-N', '--nodes', type=int, default=1, help='Number of nodes (default:1)', metavar='')
    slurm.add_argument('-n', '--tasks', type=int, default=4, help='Number of tasks per node (default:4)', metavar='')
    slurm.add_argument('-m', '--mem', type=str, default="3GB",
                       help='Memory requested, include units e.g. MB/GB (default:3GB)', metavar='')

    parser.set_defaults(paired=True)
    args = parser.parse_args()

    # Generate the Slurm script
    slurm_script = generate_slurm_script(args.job_name, args.queue, args.nodes, args.tasks, args.mem,
                                         args.files, args.outname, args.gtf, args.paired, args.strandness)

    # Submit the job
    submit_job(slurm_script, args.job_name)
