#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path
from sys import exit


def generate_slurm_script(job_name, queue, nodes, tasks, memory, directory, prefix, suffix, paired, index,
                          strandness, add_options):
    if not Path(index + ".1.ht2").is_file():
        exit('Index file does not appear to exist. Check path and basename')

    # Check if the results and aligned_hisat directories exist, if not, create them
    if not Path('results').exists():
        Path('results').mkdir(parents=True, exist_ok=True)
        print(f"results/ directory created")
    else:
        print(f"results/ directory already exists.")
    if not Path('results/aligned_hisat').exists():
        Path('results/aligned_hisat').mkdir(parents=True, exist_ok=True)

    if directory[-1] != "/":
        directory = directory + "/"

    filebasename = f'{prefix}`printf "%03d" $SLURM_ARRAY_TASK_ID`'

    slurm_script = \
        f"""#!/bin/bash -x
#SBATCH -p {queue}
#SBATCH -N {nodes}
#SBATCH -n {tasks}
#SBATCH -J {job_name}
#SBATCH --mem={memory}
#SBATCH -o jobs/{job_name}-%A_%a.out
#SBATCH -e jobs/{job_name}-%A_%a.err

module load HISAT2/2.2.1-IGB-gcc-8.2.0-Python-3.7.2 
module load SAMtools/1.17-IGB-gcc-8.2.0

mkdir results/aligned_hisat/{filebasename}/

"""
    if paired:
        if strandness.upper() == "RF":
            strandness = "--rna-strandness RF"
        elif strandness.upper() == "FR":
            strandness = "--rna-strandness FR"
        elif strandness.upper() in ("NONE", "UNSTRANDED"):
            strandness = None
        else:
            exit('Invalid strandness. Must be RF, FR or Unstranded for paired-end sequencing')

        slurm_script = slurm_script + f'hisat2 -q --phred33 --dta -t -p 4 {strandness} -x {index} -1' \
                                      f' {directory}{filebasename}_1{suffix} -2 {directory}{filebasename}_2{suffix} ' \
                                      f'-S results/aligned_hisat/{filebasename}/{filebasename}.sam --summary-file' \
                                      f' results/aligned_hisat/{filebasename}/summary.txt --new-summary' \
                                      f' {add_options} \n\n'

    else:
        if strandness.upper() == "R":
            strandness = "--rna-strandness R"
        elif strandness.upper() == "F":
            strandness = "--rna-strandness F"
        elif strandness.upper() in ("NONE", "UNSTRANDED"):
            strandness = None
        else:
            exit('Invalid strandness. Must be R, F or Unstranded for single-end sequencing')

        slurm_script = slurm_script + f'hisat2 -q --phred33 --dta -t -p 4 {strandness} -x {index} ' \
                                      f'-U {directory}{filebasename}{suffix} ' \
                                      f'-S results/aligned_hisat/{filebasename}/{filebasename}.sam ' \
                                      f'--summary-file results/aligned_hisat/{filebasename}/summary.txt --new-summary' \
                                      f' {add_options}\n\n'

    slurm_script = slurm_script + f'samtools view -b results/aligned_hisat/{filebasename}/{filebasename}.sam ' \
                                  f'-o results/aligned_hisat/{filebasename}/{filebasename}.bam\n\n' \
                                  f'samtools sort --write-index results/aligned_hisat/{filebasename}/{filebasename}.bam' \
                                  f' -o results/aligned_hisat/{filebasename}/{filebasename}.sorted.bam\n\n' \
                                  f'rm results/aligned_hisat/{filebasename}/{filebasename}.sam\n' \
                                  f'rm results/aligned_hisat/{filebasename}/{filebasename}.bam\n\n'

    return slurm_script


def submit_job(slurm_script, job_name, array):
    if Path("./src").exists():
        srcpath = "./src/"
    else:
        srcpath = "./"

    # Write the Slurm script to a file
    with open(f"{srcpath}{job_name}.sh", "w") as file:
        file.write(slurm_script)

    # Submit the job using sbatch
    try:
        subprocess.run(['sbatch', f'--array={array}', f'{srcpath}{job_name}.sh'], check=True)
        print('sbatch ' + f'--array={array} {srcpath}{job_name}.sh')
        print("Job submitted successfully!")
    except subprocess.CalledProcessError as e:
        print("Error submitting job:", e)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Submit a HISAT2 array job to Slurm',
                                     usage='python hisat2.py Directory Prefix Index Array [Options]')
    parser.add_argument('directory', type=str, help='Path to directory with fastq files (required)')
    parser.add_argument('prefix', type=str, help='Portion of filename before the sample number (required)')
    parser.add_argument('index', type=str, help="Path and basename for HISAT2 reference genome index")
    parser.add_argument('array', type=str, help='Array of numbers to process, can be range e.g. 1-5 or list e.g. 1,2,'
                                                '5 (required)')
    hisat = parser.add_argument_group('HISAT2 options')
    hisat.add_argument('-s', '--suffix', type=str, default=".bbdtrim.fastq.gz",
                       help='Portion of filename after the sample number, and _1/_2 in the case of paired sequencing '
                            '(required) (default:.bbdtrim.fastq.gz)', metavar='')
    hisat.add_argument('-S', '--single-end', dest='paired', action='store_false', help='Toggle single-end sequencing '
                                                                                       '(default:paired-end)')
    hisat.add_argument('-r', '--strandness', type=str, default="RF", help='Strandness of sequenced reads, can be RF, FR'
                                                                          ' or Unstranded for paired-end sequencing or '
                                                                          'or F, R or Unstranded for single-end '
                                                                          'sequencing (Default:RF)', metavar='')
    hisat.add_argument('--add_options', type=str, default="", help='Pass one or more additional options to HISAT2, '
                                                                   'use quotes e.g. "--max-intronlen 10000"', metavar="")
    slurm = parser.add_argument_group('SLURM options')
    slurm.add_argument('-j', '--job_name', type=str, default='hisat2', help='Name of the job (default: bbduk)',
                       metavar='')
    slurm.add_argument('-q', '--queue', type=str, default='lowmem', help='Cluster queue to submit the job '
                                                                         '(default:lowmem)', metavar='')
    slurm.add_argument('-N', '--nodes', type=int, default=1, help='Number of nodes (default:1)', metavar='')
    slurm.add_argument('-n', '--tasks', type=int, default=4, help='Number of tasks per node (default:4)',
                       metavar='')
    slurm.add_argument('-m', '--mem', type=str, default="4GB", help='Memory requested, include units e.g. MB/GB '
                                                                    '(default:4GB)', metavar='')
    parser.set_defaults(paired=True)
    args = parser.parse_args()

    # Generate the Slurm script
    slurm_script = generate_slurm_script(args.job_name, args.queue, args.nodes, args.tasks, args.mem, args.directory,
                                         args.prefix, args.suffix, args.paired, args.index, args.strandness,
                                         args.add_options)

    # Submit the job
    submit_job(slurm_script, args.job_name, args.array)
