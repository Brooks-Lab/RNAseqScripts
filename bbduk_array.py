#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path


def generate_slurm_script(job_name, queue, nodes, tasks, memory, directory, prefix, ext, paired):
    if directory[-1] != "/":
        directory = directory + "/"

    filebasename = f'{directory}{prefix}`printf "%03d" $SLURM_ARRAY_TASK_ID`'

    slurm_script = \
        f"""#!/bin/bash -x
#SBATCH -p {queue}
#SBATCH -N {nodes}
#SBATCH -n {tasks}
#SBATCH -J {job_name}
#SBATCH --mem={memory}
#SBATCH -o jobs/{job_name}-%A_%a.out
#SBATCH -e jobs/{job_name}-%A_%a.err

module load BBMap/38.94-Java-1.8.0_201

"""
    if paired:
        slurm_script = slurm_script + f'bbduk.sh in1={filebasename}_#{ext} out={filebasename}_#.bbdtrim{ext} ' \
                                      f'ref=adapters k=23 ktrim=r mink=11 qtrim=r trimq=14 minlength=20 hdist=1\n\n'
    else:
        slurm_script = slurm_script + f'bbduk.sh in={filebasename}{ext} out={filebasename}.bbdtrim{ext} ref=adapters ' \
                                    f'k=23 ktrim=r mink=11 qtrim=r trimq=14 minlength=20 hdist=1'

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
    parser = argparse.ArgumentParser(description='Submit an array job to Slurm', usage='python bbduk.py Directory '
                                                                                       'Prefix Array [Options] '
                                                                                       '')
    parser.add_argument('directory', type=str, help='Path to directory with fastq files (required)')
    parser.add_argument('prefix', type=str, help='Portion of filename before the sample number (required)')
    parser.add_argument('array', type=str, help='Array of numbers to process, can be range e.g. 1-5 or list e.g. 1,2,'
                                                '5 (required)')
    bbduk = parser.add_argument_group("BBDuk options")
    bbduk.add_argument('-e', '--ext', type=str, default=".fastq.gz", help='The fastq file extension '
                                                                          '(default:.fastq.gz)', metavar='')
    bbduk.add_argument('-S', '--single-end', dest='paired', action='store_false', help='Toggle single-end sequencing '
                                                                                       '(default:paired-end)')
    slurm = parser.add_argument_group("SLURM options")
    slurm.add_argument('-j', '--job_name', type=str, default='bbduk', help='Name of the job (default: bbduk)',
                       metavar='')
    slurm.add_argument('-q', '--queue', type=str, default='lowmem', help='Cluster queue to submit the job '
                                                                         '(default:lowmem)', metavar='')
    slurm.add_argument('-N', '--nodes', type=int, default=1, help='Number of nodes (default:1)', metavar='')
    slurm.add_argument('-n', '--tasks', type=int, default=1, help='Number of tasks per node (default:1)',
                       metavar='')
    slurm.add_argument('-m', '--mem', type=str, default="1GB", help='Memory requested, include units e.g. MB/GB '
                                                                    '(default:1GB)', metavar='')

    parser.set_defaults(paired=True)
    args = parser.parse_args()

    # Generate the Slurm script
    slurm_script = generate_slurm_script(args.job_name, args.queue, args.nodes, args.tasks, args.mem,
                                         args.directory, args.prefix, args.ext, args.paired)

    # Submit the job
    submit_job(slurm_script, args.job_name, args.array)
