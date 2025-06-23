import argparse
import subprocess
import sys

from time import sleep
import random

def partition(runtime_minutes, partition):
	if partition == 'priority':
		return partition
	if runtime_minutes <= 720: # 12 hours
		return 'short'
	elif runtime_minutes <= 7200: # 5 days
		return 'medium'
	else:
		return 'long'

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Submit a command to SLURM with resubmissions to handle timeout")
	parser.add_argument('-s', '--max_resubmit', type=int, help='Maximum number of times to submit command to SLURM', default=3)
	parser.add_argument('-l', '--min_delay', type=float, help='Minimum base delay between queries in seconds', default=60)
	parser.add_argument('-m', '--max_delay', type=float, help='Maximum base delay between queries in seconds', default=120)

	# sbatch parameters
	parser.add_argument('-J', '--job_name', required=True)
	parser.add_argument('-D', '--working_directory', default='./')
	parser.add_argument('-o', '--stdout', required=True)
	parser.add_argument('-e', '--stderr', required=True)
	parser.add_argument('-t', '--runtime', type=int, required=True)
	parser.add_argument('-p', '--partition')
	parser.add_argument('-n', '--cpus', type=int, default=1)

	memory = parser.add_mutually_exclusive_group(required=True)
	memory.add_argument('--mem_per_cpu', type=int, help='memory per cpu requested in mb')
	memory.add_argument('--mem', type=int, help='memory requested in mb')

	parser.add_argument('script', help='script for SLURM to run')

	args = parser.parse_args()

	submit_args = ['sbatch']
	submit_args += ['-J', args.job_name]
	submit_args += ['-D', args.working_directory]
	submit_args += ['-o', args.stdout]
	submit_args += ['-e', args.stderr]
	submit_args += ['-t', str(args.runtime)]
	submit_args += ['-p', partition(args.runtime, args.partition)]
	submit_args += ['-n', str(args.cpus)]
	if args.mem_per_cpu:
		submit_args += ['--mem-per-cpu', str(args.mem_per_cpu)]
	else:
		submit_args += ['--mem', str(args.mem)]
	submit_args += ['--constraint', 'groups']
	submit_args += ['--qos', 'ded_reich']
	submit_args += ['--account', 'reich']
	submit_args += ['--wrap', f"/usr/bin/env bash {args.script}"]
	# print(" ".join(submit_args))
	for n in range(args.max_resubmit):
		if n > 0:
			delay = n * random.uniform(args.min_delay, args.max_delay)
			sleep(delay)
		# Suspect that Cromwell is looking for one instance of the job number in stdout
		# stdout is uncaptured so Cromwell gets it
		result = subprocess.run(submit_args, stderr=subprocess.PIPE, universal_newlines=True)
		# any successful query will tell us whether job is alive
		if result.returncode == 0:
			sys.exit(0)
		else:
			# print(f'squeue failed with return code {result.returncode}', file=sys.stderr)
			if "sbatch: error: Batch job submission failed: Socket timed out on send/recv operation" in result.stderr: # Job was likely submitted despite error from SLURM
				sacct_cmd = ['sacct', '-u', "os.getenv('USER')", '--format=JobID,WorkDir', '-P', '-X']
				sacct_result = subprocess.run(sacct_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True) # run sacct in parsable mode to report working directories of running jobs
				for line in sacct_result.stdout:
					feilds = line.strip().split('|')
					if feilds[1] == args.working_directory:  # search for a job that matches our expected working directory...
						print(f'Submitted batch job {feilds[0]}', file=f'{args.stdout}.submit') # ...if we find one, fake a SLURM-style job submission confirmation in stdout.submit
						print(result.stderr, file=f'{args.working_directory}/sacct.submit') # log that we found the errant job via sacct in case that matters later
						print(f'Job {feilds[0]} found via sacct', file=f'{args.working_directory}/sacct.submit')
						sys.exit(0)
			# if 'Invalid job id specified' in result.stderr: # job is not alive
			# 	sys.exit(result.returncode) # reporting failure appears get Cromwell to resubmit
	print(result.stderr, file=sys.stderr)
	sys.exit(-1)
