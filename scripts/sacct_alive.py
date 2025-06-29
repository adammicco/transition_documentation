import argparse
import datetime
import subprocess
import sys
import os
import re
from time import sleep
from math import ceil
import random
from restart_job import restart_cromwell_SLURM_job
from submit import partition

STATE_RUNNING = [ # states where we act as if the job is running
	'PENDING',
	'RUNNING',
	'REQUEUED',
	'RESIZING',
	'REVOKED',
	'SUSPENDED'
]

STATE_RESUBMIT = [ # states where we resubmit the job with no changes
	'BOOT_FAIL',
	'DEADLINE',
	'NODE_FAIL',
	'PREEMPTED'
]

STATE_FINAL = [ # states we consider final and let cromwell parse the rc/stderr
	'CANCELLED',
	'COMPLETED', 
	'FAILED',
	'OUT_OF_MEMORY', # I think we expect cromwell to catch this error and resubmit...
]
	
STATE_TIMEOUT = [ # states where we up the time limit and resubmit
	'TIMEOUT'
]

# Running or queued jobs appear in squeue output
# Failed jobs do not appear
def increment_time(workDir):
	with open(f'{workDir}/execution/script.submit', 'r') as f:
		filedata = f.read()
		time_partition_list = re.findall(r'(-t [0-9]+ -p [a-z]+)',filedata,re.MULTILINE)[0].strip().split(' ')
	curr_time = time_partition_list[1]
	curr_partition = time_partition_list[3]
	new_time = ceil(int(curr_time) * args.time_increment)
	new_time_partition_list = ['-t', str(new_time), '-p', partition(new_time, curr_partition)]
	filedata = filedata.replace(' '.join(time_partition_list), ' '.join(new_time_partition_list))
	with open(f'{workDir}/execution/script.submit', 'w') as f:
		f.write(filedata)
	try:
		os.remove(f'{workDir}/execution/rc')
	except (FileNotFoundError):
		pass
	pass

def success(job, rc):
	sys.exit(rc)

def log_message(message):
	timed_message = f'{datetime.datetime.now()}\t{message}'
	with open('sacct_alive_log', 'a') as f:
		print(timed_message, file=f, flush=True)
		print(timed_message) 

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Check whether a SLURM job is alive, submitting multiple queries if one fails. Returns 0 upon squeue successfully running at least once, -1 upon check failure")
	parser.add_argument('-n', '--max_queries', type=int, help='Maximum number of queries to submit to check job status', default=4)
	parser.add_argument('-l', '--min_delay', type=float, help='Minimum base delay between queries in seconds', default=150)
	parser.add_argument('-m', '--max_delay', type=float, help='Maximum base delay between queries in seconds', default=300)
	parser.add_argument('-t', '--time_increment', type=float, help='Multiplication factor for runtime resubmission.', default=1.50)
	parser.add_argument('slurm_job_number', type=int, help='SLURM job number to check whether alive(running)')
	args = parser.parse_args()

	# check if slurm_job_number is the same as the job number in the stdout.submit script. If not, replace slurm_job_number with the job id from stdout.submit but print both job numbers in the stdout for when the job is running or resubmitted (will need an edit to restart_cromwell_SLURM_job())

	sacct_cmd = ['sacct', '-j', str(args.slurm_job_number), '-u', os.getenv('USER'), '--format=JobID,WorkDir,State', '-P', '-X', '-n']

	for n in range(args.max_queries):
		alias_job_txt =''
		if n > 0:
			delay = n * random.uniform(args.min_delay, args.max_delay)
			sleep(delay)
		# run sacct to get JobID|WorkDir|State
		sacct_result = subprocess.run(sacct_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)

		# sacct succeeds w/ no communication errors
		if sacct_result.returncode == 0:
			try:
				retJobId, workDir, state = sacct_result.stdout.splitlines()[0].split('|')
				state = state.split()[0]
				# print(state)
			except (IndexError): # sacct succeeds but the record of the job is gone
				log_message(f'no sacct information for job {args.slurm_job_number}')
				sys.exit(-1) # job is not running but we are not able to get workdir
			
			with open(f'{workDir}/execution/stdout.submit') as f:
				filedata = f.readlines()
				stdout_submit_jobid = re.findall(r'([0-9]+$)',filedata[0],re.MULTILINE)[0].strip()
				if stdout_submit_jobid != str(args.slurm_job_number):
					log_message(f'Called job ID was {args.slurm_job_number} and is now {stdout_submit_jobid}')
					sacct_cmd = ['sacct', '-j', str(stdout_submit_jobid), '-u', os.getenv('USER'), '--format=JobID,WorkDir,State', '-P', '-X', '-n']
					sacct_result = subprocess.run(sacct_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
					if sacct_result.returncode == 0:
						try:
							retJobId, workDir, state = sacct_result.stdout.splitlines()[0].split('|')
							state = state.split()[0]
							alias_job_txt = f"in place of restarted job {str(args.slurm_job_number)}"
							# print(state)
						except (IndexError): # sacct succeeds but the record of the job is gone
							log_message(f'no sacct information for job {args.slurm_job_number}')
							sys.exit(-1) # job is not running but we are not able to get workdir
					else:
						continue
			if state in STATE_FINAL:
				success(args.slurm_job_number, 0)
				# sys.exit(0) # job is not running and is in a final state
			elif state in STATE_RUNNING:
				log_message('\t'.join([retJobId, workDir, state, alias_job_txt]))
				success(args.slurm_job_number, 0)
				# sys.exit(0) # job is running -- all is good
			elif state in STATE_RESUBMIT:
				restart_cromwell_SLURM_job(f'{workDir}/execution')
				log_message('\t'.join([retJobId, workDir, state, alias_job_txt]))
				success(args.slurm_job_number, 0)
				# sys.exit(0)
			elif state in STATE_TIMEOUT:
				increment_time(workDir)
				restart_cromwell_SLURM_job(f'{workDir}/execution')
				log_message('\t'.join([retJobId, workDir, state, alias_job_txt.replace("restarted", "timed-out")]))
				success(args.slurm_job_number, 0)
				# sys.exit(0)
	# If we get to this point, all sacct checks have resulted in an error
	log_message(f'sacct had return code {sacct_result.returncode}: {sacct_result.stderr}', file=sys.stderr)
	sys.exit(-1)
