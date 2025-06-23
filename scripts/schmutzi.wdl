version 1.0

workflow schmutzi {
	input {
		File input_driver # Input driver file containing a list of mtDNA BAM files and optional parameters in TSV format (ID, bam, library_type?, length_deam_custom?) w/ header

		File mt_reference_rsrs_in = "/n/groups/reich/matt/pipeline/static/mtdna_rsrs.fa"

		File contDeam_pl
		File schmutzi_pl
		String freqs_path

		# may need to adjust these
	# 	Map[String, Int] length_deam_map = {
	# 		'ds.minus' : 40,
	# 		'ds.half' : 30,
	# 		'ds.plus' : 20,
	# 		'ss.minus' : 5,
	# 		'ss.half' : 2,
	# 		'ss.plus' : 1
	# 	}
	# }
		# diagnostic -- trying 2 bp across the board for all library types
		Map[String, Int] length_deam_map = {
			'ds.minus' : 2,
			'ds.half' : 2,
			'ds.plus' : 2,
			'ss.minus' : 2,
			'ss.half' : 2,
			'ss.plus' : 2
		}
	}


	# TODO accept input tsv file of bam \t library_type \t length_deam_custom (optional)

	# refactor so that calmd and contDeam run all bams in one batch since they are computationally cheap\
	# scatter schmutzi on the output of contDeam so that there are 4 schmutzi jobs running in parallel for each bam
	# scatter (bam in mt_bams) {
	# 	call calmd {
	# 		input:
	# 			bam = bam,
	# 			mt_reference = mt_reference_rsrs_in
	# 	}

	# 	call contDeam {
	# 		input:
	# 			bam = calmd.calmd_bam,
	# 			mt_reference = mt_reference_rsrs_in,
	# 			contDeam_pl = contDeam_pl,
	# 			length_deam_map = length_deam_map
	# 	}

	call calmd_contDeam {
		input:
			driver_file = input_driver,
			mt_reference = mt_reference_rsrs_in,
			contDeam_pl = contDeam_pl,
			length_deam_map = length_deam_map
	}

	scatter(paths_line in calmd_contDeam.paths_out) {
	# run schmutzi using --notusepredC option without --uselength option
		call schmutzi as schmutzi_notusepredC {
			input:
				paths_line = paths_line,
				mt_reference = mt_reference_rsrs_in,
				schmutzi_pl = schmutzi_pl,
				output_prefix = "notusepredC",
				freqs_path = freqs_path,
				notusepredC = true,
				uselength = false
		}

	# run schmutzi not using --notusepredC option without --uselength option
		call schmutzi as schmutzi_usepredC {
			input:
				paths_line = paths_line,
				mt_reference = mt_reference_rsrs_in,
				schmutzi_pl = schmutzi_pl,
				output_prefix = "usepredC",
				freqs_path = freqs_path,
				notusepredC = false,
				uselength = false
	}

	# run schmutzi using --notusepredC option with --uselength option
		call schmutzi as schmutzi_notusepredC_uselength {
			input:
				paths_line = paths_line,
				mt_reference = mt_reference_rsrs_in,
				schmutzi_pl = schmutzi_pl,
				output_prefix = "notusepredC_uselength",
				freqs_path = freqs_path,
				notusepredC = true,
				uselength = true
		}

	# run schmutzi not using --notusepredC option with --uselength option
		call schmutzi as schmutzi_usepredC_uselength {
			input:
				paths_line = paths_line,
				mt_reference = mt_reference_rsrs_in,
				schmutzi_pl = schmutzi_pl,
				output_prefix = "usepredC_uselength",
				freqs_path = freqs_path,
				notusepredC = false,
				uselength = true
		}
	}
}

# task calmd {
#     input {
#         # Array[File] bams
# 		File input_driver

#         File mt_reference
# 		File mt_reference_fai = mt_reference + ".fai"
# 		File mt_reference_amb = mt_reference + ".amb"
# 		File mt_reference_ann = mt_reference + ".ann"
# 		File mt_reference_bwt = mt_reference + ".bwt"
# 		File mt_reference_pac = mt_reference + ".pac"
# 		File mt_reference_sa = mt_reference + ".sa" 
#     }

# 	Array[File] bam_files = schmutzi_inputs.map(input => input.bam)

#     command <<<
# 		set -e
# 		for bam in ~{sep=' ' bam_files}
# 			do
# 				samtools calmd -bAr ~{bam} ~{mt_reference} > $(basename ~{bam} .bam).calmd.bam
# 				samtools index $(basename ~{bam} .bam).calmd.bam
# 			done
#     >>>

#     output {
# 		Array[File] calmd_bams = glob("*.calmd.bam")
#     }

#     runtime {
# 		# TODO: optimize this based on the number/size of bams being processed. Also, measure the memory usage of samtools calmd to determine a more appropriate memory allocation.
# 		cpus: 1
# 		requested_memory_mb_per_core: 4000
# 		runtime_minutes: 5
#     }
# }

task calmd_contDeam {
	input {
		File driver_file

		File mt_reference
		File mt_reference_fai = mt_reference + ".fai"
		File mt_reference_amb = mt_reference + ".amb"
		File mt_reference_ann = mt_reference + ".ann"
		File mt_reference_bwt = mt_reference + ".bwt"
		File mt_reference_pac = mt_reference + ".pac"
		File mt_reference_sa = mt_reference + ".sa"

		Map[String, Int] length_deam_map

		File contDeam_pl
		
		Int processes = 8
	}

	Int num_files = length(read_lines(driver_file))
	Int minutes = 10

	command <<<
		set -e
		python3 <<CODE
		import subprocess
		import json
		from multiprocessing import Pool
		from os.path import basename, exists
		from pathlib import Path
		import re
		from collections import Counter
		from time import sleep



		FIELD_MAPPING = {
			"ID"                    : 0,
			"bam"                   : 1,
			"histogram"             : 2,
			"library_type"          : 3,
			"length_deam_custom"    : 4,
		}

		# variables passed from WDL
		DRIVER_PATH = "~{driver_file}"
		MT_REF = "~{mt_reference}"
		LENGTH_DEAM_MAP = json.loads(open("~{write_json(length_deam_map)}", 'r').read())
		CONTDEAM_PL = "~{contDeam_pl}"
		PROCESSES = int("~{processes}")
		# =========================

		def library_type_from_histogram(histogram):
			library_types = []
			with open(histogram, 'r') as hist_file:
				while True:
					line = hist_file.readline().strip()
					if not line:
						break
					match = re.match(r"^##\s{1}histogram\s{1}\d+\s+:\s+([a-zA-Z]*:[a-zA-Z]*)", line)
					if not match:
						continue
					library_types.append(match.group(1).replace(":",".").replace("single", "ss").replace("double", "ds"))
			return Counter(library_types).most_common(1)[0][0] # return the most common library type

		def run_calmd_contDeam(id, bam, histogram, library_type, length_deam_custom):
			"""Run the calmd and contDeam process for a single BAM file."""
			try:
				# Run samtools calmd
				calmd_bam = f"{basename(bam).replace('.bam', '')}.calmd.bam"
				subprocess.run(
					["samtools", "calmd", "-bAr", bam, MT_REF],
					check=True,
					stdout=open(calmd_bam, 'wb')
				)
				subprocess.run(
					["samtools", "index", calmd_bam],
					check=True
				)

				# logic for selecting the appropriate length_deam:
				# length_deam_custom > library_type w/ mapping > histogram majority vote > ds.half
				# for library_type:
				# library_type > histogram majority vote > ds.half

				# Determine library type if not provided
				print(id, bam, histogram, library_type, length_deam_custom)
				if not library_type:
					print(f"Library type not provided for {id}, determining from histogram...")
					if not histogram:
						print(f"Warning: No histogram provided for {id}, defaulting to ds.half")
						library_type = "ds.half"
					else:
						print(f"Determining library type from histogram for {id}...")
						library_type = library_type_from_histogram(histogram)

				# Determine length_deam
				if length_deam_custom:
					print(f"Custom length_deam provided for {id}: {length_deam_custom}")
					length_deam = length_deam_custom
				else:
					print(f"Using library type for length_deam determination for {id}: {library_type}")
					length_deam = LENGTH_DEAM_MAP[library_type]

				library_type = library_type.split('.')[0].replace("ss", "single").replace("ds", "double")  # Normalize library type for output

				out_dir = f"contDeam/{id}_contDeam"
				Path.mkdir(Path(out_dir), exist_ok=True, parents=True)

				subprocess.run(
					[CONTDEAM_PL, "--lengthDeam", str(length_deam), "--library", library_type, "--out", f"{out_dir}/{id}_contDeam", MT_REF, calmd_bam]
				)

				return [str(Path.absolute(Path(calmd_bam))), str(Path.absolute(Path(f"{out_dir}/{id}_contDeam")))]
			except Exception as e:
				print(f"Error processing {id}: {e}")

		def get_from_list_or_return_none(lst, index):
			"""Helper function to safely get an item from a list or return None if index is out of range."""
			try:
				return lst[index]
			except IndexError:
				return None

		with open(DRIVER_PATH, 'r') as driver_file:
			driver_lines = [x.split() for x in driver_file.readlines()] # do not strip whitespace since there are optional fields

		pool = Pool(processes=PROCESSES)
		results = []

		for line in driver_lines:
			# Skip empty lines
			if not line:
				continue
			
			id = get_from_list_or_return_none(line, FIELD_MAPPING["ID"])
			bam = get_from_list_or_return_none(line, FIELD_MAPPING["bam"])
			histogram = get_from_list_or_return_none(line, FIELD_MAPPING["histogram"])
			library_type = get_from_list_or_return_none(line, FIELD_MAPPING["library_type"])
			length_deam_custom = get_from_list_or_return_none(line, FIELD_MAPPING["length_deam_custom"])

			results.append(pool.apply_async(run_calmd_contDeam, args=(id, bam, histogram, library_type, length_deam_custom)))

		pool.close()
		paths_out = [ result.get() for result in results]  # Wait for all processes to complete
		sleep(10)
		open("calmd_paths.txt", 'w').write("\n".join(["\t".join(x) for x in paths_out]) + '\n')
		CODE
	>>>

	output {
		# Array[File] calmd_bams = glob("*.calmd.bam")
		# Array[] contDeam_output_dir = glob("contDeam/*")
		Array[String] paths_out = read_lines("calmd_paths.txt")
	}

	runtime {
		# TODO: may need to optimize this
		cpus: if num_files < processes then num_files else processes
		requested_memory_mb_per_core: 1000
		runtime_minutes: minutes
		queue: if minutes > 720 then "medium" else "short"
	}
}

task schmutzi {
	input {
		String paths_line
		# File bam
		# File bam_bai = bam + ".bai"
		# String contDeam_output_dir # may need to change this to a String if Cromwell doesn't like a directory as input for a File type
		File schmutzi_pl
		String output_prefix
		String freqs_path
		Boolean notusepredC
		Boolean uselength

        File mt_reference
		File mt_reference_fai = mt_reference + ".fai"
		File mt_reference_amb = mt_reference + ".amb"
		File mt_reference_ann = mt_reference + ".ann"
		File mt_reference_bwt = mt_reference + ".bwt"
		File mt_reference_pac = mt_reference + ".pac"
		File mt_reference_sa = mt_reference + ".sa" 
	}

	command <<<
		set -e
		module load R/3.5.1

		bam=$(echo ~{paths_line} | awk '{print $1}') # get the first field from the paths_line input
		bam_bai=$(echo ${bam}).bai

		bam_base=$(basename ${bam})

		ln -s ${bam} .
		ln -s ${bam_bai} .

		contDeam_output_prefix=$(echo ~{paths_line} | awk '{print $2}') # get the second field from the paths_line input
		contDeam_base=$(basename ${contDeam_output_prefix})

		mkdir -p contDeam
		for file in ${contDeam_output_prefix}*
			do
			if [ -f contDeam/$(basename $file) ]; then
				rm contDeam/$(basename $file)
			fi
			ln -s ${file} contDeam/$(basename $file)
		done

		~{schmutzi_pl} -t 8 --estdeam ~{if notusepredC then '--notusepredC' else ''} ~{if uselength then '--uselength' else ''} --ref ~{mt_reference} --out ~{output_prefix} contDeam/${contDeam_base} ~{freqs_path} ${bam_base}
	>>>

	output {
		File schmutzi_output = glob("*_final*")
	}

	runtime {
		# TODO: optimize this
		cpus: 8
		requested_memory_mb_per_core: 4000
		runtime_minutes: 720
	}
}
