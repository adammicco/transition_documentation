task pulldown{
	File bamlist
	File report
	File pulldown_executable
	File python_pulldown
	Boolean ignore_experiment
	
	String release_parent_directory
	String label
	String snp_set

	command{
		python3 ${python_pulldown} --pulldown_executable ${pulldown_executable} --pulldown_label ${label} --release_directory ${release_parent_directory} ${if ignore_experiment then "--ignore_experiment" else ""} ${if (snp_set != "") then "--snp_set " + snp_set else ""} ${bamlist} ${report}
	}
	output{
		File pulldown_ind = "${label}.combined.ind"
		File pulldown_snp = "${label}.combined.snp"
		File pulldown_geno = "${label}.combined.geno"
		File dblist = "${label}.dblist"
	}
	runtime{
		cpus: 2
		requested_memory_mb_per_core: 8000
	}
}

task generate_fstat_input{
	File comparison_set_1_ind
	File comparison_set_1_snp
	File comparison_set_1_geno

	File comparison_set_2_ind
	File comparison_set_2_snp
	File comparison_set_2_geno

	File outgroup_ind
	File outgroup_snp
	File outgroup_geno

	File pulldown_ind
	File pulldown_snp
	File pulldown_geno

	String label

	File mergemany_executable
	
	
	command <<<
		# Assign pop labels to match _i and _d samples and symlink geno, snp files
		awk '{print $1, $2, $1}' ${pulldown_ind} > combined.ind
		ln -s ${pulldown_snp} combined.snp
		ln -s ${pulldown_geno} combined.geno

		ln -s ${comparison_set_1_ind} $(basename ${comparison_set_1_ind})
		ln -s ${comparison_set_1_snp} $(basename ${comparison_set_1_snp})
		ln -s ${comparison_set_1_geno} $(basename ${comparison_set_1_geno})

		ln -s ${comparison_set_2_ind} $(basename ${comparison_set_2_ind})
		ln -s ${comparison_set_2_snp} $(basename ${comparison_set_2_snp})
		ln -s ${comparison_set_2_geno} $(basename ${comparison_set_2_geno})

		ln -s ${outgroup_ind} $(basename ${outgroup_ind})
		ln -s ${outgroup_snp} $(basename ${outgroup_snp})
		ln -s ${outgroup_geno} $(basename ${outgroup_geno})

		# Generate merged sample set for F4 input
		${mergemany_executable} combined $(basename ${comparison_set_1_ind} .ind) $(basename ${comparison_set_2_ind} .ind) $(basename ${outgroup_ind} .ind) ${label}_fstat_set

		# Generate list of target libraries.
		comm -12 <(grep '_d[[:blank:]]' ${pulldown_ind} | awk '{print $3}' | sort) <(grep '_i[[:blank:]]' ${pulldown_ind} | awk '{print $3}' | sort) > ${label}_fstat_target_libraries

		# Build list of tests
		comparison1=$(echo $(basename ${comparison_set_1_ind} .ind) $(basename ${outgroup_ind} .ind))
		comparison2=$(echo $(basename ${comparison_set_2_ind} .ind) $(basename ${outgroup_ind} .ind))
		# awk -v comparison="$comparison" '{print $1"_i",$1"_d",comparison}' ${label}_fstat_target_libraries > ${label}_f4_tests
		while read -r lib
			do
				echo $lib'_i '$lib'_d '$comparison1' '$lib'_i '$lib'_d '$comparison2 > $lib'_f4_test'
			done < ${label}_fstat_target_libraries
	>>>
	output{
		File fstat_set_ind = "${label}_fstat_set.ind"
		File fstat_set_snp = "${label}_fstat_set.snp"
		File fstat_set_geno = "${label}_fstat_set.geno"
		Array[File] fstat_pop_list = glob("*_f4_test")
		Array[String] fstat_target_libraries = read_lines("${label}_fstat_target_libraries")
	}
	runtime{
		cpus: 1
		requested_memory_mb_per_core: 36000
	}
}

task run_damage_inverse_fstats{
	File fstat_set_ind
	File fstat_set_snp
	File fstat_set_geno

	File fstat_tests

	File fstat_executable

	command{
		# Build par file
		echo "indivname:	${fstat_set_ind}" >> f4_par
		echo "snpname:	${fstat_set_snp}" >> f4_par
		echo "genotypename:	${fstat_set_geno}" >> f4_par
		echo "popfilename:	${fstat_tests}" >> f4_par
		echo "allsnps:	YES" >> f4_par

		# Run f4
		${fstat_executable} -p f4_par
		
	}

	output {
		File stdout = "stdout"
	}

	runtime{
		cpus: 2
		requested_memory_mb_per_core: 16000
	}
}

task concat_fstat_results {
	File fstat_stdout

	command <<<
		for file in {${sep=',' fstat_stdout}
			do
				grep -hE 'result|badoctet' $file | awk 'name=substr($2, 1, length($2)-2) {print name,"f4_test","i:d:"$4":"$5"::i:d:"$9":"$10,"f4diff",$11,"std_error",$12,"z_score",$13,"snps?",$14}' > concat
			done
	>>>

	output {
		File results = "concat"
	}

	runtime {
		cpus: 1
		requested_memory_mb_per_core: 200
		runtime_minutes: 5
	}
}