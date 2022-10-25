#!/usr/bin/env nextflow
// This version is for deployment on a range of executors.
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.Files

def user2containerPath(base_dir, user_path, container_base_dir) {
	// Throws an error if user_path is not in the subdirectory tree of base_dir
	try{
		relative_user_path = base_dir.relativize(Paths.get(user_path))
	} catch(Exception e) {
		println "Error converting from user_path to container_path"
		println "Check user_path is in the subdirectory tree of base_dir"
		println "base_dir:  " + base_dir
		println "user_path: " + user_path
		throw e
	}
	Paths.get(container_base_dir, relative_user_path.toString()).toString()
}

// Determine mountpoint and container command
def deployment() {
	// Determine container type and file location
	params.container_type = "none"
	params.container_file = ""
	if(params.container_file=="" && params.container_type.toString().toLowerCase()!="none") throw new Exception("container_file must be specified for container_type!='none'")
	if(params.container_file.toString().toLowerCase()=="singularity" && !Files.exists(Paths.get(params.container_file))) throw new Exception("container_file ${params.container_file} does not exist")

	// Identify where the user file system is mounted to the container file system
	if(params.container_type.toString().toLowerCase()=="none") {
		params.container_mount = params.base_dir
	} else {
		params.container_mount = "/home/jovyan"
	}
	params.container_args = ""

	// Check io files exist
	if(!Files.exists(Paths.get(params.base_dir))) throw new Exception("base_dir ${params.base_dir} does not exist")
	if(!Files.exists(Paths.get(params.id_file))) throw new Exception("id_file ${params.id_file} does not exist")
	if(!Files.exists(Paths.get(params.ref_fa))) throw new Exception("ref_fa ${params.ref_fa} does not exist")
	if(!Files.exists(Paths.get(params.ref_gb))) throw new Exception("ref_gb ${params.ref_gb} does not exist")

	// Convert io files from user file system to container file system
	base_dir = Paths.get(params.base_dir)
	params.container_analysis_dir = user2containerPath(base_dir, params.analysis_dir, params.container_mount)
	params.container_id_file = user2containerPath(base_dir, params.id_file, params.container_mount)
	if(params.covariate_file=="") {
		params.container_covariate_file = ""
	} else {
		if(!Files.exists(Paths.get(params.covariate_file))) throw new Exception("covariate_file ${params.covariate_file} does not exist")
		params.container_covariate_file = user2containerPath(base_dir, params.covariate_file, params.container_mount)
	}
	params.container_ref_fa = user2containerPath(base_dir, params.ref_fa, params.container_mount)
	params.container_ref_gb = user2containerPath(base_dir, params.ref_gb, params.container_mount)

	// Convert user-specified software_file from user file system to container file system
	assert !binding.hasVariable('params.default_software_file')  // Do not allow override !
	assert !binding.hasVariable('params.default_script_dir')  // Do not allow override !
	params.default_software_file = "/usr/share/kmer_pipeline/example/pipeline_software_location.txt"
	params.default_script_dir = "/usr/local/bin"
	params.software_file = params.default_software_file
	if(params.software_file.toString().toLowerCase()==params.default_software_file) {
		// If the default software file is retained, assume the path is internal to the container
		params.container_software_file = params.software_file
	} else {
		// Otherwise assume it is on the user file system and convert path
		if(!Files.exists(Paths.get(params.software_file))) throw new Exception("software_file ${params.software_file} does not exist")
		params.container_software_file = user2containerPath(base_dir, params.software_file, params.container_mount)
	}

	// Implied container deployment command (can be overridden explicitly)
	container_type = params.container_type.toString().toLowerCase()
	switch(container_type) {
	case "none":
		params.container_cmd = ""
		break;
	case "singularity":
		assert params.container_file != ""
		params.container_cmd = "singularity exec --containall --cleanenv --home ${params.base_dir}:${params.container_mount} ${params.container_args} ${params.container_file}"
		break;
	case "docker":
		assert params.container_file != ""
		params.container_cmd = "docker run --rm -v ${params.base_dir}:${params.container_mount} ${params.container_args} ${params.container_file}"
		break;
	default:
		throw new Exception("Container type '${params.container}' not recognised")
	}
}

// Read the software location file
// *** Assume the software_file itself is in the user file system ***
// *** BUT all paths in software_file are in the container file system ***
def read_container_script_dir() {
	// By default the software file is stored within the container - next line avoids reading it from outside the container
	if(params.software_file.toString().toLowerCase()==params.default_software_file) return params.default_script_dir
	// Otherwise the software file is specified in the user file system, so it can be read directly
	def software_file = new File(params.software_file)
	assert software_file.readLines().head().split('\t')*.toLowerCase() == ['name','path']
	def rows = software_file.readLines().tail()*.split('\t')
	rows.collectEntries( {[it[0],it[1]]} )["scriptpath"]
}

// Read the name of the reference genome
def read_ref_name() {
	cmd = "${params.container_cmd} ${params.container_script_dir}/get_ref_name.Rscript ${params.container_ref_fa}"
	proc = cmd.execute()
	sout = new StringBuilder()
	ref_name_read_error = new StringBuilder()
	proc.waitForProcessOutput(sout, ref_name_read_error)
	assert ref_name_read_error.toString()==""
	assert sout.toString()!=""
	params.ref_name = sout.toString().trim()
}

// Get sample size (irrespective of phenotype)
def get_n() {
	cmd = "${params.container_cmd} wc -l ${params.container_id_file}"
	proc = cmd.execute()
	outputStream = new StringBuffer();
	proc.waitForProcessOutput(outputStream, System.err)
	Scanner scn = new Scanner(outputStream.toString())
	assert(scn.hasNextInt())
	return scn.nextInt()-1
}

// Processes
process check_params {
output:
	val true, emit: done
shell:
	'''
	echo "Prelim: Checking parameters ..."
	echo "Not yet implemented"
	'''
}

process create_analysis_dir {
input:
	val ready
output:
	val true, emit: done
shell:
	'''
	echo "Prelim: Creating analysis directories ..."
	echo "!{params.container_analysis_dir}"
	mkdir -p !{params.analysis_dir}
	mkdir -p !{params.workdir}
	mkdir -p !{params.logdir}
	'''
}

process create_analysis_file {
input:
	val ready
output:
	val true, emit: done
shell:
	'''
	echo "Prelim: Defining analysis file ..."
	echo "kmertype\tkmerlength" > !{params.analysis_file}
	echo "!{params.kmer_type}\t!{params.kmer_length}" >> !{params.analysis_file}
	'''
}

// Define required input files
def read_id_file() {
	def id_file = new File(params.id_file)
	assert id_file.readLines().head().split('\t')*.toLowerCase() == ['id','paths','pheno']
	def rows = id_file.readLines().tail()*.split('\t')
	id = rows.collect( { row -> row[0] })
	paths = rows.collect( { row -> row[1] })
	pheno = rows.collect( { row -> row[2] })
	ret = [id: id, paths: paths, pheno: pheno]
}

process filecheck_prelim {
input:
	val ready
	path inputFASTAs
	path ref_fa
	path ref_gb
output:
	val true, emit: done
shell:
	'''
	echo "Fchk 0: Initial filecheck"
	echo "inputFASTAs: !{inputFASTAs}"
	echo "ref_fa: !{ref_fa}"
	echo "ref_gb: !{ref_gb}"
	'''
}

// Step 1: Count kmers
// min(n,p)-fold parallelization
process countkmers {
input:
	val ready
	val sampid
output:
	val true, emit: done
	val sampid
shell:
if(!params.skip1)
	'''
	echo "Step 1: Counting kmers"
	echo "sampid: !{sampid}"
	ln -sfr $(pwd) !{params.workdir}/countkmers.!{sampid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/countkmers.!{sampid}.log
	!{params.container_cmd} !{params.container_script_dir}/countkmers.Rscript \
		!{sampid} \
		!{params.container_id_file} \
		!{params.container_analysis_dir} \
		!{params.output_prefix} \
		!{params.container_software_file} \
		!{params.container_analysis_file}
	'''
else
	'''
	echo "Skipping Step 1: Counting kmers"
	'''
}

def outfiles_countkmers(id_list) {
	[kmercounts:
		id_list['id'].collect({ lab -> "${params.kmer_type}kmer${params.kmer_length}/${lab}.kmer${params.kmer_length}.txt.gz" }),
	kmertotals:
		id_list['id'].collect({ lab -> "${params.kmer_type}kmer${params.kmer_length}/${lab}.kmer${params.kmer_length}.total.txt" }),
	kmers_filepaths:
		"${params.output_prefix}_${params.kmer_type}${params.kmer_length}_kmers_filepaths.txt"
	]
}

process filecheck_countkmers {
input:
	val ready
	val max_sampid
	path kmercounts
	path kmertotals
	path kmers_filepaths
output:
	val true, emit: done
shell:
if(!params.skip1 | !params.skip2)
	'''
	echo "Fchk 1: Counting kmers"
	echo "max_sampid: !{max_sampid}"
	echo "kmercounts: !{kmercounts}"
	echo "kmertotals: !{kmertotals}"
	echo "kmers_filepaths: !{kmers_filepaths}"
	ln -sfr $(pwd) !{params.workdir}/filecheck_countkmers
	ln -sfr $(pwd)/.command.log !{params.logdir}/filecheck_countkmers.log
	'''
else
	'''
	echo "Skipping Fchk 1: Counting kmers"
	'''
}

def outfiles_countkmers_protein(id_list) {
	[kmercounts:
		id_list['id'].collect({ lab -> "${params.kmer_type}kmer${params.kmer_length}/${lab}.kmer${params.kmer_length}.txt.gz" }),
	kmertotals:
		id_list['id'].collect({ lab -> "${params.kmer_type}kmer${params.kmer_length}/${lab}.kmer${params.kmer_length}.total.txt" }),
	kmers_filepaths:
		"${params.output_prefix}_${params.kmer_type}${params.kmer_length}_kmers_filepaths.txt",
	reading_frames:
		id_list['id'].collect({ lab -> "translated_contigs/${lab}_translated_all_reading_frames.fa.gz" })
	]
}

process filecheck_countkmers_protein {
input:
	val ready
	val max_sampid
	path kmercounts
	path kmertotals
	path reading_frames
	path kmers_filepaths
output:
	val true, emit: done
shell:
if(!params.skip1 | !params.skip2)
	'''
	echo "Fchk 1: Counting kmers (proteins)"
	echo "max_sampid: !{max_sampid}"
	echo "kmercounts: !{kmercounts}"
	echo "kmertotals: !{kmertotals}"
	echo "reading_frames: !{reading_frames}"
	echo "kmers_filepaths: !{kmers_filepaths}"
	ln -sfr $(pwd) !{params.workdir}/filecheck_countkmers_protein
	ln -sfr $(pwd)/.command.log !{params.logdir}/filecheck_countkmers_protein.log
	'''
else
	'''
	echo "Skipping Fchk 1: Counting kmers (proteins)"
	'''
}

// Step 2: Create unique kmer list
//   Parallel pyramid (max p-fold)
process createfullkmerlist {
input:
	val ready
	val taskid
output:
	val true, emit: done
	val taskid
shell:
if(!params.skip2)
	'''
	echo "Step 2: Creating unique kmer list"
	echo "taskid: !{taskid}"
	ln -sfr $(pwd) !{params.workdir}/createfullkmerlist.!{taskid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/createfullkmerlist.!{taskid}.log
	!{params.container_cmd} !{params.container_script_dir}/createfullkmerlist.Rscript \
		!{taskid} \
		!{params.n} \
		!{params.p} \
		!{params.output_prefix} \
		!{params.container_analysis_dir} \
		!{params.container_id_file} \
		!{params.kmer_type} \
		!{params.kmer_length} \
		!{params.container_software_file}
	'''
else
	'''
	echo "Skipping Step 2: Creating unique kmer list"
	'''
}

def outfiles_createfullkmerlist() {
	[kmermerge: "${params.kmerFilePrefix}.kmermerge.txt.gz"]
}


process filecheck_createfullkmerlist {
input:
	val ready
	path kmermerge
output:
	val true, emit: done
shell:
if(!params.skip2 | !params.skip3)
	'''
	echo "Fchk 2: Creating unique kmer list"
	echo "kmermerge: !{kmermerge}"
	ln -sfr $(pwd) !{params.workdir}/filecheck_createfullkmerlist
	ln -sfr $(pwd)/.command.log !{params.logdir}/filecheck_createfullkmerlist.log
	'''
else
	'''
	echo "Skipping Fchk 2: Creating unique kmer list"
	'''
}

// Step 3: Create kmer presence/absence patterns and kinship matrix
//   Parallel pyramid (max p-fold)
process stringlist2patternandkinship {
input:
	val ready
	val taskid
output:
	val true, emit: done
	val taskid
shell:
if(!params.skip3)
	'''
	echo "Step 3: Creating kmer presence/absence patterns and kinship matrix"
	ln -sfr $(pwd) !{params.workdir}/stringlist2patternandkinship.!{taskid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/stringlist2patternandkinship.!{taskid}.log
	!{params.container_cmd} !{params.container_script_dir}/stringlist2patternandkinship.Rscript \
		!{taskid} \
		!{params.p} \
		!{params.container_id_file} \
		!{params.kmerFilePrefix}.kmermerge.txt.gz \
		!{params.kmerFilePrefix}_kmers_filepaths.txt \
		!{params.container_analysis_dir} \
		!{params.output_prefix} \
		!{params.kmer_type} \
		!{params.container_software_file} \
		!{params.kmer_length} \
		!{params.min_count}
	'''
else
	'''
	echo "Skipping Step 3: Creating kmer presence/absence patterns and kinship matrix"
	'''
}

def outfiles_stringlist2patternandkinship() {
	[patternKey: "${params.kmerFilePrefix}.patternmerge.patternKey.txt.gz",
	patternIndex: "${params.kmerFilePrefix}.patternmerge.patternIndex.txt.gz",
	presenceCount: "${params.kmerFilePrefix}.patternmerge.presenceCount.txt.gz",
	patternKeySize: "${params.kmerFilePrefix}.patternmerge.patternKeySize.txt",
	kinship: "${params.kmerFilePrefix}.kinship.txt.gz",
	kinshipWeight: "${params.kmerFilePrefix}.kinshipWeight.txt"]
}

process filecheck_stringlist2patternandkinship {
input:
	val ready
	path patternKey
	path patternIndex
	path presenceCount
	path patternKeySize
	path kinship
	path kinshipWeight
output:
	val true, emit: done
shell:
if(!params.skip3 | !params.skip4)
	'''
	echo "Fchk 3: Creating kmer presence/absence patterns and kinship matrix"
	echo "patternKey: !{patternKey}"
	echo "patternIndex: !{patternIndex}"
	echo "presenceCount: !{presenceCount}"
	echo "patternKeySize: !{patternKeySize}"
	echo "kinship: !{kinship}"
	echo "kinshipWeight: !{kinshipWeight}"
	ln -sfr $(pwd) !{params.workdir}/filecheck_stringlist2patternandkinship
	ln -sfr $(pwd)/.command.log !{params.logdir}/filecheck_stringlist2patternandkinship.log
	'''
else
	'''
	echo "Skipping Fchk 3: Creating kmer presence/absence patterns and kinship matrix"
	'''
}

// Step 4: Run GEMMA
//   p-fold parallelization
process rungemma {
//	publishDir "${params.container_analysis_dir}", mode: 'rellink'
//	stageInMode 'rellink'
input:
	val ready
	val taskid
output:
	val true, emit: done
	val taskid
shell:
if(!params.skip4 && params.container_covariate_file=="")
	'''
	echo "Step 4: Running GEMMA"
	ln -sfr $(pwd) !{params.workdir}/rungemma.!{taskid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/rungemma.!{taskid}.log
	!{params.container_cmd} !{params.container_script_dir}/rungemma.Rscript \
		!{taskid} \
		!{params.p} \
		!{params.kmerFilePrefix} \
		!{params.container_id_file} \
		!{params.output_prefix} \
		!{params.container_analysis_dir} \
		!{params.kmer_type} \
		!{params.kmer_length} \
		!{params.container_software_file}
	cp --remove-destination $(pwd)/.command.log !{params.logdir}/rungemma.!{taskid}.log
	'''
else if(!params.skip4 && params.container_covariate_file!="")
	'''
	echo "Step 4: Running GEMMA"
	ln -sfr $(pwd) !{params.workdir}/rungemma.!{taskid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/rungemma.!{taskid}.log
	!{params.container_cmd} !{params.container_script_dir}/rungemma.Rscript \
		!{taskid} \
		!{params.p} \
		!{params.kmerFilePrefix} \
		!{params.container_id_file} \
		!{params.output_prefix} \
		!{params.container_analysis_dir} \
		!{params.kmer_type} \
		!{params.kmer_length} \
		!{params.container_software_file} \
		!{params.container_covariate_file
	cp --remove-destination $(pwd)/.command.log !{params.logdir}/rungemma.!{taskid}.log
	'''
else
	'''
	echo "Skipping Step 4: Running GEMMA"
	'''
}

def outfiles_rungemma() {
	gemmaout = "${params.container_analysis_dir}/${params.kmer_type}kmer${params.kmer_length}_gemma/output/${params.output_prefix}_${params.kmer_type}${params.kmer_length}"
	[assoc: "${gemmaout}.*.assoc.txt.gz",
	pval: "${gemmaout}.*.pval.txt.gz",
	log: "${gemmaout}.*.log.txt.gz"]
}

/*process proc_outfiles_rungemma {
input:
	val ready
output:
	val true, emit: done
	env assoc
	env pval
	env log
shell:
	'''
	gemmaout="!{params.container_analysis_dir}/!{params.kmer_type}kmer!{params.kmer_length}_gemma/output/!{params.output_prefix}_!{params.kmer_type}!{params.kmer_length}"
	assoc="${gemmaout}.*.assoc.txt.gz"
	pval="${gemmaout}.*.pval.txt.gz"
	log="${gemmaout}.*.log.txt.gz"
	'''
}*/

/*process proc_outfiles_rungemma {
input:
	val ready
output:
	val true, emit: done
path "${params.container_analysis_dir}/${params.kmer_type}kmer${params.kmer_length}_gemma/output/${params.output_prefix}_${params.kmer_type}${params.kmer_length}.*.assoc.txt.gz", emit: assoc
path "${params.container_analysis_dir}/${params.kmer_type}kmer${params.kmer_length}_gemma/output/${params.output_prefix}_${params.kmer_type}${params.kmer_length}.*.pval.txt.gz", emit: pval
path "${params.container_analysis_dir}/${params.kmer_type}kmer${params.kmer_length}_gemma/output/${params.output_prefix}_${params.kmer_type}${params.kmer_length}.*.log.txt.gz", emit: log
shell:
	'''
	gemmaout="!{params.container_analysis_dir}/!{params.kmer_type}kmer!{params.kmer_length}_gemma/output/!{params.output_prefix}_!{params.kmer_type}!{params.kmer_length}"
	assoc="${gemmaout}.*.assoc.txt.gz"
	pval="${gemmaout}.*.pval.txt.gz"
	log="${gemmaout}.*.log.txt.gz"
	echo "assoc: ${assoc}"
	echo "pval: ${pval}"
	echo "log: ${log}"
	'''
}*/


/*process filecheck_rungemma {
input:
	val ready
//	path assoc name "${params.container_analysis_dir}/${params.kmer_type}kmer${params.kmer_length}_gemma/output/${params.output_prefix}_${params.kmer_type}${params.kmer_length}.*.assoc.txt.gz"
//	path pval name "${params.container_analysis_dir}/${params.kmer_type}kmer${params.kmer_length}_gemma/output/${params.output_prefix}_${params.kmer_type}${params.kmer_length}.*.pval.txt.gz"
//	path log name "${params.container_analysis_dir}/${params.kmer_type}kmer${params.kmer_length}_gemma/output/${params.output_prefix}_${params.kmer_type}${params.kmer_length}.*.log.txt.gz"
	path assoc
	path pval
	path log
output:
	val true, emit: done
shell:
if(!params.skip4 | !params.skip5)
	'''
	echo "Fchk 4: Running GEMMA"
	echo "assoc: !{assoc}"
	echo "pval: !{pval}"
	echo "log: !{log}"
	ln -sfr $(pwd) !{params.workdir}/filecheck_rungemma
	ln -sfr $(pwd)/.command.log !{params.logdir}/filecheck_rungemma.log
	'''
else
	'''
	echo "Skipping Fchk 4: Running GEMMA"
	'''
 }*/

// Step 5: Run contig alignment only (no merging)
//   p-fold parallelization
process kmercontigalign {
input:
	val ready
	val taskid
output:
	val true, emit: done
	val taskid
shell:
if(!params.skip5)
	'''
	echo "Step 5: Running contig alignment"
	ln -sfr $(pwd) !{params.workdir}/kmercontigalign.!{taskid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/kmercontigalign.!{taskid}.log
	!{params.container_cmd} !{params.container_script_dir}/kmercontigalignonly.Rscript \
		!{taskid} \
		!{params.n} \
		!{params.output_prefix} \
		!{params.container_analysis_dir} \
		!{params.container_id_file} \
		!{params.container_ref_fa} \
		!{params.container_ref_gb} \
		!{params.kmer_type} \
		!{params.kmer_length} \
		!{params.nucmerident} \
		!{params.kmerFilePrefix}.kmermerge.txt.gz \
		!{params.container_software_file}
	'''
else
	'''
	echo "Skipping Step 5: Running contig alignment"
	'''
}

// Step 5A: Merge contig alignments
//   p5-fold parallelization
process kmercontigalignmerge {
input:
	val ready
	val taskid
output:
	val true, emit: done
	val taskid
shell:
if(!params.skip5)
	'''
	echo "Step 5A: Merging contig alignments"
	ln -sfr $(pwd) !{params.workdir}/kmercontigalignmerge.!{taskid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/kmercontigalignmerge.!{taskid}.log
	!{params.container_cmd} !{params.container_script_dir}/kmercontigalignmerge.Rscript \
		!{taskid} \
		!{params.n} \
		!{params.p5} \
		!{params.output_prefix} \
		!{params.container_analysis_dir} \
		!{params.kmergenecombination} \
		!{params.kmer_type} \
		!{params.kmer_length} \
		!{params.container_ref_fa} \
		!{params.nucmerident} \
		!{params.container_software_file}
	'''
else
	'''
	echo "Skipping Step 5A: Merging contig alignments"
	'''
}

def outfiles_kmercontigalign() {
	kmercontigalignout = "${params.container_analysis_dir}/${params.kmer_type}kmer${params.kmer_length}_kmergenealign/${params.output_prefix}_${params.kmer_type}${params.kmer_length}"
	[geneIdNameLookup: "${kmercontigalignout}_${params.ref_name}_gene_id_name_lookup.txt",
	kmerListGeneIDs: "${kmercontigalignout}_${params.ref_name}_t${params.nucmerident}_*_nucmeralign_kmer_list_gene_IDs.txt.gz",
	kmerAlignMerge: "${params.container_analysis_dir}/${params.output_prefix}_${params.kmer_type}${params.kmer_length}.${params.ref_name}_t${params.nucmerident}.kmeralignmerge.txt.gz",
	kmerAlignMergeCount: "${params.container_analysis_dir}/${params.output_prefix}_${params.kmer_type}${params.kmer_length}.${params.ref_name}_t${params.nucmerident}.kmeralignmerge.count.txt.gz"]
	// Not clear from the documentation if the last two are needed downstream
}

/*process filecheck_kmercontigalign {
input:
	val ready
	path geneIdNameLookup
	path kmerListGeneIDs
output:
	val true, emit: done
shell:
if(!params.skip5 | !params.skip6)
	'''
	echo "Fchk 5: Running contig alignment"
	echo "geneIdNameLookup: !{geneIdNameLookup}"
	echo "kmerListGeneIDs: !{kmerListGeneIDs}"
	ln -sfr $(pwd) !{params.workdir}/filecheck_kmercontigalign
	ln -sfr $(pwd)/.command.log !{params.logdir}/filecheck_kmercontigalign.log
	'''
else
	'''
	echo "Skipping Fchk 5: Running contig alignment"
	'''
}*/

// Step 6: Plot figures using contig alignment positions
//   One core
process plotManhattan {
input:
	val ready_rungemma
	val ready_kmercontigalign
output:
	val true, emit: done
shell:
if(!params.skip6)
	'''
	echo "Step 6: Plotting figures using contig alignment positions"
	ln -sfr $(pwd) !{params.workdir}/plotManhattan
	ln -sfr $(pwd)/.command.log !{params.logdir}/plotManhattan.log
	!{params.container_cmd} !{params.container_script_dir}/plotManhattan.Rscript \
		!{params.output_prefix} \
		!{params.container_analysis_dir} \
		!{params.kmerFilePrefix} \
		!{params.container_ref_gb} \
		!{params.container_ref_fa} \
		!{params.gene_lookup_file} \
		!{params.container_id_file} \
		!{params.nucmerident} \
		!{params.min_count} \
		!{params.kmer_type} \
		!{params.kmer_length} \
		!{params.minor_allele_threshold} \
		!{params.container_software_file} \
		!{params.blastident} \
		!{params.ntopgenes}
	cp --remove-destination $(pwd)/.command.log !{params.logdir}/plotManhattan.log
	'''
	/* Temporarily removed since default values cannot be explicitly specified:\
	!{params.annotateGeneFile} \
	!{params.override_signif}*/
else
	'''
	echo "Skipping Step 6: Plotting figures using contig alignment positions"
	'''
}

// Step 5B: Run bowtie2 (nucleotide kmers only) -- fast!
//   One core
process runbowtie {
input:
	val ready
output:
	val true, emit: done
shell:
	'''
	echo "Step 5B: Running bowtie2 (nucleotide kmers only)"
	ln -sfr $(pwd) !{params.workdir}/runbowtie
	ln -sfr $(pwd)/.command.log !{params.logdir}/runbowtie.log
	!{params.container_cmd} !{params.container_script_dir}/runbowtie.Rscript \
		!{params.output_prefix} \
		!{params.container_analysis_dir} \
		!{params.kmerFilePrefix} \
		!{params.container_ref_fa} \
		!{params.kmer_type} \
		!{params.kmer_length} \
		!{params.container_software_file} \
		!{params.bowtie_parameters} \
		!{params.samtools_filter}
	'''
}

// Step 6B: Plot figures using bowtie2 mapping positions (nucleotide kmers only)
//   One core
process plotManhattanbowtie {
input:
	val ready
output:
	val true, emit: done
shell:
	'''
	echo "Step 6B: Plotting figures using bowtie2 mapping positions (nucleotide kmers only)"
	ln -sfr $(pwd) !{params.workdir}/plotManhattanbowtie
	ln -sfr $(pwd)/.command.log !{params.logdir}/plotManhattanbowtie.log
	!{params.container_cmd} !{params.container_script_dir}/plotManhattanbowtie.Rscript \
		!{params.output_prefix} \
		!{params.container_analysis_dir} \
		!{params.kmerFilePrefix} \
		!{params.container_ref_gb} \
		!{params.container_ref_fa} \
		!{params.container_id_file} \
		!{params.kmer_type} \
		!{params.kmer_length} \
		!{params.minor_allele_threshold} \
		!{params.samtools_filter} \
		!{params.container_software_file} \
		!{params.blastident} \
		!{params.ntopgenes}
	cp --remove-destination $(pwd)/.command.log !{params.logdir}/plotManhattanbowtie.log
	'''
	/* Temporarily removed since default values cannot be explicitly specified \
	!{params.annotateGeneFile} \
	!{params.override_signif} */
}

// Step 7: Generate HTML report
//   One process
process genReport {
input:
	val ready
output:
	val true, emit: done
shell:
if(!params.skip7)
	'''
	echo "Step 7: Generating HTML report"
	ln -sfr $(pwd) !{params.workdir}/genReport
	ln -sfr $(pwd)/.command.log !{params.logdir}/genReport.log
	!{params.container_cmd} !{params.container_script_dir}/gen-report.Rscript \
		!{params.output_prefix} \
		!{params.kmer_type} \
		!{params.kmer_length} \
		!{params.ref_name} \
		!{params.container_ref_gb} \
		!{params.minor_allele_threshold} \
		!{params.nucmerident} \
		!{params.min_count} \
		!{params.ntopgenes} \
		!{params.container_script_dir} \
		!{params.container_analysis_dir} \
		!{params.container_logdir}
	'''
else
	'''
	echo "Skipping Step 7: Generating HTML report"
	'''
}

// Step 7B: Generate HTML gene report
//   Linear parallelization
process genGeneReport {
input:
	val ready
	val hitnum
output:
	val true, emit: done
shell:
if(!params.skip7)
	'''
	echo "Step 7B: Generating HTML gene report"
	ln -sfr $(pwd) !{params.workdir}/genGeneReport.!{hitnum}
	ln -sfr $(pwd)/.command.log !{params.logdir}/genGeneReport.!{hitnum}.log
	!{params.container_cmd} !{params.container_script_dir}/gen-gene-report.Rscript \
		!{hitnum} \
		!{params.output_prefix} \
		!{params.kmer_type} \
		!{params.kmer_length} \
		!{params.ref_name} \
		!{params.container_ref_gb} \
		!{params.minor_allele_threshold} \
		!{params.nucmerident} \
		!{params.min_count} \
		!{params.container_script_dir} \
		!{params.container_analysis_dir} \
		!{params.container_logdir}
	'''
else
	'''
	echo "Skipping Step 7B: Generating HTML gene report"
	'''
}

// Step 7B: Generate HTML protein report
//   Linear parallelization
process genProteinReport {
input:
	val ready
	val hitnum
output:
	val true, emit: done
shell:
if(!params.skip7)
	'''
	echo "Step 7B: Generating HTML protein report"
	ln -sfr $(pwd) !{params.workdir}/genProteinReport.!{hitnum}
	ln -sfr $(pwd)/.command.log !{params.logdir}/genProteinReport.!{hitnum}.log
	!{params.container_cmd} !{params.container_script_dir}/gen-protein-report.Rscript \
		!{hitnum} \
		!{params.output_prefix} \
		!{params.kmer_type} \
		!{params.kmer_length} \
		!{params.ref_name} \
		!{params.container_ref_gb} \
		!{params.minor_allele_threshold} \
		!{params.nucmerident} \
		!{params.min_count} \
		!{params.container_script_dir} \
		!{params.container_analysis_dir} \
		!{params.container_logdir}
	'''
else
	'''
	echo "Skipping Step 7B: Generating HTML protein report"
	'''
}

// Step 7C: Generate HTML unmapped reads report
//   One process
process genUnmappedReport {
input:
	val ready
output:
	val true, emit: done
shell:
if(!params.skip7)
	'''
	echo "Step 7C: Generating HTML unmapped report"
	ln -sfr $(pwd) !{params.workdir}/genUnmappedReport
	ln -sfr $(pwd)/.command.log !{params.logdir}/genUnmappedReport.log
	!{params.container_cmd} !{params.container_script_dir}/gen-unmapped-report.Rscript \
		!{params.output_prefix} \
		!{params.kmer_type} \
		!{params.kmer_length} \
		!{params.ref_name} \
		!{params.container_ref_gb} \
		!{params.minor_allele_threshold} \
		!{params.nucmerident} \
		!{params.min_count} \
		!{params.container_script_dir} \
		!{params.container_analysis_dir} \
		!{params.container_logdir}
	'''
else
	'''
	echo "Skipping Step 7C: Generating HTML unmapped report"
	'''
}

// Print explicitly specified and default parameters.
// Default parameters can be overriden by specifying them in nextflow.config
// Unspecified parameters with no default will throw an error here.
println 'Parameters'
// Output files (no defaults)
println 'Output files'
println 'base_dir:                ' + params.base_dir
println 'output_prefix:           ' + params.output_prefix
println 'analysis_dir:            ' + params.analysis_dir
println ''
// Analysis options
println 'Analysis options'
println 'kmer_type:               ' + params.kmer_type		// No default
println 'kmer_length:             ' + params.kmer_length	// No default
params.ntopgenes = 20
println 'ntopgenes:               ' + params.ntopgenes
params.minor_allele_threshold = 0.01
println 'minor_allele_threshold:  ' + params.minor_allele_threshold
params.min_count = 1
println 'min_count:               ' + params.min_count
params.nucmerident = 90
println 'nucmerident:             ' + params.nucmerident
params.bowtie_parameters = "--very-sensitive"
println 'bowtie_parameters:       ' + params.bowtie_parameters
params.samtools_filter = 10
println 'samtools_filter:         ' + params.samtools_filter
params.blastident = 70
println 'blastident:              ' + params.blastident
//Neither yet implemented because of problem explicitly specifying NULL annotateGeneFile:
//params.annotateGeneFile = "NULL"
//println 'annotateGeneFile:        ' + params.annotateGeneFile
//params.override_signif = "FALSE"
//println 'override_signif:         ' + params.override_signif
println ''
// Input files
println 'Input files'
println 'id_file:                 ' + params.id_file			// No default
params.covariate_file = ""
println 'covariate_file:          ' + params.covariate_file
println ''
// Species-specific reference genome FASTA and genbank files (no defaults)
println 'Species-specific reference genome FASTA and genbank files'
println 'ref_fa:                  ' + params.ref_fa
println 'ref_gb:                  ' + params.ref_gb
println ''
// Deployment
deployment()
println 'Deployment'
println 'maxp:                    ' + params.maxp
println 'container_type:          ' + params.container_type
println 'container_file:          ' + params.container_file
println 'container_cmd:           ' + params.container_cmd
println 'container_args:          ' + params.container_args
println 'container_mount:         ' + params.container_mount
println 'software_file:           ' + params.software_file
println ''
// Workflow parameters
println 'Workflow parameters'
params.skip1 = false
println 'skip1:                   ' + params.skip1
params.skip2 = false
println 'skip2:                   ' + params.skip2
params.skip3 = false
println 'skip3:                   ' + params.skip3
params.skip4 = false
println 'skip4:                   ' + params.skip4
params.skip5 = false
println 'skip5:                   ' + params.skip5
params.skip6 = false
println 'skip6:                   ' + params.skip6
params.skip7 = false
println 'skip7:                   ' + params.skip7
println ''
// Implied parameters constructed from explicit parameters. Can be overridden by specifying them in nextflow.config
println 'Implied parameters'
params.container_script_dir = read_container_script_dir()
println 'container_script_dir:    ' + params.container_script_dir
params.analysis_file = params.analysis_dir + "/" + params.output_prefix + "_" + params.kmer_type + params.kmer_length + ".analysis_file.txt"
println 'analysis_file:           ' + params.analysis_file
params.container_analysis_file = params.container_analysis_dir + "/" + params.output_prefix + "_" + params.kmer_type + params.kmer_length + ".analysis_file.txt"
println 'container_analysis_file: ' + params.container_analysis_file
params.kmerFilePrefix = params.container_analysis_dir + "/" + params.output_prefix + "_" + params.kmer_type + params.kmer_length
println 'kmerFilePrefix:          ' + params.kmerFilePrefix
read_ref_name()
println 'ref_name:                ' + params.ref_name
params.gene_lookup_file = params.container_analysis_dir + "/" + params.kmer_type + "kmer" + params.kmer_length + "_kmergenealign/" + params.output_prefix + "_" + params.kmer_type + params.kmer_length + "_" + params.ref_name + "_gene_id_name_lookup.txt"
println 'gene_lookup_file:        ' + params.gene_lookup_file
params.kmergenecombination = params.container_analysis_dir + "/" + params.kmer_type + "kmer" + params.kmer_length + "_kmergenealign/" + params.output_prefix + "_" + params.kmer_type + params.kmer_length + "_" + params.ref_name + "_kmergenecombination_filepaths.txt"
println 'kmergenecombination:     ' + params.kmergenecombination
params.logdir = params.analysis_dir + "/log." + params.output_prefix + "_" + params.kmer_type + params.kmer_length
println 'logdir:                  ' + params.logdir
params.container_logdir = params.container_analysis_dir + "/log." + params.output_prefix + "_" + params.kmer_type + params.kmer_length
println 'container_software_file: ' + params.container_software_file
println 'container_analysis_dir:  ' + params.container_analysis_dir
println 'container_id_file:       ' + params.container_id_file
println 'container_covariate_file:' + params.container_covariate_file
println 'container_ref_fa:        ' + params.container_ref_fa
println 'container_ref_gb:        ' + params.container_ref_gb
println 'container_logdir:        ' + params.container_logdir
params.workdir = params.analysis_dir + "/work." + params.output_prefix + "_" + params.kmer_type + params.kmer_length
println 'workdir:                 ' + params.workdir
assert !binding.hasVariable('params.n')  // Do not allow override !
params.n = get_n()
println 'n:                       ' + params.n
params.p = (int)Math.max(1, Math.min(Math.ceil(params.n/2), params.maxp))
println 'p:                       ' + params.p
params.p5 = (int)Math.max(1, Math.min(Math.ceil(params.n/5), params.maxp))
println 'p5:                      ' + params.p5

workflow {
	// Prelim: Checking parameters
	check_params()
	
	// Prelim: Creating analysis directories
	create_analysis_dir(check_params.out.done)
	
	// Prelim: Defining analysis file
	create_analysis_file(create_analysis_dir.out.done)

	// Read main input file
	id_list = read_id_file()

	if(params.kmer_type.toString().toLowerCase()=="nucleotide") {
		// Step 1: Counting kmers
		// n-fold parallelization
		countkmers(create_analysis_file.out.done, Channel.of(1..params.n))

		// Step 2: Creating unique kmer list
		// Parallel pyramid (p-fold)
		createfullkmerlist(countkmers.out.done.collect(), Channel.of(1..params.p))

		// Step 3: Creating kmer presence/absence patterns and kinship matrix
		// Parallel pyramid (maxp-fold)
		stringlist2patternandkinship(createfullkmerlist.out.done.collect(), Channel.of(1..params.maxp))

		// Step 4: Running GEMMA
		// maxp-fold parallelization
		rungemma(stringlist2patternandkinship.out.done.collect(), Channel.of(1..params.maxp))

		// Step 5: Running contig alignment (no merging)
		// n-fold parallelization
		// Could branch from step 2 (not 4)
		//kmercontigalign(createfullkmerlist.out.done.collect(), Channel.of(1..params.n))
		kmercontigalign(rungemma.out.done.collect(), Channel.of(1..params.n))

		// Step 5A: Merge contig alignments
		// p5-fold parallelization
		kmercontigalignmerge(kmercontigalign.out.done.collect(), Channel.of(1..params.p5))

		// Step 6: Plotting figures using contig alignment positions
		// One core
		plotManhattan(rungemma.out.done.collect(), kmercontigalignmerge.out.done.collect())

		// Step 5B: Running bowtie2 (nucleotide kmers only)
		// One core
		//runbowtie(plotManhattan.out.done)

		// Step 6B: Plotting figures using bowtie2 mapping positions (nucleotide kmers only)
		// One core
		//plotManhattanbowtie(runbowtie.out.done)
		
		// Step 7: Generate HTML report
		//   One process
		genReport(plotManhattan.out.done)

		// Step 7B: Generate HTML gene report
		//   Linear parallelization
		genGeneReport(genReport.out.done, Channel.of(1..params.ntopgenes))

		// Step 7C: Generate HTML unmapped reads report
		//   One process
		genUnmappedReport(genGeneReport.out.done.collect())

	} else if(params.kmer_type.toString().toLowerCase()=="protein") {
		// Step 1: Counting kmers
		// n-fold parallelization
		countkmers(create_analysis_file.out.done, Channel.of(1..params.n))

		// Step 2: Creating unique kmer list
		// Parallel pyramid (p-fold)
		createfullkmerlist(countkmers.out.done.collect(), Channel.of(1..params.p))

		// Step 3: Creating kmer presence/absence patterns and kinship matrix
		// Parallel pyramid (maxp-fold)
		stringlist2patternandkinship(createfullkmerlist.out.done.collect(), Channel.of(1..params.maxp))

		// Step 4: Running GEMMA
		// maxp-fold parallelization
		rungemma(stringlist2patternandkinship.out.done.collect(), Channel.of(1..params.maxp))

		// Step 5: Running contig alignment (no merging)
		// n-fold parallelization
		// Could branch from step 2 (not 4)
		//kmercontigalign(createfullkmerlist.out.done.collect(), Channel.of(1..params.n))
		kmercontigalign(rungemma.out.done.collect(), Channel.of(1..params.n))

		// Step 5A: Merge contig alignments
		// p5-fold parallelization
		kmercontigalignmerge(kmercontigalign.out.done.collect(), Channel.of(1..params.p5))

		// Step 6: Plotting figures using contig alignment positions
		// One core
		plotManhattan(rungemma.out.done.collect(), kmercontigalignmerge.out.done.collect())
		
		// Step 7: Generate HTML report
		//   One process
		genReport(plotManhattan.out.done)

		// Step 7B: Generate HTML protein report
		//   Linear parallelization
		genProteinReport(genReport.out.done, Channel.of(1..params.ntopgenes))

		// Step 7C: Generate HTML unmapped reads report
		//   One process
		genUnmappedReport(genProteinReport.out.done.collect())

	} else {
		throw new Exception("params.kmer_type must be nucleotide or protein")
	}

	println 'Processes registered'
}
