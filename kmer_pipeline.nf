#!/usr/bin/env nextflow
// This version is for deployment on a range of executors.
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.Files
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import groovy.json.JsonOutput
import groovy.json.JsonSlurper

// The parameters the user set (config files and command line), before any default is assigned
USER_KEYS = new TreeSet(params.keySet())

// Real path of p, following symlinks; the part that does not exist yet (e.g. analysis_dir) is kept as written
def resolvePath(p) {
	def path = Paths.get(p).toAbsolutePath().normalize()
	def tail = []
	while(path != null && !Files.exists(path)) {
		tail.add(0, path.getFileName().toString())
		path = path.getParent()
	}
	if(path == null) return Paths.get(p).toAbsolutePath().normalize()
	def real = path.toRealPath()
	tail.each { real = real.resolve(it) }
	real
}

def user2containerPath(base_dir, user_path, container_base_dir) {
	// Throws an error if user_path is not in the subdirectory tree of base_dir. A path beneath base_dir as written
	// is used as written (a symlink there may point elsewhere, bound into the container by container_args).
	// Otherwise both are resolved first, so that two spellings of the same place (for example through a
	// symlinked directory) are treated alike.
	def written_base = base_dir.toAbsolutePath().normalize()
	def written_user = Paths.get(user_path).toAbsolutePath().normalize()
	if(written_user.startsWith(written_base))
		return Paths.get(container_base_dir, written_base.relativize(written_user).toString()).toString()
	def real_base = resolvePath(base_dir.toString())
	def real_user = resolvePath(user_path)
	if(!real_user.startsWith(real_base)) {
		println "Error converting from user_path to container_path"
		println "Every input and output must lie beneath base_dir, which is the only directory mounted in the container"
		println "base_dir:  " + base_dir + (real_base.toString() != base_dir.toString() ? " (" + real_base + ")" : "")
		println "user_path: " + user_path + (real_user.toString() != user_path ? " (" + real_user + ")" : "")
		throw new Exception("${user_path} is not beneath base_dir ${base_dir}")
	}
	Paths.get(container_base_dir, real_base.relativize(real_user).toString()).toString()
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

	// Analysis, work and log directories (created by deployment_write(), after the checks)
	params.workdir = params.analysis_dir + "/work." + params.output_prefix + "_" + params.kmer_type + params.kmer_length
	params.logdir = params.analysis_dir + "/log." + params.output_prefix + "_" + params.kmer_type + params.kmer_length

	// Read main input file
	id_list = read_id_file()

	// Convert io files from user file system to container file system
	base_dir = Paths.get(params.base_dir)
	params.container_analysis_dir = user2containerPath(base_dir, params.analysis_dir, params.container_mount)
	// Each param is assigned once: Nextflow ignores later assignments to a param
	if(params.pheno_file!="" && !Files.exists(Paths.get(params.pheno_file))) throw new Exception("pheno_file ${params.pheno_file} does not exist")
	params.container_pheno_file = params.pheno_file=="" ? "" : user2containerPath(base_dir, params.pheno_file, params.container_mount)
	if(params.precomputed_dir!="" && !Files.isDirectory(Paths.get(params.precomputed_dir))) throw new Exception("precomputed_dir ${params.precomputed_dir} is not a folder")
	params.container_precomputed_dir = params.precomputed_dir=="" ? "" : user2containerPath(base_dir, params.precomputed_dir, params.container_mount)
	if(params.covariate_file=="") {
		params.container_covariate_file = ""
	} else {
		if(!Files.exists(Paths.get(params.covariate_file))) throw new Exception("covariate_file ${params.covariate_file} does not exist")
		params.container_covariate_file = user2containerPath(base_dir, params.covariate_file, params.container_mount)
	}
	if(params.containsKey('annotateGeneFile')) {
		if(!Files.exists(Paths.get(params.annotateGeneFile))) throw new Exception("annotateGeneFile ${params.annotateGeneFile} does not exist")
		params.container_annotateGeneFile = user2containerPath(base_dir, params.annotateGeneFile, params.container_mount)
	} else {
		params.container_annotateGeneFile = ""
	}
	params.container_ref_fa = user2containerPath(base_dir, params.ref_fa, params.container_mount)
	params.container_ref_gb = user2containerPath(base_dir, params.ref_gb, params.container_mount)

	// Convert paths in id_list from user file system to container file system
	id_list["containerpaths"] = id_list["paths"].collect( { user_dir -> user2containerPath(base_dir, user_dir, params.container_mount) })

	// Write the new id_list
	userfs_container_id_file = params.analysis_dir + "/" + params.output_prefix + "_" + params.kmer_type + params.kmer_length + ".container_id_file.txt"
	params.container_id_file = user2containerPath(base_dir, userfs_container_id_file, params.container_mount)
	params.container_user_id_file = user2containerPath(base_dir, params.id_file, params.container_mount)
	DEPLOYMENT_FILES = [id_list: id_list, container_id_file: userfs_container_id_file]

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
// *** EXCEPT for the default value which is interpreted as in the container file system ***
// The files and folders deployment() would have written, written once the checks have passed
def deployment_write() {
	Files.createDirectories(Paths.get(params.analysis_dir))
	Files.createDirectories(Paths.get(params.workdir))
	Files.createDirectories(Paths.get(params.logdir))
	write_container_id_file(DEPLOYMENT_FILES.id_list, DEPLOYMENT_FILES.container_id_file)
	create_analysis_file()
}

// D5: every assembly in id_file exists, can be read, and starts as a FASTA file (gzipped or not).
// Checked here, on the user's file system, where the paths in id_file are.
def check_assemblies(id_list) {
	def problems = []
	id_list['paths'].eachWithIndex { path, k ->
		def f = new File(path)
		if(!f.exists()) problems << "${id_list['id'][k]}: ${path} does not exist"
		else if(!f.canRead()) problems << "${id_list['id'][k]}: ${path} cannot be read"
		else if(f.length() == 0) problems << "${id_list['id'][k]}: ${path} is empty"
		else {
			try {
				def stream = f.newInputStream()
				if(path.endsWith(".gz")) stream = new java.util.zip.GZIPInputStream(stream)
				def first = -1
				try {
					while((first = stream.read()) != -1 && Character.isWhitespace(first as char)) {}
				} finally { stream.close() }
				if(first != ('>' as char) as int) problems << "${id_list['id'][k]}: ${path} is not a FASTA file (it should start with >)"
			} catch(Exception e) {
				problems << "${id_list['id'][k]}: ${path} cannot be read (${e.message})"
			}
		}
	}
	if(problems.size() > 20) problems = problems.take(20) + ["and ${problems.size() - 20} more assembly problems"]
	return problems.collect { "assembly ${it}".toString() }
}

// A true/false parameter: a Boolean, or the text true or false in any case
def parse_bool(name, value) {
	if(value instanceof Boolean) return value
	def v = value.toString().trim().toLowerCase()
	if(v == "true") return true
	if(v == "false") return false
	throw new Exception("${name} must be true or false, not '${value}'")
}

// Run preflight.py in the container: parameter, overwrite, reuse and -resume checks (and, with
// overwrite = true, removal of the outputs of an earlier run). Stops the workflow on any error.
// With "--finish", "finished" or "failed", records the end of the run in the run manifest.
def preflight(List extra = []) {
	if(!extra) {  // assemblies first: preflight.py may delete files (overwrite = true)
		def problems = check_assemblies(DEPLOYMENT_FILES.id_list)
		if(problems) {
			problems.each { println "Error: ${it}" }
			throw new Exception("stopped before running anything: see the errors above")
		}
	}
	def run_steps = (1..7).findAll { !SKIP[it] }.join(",")
	def result_params = [kmer_type: params.kmer_type, kmer_length: params.kmer_length.toString(),
		kmer_min_count: params.kmer_min_count.toString(), plot_min_genomes: params.plot_min_genomes.toString(),
		minor_allele_threshold: params.minor_allele_threshold.toString(), nucmerident: params.nucmerident.toString(),
		ntopgenes: params.ntopgenes.toString(), blastident: params.blastident.toString(), maxp: params.maxp.toString(),
		min_contig_length: (params.containsKey('min_contig_length') ? params.min_contig_length : 0).toString(),
		output_prefix: params.output_prefix, run_steps: run_steps, precomputed_dir: params.container_precomputed_dir,
		precomputed_prefix: params.precomputed_prefix]
	def input_files = [id_file: params.container_user_id_file, covariate_file: params.container_covariate_file,
		pheno_file: params.container_pheno_file,
		ref_fa: params.container_ref_fa, ref_gb: params.container_ref_gb]
	def cmd = params.container_cmd.tokenize() + ["${params.container_script_dir}/preflight.py".toString(),
		"--analysis-dir", params.container_analysis_dir, "--output-prefix", params.output_prefix,
		"--kmer-type", params.kmer_type, "--kmer-length", params.kmer_length.toString(),
		"--session-id", workflow.sessionId.toString(), "--run-steps", run_steps,
		"--overwrite", OVERWRITE.toString(), "--resume", workflow.resume.toString(),
		"--pid", ProcessHandle.current().pid().toString(), "--user-params", USER_KEYS.join(","),
		"--params-json", JsonOutput.toJson(result_params), "--input-files", JsonOutput.toJson(input_files),
		"--id-file", params.container_user_id_file, "--covariate-file", params.container_covariate_file,
		"--pheno-file", params.container_pheno_file, "--precomputed-dir", params.container_precomputed_dir,
		"--precomputed-prefix", params.precomputed_prefix, "--ref-name", params.ref_name.toString()] + extra
	def proc = cmd.execute()
	def sout = new StringBuilder(), serr = new StringBuilder()
	proc.waitForProcessOutput(sout, serr)
	if(serr.toString().trim()) println serr.toString().trim()
	if(proc.exitValue() != 0) throw new Exception("preflight.py failed (exit status ${proc.exitValue()})")
	def result = new JsonSlurper().parseText(sout.toString().trim().readLines().last())
	result.warnings.each { println "Warning: ${it}" }
	if(result.deleted) {
		println "Deleted ${result.deleted.size()} files of an earlier run (overwrite = true), for example:"
		result.deleted.take(20).each { println "  ${it}" }
	}
	if(result.errors) {
		result.errors.each { println "Error: ${it}" }
		throw new Exception("stopped before running anything: see the errors above")
	}
}

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
	cmd = "${params.container_cmd} ${params.container_script_dir}/get_ref_name.py --fasta-file ${params.container_ref_fa}"
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

// Read main input file
def read_id_file() {
	def id_file = new File(params.id_file)
	assert id_file.readLines().head().split('\t')*.toLowerCase() == ['id','paths','pheno']
	def rows = id_file.readLines().tail()*.split('\t')
	id = rows.collect( { row -> row[0] })
	paths = rows.collect( { row -> row[1] })
	pheno = rows.collect( { row -> row[2] })
	ret = [id: id, paths: paths, pheno: pheno]
}

// Write a copy of the main input file with paths relative to the container file system
def write_container_id_file(id_list, filename) {
	nlines = id_list['id'].size()
	assert id_list['paths'].size()==nlines
	assert id_list['pheno'].size()==nlines
	assert id_list['containerpaths'].size()==nlines
	if(nlines>0) {
		def id_file = new PrintWriter(new FileWriter(filename))
		id_file.printf("%s\t%s\t%s\n", "id", "paths", "pheno")
		for(int i=0; i<nlines; i++) {
			id_file.printf("%s\t%s\t%s\n", id_list["id"][i], id_list["containerpaths"][i], id_list["pheno"][i])
		}
		id_file.close()
	}
}

def create_analysis_file() {
	file = new PrintWriter(new FileWriter(params.analysis_file))
	file.printf("%s\t%s\n", "kmertype", "kmerlength")
	file.printf("%s\t%s\n", params.kmer_type, params.kmer_length)
	file.close()
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
	val sampid
output:
	val true, emit: done
	val sampid
shell:
if(!SKIP[1])
	'''
	echo "Step 1: Counting kmers"
	echo "sampid: !{sampid}"
	ln -sfr $(pwd) !{params.workdir}/countkmers.!{sampid} 2>/dev/null || ln -sf $(pwd) !{params.workdir}/countkmers.!{sampid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/countkmers.!{sampid}.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/countkmers.!{sampid}.log
	!{params.container_cmd} !{params.container_script_dir}/countkmers.py \
		--task-id !{sampid} \
		--id-file !{params.container_id_file} \
		--analysis-dir !{params.container_analysis_dir} \
		--output-prefix !{params.output_prefix} \
		--software-file !{params.container_software_file} \
		--analyses-list !{params.container_analysis_file} \
		!{MINLEN_ARG}
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
if(!SKIP[1] | !SKIP[2])
	'''
	echo "Fchk 1: Counting kmers"
	echo "max_sampid: !{max_sampid}"
	echo "kmercounts: !{kmercounts}"
	echo "kmertotals: !{kmertotals}"
	echo "kmers_filepaths: !{kmers_filepaths}"
	ln -sfr $(pwd) !{params.workdir}/filecheck_countkmers 2>/dev/null || ln -sf $(pwd) !{params.workdir}/filecheck_countkmers
	ln -sfr $(pwd)/.command.log !{params.logdir}/filecheck_countkmers.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/filecheck_countkmers.log
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
if(!SKIP[1] | !SKIP[2])
	'''
	echo "Fchk 1: Counting kmers (proteins)"
	echo "max_sampid: !{max_sampid}"
	echo "kmercounts: !{kmercounts}"
	echo "kmertotals: !{kmertotals}"
	echo "reading_frames: !{reading_frames}"
	echo "kmers_filepaths: !{kmers_filepaths}"
	ln -sfr $(pwd) !{params.workdir}/filecheck_countkmers_protein 2>/dev/null || ln -sf $(pwd) !{params.workdir}/filecheck_countkmers_protein
	ln -sfr $(pwd)/.command.log !{params.logdir}/filecheck_countkmers_protein.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/filecheck_countkmers_protein.log
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
if(!SKIP[2])
	'''
	echo "Step 2: Creating unique kmer list"
	echo "taskid: !{taskid}"
	ln -sfr $(pwd) !{params.workdir}/createfullkmerlist.!{taskid} 2>/dev/null || ln -sf $(pwd) !{params.workdir}/createfullkmerlist.!{taskid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/createfullkmerlist.!{taskid}.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/createfullkmerlist.!{taskid}.log
	!{params.container_cmd} !{params.container_script_dir}/createfullkmerlist.py \
		--task-id !{taskid} \
		--n !{params.n} \
		--p !{params.p} \
		--output-prefix !{params.output_prefix} \
		--analysis-dir !{params.container_analysis_dir} \
		--id-file !{params.container_id_file} \
		--kmer-type !{params.kmer_type} \
		--kmer-length !{params.kmer_length} \
		--software-file !{params.container_software_file} \
		--merge-wait-minutes !{params.merge_wait_minutes}
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
if(!SKIP[2] | !SKIP[3])
	'''
	echo "Fchk 2: Creating unique kmer list"
	echo "kmermerge: !{kmermerge}"
	ln -sfr $(pwd) !{params.workdir}/filecheck_createfullkmerlist 2>/dev/null || ln -sf $(pwd) !{params.workdir}/filecheck_createfullkmerlist
	ln -sfr $(pwd)/.command.log !{params.logdir}/filecheck_createfullkmerlist.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/filecheck_createfullkmerlist.log
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
if(!SKIP[3])
	'''
	echo "Step 3: Creating kmer presence/absence patterns and kinship matrix"
	ln -sfr $(pwd) !{params.workdir}/stringlist2patternandkinship.!{taskid} 2>/dev/null || ln -sf $(pwd) !{params.workdir}/stringlist2patternandkinship.!{taskid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/stringlist2patternandkinship.!{taskid}.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/stringlist2patternandkinship.!{taskid}.log
	!{params.container_cmd} !{params.container_script_dir}/stringlist2patternandkinship.py \
		--task-id !{taskid} \
		--p !{params.p} \
		--id-file !{params.container_id_file} \
		--fullkmerlistfile !{params.kmerFilePrefix}.kmermerge.txt.gz \
		--kmercountslistfile !{params.kmerFilePrefix}_kmers_filepaths.txt \
		--analysis-dir !{params.container_analysis_dir} \
		--output-prefix !{params.output_prefix} \
		--kmertype !{params.kmer_type} \
		--software-file !{params.container_software_file} \
		--kmer-length !{params.kmer_length} \
		--kmer-min-count !{params.kmer_min_count} \
		--merge-wait-minutes !{params.merge_wait_minutes}
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
if(!SKIP[3] | !SKIP[4])
	'''
	echo "Fchk 3: Creating kmer presence/absence patterns and kinship matrix"
	echo "patternKey: !{patternKey}"
	echo "patternIndex: !{patternIndex}"
	echo "presenceCount: !{presenceCount}"
	echo "patternKeySize: !{patternKeySize}"
	echo "kinship: !{kinship}"
	echo "kinshipWeight: !{kinshipWeight}"
	ln -sfr $(pwd) !{params.workdir}/filecheck_stringlist2patternandkinship 2>/dev/null || ln -sf $(pwd) !{params.workdir}/filecheck_stringlist2patternandkinship
	ln -sfr $(pwd)/.command.log !{params.logdir}/filecheck_stringlist2patternandkinship.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/filecheck_stringlist2patternandkinship.log
	'''
else
	'''
	echo "Skipping Fchk 3: Creating kmer presence/absence patterns and kinship matrix"
	'''
}

// Step 4: Run GEMMA
//   p-fold parallelization
// Step 4 preparation, once before the GEMMA tasks: the genomes GEMMA analyses (finite phenotype,
// complete covariates), GEMMA's phenotype file, the patterns' presence counts among those genomes
// and the kinship matrix decompressed once. Not cached: it must match the current phenotype.
process prepareGemma {
cache false
input:
	val ready
output:
	val true, emit: done
shell:
if(!SKIP[4])
	'''
	echo "Step 4: Preparing GEMMA"
	ln -sfr $(pwd) !{params.workdir}/prepareGemma 2>/dev/null || ln -sf $(pwd) !{params.workdir}/prepareGemma
	ln -sfr $(pwd)/.command.log !{params.logdir}/prepareGemma.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/prepareGemma.log
	!{params.container_cmd} !{params.container_script_dir}/prepare_gemma.py \
		--kmerfile-prefix !{params.kmerFilePrefix} \
		--id-file !{params.container_id_file} \
		--analysis-dir !{params.container_analysis_dir} \
		--output-prefix !{params.output_prefix} \
		--kmer-type !{params.kmer_type} \
		--kmer-length !{params.kmer_length} \
		!{COVARIATE_ARG} !{PHENO_ARG}
	rm -f !{params.logdir}/prepareGemma.log && cp $(pwd)/.command.log !{params.logdir}/prepareGemma.log
	'''
else
	'''
	echo "Skipping Step 4: Preparing GEMMA"
	'''
}

// After the last GEMMA task: remove the decompressed kinship matrix
process cleanupGemma {
cache false
input:
	val ready
output:
	val true, emit: done
shell:
if(!SKIP[4])
	'''
	!{params.container_cmd} !{params.container_script_dir}/prepare_gemma.py --cleanup \
		--kmerfile-prefix !{params.kmerFilePrefix} \
		--id-file !{params.container_id_file} \
		--analysis-dir !{params.container_analysis_dir} \
		--output-prefix !{params.output_prefix} \
		--kmer-type !{params.kmer_type} \
		--kmer-length !{params.kmer_length}
	'''
else
	'''
	echo "Skipping Step 4: cleaning up"
	'''
}

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
if(!SKIP[4] && params.container_covariate_file=="")
	'''
	echo "Step 4: Running GEMMA"
	ln -sfr $(pwd) !{params.workdir}/rungemma.!{taskid} 2>/dev/null || ln -sf $(pwd) !{params.workdir}/rungemma.!{taskid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/rungemma.!{taskid}.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/rungemma.!{taskid}.log
	!{params.container_cmd} !{params.container_script_dir}/rungemma.py \
		--prepared \
		--task-id !{taskid} \
		--p !{params.p} \
		--kmerfile-prefix !{params.kmerFilePrefix} \
		--id-file !{params.container_id_file} \
		--output-prefix !{params.output_prefix} \
		--analysis-dir !{params.container_analysis_dir} \
		--kmertype !{params.kmer_type} \
		--kmer-length !{params.kmer_length} \
		--software-file !{params.container_software_file}
	rm -f !{params.logdir}/rungemma.!{taskid}.log && cp $(pwd)/.command.log !{params.logdir}/rungemma.!{taskid}.log
	'''
else if(!SKIP[4] && params.container_covariate_file!="")
	'''
	echo "Step 4: Running GEMMA"
	ln -sfr $(pwd) !{params.workdir}/rungemma.!{taskid} 2>/dev/null || ln -sf $(pwd) !{params.workdir}/rungemma.!{taskid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/rungemma.!{taskid}.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/rungemma.!{taskid}.log
	!{params.container_cmd} !{params.container_script_dir}/rungemma.py \
		--prepared \
		--task-id !{taskid} \
		--p !{params.p} \
		--kmerfile-prefix !{params.kmerFilePrefix} \
		--id-file !{params.container_id_file} \
		--output-prefix !{params.output_prefix} \
		--analysis-dir !{params.container_analysis_dir} \
		--kmertype !{params.kmer_type} \
		--kmer-length !{params.kmer_length} \
		--software-file !{params.container_software_file} \
		--covariate-file !{params.container_covariate_file}
	rm -f !{params.logdir}/rungemma.!{taskid}.log && cp $(pwd)/.command.log !{params.logdir}/rungemma.!{taskid}.log
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
if(!SKIP[4] | !SKIP[5])
	'''
	echo "Fchk 4: Running GEMMA"
	echo "assoc: !{assoc}"
	echo "pval: !{pval}"
	echo "log: !{log}"
	ln -sfr $(pwd) !{params.workdir}/filecheck_rungemma 2>/dev/null || ln -sf $(pwd) !{params.workdir}/filecheck_rungemma
	ln -sfr $(pwd)/.command.log !{params.logdir}/filecheck_rungemma.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/filecheck_rungemma.log
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
if(!SKIP[5])
	'''
	echo "Step 5: Running contig alignment"
	ln -sfr $(pwd) !{params.workdir}/kmercontigalign.!{taskid} 2>/dev/null || ln -sf $(pwd) !{params.workdir}/kmercontigalign.!{taskid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/kmercontigalign.!{taskid}.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/kmercontigalign.!{taskid}.log
	!{params.container_cmd} !{params.container_script_dir}/kmercontigalignonly.py \
		--task-id !{taskid} \
		--n !{params.n} \
		--output-prefix !{params.output_prefix} \
		--output-dir !{params.container_analysis_dir} \
		--id-file !{params.container_id_file} \
		--ref-fa !{params.container_ref_fa} \
		--ref-gb !{params.container_ref_gb} \
		--kmer-type !{params.kmer_type} \
		--kmer-length !{params.kmer_length} \
		--nucmerident !{params.nucmerident} \
		--kmerseqfile !{params.kmerFilePrefix}.kmermerge.txt.gz \
		--software-file !{params.container_software_file}
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
if(!SKIP[5])
	'''
	echo "Step 5A: Merging contig alignments"
	ln -sfr $(pwd) !{params.workdir}/kmercontigalignmerge.!{taskid} 2>/dev/null || ln -sf $(pwd) !{params.workdir}/kmercontigalignmerge.!{taskid}
	ln -sfr $(pwd)/.command.log !{params.logdir}/kmercontigalignmerge.!{taskid}.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/kmercontigalignmerge.!{taskid}.log
	!{params.container_cmd} !{params.container_script_dir}/kmercontigalignmerge.py \
		--task-id !{taskid} \
		--n !{params.n} \
		--p !{params.p5} \
		--output-prefix !{params.output_prefix} \
		--analysis-dir !{params.container_analysis_dir} \
		--input-files !{params.kmergenecombination} \
		--kmer-type !{params.kmer_type} \
		--kmer-length !{params.kmer_length} \
		--ref-fa !{params.container_ref_fa} \
		--nucmerident !{params.nucmerident} \
		--software-file !{params.container_software_file} \
		--merge-wait-minutes !{params.merge_wait_minutes}
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
if(!SKIP[5] | !SKIP[6])
	'''
	echo "Fchk 5: Running contig alignment"
	echo "geneIdNameLookup: !{geneIdNameLookup}"
	echo "kmerListGeneIDs: !{kmerListGeneIDs}"
	ln -sfr $(pwd) !{params.workdir}/filecheck_kmercontigalign 2>/dev/null || ln -sf $(pwd) !{params.workdir}/filecheck_kmercontigalign
	ln -sfr $(pwd)/.command.log !{params.logdir}/filecheck_kmercontigalign.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/filecheck_kmercontigalign.log
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
if(!SKIP[6])
	'''
	echo "Step 6: Plotting figures using contig alignment positions"
	ln -sfr $(pwd) !{params.workdir}/plotManhattan 2>/dev/null || ln -sf $(pwd) !{params.workdir}/plotManhattan
	ln -sfr $(pwd)/.command.log !{params.logdir}/plotManhattan.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/plotManhattan.log
	!{params.container_cmd} !{params.container_script_dir}/plotManhattan.py \
		--output-prefix !{params.output_prefix} \
		--analysis-dir !{params.container_analysis_dir} \
		--kmerfile-prefix !{params.kmerFilePrefix} \
		--ref-gb !{params.container_ref_gb} \
		--ref-fa !{params.container_ref_fa} \
		--gene-lookup-file !{params.gene_lookup_file} \
		--id-file !{params.container_id_file} \
		--nucmerident !{params.nucmerident} \
		--plot-min-genomes !{params.plot_min_genomes} \
		--kmer-type !{params.kmer_type} \
		--kmer-length !{params.kmer_length} \
		--minor-allele-threshold !{params.minor_allele_threshold} \
		--software-file !{params.container_software_file} \
		--blastident !{params.blastident} \
		--ngenes !{params.ntopgenes} \
		!{COVARIATE_ARG} !{ANNOTATE_ARG}
	rm -f !{params.logdir}/plotManhattan.log && cp $(pwd)/.command.log !{params.logdir}/plotManhattan.log
	'''
else
	'''
	echo "Skipping Step 6: Plotting figures using contig alignment positions"
	'''
}

// Step 6, figures: draw the figures in R from the data plotManhattan wrote
//   One process; runs with Step 6 (skip6)
process plotFigures {
input:
	val ready
output:
	val true, emit: done
shell:
if(!SKIP[6])
	'''
	echo "Step 6, figures: Drawing figures in R"
	# Recorded here so that changing a figure-only parameter reruns this task under -resume
	echo "Figure parameters: ntopgenes=!{params.ntopgenes} blastident=!{params.blastident} annotateGeneFile=!{params.containsKey('annotateGeneFile') ? params.annotateGeneFile : ''}"
	ln -sfr $(pwd) !{params.workdir}/plotFigures 2>/dev/null || ln -sf $(pwd) !{params.workdir}/plotFigures
	ln -sfr $(pwd)/.command.log !{params.logdir}/plotFigures.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/plotFigures.log
	!{params.container_cmd} Rscript --vanilla !{params.container_script_dir}/Rscript_launcher.R \
		!{params.container_script_dir}/plot_figures.R \
		--data-dir !{params.container_analysis_dir}/!{params.kmer_type}kmer!{params.kmer_length}_kmergenealign_figures/figure_data
	rm -f !{params.logdir}/plotFigures.log && cp $(pwd)/.command.log !{params.logdir}/plotFigures.log
	'''
else
	'''
	echo "Skipping Step 6, figures: Drawing figures in R"
	'''
}

// Step 6B, figures: draw the bowtie2-mapping figures in R (nucleotide kmers only)
process plotFiguresbowtie {
input:
	val ready
output:
	val true, emit: done
shell:
	'''
	echo "Step 6B, figures: Drawing bowtie2-mapping figures in R"
	ln -sfr $(pwd) !{params.workdir}/plotFiguresbowtie 2>/dev/null || ln -sf $(pwd) !{params.workdir}/plotFiguresbowtie
	ln -sfr $(pwd)/.command.log !{params.logdir}/plotFiguresbowtie.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/plotFiguresbowtie.log
	!{params.container_cmd} Rscript --vanilla !{params.container_script_dir}/Rscript_launcher.R \
		!{params.container_script_dir}/plot_figures.R \
		--data-dir !{params.container_analysis_dir}/!{params.kmer_type}kmer!{params.kmer_length}_bowtie2mapping_figures/figure_data
	rm -f !{params.logdir}/plotFiguresbowtie.log && cp $(pwd)/.command.log !{params.logdir}/plotFiguresbowtie.log
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
	ln -sfr $(pwd) !{params.workdir}/runbowtie 2>/dev/null || ln -sf $(pwd) !{params.workdir}/runbowtie
	ln -sfr $(pwd)/.command.log !{params.logdir}/runbowtie.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/runbowtie.log
	!{params.container_cmd} !{params.container_script_dir}/runbowtie.py \
		--output-prefix !{params.output_prefix} \
		--analysis-dir !{params.container_analysis_dir} \
		--kmerfile-prefix !{params.kmerFilePrefix} \
		--ref-fa !{params.container_ref_fa} \
		--kmer-type !{params.kmer_type} \
		--kmer-length !{params.kmer_length} \
		--software-file !{params.container_software_file} \
		--bowtie-parameters=!{params.bowtie_parameters} \
		--samtools-filter !{params.samtools_filter}
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
	ln -sfr $(pwd) !{params.workdir}/plotManhattanbowtie 2>/dev/null || ln -sf $(pwd) !{params.workdir}/plotManhattanbowtie
	ln -sfr $(pwd)/.command.log !{params.logdir}/plotManhattanbowtie.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/plotManhattanbowtie.log
	!{params.container_cmd} !{params.container_script_dir}/plotManhattanbowtie.py \
		--output-prefix !{params.output_prefix} \
		--analysis-dir !{params.container_analysis_dir} \
		--kmerfile-prefix !{params.kmerFilePrefix} \
		--ref-gb !{params.container_ref_gb} \
		--ref-fa !{params.container_ref_fa} \
		--id-file !{params.container_id_file} \
		--kmer-type !{params.kmer_type} \
		--kmer-length !{params.kmer_length} \
		--minor-allele-threshold !{params.minor_allele_threshold} \
		--samtools-filter !{params.samtools_filter} \
		--software-file !{params.container_software_file} \
		--blastident !{params.blastident} \
		--ngenes !{params.ntopgenes} \
		!{ANNOTATE_ARG}
	rm -f !{params.logdir}/plotManhattanbowtie.log && cp $(pwd)/.command.log !{params.logdir}/plotManhattanbowtie.log
	'''
}

// Step 7: Generate HTML report
//   One process
process genReport {
input:
	val ready
output:
	val true, emit: done
shell:
if(!SKIP[7])
	'''
	echo "Step 7: Generating HTML report"
	# Recorded here so that changing a figure-only parameter reruns this task under -resume
	echo "Figure parameters: ntopgenes=!{params.ntopgenes} blastident=!{params.blastident} annotateGeneFile=!{params.containsKey('annotateGeneFile') ? params.annotateGeneFile : ''}"
	ln -sfr $(pwd) !{params.workdir}/genReport 2>/dev/null || ln -sf $(pwd) !{params.workdir}/genReport
	ln -sfr $(pwd)/.command.log !{params.logdir}/genReport.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/genReport.log
	!{params.container_cmd} !{params.container_script_dir}/gen-report.py \
		--prefix !{params.output_prefix} \
		--anatype !{params.kmer_type} \
		--k !{params.kmer_length} \
		--refname !{params.ref_name} \
		--ref-gb !{params.container_ref_gb} \
		--maf !{params.minor_allele_threshold} \
		--alignident !{params.nucmerident} \
		--plot-min-genomes !{params.plot_min_genomes} \
		--ngenes !{params.ntopgenes} \
		--srcdir !{params.container_script_dir} \
		--outdir !{params.container_analysis_dir} \
		--logdir !{params.container_logdir}
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
if(!SKIP[7])
	'''
	echo "Step 7B: Generating HTML gene report"
	# Recorded here so that changing a figure-only parameter reruns this task under -resume
	echo "Figure parameters: ntopgenes=!{params.ntopgenes} blastident=!{params.blastident} annotateGeneFile=!{params.containsKey('annotateGeneFile') ? params.annotateGeneFile : ''}"
	ln -sfr $(pwd) !{params.workdir}/genGeneReport.!{hitnum} 2>/dev/null || ln -sf $(pwd) !{params.workdir}/genGeneReport.!{hitnum}
	ln -sfr $(pwd)/.command.log !{params.logdir}/genGeneReport.!{hitnum}.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/genGeneReport.!{hitnum}.log
	!{params.container_cmd} !{params.container_script_dir}/gen-gene-report.py \
		--hit-num !{hitnum} \
		--prefix !{params.output_prefix} \
		--anatype !{params.kmer_type} \
		--k !{params.kmer_length} \
		--refname !{params.ref_name} \
		--ref-gb !{params.container_ref_gb} \
		--maf !{params.minor_allele_threshold} \
		--alignident !{params.nucmerident} \
		--plot-min-genomes !{params.plot_min_genomes} \
		--srcdir !{params.container_script_dir} \
		--outdir !{params.container_analysis_dir} \
		--logdir !{params.container_logdir}
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
if(!SKIP[7])
	'''
	echo "Step 7B: Generating HTML protein report"
	# Recorded here so that changing a figure-only parameter reruns this task under -resume
	echo "Figure parameters: ntopgenes=!{params.ntopgenes} blastident=!{params.blastident} annotateGeneFile=!{params.containsKey('annotateGeneFile') ? params.annotateGeneFile : ''}"
	ln -sfr $(pwd) !{params.workdir}/genProteinReport.!{hitnum} 2>/dev/null || ln -sf $(pwd) !{params.workdir}/genProteinReport.!{hitnum}
	ln -sfr $(pwd)/.command.log !{params.logdir}/genProteinReport.!{hitnum}.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/genProteinReport.!{hitnum}.log
	!{params.container_cmd} !{params.container_script_dir}/gen-protein-report.py \
		--hit-num !{hitnum} \
		--prefix !{params.output_prefix} \
		--anatype !{params.kmer_type} \
		--k !{params.kmer_length} \
		--refname !{params.ref_name} \
		--ref-gb !{params.container_ref_gb} \
		--maf !{params.minor_allele_threshold} \
		--alignident !{params.nucmerident} \
		--plot-min-genomes !{params.plot_min_genomes} \
		--srcdir !{params.container_script_dir} \
		--outdir !{params.container_analysis_dir} \
		--logdir !{params.container_logdir}
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
if(!SKIP[7])
	'''
	echo "Step 7C: Generating HTML unmapped report"
	# Recorded here so that changing a figure-only parameter reruns this task under -resume
	echo "Figure parameters: ntopgenes=!{params.ntopgenes} blastident=!{params.blastident} annotateGeneFile=!{params.containsKey('annotateGeneFile') ? params.annotateGeneFile : ''}"
	ln -sfr $(pwd) !{params.workdir}/genUnmappedReport 2>/dev/null || ln -sf $(pwd) !{params.workdir}/genUnmappedReport
	ln -sfr $(pwd)/.command.log !{params.logdir}/genUnmappedReport.log 2>/dev/null || ln -sf $(pwd)/.command.log !{params.logdir}/genUnmappedReport.log
	!{params.container_cmd} !{params.container_script_dir}/gen-unmapped-report.py \
		--prefix !{params.output_prefix} \
		--anatype !{params.kmer_type} \
		--k !{params.kmer_length} \
		--refname !{params.ref_name} \
		--ref-gb !{params.container_ref_gb} \
		--maf !{params.minor_allele_threshold} \
		--alignident !{params.nucmerident} \
		--plot-min-genomes !{params.plot_min_genomes} \
		--srcdir !{params.container_script_dir} \
		--outdir !{params.container_analysis_dir} \
		--logdir !{params.container_logdir}
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
// D4: min_count is split into kmer_min_count (copies of a k-mer in a genome for it to count as
// present, step 3) and plot_min_genomes (genomes a k-mer/gene combination must be seen in to be
// plotted, steps 6-7). A legacy min_count sets both, as before. Checked before any default is set.
if(params.containsKey('min_count')) {
	if(params.containsKey('kmer_min_count') || params.containsKey('plot_min_genomes'))
		throw new Exception("min_count is replaced by kmer_min_count and plot_min_genomes: set only the new parameters")
	params.kmer_min_count = params.min_count
	params.plot_min_genomes = params.min_count
	println "Warning: min_count is replaced by kmer_min_count (copies of a k-mer in a genome for it to count as present) and plot_min_genomes (genomes a k-mer/gene combination must be seen in to be plotted); min_count = ${params.min_count} sets both" + (params.min_count.toString().toInteger() > 1 ? ", so a k-mer must also occur ${params.min_count} times in a genome to count as present" : "")
} else {
	params.kmer_min_count = 1
	params.plot_min_genomes = 1
}
println 'kmer_min_count:          ' + params.kmer_min_count
println 'min_contig_length:       ' + (params.containsKey('min_contig_length') ? params.min_contig_length : 0)
println 'plot_min_genomes:        ' + params.plot_min_genomes
params.nucmerident = 90
println 'nucmerident:             ' + params.nucmerident
params.bowtie_parameters = "--very-sensitive"
println 'bowtie_parameters:       ' + params.bowtie_parameters
params.samtools_filter = 10
println 'samtools_filter:         ' + params.samtools_filter
params.blastident = 70
println 'blastident:              ' + params.blastident
// N7: how long a merging task (steps 2, 3, 5A) waits for files written by other tasks
params.merge_wait_minutes = 100
println 'merge_wait_minutes:      ' + params.merge_wait_minutes
// annotateGeneFile has no default: set it to draw close-ups for chosen genes instead of the top ntopgenes
params.override_signif = "FALSE"
println 'override_signif:         ' + params.override_signif
if(params.containsKey('annotateGeneFile')) println 'annotateGeneFile:        ' + params.annotateGeneFile
println ''
// Input files
println 'Input files'
println 'id_file:                 ' + params.id_file			// No default
params.covariate_file = ""
// N5: phenotypes from another file, matched to id_file by ID; steps 1-3 and 5 from an earlier analysis
params.pheno_file = ""
println 'pheno_file:              ' + params.pheno_file
params.precomputed_dir = ""
println 'precomputed_dir:         ' + params.precomputed_dir
params.precomputed_prefix = params.output_prefix
println 'precomputed_prefix:      ' + params.precomputed_prefix
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
// skip1-7 and overwrite, as Booleans used everywhere below (Nextflow keeps a user's value as given,
// so a quoted "false" stays text: each is converted once here)
SKIP = [:]
(1..7).each { k -> SKIP[k] = parse_bool("skip${k}", params.containsKey("skip${k}".toString()) ? params["skip${k}".toString()] : false) }
(1..7).each { k -> println "skip${k}:                   " + SKIP[k] }
OVERWRITE = parse_bool("overwrite", params.containsKey('overwrite') ? params.overwrite : false)
println 'overwrite:               ' + OVERWRITE
// N5: with precomputed_dir, steps 1-3 and 5 come from it, so they are skipped
if(params.container_precomputed_dir) {
	[1, 2, 3, 5].each { k ->
		if(params.containsKey("skip${k}".toString()) && !SKIP[k])
			throw new Exception("skip${k} = false: with precomputed_dir, steps 1-3 and 5 are read from it, not run")
		SKIP[k] = true
	}
	println 'steps 1-3 and 5 skipped: read from precomputed_dir'
}
// The folder and prefix of the outputs of steps 1-3 and 5: this analysis's, or precomputed_dir's
INPUT_DIR = params.container_precomputed_dir ?: params.container_analysis_dir
INPUT_PREFIX = params.container_precomputed_dir ? params.precomputed_prefix : params.output_prefix
println ''
// Implied parameters constructed from explicit parameters. Can be overridden by specifying them in nextflow.config
println 'Implied parameters'
params.container_script_dir = read_container_script_dir()
println 'container_script_dir:    ' + params.container_script_dir
params.analysis_file = params.analysis_dir + "/" + params.output_prefix + "_" + params.kmer_type + params.kmer_length + ".analysis_file.txt"
println 'analysis_file:           ' + params.analysis_file
params.container_analysis_file = params.container_analysis_dir + "/" + params.output_prefix + "_" + params.kmer_type + params.kmer_length + ".analysis_file.txt"
println 'container_analysis_file: ' + params.container_analysis_file
params.kmerFilePrefix = INPUT_DIR + "/" + INPUT_PREFIX + "_" + params.kmer_type + params.kmer_length
println 'kmerFilePrefix:          ' + params.kmerFilePrefix
read_ref_name()
println 'ref_name:                ' + params.ref_name
params.gene_lookup_file = INPUT_DIR + "/" + params.kmer_type + "kmer" + params.kmer_length + "_kmergenealign/" + INPUT_PREFIX + "_" + params.kmer_type + params.kmer_length + "_" + params.ref_name + "_gene_id_name_lookup.txt"
println 'gene_lookup_file:        ' + params.gene_lookup_file
params.kmergenecombination = params.container_analysis_dir + "/" + params.kmer_type + "kmer" + params.kmer_length + "_kmergenealign/" + params.output_prefix + "_" + params.kmer_type + params.kmer_length + "_" + params.ref_name + "_kmergenecombination_filepaths.txt"
println 'kmergenecombination:     ' + params.kmergenecombination
println 'logdir:                  ' + params.logdir
params.container_logdir = params.container_analysis_dir + "/log." + params.output_prefix + "_" + params.kmer_type + params.kmer_length
println 'container_software_file: ' + params.container_software_file
println 'container_analysis_dir:  ' + params.container_analysis_dir
println 'container_id_file:       ' + params.container_id_file
println 'container_covariate_file:' + params.container_covariate_file
println 'container_ref_fa:        ' + params.container_ref_fa
println 'container_ref_gb:        ' + params.container_ref_gb
println 'container_logdir:        ' + params.container_logdir
println 'workdir:                 ' + params.workdir
// The covariate file option of the scripts that take one (empty without a covariate file)
COVARIATE_ARG = params.container_covariate_file ? "--covariate-file " + params.container_covariate_file : ""
// Optional: contigs shorter than min_contig_length are ignored when counting (default 0 keeps all; 10 x the k-mer length is a sensible choice)
MINLEN_ARG = params.containsKey('min_contig_length') ? "--min-contig-length " + params.min_contig_length : ""
// Optional: genes (or geneA:geneB intergenic regions) to draw close-ups for, instead of the top ntopgenes
ANNOTATE_ARG = params.container_annotateGeneFile ? "--annotate-gene-file " + params.container_annotateGeneFile + " --override-signif " + params.override_signif : ""
PHENO_ARG = params.container_pheno_file ? "--pheno-file " + params.container_pheno_file : ""
println ''
// Checks before anything is written; then the workflow's own files
preflight()
deployment_write()
println ''
assert !binding.hasVariable('params.n')  // Do not allow override !
params.n = get_n()
println 'n:                       ' + params.n
params.p = (int)Math.max(1, Math.min(Math.ceil(params.n/2), params.maxp))
println 'p:                       ' + params.p
params.p5 = (int)Math.max(1, Math.min(Math.ceil(params.n/5), params.maxp))
println 'p5:                      ' + params.p5

// Record in the run manifest whether the run finished
workflow.onComplete {
	try {
		preflight(["--finish", workflow.success ? "finished" : "failed"])
	} catch(Exception e) {
		println "Warning: could not record the end of the run in the run manifest: ${e.message}"
	}
}

// A process whose steps are all skipped is not submitted at all (each submission is a Slurm job that would only
// print "Skipping Step n"); downstream processes get an immediately available token in its place
def RUNS(List steps) { steps.any { !SKIP[it] } }

workflow {
	if(params.kmer_type.toString().toLowerCase()=="nucleotide") {
		// Step 1: Counting kmers
		// n-fold parallelization
		if(RUNS([1])) countkmers(Channel.of(1..params.n))

		// Step 2: Creating unique kmer list
		// Parallel pyramid (p-fold)
		if(RUNS([2])) createfullkmerlist((RUNS([1]) ? countkmers.out.done : Channel.value(true)).collect(), Channel.of(1..params.p))

		// Step 3: Creating kmer presence/absence patterns and kinship matrix
		// Parallel pyramid (maxp-fold)
		if(RUNS([3])) stringlist2patternandkinship((RUNS([2]) ? createfullkmerlist.out.done : Channel.value(true)).collect(), Channel.of(1..params.maxp))

		// Step 4: Running GEMMA
		// maxp-fold parallelization
		if(RUNS([4])) prepareGemma((RUNS([3]) ? stringlist2patternandkinship.out.done : Channel.value(true)).collect())
		if(RUNS([4])) rungemma((RUNS([4]) ? prepareGemma.out.done : Channel.value(true)).collect(), Channel.of(1..params.maxp))
		if(RUNS([4])) cleanupGemma((RUNS([4]) ? rungemma.out.done : Channel.value(true)).collect())

		// Step 5: Running contig alignment (no merging)
		// n-fold parallelization
		// Could branch from step 2 (not 4)
		//kmercontigalign(createfullkmerlist.out.done.collect(), Channel.of(1..params.n))
		if(RUNS([5])) kmercontigalign((RUNS([4]) ? rungemma.out.done : Channel.value(true)).collect(), Channel.of(1..params.n))

		// Step 5A: Merge contig alignments
		// p5-fold parallelization
		if(RUNS([5])) kmercontigalignmerge((RUNS([5]) ? kmercontigalign.out.done : Channel.value(true)).collect(), Channel.of(1..params.p5))

		// Step 6: Plotting figures using contig alignment positions
		// One core
		if(RUNS([6])) plotManhattan((RUNS([4]) ? cleanupGemma.out.done : Channel.value(true)).collect(), (RUNS([5]) ? kmercontigalignmerge.out.done : Channel.value(true)).collect())

		// Step 6, figures: drawn in R
		// One core
		if(RUNS([6])) plotFigures(plotManhattan.out.done)

		// Step 5B: Running bowtie2 (nucleotide kmers only)
		// One core
		//runbowtie(plotManhattan.out.done)

		// Step 6B: Plotting figures using bowtie2 mapping positions (nucleotide kmers only)
		// One core
		//plotManhattanbowtie(runbowtie.out.done)
		//plotFiguresbowtie(plotManhattanbowtie.out.done)
		
		// Step 7: Generate HTML report
		//   One process
		if(RUNS([7])) genReport((RUNS([6]) ? plotFigures.out.done : Channel.value(true)))

		// Step 7B: Generate HTML gene report
		//   Linear parallelization
		if(RUNS([7])) genGeneReport(genReport.out.done, Channel.of(1..params.ntopgenes))

		// Step 7C: Generate HTML unmapped reads report
		//   One process
		if(RUNS([7])) genUnmappedReport(genGeneReport.out.done.collect())

	} else if(params.kmer_type.toString().toLowerCase()=="protein") {
		// Step 1: Counting kmers
		// n-fold parallelization
		if(RUNS([1])) countkmers(Channel.of(1..params.n))

		// Step 2: Creating unique kmer list
		// Parallel pyramid (p-fold)
		if(RUNS([2])) createfullkmerlist((RUNS([1]) ? countkmers.out.done : Channel.value(true)).collect(), Channel.of(1..params.p))

		// Step 3: Creating kmer presence/absence patterns and kinship matrix
		// Parallel pyramid (maxp-fold)
		if(RUNS([3])) stringlist2patternandkinship((RUNS([2]) ? createfullkmerlist.out.done : Channel.value(true)).collect(), Channel.of(1..params.maxp))

		// Step 4: Running GEMMA
		// maxp-fold parallelization
		if(RUNS([4])) prepareGemma((RUNS([3]) ? stringlist2patternandkinship.out.done : Channel.value(true)).collect())
		if(RUNS([4])) rungemma((RUNS([4]) ? prepareGemma.out.done : Channel.value(true)).collect(), Channel.of(1..params.maxp))
		if(RUNS([4])) cleanupGemma((RUNS([4]) ? rungemma.out.done : Channel.value(true)).collect())

		// Step 5: Running contig alignment (no merging)
		// n-fold parallelization
		// Could branch from step 2 (not 4)
		//kmercontigalign(createfullkmerlist.out.done.collect(), Channel.of(1..params.n))
		if(RUNS([5])) kmercontigalign((RUNS([4]) ? rungemma.out.done : Channel.value(true)).collect(), Channel.of(1..params.n))

		// Step 5A: Merge contig alignments
		// p5-fold parallelization
		if(RUNS([5])) kmercontigalignmerge((RUNS([5]) ? kmercontigalign.out.done : Channel.value(true)).collect(), Channel.of(1..params.p5))

		// Step 6: Plotting figures using contig alignment positions
		// One core
		if(RUNS([6])) plotManhattan((RUNS([4]) ? cleanupGemma.out.done : Channel.value(true)).collect(), (RUNS([5]) ? kmercontigalignmerge.out.done : Channel.value(true)).collect())

		// Step 6, figures: drawn in R
		// One core
		if(RUNS([6])) plotFigures(plotManhattan.out.done)

		// Step 7: Generate HTML report
		//   One process
		if(RUNS([7])) genReport((RUNS([6]) ? plotFigures.out.done : Channel.value(true)))

		// Step 7B: Generate HTML protein report
		//   Linear parallelization
		if(RUNS([7])) genProteinReport(genReport.out.done, Channel.of(1..params.ntopgenes))

		// Step 7C: Generate HTML unmapped reads report
		//   One process
		if(RUNS([7])) genUnmappedReport(genProteinReport.out.done.collect())

	} else {
		throw new Exception("params.kmer_type must be nucleotide or protein")
	}

	println 'Processes registered'
}
