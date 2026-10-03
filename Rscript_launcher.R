# Rscript_launcher.R: run an R script so that errors can be debugged.
#   Rscript --vanilla Rscript_launcher.R SCRIPT.R [arguments...]
# The script is read with its source references kept, and must define main(args).
# On an error this prints the message, the failing call, and the file:line of each call
# leading to it, closes any graphics devices, and exits with status 1 (so Nextflow stops).
# Warnings are printed as they happen.
options(keep.source = TRUE, warn = 1)
launcher_args = commandArgs(trailingOnly = TRUE)
if(length(launcher_args) < 1) {
	cat("Usage: Rscript --vanilla Rscript_launcher.R SCRIPT.R [arguments...]\n", file = stderr())
	quit(status = 2, save = "no")
}
script = launcher_args[1]
stamp = file.path(dirname(normalizePath(script, mustWork = FALSE)), "STAGED_SHA")
cat("Running ", normalizePath(script, mustWork = FALSE), " (", if(file.exists(stamp)) readLines(stamp, n = 1, warn = FALSE) else "not staged",
	") with ", R.version.string, "\n", sep = "")

launcher_on_error = function(e) {
	calls = sys.calls()
	cat("Error: ", conditionMessage(e), "\n", sep = "", file = stderr())
	cl = conditionCall(e)
	if(!is.null(cl)) cat("In: ", paste(deparse(cl, nlines = 1), collapse = ""), "\n", sep = "", file = stderr())
	cat("Traceback (innermost last):\n", file = stderr())
	skip = "^(withCallingHandlers|launcher_on_error|source|eval|doTryCatch|tryCatch|tryCatchList|tryCatchOne|\\.handleSimpleError|h\\(simpleError|stop\\(|withVisible)"
	for(cl in calls) {
		fn = paste(deparse(cl, nlines = 1), collapse = "")
		sr = attr(cl, "srcref")
		if(is.null(sr) && grepl(skip, fn)) next
		loc = if(!is.null(sr)) paste0(basename(getSrcFilename(sr)), ":", sr[1]) else "(R internal)"
		cat("  ", loc, "  ", substr(fn, 1, 160), "\n", sep = "", file = stderr())
	}
	graphics.off()
	quit(status = 1, save = "no")
}

withCallingHandlers({
	if(!file.exists(script)) stop("Script not found: ", script)
	source(script, keep.source = TRUE)
	if(!exists("main", mode = "function")) stop("Script ", script, " does not define main(args)")
	main(launcher_args[-1])
}, error = launcher_on_error)
