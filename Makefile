%.html: %.rmd
	echo "rmarkdown::render(basename(\"$<\"))" | R --slave

report: report.html

run:
	R --slave < script.R
