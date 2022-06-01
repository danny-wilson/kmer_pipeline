BINARIES=kmerlist2pattern pattern2kinship patterncounts patternmerge sort_strings stringlist2count stringlist2pattern

move : build
	cd C++/read_dsk && mv -t ../.. $(BINARIES)

build :
	cd C++/read_dsk && $(MAKE)
