all : clean data report methods second tutorial ITHIM-discussion
tex : methods second ITHIM-discussion

report: R/ITHIM.Rmd R/ITHIM.R
	cd R; R --vanilla -e 'source("ITHIM.R")'

tutorial: R/tutorial.Rmd R/tutorial.R
	cd R; R --vanilla -e 'source("tutorial.R")'

methods: methodsPaper/manuscript.tex tex/ITHIM.bib ./methodsPaper/tex/abstract.tex ./methodsPaper/tex/conclusion.tex ./methodsPaper/tex/discussion.tex ./methodsPaper/tex/introduction.tex ./methodsPaper/tex/materialsAndMethods.tex ./methodsPaper/tex/results.tex ./methodsPaper/tex/supportingInformation.tex ./methodsPaper/figures/fig1.pdf
	cd methodsPaper; make methodsPaper.pdf

second: secondPaper/manuscript.tex tex/ITHIM.bib ./secondPaper/tex/abstract.tex ./secondPaper/tex/discussion.tex ./secondPaper/tex/introduction.tex ./secondPaper/tex/methods.tex ./secondPaper/tex/publicHealthImplications.tex ./secondPaper/tex/results.tex
	cd secondPaper; make secondPaper.pdf

data:
	wget http://nhts.ornl.gov/2009/download/Ascii.zip
	mkdir -p ./data/NHTS
	unzip Ascii.zip -d ./data/NHTS/
	mv -v Ascii.zip ./data/NHTS/
	wget http://www.bls.gov/tus/special.requests/atuscps_2015.zip
	mkdir -p ./data/ATUS/
	unzip  atuscps_2015.zip -d ./data/ATUS/
	mv -v atuscps_2015.zip ./data/ATUS
	wget http://www.bls.gov/tus/special.requests/atussum_2015.zip
	mkdir -p ./data/ATUS/
	unzip  atussum_2015.zip -d ./data/ATUS/
	mv -v atussum_2015.zip ./data/ATUS


ITHIM-discussion: ./slides/ITHIM-discussion/ITHIM-discussion.tex
	cd ./slides/ITHIM-discussion; make ITHIM-discussion.pdf

webpage: html/index.txt
	perl ~/lib/Markdown_1.0.1/Markdown.pl ./html/index.txt > ./html/index.html
	perl ~/lib/Markdown_1.0.1/Markdown.pl ./html/tableau/tableau.txt > ./html/tableau/tableau.html

clean:
	rm -rf ./data/ATUS/
	rm -rf ./data/NHTS/
	rm -rf ./R/figure/
	mkdir -p ./R/figure/
	mkdir -p ./R/data/
	mkdir -p ./data/ATUS/
	mkdir -p ./data/NHTS/

.PHONY: data

.PHONY2: clean
