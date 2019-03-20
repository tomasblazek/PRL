filename=xblaze31

${filename}.pdf:${filename}.tex
	latex ${filename}.tex
	dvips ${filename}.dvi
	ps2pdf ${filename}.ps

clean:
	rm -rf *.aux *.log *.dvi *.ps
