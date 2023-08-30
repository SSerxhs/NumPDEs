story:
	make run
	make report
	rm main.exe
	rm *.jpg
	rm *.txt
run:
	g++ -O3 -std=c++17 main.cpp -o main.exe
	./main.exe
	python3 plot.py
report:
	xelatex report.tex
	xelatex report.tex
	rm *.log
	rm *.aux