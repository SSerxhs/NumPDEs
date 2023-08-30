run: main.cpp
	python3 data_gen.py
	g++ -std=c++17 -O2 main.cpp -o main.exe
	time ./main.exe
	python3 plot.py