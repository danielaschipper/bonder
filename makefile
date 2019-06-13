

bonder : Analize.cpp fill.cpp iface.cpp LL.cpp LLi.cpp main.cpp readwfn.cpp stdafx.cpp output.cpp readwfx.cpp
	g++ -o bond Analize.cpp fill.cpp iface.cpp LL.cpp LLi.cpp main.cpp readwfn.cpp stdafx.cpp output.cpp readwfx.cpp -I .   -std=c++11 -fopenmp  -pthread -Ofast  -fsanitize=undefined
