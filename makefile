boostH=/home/ds143/.local/lib/include
boostLib=/home/ds143/.local/lib/lib

bonder : Analize.cpp fill.cpp iface.cpp LL.cpp LLi.cpp main.cpp readwfn.cpp stdafx.cpp output.cpp readwfx.cpp
	g++ -o bond Analize.cpp fill.cpp iface.cpp LL.cpp LLi.cpp main.cpp readwfn.cpp stdafx.cpp output.cpp readwfx.cpp -I . -I $(boostH)   -std=c++11  $(boostLib)/libboost_thread.a $(boostLib)/libboost_system.a -pthread -Ofast  #-fsanitize=undefined
