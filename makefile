matrices = $(wildcard paar_header/*.h)
targets = $(matrices:.h=_check)
slpprograms = slp_heuristic LinOpt

all: slp_heuristic LinOpt $(targets)

$(targets): %_check: paar.cpp  # $*.h
	g++ $< -o $@_alg1 -O3 -Wall -std=c++11 -mpopcnt -D PAAR1 -include stdint.h -include $*.h
	g++ $< -o $@_alg2 -O3 -Wall -std=c++11 -mpopcnt -D PAAR2 -include stdint.h -include $*.h

$(slpprograms): %: %.cpp
	g++ -o $@ $< -std=c++11 -O3 -Wall

clean:
	${RM} paar_header/*_check_alg1 paar_header/*_check_alg2
	${RM} $(slpprograms)

.PHONY: clean
