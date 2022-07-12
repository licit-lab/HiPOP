.PHONY: python cpp

cpp:
	cd cpp/build && make -j install
	

python: cpp
	python -m pip install ./python

