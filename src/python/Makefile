build:
	DYLD_LIBRARY_PATH=$(shell pwd)/../clib/.libs LD_LIBRARY_PATH=$(shell pwd)/../clib/.libs LIBRARY_PATH=$(shell pwd)/../clib/.libs CPLUS_INCLUDE_PATH=$(shell pwd)/../clib/ C_INCLUDE_PATH=$(shell pwd)/../clib/ python setup.py build_ext -i
test:
	PYTHONPATH=$(shell pwd) DYLD_LIBRARY_PATH=$(shell pwd)/../clib/.libs LD_LIBRARY_PATH=$(shell pwd)/../clib/.libs LIBRARY_PATH=$(shell pwd)/../clib/.libs CPLUS_INCLUDE_PATH=$(shell pwd)/../clib/ C_INCLUDE_PATH=$(shell pwd)/../clib/ py.test tests
clean:
	rm -f pygrackle/*.so pygrackle/*.c pygrackle/*.pyc pygrackle/utilities/*.pyc
	rm -rf pygrackle.egg-info
	rm -rf build

.PHONY: build test
