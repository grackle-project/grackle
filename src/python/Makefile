test:
	PYTHONPATH=$(shell pwd) LD_LIBRARY_PATH=$(shell pwd)../clib/.libs py.test tests
clean:
	rm -f pygrackle/*.so pygrackle/*.c pygrackle/*.pyc pygrackle/utilities/*.pyc
	rm -rf pygrackle.egg_info
