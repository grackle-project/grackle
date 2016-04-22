build:
	$(MAKE) -C clib
build-python:
	cd python && python setup.py build_ext -i
clean:
	$(MAKE) -C clib  clean
	$(MAKE) -C python clean
install: build
	$(MAKE) -C clib  install
test: build build-python
	$(MAKE) -C python test
