build:
	$(MAKE) -C clib
build-python:
	$(MAKE) -C python build
clean:
	$(MAKE) -C clib  clean
	$(MAKE) -C python clean
install: build
	$(MAKE) -C clib  install
test: build build-python
	$(MAKE) -C python test
