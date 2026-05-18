# -> it is the caller's responsibility to appropriately adjust the
#    PKG_CONFIG_PATH environment variable so that pkg-config can find
#    grackle.pc
# -> it is also the the caller's responsibility to setup LD_LIBRARY_PATH
#    appropriately if they want to use Grackle as a shared library.
#    (Alternative extra logic can be added to add -rpath to the linker
#    flags to accomplish the same thing)

CFLAGS = `pkg-config --cflags grackle`
LDFLAGS = `pkg-config --libs grackle`

# flags unrelated to Grackle
UNAME := $(shell uname)
ifneq ($(UNAME), Darwin)
  OTHER_LDFLAGS=-lm
endif

example_app:
	$(CC) $(CFLAGS) -c c_example.c -o c_example.o
	$(CC) c_example.o $(LDFLAGS) $(OTHER_LDFLAGS) -o example_app 
