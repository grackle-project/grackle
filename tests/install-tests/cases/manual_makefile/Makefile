# Before Grackle 3.2, CFLAGS needed to include -DCONFIG_BFLOAT_4 or
# -DCONFIG_BFLOAT_8 based on how grackle was compiled

GRACKLE_INSTALL = /home/gr-user/local
CFLAGS = -I${GRACKLE_INSTALL}/include
LDFLAGS = -L${GRACKLE_INSTALL}/lib -lgrackle

# flags unrelated to Grackle
# -> these are ONLY needed because c_example c_example.c makes
#    direct use of the standard C library's math functions
UNAME := $(shell uname)
ifneq ($(UNAME), Darwin)
  OTHER_LDFLAGS=-lm
endif

example_app:
	$(CC) $(CFLAGS) -c c_example.c -o c_example.o
	$(CC) c_example.o $(LDFLAGS) $(OTHER_LDFLAGS) -o example_app
