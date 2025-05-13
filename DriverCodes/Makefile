# Usage: make <file> (without .c)

QCC = qcc
CFLAGS = -O2 -Wall -autolink -lm

# Detect platform
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    # macOS uses framework
    OPENGLIBS = -lfb_tiny -framework OpenGL
else
    # Linux uses standard GL linking
    OPENGLIBS = -lfb_tiny -lGL
endif

LDFLAGS = -L$(BASILISK)/gl -lglutils $(OPENGLIBS)

# Default rule: takes 'make target' and builds 'target' from 'target.c'
%: %.c
	$(QCC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f *.o *.gfs *.mp4 *.ppm *.png *~
