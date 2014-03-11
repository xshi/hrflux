#BINDIR=./bin
#LIBDIR=./lib

GCC_FLAGS=-Wall -g
ROOT_FLAGS=`root-config --cflags --libs`
BOOST_FLAGS=-I/usr/include/boost141/ -L/usr/lib64/boost141 -lboost_filesystem 

OS := $(shell uname)
ifeq ($(OS),Darwin)
	BOOST_FLAGS=-I/opt/local/include/ -L/opt/local/lib -lboost_system-mt  -lboost_filesystem-mt 
endif

SUMMARY_FLAGS=$(GCC_FLAGS) $(ROOT_FLAGS) 
FLUX_FLAGS=$(GCC_FLAGS) $(ROOT_FLAGS) $(BOOST_FLAGS) 
FITSPOT_FLAGS=$(GCC_FLAGS) $(ROOT_FLAGS) $(BOOST_FLAGS) 
DRAW_FLAGS=$(GCC_FLAGS) $(ROOT_FLAGS) 


all: summary fitspot flux draw draw_effvsdc draw_effvsdc_flux \
draw_effvsdc_flux_slices draw_effvsflux_slices draw_beamposition 
	@echo "Full build successful."


summary: summary.cc 
	g++ $(SUMMARY_FLAGS) $< -o $@ 

fitspot: fitspot.cc 
	g++ $(FITSPOT_FLAGS) $< -o $@

fitspot.o: fitspot.cc 
	g++ -c -D__LIBS__ $(FITSPOT_FLAGS) $< -o $@

flux: flux.cc fitspot.o
	g++ $(FLUX_FLAGS) $^ -o $@

draw_beamposition: draw_beamposition.cc fitspot.o
	g++ $(FITSPOT_FLAGS) $^ -o $@

draw: draw.cc 
	g++ $(DRAW_FLAGS) $< -o $@ 

%: %.cc 
	g++ $(DRAW_FLAGS) $< -o $@ 

clean:
	rm -rf summary fitspot flux draw draw_effvsdc draw_effvsdc_flux draw_effvsdc_flux_slices draw_effvsflux_slices *.o *~ *.orig  *.dSYM

