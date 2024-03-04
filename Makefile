all: build/MOM6

build/MOM6: build/Makefile
	make -C build MOM6

build/Makefile: build/config.status build/Makefile.in
	cd build && ./config.status

build/config.status: build/configure ac/deps/lib/libFMS.a
	cd build \
	&& PATH="${PATH}:$(dir $(abspath ac/makedep))" \
	  ./configure -n --with-framework=fms2

build/Makefile.in: ac/Makefile.in| build
	cp ac/Makefile.in build/

build/configure: build/configure.ac build/m4
	autoreconf build

build/configure.ac: ac/configure.ac | build
	cp ac/configure.ac build/

build/m4: ac/m4 | build
	cp -r ac/m4 build/

build:
	mkdir -p build

# FMS? (should move deps/ to build/)
ac/deps/lib/libFMS.a:
	make -C ac/deps/ -j

# Cleanup
clean:
	rm -rf build
	make -C ac/deps clean
