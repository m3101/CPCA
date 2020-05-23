help:
	@echo "Usage:"
	@echo "make lib		Downloads the required libraries"
	@echo "make tests		Builds the test program (requires a 32bit ARGB file named satie.bmp at the project root (preferrably a picture of Erik Satie))"
	@echo "make builds		Builds a shared library"
lib:
	@mkdir lib
	@cd lib;\
		git clone https://github.com/m3101/bmpreader.git;\
		cd bmpreader;make builds;\
		cp ./build/libbmpreader.so ../libbmpreader.so;\
		cp ./src/bmpreader.h ../bmpreader.h;\
		cd ..;\
		rm -rf bmpreader
tests:
	@make lib
	@if [ -d build ]; then rm -rf build;fi
	@mkdir build
	gcc -L./lib/ src/cpca.c src/linal.c test/gauss.c test/gfx.c -lX11 -lm -lbmpreader -o build/gauss -Wall -Werror -g
	gcc -L./lib/ src/cpca.c src/linal.c test/singval.c test/gfx.c -lX11 -lm -lbmpreader -o build/svd -Wall -Werror -g
builds:
	@make lib
	@if [ -d build ]; then rm -rf build;fi
	@mkdir build
	gcc -o build/cpca.o -c src/cpca.c -Wall -Werror -fpic
	gcc -o build/linal.o -c src/linal.c -Wall -Werror -fpic
	gcc -shared -o build/libcpca.so build/cpca.o build/linal.o -lm
	@cp build/libcpca.so lib/libcpca.so -f
	@rm build/*.o
	gcc -L./lib/ test/test.c -Wl,-rpath=./lib test/gfx.c -lX11 -lm -lbmpreader -lcpca -o build/test -Wall -Werror -g