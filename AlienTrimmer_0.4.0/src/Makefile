
GCJ=gcj
GCJFLAGS=-fsource=1.6 -march=native -O3 -minline-all-stringops -fomit-frame-pointer -momit-leaf-frame-pointer -fstrict-aliasing -fno-store-check -fno-bounds-check -Wall
EXEC=AlienTrimmer

AlienTrimmer: AlienTrimmer.java
	$(GCJ) $(GCJFLAGS) --main=AlienTrimmer AlienTrimmer.java -o $(EXEC)

