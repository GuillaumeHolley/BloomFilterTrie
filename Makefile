WORKDIR = `pwd`

CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

INC = 
CFLAGS = -O3 -Wextra -Wall
RESINC = 
LIBDIR = 
LIB = 
LDFLAGS = 

INC_DEBUG = $(INC)
CFLAGS_DEBUG = $(CFLAGS) -Wall -g -std=c99
RESINC_DEBUG = $(RESINC)
RCFLAGS_DEBUG = $(RCFLAGS)
LIBDIR_DEBUG = $(LIBDIR)
LIB_DEBUG = $(LIB)
LDFLAGS_DEBUG = $(LDFLAGS) -lJudy
OBJDIR_DEBUG = obj/Debug
DEP_DEBUG = 
OUT_DEBUG = bin/Debug/bft

INC_RELEASE = $(INC)
CFLAGS_RELEASE = $(CFLAGS) -O3 -Wall -std=c99
RESINC_RELEASE = $(RESINC)
RCFLAGS_RELEASE = $(RCFLAGS)
LIBDIR_RELEASE = $(LIBDIR)
LIB_RELEASE = $(LIB)
LDFLAGS_RELEASE = $(LDFLAGS) -lJudy -ljemalloc
OBJDIR_RELEASE = obj/Release
DEP_RELEASE = 
OUT_RELEASE = bin/Release/bft

OBJ_DEBUG = $(OBJDIR_DEBUG)/src/insertNode.o $(OBJDIR_DEBUG)/src/CC.o $(OBJDIR_DEBUG)/src/UC.o $(OBJDIR_DEBUG)/src/UC_annotation.o $(OBJDIR_DEBUG)/src/annotation.o $(OBJDIR_DEBUG)/src/annotation_special_nodes.o $(OBJDIR_DEBUG)/src/branchingNode.o $(OBJDIR_DEBUG)/src/deleteColorsNode.o $(OBJDIR_DEBUG)/src/fasta.o $(OBJDIR_DEBUG)/src/getRSS.o $(OBJDIR_DEBUG)/src/log2.o $(OBJDIR_DEBUG)/src/main.o $(OBJDIR_DEBUG)/src/marking.o $(OBJDIR_DEBUG)/src/presenceNode.o $(OBJDIR_DEBUG)/src/printMemory.o $(OBJDIR_DEBUG)/src/quicksort.o $(OBJDIR_DEBUG)/src/replaceAnnotation.o $(OBJDIR_DEBUG)/src/retrieveAnnotation.o $(OBJDIR_DEBUG)/src/write_to_disk.o

OBJ_RELEASE = $(OBJDIR_RELEASE)/src/insertNode.o $(OBJDIR_RELEASE)/src/CC.o $(OBJDIR_RELEASE)/src/UC.o $(OBJDIR_RELEASE)/src/UC_annotation.o $(OBJDIR_RELEASE)/src/annotation.o $(OBJDIR_RELEASE)/src/annotation_special_nodes.o $(OBJDIR_RELEASE)/src/branchingNode.o $(OBJDIR_RELEASE)/src/deleteColorsNode.o $(OBJDIR_RELEASE)/src/fasta.o $(OBJDIR_RELEASE)/src/getRSS.o $(OBJDIR_RELEASE)/src/log2.o $(OBJDIR_RELEASE)/src/main.o $(OBJDIR_RELEASE)/src/marking.o $(OBJDIR_RELEASE)/src/presenceNode.o $(OBJDIR_RELEASE)/src/printMemory.o $(OBJDIR_RELEASE)/src/quicksort.o $(OBJDIR_RELEASE)/src/replaceAnnotation.o $(OBJDIR_RELEASE)/src/retrieveAnnotation.o $(OBJDIR_RELEASE)/src/write_to_disk.o

all: debug release

clean: clean_debug clean_release

before_debug: 
	test -d bin/Debug || mkdir -p bin/Debug
	test -d $(OBJDIR_DEBUG)/src || mkdir -p $(OBJDIR_DEBUG)/src

after_debug: 

debug: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG) $(DEP_DEBUG)
	$(LD) $(LIBDIR_DEBUG) -o $(OUT_DEBUG) $(OBJ_DEBUG)  $(LDFLAGS_DEBUG) $(LIB_DEBUG)

$(OBJDIR_DEBUG)/src/insertNode.o: src/insertNode.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/insertNode.c -o $(OBJDIR_DEBUG)/src/insertNode.o

$(OBJDIR_DEBUG)/src/CC.o: src/CC.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/CC.c -o $(OBJDIR_DEBUG)/src/CC.o

$(OBJDIR_DEBUG)/src/UC.o: src/UC.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/UC.c -o $(OBJDIR_DEBUG)/src/UC.o

$(OBJDIR_DEBUG)/src/UC_annotation.o: src/UC_annotation.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/UC_annotation.c -o $(OBJDIR_DEBUG)/src/UC_annotation.o

$(OBJDIR_DEBUG)/src/annotation.o: src/annotation.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/annotation.c -o $(OBJDIR_DEBUG)/src/annotation.o

$(OBJDIR_DEBUG)/src/annotation_special_nodes.o: src/annotation_special_nodes.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/annotation_special_nodes.c -o $(OBJDIR_DEBUG)/src/annotation_special_nodes.o

$(OBJDIR_DEBUG)/src/branchingNode.o: src/branchingNode.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/branchingNode.c -o $(OBJDIR_DEBUG)/src/branchingNode.o

$(OBJDIR_DEBUG)/src/deleteColorsNode.o: src/deleteColorsNode.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/deleteColorsNode.c -o $(OBJDIR_DEBUG)/src/deleteColorsNode.o

$(OBJDIR_DEBUG)/src/fasta.o: src/fasta.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/fasta.c -o $(OBJDIR_DEBUG)/src/fasta.o

$(OBJDIR_DEBUG)/src/getRSS.o: src/getRSS.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/getRSS.c -o $(OBJDIR_DEBUG)/src/getRSS.o

$(OBJDIR_DEBUG)/src/log2.o: src/log2.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/log2.c -o $(OBJDIR_DEBUG)/src/log2.o

$(OBJDIR_DEBUG)/src/main.o: src/main.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/main.c -o $(OBJDIR_DEBUG)/src/main.o

$(OBJDIR_DEBUG)/src/marking.o: src/marking.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/marking.c -o $(OBJDIR_DEBUG)/src/marking.o

$(OBJDIR_DEBUG)/src/presenceNode.o: src/presenceNode.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/presenceNode.c -o $(OBJDIR_DEBUG)/src/presenceNode.o

$(OBJDIR_DEBUG)/src/printMemory.o: src/printMemory.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/printMemory.c -o $(OBJDIR_DEBUG)/src/printMemory.o

$(OBJDIR_DEBUG)/src/quicksort.o: src/quicksort.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/quicksort.c -o $(OBJDIR_DEBUG)/src/quicksort.o

$(OBJDIR_DEBUG)/src/replaceAnnotation.o: src/replaceAnnotation.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/replaceAnnotation.c -o $(OBJDIR_DEBUG)/src/replaceAnnotation.o

$(OBJDIR_DEBUG)/src/retrieveAnnotation.o: src/retrieveAnnotation.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/retrieveAnnotation.c -o $(OBJDIR_DEBUG)/src/retrieveAnnotation.o

$(OBJDIR_DEBUG)/src/write_to_disk.o: src/write_to_disk.c
	$(CC) $(CFLAGS_DEBUG) $(INC_DEBUG) -c src/write_to_disk.c -o $(OBJDIR_DEBUG)/src/write_to_disk.o

clean_debug: 
	rm -f $(OBJ_DEBUG) $(OUT_DEBUG)
	rm -rf bin/Debug
	rm -rf $(OBJDIR_DEBUG)/src

before_release: 
	test -d bin/Release || mkdir -p bin/Release
	test -d $(OBJDIR_RELEASE)/src || mkdir -p $(OBJDIR_RELEASE)/src

after_release: 

release: before_release out_release after_release

out_release: before_release $(OBJ_RELEASE) $(DEP_RELEASE)
	$(LD) $(LIBDIR_RELEASE) -o $(OUT_RELEASE) $(OBJ_RELEASE)  $(LDFLAGS_RELEASE) $(LIB_RELEASE)

$(OBJDIR_RELEASE)/src/insertNode.o: src/insertNode.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/insertNode.c -o $(OBJDIR_RELEASE)/src/insertNode.o

$(OBJDIR_RELEASE)/src/CC.o: src/CC.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/CC.c -o $(OBJDIR_RELEASE)/src/CC.o

$(OBJDIR_RELEASE)/src/UC.o: src/UC.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/UC.c -o $(OBJDIR_RELEASE)/src/UC.o

$(OBJDIR_RELEASE)/src/UC_annotation.o: src/UC_annotation.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/UC_annotation.c -o $(OBJDIR_RELEASE)/src/UC_annotation.o

$(OBJDIR_RELEASE)/src/annotation.o: src/annotation.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/annotation.c -o $(OBJDIR_RELEASE)/src/annotation.o

$(OBJDIR_RELEASE)/src/annotation_special_nodes.o: src/annotation_special_nodes.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/annotation_special_nodes.c -o $(OBJDIR_RELEASE)/src/annotation_special_nodes.o

$(OBJDIR_RELEASE)/src/branchingNode.o: src/branchingNode.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/branchingNode.c -o $(OBJDIR_RELEASE)/src/branchingNode.o

$(OBJDIR_RELEASE)/src/deleteColorsNode.o: src/deleteColorsNode.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/deleteColorsNode.c -o $(OBJDIR_RELEASE)/src/deleteColorsNode.o

$(OBJDIR_RELEASE)/src/fasta.o: src/fasta.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/fasta.c -o $(OBJDIR_RELEASE)/src/fasta.o

$(OBJDIR_RELEASE)/src/getRSS.o: src/getRSS.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/getRSS.c -o $(OBJDIR_RELEASE)/src/getRSS.o

$(OBJDIR_RELEASE)/src/log2.o: src/log2.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/log2.c -o $(OBJDIR_RELEASE)/src/log2.o

$(OBJDIR_RELEASE)/src/main.o: src/main.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/main.c -o $(OBJDIR_RELEASE)/src/main.o

$(OBJDIR_RELEASE)/src/marking.o: src/marking.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/marking.c -o $(OBJDIR_RELEASE)/src/marking.o

$(OBJDIR_RELEASE)/src/presenceNode.o: src/presenceNode.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/presenceNode.c -o $(OBJDIR_RELEASE)/src/presenceNode.o

$(OBJDIR_RELEASE)/src/printMemory.o: src/printMemory.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/printMemory.c -o $(OBJDIR_RELEASE)/src/printMemory.o

$(OBJDIR_RELEASE)/src/quicksort.o: src/quicksort.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/quicksort.c -o $(OBJDIR_RELEASE)/src/quicksort.o

$(OBJDIR_RELEASE)/src/replaceAnnotation.o: src/replaceAnnotation.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/replaceAnnotation.c -o $(OBJDIR_RELEASE)/src/replaceAnnotation.o

$(OBJDIR_RELEASE)/src/retrieveAnnotation.o: src/retrieveAnnotation.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/retrieveAnnotation.c -o $(OBJDIR_RELEASE)/src/retrieveAnnotation.o

$(OBJDIR_RELEASE)/src/write_to_disk.o: src/write_to_disk.c
	$(CC) $(CFLAGS_RELEASE) $(INC_RELEASE) -c src/write_to_disk.c -o $(OBJDIR_RELEASE)/src/write_to_disk.o

clean_release: 
	rm -f $(OBJ_RELEASE) $(OUT_RELEASE)
	rm -rf bin/Release
	rm -rf $(OBJDIR_RELEASE)/src

.PHONY: before_debug after_debug clean_debug before_release after_release clean_release

