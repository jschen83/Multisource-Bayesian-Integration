#################################################################
#Directory
WORK         =.
SUBFUN       =.

#Compiler
COMPILE.C    =  g++ -c
COMPILE.PLUS =  g++ -c
LINK.PLUS    =  g++ -v  
#
CFLAGS       = 
LIBS         = -lm 
#################################################################
#Main-functions
MAINDRIVER = $(WORK)/main_driver.o

#Subfunctions
PRIOR      = $(WORK)/class_prior.o 
LKHD       = $(WORK)/class_lkhd.o 
POST       = $(WORK)/class_post.o 

#Sub-General-functions
BASE        = $(SUBFUN)/basefuns.o

INCLUDE   = $(WORK)/bayes.h 

#################################################################
mainOBJ=  $(MAINDRIVER) $(POST) $(LKHD)  $(PRIOR) $(BASE) 

mymain: $(INCLUDE)  $(mainOBJ) 
	$(LINK.PLUS) $(CFLAGS) $(mainOBJ) -o $(WORK)/mymain $(LIBS)

##################################################################
# Compile all other sub-functions
##################################################################
$(MAINDRIVER): $(WORK)/main_driver.cpp
	$(COMPILE.PLUS) $(CFLAGS) $(WORK)/main_driver.cpp -o $@

###################################################################
# Clean all objective functions
###################################################################
allOBJ =$(mainOBJ)
clean:
	rm $(allOBJ)

###################################################################
# Compile other sub-functions
###################################################################
$(PRIOR): $(INCLUDE)   $(WORK)/class_prior.cpp 
	$(COMPILE.PLUS) $(CFLAGS) $(WORK)/class_prior.cpp -o $@
$(LKHD): $(INCLUDE)   $(WORK)/class_lkhd.cpp 
	$(COMPILE.PLUS) $(CFLAGS) $(WORK)/class_lkhd.cpp -o $@
$(POST): $(INCLUDE)  $(WORK)/class_post.cpp 
	$(COMPILE.PLUS) $(CFLAGS) $(WORK)/class_post.cpp -o $@
$(BASE): $(INCLUDE) $(SUBFUN)/basefuns.c 
	$(COMPILE.C) $(CFLAGS) $(SUBFUN)/basefuns.c -o $@
