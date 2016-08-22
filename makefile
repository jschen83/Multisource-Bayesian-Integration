#################################################################
#Directory
WORK         =.
SUBFUN       =.

#Compiler
COMPILE.C    =  g++ -c
COMPILE.PLUS =  g++ -c
COMPILE.F    =  g77 -c 
LINK.F       =  g77 -O2  
LINK.PLUS    =  g++ -v  
FFLAGS       = 
LIBS         = -lm 
#LIBS= -lm -lgcc -lg2c
#################################################################
#Main-functions
MAINSTTR    =$(WORK)/main_sttr.o

#Subfunctions
PRIOR       =$(WORK)/class_prior.o 
LKHD        =$(WORK)/class_lkhd.o 
POST        =$(WORK)/class_post.o 


#Sub-General-functions
BASE        = $(SUBFUN)/basefuns.o

#################################################################
mainOBJ=  $(MAINSTTR)   \
          $(POST) $(LKHD)  $(PRIOR) \
          $(BASE) 

#################################################################
INCLUDE   = $(WORK)/bayes.h 

#################################################################
# Link objective functions to  obtain executive files
#################################################################
mainsttr: $(INCLUDE)  $(mainOBJ) 
	$(LINK.PLUS) $(FFLAGS) $(mainOBJ) -o $(WORK)/mainsttr $(LIBS)

##################################################################
# Compile all other sub-functions
##################################################################
$(MAINSTTR): $(WORK)/main_sttr.cpp
	$(COMPILE.PLUS) $(FFLAGS) $(WORK)/main_sttr.cpp -o $@

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
	$(COMPILE.PLUS) $(FFLAGS) $(WORK)/class_prior.cpp -o $@
$(LKHD): $(INCLUDE)   $(WORK)/class_lkhd.cpp 
	$(COMPILE.PLUS) $(FFLAGS) $(WORK)/class_lkhd.cpp -o $@
$(POST): $(INCLUDE)  $(WORK)/class_post.cpp 
	$(COMPILE.PLUS) $(FFLAGS) $(WORK)/class_post.cpp -o $@
$(BASE): $(INCLUDE) $(SUBFUN)/basefuns.c 
	$(COMPILE.C) $(FFLAGS) $(SUBFUN)/basefuns.c -o $@
