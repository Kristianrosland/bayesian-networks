#-----------------------------------------------------------------------
#-------------------------- GOBNILP MAKE FILE --------------------------
#-----------------------------------------------------------------------

# This file was created by editing the Makefile for the linear ordering
# example in SCIP

# Set the path the SCIP directory
SCIPDIR = scip

# Include default project Makefile from SCIP
include $(SCIPDIR)/make/make.project

# This missing for scip < 3
VALGRIND = false

# Declare the phony targets
.PHONY: all doxygen manual pedmanual gobnilp clean test

# This stops the intermediate latex files being deleted
.SECONDARY:

#-----------------------------------------------------------------------
#----------------------- Command Line Arguments  -----------------------
#-----------------------------------------------------------------------

# HUSH is used to optionally hide output from latex and Doxygen
# (Turned on by default)
HUSH 				=	TRUE
ifeq ($(HUSH), TRUE)
	HUSHOUT = &>/dev/null
endif

# STATIC is used to statically link the included libraries
# (Turned off by default)
STATIC			=	FALSE
ifeq ($(STATIC), TRUE)
	STATICFLAG = -static
endif

#-----------------------------------------------------------------------
#---------------------------- Declarations  ----------------------------
#-----------------------------------------------------------------------

MAINNAME			=	gobnilp
MAINSRC			=	$(wildcard $(SRCDIR)/*.c)
MAINOBJ			=	$(notdir $(MAINSRC:.c=.o))

MAIN				=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE			=	$(BINDIR)/$(MAIN)
ifeq ($(OSTYPE),mingw)
	MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME).exe
else
	MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
endif
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

DOXYGENFILE		=	Doxyfile
DOXYGENDIR		=	docs

MANUALNAME		= 	manual
PEDMANUALNAME	=	pedigree_manual
LATEXCMD			=	@pdflatex -halt-on-error -interaction=errorstopmode

#Needs to appear first so that "make" makes gobnilp
gobnilp: $(SCIP) $(MAINFILE) $(MAINSHORTLINK)

#-----------------------------------------------------------------------
#----------------------------      All      ----------------------------
#-----------------------------------------------------------------------

all: gobnilp doxygen manual pedmanual


#-----------------------------------------------------------------------
#----------------------------    Doxygen    ----------------------------
#-----------------------------------------------------------------------

doxygen: $(DOXYGENDIR)

$(DOXYGENDIR): $(DOXYGENFILE) $(MAINSRC)
		@echo "Making Doxygen files"
		@doxygen $(HUSHOUT)

#-----------------------------------------------------------------------
#----------------------------    Manuals    ----------------------------
#-----------------------------------------------------------------------

manual: $(MANUALNAME).pdf
pedmanual: $(PEDMANUALNAME).pdf

%.pdf: %.ltx %.bbl
		@echo "Making manual"
		@$(LATEXCMD) $*.ltx $(HUSHOUT)
		@$(LATEXCMD) $*.ltx $(HUSHOUT)
		@$(LATEXCMD) $*.ltx $(HUSHOUT)
		@rm -f $*.aux $*.log $*.toc

%.bbl: %.bib %_tmp.aux
		@bibtex $*_tmp $(HUSHOUT)
		@mv $*_tmp.bbl $*.bbl
		@rm -f $*_tmp.blg

%_tmp.aux: %.ltx
		@$(LATEXCMD) -jobname $*_tmp $*.ltx $(HUSHOUT)
		@rm -f $*_tmp.log $*_tmp.pdf $*_tmp.toc

#-----------------------------------------------------------------------
#----------------------------     Clean     ----------------------------
#-----------------------------------------------------------------------

clean:
		@echo "Removing object files"
		@if [ -d obj ] ; then                                            \
			for i in $$(ls obj) ; do                                      \
				rm -f  obj/$$i/*.o ;                                       \
				rmdir  obj/$$i ;                                           \
			done ;                                                        \
			rmdir obj ;                                                   \
		fi

		@echo "Removing binaries"
		@rm -f $(BINDIR)/$(MAINNAME)*
		@if [ -d $(BINDIR) ] ; then rmdir $(BINDIR) ; fi

		@echo "Removing documentation"
		@if [ -d $(DOXYGENDIR) ] ; then                                  \
			if [ -d $(DOXYGENDIR)/html ] ; then                           \
				if [ -d $(DOXYGENDIR)/html/search ] ; then                 \
					rm -f $(DOXYGENDIR)/html/search/* ;                     \
					rmdir $(DOXYGENDIR)/html/search ;                       \
				fi ;                                                       \
				rm -f $(DOXYGENDIR)/html/* ;                               \
				rmdir $(DOXYGENDIR)/html ;                                 \
			fi ;                                                          \
			rmdir $(DOXYGENDIR) ;                                         \
		fi

		@echo "Removing manual"
		@rm -f $(MANUALNAME)_tmp.aux
		@rm -f $(MANUALNAME).bbl
		@rm -f $(MANUALNAME).pdf
		@rm -f $(PEDMANUALNAME)_tmp.aux
		@rm -f $(PEDMANUALNAME).bbl
		@rm -f $(PEDMANUALNAME).pdf

#-----------------------------------------------------------------------
#----------------------------     Test      ----------------------------
#-----------------------------------------------------------------------

test:           $(MAINFILE)
		@-(cd check && ln -fs ../$(SCIPDIR)/check/evalcheck.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/check.awk);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/getlastprob.awk);
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) "example" $(LPS) $(VALGRIND);

#-----------------------------------------------------------------------
#----------------------------    GOBNILP    ----------------------------
#-----------------------------------------------------------------------

$(MAINSHORTLINK): $(MAINFILE)
		@rm -f $@
		@cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(MAINFILE): $(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
		@$(LINKCXX) $(STATICFLAG) $(MAINOBJFILES) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
                $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o: $(SRCDIR)/%.c
		@echo "-> compiling $*.c"
		@$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c -g $< $(CC_o)$@

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

#-----------------------------------------------------------------------
