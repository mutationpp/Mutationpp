DIR=/home/jb/src/mutation.1.3
THSRC=$(DIR)/sources/thermo
TRSRC=$(DIR)/sources/transport
CHSRC=$(DIR)/sources/chemistry
GESRC=$(DIR)/sources/general
NUSRC=$(DIR)/sources/numerics
EXEC=$(DIR)/exec
OBJ=$(DIR)/sources/objects
LIB=$(DIR)/sources/objects

# FORTRAN-77 compiler given below
FC_i386=g77
FC_alpha=f90
FC_sparc=f90
FC_x86_64=g77
FC_powerpc=g77
FC=$(FC_$(MACHTYPE))
FC=gfortran
echo "FC=$FC"

LD_i386=g77
LD_alpha=f90
LD_sparc=f90
LD_x86_64=g77
LD_powerpc=g77
LD=$(LD_$(MACHTYPE))
LD=gfortran
echo "LD=$LD"

# COMPILER FLAGS GIVEN BELOW
COMPILERFLAGS_i386= 
#COMPILERFLAGS_i386= -O -m0 -YEXT_NAMES=LCS
COMPILERFLAGS_alpha=
COMPILERFLAGS_sparc=
COMPILERFLAGS_x86_64=
COMPILERFLAGS_powerpc= -fbounds-check 
COMPILERFLAGS_powerpc=
FFLAGS=$(COMPILERFLAGS_$(MACHTYPE))
FFLAGS= -g -O1

#LINKER FLAGS GIVEN BELOW
LINKFLAGS_i386= 
#LINKFLAGS_i386= -O -B101 -en -m0 -YEXT_NAMES=LCS -B108 -YCFRL=1 /usr/lib/libc.a
LINKFLAGS_alpha=
LINKFLAGS_sparc=
LINKFLAGS_x86_64= 
LINKFLAGS_powerpc= -O 
LDFLAGS=$(LINKFLAGS_$(MACHTYPE))
LDFLAGS=

# DEFINES
AR=ar
ARFLAGS=rv
CP=cp -f
RM=rm -f
RNLIB = touch

OBJS=\
 $(OBJ)/array.o\
 $(OBJ)/arrhenius.o\
 $(OBJ)/cg.o\
 $(OBJ)/cholesky.o\
 $(OBJ)/collision.o\
 $(OBJ)/correction.o\
 $(OBJ)/energy.o\
 $(OBJ)/enthalpy.o\
 $(OBJ)/equilibrium.o\
 $(OBJ)/density.o\
 $(OBJ)/entropy.o\
 $(OBJ)/eta.o\
 $(OBJ)/gamma.o\
 $(OBJ)/gibbs.o\
 $(OBJ)/gmres.o\
 $(OBJ)/initialize.o\
 $(OBJ)/lambda.o\
 $(OBJ)/lambdae.o\
 $(OBJ)/lambdaint.o\
 $(OBJ)/lambdarea.o\
 $(OBJ)/length.o\
 $(OBJ)/link.o\
 $(OBJ)/lu.o\
 $(OBJ)/minres.o\
 $(OBJ)/mixproperty.o\
 $(OBJ)/operation.o\
 $(OBJ)/presentation.o\
 $(OBJ)/readinput.o\
 $(OBJ)/sigmae.o\
 $(OBJ)/stefan.o\
 $(OBJ)/total.o\
 $(OBJ)/transfer.o\

mutobj= \
 $(OBJ)/mutation.o 

1Dentryobj= \
 $(OBJ)/1Dentry.o

shockobj= \
 $(OBJ)/shock.o

interobj= \
 $(OBJ)/interface77.o

shockinoobj= \
 $(OBJ)/shockino.o

default: mutlib 1Dentry

$(OBJ)/array.o : $(GESRC)/array.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(GESRC)/array.f

$(OBJ)/arrhenius.o : $(CHSRC)/arrhenius.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(CHSRC)/arrhenius.f

$(OBJ)/cg.o : $(NUSRC)/cg.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(NUSRC)/cg.f

$(OBJ)/cholesky.o : $(NUSRC)/cholesky.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(NUSRC)/cholesky.f

$(OBJ)/collision.o : $(TRSRC)/collision.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(TRSRC)/collision.f

$(OBJ)/correction.o : $(TRSRC)/correction.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(TRSRC)/correction.f

$(OBJ)/energy.o : $(THSRC)/energy.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(THSRC)/energy.f

$(OBJ)/entropy.o : $(THSRC)/entropy.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(THSRC)/entropy.f

$(OBJ)/enthalpy.o : $(THSRC)/enthalpy.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(THSRC)/enthalpy.f

$(OBJ)/equilibrium.o : $(CHSRC)/equilibrium.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(CHSRC)/equilibrium.f

$(OBJ)/eta.o : $(TRSRC)/eta.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(TRSRC)/eta.f

$(OBJ)/density.o : $(THSRC)/density.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(THSRC)/density.f

$(OBJ)/gamma.o : $(THSRC)/gamma.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(THSRC)/gamma.f

$(OBJ)/gibbs.o : $(THSRC)/gibbs.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(THSRC)/gibbs.f

$(OBJ)/gmres.o : $(NUSRC)/gmres.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(NUSRC)/gmres.f

$(OBJ)/interface77.o : $(GESRC)/interface77.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(GESRC)/interface77.f

$(OBJ)/lambda.o : $(TRSRC)/lambda.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(TRSRC)/lambda.f

$(OBJ)/lambdae.o : $(TRSRC)/lambdae.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(TRSRC)/lambdae.f

$(OBJ)/lambdaint.o : $(TRSRC)/lambdaint.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(TRSRC)/lambdaint.f

$(OBJ)/lambdarea.o : $(TRSRC)/lambdarea.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(TRSRC)/lambdarea.f

$(OBJ)/initialize.o : $(GESRC)/initialize.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(GESRC)/initialize.f

$(OBJ)/length.o : $(GESRC)/length.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(GESRC)/length.f

$(OBJ)/link.o : $(NUSRC)/link.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(NUSRC)/link.f

$(OBJ)/mutation.o : $(GESRC)/mutation.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(GESRC)/mutation.f

$(OBJ)/lu.o : $(NUSRC)/lu.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(NUSRC)/lu.f

$(OBJ)/minres.o : $(NUSRC)/minres.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(NUSRC)/minres.f

$(OBJ)/mixproperty.o : $(THSRC)/mixproperty.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(THSRC)/mixproperty.f

$(OBJ)/operation.o : $(NUSRC)/operation.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(NUSRC)/operation.f

$(OBJ)/presentation.o : $(GESRC)/presentation.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(GESRC)/presentation.f

$(OBJ)/readinput.o : $(GESRC)/readinput.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(GESRC)/readinput.f

$(OBJ)/1Dentry.o : $(GESRC)/1Dentry.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(GESRC)/1Dentry.f

$(OBJ)/shock.o : $(GESRC)/shock.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(GESRC)/shock.f

$(OBJ)/shockino.o : $(GESRC)/shockino.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(GESRC)/shockino.f

$(OBJ)/sigmae.o : $(TRSRC)/sigmae.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(TRSRC)/sigmae.f

$(OBJ)/stefan.o : $(TRSRC)/stefan.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(TRSRC)/stefan.f

$(OBJ)/total.o : $(THSRC)/total.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(THSRC)/total.f

$(OBJ)/transfer.o : $(CHSRC)/transfer.f
	cd $(OBJ); $(FC) $(FFLAGS) -c $(CHSRC)/transfer.f

mutlib: $(OBJS)
	@echo Machine is $(MACHTYPE)
	@echo "Building library Mutation"
	@rm -f libmutation_$(MACHTYPE).a
	cd $(LIB); $(AR) ruv libmutation.a $(OBJS)
	@cd $(LIB); $(RNLIB) libmutation.a

mutation: $(OBJS) $(mutobj)
	@echo Machine is $(MACHTYPE)
	@echo "Building mutation"
	@$(LD) -o $(EXEC)/mutation $(OBJS) $(mutobj) $(LDFLAGS)

1Dentry: $(OBJS) $(1Dentryobj)
	@echo Machine is $(MACHTYPE)
	@echo "Building 1Dentry"
	@$(LD) -o $(EXEC)/1Dentry $(OBJS) $(1Dentryobj) $(LDFLAGS)

shock: $(OBJS) $(shockobj)
	@echo Machine is $(MACHTYPE)
	@echo "Building shock"
	@$(LD) -o $(EXEC)/shock $(OBJS) $(shockobj) $(LDFLAGS)

interface: $(OBJS) $(interobj)
	@echo Machine is $(MACHTYPE)
	@echo "Building interface"
	@$(LD) -o $(EXEC)/interface77 $(OBJS) $(interobj) $(LDFLAGS)

shockino: $(OBJS) $(shockinoobj)
	@echo Machine is $(MACHTYPE)
	@echo "Building shockino"
	@$(LD) -o $(EXEC)/shockino $(OBJS) $(shockinoobj) $(LDFLAGS)

clean:
	cd $(OBJ); $(RM)  *.o
	cd $(EXEC); $(RM) libmutation.a mutation 1Dentry shock interface shockino
