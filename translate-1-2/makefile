headers = config.h channel.c node.h mesh.h model.h model_complete.h \
	model_zero_inertia.h model_diffusive.h model_kinematic.h \
	model_complete_LaxFriedrichs.h model_zero_inertia_LaxFriedrichs.h  \
	model_diffusive_upwind.h model_kinematic_upwind.h \
	model_complete_upwind.h model_zero_inertia_upwind.h  \
	# model_diffusive_LaxFriedrichs.h model_kinematic_LaxFriedrichs.h

sources = main.c channel.c node.c mesh.c model.c model_complete.c \
	model_zero_inertia.c model_diffusive.c model_kinematic.c \
	model_complete_LaxFriedrichs.c model_zero_inertia_LaxFriedrichs.c  \
	model_diffusive_upwind.c model_kinematic_upwind.c \
	model_complete_upwind.c model_zero_inertia_upwind.c \
	# model_diffusive_LaxFriedrichs.c model_kinematic_LaxFriedrichs.c

objects = main.o channel.o node.o mesh.o model.o model_complete.o \
	model_zero_inertia.o model_diffusive.o model_kinematic.o \
	model_complete_LaxFriedrichs.o model_zero_inertia_LaxFriedrichs.o  \
	model_diffusive_upwind.o model_kinematic_upwind.o \
	model_complete_upwind.o model_zero_inertia_upwind.o \
	# model_diffusive_LaxFriedrichs.o model_kinematic_LaxFriedrichs.o

manuals = reference-manual.pdf swocs-manuals/english/user-manual.pdf \
	swocs-manuals/español/manual-usuario.pdf

libraries = -lm

flags = -O2 -Wall

#compiler = gcc -c $(flags)
compiler = i586-mingw32msvc-gcc -c $(flags)

#linker = gcc $(flags)
linker = i586-mingw32msvc-gcc $(flags)

#swocs = swocs
swocs = swocs.exe

#translate-1-2 = translate-1-2
translate-1-2 = translate-1-2.exe

all: $(swocs) $(translate-1-2)

$(translate-1-2): translate-1-2.c makefile
	$(linker) translate-1-2.c -o $(translate-1-2)

$(swocs): $(objects) makefile
	$(linker) $(objects) $(libraries) -o $(swocs)

channel.o: channel.c channel.h config.h makefile
	$(compiler) channel.c -o channel.o

node.o: node.c node.h config.h makefile
	$(compiler) node.c -o node.o

mesh.o: mesh.c mesh.h node.h channel.h config.h makefile
	$(compiler) mesh.c -o mesh.o

model.o: model.c model.h mesh.h node.h channel.h config.h makefile
	$(compiler) model.c -o model.o

model_complete.o: model_complete.c model_complete.h model.h node.h channel.h \
	config.h makefile
	$(compiler) model_complete.c -o model_complete.o

model_zero_inertia.o: model_zero_inertia.c model_zero_inertia.h model.h node.h \
	channel.h config.h makefile
	$(compiler) model_zero_inertia.c -o model_zero_inertia.o

model_diffusive.o: model_diffusive.c model_diffusive.h model.h node.h \
	channel.h config.h makefile
	$(compiler) model_diffusive.c -o model_diffusive.o

model_kinematic.o: model_kinematic.c model_kinematic.h model.h node.h \
	channel.h config.h makefile
	$(compiler) model_kinematic.c -o model_kinematic.o

model_complete_LaxFriedrichs.o: model_complete_LaxFriedrichs.c \
	model_complete_LaxFriedrichs.h model.h node.h channel.h config.h makefile
	$(compiler) model_complete_LaxFriedrichs.c -o model_complete_LaxFriedrichs.o

model_zero_inertia_LaxFriedrichs.o: model_zero_inertia_LaxFriedrichs.c \
	model_zero_inertia_LaxFriedrichs.h model.h node.h channel.h config.h \
	makefile
	$(compiler) model_zero_inertia_LaxFriedrichs.c \
		-o model_zero_inertia_LaxFriedrichs.o

model_diffusive_upwind.o: model_diffusive_upwind.c model_diffusive_upwind.h \
	model.h node.h channel.h config.h makefile
	$(compiler) model_diffusive_upwind.c -o model_diffusive_upwind.o

model_kinematic_upwind.o: model_kinematic_upwind.c model_kinematic_upwind.h \
	model.h node.h channel.h config.h makefile
	$(compiler) model_kinematic_upwind.c -o model_kinematic_upwind.o

model_complete_upwind.o: model_complete_upwind.c model_complete_upwind.h \
	model.h node.h channel.h config.h makefile
	$(compiler) model_complete_upwind.c -o model_complete_upwind.o

model_zero_inertia_upwind.o: model_zero_inertia_upwind.c \
	model_zero_inertia_upwind.h model.h node.h channel.h config.h makefile
	$(compiler) model_zero_inertia_upwind.c -o model_zero_inertia_upwind.o

#model_diffusive_LaxFriedrichs.o: model_diffusive_LaxFriedrichs.c \
	model_diffusive_LaxFriedrichs.h model.h node.h channel.h config.h makefile
#	$(compiler) model_diffusive_LaxFriedrichs.c \
		-o model_diffusive_LaxFriedrichs.o

#model_kinematic_LaxFriedrichs.o: model_kinematic_LaxFriedrichs.c \
	model_kinematic_LaxFriedrichs.h model.h node.h channel.h config.h makefile
#	$(compiler) model_kinematic_LaxFriedrichs.c \
		-o model_kinematic_LaxFriedrichs.o

main.o: main.c model.h mesh.h node.h channel.h config.h makefile \
	model_kinematic.h model_diffusive.h model_zero_inertia.h model_complete.h \
	model_complete_LaxFriedrichs.h model_zero_inertia_LaxFriedrichs.h \
	model_kinematic_upwind.h model_diffusive_upwind.h \
	model_zero_inertia_upwind.h model_complete_upwind.h
	$(compiler) main.c -o main.o

manuals: $(manuals)

reference-manual.pdf: $(headers) $(sources) Doxyfile makefile
	doxygen
	cd latex; make; mv refman.pdf ../reference-manual.pdf

swocs-manuals/english/user-manual.pdf: swocs-manuals/english/user-manual.tex \
	makefile
	cd swocs-manuals/english; pdflatex user-manual; pdflatex user-manual

swocs-manuals/español/manual-usuario.pdf: \
	swocs-manuals/español/manual-usuario.tex \
	makefile
	cd swocs-manuals/español; pdflatex manual-usuario; pdflatex manual-usuario

zip: $(src) makefile
	zip swocs-src-1-2 *.h *.c makefile Doxyfile README

examples:

